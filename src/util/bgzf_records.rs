//! A BGZF pipeline for the records of a BAM file, in which each worker that inflates a block
//! also reads the records of that block.
//!
//! One thread reads the raw BGZF frames. A pool of workers inflates them (and checks their
//! CRC) in parallel; each worker then takes, in block order, the *carry* the previous block's
//! worker handed over (the bytes of a record cut by the block end), walks the record chain of
//! its own block to hand the next carry on, and only then calls `work` on every record that
//! ends in its block, while the block is still in its cache. The calling thread receives the
//! per-block outputs and calls `stitch` on them in block order.
//!
//! The records, their order and their fields are exactly what [`noodles_bam::io::Reader`]
//! reads from the same stream: the records are read with that reader, over the record bytes.
//! As in that reader, a record block size of 0 ends the stream, and a stream that ends inside
//! a record is an error.

use std::fs::File;
use std::io::{self, BufReader, Read};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};

use crossbeam::channel;
use noodles_bam as bam;

/// What a worker builds from the records of one block.
pub trait BlockOutput: Default + Send {
    /// Empties the output for reuse, keeping its allocations.
    fn clear(&mut self);
}

const BGZF_HEADER_SIZE: usize = 18;
const GZ_TRAILER_SIZE: usize = 8;
const MIN_FRAME_SIZE: usize = BGZF_HEADER_SIZE + GZ_TRAILER_SIZE;

enum Frame {
    Data(Vec<u8>),
    /// the stream of frames ended (cleanly) before this block
    End,
    Error(io::Error),
}

struct Job {
    idx: u64,
    frame: Frame,
}

struct Done<O> {
    idx: u64,
    out: io::Result<O>,
    last: bool,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ChainState {
    /// at a record boundary, or inside a record whose start bytes are pending
    Live,
    /// a record block size of 0 ended the stream
    Stopped,
    /// an earlier block failed; its error ends the parse
    Failed,
}

struct Carry {
    /// the bytes of a record that starts in an earlier block and does not end there
    pending: Vec<u8>,
    state: ChainState,
}

/// Hands the carry from block to block, in block order. Only one carry exists at a time (the
/// carry into block `i + 1` is made by the worker that took the carry into block `i`), so one
/// slot holds it; a waiting worker spins briefly on `turn` and then sleeps.
struct Chain {
    turn: AtomicU64,
    slot: Mutex<Option<Carry>>,
    cv: Condvar,
    sleepers: AtomicUsize,
    abort: AtomicBool,
}

/// How many times a worker polls `Chain::turn` before it sleeps (4096 and 65536 parse alike;
/// 512 is slower).
const SPINS: u32 = 4096;

impl Chain {
    fn new() -> Self {
        Self {
            turn: AtomicU64::new(0),
            slot: Mutex::new(Some(Carry {
                pending: Vec::new(),
                state: ChainState::Live,
            })),
            cv: Condvar::new(),
            sleepers: AtomicUsize::new(0),
            abort: AtomicBool::new(false),
        }
    }

    /// Waits for the carry into block `idx`; `None` once the pipeline is aborted.
    fn take(&self, idx: u64) -> Option<Carry> {
        let mut spins = 0_u32;
        loop {
            if self.turn.load(Ordering::Acquire) == idx {
                break;
            }
            if self.abort.load(Ordering::Acquire) {
                return None;
            }
            spins += 1;
            if spins < SPINS {
                std::hint::spin_loop();
            } else {
                let mut g = self.slot.lock().unwrap();
                self.sleepers.fetch_add(1, Ordering::SeqCst);
                while self.turn.load(Ordering::SeqCst) != idx && !self.abort.load(Ordering::SeqCst)
                {
                    g = self.cv.wait(g).unwrap();
                }
                self.sleepers.fetch_sub(1, Ordering::SeqCst);
                if self.turn.load(Ordering::SeqCst) != idx {
                    return None;
                }
                return Some(g.take().expect("the carry for this block"));
            }
        }
        Some(
            self.slot
                .lock()
                .unwrap()
                .take()
                .expect("the carry for this block"),
        )
    }

    /// Hands `carry` on to block `idx`.
    fn put(&self, idx: u64, carry: Carry) {
        {
            let mut g = self.slot.lock().unwrap();
            *g = Some(carry);
            self.turn.store(idx, Ordering::SeqCst);
        }
        if self.sleepers.load(Ordering::SeqCst) > 0 {
            self.cv.notify_all();
        }
    }

    fn abort(&self) {
        self.abort.store(true, Ordering::SeqCst);
        let _g = self.slot.lock().unwrap();
        self.cv.notify_all();
    }
}

/// Reads one BGZF frame into `buf`, as noodles' BGZF readers do: `false` at the end of the
/// stream (including a partial frame header there).
fn read_frame_into<R: Read>(reader: &mut R, buf: &mut Vec<u8>) -> io::Result<bool> {
    buf.resize(BGZF_HEADER_SIZE, 0);
    match reader.read_exact(buf) {
        Ok(()) => {}
        Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(false),
        Err(e) => return Err(e),
    }
    let bsize = u16::from_le_bytes([buf[16], buf[17]]);
    let block_size = usize::from(bsize) + 1;
    if block_size < MIN_FRAME_SIZE {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid frame size",
        ));
    }
    buf.resize(block_size, 0);
    reader.read_exact(&mut buf[BGZF_HEADER_SIZE..])?;
    Ok(true)
}

fn is_valid_header(src: &[u8]) -> bool {
    src[0..2] == [0x1f, 0x8b]
        && src[2] == 0x08
        && src[3] == 0x04
        && src[10..12] == [0x06, 0x00]
        && src[12..14] == *b"BC"
        && src[14..16] == [0x02, 0x00]
}

/// Inflates one frame into `data` (which keeps its length across blocks, so it is not cleared
/// for each block) and checks its CRC; returns the inflated length.
fn inflate_frame(
    frame: &[u8],
    inflater: &mut flate2::Decompress,
    data: &mut Vec<u8>,
) -> io::Result<usize> {
    if frame.len() < MIN_FRAME_SIZE || !is_valid_header(&frame[..BGZF_HEADER_SIZE]) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF header",
        ));
    }
    let end = frame.len() - GZ_TRAILER_SIZE;
    let cdata = &frame[BGZF_HEADER_SIZE..end];
    let crc32 = u32::from_le_bytes(frame[end..end + 4].try_into().unwrap());
    let isize = u32::from_le_bytes(frame[end + 4..].try_into().unwrap()) as usize;
    if data.len() < isize {
        data.resize(isize, 0);
    }
    let dst = &mut data[..isize];
    if isize > 0 {
        inflater.reset(false);
        inflater
            .decompress(cdata, dst, flate2::FlushDecompress::Finish)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if inflater.total_out() != isize as u64 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "failed to fill whole buffer",
            ));
        }
    }
    let mut crc = flate2::Crc::new();
    crc.update(dst);
    if crc.sum() != crc32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "block data checksum mismatch",
        ));
    }
    Ok(isize)
}

fn block_size_at(src: &[u8]) -> usize {
    u32::from_le_bytes(src[..4].try_into().unwrap()) as usize
}

struct Worker {
    inflater: flate2::Decompress,
    spare: Vec<u8>,
    rec: bam::Record,
    first_skip: usize,
}

impl Worker {
    /// Takes block `idx` (inflated into `data[..n]`) through the record chain with the carry
    /// into it, hands the carry on, then reads the block's records into `out`.
    #[allow(clippy::too_many_arguments)]
    fn process<O, W>(
        &mut self,
        idx: u64,
        last: bool,
        inflated: io::Result<usize>,
        data: &[u8],
        mut carry: Carry,
        mut out: O,
        chain: &Chain,
        work: &W,
    ) -> Done<O>
    where
        O: BlockOutput,
        W: Fn(&bam::Record, &mut O),
    {
        if carry.state == ChainState::Failed {
            chain.put(idx + 1, carry);
            out.clear();
            return Done {
                idx,
                out: Ok(out),
                last,
            };
        }
        let n = match inflated {
            Ok(n) => n,
            Err(e) => {
                carry.state = ChainState::Failed;
                chain.put(idx + 1, carry);
                return Done {
                    idx,
                    out: Err(e),
                    last,
                };
            }
        };
        if last {
            let res = if carry.state == ChainState::Live && !carry.pending.is_empty() {
                Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "failed to fill whole buffer",
                ))
            } else {
                out.clear();
                Ok(out)
            };
            chain.put(idx + 1, carry);
            return Done {
                idx,
                out: res,
                last,
            };
        }
        let d = &data[..n];
        let mut off = if idx == 0 { self.first_skip.min(n) } else { 0 };
        // the record cut by the previous block end, once complete
        let mut straddle: Option<Vec<u8>> = None;
        if carry.state == ChainState::Live && !carry.pending.is_empty() {
            let p = &mut carry.pending;
            if p.len() < 4 {
                let take = (4 - p.len()).min(n - off);
                p.extend_from_slice(&d[off..off + take]);
                off += take;
            }
            if p.len() >= 4 {
                let bs = block_size_at(p);
                if bs == 0 {
                    carry.state = ChainState::Stopped;
                    p.clear();
                } else {
                    let take = (4 + bs - p.len()).min(n - off);
                    p.extend_from_slice(&d[off..off + take]);
                    off += take;
                    if p.len() == 4 + bs {
                        self.spare.clear();
                        straddle = Some(std::mem::replace(p, std::mem::take(&mut self.spare)));
                    }
                }
            }
        }
        // the records that start and end in this block are d[start..end]
        let start = off;
        let mut end = off;
        if carry.state == ChainState::Live && carry.pending.is_empty() {
            loop {
                if n - end < 4 {
                    break;
                }
                let bs = block_size_at(&d[end..]);
                if bs == 0 {
                    carry.state = ChainState::Stopped;
                    break;
                }
                if n - end - 4 < bs {
                    break;
                }
                end += 4 + bs;
            }
            if carry.state == ChainState::Live {
                carry.pending.extend_from_slice(&d[end..]);
            }
        }
        chain.put(idx + 1, carry);

        // the records, in parallel with the other workers
        out.clear();
        let mut res = Ok(());
        if let Some(s) = straddle {
            res = read_records(&s, &mut self.rec, &mut out, work);
            self.spare = s;
        }
        if res.is_ok() {
            res = read_records(&d[start..end], &mut self.rec, &mut out, work);
        }
        Done {
            idx,
            out: res.map(|_| out),
            last: false,
        }
    }
}

/// Calls `work` on each record of `src`, a whole number of length-prefixed records.
fn read_records<O, W>(src: &[u8], rec: &mut bam::Record, out: &mut O, work: &W) -> io::Result<()>
where
    W: Fn(&bam::Record, &mut O),
{
    let mut reader = bam::io::Reader::from(src);
    while reader.read_record(rec)? != 0 {
        work(rec, out);
    }
    Ok(())
}

/// Reads the BAM records of `input` (a BGZF stream positioned at the start of a block, whose
/// first `first_skip` inflated bytes are skipped), calling `work` on each record on one of
/// `workers` worker threads, into the output of the record's block, and `stitch` on the
/// outputs on the calling thread, in block order. `stitch` returning `Ok(false)` stops the
/// parse early (with `Ok`); an error from it or from the stream stops it with that error.
pub fn for_each_block<O, W, S>(
    input: File,
    first_skip: usize,
    workers: usize,
    work: W,
    mut stitch: S,
) -> anyhow::Result<()>
where
    O: BlockOutput,
    W: Fn(&bam::Record, &mut O) + Sync,
    S: FnMut(&mut O) -> anyhow::Result<bool>,
{
    let workers = workers.max(1);
    // blocks in flight (read, not yet stitched); bounds memory and the reorder ring
    let in_flight = 4 * workers + 8;

    let chain = Chain::new();
    let (job_tx, job_rx) = channel::bounded::<Job>(2 * workers);
    let (done_tx, done_rx) = channel::unbounded::<Done<O>>();
    let (credit_tx, credit_rx) = channel::bounded::<()>(in_flight);
    for _ in 0..in_flight {
        credit_tx.send(()).unwrap();
    }
    let (free_frame_tx, free_frame_rx) = channel::unbounded::<Vec<u8>>();
    let (free_out_tx, free_out_rx) = channel::unbounded::<O>();

    std::thread::scope(|scope| -> anyhow::Result<()> {
        // the frame reader
        scope.spawn(move || {
            let mut reader = BufReader::with_capacity(1 << 20, input);
            let mut idx = 0_u64;
            loop {
                if credit_rx.recv().is_err() {
                    return;
                }
                let mut buf = free_frame_rx.try_recv().unwrap_or_default();
                let (frame, last) = match read_frame_into(&mut reader, &mut buf) {
                    Ok(true) => (Frame::Data(buf), false),
                    Ok(false) => (Frame::End, true),
                    Err(e) => (Frame::Error(e), true),
                };
                if job_tx.send(Job { idx, frame }).is_err() || last {
                    return;
                }
                idx += 1;
            }
        });

        for _ in 0..workers {
            let job_rx = job_rx.clone();
            let done_tx = done_tx.clone();
            let free_frame_tx = free_frame_tx.clone();
            let free_out_rx = free_out_rx.clone();
            let chain = &chain;
            let work = &work;
            scope.spawn(move || {
                let mut w = Worker {
                    inflater: flate2::Decompress::new(false),
                    spare: Vec::new(),
                    rec: bam::Record::default(),
                    first_skip,
                };
                let mut data = Vec::<u8>::new();
                while let Ok(Job { idx, frame }) = job_rx.recv() {
                    // inflate, in parallel with the other workers
                    let last = !matches!(frame, Frame::Data(_));
                    let inflated = match frame {
                        Frame::Data(buf) => {
                            let r = inflate_frame(&buf, &mut w.inflater, &mut data);
                            let _ = free_frame_tx.send(buf);
                            r
                        }
                        Frame::End => Ok(0),
                        Frame::Error(e) => Err(e),
                    };
                    let Some(carry) = chain.take(idx) else {
                        return;
                    };
                    let out = free_out_rx.try_recv().unwrap_or_default();
                    let done = w.process(idx, last, inflated, &data, carry, out, chain, work);
                    let _ = done_tx.send(done);
                }
            });
        }
        drop(job_rx);
        drop(done_tx);
        drop(free_frame_tx);
        drop(free_out_rx);

        // the stitcher: this thread, in block order
        let result = (|| -> anyhow::Result<()> {
            let mut ring: Vec<Option<Done<O>>> = (0..in_flight).map(|_| None).collect();
            let mut next = 0_u64;
            loop {
                let slot = (next % in_flight as u64) as usize;
                while ring[slot].is_none() {
                    let d = done_rx
                        .recv()
                        .map_err(|_| anyhow::anyhow!("the BAM block workers stopped early"))?;
                    let s = (d.idx % in_flight as u64) as usize;
                    ring[s] = Some(d);
                }
                let d = ring[slot].take().unwrap();
                debug_assert_eq!(d.idx, next);
                let _ = credit_tx.send(());
                let mut out = d.out?;
                if d.last {
                    return Ok(());
                }
                if !stitch(&mut out)? {
                    return Ok(());
                }
                let _ = free_out_tx.send(out);
                next += 1;
            }
        })();
        // stop the reader and the workers (they have finished when the stream ended)
        drop(credit_tx);
        drop(done_rx);
        chain.abort();
        result
    })
}
