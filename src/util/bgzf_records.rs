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
//! reads from the same stream: each record is read in place from the inflated block (a
//! [`RawRecord`] over its bytes, which reads every field as `noodles_bam::Record` reads it),
//! and only a record cut by a block end is first joined into a buffer. As in that reader, a
//! record block size of 0 ends the stream, and a stream that ends inside a record is an
//! error.

use std::fs::File;
use std::io::{self, BufReader, Read};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};

use crossbeam::channel;

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
/// for each block) and checks its CRC; returns the inflated length. The raw-deflate payload is
/// inflated in one call through libdeflate, reusing the worker's decompressor.
fn inflate_frame(
    frame: &[u8],
    inflater: &mut libdeflater::Decompressor,
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
        let n = inflater
            .deflate_decompress(cdata, dst)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if n != isize {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "failed to fill whole buffer",
            ));
        }
    }
    if libdeflater::crc32(dst) != crc32 {
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
    inflater: libdeflater::Decompressor,
    spare: Vec<u8>,
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
        W: Fn(&RawRecord<'_>, &mut O),
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
            res = read_records(&s, &mut out, work);
            self.spare = s;
        }
        if res.is_ok() {
            res = read_records(&d[start..end], &mut out, work);
        }
        Done {
            idx,
            out: res.map(|_| out),
            last: false,
        }
    }
}

/// Calls `work` on each record of `src`, a whole number of length-prefixed records, read in
/// place, with the outcomes of `noodles_bam::io::Reader::read_record` over the same bytes: the
/// end of `src` or a block size of 0 ends the records, a cut length prefix or record is an
/// error.
fn read_records<O, W>(mut src: &[u8], out: &mut O, work: &W) -> io::Result<()>
where
    W: Fn(&RawRecord<'_>, &mut O),
{
    loop {
        if src.is_empty() {
            return Ok(());
        }
        if src.len() < 4 {
            return Err(cut_record());
        }
        let bs = block_size_at(src);
        if bs == 0 {
            return Ok(());
        }
        let rest = &src[4..];
        if rest.len() < bs {
            return Err(cut_record());
        }
        let (buf, next) = rest.split_at(bs);
        work(&RawRecord::new(buf)?, out);
        src = next;
    }
}

fn cut_record() -> io::Error {
    io::Error::new(io::ErrorKind::UnexpectedEof, "failed to fill whole buffer")
}

/// One BAM record, read in place from its bytes (those after its block size): each field is
/// read as `noodles_bam::Record` reads it from the same bytes, and only when asked for.
#[derive(Clone, Copy)]
pub struct RawRecord<'a> {
    buf: &'a [u8],
    name_end: usize,
    cigar_end: usize,
    sequence_end: usize,
    quality_scores_end: usize,
}

impl<'a> RawRecord<'a> {
    const NAME_START: usize = 32;

    /// Indexes the record's fields, failing as `noodles_bam` does on a record too short for
    /// its fixed fields or for the lengths they give.
    #[inline]
    pub fn new(buf: &'a [u8]) -> io::Result<Self> {
        if buf.len() < Self::NAME_START {
            return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
        }
        let name_end = Self::NAME_START + usize::from(buf[8]);
        let cigar_op_count = usize::from(u16::from_le_bytes([buf[12], buf[13]]));
        let cigar_end = name_end + 4 * cigar_op_count;
        let base_count = usize::try_from(u32::from_le_bytes(buf[16..20].try_into().unwrap()))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let sequence_end = cigar_end + base_count.div_ceil(2);
        let quality_scores_end = sequence_end + base_count;
        if buf.len() < quality_scores_end {
            return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
        }
        Ok(Self {
            buf,
            name_end,
            cigar_end,
            sequence_end,
            quality_scores_end,
        })
    }

    #[inline]
    fn i32_at(&self, i: usize) -> i32 {
        i32::from_le_bytes(self.buf[i..i + 4].try_into().unwrap())
    }

    /// The reference sequence id as stored (`None`: -1, no reference).
    #[inline]
    pub fn reference_sequence_id(&self) -> Option<i32> {
        match self.i32_at(0) {
            -1 => None,
            n => Some(n),
        }
    }

    /// The zero-based alignment start as stored (`None`: -1, no position).
    #[inline]
    pub fn alignment_start(&self) -> Option<i32> {
        match self.i32_at(4) {
            -1 => None,
            n => Some(n),
        }
    }

    #[inline]
    pub fn flags(&self) -> u16 {
        u16::from_le_bytes([self.buf[14], self.buf[15]])
    }

    /// The read name without its NUL (`None`: the missing name `*`).
    #[inline]
    pub fn name(&self) -> Option<&'a [u8]> {
        match &self.buf[Self::NAME_START..self.name_end] {
            b"*\0" => None,
            n => Some(n.strip_suffix(&[0]).unwrap_or(n)),
        }
    }

    /// The number of bases (`l_seq`).
    #[inline]
    pub fn sequence_len(&self) -> usize {
        self.quality_scores_end - self.sequence_end
    }

    /// The 4-bit code of base `i` (< `sequence_len()`).
    #[inline]
    pub fn base_code(&self, i: usize) -> u8 {
        let b = self.buf[self.cigar_end + i / 2];
        if i.is_multiple_of(2) {
            b >> 4
        } else {
            b & 0x0f
        }
    }

    /// The packed CIGAR ops (`len << 4 | kind`, 4 bytes each), as `noodles_bam::Record::cigar`
    /// gives them: a record whose stored CIGAR is the placeholder `<l_seq>S<n>N` and whose data
    /// holds a `CG:B` array has that array's bytes.
    #[inline]
    pub fn cigar(&self) -> &'a [u8] {
        const SKIP: u32 = 3;
        const SOFT_CLIP: u32 = 4;
        let src = &self.buf[self.name_end..self.cigar_end];
        if src.len() == 8 {
            let op_1 = u32::from_le_bytes(src[..4].try_into().unwrap());
            let op_2 = u32::from_le_bytes(src[4..].try_into().unwrap());
            if op_1 & 0x0f == SOFT_CLIP
                && (op_1 >> 4) as usize == self.sequence_len()
                && op_2 & 0x0f == SKIP
                && let Ok(Some(buf)) = raw_cigar_tag(self.data())
            {
                return buf;
            }
        }
        src
    }

    /// The optional fields' bytes.
    #[inline]
    pub fn data(&self) -> &'a [u8] {
        &self.buf[self.quality_scores_end..]
    }

    /// The integer value of tag `tag`, as `noodles_bam::record::Data::get` and `Value::as_int`
    /// give it: `None` when no field has the tag, `Some(Ok(None))` when its value is not an
    /// integer, and an error when a field up to it cannot be decoded.
    #[inline]
    pub fn int_tag(&self, tag: [u8; 2]) -> Option<io::Result<Option<i64>>> {
        let mut src = self.data();
        while !src.is_empty() {
            match decode_field(&mut src, tag) {
                Ok(Some(v)) => return Some(Ok(v)),
                Ok(None) => {}
                Err(e) => return Some(Err(e)),
            }
        }
        None
    }
}

fn take_bytes<'a>(src: &mut &'a [u8], n: usize) -> io::Result<&'a [u8]> {
    if src.len() < n {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }
    let (a, b) = src.split_at(n);
    *src = b;
    Ok(a)
}

/// The raw bytes of a `B` array value (after its type byte), as `noodles_bam` decodes them.
fn raw_array<'a>(src: &mut &'a [u8]) -> io::Result<&'a [u8]> {
    let width = match take_bytes(src, 1)?[0] {
        b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid subtype",
            ));
        }
    };
    let n = u32::from_le_bytes(take_bytes(src, 4)?.try_into().unwrap());
    let n = usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    take_bytes(src, n * width)
}

/// Skips a value of type `ty` as `noodles_bam`'s `decode_value` reads it (with its errors),
/// returning it as `Value::as_int` gives it.
#[inline]
fn decode_value(src: &mut &[u8], ty: u8) -> io::Result<Option<i64>> {
    Ok(match ty {
        b'A' => {
            take_bytes(src, 1)?;
            None
        }
        b'c' => Some(i64::from(take_bytes(src, 1)?[0] as i8)),
        b'C' => Some(i64::from(take_bytes(src, 1)?[0])),
        b's' => Some(i64::from(i16::from_le_bytes(
            take_bytes(src, 2)?.try_into().unwrap(),
        ))),
        b'S' => Some(i64::from(u16::from_le_bytes(
            take_bytes(src, 2)?.try_into().unwrap(),
        ))),
        b'i' => Some(i64::from(i32::from_le_bytes(
            take_bytes(src, 4)?.try_into().unwrap(),
        ))),
        b'I' => Some(i64::from(u32::from_le_bytes(
            take_bytes(src, 4)?.try_into().unwrap(),
        ))),
        b'f' => {
            take_bytes(src, 4)?;
            None
        }
        b'Z' | b'H' => {
            let len = src.iter().position(|&b| b == 0).ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "string not NUL terminated")
            })?;
            *src = &src[len + 1..];
            None
        }
        b'B' => {
            raw_array(src)?;
            None
        }
        _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid type")),
    })
}

/// Decodes one field as `noodles_bam`'s `decode_field` does; when its tag is `want`, returns
/// its value as `Value::as_int` gives it.
#[inline]
fn decode_field(src: &mut &[u8], want: [u8; 2]) -> io::Result<Option<Option<i64>>> {
    let tag = take_bytes(src, 2)?;
    let ty = take_bytes(src, 1)?[0];
    let v = decode_value(src, ty)?;
    Ok((tag == want).then_some(v))
}

/// The raw bytes of the first `CG:B` array in `src`, as `noodles_bam`'s `get_raw_cigar` finds
/// them (an error when a field before it cannot be decoded).
fn raw_cigar_tag(mut src: &[u8]) -> io::Result<Option<&[u8]>> {
    while !src.is_empty() {
        let tag = take_bytes(&mut src, 2)?;
        let ty = take_bytes(&mut src, 1)?[0];
        if ty == b'B' {
            let buf = raw_array(&mut src)?;
            if tag == b"CG" {
                return Ok(Some(buf));
            }
        } else {
            decode_value(&mut src, ty)?;
        }
    }
    Ok(None)
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
    W: Fn(&RawRecord<'_>, &mut O) + Sync,
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
                    inflater: libdeflater::Decompressor::new(),
                    spare: Vec::new(),
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
