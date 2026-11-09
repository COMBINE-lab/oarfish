//! A BGZF pipeline for the records of a BAM file, in which each worker that inflates a block
//! also reads the records of that block.
//!
//! One thread reads the raw BGZF frames. A pool of workers inflates them (and checks their
//! CRC) in parallel. The record chain runs in block order: the *carry* into a block (the bytes
//! of a record cut by the previous block end) tells where its first record starts. A worker
//! that finishes inflating a block whose carry is waiting walks the block's record chain to
//! make the next carry, then calls `work` on every record that ends in its block, while the
//! block is still in its cache. A worker whose carry has not arrived parks its block and takes
//! the next job instead of waiting; the carry holder walks the parked blocks that follow its
//! own and hands each back to its worker, which then reads its records. The calling thread
//! receives the per-block outputs and calls `stitch` on them in block order.
//!
//! The records, their order and their fields are exactly what [`noodles_bam::io::Reader`]
//! reads from the same stream: each record is read in place from the inflated block (a
//! [`RawRecord`] over its bytes, which reads every field as `noodles_bam::Record` reads it),
//! and only a record cut by a block end is first joined into a buffer. As in that reader, a
//! record block size of 0 ends the stream, and a stream that ends inside a record is an
//! error.

use std::fs::File;
use std::io::{self, BufReader, Read};
use std::sync::Mutex;

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

/// An inflated block whose carry had not arrived when its worker finished inflating it: the
/// worker leaves it here and takes the next job, and whoever holds the carry walks it.
struct Parked {
    idx: u64,
    /// the worker that inflated it, and reads its records once it is walked
    owner: usize,
    last: bool,
    inflated: io::Result<usize>,
    data: Vec<u8>,
    /// the block's record chain from a guessed first record start, walked while the block
    /// was in its worker's cache
    spec: Spec,
}

/// A record chain walked from a guessed start: `starts` holds the offsets the walk visits,
/// ascending, the last being where it ends (with `stopped`: a block size of 0 there). Empty
/// when no start was guessed. The chain from an offset is fixed by the bytes, so a walk from
/// the true start that reaches any of `starts` ends where this one does.
#[derive(Default)]
struct Spec {
    starts: Vec<u32>,
    stopped: bool,
}

/// The fixed fields of a BAM record (after its block size) as far as they bound its length.
const FIXED_FIELDS: usize = 32;

/// Whether `d[o..]` plausibly starts a record: a block size large enough for the fields its
/// fixed part gives, positions and reference ids of at least -1, and a NUL-terminated
/// printable read name. Only picks where [`speculate`] starts; nothing depends on it being
/// right.
fn plausible_record(d: &[u8], o: usize) -> Option<usize> {
    let r = d.get(o..)?;
    if r.len() < 4 + FIXED_FIELDS {
        return None;
    }
    let bs = block_size_at(r);
    if !(FIXED_FIELDS..1 << 24).contains(&bs) {
        return None;
    }
    let f = &r[4..];
    let i32_at = |i: usize| i32::from_le_bytes(f[i..i + 4].try_into().unwrap());
    if i32_at(0) < -1 || i32_at(4) < -1 || i32_at(20) < -1 || i32_at(24) < -1 {
        return None;
    }
    let l_name = usize::from(f[8]);
    let n_cigar = usize::from(u16::from_le_bytes([f[12], f[13]]));
    let l_seq = u32::from_le_bytes(f[16..20].try_into().unwrap()) as usize;
    if l_name == 0
        || FIXED_FIELDS + l_name + 4 * n_cigar + l_seq.div_ceil(2) + l_seq > bs
        || f.len() < FIXED_FIELDS + l_name
    {
        return None;
    }
    let name = &f[FIXED_FIELDS..FIXED_FIELDS + l_name];
    if name[l_name - 1] != 0 || !name[..l_name - 1].iter().all(|b| (b'!'..=b'~').contains(b)) {
        return None;
    }
    Some(4 + bs)
}

/// Guesses where the first record that starts in `d` starts (the first offset that plausibly
/// starts a record followed by another, or by the block end), and walks the record chain
/// from there into `spec`, as [`Worker::walk`] walks it.
fn speculate(d: &[u8], spec: &mut Spec) {
    spec.starts.clear();
    spec.stopped = false;
    let n = d.len();
    let Some(mut end) = (0..n).find(|&o| match plausible_record(d, o) {
        Some(len) => o + len + 4 + FIXED_FIELDS > n || plausible_record(d, o + len).is_some(),
        None => false,
    }) else {
        return;
    };
    loop {
        spec.starts.push(end as u32);
        if n - end < 4 {
            break;
        }
        let bs = block_size_at(&d[end..]);
        if bs == 0 {
            spec.stopped = true;
            break;
        }
        if n - end - 4 < bs {
            break;
        }
        end += 4 + bs;
    }
}

/// What reading a walked block's records takes.
enum Plan {
    /// no records: an empty output, or the block's error
    Empty(io::Result<()>),
    /// the record cut by the previous block end (once complete), then `data[start..end]`
    Read {
        straddle: Option<Vec<u8>>,
        start: usize,
        end: usize,
    },
}

/// Work handed to a worker by another.
enum Task {
    /// the carry holder walked this worker's parked block `idx`: read its records
    Read {
        idx: u64,
        last: bool,
        data: Vec<u8>,
        plan: Plan,
        spec: Spec,
    },
    /// the parse has stopped
    Abort,
}

struct ChainCore {
    /// the block the carry goes into next
    next: u64,
    /// the carry into `next`, while no worker holds it
    carry: Option<Carry>,
    /// parked blocks, at `idx % parked.len()` (fewer blocks than that are in flight)
    parked: Vec<Option<Parked>>,
    aborted: bool,
}

/// Hands the carry from block to block, in block order. Only one carry exists at a time. A
/// worker that finishes inflating block `next` takes the carry and walks its block; one that
/// finishes a later block parks it and goes on. The carry holder then walks the parked blocks
/// that follow its own, in order, hands each one's records back to its worker (as a
/// [`Task::Read`]), and leaves the carry at the first block not yet inflated. No worker waits
/// for the carry.
struct Chain {
    core: Mutex<ChainCore>,
    tasks: Vec<channel::Sender<Task>>,
}

impl Chain {
    fn new(ring: usize, tasks: Vec<channel::Sender<Task>>) -> Self {
        Self {
            core: Mutex::new(ChainCore {
                next: 0,
                carry: Some(Carry {
                    pending: Vec::new(),
                    state: ChainState::Live,
                }),
                parked: (0..ring).map(|_| None).collect(),
                aborted: false,
            }),
            tasks,
        }
    }

    fn abort(&self) {
        self.core.lock().unwrap().aborted = true;
        for t in &self.tasks {
            let _ = t.send(Task::Abort);
        }
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
    /// inflate buffers handed back from parked blocks
    bufs: Vec<Vec<u8>>,
    /// this worker's blocks that are parked or walked but not yet read
    parked: usize,
    /// chain guesses handed back from parked blocks
    specs: Vec<Spec>,
}

impl Worker {
    /// Takes block `idx` (inflated into `data`) through the record chain if the carry into it
    /// has arrived, together with the parked blocks after it, and reads its records; parks it
    /// otherwise. `false` once the pipeline is aborted.
    #[allow(clippy::too_many_arguments)]
    fn on_inflated<O, W>(
        &mut self,
        me: usize,
        idx: u64,
        last: bool,
        inflated: io::Result<usize>,
        data: &mut Vec<u8>,
        chain: &Chain,
        work: &W,
        done_tx: &channel::Sender<Done<O>>,
        free_out_rx: &channel::Receiver<O>,
    ) -> bool
    where
        O: BlockOutput,
        W: Fn(&RawRecord<'_>, &mut O),
    {
        let mut spec: Option<Spec> = None;
        let carry = loop {
            let mut core = chain.core.lock().unwrap();
            if core.aborted {
                return false;
            }
            if core.next == idx {
                break core.carry.take().expect("the carry for this block");
            }
            debug_assert!(core.next < idx);
            let Some(sp) = spec.take() else {
                // guess the block's chain while it is in this cache, then look again
                drop(core);
                let mut sp = self.specs.pop().unwrap_or_default();
                if let Ok(n) = inflated
                    && !last
                {
                    speculate(&data[..n], &mut sp);
                }
                spec = Some(sp);
                continue;
            };
            let ring = core.parked.len() as u64;
            let slot = &mut core.parked[(idx % ring) as usize];
            debug_assert!(slot.is_none());
            let buf = self.bufs.pop().unwrap_or_default();
            *slot = Some(Parked {
                idx,
                owner: me,
                last,
                inflated,
                data: std::mem::replace(data, buf),
                spec: sp,
            });
            self.parked += 1;
            return true;
        };
        if let Some(sp) = spec {
            self.specs.push(sp);
        }
        let (mut carry, plan) = self.walk(idx, last, inflated, data, carry, None);
        // the parked blocks after this one, while they follow on
        let mut j = idx + 1;
        loop {
            let p = {
                let mut core = chain.core.lock().unwrap();
                let ring = core.parked.len() as u64;
                let slot = &mut core.parked[(j % ring) as usize];
                if slot.as_ref().is_some_and(|p| p.idx == j) {
                    slot.take().unwrap()
                } else {
                    core.next = j;
                    core.carry = Some(carry);
                    break;
                }
            };
            let (c, plan) = self.walk(j, p.last, p.inflated, &p.data, carry, Some(&p.spec));
            carry = c;
            let _ = chain.tasks[p.owner].send(Task::Read {
                idx: j,
                last: p.last,
                data: p.data,
                plan,
                spec: p.spec,
            });
            j += 1;
        }
        let out = free_out_rx.try_recv().unwrap_or_default();
        let _ = done_tx.send(self.read(idx, last, plan, data, out, work));
        true
    }

    /// Walks block `idx` (inflated into `data[..n]`) through the record chain with the carry
    /// into it: returns the carry into the next block, and what reading its records takes.
    /// Once the walk reaches an offset of `spec`'s chain, it ends where that chain ends.
    fn walk(
        &mut self,
        idx: u64,
        last: bool,
        inflated: io::Result<usize>,
        data: &[u8],
        mut carry: Carry,
        spec: Option<&Spec>,
    ) -> (Carry, Plan) {
        if carry.state == ChainState::Failed {
            return (carry, Plan::Empty(Ok(())));
        }
        let n = match inflated {
            Ok(n) => n,
            Err(e) => {
                carry.state = ChainState::Failed;
                return (carry, Plan::Empty(Err(e)));
            }
        };
        if last {
            let res = if carry.state == ChainState::Live && !carry.pending.is_empty() {
                Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "failed to fill whole buffer",
                ))
            } else {
                Ok(())
            };
            return (carry, Plan::Empty(res));
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
            let starts = spec.map_or(&[][..], |s| &s.starts[..]);
            let mut k = 0;
            loop {
                while k < starts.len() && (starts[k] as usize) < end {
                    k += 1;
                }
                if k < starts.len() && starts[k] as usize == end {
                    end = *starts.last().unwrap() as usize;
                    if spec.is_some_and(|s| s.stopped) {
                        carry.state = ChainState::Stopped;
                    }
                    break;
                }
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
        (
            carry,
            Plan::Read {
                straddle,
                start,
                end,
            },
        )
    }

    /// Reads the records of walked block `idx` into `out`.
    fn read<O, W>(
        &mut self,
        idx: u64,
        last: bool,
        plan: Plan,
        data: &[u8],
        mut out: O,
        work: &W,
    ) -> Done<O>
    where
        O: BlockOutput,
        W: Fn(&RawRecord<'_>, &mut O),
    {
        out.clear();
        match plan {
            Plan::Empty(res) => Done {
                idx,
                out: res.map(|_| out),
                last,
            },
            Plan::Read {
                straddle,
                start,
                end,
            } => {
                let mut res = Ok(());
                if let Some(s) = straddle {
                    res = read_records(&s, &mut out, work);
                    self.spare = s;
                }
                if res.is_ok() {
                    res = read_records(&data[start..end], &mut out, work);
                }
                Done {
                    idx,
                    out: res.map(|_| out),
                    last: false,
                }
            }
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

    let (task_txs, task_rxs): (Vec<_>, Vec<_>) =
        (0..workers).map(|_| channel::unbounded::<Task>()).unzip();
    let chain = Chain::new(in_flight, task_txs);
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

        for (me, task_rx) in task_rxs.into_iter().enumerate() {
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
                    bufs: Vec::new(),
                    parked: 0,
                    specs: Vec::new(),
                };
                let mut data = Vec::<u8>::new();
                let mut jobs_open = true;
                loop {
                    // reads handed back come first: their blocks hold up the stitcher
                    let task = match task_rx.try_recv() {
                        Ok(t) => Some(t),
                        Err(_) if !jobs_open => {
                            if w.parked == 0 {
                                return;
                            }
                            task_rx.recv().ok()
                        }
                        Err(_) => channel::select! {
                            recv(task_rx) -> t => t.ok(),
                            recv(job_rx) -> job => {
                                let Ok(Job { idx, frame }) = job else {
                                    jobs_open = false;
                                    continue;
                                };
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
                                if !w.on_inflated(
                                    me, idx, last, inflated, &mut data, chain, work, &done_tx,
                                    &free_out_rx,
                                ) {
                                    return;
                                }
                                None
                            }
                        },
                    };
                    match task {
                        Some(Task::Read {
                            idx,
                            last,
                            data: buf,
                            plan,
                            spec,
                        }) => {
                            let out = free_out_rx.try_recv().unwrap_or_default();
                            let _ = done_tx.send(w.read(idx, last, plan, &buf, out, work));
                            w.bufs.push(buf);
                            w.specs.push(spec);
                            w.parked -= 1;
                        }
                        Some(Task::Abort) => return,
                        None => {}
                    }
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
