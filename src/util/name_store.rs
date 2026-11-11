//! The read names kept for `--write-assignment-probs`, in the order their reads were added.
//!
//! Names stay in memory until they reach [`SPILL_AFTER_BYTES`]; after that, every name (those
//! held so far and all later ones) goes to an anonymous, LZ4-compressed temporary file that is
//! read back once, in order, when the `.prob` file is written.
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

/// How many bytes of names to hold in memory before moving them to a temporary file.
const SPILL_AFTER_BYTES: usize = 256 * 1024 * 1024;

/// An append-only list of read names that moves to disk once it grows large.
pub struct NameStore {
    mem: Vec<String>,
    mem_bytes: usize,
    spill: Option<lz4::Encoder<BufWriter<File>>>,
    spill_after: usize,
}

impl Default for NameStore {
    fn default() -> Self {
        Self::new()
    }
}

impl NameStore {
    pub fn new() -> Self {
        Self::with_spill_after(SPILL_AFTER_BYTES)
    }

    fn with_spill_after(spill_after: usize) -> Self {
        Self {
            mem: Vec::new(),
            mem_bytes: 0,
            spill: None,
            spill_after,
        }
    }

    /// Appends `name`.
    pub fn push(&mut self, name: String) -> io::Result<()> {
        if let Some(enc) = self.spill.as_mut() {
            return write_name(enc, &name);
        }
        self.mem_bytes += name.len();
        self.mem.push(name);
        if self.mem_bytes > self.spill_after {
            let mut enc = lz4::EncoderBuilder::new()
                .level(1)
                .build(BufWriter::with_capacity(1 << 20, tempfile::tempfile()?))?;
            for n in self.mem.drain(..) {
                write_name(&mut enc, &n)?;
            }
            self.mem = Vec::new();
            self.spill = Some(enc);
        }
        Ok(())
    }

    /// The names, in the order they were added.
    pub fn into_iter(self) -> Box<dyn Iterator<Item = io::Result<String>>> {
        match self.spill {
            None => Box::new(self.mem.into_iter().map(Ok)),
            Some(enc) => match reopen(enc) {
                Ok(r) => Box::new(SpillIter { r: Some(r) }),
                Err(e) => Box::new(std::iter::once(Err(e))),
            },
        }
    }
}

/// Writes `name` as its length (u32, little-endian) followed by its bytes.
fn write_name<W: Write>(w: &mut W, name: &str) -> io::Result<()> {
    let len = u32::try_from(name.len())
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "read name too long"))?;
    w.write_all(&len.to_le_bytes())?;
    w.write_all(name.as_bytes())
}

/// Finishes the LZ4 stream and returns a reader positioned at its start.
fn reopen(enc: lz4::Encoder<BufWriter<File>>) -> io::Result<lz4::Decoder<BufReader<File>>> {
    let (w, res) = enc.finish();
    res?;
    let mut f = w.into_inner().map_err(|e| e.into_error())?;
    f.seek(SeekFrom::Start(0))?;
    lz4::Decoder::new(BufReader::with_capacity(1 << 20, f))
}

struct SpillIter {
    r: Option<lz4::Decoder<BufReader<File>>>,
}

impl Iterator for SpillIter {
    type Item = io::Result<String>;

    fn next(&mut self) -> Option<Self::Item> {
        let r = self.r.as_mut()?;
        let mut len = [0u8; 4];
        // A clean end of stream is the end of the names.
        match r.read(&mut len[..1]) {
            Ok(0) => {
                self.r = None;
                return None;
            }
            Ok(_) => {}
            Err(e) => {
                self.r = None;
                return Some(Err(e));
            }
        }
        let res = r.read_exact(&mut len[1..]).and_then(|_| {
            let mut buf = vec![0u8; u32::from_le_bytes(len) as usize];
            r.read_exact(&mut buf)?;
            String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        });
        if res.is_err() {
            self.r = None;
        }
        Some(res)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(n: usize) -> Vec<String> {
        (0..n)
            .map(|i| match i % 3 {
                0 => String::new(),
                1 => format!("read_{i}\0"),
                _ => format!("m{i}/ccs/{}", "x".repeat(i % 50)),
            })
            .collect()
    }

    fn round_trip(spill_after: usize, n: usize) {
        let want = names(n);
        let mut s = NameStore::with_spill_after(spill_after);
        for w in &want {
            s.push(w.clone()).unwrap();
        }
        let bytes: usize = want.iter().map(String::len).sum();
        assert_eq!(s.spill.is_some(), bytes > spill_after);
        let got: Vec<String> = s.into_iter().collect::<io::Result<_>>().unwrap();
        assert_eq!(got, want);
    }

    #[test]
    fn in_memory() {
        round_trip(usize::MAX, 1000);
        round_trip(usize::MAX, 0);
    }

    #[test]
    fn spilled_from_the_start() {
        round_trip(0, 2);
        round_trip(0, 200_000);
    }

    #[test]
    fn spilled_part_way() {
        round_trip(10_000, 200_000);
    }
}
