use anyhow::{Context, bail};
use seqcol_rs;
use std::io::{Read, Seek, Write};
use std::str;
use tracing::{debug, info, warn};

use crate::NamedDigestVec;

const DIGEST_VERSION: u8 = 4;

pub fn get_digest_from_fasta(
    fpath: &std::path::Path,
) -> std::thread::JoinHandle<anyhow::Result<seqcol_rs::DigestResult>> {
    let fpath_clone = fpath.to_path_buf();
    std::thread::spawn(|| {
        info!("generating reference digest for {}", fpath_clone.display());
        let mut seqcol_obj = seqcol_rs::SeqCol::try_from_fasta_file(fpath_clone).unwrap();
        let digest = seqcol_obj.digest(seqcol_rs::DigestConfig {
            level: seqcol_rs::DigestLevel::Level1,
            additional_attr: vec![seqcol_rs::KnownAttr::SortedNameLengthPairs],
        });
        info!("done");
        digest
    })
}

/// Append the oarfish reference-signature footer
/// (`[json][len:u64][ver:u8]["OARFISHSIG"]`) to an on-disk index file. This is
/// format-neutral: it works for a legacy `.mmi` and for a rammap `RMMI`
/// file, because both index loaders stop at the end of their own structure and
/// ignore trailing bytes (the `.mmi` loader has no EOF check; rammap's
/// `bincode::deserialize_from` stops at the struct boundary). The footer is read
/// back by [`read_digest_footer`], which seeks from the end of the file.
pub(crate) fn append_digest_footer(
    idx_file: &str,
    digest: &NamedDigestVec,
) -> anyhow::Result<()> {
    if std::fs::exists(idx_file)? {
        let mut writer = std::fs::OpenOptions::new()
            .read(true)
            .append(true)
            .open(idx_file)?;
        let digest_json = digest.to_json();
        let json_str = serde_json::to_string(&digest_json)?;

        let str_len = json_str.len();
        let len_bytes = str_len.to_le_bytes();

        let version = DIGEST_VERSION;
        let ver_bytes = version.to_le_bytes();

        const OARFISH_FOOTER_MAGIC: &str = "OARFISHSIG";
        writer.write_all(json_str.as_bytes())?;
        writer.write_all(&len_bytes)?;
        writer.write_all(&ver_bytes)?;
        writer.write_all(OARFISH_FOOTER_MAGIC.as_bytes())?;
        drop(writer);
        Ok(())
    } else {
        anyhow::bail!("The file {} did not exist as an index", idx_file);
    }
}

fn get_digest_from_value(value: &mut serde_json::Value) -> anyhow::Result<seqcol_rs::DigestResult> {
    let seq_col_digests_value = value
        .get("seqcol_digest")
        .expect("a seqcol_digest entry must exist!");
    let seq_col_digests = seq_col_digests_value
        .as_object()
        .expect("seqcol_digest should be an object");

    // ensure we have the appropriate values
    let _lengths = seq_col_digests
        .get("lengths")
        .expect("lengths field should exist")
        .is_string();
    let _names = seq_col_digests
        .get("names")
        .expect("names field should exist")
        .is_string();
    let _seqs = seq_col_digests
        .get("sequences")
        .expect("sequences field should exist")
        .is_string();
    let _sorted_name_length_pairs = seq_col_digests
        .get("sorted_name_length_pairs")
        .expect("sorted_name_length_pairs field should exist")
        .is_string();

    let sha_digests = value["sha256_digests"]
        .as_object()
        .expect("should be an object");
    let sha256_names = sha_digests["sha256_names"]
        .as_str()
        .expect("should be a string");
    let sha256_seqs = sha_digests["sha256_seqs"]
        .as_str()
        .expect("should be a string");
    Ok(seqcol_rs::DigestResult {
        sq_digest: seqcol_rs::DigestLevelResult::Level1(seqcol_rs::Level1Digest {
            digests: seq_col_digests_value.clone(),
        }),
        sha256_seqs: Some(sha256_seqs.to_owned()),
        sha256_names: Some(sha256_names.to_owned()),
    })
}

/// Reads the oarfish reference-signature footer from an on-disk index file
/// (legacy `.mmi` or rammap `RMMI`), *if* that index was originally created
/// with oarfish. Otherwise, it returns an error. The footer is located by
/// seeking from the end of the file, so it is independent of the index body
/// format.
pub(crate) fn read_digest_footer(idx_file: &str) -> anyhow::Result<NamedDigestVec> {
    if std::fs::exists(idx_file)? {
        let mut file = std::fs::OpenOptions::new().read(true).open(idx_file)?;

        const OARFISH_FOOTER_MAGIC: &str = "OARFISHSIG";
        let magic_len = OARFISH_FOOTER_MAGIC.len() as i64;
        file.seek(std::io::SeekFrom::End(-magic_len))?;

        let mut buf = vec![0u8; magic_len as usize];
        file.read_exact(&mut buf)?;

        if str::from_utf8(&buf)? == OARFISH_FOOTER_MAGIC {
            info!("succesfully detected oarfish magic footer in index");
            let ver_offset = magic_len + 1;
            file.seek(std::io::SeekFrom::End(-ver_offset))?;
            let mut ver_buf: [u8; 1] = [0];
            file.read_exact(&mut ver_buf)?;
            let stored_ver = u8::from_le_bytes(ver_buf);
            if stored_ver < DIGEST_VERSION {
                warn!(
                    "the provided index has oarfish signature version {}, but the current version is {}.",
                    u32::from(stored_ver),
                    u32::from(DIGEST_VERSION)
                );
                warn!("please consider re-creating the index.");
            }

            let str_size_offset = ver_offset + 8;
            file.seek(std::io::SeekFrom::End(-str_size_offset))?;
            let mut sig_len_buf: [u8; 8] = [0, 0, 0, 0, 0, 0, 0, 0];
            file.read_exact(&mut sig_len_buf)?;
            let sig_len = usize::from_le_bytes(sig_len_buf);

            let sig_offset = magic_len + 1 + 8 + (sig_len as i64);
            buf.resize(sig_len, 0);
            file.seek(std::io::SeekFrom::End(-sig_offset))?;
            file.read_exact(&mut buf)?;
            let mut value: serde_json::Value = serde_json::from_slice(&buf)?;
            debug!("{}", value);

            let mut ndv = NamedDigestVec::new();

            if let Some(_seq_col_digest_value) = value.get("seqcol_digest") {
                warn!(
                    "Found an unannotated seqcol digest entry from an old version of the index! Recording it as \"deprecated\""
                );
                let digest = get_digest_from_value(&mut value)?;
                ndv.push(("deprecated".to_string(), digest));
            } else if let serde_json::Value::Object(omap) = &mut value {
                // we have to iterate
                for (k, v) in omap.iter_mut() {
                    let digest = get_digest_from_value(v)?;
                    ndv.push((k.to_string(), digest));
                }
            } else {
                bail!("The top-level value should be an object!");
            }
            Ok(ndv)
        } else {
            bail!("index did not have an oarfish footer!");
        }
    } else {
        bail!("could not find file {}", idx_file);
    }
}

/// A reference signature that is either already computed or still being
/// computed on a scoped thread (see [`DeferredDigest::spawn_from_header`]).
/// The signature is only needed when `meta_info.json` is assembled, so the BAM
/// path computes it while the alignments are parsed and quantified, and
/// [`DeferredDigest::resolve`] joins it just before it is written.
pub(crate) enum DeferredDigest<'scope> {
    Ready(NamedDigestVec),
    Pending {
        name: String,
        job: std::thread::ScopedJoinHandle<'scope, anyhow::Result<seqcol_rs::DigestResult>>,
    },
}

impl From<NamedDigestVec> for DeferredDigest<'_> {
    fn from(d: NamedDigestVec) -> Self {
        DeferredDigest::Ready(d)
    }
}

impl<'scope> DeferredDigest<'scope> {
    /// Start [`digest_from_header`] on a thread of `scope`; the result is
    /// recorded under `name` once resolved.
    pub(crate) fn spawn_from_header<'env>(
        scope: &'scope std::thread::Scope<'scope, 'env>,
        header: &'env noodles_sam::header::Header,
        name: &str,
    ) -> Self {
        DeferredDigest::Pending {
            name: name.to_string(),
            job: scope.spawn(move || digest_from_header(header)),
        }
    }

    /// The computed signature, waiting for the digest thread if it is still running.
    pub(crate) fn resolve(self) -> anyhow::Result<NamedDigestVec> {
        match self {
            DeferredDigest::Ready(d) => Ok(d),
            DeferredDigest::Pending { name, job } => {
                let d = job
                    .join()
                    .map_err(|_| anyhow::anyhow!("the seqcol digest thread panicked"))??;
                Ok(vec![(name, d)].into())
            }
        }
    }
}

pub(crate) fn digest_from_header(
    header: &noodles_sam::header::Header,
) -> anyhow::Result<seqcol_rs::DigestResult> {
    info!("calculating seqcol digest");
    let d = level1_digest_from_names_lengths(
        header
            .reference_sequences()
            .iter()
            .map(|(k, v)| (k.as_slice(), usize::from(v.length()))),
    )
    .context("failed to compute the seqcol digest for the information from the alignment header")?;
    info!("done calculating seqcol digest");
    Ok(d)
}

/// Feeds what is written into a hash.
struct HashWriter<'a, D>(&'a mut D);

impl<D: sha2::Digest> Write for HashWriter<'_, D> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.update(buf);
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// seqcol's `sha512t24u` of everything hashed into `h`: the first 24 bytes of the SHA-512,
/// base64url-encoded (32 characters).
fn sha512t24u(h: sha2::Sha512) -> [u8; 32] {
    use base64::Engine as _;
    use sha2::Digest as _;
    let d = h.finalize();
    let mut out = [0u8; 32];
    let n = base64::prelude::BASE64_URL_SAFE
        .encode_slice(&d[..24], &mut out)
        .expect("24 bytes encode to 32 base64 characters");
    debug_assert_eq!(n, 32);
    out
}

/// The Level1 seqcol digest of a SAM/BAM header's reference names and lengths with the
/// `sorted_name_length_pairs` attribute: the same [`seqcol_rs::DigestResult`] as
/// `SeqCol::from_sam_header(..).digest(Level1, [SortedNameLengthPairs])`, computed in one
/// streaming pass. The canonical JSON of the `names` and `lengths` arrays is hashed as it is
/// written instead of being built as `serde_json` values and strings, and each name/length
/// pair keeps only its 32-byte digest, so the memory is 32 bytes per sequence (seqcol builds
/// several copies of every name, about 0.35 KB per sequence).
fn level1_digest_from_names_lengths<'a, I>(it: I) -> anyhow::Result<seqcol_rs::DigestResult>
where
    I: IntoIterator<Item = (&'a [u8], usize)>,
{
    use sha2::Digest as _;
    let it = it.into_iter();
    let mut sha256_names = sha2::Sha256::new();
    let mut names = sha2::Sha512::new();
    let mut lengths = sha2::Sha512::new();
    let mut pairs: Vec<[u8; 32]> = Vec::with_capacity(it.size_hint().0);
    names.update(b"[");
    lengths.update(b"[");
    for (i, (n, l)) in it.enumerate() {
        sha256_names.update(n);
        // seqcol requires UTF-8 names (`from_sam_header` unwraps the conversion).
        let name = std::str::from_utf8(n).unwrap();
        if i > 0 {
            names.update(b",");
            lengths.update(b",");
        }
        serde_json::to_writer(HashWriter(&mut names), name)?;
        write!(HashWriter(&mut lengths), "{}", l as u64)?;
        // The canonical JSON of {"name": n, "length": l} (keys sorted).
        let mut pair = sha2::Sha512::new();
        write!(HashWriter(&mut pair), "{{\"length\":{},\"name\":", l as u64)?;
        serde_json::to_writer(HashWriter(&mut pair), name)?;
        pair.update(b"}");
        pairs.push(sha512t24u(pair));
    }
    names.update(b"]");
    lengths.update(b"]");
    // The pair digests, sorted as strings, as a JSON array of strings (base64url characters
    // need no escaping).
    pairs.sort_unstable();
    let mut snlp = sha2::Sha512::new();
    snlp.update(b"[");
    for (i, d) in pairs.iter().enumerate() {
        if i > 0 {
            snlp.update(b",");
        }
        snlp.update(b"\"");
        snlp.update(d);
        snlp.update(b"\"");
    }
    snlp.update(b"]");
    drop(pairs);

    let as_value = |d: [u8; 32]| {
        serde_json::Value::String(
            String::from_utf8(d.to_vec()).expect("base64url digests are ASCII"),
        )
    };
    let mut res = serde_json::json!({});
    res["names"] = as_value(sha512t24u(names));
    res["lengths"] = as_value(sha512t24u(lengths));
    res["sorted_name_length_pairs"] = as_value(sha512t24u(snlp));
    Ok(seqcol_rs::DigestResult {
        sq_digest: seqcol_rs::DigestLevelResult::Level1(seqcol_rs::Level1Digest { digests: res }),
        sha256_names: Some(hex::encode(sha256_names.finalize())),
        sha256_seqs: None,
    })
}

/// Computes a full Level1 SeqCol `DigestResult` from a loaded rammap aligner by
/// reconstructing each reference sequence from the (4-bit packed) index. It
/// yields the same sequence-content digest as [`get_digest_from_fasta`] would on
/// the original FASTA, and is used as the fallback signature when a rammap index
/// has no oarfish footer. The index must still retain its sequences
/// (`index.has_sequences()`); a sequence-stripped index cannot produce a
/// sequence digest and the caller should fall back to [`digest_from_header`].
pub(crate) fn digest_from_rammap_index(
    aligner: &rammap::api::Aligner,
) -> anyhow::Result<seqcol_rs::DigestResult> {
    use anyhow::ensure;
    // nt4 encoding: 0=A, 1=C, 2=G, 3=T, 4=N (rammap `Index::NT4_TO_ASCII`).
    const NT4_TO_ASCII: [u8; 5] = *b"ACGTN";
    let idx = aligner.index();
    ensure!(
        idx.has_sequences(),
        "rammap index does not retain reference sequences; cannot compute a sequence digest"
    );
    info!("calculating seqcol digest from the rammap index");
    // Lazily reconstruct (name, ASCII-sequence) pairs so we never hold all of
    // the reference sequences in memory at once.
    let name_seq_iter = (0..idx.seqs.len()).map(|rid| {
        let ts = &idx.seqs[rid];
        let nt4 = idx.get_region_nt4(rid, 0, ts.len);
        let seq: String = nt4
            .iter()
            .map(|&b| NT4_TO_ASCII[b as usize] as char)
            .collect();
        (ts.name.clone(), seq)
    });
    let mut sq = seqcol_rs::SeqCol::try_from_name_seq_iter(name_seq_iter)?;
    let d = sq
        .digest(seqcol_rs::DigestConfig {
            level: seqcol_rs::DigestLevel::Level1,
            additional_attr: vec![seqcol_rs::KnownAttr::SortedNameLengthPairs],
        })
        .context("failed to compute the seqcol digest from the rammap index")?;
    info!("done calculating seqcol digest");
    Ok(d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Round-trip: build a rammap index, write it with an appended oarfish
    /// footer, reload it, and confirm (a) the loader ignores the trailing footer
    /// (sequences load correctly), (b) the footer is readable from EOF, and
    /// (c) the index-derived digest matches the FASTA-derived digest. This guards
    /// against a future rammap-core change to trailing-byte handling.
    #[test]
    fn rammap_index_footer_roundtrip() {
        let dir = std::env::temp_dir();
        let stamp = std::process::id();
        let fasta_path = dir.join(format!("oarfish_sig_test_{stamp}.fa"));
        let idx_path = dir.join(format!("oarfish_sig_test_{stamp}.rmmi"));

        let seqs: [(&str, &str); 3] = [
            ("seqA", "ACGTACGTACGTACGTACGTACGTACGTACGT"),
            ("seqB", "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA"),
            ("seqC", "ACGTNNNNACGTACGTACGTACGTNNNNACGT"),
        ];
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            for (n, s) in &seqs {
                writeln!(f, ">{n}\n{s}").unwrap();
            }
        }

        // FASTA-derived reference signature (the source of truth).
        let fasta_digest = get_digest_from_fasta(&fasta_path).join().unwrap().unwrap();
        let fasta_digest_json = serde_json::to_string(&fasta_digest.to_json()).unwrap();
        let mut ndv = NamedDigestVec::new();
        ndv.push(("annotated_transcripts_digest".to_string(), fasta_digest));

        // Build + persist the rammap index, then append the footer.
        let aligner = rammap::api::Aligner::from_fasta(
            fasta_path.to_str().unwrap(),
            rammap::api::Preset::MapOnt,
        )
        .unwrap();
        let idx_str = idx_path.to_str().unwrap();
        aligner.save_index(idx_str).unwrap();
        append_digest_footer(idx_str, &ndv).unwrap();

        // (a) Reload through rammap's own loader: trailing footer is ignored and
        // all sequences are present.
        let reloaded = rammap::api::Aligner::from_index(idx_str, rammap::api::Preset::MapOnt).unwrap();
        assert_eq!(reloaded.index().seqs.len(), seqs.len());
        assert!(reloaded.index().has_sequences());

        // (b) Footer round-trips.
        let read_back = read_digest_footer(idx_str).unwrap();
        let read_json = serde_json::to_string(&read_back.to_json()).unwrap();
        let orig_json = serde_json::to_string(&ndv.to_json()).unwrap();
        assert_eq!(read_json, orig_json, "footer digest changed across round-trip");

        // (c) Index-derived digest equals the FASTA-derived digest.
        let idx_digest = digest_from_rammap_index(&reloaded).unwrap();
        assert_eq!(
            serde_json::to_string(&idx_digest.to_json()).unwrap(),
            fasta_digest_json,
            "index-derived sequence digest disagrees with FASTA-derived digest"
        );

        let _ = std::fs::remove_file(&fasta_path);
        let _ = std::fs::remove_file(&idx_path);
    }

    /// Round-trip for a spearfish index written as `--index-out` writes it: built to a file,
    /// with the oarfish footer appended. Confirms (a) spearfish opens it despite the trailing
    /// footer and it holds the same transcripts, (b) it maps a read as the in-memory index
    /// does, and (c) the footer reads back. This guards against a future spearfish change to
    /// trailing-byte handling.
    #[test]
    fn spearfish_index_footer_roundtrip() {
        let dir = std::env::temp_dir();
        let stamp = std::process::id();
        let fasta_path = dir.join(format!("oarfish_sig_test_{stamp}_sf.fa"));
        let idx_path = dir.join(format!("oarfish_sig_test_{stamp}.sfi"));

        // three pseudo-random 3 kb transcripts (xorshift, fixed seed)
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut base = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state >> 33) as usize % 4]
        };
        let txs: Vec<(String, String, Vec<u8>)> = (0..3)
            .map(|i| {
                (
                    format!("tx{i}"),
                    format!("tx{i}"),
                    (0..3000).map(|_| base()).collect(),
                )
            })
            .collect();
        let read = txs[1].2[500..2000].to_vec();
        {
            let mut f = std::fs::File::create(&fasta_path).unwrap();
            for (n, _, s) in &txs {
                writeln!(f, ">{n}\n{}", std::str::from_utf8(s).unwrap()).unwrap();
            }
        }

        let opts = spearfish::IndexOptions::default();
        let (mem, _) = spearfish::Index::build(txs.clone(), &opts).unwrap();
        spearfish::Index::build_to_file(txs, &opts, &idx_path).unwrap();
        let idx_str = idx_path.to_str().unwrap();
        let mut ndv = NamedDigestVec::new();
        ndv.push((
            "annotated_transcripts_digest".to_string(),
            get_digest_from_fasta(&fasta_path).join().unwrap().unwrap(),
        ));
        append_digest_footer(idx_str, &ndv).unwrap();

        // (a) spearfish opens the file with its footer; same transcripts.
        assert!(crate::util::aligner::is_spearfish_index(&idx_path));
        let disk = spearfish::Index::open(&idx_path).unwrap();
        assert_eq!(disk.n_transcripts(), mem.n_transcripts());
        for t in 0..mem.n_transcripts() as u32 {
            assert_eq!(disk.tx_name(t), mem.tx_name(t));
            assert_eq!(disk.tx_len(t), mem.tx_len(t));
        }

        // (b) the read maps to the same hits through either index.
        let map_opts = spearfish::MapOptions::preset(spearfish::Preset::Ont);
        let hits = |idx| {
            spearfish::Mapper::new(std::sync::Arc::new(idx), &map_opts)
                .unwrap()
                .map(&read)
                .hits
        };
        let (mem_hits, disk_hits) = (hits(mem), hits(disk));
        assert!(
            mem_hits.iter().any(|h| h.tx == 1),
            "read did not map to its transcript"
        );
        assert_eq!(disk_hits, mem_hits);

        // (c) the footer round-trips.
        let read_back = read_digest_footer(idx_str).unwrap();
        assert_eq!(
            serde_json::to_string(&read_back.to_json()).unwrap(),
            serde_json::to_string(&ndv.to_json()).unwrap()
        );

        let _ = std::fs::remove_file(&fasta_path);
        let _ = std::fs::remove_file(&idx_path);
        let _ = std::fs::remove_file(spearfish::provenance::sidecar(&idx_path));
    }
}
