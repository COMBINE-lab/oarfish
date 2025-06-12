use anyhow::{Context, bail};
use minimap2_sys::MmIdx;
use seqcol_rs;
use std::io::{Read, Seek, Write};
use std::str;
use std::sync::Arc;
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

pub(crate) fn append_digest_to_mm2_index(
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
        anyhow::bail!("The file {} did not exist as a minimap2 index", idx_file);
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

/// Computes a `DigestResult` from a provided minimap2 index
/// file, *if* that index was originally created with oarfish.
/// Otherwise, it return an error.
pub(crate) fn read_digest_from_mm2_index(idx_file: &str) -> anyhow::Result<NamedDigestVec> {
    if std::fs::exists(idx_file)? {
        let mut file = std::fs::OpenOptions::new().read(true).open(idx_file)?;

        const OARFISH_FOOTER_MAGIC: &str = "OARFISHSIG";
        let magic_len = OARFISH_FOOTER_MAGIC.len() as i64;
        file.seek(std::io::SeekFrom::End(-magic_len))?;

        let mut buf = vec![0u8; magic_len as usize];
        file.read_exact(&mut buf)?;

        if str::from_utf8(&buf)? == OARFISH_FOOTER_MAGIC {
            info!("succesfully detected oarfish magic footer in minimap2 index");
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
            bail!("minimap2 index did not have an oarfish footer!");
        }
    } else {
        bail!("could not find file {}", idx_file);
    }
}

pub(crate) fn digest_from_header(
    header: &noodles_sam::header::Header,
) -> anyhow::Result<seqcol_rs::DigestResult> {
    let seqcol_digest = {
        info!("calculating seqcol digest");
        let mut sc = seqcol_rs::SeqCol::from_sam_header(
            header
                .reference_sequences()
                .iter()
                .map(|(k, v)| (k.as_slice(), v.length().into())),
        );
        let d = sc
            .digest(seqcol_rs::DigestConfig {
                level: seqcol_rs::DigestLevel::Level1,
                additional_attr: vec![seqcol_rs::KnownAttr::SortedNameLengthPairs],
            })
            .context(
                "failed to compute the seqcol digest for the information from the alignment header",
            )?;
        info!("done calculating seqcol digest");
        d
    };
    Ok(seqcol_digest)
}

pub(crate) fn digest_from_index(mmi: &Arc<MmIdx>) -> anyhow::Result<seqcol_rs::DigestResult> {
    let idx_iter = crate::util::mm_utils::MMIdxNameSeqIter::from_idx(mmi);
    let mut sq = seqcol_rs::SeqCol::try_from_name_seq_iter(idx_iter)?;
    let d = sq
        .digest(seqcol_rs::DigestConfig {
            level: seqcol_rs::DigestLevel::Level1,
            additional_attr: vec![seqcol_rs::KnownAttr::SortedNameLengthPairs],
        })
        .context("failed to compute the seqcol digest for the information from the index")?;
    info!("done calculating seqcol digest");
    Ok(d)
}
