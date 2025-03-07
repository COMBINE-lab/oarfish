use anyhow::{bail, Context};
use seqcol_rs;
use std::io::{Read, Seek, Write};
use std::str;
use tracing::{debug, info};

pub(crate) fn append_digest_to_mm2_index(
    idx_file: &str,
    digest: &seqcol_rs::DigestResult,
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

        let version = 1u8;
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

pub(crate) fn read_digest_from_mm2_index(
    idx_file: &str,
) -> anyhow::Result<seqcol_rs::DigestResult> {
    if std::fs::exists(idx_file)? {
        let mut file = std::fs::OpenOptions::new().read(true).open(idx_file)?;

        const OARFISH_FOOTER_MAGIC: &str = "OARFISHSIG";
        let magic_len = OARFISH_FOOTER_MAGIC.len() as i64;
        file.seek(std::io::SeekFrom::End(-magic_len))?;

        let mut buf = vec![0u8; magic_len as usize];
        file.read_exact(&mut buf)?;

        if str::from_utf8(&buf)? == OARFISH_FOOTER_MAGIC {
            info!("succesfully detected oarfish magic footer in minimap2 index");
            let str_size_offset = magic_len + 1 + 8;
            file.seek(std::io::SeekFrom::End(-str_size_offset))?;
            let mut sig_len_buf: [u8; 8] = [0, 0, 0, 0, 0, 0, 0, 0];
            file.read_exact(&mut sig_len_buf)?;
            let sig_len = usize::from_le_bytes(sig_len_buf);

            let sig_offset = magic_len + 1 + 8 + (sig_len as i64);
            buf.resize(sig_len, 0);
            file.seek(std::io::SeekFrom::End(-sig_offset))?;
            file.read_exact(&mut buf)?;
            let value: serde_json::Value = serde_json::from_slice(&buf)?;

            debug!("{}", value);
            let seq_col_digests = value["seqcol_digest"]
                .as_object()
                .expect("should be an object");
            let lengths = seq_col_digests["lengths"]
                .as_str()
                .expect("should be a string");
            let names = seq_col_digests["names"]
                .as_str()
                .expect("should be a string");
            let seqs = seq_col_digests["sequences"]
                .as_str()
                .expect("should have sequences");
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
                    lengths: lengths.to_owned(),
                    names: names.to_owned(),
                    sequences: Some(seqs.to_owned()),
                    sorted_name_length_pairs: None,
                }),
                sha256_seqs: Some(sha256_seqs.to_owned()),
                sha256_names: Some(sha256_names.to_owned()),
            })
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
        let sc = seqcol_rs::SeqCol::from_sam_header(
            header
                .reference_sequences()
                .iter()
                .map(|(k, v)| (k.as_slice(), v.length().into())),
        );
        let d = sc
            .digest(seqcol_rs::DigestConfig {
                level: seqcol_rs::DigestLevel::Level1,
                with_seqname_pairs: false,
            })
            .context(
                "failed to compute the seqcol digest for the information from the alignment header",
            )?;
        info!("done calculating seqcol digest");
        d
    };
    Ok(seqcol_digest)
}
