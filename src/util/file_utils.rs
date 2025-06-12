use anyhow::bail;
use niffler;
use tracing::{error, info};

#[cfg(target_family = "unix")]
use std::os::unix::fs::FileTypeExt;
use std::process::Command;

use std::io::Read; // for chain implementation
use std::path::{Path, PathBuf};

use crate::util::digest_utils;

#[allow(dead_code)]
pub(crate) enum SourceType {
    Fastx(PathBuf),
    ExistingMM2Index(PathBuf),
    ExistingOarfishIndex(PathBuf),
}

impl SourceType {
    /// determine the type of the source and return the path
    /// wrapped in the correct variant.
    pub fn from_path<P: AsRef<Path>>(p: P) -> Self {
        if is_fasta(p.as_ref()).unwrap_or(false) {
            Self::Fastx(PathBuf::from(p.as_ref()))
        } else {
            match digest_utils::read_digest_from_mm2_index(
                p.as_ref().to_str().expect("can be represented as a str"),
            ) {
                // we read a pre-computed digest from an oarfish-constructed
                // minimap2 index
                Ok(_d) => Self::ExistingOarfishIndex(PathBuf::from(p.as_ref())),
                _ => Self::ExistingMM2Index(PathBuf::from(p.as_ref())),
            }
        }
    }

    #[allow(dead_code)]
    pub fn is_raw_mm2_index(&self) -> bool {
        matches!(&self, Self::ExistingMM2Index(_))
    }

    #[allow(dead_code)]
    pub fn is_oarfish_index(&self) -> bool {
        matches!(&self, Self::ExistingOarfishIndex(_))
    }

    #[allow(dead_code)]
    pub fn is_fasta(&self) -> bool {
        matches!(&self, Self::Fastx(_))
    }
}
pub fn is_fasta(fname: &std::path::Path) -> anyhow::Result<bool> {
    match niffler::from_path(fname) {
        Ok((mut reader, _format)) => {
            let mut first_char = vec![0_u8];
            reader.read_exact(&mut first_char)?;
            drop(reader);
            Ok(first_char[0] == b'>' || first_char[0] == b'@')
        }
        _ => Ok(false),
    }
}

#[derive(Debug)]
pub(crate) struct RefSource {
    file_path: std::path::PathBuf,
    concat_handle: Option<std::thread::JoinHandle<anyhow::Result<()>>>,
}

impl RefSource {
    #[allow(dead_code)]
    pub fn file_path_buf(&self) -> PathBuf {
        self.file_path.clone()
    }
    pub fn file_path(&self) -> &Path {
        &self.file_path
    }

    pub fn join_if_needed(self) -> anyhow::Result<()> {
        if let Some(concat_handle) = self.concat_handle {
            if let Err(e) = concat_handle.join() {
                bail!(
                    "Failed to concatenate reference and novel input transcript sequences: {:#?}",
                    e
                );
            } else {
                info!("joined successfully!");
            }
            std::fs::remove_file(self.file_path)?;
        }
        Ok(())
    }
}

pub(crate) fn get_ref_source(
    annotated: Option<PathBuf>,
    novel: Option<PathBuf>,
) -> anyhow::Result<RefSource> {
    let concat_handle: Option<std::thread::JoinHandle<_>>;

    if annotated.as_ref().or(novel.as_ref()).is_none() {
        bail!("at least one of --annotated or --novel but be provided");
    }

    // The `ref_file` input argument is either a FASTA file with annotated
    // sequences, in which case we will compute the proper digest in a separate
    // thread, OR an existing minimap2 index, in which case we won't attempt
    // to treat it as a FASTA file and we will later get the digest from
    // the index.
    let input_path = if annotated.is_some() && novel.is_some() {
        let pid = std::process::id();
        let fifo_fname = format!("combine_transcripts_{}.fifo", pid);
        create_fifo_if_absent(&fifo_fname)?;

        let fifo_fname_clone = fifo_fname.to_string().clone();
        let ref_paths = annotated.clone().expect("annotated should exist");
        let novel_paths = novel.clone().expect("novel txps should exist");
        // a thread that will concatenate the reference transcripts and then the novel
        // trancsripts
        concat_handle = Some(std::thread::spawn(move || {
            let fifo_path = std::path::Path::new(fifo_fname_clone.as_str());
            let mut ff = std::fs::File::options().write(true).open(fifo_path)?;
            let ref_read = std::io::BufReader::new(std::fs::File::open(ref_paths)?);
            let novel_read = std::io::BufReader::new(std::fs::File::open(novel_paths)?);
            let mut reader = ref_read.chain(novel_read);
            info!("before copy");
            match std::io::copy(&mut reader, &mut ff) {
                Err(e) => {
                    error!("Error: {:#?}, in copying input to output", e);
                }
                Ok(nb) => {
                    info!("copied {} bytes from input to output across fifo", nb)
                }
            }
            drop(reader);
            drop(ff);
            Ok(())
        }));

        PathBuf::from(&fifo_fname)
    } else {
        concat_handle = None;
        annotated
            .clone()
            .or(novel.clone())
            .expect("either reference or novel transcripts must be provided")
    };

    Ok(RefSource {
        file_path: input_path.clone(),
        concat_handle,
    })
}

pub fn create_fifo_if_absent<P: AsRef<std::path::Path>>(path_in: P) -> anyhow::Result<()> {
    let path = path_in.as_ref();
    // check if the file already exists and IS a fifo; if so
    // then open it for writing, otherwise make it.
    let fifo_exists = if std::fs::exists(path)? {
        let minfo = std::fs::metadata(path)?;
        if cfg!(target_family = "unix") {
            if minfo.file_type().is_fifo() {
                eprintln!(
                    "The path {} already existed as is a fifo, so using that for communication.",
                    path.display()
                );
                true
            } else {
                // the file existed but wasn't a fifo
                bail!(
                    "The file {} existed already, but wasn't a fifo, so it can't be used as a named pipe. Please remove the file or provide a named pipe instead.",
                    path.display()
                );
            }
        } else {
            // the file existed but wasn't a fifo
            bail!("Named pipes are not supported on non-unix (i.e. non linux/MacOS) systems.");
        }
    } else {
        false
    };

    if !fifo_exists {
        if cfg!(target_family = "unix") {
            let status = Command::new("mkfifo").arg(path).status()?;
            if !status.success() {
                bail!("`mkfifo` command failed with exit status {:#?}", status);
            }
        } else {
            bail!("Named pipes are not supported on non-unix (i.e. non linux/MacOS) systems.");
        }
    }
    Ok(())
}
