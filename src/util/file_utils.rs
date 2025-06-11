use anyhow::bail;
#[cfg(target_family = "unix")]
use std::os::unix::fs::FileTypeExt;
use std::process::Command;

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
