use arrow2::{
    array::Array,
    chunk::Chunk,
    datatypes::Schema,
    io::parquet::write::{
        CompressionOptions, Encoding, FileWriter, RowGroupIterator, Version, WriteOptions,
        transverse,
    },
};
use std::fs::File;

/// Write a chunk of values `chunk` to a file specified at the
/// provided `path` using `scheme`. This raises an [anyhow::Error] if
/// there is an error in creating the file or writing the contents.
pub(crate) fn write_chunk_to_file(
    path: &str,
    schema: Schema,
    chunk: Chunk<Box<dyn Array>>,
) -> anyhow::Result<()> {
    let options = WriteOptions {
        write_statistics: true,
        compression: CompressionOptions::Zstd(None),
        version: Version::V2,
        data_pagesize_limit: None,
    };

    // Create a new empty file
    let file = File::create(path)?;

    let mut writer = FileWriter::try_new(file, schema.clone(), options)?;

    let iter = vec![Ok(chunk)];
    let encodings = schema
        .fields
        .iter()
        .map(|f| transverse(&f.data_type, |_| Encoding::Plain))
        .collect::<Vec<Vec<Encoding>>>();
    let row_groups = RowGroupIterator::try_new(iter.into_iter(), &schema, options, encodings)?;
    for group in row_groups {
        writer.write(group?)?;
    }
    let _size = writer.end(None)?;
    Ok(())
}
