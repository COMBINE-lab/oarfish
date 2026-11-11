use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{Field, Schema};
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, Encoding, ZstdLevel};
use parquet::file::properties::{EnabledStatistics, WriterProperties, WriterVersion};
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

/// Write `columns` (name, values) as one table to the Parquet file `path`:
/// non-nullable columns, plain-encoded and zstd-compressed, with statistics.
/// This raises an [anyhow::Error] if there is an error in creating the file or
/// writing the contents.
pub(crate) fn write_columns(path: &Path, columns: Vec<(String, ArrayRef)>) -> anyhow::Result<()> {
    let fields: Vec<Field> = columns
        .iter()
        .map(|(name, a)| Field::new(name, a.data_type().clone(), false))
        .collect();
    let schema = Arc::new(Schema::new(fields));
    let batch = RecordBatch::try_new(
        schema.clone(),
        columns.into_iter().map(|(_, a)| a).collect(),
    )?;

    let props = WriterProperties::builder()
        .set_writer_version(WriterVersion::PARQUET_2_0)
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .set_dictionary_enabled(false)
        .set_encoding(Encoding::PLAIN)
        .set_statistics_enabled(EnabledStatistics::Page)
        .build();

    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow_array::{Array, Float64Array};
    use arrow_schema::DataType;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    #[test]
    fn columns_read_back_unchanged() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.infreps.pq");
        let a: Vec<f64> = (0..100_000).map(|i| i as f64 * 0.25).collect();
        let b: Vec<f64> = (0..100_000).map(|i| (i % 7) as f64).collect();
        write_columns(
            &path,
            vec![
                (
                    "bootstrap.0".into(),
                    Arc::new(Float64Array::from(a.clone())),
                ),
                (
                    "bootstrap.1".into(),
                    Arc::new(Float64Array::from(b.clone())),
                ),
            ],
        )
        .unwrap();

        let reader = ParquetRecordBatchReaderBuilder::try_new(File::open(&path).unwrap())
            .unwrap()
            .build()
            .unwrap();
        let (mut got_a, mut got_b) = (vec![], vec![]);
        for batch in reader {
            let batch = batch.unwrap();
            let s = batch.schema();
            assert_eq!(s.field(0).name(), "bootstrap.0");
            assert_eq!(s.field(1).name(), "bootstrap.1");
            assert_eq!(s.field(0).data_type(), &DataType::Float64);
            assert!(s.fields().iter().all(|f| !f.is_nullable()));
            let col = |i: usize| {
                batch
                    .column(i)
                    .as_any()
                    .downcast_ref::<Float64Array>()
                    .unwrap()
                    .values()
                    .to_vec()
            };
            got_a.extend(col(0));
            got_b.extend(col(1));
        }
        assert_eq!((got_a, got_b), (a, b));
    }
}
