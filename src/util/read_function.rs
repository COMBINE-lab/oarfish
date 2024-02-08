use crate::util::oarfish_types::ShortReadRecord;
use csv::ReaderBuilder;
use ndarray::Array2;
use std::fs::File;


/// Read the short read quantification from the file `short_read_path`
pub fn read_short_quant_vec(short_read_path: &str, txps_name: &[String]) -> anyhow::Result<Vec<f64>> {
    // try to open the short read file
    let file = File::open(short_read_path)?;

    // read the data as a tab-delimited TSV file.
    let mut rdr = ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(file);

    // deserialize CSV records into a Vec of ShortReadRecords
    let records: Vec<ShortReadRecord> = rdr
        .deserialize()
        .collect::<Result<Vec<ShortReadRecord>, csv::Error>>()
        .unwrap_or_else(|err| {
            eprintln!("Failed to deserialize CSV records: {}", err);
            std::process::exit(1);
        })
        .into_iter()
        .collect();

    // Convert the Vec of Records into an ndarray::Array2
    let mut data_array = Array2::from_shape_vec(
        (records.len(), 2),
        records
            .iter()
            .flat_map(|r| vec![r.name.clone(), r.num_reads.to_string()])
            .collect(),
    )
    .unwrap();

    //obtain the first column of the data_array and check if all the elements exist in txps_name
    let first_column: Vec<String> = data_array.column(0).to_owned().to_vec();
    let all_elements_exist = first_column
        .iter()
        .all(|element| txps_name.contains(element));

    if all_elements_exist && (first_column.len() != txps_name.len()) {
        let txps_not_in_first_column: Vec<String> = txps_name
            .iter()
            .filter(|&x| !first_column.contains(x))
            .cloned()
            .collect();

        // Create new records for the missing elements and append them to the records vector
        let missing_records: Vec<ShortReadRecord> = txps_not_in_first_column
            .iter()
            .map(|name| ShortReadRecord::empty(&name))
            .collect();

        let mut updated_records = records.clone();
        updated_records.extend(missing_records);

        // Update the data_array with the new records
        data_array = Array2::from_shape_vec(
            (updated_records.len(), 2),
            updated_records
                .iter()
                .flat_map(|r| vec![r.name.clone(), r.num_reads.to_string()])
                .collect(),
        )
        .unwrap();
    } else if !all_elements_exist {
        panic!("All the transcripts in the short read quants do not exist in the header file of the BAM file.");
    }

    //define the indices vector
    let first_column: Vec<String> = data_array.column(0).to_owned().to_vec();
    let indices: Vec<usize> = txps_name
        .iter()
        .map(|x| first_column.iter().position(|y| y == x).unwrap_or_default())
        .collect();

    // Rearrange rows based on an arbitrary indices vector
    data_array = data_array.select(ndarray::Axis(0), &indices);

    let second_column: ndarray::Array<_, _> = data_array.column(1).to_owned();
    let second_column_vec: Vec<f64> = second_column
        .iter()
        .map(|count| count.parse::<f64>().expect("Failed to parse as f64"))
        .collect();

    Ok(second_column_vec)
}
