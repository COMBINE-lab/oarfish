
use crate::variables::{TranscriptInfo, Record, EMInfo};
use std::{
    fs::OpenOptions,
    fs::File,
    io::{self, BufWriter, Write},
};
use csv::ReaderBuilder;
use ndarray::Array2;
use noodles_sam as sam;
use sam::record::data::field::tag;

pub fn bin_transcript_decision_rule(t: &TranscriptInfo, num_bins: &u32) -> (Vec<u32>, Vec<f32>, usize, Vec<f64>) {

    let mut num_discarded_read: usize = 0;
    let transcript_len = t.len.get() as u32; //transcript length
    let nbins = *num_bins;
    let bin_size = transcript_len as f32 / nbins as f32;
    let bins: Vec<std::ops::Range<f32>> = (0..nbins)
        .map(|i| if i != (nbins - 1) {i as f32 * bin_size..(i + 1) as f32 * bin_size} else {i as f32 * bin_size..transcript_len as f32})
        .collect();

    let mut bin_counts: Vec<u32> = vec![0; bins.len()];
    let bin_lengths: Vec<f32> = bins.iter().map(|range| range.end - range.start).collect();
    let mut bin_coverage: Vec<f64> = vec![0.0; bins.len()];

    let dr_threshold = (bins[0].end - bins[0].start)/2.0;


    for read in t.ranges.iter(){

        let mut discarded_read_flag: bool = true;
        let mut coverage_temp = 0.0;

        for (i, bin) in bins.iter().enumerate(){
            if read.start as f32 <= bin.start && read.end as f32 >= bin.end{
                bin_counts[i] = bin_counts[i] + 1;
                discarded_read_flag = false;
                bin_coverage[i] = 1.0;
            } 
            else if (read.start as f32 >= bin.start) && (read.end as f32 >= bin.start) && ((read.start as f32) < bin.end) && ((read.end as f32) < bin.end)  {
                if (read.end as f32 - read.start as f32) >= dr_threshold {
                    bin_counts[i] = bin_counts[i] + 1;
                    discarded_read_flag = false;

                    coverage_temp = (read.end - read.start) as f64/ (bin.end - bin.start) as f64;
                    bin_coverage[i] = if bin_coverage[i] < coverage_temp {coverage_temp} else {bin_coverage[i]};
                }
            }
            else if read.start as f32 >= bin.start && (read.start as f32) < bin.end  {
                if (bin.end as f32 - read.start as f32) >= dr_threshold {
                    bin_counts[i] = bin_counts[i] + 1;
                    discarded_read_flag = false;

                    coverage_temp = (bin.end - read.start as f32) as f64/ (bin.end - bin.start) as f64;
                    bin_coverage[i] = if bin_coverage[i] < coverage_temp {coverage_temp} else {bin_coverage[i]};
                }
            }
            else if read.end as f32 > bin.start && read.end as f32 <= bin.end  {
                if (read.end as f32 - bin.start) >= dr_threshold {
                    bin_counts[i] = bin_counts[i] + 1;
                    discarded_read_flag = false;

                    coverage_temp = (read.end as f32 - bin.start) as f64/ (bin.end - bin.start) as f64;
                    bin_coverage[i] = if bin_coverage[i] < coverage_temp {coverage_temp} else {bin_coverage[i]};
                }
            }
            
        }

        if discarded_read_flag {
            num_discarded_read = num_discarded_read + 1;
        }
    }

 (bin_counts, bin_lengths, num_discarded_read, bin_coverage)

}


pub fn bin_transcript_normalize_counts(t: &TranscriptInfo, num_bins: &u32) -> (Vec<f32>, Vec<f32>, usize, Vec<f64>) {

    let mut num_discarded_read: usize = 0;
    let transcript_len = t.len.get() as u32; //transcript length
    let nbins = *num_bins;
    let bin_size = transcript_len as f32 / nbins as f32;
    let bins: Vec<std::ops::Range<f32>> = (0..nbins)
    .map(|i| if i != (nbins - 1) {i as f32 * bin_size..(i + 1) as f32 * bin_size} else {i as f32 * bin_size..transcript_len as f32})
        .collect();

    let mut bin_counts: Vec<f32> = vec![0.0; bins.len()];
    let bin_lengths: Vec<f32> = bins.iter().map(|range| range.end - range.start).collect();
    let mut bin_coverage: Vec<f64> = vec![0.0; bins.len()];


    for read in t.ranges.iter(){

        let mut discarded_read_flag: bool = true;
        let mut coverage_temp = 0.0;

        for (i, bin) in bins.iter().enumerate(){
            if (read.start as f32) <= bin.start && (read.end as f32) >= bin.end{
                bin_counts[i] = bin_counts[i] + 1.0;
                discarded_read_flag = false;
                bin_coverage[i] = 1.0;
            } 
            else if (read.start as f32) >= bin.start && (read.end as f32) >= bin.start && (read.start as f32) < bin.end && (read.end as f32) < bin.end  {
                    bin_counts[i] = bin_counts[i] + ((read.end - read.start) as f32 / (bin.end - bin.start)as f32);
                    discarded_read_flag = false;

                    coverage_temp = (read.end - read.start) as f64/ (bin.end - bin.start) as f64;
                    bin_coverage[i] = if bin_coverage[i] < coverage_temp {coverage_temp} else {bin_coverage[i]};
            }
            else if (read.start as f32) >= bin.start && (read.start as f32) < bin.end  {
                    bin_counts[i] = bin_counts[i] + ((bin.end - read.start as f32) as f32/ (bin.end - bin.start) as f32);
                    discarded_read_flag = false;

                    coverage_temp = (bin.end - read.start as f32) as f64/ (bin.end - bin.start) as f64;
                    bin_coverage[i] = if bin_coverage[i] < coverage_temp {coverage_temp} else {bin_coverage[i]};

            }
            else if (read.end as f32) > bin.start && (read.end as f32) <= bin.end  {
                    bin_counts[i] = bin_counts[i] + ((read.end as f32 - bin.start) as f32 / (bin.end - bin.start) as f32);
                    discarded_read_flag = false;

                    coverage_temp = (read.end as f32 - bin.start) as f64/ (bin.end - bin.start) as f64;
                    bin_coverage[i] = if bin_coverage[i] < coverage_temp {coverage_temp} else {bin_coverage[i]};
                }
        }

        if discarded_read_flag {
            num_discarded_read = num_discarded_read + 1;
        }
            
    }

 (bin_counts, bin_lengths, num_discarded_read, bin_coverage)

}

pub fn factorial_ln(n: u32) -> f64 {
    if n == 0 {
        0.0
    } else {
        (1..=n).map(|x| (x as f64).ln()).sum()
    }
}

//this part is taken from dev branch
pub fn write_out_count(
    output: String,
    header: &noodles_sam::header::Header,
    txps: &Vec<TranscriptInfo>,
    counts: &[f64],
) -> io::Result<()> {
    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    writeln!(writer, "tname\tcoverage\tlen\tnum_reads").expect("Couldn't write to output file.");
    // loop over the transcripts in the header and fill in the relevant
    // information here.
            
    for (i, (rseq, rmap)) in header.reference_sequences().iter().enumerate() {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            rseq,
            txps[i].coverage,
            rmap.length(),
            counts[i]
        )
        .expect("Couldn't write to output file.");
    }

    Ok(())
}


pub fn write_out_probs(
    cdf_output: String, 
    cov_prob_output: String,
    emi: &EMInfo,
    txps_name: &Vec<String>,
) -> io::Result<()> {
    let write_cdf = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(cdf_output)
        .expect("Couldn't create output file");
    let mut writer_cdf = BufWriter::new(write_cdf);

    writeln!(writer_cdf, "Txps_Name\tCDF_Values").expect("Couldn't write to output file.");
    for (i, txp) in txps_name.iter().enumerate() {
        let cdf_values: String = emi.txp_info[i].coverage_prob.iter()
            .map(|value| value.to_string())
            .collect::<Vec<String>>()
            .join("\t");

        writeln!(
            writer_cdf,
            "{}\t{}",
            *txp,
            cdf_values,
        )
        .expect("Couldn't write to output file.");
    }



    let write_coverage_prob = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(cov_prob_output)
        .expect("Couldn't create output file");
    let mut writer_coverage_prob = BufWriter::new(write_coverage_prob);

    writeln!(writer_coverage_prob, "Read\tTxp\tread_start\tread_end\tAS\tAS_prob\tcoverage_prob").expect("Couldn't write to output file.");
    for (alns, probs, cov_probs) in emi.eq_map.iter() {
        for (a, (p, cov_p)) in alns.iter().zip(probs.iter().zip(cov_probs.iter())) {
            let target_id = a.reference_sequence_id().unwrap();
            let txp_n: String = txps_name[target_id].clone();
            let mut rstring: String = "None".to_string();
            if let Some(rname) = a.read_name() {
                rstring = <noodles_sam::record::read_name::ReadName as AsRef<str>>::as_ref(rname).to_owned();
            }
            let read_start = a.alignment_start().unwrap().get() as u32;
            let read_end = a.alignment_end().unwrap().get() as u32;
            let as_value: &sam::record::data::field::Value = a.data().get(&tag::ALIGNMENT_SCORE).expect("could not get value");
            let as_val = as_value.as_int().unwrap() as i32;
            let as_prob = *p;

            writeln!(
                writer_coverage_prob,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                rstring,
                txp_n,
                read_start,
                read_end,
                as_val,
                as_prob,
                *cov_p
            )
            .expect("Couldn't write to output file.");
        }
    }
    Ok(())
}


//this part is taken from dev branch
pub fn write_out_stat(
    output: String,
    header: &noodles_sam::header::Header,
    num_alignments: &[usize],
    num_discarded_alignments: &[usize],
    num_reads: &usize,
    num_discarded_reads_decision_rule: &usize,
    num_discarded_reads_em: &usize,
    counts: &[f64],
) -> io::Result<()> {
    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);

    writeln!(writer, "tname\tnum_alignments\tnum_discarded_alignments\tnum_reads\tnum_discarded_reads_DR\tnum_discarded_reads_EM\tnum_NOTdiscarded_reads\tcount").expect("Couldn't write to output file.");
    
    let num_notdiscarded_reads = if *num_reads > (num_discarded_reads_decision_rule + num_discarded_reads_em) {num_reads - (num_discarded_reads_decision_rule + num_discarded_reads_em)} else {0};

    writeln!(writer, "overall\t{}\t{}\t{}\t{}\t{}\t{}\t{}", 
            num_alignments.iter().cloned().sum::<usize>(), 
            num_discarded_alignments.iter().cloned().sum::<usize>(),  
            num_reads,
            num_discarded_reads_decision_rule,
            num_discarded_reads_em,
            num_notdiscarded_reads,
            counts.iter().cloned().sum::<f64>()).expect("Couldn't write to output file.");

    writeln!(writer, "tname\tnum_alignments\tnum_discarded_alignments\tcount").expect("Couldn't write to output file.");
    // loop over the transcripts in the header and fill in the relevant
    // information here.
            
    for (i, (_rseq, _rmap)) in header.reference_sequences().iter().enumerate() {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            _rseq,
            num_alignments[i],
            num_discarded_alignments[i],
            counts[i]
        )
        .expect("Couldn't write to output file.");
    }

    Ok(())
}




pub fn short_quant_vec(short_read_path: String, txps_name: Vec<String>) -> Vec<f64> {

    let file = File::open(short_read_path).expect("Failed to open file");

    // Read data from a CSV file (replace with your data source)
    let mut rdr = ReaderBuilder::new().has_headers(true).delimiter(b'\t').from_reader(file);

    // Deserialize CSV records into a Vec of Records
    let records: Vec<Record> = rdr.deserialize().collect::<Result<Vec<Record>, csv::Error>>()
                               .unwrap_or_else(|err| {
                                   eprintln!("Failed to deserialize CSV records: {}", err);
                                   std::process::exit(1);
                               })
                               .into_iter()
                               .collect();

    // Convert the Vec of Records into an ndarray::Array2
    let mut data_array = Array2::from_shape_vec((records.len(), 2), records.iter().flat_map(|r| vec![r.Name.clone(), r.NumReads.to_string()]).collect()).unwrap();

    //obtain the first column of the data_array and check if all the elements exist in txps_name
    let first_column: Vec<String> = data_array.column(0).to_owned().to_vec();
    let all_elements_exist = first_column.iter().all(|element| txps_name.contains(element));

    if all_elements_exist && (first_column.len() != txps_name.len()){

        let txps_not_in_first_column: Vec<String> = txps_name
        .iter()
        .filter(|&x| !first_column.contains(x))
        .cloned()
        .collect();

        
        // Create new records for the missing elements and append them to the records vector
        let missing_records: Vec<Record> = txps_not_in_first_column
            .iter()
            .map(|name| Record {
                Name: name.clone(),
                Length: 0,
                EffectiveLength: 0.0,
                TPM: 0.0,
                NumReads: 0.0,
            })
            .collect();

        let mut updated_records = records.clone();
        updated_records.extend(missing_records);

        // Update the data_array with the new records
        data_array = Array2::from_shape_vec((updated_records.len(), 2), updated_records.iter().flat_map(|r| vec![r.Name.clone(), r.NumReads.to_string()]).collect()).unwrap();

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

    let second_column: ndarray::Array<_,_> = data_array.column(1).to_owned();
    let second_column_vec: Vec<f64> = second_column.iter().map(|count| count.parse::<f64>().expect("Failed to parse as f64")).collect();

    second_column_vec

}