
use crate::variables::TranscriptInfo;
use std::{
    fs::OpenOptions,
    io::{self, BufWriter, Write},
};

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
pub fn write_output(
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

pub fn write_read_coverage(
    output: String, 
    read_coverage_probs: &Vec<Vec<f64>>,
) -> io::Result<()> {
    let write = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(output)
        .expect("Couldn't create output file");
    let mut writer = BufWriter::new(write);
    //writeln!(writer, "tname\tnum_alignments\tnum_discarded_alignments\tcount").expect("Couldn't write to output file.");

    for prob in read_coverage_probs.iter() {
        for (i, p) in prob.iter().enumerate(){
            write!(
                writer,
                "{}\t",
                p,
            ).expect("Couldn't write to output file.");
        }
        writeln!(writer).expect("Couldn't write to output file.");
    }

    Ok(())
}