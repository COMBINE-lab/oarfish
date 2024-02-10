use crate::util::oarfish_types::TranscriptInfo;

#[allow(dead_code)]
pub fn bin_transcript_normalize_counts(
    t: &TranscriptInfo,
    num_bins: &u32,
) -> (Vec<f32>, Vec<f32>, usize, Vec<f64>) {
    let num_discarded_read: usize = 0;
    let transcript_len = t.len.get() as u32; //transcript length
    let nbins = *num_bins;
    let bin_size = transcript_len as f32 / nbins as f32;

    // create a vector of bins, where each bin is an
    // std::ops::Range<f32> struct
    let bins: Vec<std::ops::Range<f32>> = (0..nbins)
        .map(|i| {
            if i != (nbins - 1) {
                i as f32 * bin_size..(i + 1) as f32 * bin_size
            } else {
                i as f32 * bin_size..(transcript_len + 1) as f32
            }
        })
        .collect();

    let bin_counts: Vec<f32> = vec![0.0; bins.len()];
    let bin_coverage: Vec<f64> = vec![0.0; bins.len()];
    let bin_lengths: Vec<f32> = bins.iter().map(|range| range.end - range.start).collect();
    /*
    for read in t.ranges.iter() {
        let mut discarded_read_flag: bool = true;

        for (i, bin) in bins.iter().enumerate() {
            let mut bin_inc: f64 = 0.0;
            if (read.start as f32) <= bin.start && (read.end as f32) >= bin.end {
                // read contains the entire bin
                bin_inc = 1.0;
                discarded_read_flag = false;
            } else if (read.start as f32) > bin.start
                && (read.end as f32) > bin.start
                && (read.start as f32) <= bin.end
                && (read.end as f32) <= bin.end
            {
                // read is contained within the bin
                bin_inc = (read.end - read.start) as f64 / (bin.end - bin.start) as f64;
                discarded_read_flag = false;
            } else if (read.start as f32) >= bin.start && (read.start as f32) < bin.end {
                // read starts in this bin, but doesn't end in it
                bin_inc = (bin.end - read.start as f32) as f64 / (bin.end - bin.start) as f64;
                discarded_read_flag = false;
            } else if (read.end as f32) > bin.start && (read.end as f32) <= bin.end {
                // read ends in this bin, but doesn't start in it
                bin_inc = (read.end as f32 - bin.start) as f64 / (bin.end - bin.start) as f64;
                discarded_read_flag = false;
            }
            bin_counts[i] += bin_inc as f32;
            bin_coverage[i] = bin_inc;
        }

        if discarded_read_flag {
            num_discarded_read += 1;
        }
    }
    */

    (bin_counts, bin_lengths, num_discarded_read, bin_coverage)
}
