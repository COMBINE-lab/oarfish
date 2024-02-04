use crate::util::oarfish_types::TranscriptInfo;

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
        let mut coverage_temp;

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
