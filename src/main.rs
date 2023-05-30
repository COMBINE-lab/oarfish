use clap::Parser;

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufReader}, num::NonZeroUsize,
};

use bio_types::annot::loc::Loc;
use bio_types::annot::spliced::Spliced;
use bio_types::strand::Strand;
use noodles_gtf as gtf;
use noodles_gtf::record::Strand as NoodlesStrand;
use noodles_bam as bam;

use bio_types::annot::contig::Contig;

use coitrees::{COITree, IntervalNode};

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[clap(short, long, value_parser)]
    alignments: String,
}

#[derive(Debug, Eq, PartialEq)]
struct TranscriptInfo {
    len: NonZeroUsize,
    num_compat: u32,
}

impl TranscriptInfo {
    fn new() -> Self {
        Self { len: NonZeroUsize::new(0).unwrap(), num_compat: 0 }
    }
    fn with_len(len: NonZeroUsize) -> Self {
        Self { len, num_compat: 0 }
    }
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    
    let mut reader = File::open(args.alignments)
        .map(BufReader::new)
        .map(bam::Reader::new)?;

    let header = reader.read_header()?;

    for (prog, pmap) in header.programs().iter()  {
        println!("program: {}", prog);
    }

    let mut txps: Vec<TranscriptInfo> = Vec::with_capacity(header.reference_sequences().len());

    for (rseq, rmap) in header.reference_sequences().iter()  {
        // println!("ref: {}, rmap : {:?}", rseq, rmap.length());
        txps.push(TranscriptInfo::with_len(rmap.length()));
    }

    //let mut rmap = HashMap<usize, ::new();
    //
    let mut prev_read = String::new();
    let mut records_for_read = vec![];

    for result in reader.records(&header) {
        let record = result?;
        let record_copy = record.clone();
        if let Some(rname) = record.read_name() {
            let rstring: String = <noodles_sam::record::read_name::ReadName as AsRef<str>>::as_ref(rname).to_owned();
            if prev_read == rstring {
                records_for_read.push(record_copy);
                if let Some(ref_id) = record.reference_sequence_id() {
                    txps[ref_id].num_compat += 1;
                }
            } else {
                if !prev_read.is_empty() {
                    //println!("the previous read had {} mappings", records_for_read.len());
                    records_for_read.clear();
                }
                prev_read = rstring;
                records_for_read.push(record_copy);
            }
        }
    }

    let mut num_alive = 0;
    for txp in txps.iter() {
       if txp.num_compat > 0 {
            num_alive += 1;
        } 
    }
   
    println!("Number of transcripts with > 0 compatible reads : {}", num_alive);

    Ok(())
}

#[allow(unused)]
fn main_old() -> io::Result<()> {
    let args = Args::parse();

    let mut reader = File::open(args.alignments)
        .map(BufReader::new)
        .map(gtf::Reader::new)?;
    let mut evec = Vec::new();
    let mut tvec = Vec::new();
    let mut tmap = HashMap::new();

    for result in reader.records() {
        let record = result?;
        match record.ty() {
            "exon" => {
                let s: isize = (usize::from(record.start()) as isize) - 1;
                let e: isize = usize::from(record.end()) as isize;
                let l: usize = (e - s).try_into().unwrap();
                let mut t = String::new();
                for e in record.attributes().iter() {
                    if e.key() == "transcript_id" {
                        t = e.value().to_owned();
                    }
                }

                let ni = tmap.len();
                let tid = *tmap.entry(t.clone()).or_insert(ni);

                // if this is what we just inserted
                if ni == tid {
                    tvec.push(Spliced::new(0, 1, 1, Strand::Forward));
                }

                let strand = match record.strand().unwrap() {
                    NoodlesStrand::Forward => Strand::Forward,
                    NoodlesStrand::Reverse => Strand::Reverse,
                };
                let c = Contig::new(tid, s, l, strand);
                evec.push(c);
            }
            "transcript" => {
                let mut t = String::new();
                for e in record.attributes().iter() {
                    if e.key() == "transcript_id" {
                        t = e.value().to_owned();
                    }
                }
                let ni = tmap.len();
                let tid = *tmap.entry(t.clone()).or_insert(ni);

                // if this is what we just inserted
                if ni == tid {
                    tvec.push(Spliced::new(0, 1, 1, Strand::Forward));
                }
            }
            _ => {}
        }
    }

    let mut txp_to_exon = HashMap::new();

    let mut l = 0;
    let mut max_len = 0;
    for (i, e) in evec.iter().enumerate() {
        let mut v = txp_to_exon.entry(e.refid()).or_insert(vec![]);
        v.push(i);
        l = v.len();
        if l > max_len {
            max_len = l;
        }
    }

    let mut txp_features: HashMap<usize, _> = HashMap::new();

    for (k, v) in txp_to_exon.iter_mut() {
        let strand = evec[v[0]].strand();
        v.sort_unstable_by_key(|x| evec[*x as usize].start());
        let s = evec[v[0]].start();
        let starts: Vec<usize> = v
            .iter()
            .map(|e| (evec[*e as usize].start() - s) as usize)
            .collect();
        let lens: Vec<usize> = v.iter().map(|e| evec[*e as usize].length()).collect();
        println!("lens = {:?}, starts = {:?}", lens, starts);
        txp_features.insert(
            **k,
            Spliced::with_lengths_starts(k, s, &lens, &starts, strand).unwrap(),
        );
    }

    let interval_vec: Vec<IntervalNode::<usize, usize>> = evec.iter().enumerate().map(|(i, e)|
                                                                                      IntervalNode::new(e.start() as i32,
                                                                                      (e.start() + e.length() as isize) as i32,
                                                                                      i)).collect();
    let ct = COITree::new(interval_vec);

    println!("parsed {} exons", evec.len());
    println!("parsed {} transcripts", tvec.len());
    println!("max exon transcript had {} exons", max_len);
    Ok(())
}
