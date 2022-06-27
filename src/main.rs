use clap::Parser;

use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufReader},
};

use bio_types::strand::Strand;
use bio_types::annot::spliced::Spliced;
use noodles_gtf as gtf;
use noodles_gtf::record::Strand as NoodlesStrand;

use bio_types::annot::contig::Contig;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Name of the person to greet
    #[clap(short, long, value_parser)]
    name: String,

    /// Number of times to greet
    #[clap(short, long, value_parser, default_value_t = 1)]
    count: u8,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    let mut reader = File::open(args.name)
        .map(BufReader::new)
        .map(gtf::Reader::new)?;
    let mut evec = Vec::new();
    let mut tvec = Vec::new();
    let mut tmap = HashMap::new();

    for result in reader.records() {
        let record = result?;
        match record.ty() {
            "exon" => {
                let s: isize = usize::from(record.start()) as isize;
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

    println!("parsed {} exons", evec.len());
    println!("parsed {} transcripts", tvec.len());
    Ok(())
}
