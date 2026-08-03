//! juncprobe -- fast tests for the splice-junction endpoint signal.
//!
//! phase0   : likelihood ratio that an alignment boundary coincides with an
//!            INTERNAL splice junction of the transcript, split by whether that
//!            transcript actually generated the read (simulator read names).
//! reweight : Option 1. Take oarfish's per-read posteriors (.prob), down-weight
//!            junction-terminated alignments, renormalise within the read, and
//!            measure whether posterior mass moves toward the true source.
//!            No EM -- this asks only whether the signal targets real errors.

use anyhow::{anyhow, Context, Result};
use clap::{Parser, Subcommand};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    Phase0 {
        #[arg(long)]
        bam: String,
        #[arg(long)]
        junctions: String,
        #[arg(long, default_value_t = 16)]
        threads: usize,
    },
    Reweight {
        #[arg(long)]
        bam: String,
        #[arg(long)]
        junctions: String,
        #[arg(long)]
        prob: String,
        #[arg(long, default_value_t = 16)]
        threads: usize,
        /// tolerance windows (nt) to evaluate
        #[arg(long, value_delimiter = ',', default_values_t = vec![0u32,2,5,10])]
        windows: Vec<u32>,
        /// log-penalties to evaluate
        #[arg(long, value_delimiter = ',', default_values_t = vec![0.5f64,1.0,1.5,2.3])]
        lambdas: Vec<f64>,
    },
    /// Information ceiling: for reads the EM misassigns, is the true source's
    /// alignment even distinguishable from the winner's by alignment score?
    Ceiling {
        #[arg(long)]
        bam: String,
        #[arg(long)]
        junctions: String,
        #[arg(long)]
        prob: String,
        #[arg(long, default_value_t = 16)]
        threads: usize,
    },
    /// Option 2: rebuild the conditionals from the .prob, apply the junction
    /// term, re-run EM offline, and score against transcript-level truth.
    Em {
        #[arg(long)]
        bam: String,
        #[arg(long)]
        junctions: String,
        #[arg(long)]
        prob: String,
        /// oarfish .quant, used to recover conditionals as p_rt / theta_t
        #[arg(long)]
        quant: String,
        /// truth: <transcript>\t<count>, optional header
        #[arg(long)]
        truth: String,
        #[arg(long, default_value_t = 16)]
        threads: usize,
        #[arg(long, value_delimiter = ',', default_values_t = vec![5u32])]
        windows: Vec<u32>,
        #[arg(long, value_delimiter = ',', default_values_t = vec![0.0f64, 6.0])]
        lambdas: Vec<f64>,
        #[arg(long, default_value_t = 1000)]
        max_iter: usize,
    },
}

/// strip a trailing `.<digits>` accession version
fn strip(n: &str) -> &str {
    if let Some(i) = n.rfind('.') {
        if i + 1 < n.len() && n[i + 1..].bytes().all(|b| b.is_ascii_digit()) {
            return &n[..i];
        }
    }
    n
}

/// NanoSim/TKSM read names encode the source transcript:
/// `NM_015477_2260_aligned_7818918_F_54_99_25` -> `NM_015477`
fn true_source(name: &str) -> Option<&str> {
    for sep in ["_aligned_", "_unaligned_", "_perfect_"] {
        if let Some(i) = name.find(sep) {
            let pre = &name[..i];
            return pre.rfind('_').map(|j| &pre[..j]);
        }
    }
    None
}

fn load_junctions(path: &str) -> Result<FxHashMap<String, Vec<u32>>> {
    let mut m = FxHashMap::default();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let mut it = line.split('\t');
        let tid = it.next().ok_or_else(|| anyhow!("bad junctions line"))?;
        let _len = it.next();
        let js = it.next().unwrap_or("");
        let v: Vec<u32> = if js.is_empty() {
            Vec::new()
        } else {
            js.split(',').filter_map(|x| x.parse().ok()).collect()
        };
        // keys are stored unversioned so they match the stripped BAM ref names
        m.insert(strip(tid).to_string(), v);
    }
    Ok(m)
}

/// nearest distance from `x` to any junction in the sorted slice
#[inline]
fn nearest(js: &[u32], x: i64) -> i64 {
    if js.is_empty() {
        return i64::MAX;
    }
    let i = js.partition_point(|&j| (j as i64) < x);
    let mut best = i64::MAX;
    if i < js.len() {
        best = best.min((js[i] as i64 - x).abs());
    }
    if i > 0 {
        best = best.min((js[i - 1] as i64 - x).abs());
    }
    best
}

struct Bam {
    reader: bam::io::Reader<bgzf::io::MultithreadedReader<File>>,
    /// per BAM reference id -> junctions (empty if unknown/unusable)
    refjunc: Vec<Vec<u32>>,
    refname: Vec<String>,
}

fn open_bam(path: &str, threads: usize, junc: &FxHashMap<String, Vec<u32>>) -> Result<Bam> {
    let f = File::open(path).with_context(|| format!("open {path}"))?;
    let wc = std::num::NonZeroUsize::new(threads.max(1)).unwrap();
    let dec = bgzf::io::MultithreadedReader::with_worker_count(wc, f);
    let mut reader = bam::io::Reader::from(dec);
    let header = reader.read_header()?;
    let mut refjunc = Vec::new();
    let mut refname = Vec::new();
    for (name, _) in header.reference_sequences() {
        let n = String::from_utf8_lossy(name.as_ref()).to_string();
        let key = strip(&n).to_string();
        refjunc.push(junc.get(&key).cloned().unwrap_or_default());
        refname.push(key);
    }
    Ok(Bam {
        reader,
        refjunc,
        refname,
    })
}

/// reference bases consumed by the CIGAR
#[inline]
fn ref_span(rec: &bam::Record) -> i64 {
    use noodles_sam::alignment::record::cigar::op::Kind;
    let mut n = 0i64;
    for op in rec.cigar().iter().flatten() {
        match op.kind() {
            Kind::Match
            | Kind::Deletion
            | Kind::Skip
            | Kind::SequenceMatch
            | Kind::SequenceMismatch => n += op.len() as i64,
            _ => {}
        }
    }
    n
}

fn phase0(bam: &str, junctions: &str, threads: usize) -> Result<()> {
    let junc = load_junctions(junctions)?;
    let mut b = open_bam(bam, threads, &junc)?;
    let windows = [0u32, 2, 5, 10, 20];
    // [true/decoy][window]
    let mut hit3 = [[0u64; 5]; 2];
    let mut hit5 = [[0u64; 5]; 2];
    let mut tot = [0u64; 2];
    let mut rec = bam::Record::default();
    let mut n = 0u64;
    while b.reader.read_record(&mut rec)? != 0 {
        n += 1;
        let Some(rid) = rec.reference_sequence_id().transpose()? else {
            continue;
        };
        let js = &b.refjunc[rid];
        if js.is_empty() {
            continue;
        }
        let Some(pos) = rec.alignment_start().transpose()? else {
            continue;
        };
        let a = usize::from(pos) as i64 - 1; // 0-based
        let e = a + ref_span(&rec);
        let name = String::from_utf8_lossy(rec.name().unwrap_or_default().as_ref()).to_string();
        let Some(src) = true_source(&name) else {
            continue;
        };
        let k = if src == b.refname[rid] { 0 } else { 1 };
        tot[k] += 1;
        let d3 = nearest(js, e);
        let d5 = nearest(js, a);
        for (wi, &w) in windows.iter().enumerate() {
            if d3 <= w as i64 {
                hit3[k][wi] += 1;
            }
            if d5 <= w as i64 {
                hit5[k][wi] += 1;
            }
        }
    }
    println!(
        "records: {n}   used: true={} decoy={}",
        tot[0], tot[1]
    );
    for (lbl, h) in [("3' END", &hit3), ("5' START", &hit5)] {
        println!("\n{lbl} within w of an INTERNAL junction");
        println!("{:<6} {:<12} {:<12} {}", "w", "P(hit|true)", "P(hit|decoy)", "LR");
        for (wi, &w) in windows.iter().enumerate() {
            let pt = h[0][wi] as f64 / tot[0].max(1) as f64;
            let pd = h[1][wi] as f64 / tot[1].max(1) as f64;
            println!(
                "{:<6} {:<12.4} {:<12.4} {:.2}x",
                w,
                pt,
                pd,
                if pt > 0.0 { pd / pt } else { f64::INFINITY }
            );
        }
    }
    Ok(())
}

fn reweight(
    bam: &str,
    junctions: &str,
    prob: &str,
    threads: usize,
    windows: &[u32],
    lambdas: &[f64],
) -> Result<()> {
    let junc = load_junctions(junctions)?;
    let mut b = open_bam(bam, threads, &junc)?;

    // pass 1: read name -> [(ref id, min endpoint distance to internal junction)]
    let mut per_read: FxHashMap<Box<[u8]>, Vec<(u32, u32)>> = FxHashMap::default();
    let mut rec = bam::Record::default();
    let mut n = 0u64;
    while b.reader.read_record(&mut rec)? != 0 {
        n += 1;
        let Some(rid) = rec.reference_sequence_id().transpose()? else {
            continue;
        };
        let Some(pos) = rec.alignment_start().transpose()? else {
            continue;
        };
        let js = &b.refjunc[rid];
        let d = if js.is_empty() {
            u32::MAX
        } else {
            let a = usize::from(pos) as i64 - 1;
            let e = a + ref_span(&rec);
            nearest(js, a).min(nearest(js, e)).min(u32::MAX as i64) as u32
        };
        let name = rec.name().unwrap_or_default();
        let key: Box<[u8]> = AsRef::<[u8]>::as_ref(name).to_vec().into_boxed_slice();
        per_read.entry(key).or_default().push((rid as u32, d));
    }
    eprintln!("BAM: {n} records, {} reads", per_read.len());

    // pass 2: .prob -- header, transcript names, then per-read posteriors
    let f = File::open(prob).with_context(|| format!("open {prob}"))?;
    let mut r = BufReader::with_capacity(1 << 22, f);
    let mut line = String::new();
    r.read_line(&mut line)?;
    let mut it = line.split_whitespace();
    let n_txp: usize = it.next().ok_or_else(|| anyhow!("bad .prob header"))?.parse()?;
    let n_read: usize = it.next().unwrap_or("0").parse().unwrap_or(0);
    // prob transcript index -> BAM ref id
    let mut bam_ref_of: Vec<u32> = Vec::with_capacity(n_txp);
    let mut name_to_ref: FxHashMap<&str, u32> = FxHashMap::default();
    for (i, nm) in b.refname.iter().enumerate() {
        name_to_ref.insert(nm.as_str(), i as u32);
    }
    let mut names = Vec::with_capacity(n_txp);
    for _ in 0..n_txp {
        line.clear();
        r.read_line(&mut line)?;
        let full = line.trim_end();
        let base = full.split('|').next().unwrap_or(full);
        names.push(strip(base).to_string());
    }
    let mut unmapped = 0usize;
    for nm in &names {
        match name_to_ref.get(nm.as_str()) {
            Some(&i) => bam_ref_of.push(i),
            None => {
                bam_ref_of.push(u32::MAX);
                unmapped += 1;
            }
        }
    }
    eprintln!(
        ".prob: {n_txp} transcripts ({unmapped} not resolvable to BAM refs), {n_read} reads"
    );

    // accumulate: baseline and reweighted posterior mass on the true source
    let nw = windows.len();
    let nl = lambdas.len();
    let mut base_mass = 0f64;
    let mut new_mass = vec![0f64; nw * nl];
    let mut reads_scored = 0u64;
    let mut reads_multi = 0u64;
    let mut cand = Vec::<(u32, f64)>::new();

    line.clear();
    while r.read_line(&mut line)? != 0 {
        let s = line.trim_end();
        if s.is_empty() {
            line.clear();
            continue;
        }
        let mut f = s.split('\t');
        let rname = f.next().unwrap_or("");
        let ncand: usize = f.next().unwrap_or("0").parse().unwrap_or(0);
        cand.clear();
        for _ in 0..ncand {
            let (Some(i), Some(p)) = (f.next(), f.next()) else {
                break;
            };
            if let (Ok(i), Ok(p)) = (i.parse::<usize>(), p.parse::<f64>()) {
                if i < bam_ref_of.len() {
                    cand.push((bam_ref_of[i], p));
                }
            }
        }
        let Some(src) = true_source(rname) else {
            line.clear();
            continue;
        };
        let Some(&src_ref) = name_to_ref.get(src) else {
            line.clear();
            continue;
        };
        let hits = per_read.get(rname.as_bytes());
        reads_scored += 1;
        if cand.len() > 1 {
            reads_multi += 1;
        }
        let tot: f64 = cand.iter().map(|c| c.1).sum();
        if tot <= 0.0 {
            line.clear();
            continue;
        }
        base_mass += cand
            .iter()
            .filter(|c| c.0 == src_ref)
            .map(|c| c.1)
            .sum::<f64>()
            / tot;

        for (wi, &w) in windows.iter().enumerate() {
            for (li, &lam) in lambdas.iter().enumerate() {
                let mut sum = 0f64;
                let mut on_src = 0f64;
                for &(rid, p) in &cand {
                    let mut q = p;
                    if let Some(hs) = hits {
                        if let Some(&(_, d)) = hs.iter().find(|(r2, _)| *r2 == rid) {
                            if d as u64 <= w as u64 {
                                q *= (-lam).exp();
                            }
                        }
                    }
                    sum += q;
                    if rid == src_ref {
                        on_src += q;
                    }
                }
                if sum > 0.0 {
                    new_mass[wi * nl + li] += on_src / sum;
                }
            }
        }
        line.clear();
    }

    println!(
        "reads scored: {reads_scored} (multi-candidate: {reads_multi})"
    );
    println!(
        "baseline mean posterior mass on TRUE source: {:.6}",
        base_mass / reads_scored.max(1) as f64
    );
    println!("\nreweighted (higher is better):");
    print!("{:<8}", "w\\lam");
    for l in lambdas {
        print!("{:<12.2}", l);
    }
    println!();
    for (wi, &w) in windows.iter().enumerate() {
        print!("{:<8}", w);
        for (li, _) in lambdas.iter().enumerate() {
            let v = new_mass[wi * nl + li] / reads_scored.max(1) as f64;
            print!("{:<+12.6}", v - base_mass / reads_scored.max(1) as f64);
        }
        println!();
    }
    println!("\n(cells are DELTA vs baseline mean mass on the true source)");
    Ok(())
}

/// average ranks, ties shared
fn ranks(v: &[f64]) -> Vec<f64> {
    let n = v.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_unstable_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap());
    let mut out = vec![0f64; n];
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && v[idx[j]] == v[idx[i]] {
            j += 1;
        }
        let r = (i + j - 1) as f64 / 2.0 + 1.0;
        for k in i..j {
            out[idx[k]] = r;
        }
        i = j;
    }
    out
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut cov = 0.0;
    let mut vx = 0.0;
    let mut vy = 0.0;
    for i in 0..x.len() {
        let a = x[i] - mx;
        let b = y[i] - my;
        cov += a * b;
        vx += a * a;
        vy += b * b;
    }
    if vx > 0.0 && vy > 0.0 {
        cov / (vx.sqrt() * vy.sqrt())
    } else {
        f64::NAN
    }
}

fn read_truth_map(path: &str) -> Result<FxHashMap<String, f64>> {
    let mut m = FxHashMap::default();
    for (i, line) in BufReader::new(File::open(path)?).lines().enumerate() {
        let line = line?;
        let mut it = line.split_whitespace();
        let (Some(t), Some(c)) = (it.next(), it.next()) else {
            continue;
        };
        match c.parse::<f64>() {
            Ok(v) => {
                *m.entry(strip(t).to_string()).or_insert(0.0) += v;
            }
            Err(_) if i == 0 => continue, // header
            Err(e) => return Err(anyhow!("{path}:{}: {e}", i + 1)),
        }
    }
    Ok(m)
}

#[allow(clippy::too_many_arguments)]
fn em_mode(
    bam: &str,
    junctions: &str,
    prob: &str,
    quant: &str,
    truth_path: &str,
    threads: usize,
    windows: &[u32],
    lambdas: &[f64],
    max_iter: usize,
) -> Result<()> {
    let junc = load_junctions(junctions)?;
    let mut b = open_bam(bam, threads, &junc)?;

    // pass 1: per read, the min endpoint-to-internal-junction distance per ref
    let mut per_read: FxHashMap<Box<[u8]>, Vec<(u32, u32)>> = FxHashMap::default();
    let mut rec = bam::Record::default();
    while b.reader.read_record(&mut rec)? != 0 {
        let Some(rid) = rec.reference_sequence_id().transpose()? else {
            continue;
        };
        let Some(pos) = rec.alignment_start().transpose()? else {
            continue;
        };
        let js = &b.refjunc[rid];
        let d = if js.is_empty() {
            u32::MAX
        } else {
            let a = usize::from(pos) as i64 - 1;
            let e = a + ref_span(&rec);
            nearest(js, a).min(nearest(js, e)).min(u32::MAX as i64) as u32
        };
        let key: Box<[u8]> = AsRef::<[u8]>::as_ref(rec.name().unwrap_or_default())
            .to_vec()
            .into_boxed_slice();
        let v = per_read.entry(key).or_default();
        match v.iter_mut().find(|(r, _)| *r == rid as u32) {
            Some(slot) => slot.1 = slot.1.min(d),
            None => v.push((rid as u32, d)),
        }
    }
    eprintln!("BAM pass done: {} reads", per_read.len());

    // theta from the .quant, keyed by stripped transcript name
    let mut theta_by_name: FxHashMap<String, f64> = FxHashMap::default();
    {
        let mut rdr = BufReader::new(File::open(quant)?);
        let mut line = String::new();
        rdr.read_line(&mut line)?; // header
        line.clear();
        while rdr.read_line(&mut line)? != 0 {
            let s = line.trim_end();
            let mut f = s.split('\t');
            if let (Some(nm), Some(_len), Some(nr)) = (f.next(), f.next(), f.next()) {
                let base = nm.split('|').next().unwrap_or(nm);
                *theta_by_name.entry(strip(base).to_string()).or_insert(0.0) +=
                    nr.parse::<f64>().unwrap_or(0.0);
            }
            line.clear();
        }
    }

    // .prob
    let f = File::open(prob)?;
    let mut r = BufReader::with_capacity(1 << 22, f);
    let mut line = String::new();
    r.read_line(&mut line)?;
    let n_txp: usize = line
        .split_whitespace()
        .next()
        .ok_or_else(|| anyhow!("bad header"))?
        .parse()?;
    let mut names = Vec::with_capacity(n_txp);
    for _ in 0..n_txp {
        line.clear();
        r.read_line(&mut line)?;
        let full = line.trim_end();
        let base = full.split('|').next().unwrap_or(full);
        names.push(strip(base).to_string());
    }
    let mut name_to_ref: FxHashMap<&str, u32> = FxHashMap::default();
    for (i, nm) in b.refname.iter().enumerate() {
        name_to_ref.insert(nm.as_str(), i as u32);
    }
    let bam_ref_of: Vec<u32> = names
        .iter()
        .map(|n| *name_to_ref.get(n.as_str()).unwrap_or(&u32::MAX))
        .collect();
    let theta0: Vec<f64> = names
        .iter()
        .map(|n| *theta_by_name.get(n.as_str()).unwrap_or(&0.0))
        .collect();

    // flattened conditionals: c_rt propto p_rt / theta_t
    let mut cand_t: Vec<u32> = Vec::new();
    let mut cand_c: Vec<f64> = Vec::new();
    let mut cand_d: Vec<u32> = Vec::new();
    let mut offs: Vec<u64> = vec![0];
    line.clear();
    while r.read_line(&mut line)? != 0 {
        let s = line.trim_end();
        if s.is_empty() {
            line.clear();
            continue;
        }
        let mut f = s.split('\t');
        let rname = f.next().unwrap_or("");
        let nc: usize = f.next().unwrap_or("0").parse().unwrap_or(0);
        let hits = per_read.get(rname.as_bytes());
        let mut added = 0u64;
        for _ in 0..nc {
            let (Some(i), Some(p)) = (f.next(), f.next()) else {
                break;
            };
            let (Ok(i), Ok(p)) = (i.parse::<usize>(), p.parse::<f64>()) else {
                continue;
            };
            if i >= theta0.len() || theta0[i] <= 0.0 || p <= 0.0 {
                continue;
            }
            let rid = bam_ref_of[i];
            let d = hits
                .and_then(|h| h.iter().find(|(x, _)| *x == rid).map(|(_, d)| *d))
                .unwrap_or(u32::MAX);
            cand_t.push(i as u32);
            cand_c.push(p / theta0[i]);
            cand_d.push(d);
            added += 1;
        }
        if added > 0 {
            offs.push(offs.last().unwrap() + added);
        }
        line.clear();
    }
    let n_reads = offs.len() - 1;
    eprintln!(
        "conditionals: {} reads, {} (read,txp) pairs",
        n_reads,
        cand_t.len()
    );

    let truth = read_truth_map(truth_path)?;
    let tvec: Vec<f64> = names
        .iter()
        .map(|n| *truth.get(n.as_str()).unwrap_or(&0.0))
        .collect();
    let rt = ranks(&tvec);

    println!("{:<8} {:<8} {:<10} {:<10}", "w", "lambda", "spearman", "mard");
    for &w in windows {
        for &lam in lambdas {
            let pen = (-lam).exp();
            let mut theta = vec![1.0f64 / n_txp as f64; n_txp];
            let mut next = vec![0.0f64; n_txp];
            for _ in 0..max_iter {
                next.iter_mut().for_each(|v| *v = 0.0);
                for rd in 0..n_reads {
                    let (s, e) = (offs[rd] as usize, offs[rd + 1] as usize);
                    let mut tot = 0.0;
                    for k in s..e {
                        let t = cand_t[k] as usize;
                        let mut c = cand_c[k] * theta[t];
                        if cand_d[k] <= w {
                            c *= pen;
                        }
                        tot += c;
                    }
                    if tot <= 0.0 {
                        continue;
                    }
                    for k in s..e {
                        let t = cand_t[k] as usize;
                        let mut c = cand_c[k] * theta[t];
                        if cand_d[k] <= w {
                            c *= pen;
                        }
                        next[t] += c / tot;
                    }
                }
                std::mem::swap(&mut theta, &mut next);
            }
            let re = ranks(&theta);
            let sp = pearson(&rt, &re);
            let ts: f64 = tvec.iter().sum();
            let es: f64 = theta.iter().sum();
            let scale = if es > 0.0 { ts / es } else { 0.0 };
            let mut mard = 0.0;
            for i in 0..n_txp {
                let (a, bb) = (tvec[i], theta[i] * scale);
                let d = a.abs() + bb.abs();
                if d > 0.0 {
                    mard += (a - bb).abs() / d;
                }
            }
            mard /= n_txp as f64;
            println!("{:<8} {:<8.2} {:<10.4} {:<10.4}", w, lam, sp, mard);
        }
    }
    Ok(())
}


/// Information ceiling. For each read, compare the alignment of the TRUE source
/// against the alignment of the transcript the EM actually favoured. If the two
/// alignments score identically, no per-alignment feature -- junction-based or
/// otherwise -- can separate them: the reads are genuinely uninformative and only
/// abundance can resolve them. This bounds what ANY annotation-derived signal
/// could achieve, independent of its formulation.
fn ceiling(bam: &str, junctions: &str, prob: &str, threads: usize) -> Result<()> {
    let junc = load_junctions(junctions)?;
    let mut b = open_bam(bam, threads, &junc)?;

    // read -> per-ref (alignment score, junction-terminated flag)
    let mut per_read: FxHashMap<Box<[u8]>, Vec<(u32, i32, bool)>> = FxHashMap::default();
    let mut rec = bam::Record::default();
    while b.reader.read_record(&mut rec)? != 0 {
        let Some(rid) = rec.reference_sequence_id().transpose()? else { continue };
        let Some(pos) = rec.alignment_start().transpose()? else { continue };
        let js = &b.refjunc[rid];
        let a = usize::from(pos) as i64 - 1;
        let e = a + ref_span(&rec);
        let hit = !js.is_empty() && (nearest(js, a) <= 5 || nearest(js, e) <= 5);
        // alignment score from the AS tag
        let mut score = i32::MIN;
        if let Some(Ok(v)) = rec.data().get(b"AS") {
            score = match v {
                noodles_sam::alignment::record::data::field::Value::Int32(x) => x,
                noodles_sam::alignment::record::data::field::Value::Int16(x) => x as i32,
                noodles_sam::alignment::record::data::field::Value::Int8(x) => x as i32,
                noodles_sam::alignment::record::data::field::Value::UInt32(x) => x as i32,
                noodles_sam::alignment::record::data::field::Value::UInt16(x) => x as i32,
                noodles_sam::alignment::record::data::field::Value::UInt8(x) => x as i32,
                _ => i32::MIN,
            };
        }
        let key: Box<[u8]> = AsRef::<[u8]>::as_ref(rec.name().unwrap_or_default())
            .to_vec().into_boxed_slice();
        let v = per_read.entry(key).or_default();
        match v.iter_mut().find(|(r, _, _)| *r == rid as u32) {
            Some(slot) => { if score > slot.1 { slot.1 = score; } slot.2 |= hit; }
            None => v.push((rid as u32, score, hit)),
        }
    }
    eprintln!("BAM: {} reads", per_read.len());

    let f = File::open(prob)?;
    let mut r = BufReader::with_capacity(1 << 22, f);
    let mut line = String::new();
    r.read_line(&mut line)?;
    let n_txp: usize = line.split_whitespace().next().unwrap().parse()?;
    let mut names = Vec::with_capacity(n_txp);
    for _ in 0..n_txp {
        line.clear();
        r.read_line(&mut line)?;
        let full = line.trim_end();
        names.push(strip(full.split('|').next().unwrap_or(full)).to_string());
    }
    let mut name_to_ref: FxHashMap<&str, u32> = FxHashMap::default();
    for (i, nm) in b.refname.iter().enumerate() { name_to_ref.insert(nm.as_str(), i as u32); }
    let bam_ref_of: Vec<u32> = names.iter()
        .map(|n| *name_to_ref.get(n.as_str()).unwrap_or(&u32::MAX)).collect();

    // classify misassigned reads
    let (mut n_multi, mut n_mis) = (0u64, 0u64);
    let (mut tie, mut truth_better, mut truth_worse) = (0u64, 0u64, 0u64);
    let (mut tie_junc, mut better_junc) = (0u64, 0u64);
    let mut cand = Vec::<(u32, f64)>::new();
    line.clear();
    while r.read_line(&mut line)? != 0 {
        let s = line.trim_end();
        if s.is_empty() { line.clear(); continue; }
        let mut f = s.split('\t');
        let rname = f.next().unwrap_or("");
        let nc: usize = f.next().unwrap_or("0").parse().unwrap_or(0);
        cand.clear();
        for _ in 0..nc {
            let (Some(i), Some(p)) = (f.next(), f.next()) else { break };
            if let (Ok(i), Ok(p)) = (i.parse::<usize>(), p.parse::<f64>()) {
                if i < bam_ref_of.len() { cand.push((bam_ref_of[i], p)); }
            }
        }
        if cand.len() < 2 { line.clear(); continue; }
        n_multi += 1;
        let Some(src) = true_source(rname) else { line.clear(); continue };
        let Some(&src_ref) = name_to_ref.get(src) else { line.clear(); continue };
        let tot: f64 = cand.iter().map(|c| c.1).sum();
        let on_src: f64 = cand.iter().filter(|c| c.0 == src_ref).map(|c| c.1).sum();
        if tot <= 0.0 || on_src / tot >= 0.5 { line.clear(); continue; }
        n_mis += 1;
        // winner = highest posterior candidate that is not the true source
        let mut win = (u32::MAX, -1.0f64);
        for &(rid, p) in &cand { if rid != src_ref && p > win.1 { win = (rid, p); } }
        let hits = per_read.get(rname.as_bytes());
        let sc = |rid: u32| hits.and_then(|h| h.iter().find(|(x,_,_)| *x==rid).map(|(_,s,j)| (*s,*j)));
        let (Some((s_src, j_src)), Some((s_win, j_win))) = (sc(src_ref), sc(win.0)) else { line.clear(); continue };
        if s_src == i32::MIN || s_win == i32::MIN { line.clear(); continue; }
        if s_src == s_win {
            tie += 1;
            if j_win && !j_src { tie_junc += 1; }
        } else if s_src > s_win {
            truth_better += 1;
            if j_win && !j_src { better_junc += 1; }
        } else { truth_worse += 1; }
        line.clear();
    }
    println!("multi-candidate reads: {n_multi}");
    println!("misassigned (true source has <50% posterior): {n_mis}\n");
    let d = n_mis.max(1) as f64;
    println!("Of misassigned reads, comparing the TRUE source's alignment to the winner's:");
    println!("  {:<44} {:>9}  {:>6.2}%", "TIED score (information-theoretically stuck)", tie, 100.0*tie as f64/d);
    println!("  {:<44} {:>9}  {:>6.2}%", "true source scores HIGHER (recoverable)", truth_better, 100.0*truth_better as f64/d);
    println!("  {:<44} {:>9}  {:>6.2}%", "true source scores LOWER (read is corrupted)", truth_worse, 100.0*truth_worse as f64/d);
    println!("\n  of the TIED cases, junction signal favours truth: {} ({:.2}% of all misassigned)",
             tie_junc, 100.0*tie_junc as f64/d);
    println!("  of the HIGHER cases, junction signal favours truth: {} ({:.2}% of all misassigned)",
             better_junc, 100.0*better_junc as f64/d);
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Phase0 {
            bam,
            junctions,
            threads,
        } => phase0(&bam, &junctions, threads),
        Cmd::Reweight {
            bam,
            junctions,
            prob,
            threads,
            windows,
            lambdas,
        } => reweight(&bam, &junctions, &prob, threads, &windows, &lambdas),
        Cmd::Ceiling {
            bam,
            junctions,
            prob,
            threads,
        } => ceiling(&bam, &junctions, &prob, threads),
        Cmd::Em {
            bam,
            junctions,
            prob,
            quant,
            truth,
            threads,
            windows,
            lambdas,
            max_iter,
        } => em_mode(
            &bam, &junctions, &prob, &quant, &truth, threads, &windows, &lambdas, max_iter,
        ),
    }
}
