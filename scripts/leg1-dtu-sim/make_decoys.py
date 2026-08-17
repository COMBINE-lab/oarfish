#!/usr/bin/env python3
"""Genome-scale SIRV-O analog: fabricate plausible-but-false decoy isoforms
at expressed loci by perturbing real expressed transcripts' splice structure.

For each selected gene (origin transcript tpm >= --min-tpm), one decoy is
created from its most-expressed multi-exon transcript by one perturbation,
drawn uniformly from:
  skip    - drop one random internal exon           (needs >= 3 exons)
  alt3    - drop the last exon (shorter 3' end)     (needs >= 2 exons)
  alt5    - drop the first exon (shorter 5' end)    (needs >= 2 exons)
  extend  - extend one internal exon 60-200 nt into the following intron

Outputs: decoys.gtf (records), decoys.fasta (spliced sequences for mapping),
decoys.tsv (decoy_id, origin transcript, gene, type). Deterministic given
--seed.
"""
import argparse
import random
import re
from collections import defaultdict

TID = re.compile(r'transcript_id "([^"]+)"')
GID = re.compile(r'gene_id "([^"]+)"')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--genome', required=True)
    ap.add_argument('--abundance', required=True, help='target_id\ttpm\tcell')
    ap.add_argument('--n-decoys', type=int, default=20000)
    ap.add_argument('--min-tpm', type=float, default=2.0)
    ap.add_argument('--seed', type=int, default=20260818)
    ap.add_argument('--out-prefix', required=True)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    tpm = {}
    with open(a.abundance) as fh:
        next(fh)
        for line in fh:
            p = line.split('\t')
            if len(p) >= 2:
                tpm[p[0]] = float(p[1])

    # transcript structures
    tx = {}
    for line in open(a.gtf):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        mt, mg = TID.search(f[8]), GID.search(f[8])
        if not (mt and mg):
            continue
        t = mt.group(1)
        rec = tx.setdefault(t, [f[0], f[6], mg.group(1), []])
        rec[3].append((int(f[3]), int(f[4])))
    for t in tx:
        tx[t][3].sort()

    # per gene: best expressed multi-exon transcript
    best = {}
    for t, (chrom, strand, gene, exons) in tx.items():
        v = tpm.get(t, 0.0)
        if v < a.min_tpm or len(exons) < 2:
            continue
        if gene not in best or v > best[gene][1]:
            best[gene] = (t, v)
    genes = sorted(best)
    rng.shuffle(genes)
    genes = genes[:a.n_decoys]
    print(f'eligible genes: {len(best)}; making decoys for {len(genes)}')

    # genome (memory: ~3.2 GB, fine)
    seqs = {}
    name = None
    chunks = defaultdict(list)
    for line in open(a.genome):
        if line.startswith('>'):
            name = line[1:].split()[0]
        else:
            chunks[name].append(line.strip())
    for k, v in chunks.items():
        seqs[k] = ''.join(v)
    del chunks
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')

    gtf_out = open(f'{a.out_prefix}.gtf', 'w')
    fa_out = open(f'{a.out_prefix}.fasta', 'w')
    man = open(f'{a.out_prefix}.tsv', 'w')
    man.write('decoy\torigin\tgene\ttype\n')
    made = 0
    kinds = ['skip', 'alt3', 'alt5', 'extend']
    for i, g in enumerate(genes):
        origin, _ = best[g]
        chrom, strand, gene, exons = tx[origin]
        if chrom not in seqs:
            continue
        options = [k for k in kinds
                   if (k == 'skip' and len(exons) >= 3)
                   or (k in ('alt3', 'alt5') and len(exons) >= 2)
                   or (k == 'extend' and len(exons) >= 2)]
        if not options:
            continue
        kind = rng.choice(options)
        ex = list(exons)
        if kind == 'skip':
            ex.pop(rng.randrange(1, len(ex) - 1))
        elif kind == 'alt3':
            # genomic order; 3' end depends on strand
            if strand == '+':
                ex = ex[:-1]
            else:
                ex = ex[1:]
        elif kind == 'alt5':
            if strand == '+':
                ex = ex[1:]
            else:
                ex = ex[:-1]
        else:  # extend an internal exon into the following intron
            j = rng.randrange(0, len(ex) - 1)
            gap = ex[j + 1][0] - ex[j][1] - 1
            if gap < 80:
                continue
            add = min(rng.randrange(60, 201), gap - 20)
            ex[j] = (ex[j][0], ex[j][1] + add)
        did = f'DECOY{made:05d}_{kind}'
        seq = ''.join(seqs[chrom][s - 1:e] for s, e in ex)
        if strand == '-':
            seq = seq.translate(comp)[::-1]
        if len(seq) < 200:
            continue
        attrs = f'gene_id "{gene}"; transcript_id "{did}"; decoy_origin "{origin}";'
        lo, hi = ex[0][0], ex[-1][1]
        gtf_out.write(f'{chrom}\tdecoy\ttranscript\t{lo}\t{hi}\t.\t{strand}\t.\t{attrs}\n')
        for s, e in ex:
            gtf_out.write(f'{chrom}\tdecoy\texon\t{s}\t{e}\t.\t{strand}\t.\t{attrs}\n')
        fa_out.write(f'>{did}\n{seq}\n')
        man.write(f'{did}\t{origin}\t{gene}\t{kind}\n')
        made += 1
    print(f'made {made} decoys')


if __name__ == '__main__':
    main()
