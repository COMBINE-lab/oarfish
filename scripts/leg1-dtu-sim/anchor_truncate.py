#!/usr/bin/env python3
"""dRNA-faithful truncation for TKSM mdf files: trim intervals from the
chain START (transcript 5' end), preserving the 3' end — dRNA sequences from
the poly(A) and loses the 5' side. Target read length ~ lognormal(mu,sigma),
clamped to [50, molecule length]."""
import argparse, math, random, sys
ap = argparse.ArgumentParser()
ap.add_argument('--input', required=True)
ap.add_argument('--output', required=True)
ap.add_argument('--mu', type=float, default=6.4)
ap.add_argument('--sigma', type=float, default=0.7)
ap.add_argument('--seed', type=int, default=1)
a = ap.parse_args()
rng = random.Random(a.seed)
out = open(a.output, 'w')
def flush(hdr, ivs):
    if hdr is None: return
    L = sum(e - s + 1 for _, s, e, _ in ivs)
    T = max(50, min(L, int(math.exp(rng.gauss(a.mu, a.sigma)))))
    drop = L - T
    kept = []
    for chrom, s, e, rest in ivs:
        ln = e - s + 1
        if drop >= ln:
            drop -= ln
            continue
        if drop > 0:
            # trim from the transcript-5' side of this interval:
            # '+' strand exon: 5' side = low coords; '-' strand: high coords
            if rest.startswith('-'):
                e -= drop
            else:
                s += drop
            drop = 0
        kept.append((chrom, s, e, rest))
    out.write(hdr)
    for chrom, s, e, rest in kept:
        out.write(f'{chrom}\t{s}\t{e}\t{rest}')
hdr = None; ivs = []
for line in open(a.input):
    if line.startswith('+'):
        flush(hdr, ivs); hdr = line; ivs = []
    elif line.strip():
        f = line.split('\t', 3)
        if len(f) >= 4: ivs.append((f[0], int(f[1]), int(f[2]), f[3]))
flush(hdr, ivs)
