#!/usr/bin/env python3
"""Generate the Leg-1 design: per-sample abundance TSVs + truth tables.

See spec.md (pre-registered). Deterministic given --seed.
"""
import argparse
import gzip
import math
import random
import re
from collections import defaultdict
from pathlib import Path

TID = re.compile(r'transcript_id "([^"]+)"')
GID = re.compile(r'gene_id "([^"]+)"')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base-quant', required=True)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=20260817)
    ap.add_argument('--n-dte', type=int, default=2000)
    ap.add_argument('--n-presence', type=int, default=400)
    ap.add_argument('--n-dtu', type=int, default=400)
    ap.add_argument('--sigma-log', type=float, default=0.2)
    ap.add_argument('--reps', type=int, default=4)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = Path(a.out)
    (out / 'samples').mkdir(parents=True, exist_ok=True)

    # transcript -> gene from the GTF (versioned ids, matching --use-whole-id)
    t2g = {}
    for line in open(a.gtf):
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'transcript':
            continue
        mt, mg = TID.search(f[8]), GID.search(f[8])
        if mt and mg:
            t2g[mt.group(1)] = mg.group(1)

    # base profile restricted to the annotation
    base = {}
    with open(a.base_quant) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2 and p[0] in t2g:
                v = float(p[1])
                if v > 0:
                    base[p[0]] = v
    absent = sorted(t for t in t2g if t not in base)
    print(f'annotation transcripts: {len(t2g)}; base expressed: {len(base)}; absent pool: {len(absent)}')

    taken = set()

    def sample_from(pool, n):
        pool = [t for t in pool if t not in taken]
        chosen = rng.sample(pool, n)
        taken.update(chosen)
        return chosen

    # DTE
    dte_pool = sorted(t for t, v in base.items() if v >= 2)
    dte = sample_from(dte_pool, a.n_dte)
    rng.shuffle(dte)
    half = a.n_dte // 2
    lfc = {}
    mags = [0.5, 1.0, 2.0]
    for i, t in enumerate(dte[:half]):
        lfc[t] = mags[i % 3]
    for i, t in enumerate(dte[half:]):
        lfc[t] = -mags[i % 3]

    # presence-off
    off_pool = sorted(t for t, v in base.items() if 2 <= v <= 50)
    p_off = sample_from(off_pool, a.n_presence)
    # presence-on: absent transcripts, tpm resampled from the off set's values
    p_on = sample_from(absent, a.n_presence)
    on_tpm = {t: base[rng.choice(p_off)] for t in p_on}

    # DTU: genes with >=2 isoforms at tpm>=2, none of whose isoforms is taken
    g2t = defaultdict(list)
    for t, v in base.items():
        if v >= 2:
            g2t[t2g[t]].append(t)
    dtu_pool = sorted(
        g for g, ts in g2t.items()
        if len(ts) >= 2 and not any(t in taken for t in ts)
    )
    dtu_genes = rng.sample(dtu_pool, a.n_dtu)
    dtu_swap = {}
    for g in dtu_genes:
        ts = sorted(g2t[g], key=lambda t: -base[t])
        t1, t2 = ts[0], ts[1]
        dtu_swap[t1], dtu_swap[t2] = t2, t1
        taken.update(ts[:2])

    # condition means
    def cond_tpm(t, cond):
        v = base.get(t, 0.0)
        if cond == 'B':
            if t in lfc:
                v *= 2 ** lfc[t]
            if t in p_off:
                v = 0.0
            if t in on_tpm:
                v = on_tpm[t]
            if t in dtu_swap:
                v = base.get(dtu_swap[t], 0.0)
        return v

    universe = sorted(set(base) | set(p_on))
    with open(out / 'truth_sets.tsv', 'w') as w:
        w.write('transcript\tgene\tset\tlog2fc\ttpm_A\ttpm_B\n')
        for t in universe:
            s = 'null'
            l = 0.0
            if t in lfc:
                s, l = 'dte', lfc[t]
            elif t in p_off:
                s = 'presence_off_B'
            elif t in on_tpm:
                s = 'presence_on_B'
            elif t in dtu_swap:
                s = 'dtu'
            w.write(f'{t}\t{t2g.get(t, ".")}\t{s}\t{l}\t{cond_tpm(t, "A"):.4f}\t{cond_tpm(t, "B"):.4f}\n')

    # per-sample jittered abundances
    manifest = open(out / 'manifest.tsv', 'w')
    manifest.write('sample\tcondition\tseed\tabundance\n')
    for ci, cond in enumerate(('A', 'B')):
        for r in range(1, a.reps + 1):
            sample = f'{cond}{r}'
            sseed = a.seed + ci * a.reps + r
            jrng = random.Random(sseed * 7919)
            path = out / 'samples' / f'{sample}.tsv'
            with open(path, 'w') as w:
                w.write('target_id\ttpm\tcell\n')
                for t in universe:
                    v = cond_tpm(t, cond)
                    if v <= 0:
                        continue
                    v *= math.exp(jrng.gauss(0.0, a.sigma_log))
                    w.write(f'{t}\t{v:.6f}\t.\n')
            manifest.write(f'{sample}\t{cond}\t{sseed}\t{path}\n')
    manifest.close()
    print(f'sets: dte={len(dte)} presence_off={len(p_off)} presence_on={len(p_on)} dtu_genes={len(dtu_genes)}')
    print(f'wrote {2 * a.reps} sample abundance files under {out}/samples')


if __name__ == '__main__':
    main()
