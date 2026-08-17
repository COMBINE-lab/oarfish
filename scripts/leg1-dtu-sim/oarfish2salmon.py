#!/usr/bin/env python3
"""Convert an oarfish output (<prefix>.quant + <prefix>.infreps.pq) into a
salmon-format quantification directory readable by tximport
(dropInfReps=FALSE), fishpond/swish, and edgeR::catchSalmon.

Layout written:
  <out>/quant.sf
  <out>/cmd_info.json
  <out>/aux_info/meta_info.json
  <out>/aux_info/bootstrap/bootstraps.gz   (num_boot x num_txp float64, LE)
  <out>/aux_info/bootstrap/names.tsv.gz

EffectiveLength is set equal to Length (documented shortcut: oarfish does not
compute effective lengths; length-derived biases were measured
rank-preserving on this data).
"""
import argparse
import gzip
import json
import struct
from pathlib import Path

import pyarrow.parquet as pq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--prefix', required=True, help='oarfish output prefix')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    out = Path(a.out)
    (out / 'aux_info' / 'bootstrap').mkdir(parents=True, exist_ok=True)

    names, lens, counts = [], [], []
    with open(f'{a.prefix}.quant') as fh:
        header = fh.readline()
        assert header.startswith('tname')
        for line in fh:
            p = line.rstrip('\n').split('\t')
            names.append(p[0])
            lens.append(int(p[1]))
            counts.append(float(p[2]))
    n = len(names)
    denom = sum(c / l for c, l in zip(counts, lens) if l > 0) or 1.0
    with open(out / 'quant.sf', 'w') as w:
        w.write('Name\tLength\tEffectiveLength\tTPM\tNumReads\n')
        for nm, l, c in zip(names, lens, counts):
            tpm = (c / l) / denom * 1e6 if l > 0 else 0.0
            w.write(f'{nm}\t{l}\t{l:.3f}\t{tpm:.6f}\t{c:.3f}\n')

    t = pq.read_table(f'{a.prefix}.infreps.pq')
    boot_cols = [c for c in t.column_names if c.startswith('bootstrap.')]
    boot_cols.sort(key=lambda c: int(c.split('.')[1]))
    nb = len(boot_cols)
    with gzip.open(out / 'aux_info' / 'bootstrap' / 'bootstraps.gz', 'wb') as w:
        for c in boot_cols:
            v = t.column(c).to_pylist()
            # bootstrap vectors may carry novel states past the annotated
            # transcripts; truncate to the quant.sf universe
            w.write(struct.pack(f'<{n}d', *v[:n]))
    with gzip.open(out / 'aux_info' / 'bootstrap' / 'names.tsv.gz', 'wt') as w:
        w.write('\t'.join(names))
    meta = {
        'salmon_version': '1.10.0',
        'samp_type': 'bootstrap',
        'num_bootstraps': nb,
        'num_targets': n,
        'num_valid_targets': n,
        'serialized_eq_classes': False,
        'note': 'converted from oarfish output by oarfish2salmon.py',
    }
    json.dump(meta, open(out / 'aux_info' / 'meta_info.json', 'w'), indent=1)
    json.dump({'salmon_program': 'oarfish-convert', 'auxDir': 'aux_info'},
              open(out / 'cmd_info.json', 'w'), indent=1)
    print(f'{a.prefix} -> {out}: {n} targets, {nb} bootstraps')


if __name__ == '__main__':
    main()
