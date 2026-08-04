#!/usr/bin/env python3
"""Score annotation-omission detection against the holdout manifest.

A flagged locus is a true positive if any of its member transcripts belongs to a
gene that lost an expressed isoform in the holdout -- the definition used by
docs/annotation-omission-evaluation-2026-07-25.md.
"""

import argparse
import csv
import re
import sys

TID = re.compile(r'transcript_id "([^"]+)"')
GID = re.compile(r'gene_id "([^"]+)"')


def norm(t):
    i = t.rfind('.')
    return t[:i] if i > 0 and t[i + 1:].isdigit() else t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--unexplained', required=True)
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--gtf', required=True, help='FULL annotation, for transcript->gene')
    a = ap.parse_args()

    lost_genes, deleted = set(), set()
    for line in open(a.manifest):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if f[0] == 'G':
            lost_genes.add(f[1])
        elif f[0] == 'T':
            deleted.add(f[1])

    t2g = {}
    for line in open(a.gtf):
        if line.startswith('#'):
            continue
        mt, mg = TID.search(line), GID.search(line)
        if mt and mg:
            t2g.setdefault(norm(mt.group(1)), mg.group(1))

    rows = list(csv.DictReader(open(a.unexplained), delimiter='\t'))
    sys.stderr.write('flagged loci: %d ; genes that lost an expressed isoform: %d\n'
                     % (len(rows), len(lost_genes)))

    scored = []
    for r in rows:
        txps = [norm(x) for x in r['transcripts'].split(',') if x and x != '...']
        genes = {t2g[t] for t in txps if t in t2g}
        tp = bool(genes & lost_genes)
        scored.append((float(r.get('flagged_reads', 0) or 0), tp, genes))

    print('%-20s %-8s %-8s %-8s %-10s' % ('min flagged reads', 'loci', 'TP', 'FP', 'precision'))
    for thr in (1, 2, 3, 5, 10, 50):
        sel = [s for s in scored if s[0] >= thr]
        if not sel:
            print('%-20d %-8d %-8s %-8s %-10s' % (thr, 0, '-', '-', '-'))
            continue
        tp = sum(1 for s in sel if s[1])
        print('%-20d %-8d %-8d %-8d %-10.3f' % (thr, len(sel), tp, len(sel) - tp, tp / len(sel)))

    found = set()
    for s in scored:
        found |= (s[2] & lost_genes)
    print()
    print('genes correctly identified as incomplete: %d / %d  (recall %.4f)'
          % (len(found), len(lost_genes), len(found) / max(len(lost_genes), 1)))


if __name__ == '__main__':
    main()
