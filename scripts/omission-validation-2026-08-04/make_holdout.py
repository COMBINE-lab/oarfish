#!/usr/bin/env python3
"""Build a TranSigner-style holdout annotation.

Transcripts are deleted from the set the reads were actually simulated *from*
(the expressed set), not from the catalogue at large: deleting unexpressed
transcripts orphans no reads and measures nothing.

Emits the holdout GTF plus a manifest of which transcripts/genes were removed,
so detection precision can be scored against it.
"""

import argparse
import random
import re
import sys

TID = re.compile(r'transcript_id "([^"]+)"')
GID = re.compile(r'gene_id "([^"]+)"')


def norm(t):
    i = t.rfind('.')
    return t[:i] if i > 0 and t[i + 1:].isdigit() else t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--truth', required=True, help='transcript<TAB>count')
    ap.add_argument('--fraction', type=float, required=True)
    ap.add_argument('--seed', type=int, default=20260804)
    ap.add_argument('--out-gtf', required=True)
    ap.add_argument('--out-manifest', required=True)
    a = ap.parse_args()

    expressed = set()
    for line in open(a.truth):
        f = line.split()
        if len(f) >= 2:
            try:
                if float(f[1]) > 0:
                    expressed.add(norm(f[0]))
            except ValueError:
                continue
    sys.stderr.write('expressed transcripts in truth: %d\n' % len(expressed))

    # transcript -> gene, restricted to what the annotation actually contains
    t2g = {}
    for line in open(a.gtf):
        if line.startswith('#'):
            continue
        mt, mg = TID.search(line), GID.search(line)
        if mt and mg:
            t2g.setdefault(norm(mt.group(1)), mg.group(1))
    sys.stderr.write('transcripts in annotation: %d\n' % len(t2g))

    candidates = sorted(expressed & set(t2g))
    sys.stderr.write('expressed AND annotated (deletable): %d\n' % len(candidates))
    rng = random.Random(a.seed)
    n = int(round(len(candidates) * a.fraction))
    deleted = set(rng.sample(candidates, n))
    sys.stderr.write('deleting %d (%.0f%%)\n' % (len(deleted), 100 * a.fraction))

    kept = wrote = 0
    with open(a.out_gtf, 'w') as out:
        for line in open(a.gtf):
            if line.startswith('#'):
                out.write(line)
                continue
            mt = TID.search(line)
            if mt and norm(mt.group(1)) in deleted:
                kept += 1
                continue
            out.write(line)
            wrote += 1
    sys.stderr.write('gtf lines written %d, dropped %d\n' % (wrote, kept))

    # genes that lost at least one expressed isoform
    lost_genes = {t2g[t] for t in deleted}
    with open(a.out_manifest, 'w') as m:
        m.write('#deleted_transcripts\t%d\n' % len(deleted))
        m.write('#genes_losing_expressed_isoform\t%d\n' % len(lost_genes))
        for t in sorted(deleted):
            m.write('T\t%s\t%s\n' % (t, t2g[t]))
        for g in sorted(lost_genes):
            m.write('G\t%s\n' % g)
    sys.stderr.write('genes losing an expressed isoform: %d\n' % len(lost_genes))


if __name__ == '__main__':
    main()
