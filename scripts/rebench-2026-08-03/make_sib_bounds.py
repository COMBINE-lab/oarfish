#!/usr/bin/env python3
"""Build the sibling-terminal boundary table for the F4.1 zone-separation
measurement.

For each transcript t, list the transcript-coordinate offsets of its internal
exon boundaries (same convention as junctions.tsv: cumulative exon length in
transcript orientation, so an offset J means transcript positions [0,J) come
from exons up to the junction) that genomically coincide (within --tol nt)
with ANOTHER same-chrom, same-strand transcript's terminal exon boundary
(transcript start or end). A read whose alignment terminates at such an
offset is better explained by the sibling terminating there.

Also re-derives junctions.tsv from the GTF and verifies the offset convention
against an existing junctions.tsv (--check) before writing anything.
"""
import argparse
import re
import sys
from collections import defaultdict

TID = re.compile(r'transcript_id "([^"]+)"')


def norm(t):
    i = t.rfind('.')
    return t[:i] if i > 0 and t[i + 1:].isdigit() else t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--check', help='existing junctions.tsv to verify convention')
    ap.add_argument('--tol', type=int, default=5)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    # transcript -> (chrom, strand, [exon (start,end) 1-based inclusive, genomic order])
    tx = {}
    for line in open(a.gtf):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        m = TID.search(f[8])
        if not m:
            continue
        t = m.group(1)
        rec = tx.setdefault(t, (f[0], f[6], []))
        rec[2].append((int(f[3]), int(f[4])))
    for t in tx:
        tx[t][2].sort()

    # verify junction-offset convention against the existing table
    if a.check:
        existing = {}
        for line in open(a.check):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3 and p[2]:
                existing[norm(p[0])] = [int(x) for x in p[2].split(',')]
        ok = bad = 0
        for t, (chrom, strand, exons) in tx.items():
            key = norm(t)
            if key not in existing or len(exons) < 2:
                continue
            lens = [e - s + 1 for s, e in exons]
            if strand == '-':
                lens = lens[::-1]
            offs = []
            c = 0
            for l in lens[:-1]:
                c += l
                offs.append(c)
            if offs == existing[key]:
                ok += 1
            else:
                bad += 1
                if bad <= 3:
                    sys.stderr.write(f'MISMATCH {t}: derived {offs[:5]} vs table {existing[key][:5]}\n')
        sys.stderr.write(f'convention check: {ok} match, {bad} mismatch\n')
        if bad > ok * 0.05:
            sys.stderr.write('convention check FAILED; not writing output\n')
            sys.exit(1)

    # terminal boundaries: (chrom, strand) -> sorted genomic positions of
    # transcript termini (both 5' and 3' genomic ends)
    term = defaultdict(list)
    for t, (chrom, strand, exons) in tx.items():
        term[(chrom, strand)].append((exons[0][0], t))
        term[(chrom, strand)].append((exons[-1][1], t))
    for k in term:
        term[k].sort()

    import bisect
    def has_other_terminus(chrom, strand, pos, self_t):
        arr = term[(chrom, strand)]
        i = bisect.bisect_left(arr, (pos - a.tol, ''))
        while i < len(arr) and arr[i][0] <= pos + a.tol:
            if arr[i][1] != self_t:
                return True
            i += 1
        return False

    n_tx = n_with = n_off = 0
    with open(a.out, 'w') as w:
        for t, (chrom, strand, exons) in tx.items():
            if len(exons) < 2:
                continue
            n_tx += 1
            lens = [e - s + 1 for s, e in exons]
            g_ends = [e for _, e in exons[:-1]]     # donor side (genomic)
            g_starts = [s for s, _ in exons[1:]]    # acceptor side (genomic)
            if strand == '-':
                lens = lens[::-1]
                # junction j (transcript order) corresponds to genomic junction
                # n-2-j between exons (reversed order)
                g_ends = g_ends[::-1]
                g_starts = g_starts[::-1]
            offs = []
            c = 0
            for j, l in enumerate(lens[:-1]):
                c += l
                # the two genomic positions flanking this junction
                if strand == '+':
                    gpos = (g_ends[j], g_starts[j])
                else:
                    gpos = (g_starts[j], g_ends[j])
                if any(has_other_terminus(chrom, strand, g, t) for g in gpos):
                    offs.append(c)
            if offs:
                n_with += 1
                n_off += len(offs)
                w.write(f'{t}\t{",".join(map(str, offs))}\n')
    sys.stderr.write(
        f'multi-exon transcripts: {n_tx}; with sibling-terminal internal '
        f'boundaries: {n_with}; boundaries: {n_off}\n'
    )


if __name__ == '__main__':
    main()
