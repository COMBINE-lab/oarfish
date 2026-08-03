#!/usr/bin/env python3
"""Phase 0b: is 'fraction of alignments terminating at an internal junction'
a per-transcript feature that separates truly-expressed from truly-absent
transcripts? If so it is a candidate feature for the presence/absence
component, which carries the +0.0261 detection headroom."""

import bisect
import gzip
import pickle
import struct
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, '.')
from metrics import read_truth

BAM, TRUTH, CAP = sys.argv[1], sys.argv[2], int(sys.argv[3])
W = 5  # tolerance, nt


def strip(n):
    i = n.rfind('.')
    return n[:i] if i > 0 and n[i + 1:].isdigit() else n


junc = pickle.load(open('junc/junctions.pkl', 'rb'))
J = {strip(k): v for k, v in junc.items()}
truth = read_truth(TRUTH, 'counts')

f = gzip.open(BAM, 'rb')
assert f.read(4) == b'BAM\x01'
lt, = struct.unpack('<i', f.read(4))
f.read(lt)
nref, = struct.unpack('<i', f.read(4))
refs = []
for _ in range(nref):
    ln, = struct.unpack('<i', f.read(4))
    nm = f.read(ln)[:-1].decode()
    f.read(4)
    refs.append(strip(nm))

CONSUME = {0, 2, 3, 7, 8}
n_aln = defaultdict(int)
n_hit = defaultdict(int)
n = 0
while n < CAP:
    b = f.read(4)
    if len(b) < 4:
        break
    bs, = struct.unpack('<i', b)
    rec = f.read(bs)
    if len(rec) < bs:
        break
    refid, pos = struct.unpack('<ii', rec[0:8])
    lrn = rec[8]
    ncig, = struct.unpack('<H', rec[12:14])
    flag, = struct.unpack('<H', rec[14:16])
    n += 1
    if refid < 0 or flag & 0x4:
        continue
    t = refs[refid]
    if t not in J:
        continue
    L, js = J[t]
    if not js:
        continue
    o = 32 + lrn
    span = 0
    for i in range(ncig):
        v, = struct.unpack('<I', rec[o + 4 * i:o + 4 * i + 4])
        if (v & 0xf) in CONSUME:
            span += v >> 4
    a, e = pos, pos + span
    n_aln[t] += 1
    ok = False
    for x in (a, e):
        i = bisect.bisect_left(js, x)
        for c in (i - 1, i):
            if 0 <= c < len(js) and abs(js[c] - x) <= W:
                ok = True
                break
        if ok:
            break
    if ok:
        n_hit[t] += 1

print('records scanned: %d ; transcripts with >=1 alignment: %d' % (n, len(n_aln)))
for MIN in (10, 50):
    ex, ab = [], []
    for t, c in n_aln.items():
        if c < MIN:
            continue
        fr = n_hit[t] / c
        (ex if truth.get(t, 0.0) > 0 else ab).append(fr)
    ex, ab = np.array(ex), np.array(ab)
    print()
    print('--- transcripts with >=%d alignments ---' % MIN)
    print('  truly EXPRESSED: n=%-7d mean junc-terminated frac = %.4f  median %.4f'
          % (len(ex), ex.mean(), np.median(ex)))
    print('  truly ABSENT   : n=%-7d mean junc-terminated frac = %.4f  median %.4f'
          % (len(ab), ab.mean(), np.median(ab)))
    if ex.mean() > 0:
        print('  ratio of means : %.2fx' % (ab.mean() / ex.mean()))
    for thr in (0.2, 0.4, 0.6):
        pe, pa = (ex > thr).mean(), (ab > thr).mean()
        print('   frac>%.1f : expressed %.3f  absent %.3f   enrichment %s'
              % (thr, pe, pa, ('%.2fx' % (pa / pe)) if pe > 0 else 'inf'))
