#!/usr/bin/env python3
"""Score DE results against the Leg-1 truth sets: per-set power and empirical
FDR at a q/FDR threshold.

Truth semantics at transcript level: dte and presence_on/off transcripts are
truly differential; each dtu gene's two swapped isoforms are truly
differential (opposite directions); everything else expressed in either
condition is null (modulo replicate noise, which is shared).
"""
import argparse
import csv
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--truth', required=True, help='design/truth_sets.tsv')
    ap.add_argument('--de', required=True, help='DE table')
    ap.add_argument('--qcol', default='qvalue')
    ap.add_argument('--tcol', default='transcript')
    ap.add_argument('--thresh', type=float, default=0.05)
    ap.add_argument('--label', default='')
    a = ap.parse_args()

    truth = {}
    tpm = {}
    for r in csv.DictReader(open(a.truth), delimiter='\t'):
        truth[r['transcript']] = r['set']
        tpm[r['transcript']] = (float(r['tpm_A']), float(r['tpm_B']))

    def strip(n):
        n = n.split('|')[0]
        return n

    called = set()
    tested = set()
    for r in csv.DictReader(open(a.de), delimiter='\t'):
        t = strip(r[a.tcol])
        tested.add(t)
        try:
            q = float(r[a.qcol])
        except ValueError:
            continue
        if q <= a.thresh:
            called.add(t)

    sets = defaultdict(lambda: [0, 0])  # set -> [tested, called]
    for t in tested:
        s = truth.get(t, 'null')
        if s == 'dtu':
            # only the two swapped isoforms differ; truth file marks exactly those
            pass
        sets[s][0] += 1
        if t in called:
            sets[s][1] += 1

    true_sets = ('dte', 'presence_off_B', 'presence_on_B', 'dtu')
    n_true_called = sum(sets[s][1] for s in true_sets)
    n_null_called = sets['null'][1]
    total_called = n_true_called + n_null_called
    print(f'== {a.label or a.de} (thresh {a.thresh}) ==')
    print(f'{"set":<16}{"tested":>8}{"called":>8}{"power":>8}')
    for s in true_sets:
        te, ca = sets[s]
        print(f'{s:<16}{te:>8}{ca:>8}{ca / max(te, 1):>8.3f}')
    te, ca = sets['null']
    print(f'{"null":<16}{te:>8}{ca:>8}{"-":>8}')
    print(f'empirical FDR: {n_null_called}/{total_called} = '
          f'{n_null_called / max(total_called, 1):.4f}')


if __name__ == '__main__':
    main()
