#!/usr/bin/env python3
"""F4.1 joint analysis: fit a small logistic model on the NA12878-cdna
overlap zone, evaluate FROZEN on the other five samples' zones, and check
whether the score separates the TKSM presence-suppressed-real set."""
import csv
import math
import sys

import numpy as np

BASE = '/scratch1/rob/long-read-ecosystem/rebench-2026-08-03/p0-validation'
SAMPLES = [
    'nanosim-NA12878-cdna', 'nanosim-H9-cdna', 'nanosim-NA12878-drna',
    'nanosim-H9-drna', 'tksm-RSII', 'tksm-SQ2',
]

tlen = {}
for line in open('/scratch1/rob/long-read-ecosystem/rebench-2026-08-03/junc/junctions.tsv'):
    p = line.split('\t')
    t = p[0]
    i = t.rfind('.')
    if i > 0 and t[i + 1:].isdigit():
        t = t[:i]
    tlen[t] = int(p[1])


def load(sample, lo=0.0, hi=5.0):
    path = {
        'nanosim-NA12878-cdna': f'{BASE}/zonesep_cdna.tsv',
    }.get(sample, f'{BASE}/zonesep_{sample}.tsv')
    X, y, names = [], [], []
    for r in csv.DictReader(open(path), delimiter='\t'):
        est = float(r['est'])
        if not (lo < est <= hi):
            continue
        L = tlen.get(r['tname'], 1000)
        X.append([
            float(r['s3_left']),
            float(r['s3_right']),
            math.log(max(L, 50)),
            float(r['mean_cands']),
            math.log(max(float(r['soft_reads']), 1e-3)),
        ])
        y.append(1.0 if float(r['truth']) > 0 else 0.0)
        names.append(r['tname'])
    return np.array(X), np.array(y), names


def fit_logistic(X, y, iters=3000, lr=0.05):
    mu, sd = X.mean(0), X.std(0) + 1e-9
    Z = (X - mu) / sd
    w = np.zeros(Z.shape[1] + 1)
    A = np.hstack([Z, np.ones((len(Z), 1))])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-A @ w))
        w -= lr * (A.T @ (p - y)) / len(y)
    return mu, sd, w


def score(X, mu, sd, w):
    A = np.hstack([(X - mu) / sd, np.ones((len(X), 1))])
    return A @ w


def auc(s, y):
    order = np.argsort(s)
    r = np.empty(len(s))
    sr = s[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and sr[j + 1] == sr[i]:
            j += 1
        r[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    n1 = y.sum()
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float('nan')
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


Xtr, ytr, _ = load('nanosim-NA12878-cdna')
mu, sd, w = fit_logistic(Xtr, ytr)
feat = ['s3_left', 's3_right', 'log_tlen', 'mean_cands', 'log_soft']
print('model (fit on NA12878-cdna zone, est<=5):',
      {f: round(float(v), 3) for f, v in zip(feat, w[:-1])})
print(f"\n{'sample':<24}{'n':>7}{'%real':>7}{'AUC full':>10}{'AUC est<=1':>12}")
for s in SAMPLES:
    try:
        X, y, _ = load(s)
    except FileNotFoundError:
        print(f'{s:<24}   (missing)')
        continue
    a = auc(score(X, mu, sd, w), y)
    X1, y1, _ = load(s, hi=1.0)
    a1 = auc(score(X1, mu, sd, w), y1) if len(X1) else float('nan')
    tag = ' (train)' if s == 'nanosim-NA12878-cdna' else ''
    print(f'{s:<24}{len(X):>7}{100*y.mean():>6.0f}%{a:>10.3f}{a1:>12.3f}{tag}')

# TKSM rescue: does the frozen score separate presence-suppressed real vs zero?
print('\nTKSM SQ2 presence-suppressed set (q<0.5), frozen score:')


def norm(n):
    n = n.split('|')[0]
    i = n.rfind('.')
    return n[:i] if i > 0 and n[i + 1:].isdigit() else n


truth = {}
f = open('/scratch1/rob/long-read-ecosystem/oarfish-evaluation-data/sim-panel/tksm/ground_truth/SQ2_ground_truth.csv')
rd = csv.reader(f, delimiter='\t')
next(rd)
for r in rd:
    truth[norm(r[0])] = float(r[1])
qsupp = set()
for r in csv.DictReader(open(f'{BASE}/presence/logistic/tksm-SQ2.presence.tsv'), delimiter='\t'):
    if float(r['presence_prob']) < 0.5 and float(r['est_reads']) >= 0 :
        qsupp.add(norm(r['tname']))
X, y, names = load('tksm-SQ2', lo=0.0, hi=1e18)
mask = np.array([n in qsupp for n in names])
if mask.sum():
    Xs, ys = X[mask], np.array([1.0 if truth.get(n, 0) > 0 else 0.0 for n, m in zip(names, mask) if m])
    print(f'  n={mask.sum()} suppressed-with-stats, real={100*ys.mean():.1f}%, AUC={auc(score(Xs, mu, sd, w), ys):.3f}')
else:
    print('  no overlap')
