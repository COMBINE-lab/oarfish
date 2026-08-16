#!/usr/bin/env python3
"""Re-measure the presence/absence headroom (perfect detection / perfect
abundance oracles) on converged-EM quants, per Panel B sample.

Usage: headroom.py <manifest.tsv> <quant_dir> (expects <quant_dir>/<sample>.quant)
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from metrics import evaluate, read_quant, read_truth

manifest, qdir = Path(sys.argv[1]), Path(sys.argv[2])
rows = list(csv.DictReader(open(manifest), delimiter="\t"))
print("sample\tbaseline\tperfect_detection\tperfect_abundance\tfp_txps\tfp_mass_pct")
agg = {"b": [], "d": [], "a": []}
for s in rows:
    truth = read_truth(s["truth"], s["truth_type"])
    qp = qdir / (s["sample"] + ".quant")
    if not qp.exists():
        continue
    est = read_quant(qp)
    base = evaluate(truth, est, rescale=False)["spearman"]
    # perfect detection: truth-zeros forced to zero
    det = {k: (v if truth.get(k, 0.0) > 0 else 0.0) for k, v in est.items()}
    d = evaluate(truth, det, rescale=False)["spearman"]
    # perfect abundance: expressed set to truth (zeros keep estimate)
    ab = dict(est)
    for k, v in truth.items():
        if v > 0:
            ab[k] = v
    a = evaluate(truth, ab, rescale=False)["spearman"]
    fp = [(k, v) for k, v in est.items() if v > 0 and truth.get(k, 0.0) == 0]
    tot = sum(est.values()) or 1.0
    print(
        f"{s['sample']}\t{base:.4f}\t{d:.4f}\t{a:.4f}\t{len(fp)}\t"
        f"{100.0 * sum(v for _, v in fp) / tot:.3f}"
    )
    agg["b"].append(base)
    agg["d"].append(d)
    agg["a"].append(a)
n = len(agg["b"]) or 1
print(
    f"MEAN\t{sum(agg['b'])/n:.4f}\t{sum(agg['d'])/n:.4f}\t{sum(agg['a'])/n:.4f}\t\t"
)
