#!/usr/bin/env python3
"""Merge every results TSV into the single machine-readable deliverable."""

import csv
import sys
from pathlib import Path

SRC = [
    ("panelA/results.tsv", "A", "demote"),
    ("panelA_sirv/results.tsv", "A-sirv", "demote"),
    ("panelB/results.tsv", "B", "demote"),
    ("control/results_A.tsv", "A", "main"),
    ("control/results_B.tsv", "B", "main"),
]

FIELDS = [
    "panel", "sample", "arm", "binary", "technology", "truth_type",
    "status", "wall_s", "peak_rss_gb", "exit_status",
    "n_universe", "ref_md5", "truth_orphans", "n_truth", "n_shared",
    "spearman", "kendall", "pearson_log1p", "ccc_log1p", "rmse", "mard",
    "spearman_inner", "mard_inner", "estimate_cv", "truth_is_constant",
    "rescaled",
    "kernel", "rank_blend_active", "learned_censor_scale_nt",
    "coverage_s", "em_s",
]


def main():
    base = Path(__file__).resolve().parent
    out_path = base / "coverage-rebenchmark-2026-08-03.tsv"
    rows = []
    for rel, panel, binary in SRC:
        p = base / rel
        if not p.exists():
            sys.stderr.write("missing (skipped): %s\n" % rel)
            continue
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                r["panel"] = panel
                r["binary"] = binary
                rows.append(r)

    rows.sort(key=lambda r: (r["panel"], r["sample"], r["binary"], r["arm"]))
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("wrote %s (%d rows)" % (out_path, len(rows)))


if __name__ == "__main__":
    main()
