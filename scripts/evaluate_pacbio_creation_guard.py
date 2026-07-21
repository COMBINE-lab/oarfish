#!/usr/bin/env python3
"""Post-hoc grid for a one-sided coverage-driven abundance creation guard."""

import argparse
import csv
import gzip
import json
from pathlib import Path

from evaluate_quant import moments, ranks, read_quant, transcript_key


def illumina_truth(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        return {transcript_key(row["Name"].split("|", 1)[0], True): float(row["NumReads"])
                for row in csv.DictReader(handle, delimiter="\t")}


def count_truth(path):
    result = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if fields:
                result[transcript_key(fields[0], True)] = float(fields[1])
    return result


def metrics(x, y):
    import math
    if sum(y):
        scale = sum(x) / sum(y)
        y = [value * scale for value in y]
    pearson, ccc = moments(x, y)
    spearman, _ = moments(ranks(x), ranks(y))
    rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(x, y)) / len(x))
    mard = sum(abs(a - b) / (abs(a) + abs(b)) if a or b else 0.0
               for a, b in zip(x, y)) / len(x)
    return dict(pearson=pearson, spearman=spearman, ccc=ccc, rmse=rmse, mard=mard)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--truth-type", choices=("illumina", "counts"), required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--physical", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    truth = illumina_truth(args.truth) if args.truth_type == "illumina" else count_truth(args.truth)
    baseline = read_quant(args.baseline, True, "|")
    physical = read_quant(args.physical, True, "|")
    names = sorted(set(truth) & set(baseline) & set(physical))
    x = [truth[name] for name in names]
    b = [baseline[name] for name in names]
    p = [physical[name] for name in names]
    total = sum(b)
    rows = []
    for fold in (1, 2, 4, 9, 19, 99):
        for additive_cpm in (0.1, 0.5, 1, 2, 5, 10):
            additive = additive_cpm * total / 1_000_000
            guarded = [min(after, before * (fold + 1) + additive)
                       if after > before else after for before, after in zip(b, p)]
            rows.append(dict(sample=args.sample, fold=fold, additive_cpm=additive_cpm,
                             **metrics(x, guarded)))
    rows.extend((dict(sample=args.sample, fold=label, additive_cpm="-", **metrics(x, values))
                 for label, values in (("baseline", b), ("physical", p))))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
