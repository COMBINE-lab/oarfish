#!/usr/bin/env python3
"""Evaluate E0-to-E2 SIRV fold-change recovery at matched read depths."""

import argparse
import csv
import math
from pathlib import Path

from evaluate_quant import moments, ranks, read_quant, transcript_key


def concentrations(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            transcript_key(row["transcript"], True): (float(row["E0"]), float(row["E2"]))
            for row in csv.DictReader(handle, delimiter="\t")
        }


def normalized_quant(path, names):
    values = read_quant(path, True, "|")
    result = {name: values.get(name, 0.0) for name in names}
    total = sum(result.values())
    return {name: value / total for name, value in result.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("concentrations", type=Path)
    parser.add_argument("quant_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--full-quant-dir", type=Path)
    args = parser.parse_args()
    truth = concentrations(args.concentrations)
    fields = ("depth", "model", "transcripts", "fc_pearson", "fc_spearman", "fc_rmse", "direction_accuracy")
    with args.output.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fields, delimiter="\t")
        writer.writeheader()
        for depth in (10000, 25000, 50000):
            for model in ("none", "auto", "abundance-blend"):
                e0 = normalized_quant(args.quant_dir / f"sirv-e0-drna-{depth // 1000}k-{model}.quant", truth)
                e2 = normalized_quant(args.quant_dir / f"sirv-e2-drna-{depth // 1000}k-{model}.quant", truth)
                floor = min(value for values in (e0, e2) for value in values.values() if value > 0) / 2
                expected = [math.log2(e2_t / e0_t) for e0_t, e2_t in truth.values()]
                observed = [math.log2((e2[name] + floor) / (e0[name] + floor)) for name in truth]
                pearson, _ = moments(expected, observed)
                spearman, _ = moments(ranks(expected), ranks(observed))
                rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(expected, observed)) / len(expected))
                direction = sum((a == 0 and abs(b) < 0.5) or (a * b > 0) for a, b in zip(expected, observed)) / len(expected)
                writer.writerow(dict(zip(fields, (depth, model, len(truth), pearson, spearman, rmse, direction))))
        if args.full_quant_dir:
            for model in ("none", "auto", "abundance-blend"):
                e0 = normalized_quant(args.full_quant_dir / f"sirv-e0-drna-full-{model}.quant", truth)
                e2 = normalized_quant(args.full_quant_dir / f"sirv-e2-drna-full-{model}.quant", truth)
                floor = min(value for values in (e0, e2) for value in values.values() if value > 0) / 2
                expected = [math.log2(e2_t / e0_t) for e0_t, e2_t in truth.values()]
                observed = [math.log2((e2[name] + floor) / (e0[name] + floor)) for name in truth]
                pearson, _ = moments(expected, observed)
                spearman, _ = moments(ranks(expected), ranks(observed))
                rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(expected, observed)) / len(expected))
                direction = sum((a == 0 and abs(b) < 0.5) or (a * b > 0) for a, b in zip(expected, observed)) / len(expected)
                writer.writerow(dict(zip(fields, ("full", model, len(truth), pearson, spearman, rmse, direction))))


if __name__ == "__main__":
    main()
