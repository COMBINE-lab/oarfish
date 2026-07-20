#!/usr/bin/env python3
"""Reproducible truth-versus-oarfish accuracy summary.

Truth is a two-column, headerless transcript/count file.  Quantification is an
oarfish TSV containing `tname` and `num_reads`.  The JSON schema is intentionally
stable so staged coverage reports can concatenate results without scraping logs.
"""

import argparse
import csv
import json
import math


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + j - 1) / 2 + 1
        for k in range(i, j):
            out[order[k]] = rank
        i = j
    return out


def moments(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    vx = sum((v - mx) ** 2 for v in x) / len(x)
    vy = sum((v - my) ** 2 for v in y) / len(y)
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y)) / len(x)
    pearson = cov / math.sqrt(vx * vy) if vx > 0 and vy > 0 else 0.0
    ccc_denom = vx + vy + (mx - my) ** 2
    ccc = 2 * cov / ccc_denom if ccc_denom > 0 else 0.0
    return pearson, ccc


def transcript_key(name, strip_version, delimiter=None):
    if delimiter:
        name = name.split(delimiter, 1)[0]
    return name.rsplit(".", 1)[0] if strip_version and name.rsplit(".", 1)[-1].isdigit() else name


def read_truth(path, strip_version, delimiter):
    result = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if fields:
                key = transcript_key(fields[0], strip_version, delimiter)
                result[key] = result.get(key, 0.0) + float(fields[1])
    return result


def read_quant(path, strip_version, delimiter):
    with open(path, encoding="utf-8", newline="") as handle:
        result = {}
        for row in csv.DictReader(handle, delimiter="\t"):
            key = transcript_key(row["tname"], strip_version, delimiter)
            result[key] = result.get(key, 0.0) + float(row["num_reads"])
        return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("truth")
    parser.add_argument("quant")
    parser.add_argument("--label", default="")
    parser.add_argument("--strip-version", action="store_true",
                        help="match accessions after removing a terminal numeric version")
    parser.add_argument("--id-delimiter",
                        help="compare only the portion of each identifier before this delimiter")
    parser.add_argument("--raw-counts", action="store_true",
                        help="do not rescale estimates to the truth library size")
    parser.add_argument("--intersection", action="store_true",
                        help="compare only transcript identifiers present in both inputs")
    args = parser.parse_args()
    truth = read_truth(args.truth, args.strip_version, args.id_delimiter)
    estimate = read_quant(args.quant, args.strip_version, args.id_delimiter)
    names = sorted(set(truth) & set(estimate) if args.intersection
                   else set(truth) | set(estimate))
    x = [truth.get(n, 0.0) for n in names]
    y = [estimate.get(n, 0.0) for n in names]
    original_estimated_total = sum(y)
    if not args.raw_counts and original_estimated_total > 0:
        scale = sum(x) / original_estimated_total
        y = [value * scale for value in y]
    pearson, ccc = moments(x, y)
    spearman, _ = moments(ranks(x), ranks(y))
    rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(x, y)) / len(x))
    mard = sum(abs(a - b) / (abs(a) + abs(b)) if a or b else 0.0 for a, b in zip(x, y)) / len(x)
    print(json.dumps({
        "label": args.label, "transcripts": len(names), "truth_total": sum(x),
        "estimated_total": original_estimated_total, "comparison_scale": sum(x) / original_estimated_total,
        "pearson": pearson, "spearman": spearman,
        "ccc": ccc, "rmse": rmse, "mard": mard,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
