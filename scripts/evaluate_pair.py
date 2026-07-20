#!/usr/bin/env python3
"""Compare two oarfish estimates when transcript-level truth is unavailable."""

import argparse
import csv
import json
import math


def read_quant(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return {row["tname"]: float(row["num_reads"])
                for row in csv.DictReader(handle, delimiter="\t")}


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + j - 1) / 2 + 1
        for k in range(i, j):
            result[order[k]] = rank
        i = j
    return result


def correlation(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    vx = sum((v - mx) ** 2 for v in x)
    vy = sum((v - my) ** 2 for v in y)
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return cov / math.sqrt(vx * vy) if vx and vy else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("first")
    parser.add_argument("second")
    parser.add_argument("--label", default="")
    parser.add_argument("--min-total-proportion", type=float, default=1e-6,
                        help="retain transcripts whose two proportions sum to at least this value")
    args = parser.parse_args()
    first, second = read_quant(args.first), read_quant(args.second)
    first_total, second_total = sum(first.values()), sum(second.values())
    names = sorted(set(first) | set(second))
    pairs = [(first.get(name, 0.0) / first_total,
              second.get(name, 0.0) / second_total) for name in names]
    pairs = [(a, b) for a, b in pairs if a + b >= args.min_total_proportion]
    x, y = zip(*pairs)
    midpoint = [(a + b) / 2 for a, b in pairs]
    jsd = 0.5 * sum(a * math.log(a / m) if a else 0.0
                    for a, m in zip(x, midpoint))
    jsd += 0.5 * sum(b * math.log(b / m) if b else 0.0
                     for b, m in zip(y, midpoint))
    print(json.dumps({
        "label": args.label,
        "transcripts": len(pairs),
        "pearson": correlation(x, y),
        "spearman": correlation(ranks(x), ranks(y)),
        "l1_distance": sum(abs(a - b) for a, b in pairs),
        "jensen_shannon": jsd,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
