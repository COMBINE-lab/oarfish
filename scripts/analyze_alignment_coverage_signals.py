#!/usr/bin/env python3
"""Analyze candidate-relative coverage evidence exported by Oarfish."""

import argparse
import csv
import itertools
import math
import statistics
import zlib
from collections import defaultdict
from pathlib import Path


def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def entropy(xs):
    total = sum(xs)
    return -sum((x / total) * math.log(x / total) for x in xs if x > 0) if total else 0.0


def js(left, right):
    sl, sr = sum(left), sum(right)
    if not sl or not sr:
        return float("nan")
    l, r = [x / sl for x in left], [x / sr for x in right]
    middle = [(x + y) / 2 for x, y in zip(l, r)]
    kl = lambda p: sum(x * math.log(x / m) for x, m in zip(p, middle) if x > 0 and m > 0)
    return 0.5 * (kl(l) + kl(r))


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(values)
    for rank, index in enumerate(order):
        out[index] = rank
    return out


def spearman(x, y):
    if len(x) < 3:
        return float("nan")
    x, y = ranks(x), ranks(y)
    mx, my = mean(x), mean(y)
    den = math.sqrt(sum((v - mx) ** 2 for v in x) * sum((v - my) ** 2 for v in y))
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / den if den else float("nan")


def arrays(row):
    split = lambda key, cast: [cast(x) for x in row[key].split(",")]
    return (row["transcripts"].split(","), split("lengths", int),
            split("starts", int), split("ends", int),
            row["strands"].split(","), split("score_probs", float),
            split("coverage_probs", float))


def length_class(length):
    return "lt750" if length < 750 else "750_1249" if length < 1250 else "1250_1999" if length < 2000 else "ge2000"


def analyze(path):
    stats = defaultdict(list)
    hist = {(lc, half): [1.0] * 400 for lc in ("lt750", "750_1249", "1250_1999", "ge2000") for half in (0, 1)}
    candidate_features = defaultdict(list)
    rows = ambiguous = high_conf = 0
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows += 1
            names, lengths, starts, ends, strands, score, cov = arrays(row)
            half = zlib.crc32(row["read"].encode("utf-8")) & 1
            for length, start, end, strand, cp in zip(lengths, starts, ends, strands, cov):
                left = max(0.0, (start - 1) / length)
                right = max(0.0, (length - end) / length)
                five, three = (right, left) if strand == "Reverse" else (left, right)
                aligned = max(0.0, (end - start + 1) / length)
                candidate_features["coverage"].append(cp)
                candidate_features["length"].append(length)
                candidate_features["five_gap"].append(five)
                candidate_features["three_gap"].append(three)
                candidate_features["aligned_fraction"].append(aligned)
                if len(names) == 1:
                    x, y = min(19, int(five * 20)), min(19, int(three * 20))
                    hist[length_class(lengths[0]), half][x * 20 + y] += 1
            if len(names) < 2:
                continue
            ambiguous += 1
            score_total, cov_total = sum(score), sum(cov)
            score_n = [x / score_total for x in score]
            cov_n = [x / cov_total for x in cov]
            combined = [x * y for x, y in zip(score_n, cov_n)]
            combined_total = sum(combined)
            combined = [x / combined_total for x in combined]
            sw, cw, bw = max(range(len(score)), key=score_n.__getitem__), max(range(len(cov)), key=cov_n.__getitem__), max(range(len(combined)), key=combined.__getitem__)
            confidence = score_n[sw]
            if confidence >= 0.9:
                high_conf += 1
                stats["coverage_agrees_high_conf"].append(cw == sw)
                stats["combined_agrees_high_conf"].append(bw == sw)
                stats["coverage_score_winner_logodds"].append(math.log(max(cov_n[sw], 1e-300)) - mean(math.log(max(cov_n[i], 1e-300)) for i in range(len(cov_n)) if i != sw))
            stats["score_entropy"].append(entropy(score_n))
            stats["coverage_entropy"].append(entropy(cov_n))
            stats["combined_entropy"].append(entropy(combined))
            stats["coverage_log_range"].append(math.log(max(cov_n) / max(min(cov_n), 1e-300)))
            stats["winner_agreement"].append(sw == cw)
    summary = {
        "sample": path.name[:-len(".coverage_signals.tsv")],
        "sampled_reads": rows, "ambiguous_reads": ambiguous,
        "high_confidence_reads": high_conf,
        "coverage_winner_agreement": mean(stats["winner_agreement"]),
        "coverage_high_conf_agreement": mean(stats["coverage_agrees_high_conf"]),
        "combined_high_conf_agreement": mean(stats["combined_agrees_high_conf"]),
        "mean_coverage_score_winner_logodds": mean(stats["coverage_score_winner_logodds"]),
        "mean_score_entropy": mean(stats["score_entropy"]),
        "mean_coverage_entropy": mean(stats["coverage_entropy"]),
        "mean_combined_entropy": mean(stats["combined_entropy"]),
        "mean_coverage_log_range": mean(stats["coverage_log_range"]),
        "split_half_js_lt750": js(hist["lt750", 0], hist["lt750", 1]),
        "split_half_js_750_1249": js(hist["750_1249", 0], hist["750_1249", 1]),
        "split_half_js_1250_1999": js(hist["1250_1999", 0], hist["1250_1999", 1]),
        "split_half_js_ge2000": js(hist["ge2000", 0], hist["ge2000", 1]),
    }
    for feature in ("length", "five_gap", "three_gap", "aligned_fraction"):
        summary[f"coverage_rho_{feature}"] = spearman(candidate_features[feature], candidate_features["coverage"])
    pooled = {lc: [a + b for a, b in zip(hist[lc, 0], hist[lc, 1])] for lc in ("lt750", "750_1249", "1250_1999", "ge2000")}
    return summary, pooled


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tables", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    results = [analyze(path) for path in args.tables]
    rows = [result[0] for result in results]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, rows[0], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    print("metric\tmean\tmin\tmax")
    for key in rows[0]:
        if key == "sample": continue
        values = [float(row[key]) for row in rows]
        print(f"{key}\t{mean(values):.6g}\t{min(values):.6g}\t{max(values):.6g}")
    print("\nlength_class\tmean_between_sample_js\tmax_between_sample_js")
    for lc in ("lt750", "750_1249", "1250_1999", "ge2000"):
        values = [js(left[1][lc], right[1][lc]) for left, right in itertools.combinations(results, 2)]
        if values:
            print(f"{lc}\t{mean(values):.6g}\t{max(values):.6g}")


if __name__ == "__main__":
    main()
