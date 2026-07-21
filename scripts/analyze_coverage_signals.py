#!/usr/bin/env python3
"""Relate observable Oarfish signals to useful Kinnex coverage corrections.

The analysis deliberately treats the truth only as an evaluation target.  All
predictors (length, inferred abundance, unique support, ambiguous support, and
the proposed coverage correction) are available from an ordinary Oarfish run.
"""

import argparse
import csv
import math
import statistics
from pathlib import Path

from evaluate_kinnex_fig2 import lengths_from_gtf, quant_counts, xlsx_truth


def mean(values):
    values = list(values)
    return sum(values) / len(values)


def rank(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        value = (i + j - 1) / 2 + 1
        for k in order[i:j]:
            out[k] = value
        i = j
    return out


def pearson(x, y):
    if len(x) < 3:
        return float("nan")
    mx, my = mean(x), mean(y)
    xx = sum((v - mx) ** 2 for v in x)
    yy = sum((v - my) ** 2 for v in y)
    return (sum((a - mx) * (b - my) for a, b in zip(x, y))
            / math.sqrt(xx * yy)) if xx and yy else float("nan")


def spearman(x, y):
    return pearson(rank(x), rank(y))


def zscore(values):
    center = mean(values)
    sd = statistics.pstdev(values)
    return [(v - center) / sd for v in values]


def log_cpm(counts, names):
    total = sum(counts.get(name, 0.0) for name in names)
    return {name: math.log2((counts.get(name, 0.0) + 2.0)
                            / (total + 4.0) * 1_000_000.0)
            for name in names}


def support(path, transcript_names):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != len(transcript_names):
        raise ValueError(f"{path}: support/quant row mismatch")
    result = {}
    for name, row in zip(transcript_names, rows):
        unique = float(row["unique_reads"])
        ambiguous = float(row["ambig_reads"])
        total = unique + ambiguous
        result[name] = (unique, ambiguous, total)
    return result


def quant_and_support(directory, sample, model):
    quant_path = directory / f"{sample}-{model}.quant"
    counts = quant_counts(quant_path)
    names = []
    with quant_path.open(encoding="utf-8", newline="") as handle:
        names = [row["tname"] for row in csv.DictReader(handle, delimiter="\t")]
    info = support(directory / f"{sample}-{model}.ambig_info.tsv", names)
    return counts, info


def summarize(rows, panel):
    print(f"\n{panel}: {len(rows)} transcript/replicate observations")
    print("signal\tSpearman(signal, desired correction)\t"
          "Spearman(signal, realized error improvement)")
    signals = ("length", "log_support", "unique_fraction", "log_abundance",
               "abs_proposed_correction", "proposed_correction")
    for signal in signals:
        x = [row[signal] for row in rows]
        desired = [row["desired_correction"] for row in rows]
        gain = [row["error_improvement"] for row in rows]
        print(f"{signal}\t{spearman(x, desired):.5f}\t{spearman(x, gain):.5f}")

    print("\nlength_bin\tn\tmean_none_abs_error\tmean_auto_abs_error\t"
          "auto_win_fraction\tcorrection_alignment")
    cuts = ((0, 750), (750, 1250), (1250, 2000), (2000, float("inf")))
    for lo, hi in cuts:
        subset = [row for row in rows if lo <= row["length"] < hi]
        if not subset:
            continue
        align = spearman([row["proposed_correction"] for row in subset],
                         [row["desired_correction"] for row in subset])
        print(f"[{lo},{hi})\t{len(subset)}\t"
              f"{mean(row['none_abs_error'] for row in subset):.5f}\t"
              f"{mean(row['auto_abs_error'] for row in subset):.5f}\t"
              f"{mean(row['error_improvement'] > 0 for row in subset):.5f}\t"
              f"{align:.5f}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("quant_dir", type=Path)
    parser.add_argument("truth_xlsx", type=Path)
    parser.add_argument("gtf", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    truth = xlsx_truth(args.truth_xlsx)
    lengths = lengths_from_gtf(args.gtf)
    names = sorted(set(truth) & set(lengths))
    cached = {}
    for day in ("day0", "day5"):
        for rep in range(1, 4):
            sample = f"{day}-rep{rep}"
            for model in ("none", "auto"):
                cached[sample, model] = quant_and_support(args.quant_dir, sample, model)

    rows = []
    # Panel D: corrections to directly comparable log2 fold changes.
    for rep in range(1, 4):
        d0, d5 = f"day0-rep{rep}", f"day5-rep{rep}"
        none0 = log_cpm(cached[d0, "none"][0], names)
        none5 = log_cpm(cached[d5, "none"][0], names)
        auto0 = log_cpm(cached[d0, "auto"][0], names)
        auto5 = log_cpm(cached[d5, "auto"][0], names)
        for name in names:
            observed_none = none0[name] - none5[name]
            observed_auto = auto0[name] - auto5[name]
            expected = math.log2(truth[name]["e1"] / truth[name]["e2"])
            unique = sum(cached[s, "none"][1][name][0] for s in (d0, d5))
            ambiguous = sum(cached[s, "none"][1][name][1] for s in (d0, d5))
            total = unique + ambiguous
            rows.append({"panel": "D", "replicate": rep, "transcript": name,
                         "length": lengths[name], "log_support": math.log1p(total),
                         "unique_fraction": unique / total if total else 0.0,
                         "log_abundance": mean((none0[name], none5[name])),
                         "proposed_correction": observed_auto - observed_none,
                         "abs_proposed_correction": abs(observed_auto - observed_none),
                         "desired_correction": expected - observed_none,
                         "none_abs_error": abs(expected - observed_none),
                         "auto_abs_error": abs(expected - observed_auto),
                         "error_improvement": abs(expected - observed_none) - abs(expected - observed_auto)})

    # Panel E: standardize within replicate because the paper compares
    # concentration and log-CPM by correlation, not on a common absolute scale.
    expected_e = zscore([math.log(truth[name]["mw"] * truth[name]["e1"])
                         for name in names])
    for rep in range(1, 4):
        sample = f"day0-rep{rep}"
        none = log_cpm(cached[sample, "none"][0], names)
        auto = log_cpm(cached[sample, "auto"][0], names)
        none_z = zscore([none[name] for name in names])
        auto_z = zscore([auto[name] for name in names])
        for i, name in enumerate(names):
            unique, ambiguous, total = cached[sample, "none"][1][name]
            rows.append({"panel": "E", "replicate": rep, "transcript": name,
                         "length": lengths[name], "log_support": math.log1p(total),
                         "unique_fraction": unique / total if total else 0.0,
                         "log_abundance": none[name],
                         "proposed_correction": auto_z[i] - none_z[i],
                         "abs_proposed_correction": abs(auto_z[i] - none_z[i]),
                         "desired_correction": expected_e[i] - none_z[i],
                         "none_abs_error": abs(expected_e[i] - none_z[i]),
                         "auto_abs_error": abs(expected_e[i] - auto_z[i]),
                         "error_improvement": abs(expected_e[i] - none_z[i]) - abs(expected_e[i] - auto_z[i])})

    for panel in ("D", "E"):
        summarize([row for row in rows if row["panel"] == panel], panel)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, rows[0], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
