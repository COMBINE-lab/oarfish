#!/usr/bin/env python3
"""Diagnose paired calibration changes from the PacBio physical endpoint model."""

import argparse
import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np


def tx_id(name):
    value = name.split("|", 1)[0]
    return value.rsplit(".", 1)[0]


def gene_id(name):
    fields = name.split("|")
    value = fields[1] if len(fields) > 1 and fields[1] else fields[0]
    return value.rsplit(".", 1)[0]


def read_truth(path):
    result = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[tx_id(row["Name"])] = float(row["NumReads"])
    return result


def read_quant(path, ambiguity_path):
    values, annotations = {}, {}
    with path.open(encoding="utf-8", newline="") as quant_handle, \
            ambiguity_path.open(encoding="utf-8", newline="") as ambiguity_handle:
        quant_rows = csv.DictReader(quant_handle, delimiter="\t")
        ambiguity_rows = csv.DictReader(ambiguity_handle, delimiter="\t")
        for quant, ambiguity in zip(quant_rows, ambiguity_rows):
            key = tx_id(quant["tname"])
            values[key] = float(quant["num_reads"])
            unique = float(ambiguity["unique_reads"])
            ambig = float(ambiguity["ambig_reads"])
            annotations[key] = {
                "gene": gene_id(quant["tname"]),
                "length": float(quant["len"]),
                "unique": unique,
                "ambiguous": ambig,
            }
    return values, annotations


def rankdata(values):
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    boundaries = np.r_[0, np.flatnonzero(np.diff(sorted_values)) + 1, len(values)]
    ranks = np.empty(len(values), dtype=float)
    for start, end in zip(boundaries[:-1], boundaries[1:]):
        ranks[order[start:end]] = (start + end + 1) / 2
    return ranks


def metrics(x, y):
    mx, my = x.mean(), y.mean()
    vx, vy = np.mean((x - mx) ** 2), np.mean((y - my) ** 2)
    covariance = np.mean((x - mx) * (y - my))
    pearson = covariance / math.sqrt(vx * vy) if vx and vy else 0.0
    ccc = 2 * covariance / (vx + vy + (mx - my) ** 2)
    rmse = math.sqrt(np.mean((x - y) ** 2))
    denominator = np.abs(x) + np.abs(y)
    mard = np.divide(np.abs(x - y), denominator, out=np.zeros_like(x), where=denominator != 0).mean()
    rx, ry = rankdata(x), rankdata(y)
    spearman = np.corrcoef(rx, ry)[0, 1]
    return {"pearson": pearson, "spearman": spearman, "ccc": ccc,
            "rmse": rmse, "mard": mard}


def sufficient_by_gene(x, y, rx, ry, groups, group_count):
    columns = np.column_stack((
        np.ones(len(x)), x, y, x * x, y * y, x * y,
        (x - y) ** 2,
        np.divide(np.abs(x - y), np.abs(x) + np.abs(y),
                  out=np.zeros_like(x), where=(np.abs(x) + np.abs(y)) != 0),
        rx, ry, rx * rx, ry * ry, rx * ry,
    ))
    return np.vstack([np.bincount(groups, weights=columns[:, i], minlength=group_count)
                      for i in range(columns.shape[1])]).T


def metrics_from_sums(sums):
    n, sx, sy, sxx, syy, sxy, sse, smard, srx, sry, srxx, sryy, srxy = sums
    mx, my, mrx, mry = sx / n, sy / n, srx / n, sry / n
    vx, vy = sxx / n - mx * mx, syy / n - my * my
    cov = sxy / n - mx * my
    rvx, rvy = srxx / n - mrx * mrx, sryy / n - mry * mry
    rcov = srxy / n - mrx * mry
    return np.array((
        cov / math.sqrt(max(vx * vy, 1e-300)),
        rcov / math.sqrt(max(rvx * rvy, 1e-300)),
        2 * cov / (vx + vy + (mx - my) ** 2),
        math.sqrt(sse / n), smard / n,
    ))


def cluster_bootstrap(x, old, new, genes, replicates, seed):
    _, group = np.unique(genes, return_inverse=True)
    group_count = group.max() + 1
    rx = rankdata(x)
    old_sums = sufficient_by_gene(x, old, rx, rankdata(old), group, group_count)
    new_sums = sufficient_by_gene(x, new, rx, rankdata(new), group, group_count)
    rng = np.random.RandomState(seed)
    differences = np.empty((replicates, 5))
    for iteration in range(replicates):
        weights = rng.multinomial(group_count, np.full(group_count, 1 / group_count))
        differences[iteration] = (metrics_from_sums(weights @ new_sums)
                                  - metrics_from_sums(weights @ old_sums))
    names = ("pearson", "spearman", "ccc", "rmse", "mard")
    return {name: {
        "median_delta": float(np.median(differences[:, index])),
        "ci95": [float(v) for v in np.percentile(differences[:, index], (2.5, 97.5))],
        "probability_improved": float(np.mean(differences[:, index] > 0)
                                      if name not in ("rmse", "mard")
                                      else np.mean(differences[:, index] < 0)),
    } for index, name in enumerate(names)}


def aggregate_genes(x, old, new, genes):
    _, group = np.unique(genes, return_inverse=True)
    count = group.max() + 1
    return (np.bincount(group, weights=x, minlength=count),
            np.bincount(group, weights=old, minlength=count),
            np.bincount(group, weights=new, minlength=count))


def subset_summary(x, old, new, masks):
    result = {}
    for label, mask in masks.items():
        if mask.sum() < 2:
            continue
        before, after = metrics(x[mask], old[mask]), metrics(x[mask], new[mask])
        result[label] = {key: {"old": before[key], "new": after[key],
                               "delta": after[key] - before[key]}
                         for key in before}
        result[label]["transcripts"] = int(mask.sum())
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--physical", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260721)
    args = parser.parse_args()

    truth = read_truth(args.truth)
    old_values, annotations = read_quant(args.baseline, args.baseline.with_suffix(".ambig_info.tsv"))
    new_values, _ = read_quant(args.physical, args.physical.with_suffix(".ambig_info.tsv"))
    names = sorted(set(truth) & set(old_values) & set(new_values))
    x = np.array([truth[name] for name in names])
    old = np.array([old_values[name] for name in names])
    new = np.array([new_values[name] for name in names])
    old *= x.sum() / old.sum()
    new *= x.sum() / new.sum()
    genes = np.array([annotations[name]["gene"] for name in names])
    lengths = np.array([annotations[name]["length"] for name in names])
    unique = np.array([annotations[name]["unique"] for name in names])
    ambiguous = np.array([annotations[name]["ambiguous"] for name in names])
    support = unique + ambiguous
    ambiguity = np.divide(ambiguous, support, out=np.zeros_like(support), where=support > 0)

    positive = x > 0
    abundance_cuts = np.percentile(x[positive], (100 / 3, 200 / 3))
    masks = {
        "truth_zero": ~positive,
        "abundance_low": positive & (x <= abundance_cuts[0]),
        "abundance_mid": positive & (x > abundance_cuts[0]) & (x <= abundance_cuts[1]),
        "abundance_high": x > abundance_cuts[1],
        "length_lt_1kb": lengths < 1000,
        "length_1_2kb": (lengths >= 1000) & (lengths < 2000),
        "length_ge_2kb": lengths >= 2000,
        "no_read_support": support == 0,
        "unique_only": (support > 0) & (ambiguity == 0),
        "partly_ambiguous": (ambiguity > 0) & (ambiguity < 1),
        "ambiguous_only": ambiguity == 1,
    }
    squared_delta = (x - new) ** 2 - (x - old) ** 2
    order = np.argsort(squared_delta)[::-1]
    total_delta = squared_delta.sum()
    top = []
    for index in order[:100]:
        top.append({"transcript": names[index], "gene": genes[index],
                    "truth": x[index], "baseline": old[index], "physical": new[index],
                    "squared_error_delta": squared_delta[index],
                    "length": lengths[index], "ambiguity_fraction": ambiguity[index]})
    contribution = {}
    for count in (1, 5, 10, 100):
        contribution[str(count)] = float(squared_delta[order[:count]].sum() / total_delta) \
            if total_delta else math.nan

    sensitivity = {}
    exclusions = {
        "exclude_largest_transcript": np.arange(len(x)) != order[0],
        "exclude_largest_gene": genes != genes[order[0]],
    }
    for label, mask in exclusions.items():
        before, after = metrics(x[mask], old[mask]), metrics(x[mask], new[mask])
        sensitivity[label] = {
            key: {"baseline": before[key], "physical": after[key],
                  "delta": after[key] - before[key]}
            for key in before
        }

    gx, gold, gnew = aggregate_genes(x, old, new, genes)
    result = {
        "sample": args.sample, "transcripts": len(names), "genes": len(gx),
        "transcript_metrics": {"baseline": metrics(x, old), "physical": metrics(x, new)},
        "gene_metrics": {"baseline": metrics(gx, gold), "physical": metrics(gx, gnew)},
        "gene_cluster_bootstrap": cluster_bootstrap(
            x, old, new, genes, args.bootstrap_replicates, args.seed),
        "strata": subset_summary(x, old, new, masks),
        "largest_rmse_regression_contributors": top,
        "top_squared_error_contribution_fraction": contribution,
        "influence_sensitivity": sensitivity,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
