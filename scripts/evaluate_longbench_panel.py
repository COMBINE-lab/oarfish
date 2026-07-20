#!/usr/bin/env python3
"""Evaluate the fixed LongBench panel against matched Illumina estimates."""

import argparse
import csv
import gzip
import json
import math
from pathlib import Path

from evaluate_quant import moments, ranks, read_quant, transcript_key


CELLS = ("H69", "H146", "H211", "H526", "H1975", "H2228", "HCC827", "SHP77")
TECHNOLOGIES = ("cdna", "drna", "pb")
MODELS = ("none", "adaptive", "auto")


def read_illumina(path):
    result = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = transcript_key(row["Name"].split("|", 1)[0], True)
            result[key] = result.get(key, 0.0) + float(row["NumReads"])
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("quant_dir", type=Path)
    parser.add_argument("matrix_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument(
        "--available-only",
        action="store_true",
        help="skip sample/technology combinations absent from a confirmation tier",
    )
    args = parser.parse_args()
    fields = [
        "cell", "technology", "model", "transcripts", "pearson", "spearman",
        "ccc", "rmse", "mard", "coverage_seconds", "em_seconds",
        "kernel", "intact_fraction", "degradation_hazard_per_kb",
        "degradation_correction_weight", "mean_degradation_strength",
        "technical_truncation_fraction", "degraded_fraction", "mean_reliability",
        "degradation_shape", "held_out_log_predictive_density", "hazard_fold_sd",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for cell in CELLS:
            truth = read_illumina(args.matrix_dir / f"{cell}.tsv.gz")
            for technology in TECHNOLOGIES:
                for model in MODELS:
                    prefix = args.quant_dir / f"{cell}-{technology}-{model}"
                    if args.available_only and not prefix.with_suffix(".quant").exists():
                        continue
                    estimate = read_quant(prefix.with_suffix(".quant"), True, "|")
                    names = sorted(set(truth) & set(estimate))
                    x = [truth[name] for name in names]
                    y = [estimate[name] for name in names]
                    if sum(y) > 0:
                        scale = sum(x) / sum(y)
                        y = [value * scale for value in y]
                    pearson, ccc = moments(x, y)
                    spearman, _ = moments(ranks(x), ranks(y))
                    rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(x, y)) / len(x))
                    mard = sum(abs(a - b) / (abs(a) + abs(b)) if a or b else 0.0
                               for a, b in zip(x, y)) / len(x)
                    with prefix.with_suffix(".meta_info.json").open(encoding="utf-8") as handle:
                        metadata = json.load(handle)
                    diagnostics = metadata.get("coverage_diagnostics", {})
                    degradation = diagnostics.get("degradation", {})
                    adaptive = diagnostics.get("adaptive", diagnostics)
                    writer.writerow({
                        "cell": cell, "technology": technology, "model": model,
                        "transcripts": len(names), "pearson": pearson,
                        "spearman": spearman, "ccc": ccc, "rmse": rmse,
                        "mard": mard,
                        "coverage_seconds": metadata["coverage_model_time"]["seconds"],
                        "em_seconds": metadata["em_time"]["seconds"],
                        "kernel": diagnostics.get("technology_kernel", ""),
                        "intact_fraction": degradation.get("intact_fraction", ""),
                        "degradation_hazard_per_kb": degradation.get(
                            "degradation_hazard_per_kb", ""),
                        "degradation_correction_weight": degradation.get(
                            "degradation_correction_weight", ""),
                        "mean_degradation_strength": degradation.get(
                            "mean_degradation_strength", ""),
                        "technical_truncation_fraction": degradation.get(
                            "technical_truncation_fraction", ""),
                        "degraded_fraction": degradation.get("degraded_fraction", ""),
                        "degradation_shape": degradation.get("degradation_shape", ""),
                        "held_out_log_predictive_density": degradation.get(
                            "held_out_log_predictive_density", ""),
                        "hazard_fold_sd": degradation.get("hazard_fold_sd", ""),
                        "mean_reliability": adaptive.get("mean_reliability", ""),
                    })


if __name__ == "__main__":
    main()
