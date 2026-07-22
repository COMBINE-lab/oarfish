#!/usr/bin/env python3
"""Compare no coverage, auto, and abundance blending on truth-bearing samples."""

import argparse
import csv
import json
import math
import re
import subprocess
from pathlib import Path

from evaluate_quant import moments, ranks, read_quant, read_truth, transcript_key


MODELS = ("none", "auto", "abundance-blend")


def read_concentrations(path, mix):
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            transcript_key(row["transcript"], True): float(row[mix])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def summarize(truth, quant_path):
    estimate = read_quant(quant_path, True, "|")
    names = sorted(truth)
    x = [truth[name] for name in names]
    y = [estimate.get(name, 0.0) for name in names]
    if sum(y):
        scale = sum(x) / sum(y)
        y = [value * scale for value in y]
    pearson, ccc = moments(x, y)
    spearman, _ = moments(ranks(x), ranks(y))
    rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(x, y)) / len(x))
    mard = sum(
        abs(a - b) / (abs(a) + abs(b)) if a or b else 0.0
        for a, b in zip(x, y)
    ) / len(x)
    mean = sum(y) / len(y)
    estimate_cv = (
        math.sqrt(sum((value - mean) ** 2 for value in y) / len(y)) / mean
        if mean else math.nan
    )
    ordered = sorted(x)
    cut1, cut2 = ordered[len(ordered) // 3], ordered[2 * len(ordered) // 3]
    strata = (x_i <= cut1 and "low" or x_i <= cut2 and "mid" or "high" for x_i in x)
    errors = {"low": [], "mid": [], "high": []}
    for a, b, stratum in zip(x, y, strata):
        errors[stratum].append(abs(a - b) / (abs(a) + abs(b)) if a or b else 0.0)
    stratum_mard = {
        key: sum(values) / len(values) if values else math.nan
        for key, values in errors.items()
    }
    return (len(names), pearson, spearman, ccc, rmse, mard, estimate_cv,
            stratum_mard["low"], stratum_mard["mid"], stratum_mard["high"])


def wall_seconds(text):
    match = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", text)
    parts = [float(value) for value in match.group(1).split(":")]
    return sum(value * 60 ** power for power, value in enumerate(reversed(parts)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--binary", type=Path, default=Path("target/release/oarfish"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS))
    parser.add_argument("--samples", nargs="+")
    parser.add_argument("--candidate-pruning", choices=("none", "dominance"), default="none")
    parser.add_argument("--dominance-bayes-factor", type=float, default=2.0)
    parser.add_argument("--censoring-model", choices=("none", "adaptive"), default="none")
    parser.add_argument("--rank-blend", choices=("none", "fixed", "auto"), default="auto")
    parser.add_argument("--rank-blend-floor", type=float, default=0.8)
    parser.add_argument("--coverage-abundance-midpoint-per-million", type=float, default=300.0)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with args.manifest.open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    if args.samples:
        selected = set(args.samples)
        samples = [sample for sample in samples if sample["sample"] in selected]

    fields = (
        "sample", "technology", "mix", "depth", "model", "transcripts",
        "pearson", "spearman", "ccc", "rmse", "mard", "estimate_cv",
        "low_mard", "mid_mard", "high_mard", "wall_seconds", "peak_rss_kb",
        "alignment_seconds", "coverage_seconds", "em_seconds", "em_evaluations",
    )
    with (args.output_dir / "results.tsv").open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fields, delimiter="\t")
        writer.writeheader()
        for sample in samples:
            truth = (
                read_truth(Path(sample["truth"]), True, None)
                if sample.get("truth_type") == "counts"
                else read_concentrations(Path(sample["truth"]), sample["mix"])
            )
            for model in args.models:
                prefix = args.output_dir / f'{sample["sample"]}-{model}'
                time_path = prefix.with_suffix(".time.txt")
                coverage_model = "none" if model == "none" else "auto"
                ablation = "abundance-blend" if model == "abundance-blend" else "full"
                command = [
                    "/usr/bin/time", "-v", "-o", str(time_path), str(args.binary),
                    "--reads", sample["reads"], "--annotated", sample["reference"],
                    "--output", str(prefix), "--seq-tech", sample["technology"],
                    "--coverage-model", coverage_model, "--coverage-ablation", ablation,
                    "--filter-group", "no-filters", "--threads", str(args.threads),
                    "--candidate-pruning", args.candidate_pruning,
                    "--dominance-bayes-factor", str(args.dominance_bayes_factor),
                    "--censoring-model", args.censoring_model,
                    "--rank-blend", args.rank_blend,
                    "--rank-blend-floor", str(args.rank_blend_floor),
                    "--coverage-abundance-midpoint-per-million",
                    str(args.coverage_abundance_midpoint_per_million),
                ]
                with prefix.with_suffix(".log").open("w", encoding="utf-8") as log:
                    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
                timing = time_path.read_text(encoding="utf-8")
                rss = int(re.search(
                    r"Maximum resident set size \(kbytes\): (\d+)", timing
                ).group(1))
                metadata = json.loads(prefix.with_suffix(".meta_info.json").read_text())
                metrics = summarize(truth, prefix.with_suffix(".quant"))
                writer.writerow(dict(zip(fields, (
                    sample["sample"], sample["technology"], sample["mix"], sample["depth"],
                    model, *metrics, wall_seconds(timing), rss,
                    metadata["alignment_time"]["seconds"],
                    metadata["coverage_model_time"]["seconds"],
                    metadata["em_time"]["seconds"], metadata["em_evaluations"],
                ))))
                out.flush()


if __name__ == "__main__":
    main()
