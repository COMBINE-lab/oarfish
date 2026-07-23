#!/usr/bin/env python3
"""Run reproducible coverage ablations with wall time, peak RSS, and accuracy."""

import argparse
import csv
import gzip
import json
import math
import re
import subprocess
from pathlib import Path

from evaluate_quant import moments, ranks, read_quant, read_truth, transcript_key

ABLATIONS = (
    "none", "full", "no-logistic", "no-endpoint", "no-support-gate",
    "no-agreement-gate", "no-bayes-cap", "quality-gate",
    "uncertainty-gate", "eq-class-training", "all-candidates",
    "abundance-blend", "pacbio-physical-endpoint",
)


def illumina_truth(path):
    result = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = transcript_key(row["Name"].split("|", 1)[0], True)
            result[key] = result.get(key, 0.0) + float(row["NumReads"])
    return result


def accuracy(truth, quant):
    estimate = read_quant(quant, True, "|")
    names = sorted(set(truth) & set(estimate))
    x, y = [truth[n] for n in names], [estimate[n] for n in names]
    if sum(y) > 0:
        scale = sum(x) / sum(y)
        y = [v * scale for v in y]
    pearson, ccc = moments(x, y)
    spearman, _ = moments(ranks(x), ranks(y))
    rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(x, y)) / len(x))
    mard = sum(abs(a - b) / (abs(a) + abs(b)) if a or b else 0.0
               for a, b in zip(x, y)) / len(x)
    return len(names), pearson, spearman, ccc, rmse, mard


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path,
                        help="TSV: sample, bam, technology, truth, truth_type")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--binary", type=Path, default=Path("target/release/oarfish"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--ablations", nargs="+", choices=ABLATIONS,
                        default=list(ABLATIONS))
    parser.add_argument("--samples", nargs="+",
                        help="run only these manifest sample names")
    parser.add_argument("--result-name", default="results.tsv")
    parser.add_argument("--endpoint-weight", type=float, default=0.5)
    parser.add_argument("--score-prob-denom", type=float,
                        help="override alignment score-to-probability temperature")
    parser.add_argument(
        "--alignment-calibration", choices=("none", "agreement", "auto"), default="auto"
    )
    parser.add_argument("--candidate-pruning", choices=("none", "dominance"), default="none")
    parser.add_argument("--dominance-bayes-factor", type=float, default=4.0)
    parser.add_argument("--censoring-model", choices=("none", "adaptive"), default="none")
    parser.add_argument("--rank-blend", choices=("none", "fixed", "auto"), default="auto")
    parser.add_argument("--rank-blend-floor", type=float, default=0.8)
    parser.add_argument("--coverage-abundance-midpoint-per-million", type=float, default=300.0)
    parser.add_argument("--em-accel", choices=("none", "squarem", "daarem"),
                        default="squarem")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(args.manifest.open(), delimiter="\t"))
    if args.samples:
        rows = [row for row in rows if row["sample"] in set(args.samples)]
    result_path = args.output_dir / args.result_name
    fields = ("sample", "technology", "ablation", "repeat", "transcripts",
              "pearson", "spearman", "ccc", "rmse", "mard", "wall_seconds",
              "peak_rss_kb", "alignment_seconds", "coverage_seconds", "em_seconds",
              "mean_reliability", "mean_quality_gate", "mean_uncertainty_gate",
              "ambiguous_training_reads", "capped_reads")
    with result_path.open("w", newline="", encoding="utf-8") as destination:
        writer = csv.DictWriter(destination, fields, delimiter="\t")
        writer.writeheader()
        for item in rows:
            truth = (illumina_truth(Path(item["truth"])) if item["truth_type"] == "illumina"
                     else read_truth(Path(item["truth"]), True, None))
            for ablation in args.ablations:
                for repeat in range(1, args.repeats + 1):
                    prefix = args.output_dir / f'{item["sample"]}-{ablation}-r{repeat}'
                    time_file = prefix.with_suffix(".time.txt")
                    model = "none" if ablation == "none" else "auto"
                    ablation_arg = "full" if ablation == "none" else ablation
                    command = ["/usr/bin/time", "-v", "-o", str(time_file),
                               str(args.binary), "--alignments", item["bam"],
                               "--output", str(prefix), "--coverage-model", model,
                               "--seq-tech", item["technology"], "--coverage-ablation",
                               ablation_arg, "--filter-group", "no-filters", "--threads",
                               str(args.threads), "--em-accel", args.em_accel,
                               "--endpoint-weight", str(args.endpoint_weight)]
                    if args.score_prob_denom is not None:
                        command.extend(["--score-prob-denom", str(args.score_prob_denom)])
                    command.extend(["--candidate-pruning", args.candidate_pruning,
                                    "--dominance-bayes-factor", str(args.dominance_bayes_factor),
                                    "--censoring-model", args.censoring_model,
                                    "--rank-blend", args.rank_blend,
                                    "--rank-blend-floor", str(args.rank_blend_floor),
                                    "--coverage-abundance-midpoint-per-million",
                                    str(args.coverage_abundance_midpoint_per_million)])
                    command.extend(["--alignment-calibration", args.alignment_calibration])
                    with prefix.with_suffix(".log").open("w") as log:
                        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
                    timing = time_file.read_text()
                    wall = re.search(r"Elapsed \(wall clock\) time.*: ([0-9:.]+)", timing).group(1)
                    parts = [float(v) for v in wall.split(":")]
                    wall_seconds = sum(v * 60 ** i for i, v in enumerate(reversed(parts)))
                    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", timing).group(1))
                    metadata = json.loads(prefix.with_suffix(".meta_info.json").read_text())
                    diag = metadata.get("coverage_diagnostics", {})
                    adaptive = diag.get("adaptive", diag)
                    metrics = accuracy(truth, prefix.with_suffix(".quant"))
                    writer.writerow(dict(zip(fields, (
                        item["sample"], item["technology"], ablation, repeat, *metrics,
                        wall_seconds, rss, metadata["alignment_time"]["seconds"],
                        metadata["coverage_model_time"]["seconds"],
                        metadata["em_time"]["seconds"], adaptive.get("mean_reliability", ""),
                        adaptive.get("mean_quality_gate", ""),
                        adaptive.get("mean_uncertainty_gate", ""),
                        adaptive.get("ambiguous_training_reads", ""),
                        adaptive.get("capped_reads", "")))))
                    destination.flush()


if __name__ == "__main__":
    main()
