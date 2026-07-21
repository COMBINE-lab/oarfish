#!/usr/bin/env python3
"""Run coverage-model variants on the public Kinnex Figure 2 BAM subset."""

import argparse
import concurrent.futures
import csv
import json
import re
import subprocess
from pathlib import Path


MODELS = {
    "none": ("none", "full"),
    "logistic": ("logistic", "full"),
    "auto": ("auto", "full"),
    "abundance-blend": ("auto", "abundance-blend"),
    "pacbio-physical-endpoint": ("auto", "pacbio-physical-endpoint"),
}


def timed_value(text, label):
    match = re.search(rf"{re.escape(label)}: (.+)", text)
    return match.group(1).strip() if match else ""


def run_one(binary, bam, sample, model, output_dir, threads, em_accel, logistic_weight,
            endpoint_weight):
    coverage_model, ablation = MODELS[model]
    prefix = output_dir / f"{sample}-{model}"
    log_path = prefix.with_suffix(".log")
    time_path = prefix.with_suffix(".time.txt")
    command = [
        "/usr/bin/time", "-v", "-o", str(time_path), str(binary),
        "--threads", str(threads), "--filter-group", "no-filters",
        "--alignments", str(bam), "--output", str(prefix),
        "--coverage-model", coverage_model, "--coverage-ablation", ablation,
        "--logistic-weight", str(logistic_weight),
        "--endpoint-weight", str(endpoint_weight),
        "--seq-tech", "pac-bio", "--em-accel", em_accel,
    ]
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    timing = time_path.read_text(encoding="utf-8")
    metadata = json.loads(prefix.with_suffix(".meta_info.json").read_text())
    return {
        "sample": sample,
        "model": model,
        "em_accel": em_accel,
        "logistic_weight": logistic_weight,
        "endpoint_weight": endpoint_weight,
        "bam": str(bam),
        "wall_time": timed_value(timing, "Elapsed (wall clock) time (h:mm:ss or m:ss)"),
        "peak_rss_kb": timed_value(timing, "Maximum resident set size (kbytes)"),
        "alignment_seconds": metadata["alignment_time"]["seconds"],
        "coverage_seconds": metadata["coverage_model_time"]["seconds"],
        "em_seconds": metadata["em_time"]["seconds"],
        "em_evaluations": metadata["em_evaluations"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--binary", type=Path, default=Path("target/release/oarfish"))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS))
    parser.add_argument("--em-accel", choices=("none", "squarem", "daarem"), default="none")
    parser.add_argument("--logistic-weight", type=float, default=1.0)
    parser.add_argument("--endpoint-weight", type=float, default=0.5)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    for day in ("day0", "day5"):
        for replicate in range(1, 4):
            sample = f"{day}-rep{replicate}"
            bam = args.data_dir / sample / f"{sample}.aligned.bam"
            for model in args.models:
                jobs.append((args.binary.resolve(), bam, sample, model,
                             args.output_dir.resolve(), args.threads, args.em_accel,
                             args.logistic_weight, args.endpoint_weight))
    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(run_one, *job) for job in jobs]
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            rows.append(row)
            print(f"completed {row['sample']} {row['model']}", flush=True)
    fields = list(rows[0])
    with (args.output_dir / "run_metrics.tsv").open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: (row["sample"], row["model"])))


if __name__ == "__main__":
    main()
