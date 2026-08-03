#!/usr/bin/env python3
"""Run the 9-arm coverage re-benchmark (2026-08-03) over a panel.

Arms are defined by docs/coverage-auto-rebenchmark-plan-2026-08-03.md.
Emits one TSV row per (arm, sample) with accuracy, wall time and peak RSS.
"""

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from metrics import (evaluate, read_concentrations, read_quant, read_truth,
                     shared_universe)

# name -> extra oarfish flags. `--filter-group no-filters` is added to every run.
ARMS = {
    "none": [],
    "logistic": ["--model-coverage"],
    "adaptive-bare": [
        "--coverage-model", "auto",
        "--alignment-calibration", "none",
        "--candidate-pruning", "none",
        "--censoring-model", "none",
        "--rank-blend", "none",
    ],
    "auto-new": ["--coverage-model", "auto"],
    "auto-old": [
        "--coverage-model", "auto",
        "--rank-blend", "auto",
        "--candidate-pruning", "auto",
    ],
    "+calib": [
        "--coverage-model", "auto",
        "--alignment-calibration", "agreement",
        "--candidate-pruning", "none",
        "--censoring-model", "none",
        "--rank-blend", "none",
    ],
    "+censor": [
        "--coverage-model", "auto",
        "--alignment-calibration", "none",
        "--candidate-pruning", "none",
        "--censoring-model", "adaptive",
        "--rank-blend", "none",
    ],
    "+prune": [
        "--coverage-model", "auto",
        "--alignment-calibration", "none",
        "--candidate-pruning", "auto",
        "--censoring-model", "none",
        "--rank-blend", "none",
    ],
    "+rank": [
        "--coverage-model", "auto",
        "--alignment-calibration", "none",
        "--candidate-pruning", "none",
        "--censoring-model", "none",
        "--rank-blend", "auto",
    ],
}

FIELDS = [
    "panel", "sample", "arm", "binary", "technology", "truth_type",
    "status", "wall_s", "peak_rss_gb", "exit_status",
    "n_union", "n_truth", "n_est", "n_shared",
    "spearman", "kendall", "pearson_log1p", "ccc_log1p", "rmse", "mard",
    "spearman_inner", "mard_inner", "rescaled",
    "estimate_cv", "truth_is_constant",
    "n_universe", "ref_md5", "truth_orphans",
    "kernel", "rank_blend_active", "learned_censor_scale_nt",
    "coverage_s", "em_s", "quant", "cmd",
]

TIME_WALL = re.compile(r"Elapsed \(wall clock\) time.*?:\s*([0-9:.]+)")
TIME_RSS = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
TIME_EXIT = re.compile(r"Exit status:\s*(\d+)")


def parse_wall(text):
    parts = text.split(":")
    seconds = float(parts[-1])
    if len(parts) > 1:
        seconds += 60 * float(parts[-2])
    if len(parts) > 2:
        seconds += 3600 * float(parts[-3])
    return seconds


def run_one(binary, sample, arm, outdir, threads, extra_flags, timing_path):
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / sample["sample"]
    # --seq-tech is required even in alignment mode: bulk.rs gates the PacBio
    # physical-endpoint kernel on it, and the original selection harness
    # (scripts/run_coverage_ablation.py) passed it too.
    if sample.get("bam"):
        source = ["-a", sample["bam"]]
    else:
        # Read mode, as the original public-manifest harness
        # (scripts/run_coverage_default_benchmark.py) ran these samples.
        # A prebuilt --index is preferred so the reference is byte-identical
        # across arms and is not rebuilt once per run.
        source = ["--reads", sample["reads"]]
        if sample.get("index"):
            source += ["--index", sample["index"]]
        else:
            source += ["--annotated", sample["reference"]]

    cmd = [
        "/usr/bin/time", "-v", "-o", str(timing_path),
        str(binary),
    ] + source + [
        "--seq-tech", sample["technology"],
        "--filter-group", "no-filters",
        "-j", str(threads),
        "-o", str(prefix),
    ] + ARMS[arm] + extra_flags

    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL,
                          stderr=subprocess.PIPE, universal_newlines=True)
    return cmd, proc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("manifest", type=Path,
                    help="TSV with sample,bam,technology,truth,truth_type[,extra]")
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--binary-label", default="demote")
    ap.add_argument("--panel", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--arms", nargs="+", default=list(ARMS))
    ap.add_argument("--samples", nargs="+")
    ap.add_argument("--rescale", choices=["yes", "no"], required=True,
                    help="scale estimate totals to truth (yes for proxy truth)")
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--keep-quant", action="store_true")
    args = ap.parse_args()

    with open(args.manifest) as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))
    if args.samples:
        samples = [s for s in samples if s["sample"] in args.samples]

    args.results.parent.mkdir(parents=True, exist_ok=True)
    new_file = not args.results.exists()
    out = open(args.results, "a", newline="")
    writer = csv.DictWriter(out, fieldnames=FIELDS, delimiter="\t",
                            extrasaction="ignore")
    if new_file:
        writer.writeheader()
        out.flush()

    truth_cache = {}
    rescale = args.rescale == "yes"

    for sample in samples:
        # SIRV rows carry no truth_type but do carry a mixture column; the
        # original harness keyed off exactly this.
        ttype = sample.get("truth_type") or ""
        if not ttype and sample.get("mix"):
            ttype = "concentrations"
        sample["truth_type"] = ttype

        key = (sample["truth"], ttype, sample.get("mix", ""))
        if key not in truth_cache:
            if ttype == "concentrations":
                truth_cache[key] = read_concentrations(sample["truth"],
                                                       sample["mix"])
            else:
                truth_cache[key] = read_truth(sample["truth"], ttype)
        truth = truth_cache[key]
        extra = sample.get("extra_flags", "") or ""
        extra_flags = extra.split()

        # Pinned per sample by the first arm that runs, then enforced for every
        # later arm: compared runs must be scored over an identical reference
        # transcriptome, so a differing .quant key set is a hard error rather
        # than a silently different denominator.
        universe = None
        ref_md5 = None

        for arm in args.arms:
            wdir = args.outdir / arm
            timing = wdir / (sample["sample"] + ".time.txt")
            wdir.mkdir(parents=True, exist_ok=True)
            started = time.time()
            cmd, proc = run_one(args.binary, sample, arm, wdir,
                                args.threads, extra_flags, timing)
            row = {
                "panel": args.panel, "sample": sample["sample"], "arm": arm,
                "binary": args.binary_label, "technology": sample["technology"],
                "truth_type": sample["truth_type"], "rescaled": args.rescale,
                "cmd": " ".join(cmd),
            }

            timing_text = timing.read_text() if timing.exists() else ""
            m = TIME_WALL.search(timing_text)
            row["wall_s"] = round(parse_wall(m.group(1)), 2) if m else ""
            m = TIME_RSS.search(timing_text)
            row["peak_rss_gb"] = round(int(m.group(1)) / 1048576, 3) if m else ""
            m = TIME_EXIT.search(timing_text)
            row["exit_status"] = m.group(1) if m else str(proc.returncode)

            quant = wdir / (sample["sample"] + ".quant")
            if proc.returncode != 0 or not quant.exists():
                row["status"] = "FAILED"
                sys.stderr.write(
                    "FAILED %s/%s rc=%d\n%s\n" % (sample["sample"], arm,
                                                  proc.returncode,
                                                  proc.stderr[-2000:]))
            else:
                est = read_quant(str(quant))

                digest = hashlib.md5(
                    "\n".join(sorted(est)).encode()).hexdigest()
                if universe is None:
                    universe = shared_universe(truth, est, sample["truth_type"])
                    ref_md5 = digest
                elif digest != ref_md5:
                    raise SystemExit(
                        "reference transcriptome differs between arms for %s: "
                        "arm %s has md5 %s, expected %s" %
                        (sample["sample"], arm, digest, ref_md5))

                row.update(evaluate(truth, est, rescale, universe=universe))
                row["n_universe"] = len(universe)
                row["ref_md5"] = digest
                row["truth_orphans"] = len(set(truth) - universe)
                row["status"] = "ok"
                row["quant"] = str(quant)

                # Record what the run actually activated, so "auto enabled four
                # things silently" cannot happen again in this table.
                meta_path = wdir / (sample["sample"] + ".meta_info.json")
                if meta_path.exists():
                    meta = json.loads(meta_path.read_text())
                    diag = meta.get("coverage_diagnostics", {}) or {}
                    row["kernel"] = diag.get("technology_kernel", "")
                    row["rank_blend_active"] = (
                        diag.get("rank_blend_selection", {}) or {}).get("active", "")
                    row["learned_censor_scale_nt"] = (
                        diag.get("censoring", {}) or {}).get("learned_scale_nt", "")
                    row["coverage_s"] = (
                        meta.get("coverage_model_time", {}) or {}).get("seconds", "")
                    row["em_s"] = (meta.get("em_time", {}) or {}).get("seconds", "")
                if not args.keep_quant:
                    quant.unlink()

            writer.writerow(row)
            out.flush()
            print("[%s] %-22s %-14s %6.1fs spearman=%s" % (
                args.panel, sample["sample"], arm, time.time() - started,
                row.get("spearman", "NA")), flush=True)

    out.close()


if __name__ == "__main__":
    main()
