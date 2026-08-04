#!/usr/bin/env python3
"""Summarize coverage-ablation results split by the *provenance of their truth*.

The 28-case ablation panel is 27 libraries scored against a matched Illumina
comparator and one scored against molecular truth. A pooled mean is therefore
~96% comparator, and on three separate occasions a candidate that improved the
comparator mean was rejected by molecular truth:

  * unique-read-only profiles (2026-07-21): better Spearman/MARD in all 24
    LongBench libraries, synthetic CCC 0.99704 -> 0.94994
  * endpoint nucleotide measure (2026-07-24): ONT dRNA Spearman +0.000515 (8/9)
    while true-transcript discrimination fell 0.043
  * score temperature (2026-07-24): comparator prefers D -> 0.5, molecular truth
    prefers D -> 20, both monotonically

This script reports the two tiers separately and gates promotion on the
truth-bearing tier, so a comparator-driven gain cannot be adopted silently.

Usage
-----
    python3 scripts/summarize_panel.py RESULTS.tsv [RESULTS.tsv ...] \
        --manifest benchmarks/coverage_ablation_manifest.tsv \
        [--manifest benchmarks/truth_tier_manifest.tsv] \
        [--control full] [--candidate ARM ...]
"""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

# Primary metrics per docs/evaluation.md; `higher_is_better` drives win counts.
METRICS = (
    ("ccc", True, True),
    ("mard", False, True),
    ("spearman", True, False),
    ("pearson", True, False),
    ("rmse", False, False),
)
PRIMARY = tuple(name for name, _, primary in METRICS if primary)
TIER_OF_TRUTH_TYPE = {"counts": "truth", "illumina": "comparator"}


def load_tiers(manifests):
    """sample -> tier, from each manifest's `truth_type` column."""
    tiers = {}
    for path in manifests:
        with open(path, encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                truth_type = (row.get("truth_type") or "").strip()
                tiers[row["sample"]] = TIER_OF_TRUTH_TYPE.get(truth_type, "comparator")
    return tiers


def load_results(paths):
    rows = []
    for path in paths:
        with open(path, encoding="utf-8", newline="") as handle:
            rows.extend(csv.DictReader(handle, delimiter="\t"))
    # Average repeats so a candidate run twice does not outvote one run once.
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["sample"], row["technology"], row["ablation"])].append(row)
    merged = {}
    for key, group in grouped.items():
        merged[key] = {
            name: statistics.mean(float(r[name]) for r in group)
            for name, _, _ in METRICS
        }
    return merged


def summarize(merged, tiers, control, candidate):
    """Per-tier deltas for one candidate arm versus the control arm.

    Also keyed by (tier, technology), because a change that applies to only one
    technology is diluted in a tier mean by the technologies it cannot affect.
    """
    per_tier = defaultdict(lambda: defaultdict(list))
    for (sample, technology, ablation), values in merged.items():
        if ablation != candidate:
            continue
        base = merged.get((sample, technology, control))
        if base is None:
            continue
        tier = tiers.get(sample, "comparator")
        for key in (tier, (tier, technology)):
            for name, _, _ in METRICS:
                per_tier[key][name].append(values[name] - base[name])
            per_tier[key]["_samples"].append(sample)
    return per_tier


def render(per_tier, candidate, control, tolerance):
    print(f"\n=== {candidate}  vs  {control} ===")
    verdict_lines = []
    for tier in ("truth", "comparator"):
        deltas = per_tier.get(tier)
        if not deltas:
            print(f"  [{tier}] no samples")
            continue
        n = len(deltas["_samples"])
        print(f"  [{tier}]  n={n}  ({', '.join(sorted(deltas['_samples'])[:4])}"
              f"{', ...' if n > 4 else ''})")
        # A sample the change cannot reach contributes an exact zero. Those
        # zeros drag the median to zero regardless of the candidate's merit, so
        # primary metrics are judged over *affected* samples only.
        affected = [i for i in range(n)
                    if any(abs(deltas[name][i]) > 1e-12 for name in PRIMARY)]
        inert = n - len(affected)
        if inert:
            print(f"    {inert}/{n} samples inert (exact zero on every primary metric)")
        for name, higher_is_better, primary in METRICS:
            values = deltas[name]
            mean = statistics.mean(values)
            wins = sum(1 for v in values if (v > 0) == higher_is_better and v != 0)
            mark = "*" if primary else " "
            line = f"    {mark} {name:<9}{mean:+12.6f}   wins {wins}/{n}"
            if primary and affected:
                sub = [values[i] for i in affected]
                line += (f"   | affected n={len(sub)}"
                         f" median {statistics.median(sub):+.6f}")
            print(line)
            if tier == "truth" and primary and affected:
                sub = [values[i] for i in affected]
                median = statistics.median(sub)
                regressed = -median if higher_is_better else median
                if regressed > tolerance:
                    verdict_lines.append(
                        f"{name} median over affected samples regresses by "
                        f"{regressed:.6f} (> {tolerance})")
        # A technology the change cannot reach contributes only zeros, which
        # drags the tier mean toward zero; show the split so that is visible.
        techs = sorted(key[1] for key in per_tier
                       if isinstance(key, tuple) and key[0] == tier)
        if len(techs) > 1:
            print(f"    by technology (CCC / Spearman, n):")
            for tech in techs:
                sub = per_tier[(tier, tech)]
                k = len(sub["_samples"])
                inert = all(abs(v) < 1e-12 for v in sub["ccc"])
                note = "   [inert - change does not reach this technology]" if inert else ""
                print(f"      {tech:<10}{statistics.mean(sub['ccc']):+12.6f}"
                      f"{statistics.mean(sub['spearman']):+12.6f}  n={k}{note}")
    print("  (* = primary metric per docs/evaluation.md)")
    if not per_tier.get("truth"):
        print("  VERDICT: INDETERMINATE - no truth-bearing samples in these results.")
        print("           Add benchmarks/truth_tier_manifest.tsv to the run.")
    elif verdict_lines:
        print("  VERDICT: REJECT - truth tier regresses:")
        for line in verdict_lines:
            print(f"           - {line}")
    else:
        print("  VERDICT: truth tier clean; judge the comparator tier on its merits.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results", nargs="+", type=Path)
    parser.add_argument("--manifest", action="append", required=True, type=Path,
                        help="may be repeated; supplies each sample's truth_type")
    parser.add_argument("--control", default="full")
    parser.add_argument("--candidate", nargs="*",
                        help="default: every non-control arm present")
    parser.add_argument("--truth-tolerance", type=float, default=0.0,
                        help="allowed mean regression of a primary metric in the "
                             "truth tier before the candidate is rejected")
    args = parser.parse_args()

    tiers = load_tiers(args.manifest)
    merged = load_results(args.results)
    arms = sorted({ablation for _, _, ablation in merged})
    if args.control not in arms:
        raise SystemExit(f"control arm {args.control!r} not present; have {arms}")
    candidates = args.candidate or [a for a in arms if a != args.control]

    counts = defaultdict(int)
    for sample, _, _ in merged:
        counts[tiers.get(sample, "comparator")] += 1
    unique = {s for s, _, _ in merged}
    print(f"arms: {', '.join(arms)}")
    print("tier composition: " + ", ".join(
        f"{t}={len({s for s in unique if tiers.get(s, 'comparator') == t})}"
        for t in ("truth", "comparator")))

    for candidate in candidates:
        render(summarize(merged, tiers, args.control, candidate),
               candidate, args.control, args.truth_tolerance)


if __name__ == "__main__":
    main()
