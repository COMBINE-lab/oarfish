#!/usr/bin/env python3
"""Aggregate the 9-arm re-benchmark and apply the plan's decision rules.

Everything is measured against `logistic`, per decision rule 1: beating `none`
is not evidence. Reports per-sample and per-family, never panel means alone.
"""

import argparse
import csv
import math
from collections import OrderedDict, defaultdict
from pathlib import Path

ARM_ORDER = ["none", "logistic", "adaptive-bare", "auto-new", "auto-old",
             "+calib", "+censor", "+prune", "+rank"]

FEATURE_ARM = OrderedDict([
    ("rank blending", "+rank"),
    ("dominance pruning", "+prune"),
    ("alignment calibration", "+calib"),
    ("censoring", "+censor"),
])


def family(row):
    """Dataset family. The regression is family-specific, so this is the unit
    decision rule 1's 'no family regressing more than 0.01 Spearman' applies to."""
    s = row["sample"]
    if row["panel"] == "B":
        if s.startswith("tksm"):
            return "B:tksm-pacbio"
        return "B:nanosim-" + ("drna" if "drna" in s else "cdna")
    if row["panel"] == "A-sirv":
        if s.startswith("independent"):
            return "A:independent-sim"
        # E0 is equal-molar: truth variance is 0, so every rank/correlation
        # metric is undefined there and only MARD / dispersion carry signal.
        return "A:sirv-E0" if "-e0-" in s else "A:sirv-E2"
    if s == "synthetic-drna":
        return "A:synthetic-drna"
    return "A:" + row["technology"]


def load(paths):
    rows = []
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("status") != "ok":
                    rows.append(r)
                    continue
                for k in ("spearman", "kendall", "pearson_log1p", "ccc_log1p",
                          "rmse", "mard", "spearman_inner", "mard_inner",
                          "wall_s", "peak_rss_gb"):
                    try:
                        r[k] = float(r[k])
                    except (TypeError, ValueError):
                        r[k] = float("nan")
                r["family"] = family(r)
                rows.append(r)
    return rows


def mean(vals):
    vals = [v for v in vals if v == v]
    return sum(vals) / len(vals) if vals else float("nan")


def index(rows):
    """(panel, sample, arm) -> row, for ok rows only."""
    return {(r["panel"], r["sample"], r["arm"]): r
            for r in rows if r.get("status") == "ok"}


def per_sample_delta(rows, metric="spearman"):
    """arm -> list of (panel, sample, family, delta vs logistic)."""
    idx = index(rows)
    out = defaultdict(list)
    for (panel, sample, arm), r in idx.items():
        base = idx.get((panel, sample, "logistic"))
        if base is None:
            continue
        out[arm].append((panel, sample, r["family"],
                         r[metric] - base[metric]))
    return out


def fmt(v, nd=4):
    if v != v:
        return "NA"
    return ("%%.%df" % nd) % v


def table_by_family(rows, metric, arms, out):
    fams = sorted({r["family"] for r in rows if r.get("status") == "ok"})
    out.append("| family | n | " + " | ".join(arms) + " |")
    out.append("|" + "---|" * (len(arms) + 2))
    for fam in fams:
        sel = [r for r in rows if r.get("status") == "ok" and r["family"] == fam]
        n = len({r["sample"] for r in sel})
        cells = []
        for arm in arms:
            cells.append(fmt(mean([r[metric] for r in sel if r["arm"] == arm])))
        out.append("| %s | %d | %s |" % (fam, n, " | ".join(cells)))
    return out


def table_per_sample(rows, metric, arms, panel):
    """Per-sample table. Decision rule: never report panel means alone, since
    the known regression is dataset-family-specific and a mean over 28 cases
    hides it entirely."""
    sel = [r for r in rows if r["panel"] == panel and r.get("status") == "ok"]
    samples = sorted({r["sample"] for r in sel})
    idx = {(r["sample"], r["arm"]): r for r in sel}
    out = ["| sample | " + " | ".join(arms) + " | best |",
           "|" + "---|" * (len(arms) + 2)]
    for s in samples:
        vals = []
        for a in arms:
            r = idx.get((s, a))
            vals.append(r[metric] if r else float("nan"))
        best = arms[max(range(len(vals)), key=lambda i: (vals[i] if vals[i] == vals[i] else -9e9))]
        cells = []
        for a, v in zip(arms, vals):
            cells.append(("**%s**" % fmt(v)) if a == best else fmt(v))
        out.append("| %s | %s | %s |" % (s, " | ".join(cells), best))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="+", type=Path)
    ap.add_argument("--metric", default="spearman")
    ap.add_argument("--per-sample", action="store_true")
    args = ap.parse_args()

    rows = load(args.results)
    ok = [r for r in rows if r.get("status") == "ok"]
    failed = [r for r in rows if r.get("status") != "ok"]
    arms = [a for a in ARM_ORDER if any(r["arm"] == a for r in ok)]

    print("runs: %d ok, %d failed" % (len(ok), len(failed)))
    for r in failed:
        print("  FAILED %s/%s/%s" % (r["panel"], r["sample"], r["arm"]))
    print()

    for panel in sorted({r["panel"] for r in ok}):
        sub = [r for r in ok if r["panel"] == panel]
        print("=== Panel %s: mean %s by family ===" % (panel, args.metric))
        out = []
        table_by_family(sub, args.metric, arms, out)
        print("\n".join(out))
        print()
        if args.per_sample:
            print("=== Panel %s: per-sample %s ===" % (panel, args.metric))
            print("\n".join(table_per_sample(ok, args.metric, arms, panel)))
            print()

    # Head-to-head vs logistic
    deltas = per_sample_delta(ok, args.metric)
    print("=== vs logistic (%s) ===" % args.metric)
    print("| arm | n | mean d | median d | wins | losses | worst family (mean d) |")
    print("|---|---|---|---|---|---|---|")
    for arm in arms:
        if arm == "logistic":
            continue
        ds = deltas.get(arm, [])
        if not ds:
            continue
        vals = sorted(d for _, _, _, d in ds)
        byfam = defaultdict(list)
        for _, _, fam, d in ds:
            byfam[fam].append(d)
        fam_means = {f: mean(v) for f, v in byfam.items()}
        worst = min(fam_means.items(), key=lambda kv: kv[1])
        med = vals[len(vals) // 2] if len(vals) % 2 else \
            (vals[len(vals) // 2 - 1] + vals[len(vals) // 2]) / 2
        print("| %s | %d | %s | %s | %d | %d | %s (%s) |" % (
            arm, len(vals), fmt(mean(vals)), fmt(med),
            sum(1 for v in vals if v > 0), sum(1 for v in vals if v < 0),
            worst[0], fmt(worst[1])))
    print()

    # Decision rule 1 / 3 evaluation per feature
    print("=== decision rule 1: beats logistic on union, no family regressing >0.01 ===")
    print("| feature | arm | union mean d | worst family d | marginal vs adaptive-bare | verdict |")
    print("|---|---|---|---|---|---|")
    idx = index(ok)
    for feat, arm in FEATURE_ARM.items():
        ds = deltas.get(arm, [])
        if not ds:
            continue
        byfam = defaultdict(list)
        for _, _, fam, d in ds:
            byfam[fam].append(d)
        fam_means = {f: mean(v) for f, v in byfam.items()}
        worst = min(fam_means.items(), key=lambda kv: kv[1])
        union = mean([d for _, _, _, d in ds])
        marg = []
        for (panel, sample, a), r in idx.items():
            if a != arm:
                continue
            bare = idx.get((panel, sample, "adaptive-bare"))
            if bare:
                marg.append(r[args.metric] - bare[args.metric])
        passes = union > 0 and worst[1] >= -0.01
        print("| %s | %s | %s | %s (%s) | %s | %s |" % (
            feat, arm, fmt(union), fmt(worst[1]), worst[0], fmt(mean(marg)),
            "DEFAULT-ON" if passes else "FAILS RULE 1"))
    print()

    # Rule 4
    ds = deltas.get("adaptive-bare", [])
    if ds:
        byfam = defaultdict(list)
        for _, _, fam, d in ds:
            byfam[fam].append(d)
        print("=== decision rule 4: does adaptive-bare beat logistic? ===")
        print("union mean d = %s" % fmt(mean([d for _, _, _, d in ds])))
        for f, v in sorted(byfam.items()):
            print("  %-22s %s" % (f, fmt(mean(v))))
        print()

    # Rule 5: panel agreement
    print("=== decision rule 5: panel agreement (mean d vs logistic) ===")
    print("| arm | Panel A | Panel B | agree? |")
    print("|---|---|---|---|")
    for arm in arms:
        if arm == "logistic":
            continue
        ds = deltas.get(arm, [])
        a = mean([d for p, _, _, d in ds if p == "A"])
        b = mean([d for p, _, _, d in ds if p == "B"])
        agree = "yes" if (a == a and b == b and (a > 0) == (b > 0)) else "NO"
        print("| %s | %s | %s | %s |" % (arm, fmt(a), fmt(b), agree))


if __name__ == "__main__":
    main()
