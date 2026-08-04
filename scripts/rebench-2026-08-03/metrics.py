#!/usr/bin/env python3
"""Accuracy metrics for the coverage-model re-benchmark (2026-08-03).

Join semantics (fixed by docs/coverage-auto-rebenchmark-plan-2026-08-03.md):
  * strip an accession version suffix (`\\.\\d+$`) from both sides
  * OUTER join the truth and estimate key sets
  * fill missing entries with 0.0

Metrics: Spearman, Kendall (tau-b), Pearson on log1p, CCC on log1p, RMSE, MARD.
"""

import csv
import gzip
import math
import re

VERSION_RE = re.compile(r"\.\d+$")


def strip_key(name, delimiter="|"):
    """Normalise a transcript identifier: take the field before `delimiter`,
    then remove a trailing accession version."""
    if delimiter and delimiter in name:
        name = name.split(delimiter, 1)[0]
    return VERSION_RE.sub("", name)


def _looks_numeric(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def read_truth(path, truth_type):
    """Read a truth table into {stripped_id: value}.

    Handles the three formats present across both panels:
      * `illumina` - gzipped salmon-style TSV with Name/NumReads, Name is
        pipe-delimited and versioned.
      * `counts`   - whitespace/TSV `id<TAB>count`. May or may not carry a
        header, and may or may not be versioned (NanoSim is headerless and
        unversioned; TKSM has a `Transcript_ID/Count` header and IS versioned).
    """
    result = {}
    if truth_type == "illumina":
        with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                key = strip_key(row["Name"])
                result[key] = result.get(key, 0.0) + float(row["NumReads"])
        return result

    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh):
            fields = line.split()
            if len(fields) < 2:
                continue
            if not _looks_numeric(fields[1]):
                # header row (e.g. TKSM's "Transcript_ID\tCount"); only tolerated
                # on the first line, otherwise it is a real parse error.
                if lineno == 0:
                    continue
                raise ValueError(f"{path}:{lineno + 1}: non-numeric count {fields[1]!r}")
            key = strip_key(fields[0])
            result[key] = result.get(key, 0.0) + float(fields[1])
    return result


def read_concentrations(path, mix):
    """SIRV spike-in truth: a concentration table with one column per mixture
    (E0/E1/E2). `mix` selects the column. Matches
    scripts/run_coverage_default_benchmark.py:read_concentrations."""
    with open(path, encoding="utf-8", newline="") as fh:
        return {strip_key(row["transcript"]): float(row[mix])
                for row in csv.DictReader(fh, delimiter="\t")}


def read_quant(path):
    """Read an oarfish .quant into {stripped_id: num_reads}."""
    result = {}
    with open(path, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = strip_key(row["tname"])
            result[key] = result.get(key, 0.0) + float(row["num_reads"])
    return result


def ranks(values):
    """Average ranks, ties shared."""
    order = sorted(range(len(values)), key=values.__getitem__)
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + j - 1) / 2 + 1
        for k in range(i, j):
            out[order[k]] = rank
        i = j
    return out


def _pearson_ccc(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    vx = sum((v - mx) ** 2 for v in x) / n
    vy = sum((v - my) ** 2 for v in y) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y)) / n
    pearson = cov / math.sqrt(vx * vy) if vx > 0 and vy > 0 else 0.0
    denom = vx + vy + (mx - my) ** 2
    ccc = 2 * cov / denom if denom > 0 else 0.0
    return pearson, ccc


def kendall_tau_b(x, y):
    """Kendall tau-b via a merge-sort O(n log n) discordant count.

    tau_b = (C - D) / sqrt((n0 - n1) * (n0 - n2))
    """
    n = len(x)
    if n < 2:
        return 0.0
    pairs = sorted(zip(x, y))
    xs = [p[0] for p in pairs]
    ys = [p[1] for p in pairs]

    def tie_count(vals):
        total = 0
        i = 0
        while i < len(vals):
            j = i + 1
            while j < len(vals) and vals[j] == vals[i]:
                j += 1
            run = j - i
            total += run * (run - 1) // 2
            i = j
        return total

    n0 = n * (n - 1) // 2
    n1 = tie_count(xs)              # ties in x
    n2 = tie_count(sorted(ys))      # ties in y

    # joint ties (tied in both x and y) must not be counted as discordant
    joint = 0
    i = 0
    while i < n:
        j = i + 1
        while j < n and xs[j] == xs[i]:
            j += 1
        joint += tie_count(sorted(ys[i:j]))
        i = j

    # count discordant pairs by counting inversions in ys (x already ascending,
    # ties in x grouped and internally sorted by y so they contribute none)
    seq = []
    i = 0
    while i < n:
        j = i + 1
        while j < n and xs[j] == xs[i]:
            j += 1
        seq.extend(sorted(ys[i:j]))
        i = j

    def sort_count(arr):
        if len(arr) < 2:
            return arr, 0
        mid = len(arr) // 2
        left, a = sort_count(arr[:mid])
        right, b = sort_count(arr[mid:])
        merged, inv = [], a + b
        i = j = 0
        while i < len(left) and j < len(right):
            if left[i] <= right[j]:
                merged.append(left[i]); i += 1
            else:
                merged.append(right[j]); j += 1
                inv += len(left) - i
        merged.extend(left[i:]); merged.extend(right[j:])
        return merged, inv

    _, discordant = sort_count(seq)
    concordant = n0 - n1 - n2 + joint - discordant
    denom = math.sqrt((n0 - n1) * (n0 - n2))
    return (concordant - discordant) / denom if denom > 0 else 0.0


def _ranks_np(v):
    """Average ranks with ties shared (numpy equivalent of `ranks`)."""
    import numpy as np
    from scipy.stats import rankdata
    return rankdata(v)


def _pearson_ccc_np(x, y):
    """Pearson and CCC. Returns NaN - not 0.0 - when either vector is constant.

    This matters for the SIRV E0 mixture, which is equal-molar: every spike-in
    has concentration exactly 1, so truth variance is 0 and every correlation
    against it is undefined. Returning 0.0 there (as the older helpers did)
    would silently average six undefined samples into each arm's panel mean and
    read as a huge regression.
    """
    import numpy as np
    mx, my = x.mean(), y.mean()
    vx, vy = x.var(), y.var()
    cov = ((x - mx) * (y - my)).mean()
    pearson = (float(cov / math.sqrt(vx * vy))
               if vx > 0 and vy > 0 else float("nan"))
    denom = vx + vy + (mx - my) ** 2
    ccc = float(2 * cov / denom) if denom > 0 else float("nan")
    return pearson, ccc


def shared_universe(truth, estimate, truth_type):
    """The reference transcriptome shared by the run and the truth.

    Every arm compared against another must be scored over exactly this set, and
    a transcript may only be filled with 0 when 0 is a real measurement rather
    than an annotation gap.

    * `illumina` truth enumerates its own full annotation (every transcript,
      including the zeros), so the shared reference is the intersection of the
      two annotations. Transcripts the Illumina reference never contained are
      excluded rather than scored as truth=0.
    * `counts` truth (simulations) is a sparse list of non-zero transcripts over
      the *same* reference the BAM was built on, so the shared reference is the
      run's own transcriptome and absent truth entries are genuine zeros.
    """
    est_keys = set(estimate)
    if truth_type in ("illumina", "concentrations"):
        # Both enumerate their own complete annotation: the Illumina quant lists
        # every transcript it measured (zeros included), and the SIRV
        # concentration table lists every spike-in in the mix. Anything outside
        # is an annotation gap, not a measured zero.
        return est_keys & set(truth)
    return est_keys


def evaluate(truth, estimate, rescale, universe=None):
    """Outer-join `truth` and `estimate` (fill 0) and compute all metrics.

    `universe`, when given, restricts scoring to a fixed reference transcriptome
    (see `shared_universe`) so that compared runs are scored over an identical
    identifier set. Without it the union of the two key sets is used.

    `rescale` multiplies the estimate so its total matches truth's. Required
    when truth is a cross-protocol proxy (Panel A matched-Illumina), where the
    two vectors are in different units; must be OFF when truth and estimate are
    both simulated read counts in the same units (Panel B). Rescaling is
    computed after the universe restriction, so it normalises over exactly the
    transcripts being scored.

    Uses numpy/scipy for the heavy vectors. `kendall_tau_b` and `_pearson_ccc`
    above are the reference implementations and agree with scipy exactly over
    200 randomised tie-heavy trials.
    """
    import numpy as np
    from scipy.stats import kendalltau

    if universe is None:
        names = sorted(set(truth) | set(estimate))
    else:
        names = sorted(universe)
    x = np.fromiter((truth.get(n, 0.0) for n in names), dtype=float, count=len(names))
    y = np.fromiter((estimate.get(n, 0.0) for n in names), dtype=float, count=len(names))

    if rescale:
        sy = y.sum()
        if sy > 0:
            y = y * (x.sum() / sy)

    lx, ly = np.log1p(x), np.log1p(y)

    pearson_log, ccc_log = _pearson_ccc_np(lx, ly)
    spearman, _ = _pearson_ccc_np(_ranks_np(x), _ranks_np(y))
    kendall = float(kendalltau(x, y)[0])
    rmse = float(np.sqrt(np.mean((x - y) ** 2)))
    denom = np.abs(x) + np.abs(y)
    with np.errstate(invalid="ignore", divide="ignore"):
        ard = np.where(denom > 0, np.abs(x - y) / denom, 0.0)
    mard = float(ard.mean())

    # Secondary view restricted to transcripts the truth actually reports as
    # non-zero, i.e. dropping the zero-fill entirely.
    shared = sorted(k for k in names if truth.get(k, 0.0) > 0)
    if shared:
        ix = np.fromiter((truth[n] for n in shared), dtype=float, count=len(shared))
        iy = np.fromiter((estimate.get(n, 0.0) for n in shared),
                         dtype=float, count=len(shared))
        if rescale and iy.sum() > 0:
            iy = iy * (ix.sum() / iy.sum())
        spearman_inner, _ = _pearson_ccc_np(_ranks_np(ix), _ranks_np(iy))
        idenom = np.abs(ix) + np.abs(iy)
        with np.errstate(invalid="ignore", divide="ignore"):
            mard_inner = float(np.where(idenom > 0, np.abs(ix - iy) / idenom, 0.0).mean())
    else:
        spearman_inner, mard_inner = float("nan"), float("nan")

    # Dispersion of the estimate. For an equal-molar mixture (SIRV E0) this is
    # the metric that actually carries signal: truth is flat, so the question is
    # how uniform the estimate is, not how it correlates.
    my = y.mean()
    estimate_cv = float(y.std() / my) if my > 0 else float("nan")

    return {
        "estimate_cv": estimate_cv,
        "truth_is_constant": int(x.var() == 0),
        "n_union": len(names),
        "n_truth": len(truth),
        "n_est": len(estimate),
        "n_shared": len(shared),
        "spearman": spearman,
        "kendall": kendall,
        "pearson_log1p": pearson_log,
        "ccc_log1p": ccc_log,
        "rmse": rmse,
        "mard": mard,
        "spearman_inner": spearman_inner,
        "mard_inner": mard_inner,
    }
