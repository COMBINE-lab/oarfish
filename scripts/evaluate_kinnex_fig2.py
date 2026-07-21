#!/usr/bin/env python3
"""Reproduce Kinnex Figure 2D/E metrics for Oarfish coverage variants."""

import argparse
import csv
import math
import statistics
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


NS = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def pearson(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    xx = sum((v - mx) ** 2 for v in x)
    yy = sum((v - my) ** 2 for v in y)
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(xx * yy)


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + end - 1) / 2 + 1
        for idx in order[start:end]:
            result[idx] = rank
        start = end
    return result


def xlsx_truth(path):
    with zipfile.ZipFile(path) as archive:
        shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        shared = ["".join(node.itertext()) for node in shared_root.findall("x:si", NS)]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    truth = {}
    for row in sheet.findall(".//x:row", NS):
        cells = {}
        for cell in row.findall("x:c", NS):
            ref = cell.attrib["r"]
            column = "".join(c for c in ref if c.isalpha())
            value = cell.find("x:v", NS)
            if value is None:
                continue
            raw = value.text
            cells[column] = shared[int(raw)] if cell.attrib.get("t") == "s" else raw
        # readxl starts at the populated B column, so the paper's ...6/...8/
        # ...11/...12 fields correspond to worksheet G/I/L/M.
        name = cells.get("G", "")
        if name.startswith("SIRV") and name[4:].isdigit():
            try:
                truth[name] = {
                    "mw": float(cells["I"]), "e1": float(cells["L"]),
                    "e2": float(cells["M"]),
                }
            except (KeyError, ValueError):
                pass
    return truth


def lengths_from_gtf(path):
    lengths = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = dict(
                part.strip().split(" ", 1) for part in fields[8].split(";")
                if part.strip() and " " in part
            )
            name = attrs["transcript_id"].strip('"')
            if name.startswith("SIRV") and name[4:].isdigit():
                lengths[name] = lengths.get(name, 0) + int(fields[4]) - int(fields[3]) + 1
    return lengths


def quant_counts(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["tname"]: float(row["num_reads"])
                for row in csv.DictReader(handle, delimiter="\t")}


def formatted_counts(path):
    result = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["transcript_id"]] = {
                sample: float(row[sample])
                for sample in ("day0-rep1", "day0-rep2", "day0-rep3",
                               "day5-rep1", "day5-rep2", "day5-rep3")
            }
    return result


def log_cpm(columns, names):
    totals = {sample: sum(counts.get(name, 0.0) for name in names)
              for sample, counts in columns.items()}
    mean_total = sum(totals.values()) / len(totals)
    result = {}
    for sample, counts in columns.items():
        scaled_prior = 2.0 * totals[sample] / mean_total
        adjusted_total = totals[sample] + 2.0 * scaled_prior
        result[sample] = {
            name: math.log2((counts.get(name, 0.0) + scaled_prior)
                            / adjusted_total * 1_000_000.0)
            for name in names
        }
    return result


def metrics(model, columns, truth, lengths, writer):
    names = sorted(set(truth) & set(lengths))
    cpm = log_cpm(columns, names)
    expected_d, observed_d, d_lengths = [], [], []
    for rep in range(1, 4):
        for name in names:
            expected_d.append(math.log2(truth[name]["e1"] / truth[name]["e2"]))
            observed_d.append(cpm[f"day0-rep{rep}"][name] - cpm[f"day5-rep{rep}"][name])
            d_lengths.append(lengths[name])
    expected_e, observed_e, e_lengths = [], [], []
    for rep in range(1, 4):
        for name in names:
            expected_e.append(math.log(truth[name]["mw"] * truth[name]["e1"]))
            observed_e.append(cpm[f"day0-rep{rep}"][name])
            e_lengths.append(lengths[name])
    for panel, expected, observed, lens in (
        ("D", expected_d, observed_d, d_lengths),
        ("E", expected_e, observed_e, e_lengths),
    ):
        for stratum, keep in (
            ("all", [True] * len(lens)),
            ("short", [length < 1250 for length in lens]),
            ("long", [length >= 1250 for length in lens]),
        ):
            x = [v for v, use in zip(expected, keep) if use]
            y = [v for v, use in zip(observed, keep) if use]
            residual_rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(x, y)) / len(x)) if panel == "D" else ""
            writer.writerow({"model": model, "panel": panel, "stratum": stratum,
                             "observations": len(x), "pearson": pearson(x, y),
                             "spearman": pearson(ranks(x), ranks(y)),
                             "rmse": residual_rmse})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("quant_dir", type=Path)
    parser.add_argument("truth_xlsx", type=Path)
    parser.add_argument("gtf", type=Path)
    parser.add_argument("published_counts", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--models", nargs="+", default=("none", "logistic", "auto", "abundance-blend"))
    args = parser.parse_args()
    truth = xlsx_truth(args.truth_xlsx)
    lengths = lengths_from_gtf(args.gtf)
    fields = ("model", "panel", "stratum", "observations", "pearson", "spearman", "rmse")
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        published = formatted_counts(args.published_counts)
        columns = {sample: {name: row[sample] for name, row in published.items()}
                   for sample in next(iter(published.values()))}
        metrics("published-logistic", columns, truth, lengths, writer)
        for model in args.models:
            columns = {
                f"{day}-rep{rep}": quant_counts(
                    args.quant_dir / f"{day}-rep{rep}-{model}.quant")
                for day in ("day0", "day5") for rep in range(1, 4)
            }
            metrics(model, columns, truth, lengths, writer)
    print(f"truth={len(truth)} SIRVs; lengths={len(lengths)} SIRVs")


if __name__ == "__main__":
    main()
