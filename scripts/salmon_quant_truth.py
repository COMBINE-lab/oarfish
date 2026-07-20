#!/usr/bin/env python3
"""Convert a Salmon-style quantification table to evaluator truth/count TSV."""

import argparse
import csv
import gzip


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("--id-delimiter", default="|")
    args = parser.parse_args()
    opener = gzip.open if args.input.endswith(".gz") else open
    with opener(args.input, "rt", encoding="utf-8", newline="") as source, open(
        args.output, "w", encoding="utf-8"
    ) as destination:
        for row in csv.DictReader(source, delimiter="\t"):
            name = row["Name"].split(args.id_delimiter, 1)[0]
            destination.write(f"{name}\t{row['NumReads']}\n")


if __name__ == "__main__":
    main()
