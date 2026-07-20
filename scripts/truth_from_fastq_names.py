#!/usr/bin/env python3
"""Count source transcripts encoded as the first FASTQ header token."""

import argparse
import gzip


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq")
    parser.add_argument("output")
    parser.add_argument("--delimiter", default="|")
    args = parser.parse_args()
    opener = gzip.open if args.fastq.endswith(".gz") else open
    counts = {}
    with opener(args.fastq, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle):
            if line_number % 4 == 0:
                name = line[1:].split()[0].split(args.delimiter, 1)[0]
                counts[name] = counts.get(name, 0) + 1
    with open(args.output, "w", encoding="utf-8") as output:
        for name in sorted(counts):
            output.write(f"{name}\t{counts[name]}\n")


if __name__ == "__main__":
    main()
