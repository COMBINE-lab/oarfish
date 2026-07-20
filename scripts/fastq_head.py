#!/usr/bin/env python3
"""Copy the first N records from a plain or gzip-compressed FASTQ."""

import argparse
import gzip
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("records", type=int)
    args = parser.parse_args()
    source = (
        sys.stdin
        if args.input == "-"
        else (gzip.open(args.input, "rt", encoding="utf-8")
              if args.input.endswith(".gz")
              else open(args.input, encoding="utf-8"))
    )
    with source, open(args.output, "w", encoding="utf-8") as destination:
        for _ in range(args.records * 4):
            line = source.readline()
            if not line:
                break
            destination.write(line)


if __name__ == "__main__":
    main()
