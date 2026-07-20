#!/usr/bin/env python3
"""Write equal-molar transcript truth for a SIRV E0 annotation."""

import argparse
import re


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf")
    parser.add_argument("output")
    args = parser.parse_args()
    names = set()
    with open(args.gtf, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) == 9 and fields[0].startswith("SIRV"):
                match = re.search(r'(?:^|;\s*)transcript_id\s+"([^"]+)"', fields[8])
                if match:
                    names.add(match.group(1))
    with open(args.output, "w", encoding="utf-8") as output:
        for name in sorted(names):
            output.write(f"{name}\t1\n")


if __name__ == "__main__":
    main()
