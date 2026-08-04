#!/usr/bin/env python3
"""Project the SIRV concentration table onto a single mix.

`sirv_concentrations.tsv` holds one column per mix (E0/E1/E2), but
`evaluate_quant.read_truth` reads a headerless two-column `name<TAB>value`
file. This writes that projection so SIRV-based samples can be scored by the
same driver as every other truth-bearing panel.

    python3 scripts/sirv_mix_truth.py sirv_concentrations.tsv E2 sirv_truth_E2.tsv
"""
import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("concentrations")
    parser.add_argument("mix", help="column to project, e.g. E0, E1 or E2")
    parser.add_argument("output")
    args = parser.parse_args()

    with open(args.concentrations, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"{args.concentrations} has no rows")
    name_field = list(rows[0])[0]
    if args.mix not in rows[0]:
        raise SystemExit(f"mix {args.mix!r} not in columns {list(rows[0])}")

    written = 0
    with open(args.output, "w", encoding="utf-8") as out:
        for row in rows:
            value = (row[args.mix] or "").strip()
            if not value:
                continue
            # A mix may omit a transcript entirely; a zero-concentration entry is
            # still informative (it should not be expressed) so it is retained.
            out.write(f"{row[name_field]}\t{float(value)}\n")
            written += 1
    print(f"wrote {written} transcripts for mix {args.mix} -> {args.output}")


if __name__ == "__main__":
    main()
