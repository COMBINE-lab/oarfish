#!/usr/bin/env python3
"""Extract transcript FASTA sequences from a genomic FASTA and GTF.

This intentionally uses only the Python standard library so evaluation
manifests do not depend on a particular gffread installation.
"""

import argparse
import gzip
import re


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")


def read_fasta(path):
    sequences = {}
    name, pieces = None, []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(pieces).upper()
                name, pieces = line[1:].split()[0], []
            else:
                pieces.append(line.strip())
    if name is not None:
        sequences[name] = "".join(pieces).upper()
    return sequences


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genome")
    parser.add_argument("gtf")
    parser.add_argument("output")
    args = parser.parse_args()

    genome = read_fasta(args.genome)
    transcripts = {}
    strands = {}
    with open_text(args.gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            match = re.search(r'(?:^|;\s*)transcript_id\s+"([^"]+)"', fields[8])
            if not match:
                continue
            transcript = match.group(1)
            transcripts.setdefault(transcript, []).append(
                (fields[0], int(fields[3]) - 1, int(fields[4]))
            )
            strands[transcript] = fields[6]

    with open(args.output, "w", encoding="utf-8") as output:
        for transcript in sorted(transcripts):
            exons = sorted(transcripts[transcript], key=lambda exon: exon[1])
            sequence = "".join(genome[seqname][start:end] for seqname, start, end in exons)
            if strands[transcript] == "-":
                sequence = reverse_complement(sequence)
            output.write(f">{transcript}\n")
            for start in range(0, len(sequence), 80):
                output.write(sequence[start:start + 80] + "\n")


if __name__ == "__main__":
    main()
