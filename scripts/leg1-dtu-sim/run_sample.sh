#!/usr/bin/env bash
# Generate one Leg-1 sample end-to-end (see spec.md).
# Usage: run_sample.sh <sample> <seed> <abundance.tsv> <n_molecules> <outdir>
set -euo pipefail

SAMPLE=$1; SEED=$2; ABUND=$3; NMOL=$4; OUT=$5
ROOT=/scratch1/rob/long-read-ecosystem
TKSM=$ROOT/tools/tksm-env/bin/tksm
GTF=$ROOT/rebench-2026-08-03/junc/refseq.gtf
GENOME=$ROOT/test_data/GCF_000001405.40_GRCh38.p14_genomic.fna
KDE=$ROOT/oarfish-evaluation-data/leg1-dtu-sim/models/SQ2_kde.json
mkdir -p "$OUT"
W=$OUT/$SAMPLE

echo "[$SAMPLE] transcribe $(date '+%T')"
$TKSM transcribe -g "$GTF" -a "$ABUND" --use-whole-id --non-coding \
    --molecule-count "$NMOL" --molecule-prefix "M${SAMPLE}_" \
    -s $((SEED * 10 + 1)) -o "$W.molecules.mdf"
echo "[$SAMPLE] polyA $(date '+%T')"
$TKSM polyA -i "$W.molecules.mdf" -o "$W.polya.mdf" \
    --normal 30,7 --min-length 10 -s $((SEED * 10 + 2))
echo "[$SAMPLE] truncate $(date '+%T')"
$TKSM truncate -i "$W.polya.mdf" -o "$W.trunc.mdf" \
    --kde-model "$KDE" -s $((SEED * 10 + 3))
echo "[$SAMPLE] shuffle $(date '+%T')"
$TKSM shuffle -i "$W.trunc.mdf" -o "$W.shuf.mdf" -s $((SEED * 10 + 4))
echo "[$SAMPLE] sequence $(date '+%T')"
$TKSM sequence -i "$W.shuf.mdf" -r "$GENOME" -o "$W.fastq" \
    -t 16 --badread-error-model pacbio2016 --badread-qscore-model pacbio2016 \
    --badread-identity 99.19,99.99,2.09
echo "[$SAMPLE] truth map + compress $(date '+%T')"
# molecule -> tid from the transcribe mdf; uuid -> molecule from fastq headers
python3 - "$W" <<'EOF'
import sys, gzip, re
w = sys.argv[1]
tid_re = re.compile(r'tid=([^;,\s]+)')
m2t = {}
for line in open(f'{w}.molecules.mdf'):
    if line.startswith('+'):
        f = line.split('\t')
        m = tid_re.search(f[2]) if len(f) > 2 else None
        if m:
            m2t[f[0][1:]] = m.group(1)
out = gzip.open(f'{w}.read_truth.tsv.gz', 'wt')
n = miss = 0
with open(f'{w}.fastq') as fq:
    for i, line in enumerate(fq):
        if i % 4:
            continue
        p = line.rstrip().split()
        mid = next((x[12:] for x in p if x.startswith('molecule_id=')), None)
        t = m2t.get(mid) if mid else None
        if t:
            out.write(f'{p[0][1:]}\t{t}\n')
            n += 1
        else:
            miss += 1
out.close()
print(f'truth map: {n} reads, {miss} unmatched')
EOF
gzip -f "$W.fastq"
rm -f "$W.molecules.mdf" "$W.polya.mdf" "$W.trunc.mdf" "$W.shuf.mdf"
echo "[$SAMPLE] done $(date '+%T')"
