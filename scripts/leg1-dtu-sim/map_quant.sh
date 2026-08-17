#!/usr/bin/env bash
# Map one Leg-1 sample to the transcriptome and quantify both arms
# (presence stack / plain), 30 presence-aware bootstraps each, then convert
# to salmon-format dirs for fishpond/edgeR.
# Usage: map_quant.sh <sample> <fastq.gz> <outdir>
set -euo pipefail
SAMPLE=$1; FQ=$2; OUT=$3
ROOT=/scratch1/rob/long-read-ecosystem
MM2=$ROOT/tools/mm2src/minimap2
ST=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
TXOME=$ROOT/test_data/GCF_000001405.40_GRCh38.p14_rna.fna
OARFISH=$ROOT/oarfish/target/release/oarfish
CONV=$ROOT/oarfish/scripts/leg1-dtu-sim/oarfish2salmon.py
PY=/nfshomes/nomad/miniconda3/bin/python
mkdir -p "$OUT"/{bam,quant,salmon-pres,salmon-plain}

BAM=$OUT/bam/$SAMPLE.bam
if [[ ! -s $BAM.done ]]; then
  $MM2 -t 32 -ax map-hifi -N 100 "$TXOME" "$FQ" 2> "$OUT/bam/$SAMPLE.mm2.log" \
    | $ST view -@ 4 -b -o "$BAM" -
  date > "$BAM.done"
fi

COMMON="--seq-tech pac-bio-hifi --filter-group no-filters -j 16 --model-coverage \
  --em-accel squarem --convergence-l1-thresh 1e-6 --num-bootstraps 30"
$OARFISH -a "$BAM" $COMMON \
  --presence-model spike-slab --presence-endpoint-alpha 1.0 \
  -o "$OUT/quant/${SAMPLE}_pres" > "$OUT/quant/${SAMPLE}_pres.log" 2>&1
$OARFISH -a "$BAM" $COMMON \
  -o "$OUT/quant/${SAMPLE}_plain" > "$OUT/quant/${SAMPLE}_plain.log" 2>&1

$PY "$CONV" --prefix "$OUT/quant/${SAMPLE}_pres" --out "$OUT/salmon-pres/$SAMPLE"
$PY "$CONV" --prefix "$OUT/quant/${SAMPLE}_plain" --out "$OUT/salmon-plain/$SAMPLE"
echo "[$SAMPLE] map+quant done $(date '+%T')"
