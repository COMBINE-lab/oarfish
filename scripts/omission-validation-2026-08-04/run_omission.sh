#!/bin/bash
# Detection validation for --model-unannotated-isoforms, from the existing
# genome BAM (NanoSim read names carry per-read truth).
#   full     : correct annotation -> false-positive floor
#   origin50 : 50% of expressed transcripts deleted -> detection precision
# 8 threads: correctness run, not a performance measurement.
set -u
cd /scratch1/rob/long-read-ecosystem/rebench-2026-08-03
export TMPDIR=/scratch2/tmp
B=/scratch1/rob/long-read-ecosystem/oarfish-demote/target-demote/release/oarfish
GBAM=/scratch1/rob/long-read-ecosystem/eval/parity/genome.bam
GFA=/scratch1/rob/long-read-ecosystem/test_data/GCF_000001405.40_GRCh38.p14_genomic.fna
FULL=/scratch1/rob/long-read-ecosystem/eval/GCF.pc_lncrna.matched.gtf

echo "--- full annotation, model ON (false-positive floor)"
$B --genome-alignments "$GBAM" --annotation "$FULL" --genome-fasta "$GFA" \
   --model-unannotated-isoforms --filter-group no-filters -j 8 \
   -o omission/full_on > omission/full_on.log 2>&1 && echo "   ok" || echo "   FAILED"

echo "--- origin50 holdout, model ON"
$B --genome-alignments "$GBAM" --annotation omission/origin50.gtf --genome-fasta "$GFA" \
   --model-unannotated-isoforms --filter-group no-filters -j 8 \
   -o omission/o50_on > omission/o50_on.log 2>&1 && echo "   ok" || echo "   FAILED"

echo "--- origin50 holdout, model OFF (baseline for error comparison)"
$B --genome-alignments "$GBAM" --annotation omission/origin50.gtf --genome-fasta "$GFA" \
   --filter-group no-filters -j 8 \
   -o omission/o50_off > omission/o50_off.log 2>&1 && echo "   ok" || echo "   FAILED"

echo OMISSION_DONE
