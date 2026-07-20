#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 4 ]]; then
    echo "usage: $0 NAME_SORTED_BAM OUTPUT_DIR [THREADS] [SEQ_TECH]" >&2
    exit 2
fi

bam=$1
output_dir=$2
threads=${3:-4}
technology=${4:-}
oarfish_bin=${OARFISH_BIN:-target/release/oarfish}
mkdir -p "$output_dir"

models=(none logistic endpoint hybrid adaptive)
if [[ -n $technology ]]; then
    models+=(auto)
fi
if [[ $technology == ont-drna ]]; then
    models+=(degradation)
fi
for model in "${models[@]}"; do
    technology_args=()
    if [[ -n $technology ]]; then
        technology_args=(--seq-tech "$technology")
    fi
    RUST_LOG=${RUST_LOG:-warn} "$oarfish_bin" \
        --alignments "$bam" \
        --output "$output_dir/$model" \
        --coverage-model "$model" \
        --filter-group no-filters \
        "${technology_args[@]}" \
        --threads "$threads"
done

for model in "${models[@]}"; do
    jq -c --arg model "$model" '{
        model: $model,
        coverage_seconds: .coverage_model_time.seconds,
        em_seconds: .em_time.seconds,
        coverage_diagnostics: .coverage_diagnostics
    }' "$output_dir/$model.meta_info.json"
done
