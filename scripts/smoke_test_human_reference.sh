#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <reference.fasta> <output-directory>" >&2
    exit 2
fi

reference_fasta=$1
output_directory=$2
dwgsim_bin=${DWGSIM_BIN:-./dwgsim}
read_pairs=${HUMAN_SMOKE_READ_PAIRS:-100}
random_seed=${HUMAN_SMOKE_SEED:-13}

for dependency in "$dwgsim_bin" gzip awk; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        echo "test-human-reference: required command not found: $dependency" >&2
        exit 1
    fi
done

if [[ ! -s $reference_fasta ]]; then
    echo "test-human-reference: reference FASTA not found: $reference_fasta" >&2
    echo "Run 'make download' first." >&2
    exit 1
fi
if [[ ! $read_pairs =~ ^[1-9][0-9]*$ ]]; then
    echo "test-human-reference: HUMAN_SMOKE_READ_PAIRS must be a positive integer" >&2
    exit 2
fi

mkdir -p "$output_directory"
output_prefix=${output_directory}/grch38-p14-smoke

echo "Simulating $read_pairs read pairs against the complete GRCh38.p14 reference..."
"$dwgsim_bin" \
    -z "$random_seed" \
    -N "$read_pairs" \
    -r 0 \
    -M 1 \
    "$reference_fasta" \
    "$output_prefix"

bfast_fastq=${output_prefix}.bfast.fastq.gz
bwa_read1_fastq=${output_prefix}.bwa.read1.fastq.gz
bwa_read2_fastq=${output_prefix}.bwa.read2.fastq.gz

for fastq in "$bfast_fastq" "$bwa_read1_fastq" "$bwa_read2_fastq"; do
    if [[ ! -s $fastq ]]; then
        echo "test-human-reference: expected FASTQ was not generated: $fastq" >&2
        exit 1
    fi
    gzip --test "$fastq"
done

check_fastq_lines() {
    local fastq=$1
    local expected_lines=$2
    local actual_lines

    actual_lines=$(gzip --decompress --stdout "$fastq" | awk 'END { print NR }')
    if [[ $actual_lines -ne $expected_lines ]]; then
        echo "test-human-reference: $fastq has $actual_lines lines; expected $expected_lines" >&2
        exit 1
    fi
}

check_fastq_lines "$bfast_fastq" "$((read_pairs * 8))"
check_fastq_lines "$bwa_read1_fastq" "$((read_pairs * 4))"
check_fastq_lines "$bwa_read2_fastq" "$((read_pairs * 4))"

echo "GRCh38 smoke test passed. Output: $output_directory"
