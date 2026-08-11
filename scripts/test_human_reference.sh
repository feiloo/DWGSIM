#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <reference.fasta> <output-directory>" >&2
    exit 2
fi

reference_fasta=$1
output_directory=$2
gatk_bin=${GATK:-gatk}
dwgsim_bin=${DWGSIM_BIN:-./dwgsim}
samtools_bin=${SAMTOOLS_BIN:-samtools/samtools}
read_pairs=${HUMAN_TEST_READ_PAIRS:-1000}
mutation_rate=${HUMAN_TEST_MUTATION_RATE:-0.000001}
random_seed=${HUMAN_TEST_SEED:-13}

require_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "test-human-reference: required command not found: $1" >&2
        exit 1
    fi
}

require_command "$gatk_bin"
require_command "$dwgsim_bin"
require_command "$samtools_bin"
require_command gzip

if [[ ! -s $reference_fasta ]]; then
    echo "test-human-reference: reference FASTA not found: $reference_fasta" >&2
    echo "Run 'make download' first." >&2
    exit 1
fi
if [[ ! $read_pairs =~ ^[1-9][0-9]*$ ]]; then
    echo "test-human-reference: HUMAN_TEST_READ_PAIRS must be a positive integer" >&2
    exit 2
fi

reference_fai=${reference_fasta}.fai
reference_dict=${reference_fasta%.*}.dict

if [[ ! -s $reference_fai ]]; then
    echo "Creating FASTA index with SAMtools..."
    "$samtools_bin" faidx "$reference_fasta"
fi
if [[ ! -s $reference_fai ]]; then
    echo "test-human-reference: SAMtools did not create $reference_fai" >&2
    exit 1
fi

if [[ ! -s $reference_dict ]]; then
    echo "Creating sequence dictionary with GATK..."
    "$gatk_bin" CreateSequenceDictionary \
        -R "$reference_fasta" \
        -O "$reference_dict"
fi
if [[ ! -s $reference_dict ]]; then
    echo "test-human-reference: GATK did not create $reference_dict" >&2
    exit 1
fi

mkdir -p "$output_directory"
output_prefix=${output_directory}/grch38-p14

echo "Simulating $read_pairs read pairs against the complete GRCh38.p14 reference..."
"$dwgsim_bin" \
    -z "$random_seed" \
    -N "$read_pairs" \
    -r "$mutation_rate" \
    "$reference_fasta" \
    "$output_prefix"

gzip --test \
    "${output_prefix}.bfast.fastq.gz" \
    "${output_prefix}.bwa.read1.fastq.gz" \
    "${output_prefix}.bwa.read2.fastq.gz"

echo "Validating generated variants and reference alleles with GATK..."
"$gatk_bin" ValidateVariants \
    -R "$reference_fasta" \
    -V "${output_prefix}.mutations.vcf" \
    --validation-type-to-exclude IDS \
    --validation-type-to-exclude ALLELES \
    --validation-type-to-exclude CHR_COUNTS

echo "Full human-reference test passed. Output: $output_directory"
