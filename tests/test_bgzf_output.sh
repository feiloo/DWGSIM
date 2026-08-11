#!/usr/bin/env bash

set -euo pipefail

dwgsim_bin=${DWGSIM_BIN:-./dwgsim}
bgzip_bin=${BGZIP_BIN:-samtools/bgzip}
reference_fasta=${1:-samtools/examples/ex1.fa}
read_pairs=${BGZF_TEST_READ_PAIRS:-50}

for dependency in "$dwgsim_bin" "$bgzip_bin" gzip awk cmp mktemp; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        echo "test-bgzf: required command not found: $dependency" >&2
        exit 1
    fi
done

if [[ ! -s $reference_fasta ]]; then
    echo "test-bgzf: reference FASTA not found: $reference_fasta" >&2
    exit 1
fi

test_root=$(mktemp -d "${TMPDIR:-/tmp}/dwgsim-bgzf-test.XXXXXX")
trap 'rm -rf "$test_root"' EXIT

verify_fastq() {
    local fastq=$1
    local expected_lines=$2
    local actual_lines

    if [[ ! -s $fastq ]]; then
        echo "test-bgzf: expected output was not generated: $fastq" >&2
        exit 1
    fi

    gzip --test "$fastq"
    "$bgzip_bin" -dc "$fastq" >/dev/null
    actual_lines=$(gzip --decompress --stdout "$fastq" |
        awk 'END { print NR }')
    if [[ $actual_lines -ne $expected_lines ]]; then
        echo "test-bgzf: $fastq has $actual_lines lines; expected $expected_lines" >&2
        exit 1
    fi
}

run_case() {
    local label=$1
    local output_mode=$2
    local threads=$3
    local case_directory=${test_root}/${label}
    local output_prefix=${case_directory}/reads

    mkdir -p "$case_directory"
    "$dwgsim_bin" \
        -z 13 \
        -t "$threads" \
        -N "$read_pairs" \
        -r 0 \
        -M 1 \
        -o "$output_mode" \
        "$reference_fasta" \
        "$output_prefix"

    case $output_mode in
        0)
            verify_fastq "${output_prefix}.bfast.fastq.gz" "$((read_pairs * 8))"
            verify_fastq "${output_prefix}.bwa.read1.fastq.gz" "$((read_pairs * 4))"
            verify_fastq "${output_prefix}.bwa.read2.fastq.gz" "$((read_pairs * 4))"
            ;;
        1)
            verify_fastq "${output_prefix}.bwa.read1.fastq.gz" "$((read_pairs * 4))"
            verify_fastq "${output_prefix}.bwa.read2.fastq.gz" "$((read_pairs * 4))"
            ;;
        2)
            verify_fastq "${output_prefix}.bfast.fastq.gz" "$((read_pairs * 8))"
            ;;
        *)
            echo "test-bgzf: unsupported output mode: $output_mode" >&2
            exit 2
            ;;
    esac
}

run_case bwa-single 1 1
run_case bwa-threaded 1 4
run_case all-threaded 0 4
run_case bfast-threaded 2 4

for read_end in read1 read2; do
    cmp \
        <(gzip --decompress --stdout \
            "${test_root}/bwa-single/reads.bwa.${read_end}.fastq.gz") \
        <(gzip --decompress --stdout \
            "${test_root}/bwa-threaded/reads.bwa.${read_end}.fastq.gz")
done

if "$dwgsim_bin" -t 0 -N 1 "$reference_fasta" \
    "${test_root}/invalid-thread-count" >/dev/null 2>&1; then
    echo "test-bgzf: -t 0 unexpectedly succeeded" >&2
    exit 1
fi

echo "BGZF FASTQ tests passed."
