#!/usr/bin/env bash

set -euo pipefail

dwgsim_bin=${DWGSIM_BIN:-./dwgsim}
bgzip_bin=${BGZIP_BIN:-samtools/bgzip}
reference_fasta=${1:-samtools/examples/ex1.fa}
read_pairs=${BGZF_TEST_READ_PAIRS:-1000}

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
    local thread_setting=$3
    local case_directory=${test_root}/${label}
    local output_prefix=${case_directory}/reads
    local -a thread_options=()

    if [[ $thread_setting != default ]]; then
        thread_options=(-t "$thread_setting")
    fi

    mkdir -p "$case_directory"
    "$dwgsim_bin" \
        -z 13 \
        "${thread_options[@]}" \
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

run_case all-single 0 1
run_case all-threaded 0 4
run_case all-default 0 default
run_case bwa-threaded 1 4
run_case bfast-threaded 2 4

for suffix in bfast.fastq.gz bwa.read1.fastq.gz bwa.read2.fastq.gz; do
    cmp "${test_root}/all-single/reads.${suffix}" \
        "${test_root}/all-threaded/reads.${suffix}"
    cmp "${test_root}/all-single/reads.${suffix}" \
        "${test_root}/all-default/reads.${suffix}"
done

if command -v getconf >/dev/null 2>&1; then
    online_threads=$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)
    if [[ $online_threads =~ ^[1-9][0-9]*$ ]]; then
        help_output=$("$dwgsim_bin" -h 2>&1 || true)
        reported_default=$(printf '%s\n' "$help_output" |
            awk '/-t INT/ { value=$0; sub(/^.*\[/, "", value); sub(/\].*$/, "", value); print value; exit }')
        if [[ $reported_default != "$online_threads" ]]; then
            echo "test-bgzf: default thread count is $reported_default; expected $online_threads" >&2
            exit 1
        fi
    fi
fi

if "$dwgsim_bin" -t 0 -N 1 "$reference_fasta" \
    "${test_root}/invalid-thread-count" >/dev/null 2>&1; then
    echo "test-bgzf: -t 0 unexpectedly succeeded" >&2
    exit 1
fi

echo "BGZF FASTQ tests passed."
