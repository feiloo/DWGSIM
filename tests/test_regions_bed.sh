#!/usr/bin/env bash

set -euo pipefail

dwgsim_bin=${DWGSIM_BIN:-./dwgsim}

for dependency in "$dwgsim_bin" awk gzip grep; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        echo "test-bed: required command not found: $dependency" >&2
        exit 1
    fi
done

test_root=$(mktemp -d "${TMPDIR:-/tmp}/dwgsim-bed-test.XXXXXX")
trap 'rm -rf "$test_root"' EXIT

reference_fasta=${test_root}/reference.fa
awk 'BEGIN {
    print ">chr1"
    for (i = 0; i < 500; i++) printf "ACGT"
    print ""
    print ">chr2"
    for (i = 0; i < 500; i++) printf "TGCA"
    print ""
}' > "$reference_fasta"

valid_bed=${test_root}/valid.bed
printf '%s\n' \
    '# comment' \
    'track name=targets' \
    'browser position chr1:1-2000' \
    $'chr1\t0\t500\tfirst-target' \
    $'chr2\t100\t900\tsecond-target' > "$valid_bed"

valid_prefix=${test_root}/valid
"$dwgsim_bin" \
    -z 13 -N 12 -1 20 -2 20 -d 40 -s 0 -y 0 -r 0 -M 1 -o 2 \
    -x "$valid_bed" "$reference_fasta" "$valid_prefix" >/dev/null 2>&1

valid_fastq=${valid_prefix}.bfast.fastq.gz
gzip --test "$valid_fastq"
actual_lines=$(gzip --decompress --stdout "$valid_fastq" | awk 'END { print NR }')
if [[ $actual_lines -ne 96 ]]; then
    echo "test-bed: valid BED produced $actual_lines FASTQ lines; expected 96" >&2
    exit 1
fi

# The last BED target need not be the last FASTA contig. The full -N
# remainder must be assigned to the last targeted contig without retrying on
# a later, untargeted contig.
first_contig_bed=${test_root}/first-contig-only.bed
printf 'chr1\t0\t500\n' > "$first_contig_bed"
first_contig_prefix=${test_root}/first-contig-only
"$dwgsim_bin" \
    -z 13 -N 3 -1 20 -2 20 -d 40 -s 0 -y 0 -r 0 -M 1 -o 2 \
    -x "$first_contig_bed" "$reference_fasta" "$first_contig_prefix" >/dev/null 2>&1
first_contig_lines=$(gzip --decompress --stdout "${first_contig_prefix}.bfast.fastq.gz" |
    awk 'END { print NR }')
if [[ $first_contig_lines -ne 24 ]]; then
    echo "test-bed: last-target allocation produced $first_contig_lines FASTQ lines; expected 24" >&2
    exit 1
fi

expect_invalid_bed() {
    local label=$1
    local contents=$2
    local expected_message=$3
    local bed=${test_root}/${label}.bed
    local stderr_file=${test_root}/${label}.stderr
    local prefix=${test_root}/${label}

    printf '%b' "$contents" > "$bed"
    if "$dwgsim_bin" \
        -z 13 -N 1 -1 20 -2 20 -d 40 -s 0 -y 0 -r 0 -M 1 -o 2 \
        -x "$bed" "$reference_fasta" "$prefix" >/dev/null 2>"$stderr_file"; then
        echo "test-bed: invalid case '$label' unexpectedly succeeded" >&2
        exit 1
    fi
    if ! grep --fixed-strings --quiet -- "$expected_message" "$stderr_file"; then
        echo "test-bed: invalid case '$label' did not report '$expected_message'" >&2
        sed -n '1,20p' "$stderr_file" >&2
        exit 1
    fi
}

expect_invalid_bed empty '# comment only\n' 'regions BED contains no intervals'
expect_invalid_bed missing-field 'chr1\t0\n' 'has fewer than 3 fields'
expect_invalid_bed non-integer 'chr1\tzero\t10\n' 'non-integer or out-of-range coordinate'
expect_invalid_bed negative 'chr1\t-1\t10\n' 'non-integer or out-of-range coordinate'
expect_invalid_bed coordinate-overflow 'chr1\t0\t4294967296\n' 'non-integer or out-of-range coordinate'
expect_invalid_bed zero-length 'chr1\t10\t10\n' 'must satisfy start < end'
expect_invalid_bed reversed 'chr1\t20\t10\n' 'must satisfy start < end'
expect_invalid_bed out-of-range 'chr1\t0\t2001\n' 'end is outside reference contig'
expect_invalid_bed unknown-contig 'chr3\t0\t10\n' 'contig not found in reference'
expect_invalid_bed unsorted-start 'chr1\t20\t30\nchr1\t10\t15\n' 'not sorted by start'
expect_invalid_bed contig-order 'chr2\t0\t10\nchr1\t0\t10\n' 'not in reference contig order'

echo "Regions BED tests passed."
