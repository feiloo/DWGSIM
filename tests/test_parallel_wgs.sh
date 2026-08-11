#!/usr/bin/env bash

set -euo pipefail

dwgsim_bin=${DWGSIM_BIN:-./dwgsim}
samtools_bin=${SAMTOOLS_BIN:-./samtools/samtools}
bgzip_bin=${BGZIP_BIN:-./samtools/bgzip}
reference_fasta=${1:-samtools/examples/ex1.fa}
read_pairs=${PARALLEL_WGS_TEST_READ_PAIRS:-20000}

for dependency in "$dwgsim_bin" "$samtools_bin" "$bgzip_bin" gzip awk cmp \
    mktemp paste tail od tr; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        echo "test-parallel-wgs: required command not found: $dependency" >&2
        exit 1
    fi
done

if [[ ! -s $reference_fasta ]]; then
    echo "test-parallel-wgs: reference FASTA not found: $reference_fasta" >&2
    exit 1
fi
if [[ ! $read_pairs =~ ^[1-9][0-9]*$ ]]; then
    echo "test-parallel-wgs: PARALLEL_WGS_TEST_READ_PAIRS must be positive" >&2
    exit 2
fi

test_root=$(mktemp -d "${TMPDIR:-/tmp}/dwgsim-parallel-wgs.XXXXXX")
trap 'rm -rf "$test_root"' EXIT
# Replicate the first bundled contig under 160 names. With 20,000 pairs this
# creates 160 canonical tasks, so the -t 128 case exercises 128 real workers
# without making the test output large.
awk '
    /^>/ {
        if (found) exit
        found = 1
        next
    }
    found { sequence = sequence $0 }
    END {
        if (length(sequence) < 650) exit 1
        for (contig = 1; contig <= 160; contig++) {
            printf ">parallel%d\n", contig
            for (position = 1; position <= length(sequence); position += 60) {
                print substr(sequence, position, 60)
            }
        }
    }
' "$reference_fasta" >"$test_root/reference.fa"
"$samtools_bin" faidx "$test_root/reference.fa"

run_case() {
    local label=$1
    local threads=$2
    local compression_level=${3:-4}
    local prefix=${test_root}/${label}
    local -a thread_options=()

    if [[ $threads != automatic ]]; then
        thread_options=(-t "$threads")
    fi

    "$dwgsim_bin" \
        -z 13 \
        "${thread_options[@]}" \
        -l "$compression_level" \
        -N "$read_pairs" \
        -1 150 \
        -2 150 \
        -y 0 \
        -r 0 \
        -M 1 \
        -o 1 \
        "$test_root/reference.fa" \
        "$prefix" >"${prefix}.log" 2>&1

    grep -q 'deterministic WGS v1' "${prefix}.log"
    grep -q 'in 160 fixed tasks' "${prefix}.log"
    [[ -s ${prefix}.dwgsim.complete ]]
    [[ -s ${prefix}.bwa.read1.fastq.gz ]]
    [[ -s ${prefix}.bwa.read2.fastq.gz ]]
    gzip --test "${prefix}.bwa.read1.fastq.gz"
    gzip --test "${prefix}.bwa.read2.fastq.gz"
    "$bgzip_bin" -dc "${prefix}.bwa.read1.fastq.gz" >/dev/null
    "$bgzip_bin" -dc "${prefix}.bwa.read2.fastq.gz" >/dev/null
}

verify_canonical_eof() {
    local fastq=$1
    local expected=1f8b08040000000000ff0600424302001b0003000000000000000000
    local actual

    actual=$(tail -c 28 "$fastq" | od -An -tx1 | tr -d ' \n')
    if [[ $actual != "$expected" ]]; then
        echo "test-parallel-wgs: non-canonical BGZF EOF in $fastq" >&2
        exit 1
    fi
}

verify_records() {
    local prefix=$1
    local expected_lines=$((read_pairs * 4))
    local read1=${prefix}.bwa.read1.fastq.gz
    local read2=${prefix}.bwa.read2.fastq.gz

    gzip --decompress --stdout "$read1" |
        awk -v expected_lines="$expected_lines" '
            NR % 4 == 2 && length($0) != 150 { exit 1 }
            NR % 4 == 0 && length($0) != 150 { exit 1 }
            END { if (NR != expected_lines) exit 1 }
        '
    gzip --decompress --stdout "$read2" |
        awk -v expected_lines="$expected_lines" '
            NR % 4 == 2 && length($0) != 150 { exit 1 }
            NR % 4 == 0 && length($0) != 150 { exit 1 }
            END { if (NR != expected_lines) exit 1 }
        '
    paste \
        <(gzip --decompress --stdout "$read1" | awk 'NR % 4 == 1') \
        <(gzip --decompress --stdout "$read2" | awk 'NR % 4 == 1') |
        awk -F '\t' -v expected="$read_pairs" '
            {
                read1 = $1
                read2 = $2
                if (read1 !~ /\/1$/ || read2 !~ /\/2$/) exit 1
                sub(/\/1$/, "", read1)
                sub(/\/2$/, "", read2)
                if (read1 != read2) exit 1
                pairs++
            }
            END { if (pairs != expected) exit 1 }
        '
}

run_case threads-1 1
run_case threads-2 2
run_case threads-8 8
run_case threads-64 64
run_case threads-128 128
run_case threads-automatic automatic
run_case repeat-8 8

verify_records "$test_root/threads-1"
verify_canonical_eof "$test_root/threads-1.bwa.read1.fastq.gz"
verify_canonical_eof "$test_root/threads-1.bwa.read2.fastq.gz"

for mate in read1 read2; do
    baseline=${test_root}/threads-1.bwa.${mate}.fastq.gz
    cmp "$baseline" "${test_root}/threads-2.bwa.${mate}.fastq.gz"
    cmp "$baseline" "${test_root}/threads-8.bwa.${mate}.fastq.gz"
    cmp "$baseline" "${test_root}/threads-64.bwa.${mate}.fastq.gz"
    cmp "$baseline" "${test_root}/threads-128.bwa.${mate}.fastq.gz"
    cmp "$baseline" "${test_root}/threads-automatic.bwa.${mate}.fastq.gz"
    cmp "$baseline" "${test_root}/repeat-8.bwa.${mate}.fastq.gz"
done

for compression_level in {1..9}; do
    if [[ $compression_level -eq 4 ]]; then
        continue
    fi
    run_case "level-${compression_level}-threads-1" 1 "$compression_level"
    run_case "level-${compression_level}-threads-8" 8 "$compression_level"
    for mate in read1 read2; do
        level_single=${test_root}/level-${compression_level}-threads-1.bwa.${mate}.fastq.gz
        level_threaded=${test_root}/level-${compression_level}-threads-8.bwa.${mate}.fastq.gz
        cmp "$level_single" "$level_threaded"
        cmp \
            <(gzip --decompress --stdout "$test_root/threads-1.bwa.${mate}.fastq.gz") \
            <(gzip --decompress --stdout "$level_single")
    done
done

if find "$test_root" -type f -name '*.partial.*' -print -quit | grep -q .; then
    echo "test-parallel-wgs: staging file remained after success" >&2
    exit 1
fi

if "$dwgsim_bin" -l 0 -N 1 "$test_root/reference.fa" \
    "$test_root/invalid-level" >/dev/null 2>&1; then
    echo "test-parallel-wgs: -l 0 unexpectedly succeeded" >&2
    exit 1
fi
if "$dwgsim_bin" -l 10 -N 1 "$test_root/reference.fa" \
    "$test_root/invalid-level" >/dev/null 2>&1; then
    echo "test-parallel-wgs: -l 10 unexpectedly succeeded" >&2
    exit 1
fi

# This contig passes the loader's non-ACGT threshold but has no 150-base read
# with zero Ns. Generation therefore fails after staging files are opened.
awk '
    BEGIN {
        print ">unplaceable"
        sequence = ""
        for (base = 0; base < 100; base++) sequence = sequence "A"
        for (base = 0; base < 900; base++) sequence = sequence "N"
        for (position = 1; position <= length(sequence); position += 60) {
            print substr(sequence, position, 60)
        }
    }
' >"$test_root/unplaceable.fa"
"$samtools_bin" faidx "$test_root/unplaceable.fa"
if "$dwgsim_bin" -z 13 -t 2 -l 4 -N 1 -1 150 -2 150 -y 0 -r 0 \
    -M 1 -o 1 "$test_root/unplaceable.fa" \
    "$test_root/failed-generation" >/dev/null 2>&1; then
    echo "test-parallel-wgs: unplaceable generation unexpectedly succeeded" >&2
    exit 1
fi
if find "$test_root" -maxdepth 1 -type f \
    \( -name 'failed-generation*.fastq.gz' \
       -o -name 'failed-generation*.complete' \
       -o -name 'failed-generation*.partial.*' \) \
    -print -quit | grep -q .; then
    echo "test-parallel-wgs: failed generation published an output" >&2
    exit 1
fi

echo "Deterministic parallel WGS tests passed."
