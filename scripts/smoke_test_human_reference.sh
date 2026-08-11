#!/usr/bin/env bash

set -euo pipefail

online_cpu_count() {
    local detected

    if command -v getconf >/dev/null 2>&1; then
        detected=$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)
        if [[ $detected =~ ^[1-9][0-9]*$ ]]; then
            printf '%s\n' "$detected"
            return
        fi
    fi
    printf '1\n'
}

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 <reference.fasta> <output-directory> <wgs|wes|wgs-filtered> [regions.bed]" >&2
    exit 2
fi

reference_fasta=$1
output_directory=$2
profile=$3
regions_bed=${4:-}
dwgsim_bin=${DWGSIM_BIN:-./dwgsim}
read_pairs=${HUMAN_SMOKE_READ_PAIRS:-100}
random_seed=${HUMAN_SMOKE_SEED:-13}
threads=${HUMAN_SMOKE_THREADS:-$(online_cpu_count)}

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
if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo "test-human-reference: HUMAN_SMOKE_THREADS must be a positive integer" >&2
    exit 2
fi

bed_options=()
case $profile in
    wgs)
        description="whole-genome"
        if [[ -n $regions_bed ]]; then
            echo "test-human-reference: wgs does not accept a BED argument" >&2
            exit 2
        fi
        ;;
    wes)
        description="whole-exome (GIAB RefSeq CDS plus configured padding)"
        bed_options=(-x "$regions_bed")
        ;;
    wgs-filtered)
        description="whole-genome excluding ENCODE ENCFF356LFX"
        bed_options=(-x "$regions_bed")
        ;;
    *)
        echo "test-human-reference: unknown profile: $profile" >&2
        exit 2
        ;;
esac

if [[ $profile != wgs && ! -s $regions_bed ]]; then
    echo "test-human-reference: regions BED not found for $profile: $regions_bed" >&2
    exit 1
fi

case_directory=${output_directory}/${profile}
mkdir -p "$case_directory"
output_prefix=${case_directory}/grch38-p14-${profile}-smoke
dwgsim_log=${case_directory}/dwgsim.log

echo "Simulating $read_pairs read pairs for $description against complete GRCh38.p14..."
if ! "$dwgsim_bin" \
    -z "$random_seed" \
    -t "$threads" \
    -N "$read_pairs" \
    -1 150 \
    -2 150 \
    -d 350 \
    -s 50 \
    -y 0 \
    -r 0 \
    -M 1 \
    "${bed_options[@]}" \
    "$reference_fasta" \
    "$output_prefix" >"$dwgsim_log" 2>&1; then
    echo "test-human-reference: DWGSIM failed for $profile; last log lines:" >&2
    tail -n 40 "$dwgsim_log" >&2
    exit 1
fi

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

random_headers=$(gzip --decompress --stdout "$bfast_fastq" |
    awk 'NR % 4 == 1 && /^@rand_/ { count++ } END { print count + 0 }')
if [[ $random_headers -ne 0 ]]; then
    echo "test-human-reference: $profile generated $random_headers off-reference random reads" >&2
    exit 1
fi

echo "GRCh38 $profile smoke test passed. Output: $case_directory"
