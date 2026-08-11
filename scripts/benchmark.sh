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

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <reference.fasta> <output-directory>" >&2
    exit 2
fi

reference_fasta=$1
output_directory=$2
dwgsim_bin=${DWGSIM_BIN:-./dwgsim}
time_bin=${TIME_BIN:-/usr/bin/time}
read_pairs=${BENCHMARK_READ_PAIRS:-250000}
read_length_1=${BENCHMARK_READ_LENGTH_1:-100}
read_length_2=${BENCHMARK_READ_LENGTH_2:-100}
random_seed=${BENCHMARK_SEED:-13}
benchmark_threads=${BENCHMARK_THREADS:-$(online_cpu_count)}

for dependency in "$dwgsim_bin" "$time_bin" date gzip awk wc tee; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        echo "benchmark: required command not found: $dependency" >&2
        exit 1
    fi
done

if [[ ! -s $reference_fasta ]]; then
    echo "benchmark: reference FASTA not found: $reference_fasta" >&2
    exit 1
fi

for setting in \
    "BENCHMARK_READ_PAIRS:$read_pairs" \
    "BENCHMARK_READ_LENGTH_1:$read_length_1" \
    "BENCHMARK_READ_LENGTH_2:$read_length_2" \
    "BENCHMARK_THREADS:$benchmark_threads"; do
    setting_name=${setting%%:*}
    setting_value=${setting#*:}
    if [[ ! $setting_value =~ ^[1-9][0-9]*$ ]]; then
        echo "benchmark: $setting_name must be a positive integer" >&2
        exit 2
    fi
done

if [[ ! $random_seed =~ ^[0-9]+$ ]]; then
    echo "benchmark: BENCHMARK_SEED must be a non-negative integer" >&2
    exit 2
fi

mkdir -p "$output_directory"
output_prefix=${output_directory}/dwgsim-benchmark
resource_metrics_file=${output_directory}/resource-metrics.txt
report_file=${output_directory}/benchmark.txt

echo "Generating $read_pairs paired ${read_length_1}+${read_length_2} bp reads with $benchmark_threads thread(s)..."
start_nanoseconds=$(date +%s%N)
LC_ALL=C "$time_bin" \
    --format='user_seconds=%U\nsystem_seconds=%S\ncpu_percent=%P\nmax_rss_kib=%M' \
    --output="$resource_metrics_file" \
    "$dwgsim_bin" \
    -z "$random_seed" \
    -t "$benchmark_threads" \
    -N "$read_pairs" \
    -1 "$read_length_1" \
    -2 "$read_length_2" \
    -r 0 \
    -M 1 \
    -o 1 \
    "$reference_fasta" \
    "$output_prefix"
end_nanoseconds=$(date +%s%N)

elapsed_seconds=$(awk \
    -v start="$start_nanoseconds" \
    -v end="$end_nanoseconds" \
    'BEGIN { printf "%.6f", (end - start) / 1000000000 }')

if ! awk -v elapsed="$elapsed_seconds" 'BEGIN { exit !(elapsed > 0) }'; then
    echo "benchmark: measured wall time was not positive" >&2
    exit 1
fi

read_1_fastq=${output_prefix}.bwa.read1.fastq.gz
read_2_fastq=${output_prefix}.bwa.read2.fastq.gz
expected_lines=$((read_pairs * 4))

for fastq in "$read_1_fastq" "$read_2_fastq"; do
    if [[ ! -s $fastq ]]; then
        echo "benchmark: expected FASTQ was not generated: $fastq" >&2
        exit 1
    fi
    gzip --test "$fastq"
    actual_lines=$(gzip --decompress --stdout "$fastq" | awk 'END { print NR }')
    if [[ $actual_lines -ne $expected_lines ]]; then
        echo "benchmark: $fastq has $actual_lines lines; expected $expected_lines" >&2
        exit 1
    fi
done

read_metric() {
    local metric_name=$1
    awk -F= -v metric_name="$metric_name" \
        '$1 == metric_name { print substr($0, length(metric_name) + 2); exit }' \
        "$resource_metrics_file"
}

user_seconds=$(read_metric user_seconds)
system_seconds=$(read_metric system_seconds)
cpu_percent=$(read_metric cpu_percent)
max_rss_kib=$(read_metric max_rss_kib)
compressed_output_bytes=$(( $(wc -c < "$read_1_fastq") + $(wc -c < "$read_2_fastq") ))

awk \
    -v reference="$reference_fasta" \
    -v pairs="$read_pairs" \
    -v length_1="$read_length_1" \
    -v length_2="$read_length_2" \
    -v seed="$random_seed" \
    -v threads="$benchmark_threads" \
    -v elapsed="$elapsed_seconds" \
    -v user="$user_seconds" \
    -v system_cpu="$system_seconds" \
    -v cpu="$cpu_percent" \
    -v max_rss_kib="$max_rss_kib" \
    -v output_bytes="$compressed_output_bytes" '
BEGIN {
    reads = pairs * 2
    bases = pairs * (length_1 + length_2)
    output_mib = output_bytes / 1048576

    print "DWGSIM performance benchmark"
    print "reference=" reference
    printf "read_pairs=%.0f\n", pairs
    printf "reads=%.0f\n", reads
    printf "read_length_1=%.0f\n", length_1
    printf "read_length_2=%.0f\n", length_2
    printf "bases=%.0f\n", bases
    printf "seed=%.0f\n", seed
    printf "threads=%.0f\n", threads
    printf "elapsed_seconds=%.6f\n", elapsed
    printf "user_seconds=%.2f\n", user
    printf "system_seconds=%.2f\n", system_cpu
    print "cpu_percent=" cpu
    printf "max_rss_kib=%.0f\n", max_rss_kib
    printf "max_rss_mib=%.2f\n", max_rss_kib / 1024
    printf "compressed_output_bytes=%.0f\n", output_bytes
    printf "compressed_output_mib=%.2f\n", output_mib
    printf "read_pairs_per_second=%.2f\n", pairs / elapsed
    printf "reads_per_second=%.2f\n", reads / elapsed
    printf "bases_per_second=%.2f\n", bases / elapsed
    printf "megabases_per_second=%.2f\n", bases / elapsed / 1000000
    printf "compressed_mib_per_second=%.2f\n", output_mib / elapsed
}' | tee "$report_file"

echo "Benchmark report: $report_file"
