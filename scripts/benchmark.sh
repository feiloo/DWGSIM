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
read_length_1=${BENCHMARK_READ_LENGTH_1:-150}
read_length_2=${BENCHMARK_READ_LENGTH_2:-150}
random_seed=${BENCHMARK_SEED:-13}
benchmark_threads=${BENCHMARK_THREADS:-$(online_cpu_count)}
compression_level=${BENCHMARK_COMPRESSION_LEVEL:-4}
estimate_coverage=${BENCHMARK_ESTIMATE_COVERAGE:-100}
measure_startup=${BENCHMARK_MEASURE_STARTUP:-0}

for dependency in "$dwgsim_bin" "$time_bin" date gzip awk wc tee tail mktemp; do
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
    "BENCHMARK_THREADS:$benchmark_threads" \
    "BENCHMARK_ESTIMATE_COVERAGE:$estimate_coverage"; do
    setting_name=${setting%%:*}
    setting_value=${setting#*:}
    if [[ ! $setting_value =~ ^[1-9][0-9]*$ ]]; then
        echo "benchmark: $setting_name must be a positive integer" >&2
        exit 2
    fi
done
if [[ ! $compression_level =~ ^[1-9]$ ]]; then
    echo "benchmark: BENCHMARK_COMPRESSION_LEVEL must be an integer from 1 to 9" >&2
    exit 2
fi

if [[ ! $random_seed =~ ^[0-9]+$ ]]; then
    echo "benchmark: BENCHMARK_SEED must be a non-negative integer" >&2
    exit 2
fi
if [[ $measure_startup != 0 && $measure_startup != 1 ]]; then
    echo "benchmark: BENCHMARK_MEASURE_STARTUP must be 0 or 1" >&2
    exit 2
fi

if [[ -s ${reference_fasta}.fai ]]; then
    reference_bases=$(awk '{ bases += $2 } END { printf "%.0f", bases }' \
        "${reference_fasta}.fai")
else
    reference_bases=$(awk '!/^>/ { gsub(/[[:space:]]/, ""); bases += length($0) }
        END { printf "%.0f", bases }' "$reference_fasta")
fi
if [[ ! $reference_bases =~ ^[1-9][0-9]*$ ]]; then
    echo "benchmark: could not determine a positive reference length" >&2
    exit 1
fi

mkdir -p "$output_directory"
output_prefix=${output_directory}/dwgsim-benchmark
resource_metrics_file=${output_directory}/resource-metrics.txt
report_file=${output_directory}/benchmark.txt
dwgsim_log=${output_directory}/dwgsim.log

echo "Generating $read_pairs paired ${read_length_1}+${read_length_2} bp reads with $benchmark_threads thread(s)..."
start_nanoseconds=$(date +%s%N)
if ! LC_ALL=C "$time_bin" \
    --format='user_seconds=%U\nsystem_seconds=%S\ncpu_percent=%P\nmax_rss_kib=%M' \
    --output="$resource_metrics_file" \
    "$dwgsim_bin" \
    -z "$random_seed" \
    -t "$benchmark_threads" \
    -l "$compression_level" \
    -N "$read_pairs" \
    -1 "$read_length_1" \
    -2 "$read_length_2" \
    -y 0 \
    -r 0 \
    -M 1 \
    -o 1 \
    "$reference_fasta" \
    "$output_prefix" >"$dwgsim_log" 2>&1; then
    echo "benchmark: DWGSIM failed; last log lines:" >&2
    tail -n 40 "$dwgsim_log" >&2
    exit 1
fi
end_nanoseconds=$(date +%s%N)

elapsed_seconds=$(awk \
    -v start="$start_nanoseconds" \
    -v end="$end_nanoseconds" \
    'BEGIN { printf "%.6f", (end - start) / 1000000000 }')

if ! awk -v elapsed="$elapsed_seconds" 'BEGIN { exit !(elapsed > 0) }'; then
    echo "benchmark: measured wall time was not positive" >&2
    exit 1
fi

startup_seconds=0
startup_pairs=0
if [[ $measure_startup == 1 ]]; then
    baseline_directory=$(mktemp -d "${TMPDIR:-/tmp}/dwgsim-benchmark-baseline.XXXXXX")
    trap 'rm -rf "$baseline_directory"' EXIT
    baseline_prefix=${baseline_directory}/baseline
    baseline_log=${baseline_directory}/dwgsim.log

    echo "Measuring one-pair full-reference startup baseline..."
    baseline_start_nanoseconds=$(date +%s%N)
    if ! "$dwgsim_bin" \
        -z "$random_seed" \
        -t "$benchmark_threads" \
        -l "$compression_level" \
        -N 1 \
        -1 "$read_length_1" \
        -2 "$read_length_2" \
        -y 0 \
        -r 0 \
        -M 1 \
        -o 1 \
        "$reference_fasta" \
        "$baseline_prefix" >"$baseline_log" 2>&1; then
        echo "benchmark: startup-baseline DWGSIM failed; last log lines:" >&2
        tail -n 40 "$baseline_log" >&2
        exit 1
    fi
    baseline_end_nanoseconds=$(date +%s%N)
    startup_seconds=$(awk \
        -v start="$baseline_start_nanoseconds" \
        -v end="$baseline_end_nanoseconds" \
        'BEGIN { printf "%.6f", (end - start) / 1000000000 }')
    startup_pairs=1
fi

read_1_fastq=${output_prefix}.bwa.read1.fastq.gz
read_2_fastq=${output_prefix}.bwa.read2.fastq.gz
expected_lines=$((read_pairs * 4))

for fastq in "$read_1_fastq" "$read_2_fastq"; do
    if [[ ! -s $fastq ]]; then
        echo "benchmark: expected FASTQ was not generated: $fastq" >&2
        exit 1
    fi
    if ! actual_lines=$(gzip --decompress --stdout "$fastq" |
        awk 'END { print NR }'); then
        echo "benchmark: compressed FASTQ integrity check failed: $fastq" >&2
        exit 1
    fi
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
    -v compression_level="$compression_level" \
    -v reference_bases="$reference_bases" \
    -v estimate_coverage="$estimate_coverage" \
    -v elapsed="$elapsed_seconds" \
    -v startup_seconds="$startup_seconds" \
    -v startup_pairs="$startup_pairs" \
    -v user="$user_seconds" \
    -v system_cpu="$system_seconds" \
    -v cpu="$cpu_percent" \
    -v max_rss_kib="$max_rss_kib" \
    -v output_bytes="$compressed_output_bytes" '
BEGIN {
    reads = pairs * 2
    pair_bases = length_1 + length_2
    bases = pairs * pair_bases
    output_mib = output_bytes / 1048576
    raw_pair_rate = pairs / elapsed
    measured_startup_pairs = startup_pairs
    measured_startup_seconds = startup_seconds
    throughput_pairs = pairs - startup_pairs
    throughput_seconds = elapsed - startup_seconds
    startup_adjustment_used = 1
    fixed_seconds = startup_seconds
    if (throughput_pairs <= 0 || throughput_seconds <= 0) {
        throughput_pairs = pairs
        throughput_seconds = elapsed
        startup_adjustment_used = 0
        fixed_seconds = 0
    }
    pair_rate = throughput_pairs / throughput_seconds
    estimated_pairs = int((reference_bases * estimate_coverage + pair_bases - 1) / pair_bases)
    estimated_seconds = fixed_seconds + estimated_pairs / pair_rate
    estimated_output_bytes = output_bytes * estimated_pairs / pairs
    rounded_seconds = int(estimated_seconds + 0.5)
    duration_days = int(rounded_seconds / 86400)
    duration_hours = int((rounded_seconds - duration_days * 86400) / 3600)
    duration_minutes = int((rounded_seconds - duration_days * 86400 - duration_hours * 3600) / 60)
    duration_seconds = rounded_seconds % 60

    print "DWGSIM performance benchmark"
    print "mode=wgs"
    print "reference=" reference
    printf "reference_bases=%.0f\n", reference_bases
    printf "read_pairs=%.0f\n", pairs
    printf "reads=%.0f\n", reads
    printf "read_length_1=%.0f\n", length_1
    printf "read_length_2=%.0f\n", length_2
    printf "bases=%.0f\n", bases
    printf "seed=%.0f\n", seed
    printf "threads=%.0f\n", threads
    printf "compression_level=%.0f\n", compression_level
    printf "elapsed_seconds=%.6f\n", elapsed
    printf "startup_baseline_pairs=%.0f\n", measured_startup_pairs
    printf "startup_baseline_seconds=%.6f\n", measured_startup_seconds
    printf "startup_adjustment_used=%.0f\n", startup_adjustment_used
    printf "startup_adjusted_elapsed_seconds=%.6f\n", throughput_seconds
    printf "user_seconds=%.2f\n", user
    printf "system_seconds=%.2f\n", system_cpu
    print "cpu_percent=" cpu
    printf "max_rss_kib=%.0f\n", max_rss_kib
    printf "max_rss_mib=%.2f\n", max_rss_kib / 1024
    printf "compressed_output_bytes=%.0f\n", output_bytes
    printf "compressed_output_mib=%.2f\n", output_mib
    printf "read_pairs_per_second=%.2f\n", raw_pair_rate
    printf "reads_per_second=%.2f\n", reads / elapsed
    printf "bases_per_second=%.2f\n", bases / elapsed
    printf "megabases_per_second=%.2f\n", bases / elapsed / 1000000
    printf "compressed_mib_per_second=%.2f\n", output_mib / elapsed
    printf "startup_adjusted_read_pairs_per_second=%.2f\n", pair_rate
    printf "startup_adjusted_reads_per_second=%.2f\n", pair_rate * 2
    printf "startup_adjusted_bases_per_second=%.2f\n", pair_rate * pair_bases
    printf "startup_adjusted_megabases_per_second=%.2f\n", pair_rate * pair_bases / 1000000
    printf "estimate_coverage=%.0f\n", estimate_coverage
    printf "estimated_wgs_read_pairs=%.0f\n", estimated_pairs
    printf "estimated_wgs_seconds=%.2f\n", estimated_seconds
    printf "estimated_wgs_hours=%.2f\n", estimated_seconds / 3600
    printf "estimated_wgs_days=%.3f\n", estimated_seconds / 86400
    printf "estimated_wgs_duration=%dd %02dh %02dm %02ds\n", duration_days, duration_hours, duration_minutes, duration_seconds
    printf "estimated_compressed_output_bytes=%.0f\n", estimated_output_bytes
    printf "estimated_compressed_output_gib=%.2f\n", estimated_output_bytes / 1073741824
}' | tee "$report_file"

echo "Benchmark report: $report_file"
