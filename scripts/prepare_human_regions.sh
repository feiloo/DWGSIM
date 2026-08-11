#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 7 ]]; then
    echo "Usage: $0 <reference.fa> <samtools> <assembly-report> <wes.bed.gz> <blacklist.bed.gz> <wes-output.bed> <filtered-wgs-output.bed>" >&2
    exit 2
fi

reference_fasta=$1
samtools_bin=$2
assembly_report=$3
wes_source=$4
blacklist_source=$5
wes_output=$6
filtered_wgs_output=$7
wes_padding=${HUMAN_WES_PADDING:-250}
reference_fai=${reference_fasta}.fai

for dependency in "$samtools_bin" awk gzip; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        echo "prepare-human-regions: required command not found: $dependency" >&2
        exit 1
    fi
done

for input in "$reference_fasta" "$assembly_report" "$wes_source" "$blacklist_source"; do
    if [[ ! -s $input ]]; then
        echo "prepare-human-regions: input not found: $input" >&2
        exit 1
    fi
done
if [[ ! $wes_padding =~ ^[0-9]+$ ]]; then
    echo "prepare-human-regions: HUMAN_WES_PADDING must be a non-negative integer" >&2
    exit 2
fi

gzip --test "$wes_source"
gzip --test "$blacklist_source"

if [[ ! -s $reference_fai || $reference_fasta -nt $reference_fai ]]; then
    echo "Indexing the GRCh38 FASTA for exact contig order and lengths..."
    "$samtools_bin" faidx "$reference_fasta"
fi

mkdir -p "$(dirname "$wes_output")" "$(dirname "$filtered_wgs_output")"
work_directory=$(mktemp -d "${TMPDIR:-/tmp}/dwgsim-human-regions.XXXXXX")
trap 'rm -rf "$work_directory"' EXIT

wes_uncompressed=${work_directory}/wes-source.bed
blacklist_uncompressed=${work_directory}/blacklist-source.bed
wes_temporary=${work_directory}/wes-derived.bed
filtered_temporary=${work_directory}/filtered-wgs-derived.bed
gzip --decompress --stdout "$wes_source" > "$wes_uncompressed"
gzip --decompress --stdout "$blacklist_source" > "$blacklist_uncompressed"

# Translate UCSC chromosome names to the RefSeq accessions used by the NCBI
# FASTA, pad capture targets, and merge overlapping or touching intervals.
awk -F '\t' -v OFS='\t' -v padding="$wes_padding" '
BEGIN {
    print "# Source: GIAB v3.6 Functional/GRCh38_refseq_cds.bed.gz"
    print "# Transform: UCSC names to GRCh38.p14 RefSeq accessions; padded " padding " bp; merged"
}
FILENAME == ARGV[1] {
    if ($0 !~ /^#/ && NF >= 10 && $7 != "na" && $10 != "na") {
        ucsc = $10
        sub(/\r$/, "", ucsc)
        alias[ucsc] = $7
        contig_lengths[$7] = $9 + 0
    }
    next
}
$0 ~ /^[[:space:]]*$/ || $0 ~ /^[[:space:]]*#/ ||
$1 == "track" || $1 == "browser" { next }
{
    if (NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 >= $3) {
        print "prepare-human-regions: invalid WES BED record: " $0 > "/dev/stderr"
        failed = 1
        exit 1
    }
    chrom = alias[$1]
    if (chrom == "") {
        print "prepare-human-regions: WES contig has no RefSeq alias: " $1 > "/dev/stderr"
        failed = 1
        exit 1
    }
    raw_start = $2 + 0
    raw_end = $3 + 0
    if (source_count[chrom] > 0 && raw_start < source_last_start[chrom]) {
        print "prepare-human-regions: WES BED is not sorted within " $1 > "/dev/stderr"
        failed = 1
        exit 1
    }
    source_count[chrom]++
    source_last_start[chrom] = raw_start
    start = raw_start - padding
    if (start < 0) start = 0
    end = raw_end + padding
    if (end > contig_lengths[chrom]) end = contig_lengths[chrom]
    if (start >= end) {
        print "prepare-human-regions: invalid padded WES interval: " $0 > "/dev/stderr"
        failed = 1
        exit 1
    }

    if (have_interval && chrom == previous_chrom && start <= previous_end) {
        if (end > previous_end) previous_end = end
    }
    else {
        flush_interval()
        previous_chrom = chrom
        previous_start = start
        previous_end = end
        have_interval = 1
    }
}
END {
    if (!failed) flush_interval()
}
function flush_interval() {
    if (have_interval) print previous_chrom, previous_start, previous_end
}
' "$assembly_report" "$wes_uncompressed" > "$wes_temporary"

# Store the small ENCODE blacklist by RefSeq accession, then walk the FASTA
# index in its exact order and emit the complement as an inclusion BED.
awk -F '\t' -v OFS='\t' '
BEGIN {
    print "# Source exclusions: ENCODE GRCh38 blacklist ENCFF356LFX"
    print "# Transform: UCSC names to GRCh38.p14 RefSeq accessions; complement across every FASTA contig"
}
FILENAME == ARGV[1] {
    if ($0 !~ /^#/ && NF >= 10 && $7 != "na" && $10 != "na") {
        ucsc = $10
        sub(/\r$/, "", ucsc)
        alias[ucsc] = $7
    }
    next
}
FILENAME == ARGV[2] {
    if ($0 ~ /^[[:space:]]*$/ || $0 ~ /^[[:space:]]*#/ ||
        $1 == "track" || $1 == "browser") next
    if (NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 >= $3) {
        print "prepare-human-regions: invalid blacklist BED record: " $0 > "/dev/stderr"
        failed = 1
        exit 1
    }
    chrom = alias[$1]
    if (chrom == "") {
        print "prepare-human-regions: blacklist contig has no RefSeq alias: " $1 > "/dev/stderr"
        failed = 1
        exit 1
    }
    if (count[chrom] > 0 && ($2 + 0) < last_start[chrom]) {
        print "prepare-human-regions: blacklist BED is not sorted within " $1 > "/dev/stderr"
        failed = 1
        exit 1
    }
    i = ++count[chrom]
    starts[chrom, i] = $2 + 0
    ends[chrom, i] = $3 + 0
    last_start[chrom] = $2 + 0
    next
}
FILENAME == ARGV[3] {
    chrom = $1
    contig_length = $2 + 0
    cursor = 0
    for (i = 1; i <= count[chrom]; i++) {
        start = starts[chrom, i]
        end = ends[chrom, i]
        if (start < 0 || contig_length < end) {
            print "prepare-human-regions: blacklist interval is outside " chrom > "/dev/stderr"
            failed = 1
            exit 1
        }
        if (end <= cursor) continue
        if (start < cursor) start = cursor
        if (cursor < start) print chrom, cursor, start
        cursor = end
    }
    if (cursor < contig_length) print chrom, cursor, contig_length
}
END {
    if (failed) exit 1
}
' "$assembly_report" "$blacklist_uncompressed" "$reference_fai" > "$filtered_temporary"

validate_bed() {
    local bed=$1
    local label=$2
    awk -F '\t' -v label="$label" '
    FILENAME == ARGV[1] {
        order[$1] = ++contig_count
        contig_lengths[$1] = $2 + 0
        next
    }
    $0 ~ /^[[:space:]]*$/ || $0 ~ /^[[:space:]]*#/ ||
    $1 == "track" || $1 == "browser" { next }
    {
        if (NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) fail("malformed interval")
        current_order = order[$1]
        if (current_order == 0) fail("contig absent from FASTA index")
        if (current_order < previous_order ||
            (current_order == previous_order && $2 < previous_start)) {
            fail("intervals are not in FASTA coordinate order")
        }
        if ($2 >= $3 || $3 > contig_lengths[$1]) fail("interval is empty or outside its contig")
        previous_order = current_order
        previous_start = $2 + 0
        records++
        bases += ($3 - $2)
    }
    END {
        if (records == 0) fail("contains no intervals")
        printf "Prepared %s: %d intervals, %.0f bases\n", label, records, bases > "/dev/stderr"
    }
    function fail(message) {
        print "prepare-human-regions: " label ": " message " at line " FNR > "/dev/stderr"
        exit 1
    }
    ' "$reference_fai" "$bed"
}

validate_bed "$wes_temporary" "WES RefSeq CDS targets"
validate_bed "$filtered_temporary" "WGS excluding ENCODE ENCFF356LFX"

mv "$wes_temporary" "$wes_output"
mv "$filtered_temporary" "$filtered_wgs_output"
echo "Human-region BEDs ready in $(dirname "$wes_output")"
