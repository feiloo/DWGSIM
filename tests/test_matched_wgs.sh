#!/usr/bin/env bash

set -euo pipefail

dwgsim_bin=${DWGSIM_BIN:-./dwgsim}
samtools_bin=${SAMTOOLS_BIN:-./samtools/samtools}
bgzip_bin=${BGZIP_BIN:-./samtools/bgzip}
reference_fasta=${1:-samtools/examples/ex1.fa}
read_pairs=${MATCHED_WGS_TEST_READ_PAIRS:-10000}

for dependency in "$dwgsim_bin" "$samtools_bin" "$bgzip_bin" gzip awk cmp \
    cp find mktemp paste tail od tr grep; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        echo "test-matched-wgs: required command not found: $dependency" >&2
        exit 1
    fi
done
if [[ ! -s $reference_fasta ]]; then
    echo "test-matched-wgs: reference FASTA not found: $reference_fasta" >&2
    exit 1
fi
if [[ ! $read_pairs =~ ^[1-9][0-9]*$ ]]; then
    echo "test-matched-wgs: MATCHED_WGS_TEST_READ_PAIRS must be positive" >&2
    exit 2
fi
if ((read_pairs < 2000)); then
    echo "test-matched-wgs: MATCHED_WGS_TEST_READ_PAIRS must be at least 2000 for the VAF check" >&2
    exit 2
fi

test_root=$(mktemp -d "${TMPDIR:-/tmp}/dwgsim-matched-wgs.XXXXXX")
trap 'rm -rf "$test_root"' EXIT

# A canonical 160-contig manifest creates 160 real tasks even for this small
# test, so the 128-worker case is a scheduling test rather than a label-only
# thread-count check.
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
            printf ">matched%d\n", contig
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
    local prefix=$test_root/$label
    local -a thread_options=()

    if [[ $threads != automatic ]]; then
        thread_options=(-t "$threads")
    fi

    "$dwgsim_bin" \
        --matched \
        --somatic-rate 0.01 \
        --tumor-vaf 0.25 \
        -z 29 \
        "${thread_options[@]}" \
        -l 4 \
        -N "$read_pairs" \
        -1 150 \
        -2 150 \
        -e 0.001 \
        -E 0.001 \
        -r 0.01 \
        -R 0.25 \
        "$test_root/reference.fa" \
        "$prefix" >"${prefix}.log" 2>&1

    grep -q 'deterministic matched v1' "${prefix}.log"
    grep -q 'in 160 fixed tasks' "${prefix}.log"
    [[ -s ${prefix}.matched.complete ]]
    [[ -s ${prefix}.germline.vcf ]]
    [[ -s ${prefix}.somatic.vcf ]]
}

run_tumor_only_case() {
    local label=$1
    local threads=$2
    local prefix=$test_root/$label
    local -a thread_options=()

    if [[ $threads != automatic ]]; then
        thread_options=(-t "$threads")
    fi

    "$dwgsim_bin" \
        --tumor-only \
        --somatic-rate 0.01 \
        --tumor-vaf 0.25 \
        -z 29 \
        "${thread_options[@]}" \
        -l 4 \
        -N "$read_pairs" \
        -1 150 \
        -2 150 \
        -e 0.001 \
        -E 0.001 \
        -r 0.01 \
        -R 0.25 \
        "$test_root/reference.fa" \
        "$prefix" >"${prefix}.log" 2>&1

    grep -q 'deterministic tumor-only v1' "${prefix}.log"
    grep -q 'in 160 fixed tasks' "${prefix}.log"
    [[ -s ${prefix}.tumor-only.complete ]]
    [[ -s ${prefix}.germline.vcf ]]
    [[ -s ${prefix}.somatic.vcf ]]
    [[ ! -e ${prefix}.normal.bwa.read1.fastq.gz ]]
    [[ ! -e ${prefix}.normal.bwa.read2.fastq.gz ]]
    [[ ! -e ${prefix}.matched.complete ]]
}

verify_eof() {
    local fastq=$1
    local expected=1f8b08040000000000ff0600424302001b0003000000000000000000
    local actual
    actual=$(tail -c 28 "$fastq" | od -An -tx1 | tr -d ' \n')
    [[ $actual == "$expected" ]]
}

verify_sample() {
    local prefix=$1
    local sample=$2
    local expected_pairs=$3
    local read1=${prefix}.${sample}.bwa.read1.fastq.gz
    local read2=${prefix}.${sample}.bwa.read2.fastq.gz
    local expected_lines=$((expected_pairs * 4))

    [[ -s $read1 && -s $read2 ]]
    gzip --test "$read1"
    gzip --test "$read2"
    "$bgzip_bin" -dc "$read1" >/dev/null
    "$bgzip_bin" -dc "$read2" >/dev/null
    verify_eof "$read1"
    verify_eof "$read2"
    gzip -dc "$read1" |
        awk -v expected="$expected_lines" -v sample="$sample" '
            NR == 1 && $0 !~ ("^@" sample "_") { exit 1 }
            NR % 4 == 2 && length($0) != 150 { exit 1 }
            NR % 4 == 0 && length($0) != 150 { exit 1 }
            END { if (NR != expected) exit 1 }
        '
    gzip -dc "$read2" |
        awk -v expected="$expected_lines" -v sample="$sample" '
            NR == 1 && $0 !~ ("^@" sample "_") { exit 1 }
            NR % 4 == 2 && length($0) != 150 { exit 1 }
            NR % 4 == 0 && length($0) != 150 { exit 1 }
            END { if (NR != expected) exit 1 }
        '
    paste \
        <(gzip -dc "$read1" | awk 'NR % 4 == 1 { sub(/\/1$/, ""); print }') \
        <(gzip -dc "$read2" | awk 'NR % 4 == 1 { sub(/\/2$/, ""); print }') |
        awk -v expected="$expected_pairs" '
            $1 != $2 { exit 1 }
            END { if (NR != expected) exit 1 }
        '
}

verify_truth() {
    local prefix=$1

    awk -F '\t' '
        /^#/ { next }
        NF != 11 || $7 != "PASS" || $8 !~ /SCOPE=GERMLINE/ { exit 1 }
        $10 != $11 { exit 1 }
        $10 !~ /^(1\|0|0\|1|1\|1)$/ { exit 1 }
        { count++ }
        END { if (count == 0) exit 1 }
    ' "${prefix}.germline.vcf"
    awk -F '\t' '
        /^#/ { next }
        NF != 11 || $7 != "PASS" || $8 !~ /^SOMATIC;/ { exit 1 }
        $8 !~ /SCOPE=SOMATIC/ || $8 !~ /EXPECTED_VAF=0.25/ { exit 1 }
        $10 != "0|0" || $11 !~ /^(1\|0|0\|1)$/ { exit 1 }
        { count++ }
        END { if (count == 0) exit 1 }
    ' "${prefix}.somatic.vcf"
    grep -q 'TYPE=SNV' "${prefix}.germline.vcf"
    grep -q 'TYPE=INS' "${prefix}.germline.vcf"
    grep -q 'TYPE=DEL' "${prefix}.germline.vcf"
}

run_case one 1
run_case eight 8
run_case eight_repeat 8
run_case one_twenty_eight 128
run_case automatic automatic
run_tumor_only_case tumor_one 1
run_tumor_only_case tumor_eight 8
run_tumor_only_case tumor_eight_repeat 8
run_tumor_only_case tumor_one_twenty_eight 128
run_tumor_only_case tumor_automatic automatic

for prefix in one eight eight_repeat one_twenty_eight automatic; do
    verify_sample "$test_root/$prefix" normal "$read_pairs"
    verify_sample "$test_root/$prefix" tumor "$read_pairs"
    verify_truth "$test_root/$prefix"
done

for prefix in tumor_one tumor_eight tumor_eight_repeat \
    tumor_one_twenty_eight tumor_automatic; do
    verify_sample "$test_root/$prefix" tumor "$read_pairs"
    verify_truth "$test_root/$prefix"
    grep -q '^format=dwgsim-deterministic-tumor-only-v1$' \
        "$test_root/$prefix.tumor-only.complete"
    grep -q '^normal_read_pairs=0$' \
        "$test_root/$prefix.tumor-only.complete"
    grep -q "^tumor_read_pairs=$read_pairs$" \
        "$test_root/$prefix.tumor-only.complete"
    if grep -q '^normal_read[12]_bytes=' \
        "$test_root/$prefix.tumor-only.complete"; then
        echo "test-matched-wgs: tumor-only manifest contains normal FASTQs" >&2
        exit 1
    fi
done

for suffix in \
    normal.bwa.read1.fastq.gz normal.bwa.read2.fastq.gz \
    tumor.bwa.read1.fastq.gz tumor.bwa.read2.fastq.gz \
    germline.vcf somatic.vcf; do
    cmp "$test_root/one.$suffix" "$test_root/eight.$suffix"
    cmp "$test_root/one.$suffix" "$test_root/eight_repeat.$suffix"
    cmp "$test_root/one.$suffix" "$test_root/one_twenty_eight.$suffix"
    cmp "$test_root/one.$suffix" "$test_root/automatic.$suffix"
done
for suffix in \
    tumor.bwa.read1.fastq.gz tumor.bwa.read2.fastq.gz \
    germline.vcf somatic.vcf; do
    cmp "$test_root/tumor_one.$suffix" "$test_root/tumor_eight.$suffix"
    cmp "$test_root/tumor_one.$suffix" \
        "$test_root/tumor_eight_repeat.$suffix"
    cmp "$test_root/tumor_one.$suffix" \
        "$test_root/tumor_one_twenty_eight.$suffix"
    cmp "$test_root/tumor_one.$suffix" \
        "$test_root/tumor_automatic.$suffix"
    cmp "$test_root/one.$suffix" "$test_root/tumor_one.$suffix"
done
if cmp -s "$test_root/one.normal.bwa.read1.fastq.gz" \
          "$test_root/one.tumor.bwa.read1.fastq.gz"; then
    echo "test-matched-wgs: normal and tumor libraries unexpectedly match" >&2
    exit 1
fi

# Independently allocated library sizes exercise tasks with an empty sample
# tail while preserving four-stream commit order.
"$dwgsim_bin" --matched --normal-pairs 33 --tumor-pairs 777 \
    --somatic-rate 0.01 --tumor-vaf 0.25 -r 0.01 -R 0 \
    -e 0 -E 0 -z 31 -t 8 "$test_root/reference.fa" \
    "$test_root/unequal" >"$test_root/unequal.log" 2>&1
verify_sample "$test_root/unequal" normal 33
verify_sample "$test_root/unequal" tumor 777
grep -q '^normal_read_pairs=33$' "$test_root/unequal.matched.complete"
grep -q '^tumor_read_pairs=777$' "$test_root/unequal.matched.complete"

# A sample-specific count is sufficient in tumor-only mode and must not create
# placeholder normal streams.
"$dwgsim_bin" --tumor-only --tumor-pairs 777 \
    --somatic-rate 0.01 --tumor-vaf 0.25 -r 0.01 -R 0 \
    -e 0 -E 0 -z 31 -t 8 "$test_root/reference.fa" \
    "$test_root/tumor-explicit" >"$test_root/tumor-explicit.log" 2>&1
verify_sample "$test_root/tumor-explicit" tumor 777
grep -q '^normal_read_pairs=0$' \
    "$test_root/tumor-explicit.tumor-only.complete"
grep -q '^tumor_read_pairs=777$' \
    "$test_root/tumor-explicit.tumor-only.complete"
[[ ! -e $test_root/tumor-explicit.normal.bwa.read1.fastq.gz ]]
[[ ! -e $test_root/tumor-explicit.normal.bwa.read2.fastq.gz ]]

# With germline events and sequencing errors disabled, somatic SNV alleles
# must be absent from normal and converge on the configured tumor VAF.
"$dwgsim_bin" --matched -N "$read_pairs" --somatic-rate 0.01 \
    --tumor-vaf 0.25 -r 0 -R 0 -e 0 -E 0 -Q 0 -z 37 -t 8 \
    "$test_root/reference.fa" "$test_root/vaf" \
    >"$test_root/vaf.log" 2>&1
"$dwgsim_bin" --tumor-only -N "$read_pairs" --somatic-rate 0.01 \
    --tumor-vaf 0.25 -r 0 -R 0 -e 0 -E 0 -Q 0 -z 37 -t 8 \
    "$test_root/reference.fa" "$test_root/vaf-tumor-only" \
    >"$test_root/vaf-tumor-only.log" 2>&1

cmp "$test_root/vaf.tumor.bwa.read1.fastq.gz" \
    "$test_root/vaf-tumor-only.tumor.bwa.read1.fastq.gz"
cmp "$test_root/vaf.tumor.bwa.read2.fastq.gz" \
    "$test_root/vaf-tumor-only.tumor.bwa.read2.fastq.gz"
cmp "$test_root/vaf.germline.vcf" \
    "$test_root/vaf-tumor-only.germline.vcf"
cmp "$test_root/vaf.somatic.vcf" \
    "$test_root/vaf-tumor-only.somatic.vcf"

measure_somatic_alleles() {
    local prefix=$1
    local sample=$2
    local fastq=${prefix}.${sample}.bwa.read1.fastq.gz

    awk -F '\t' -v sample="$sample" '
        function complement(base) {
            if (base == "A") return "T"
            if (base == "C") return "G"
            if (base == "G") return "C"
            return "A"
        }
        FNR == NR {
            if ($0 !~ /^#/ && $8 ~ /TYPE=SNV/) {
                key = $1 SUBSEP $2
                reference[key] = $4
                alternate[key] = $5
            }
            next
        }
        FNR % 4 == 1 {
            split($0, header, "_")
            contig = header[2]
            start = header[3] + 0
            strand = header[5] + 0
            next
        }
        FNR % 4 == 2 {
            sequence = $0
            for (key in alternate) {
                split(key, locus, SUBSEP)
                position = locus[2] + 0
                if (locus[1] != contig ||
                    position < start || start + 149 < position) continue
                cycle = strand ? start + 149 - position + 1 : position - start + 1
                expected_ref = strand ? complement(reference[key]) : reference[key]
                expected_alt = strand ? complement(alternate[key]) : alternate[key]
                observed = substr(sequence, cycle, 1)
                if (observed == expected_alt) alt_count++
                else if (observed == expected_ref) ref_count++
                else other_count++
            }
        }
        END {
            total = alt_count + ref_count
            if (other_count != 0 || total < 1000) exit 1
            if (sample == "normal" && alt_count != 0) exit 1
            if (sample == "tumor") {
                ratio = alt_count / total
                if (ratio < 0.20 || 0.30 < ratio) exit 1
            }
        }
    ' "${prefix}.somatic.vcf" <(gzip -dc "$fastq")
}

measure_somatic_alleles "$test_root/vaf" normal
measure_somatic_alleles "$test_root/vaf" tumor
measure_somatic_alleles "$test_root/vaf-tumor-only" tumor

# Variant modes must never fall back to legacy semantics for unsupported input.
if "$dwgsim_bin" --matched -N 10 -x /dev/null \
    "$test_root/reference.fa" "$test_root/rejected-bed" \
    >"$test_root/rejected-bed.log" 2>&1; then
    echo "test-matched-wgs: matched BED input was not rejected" >&2
    exit 1
fi
if "$dwgsim_bin" --matched -N 10 -F 0.25 \
    "$test_root/reference.fa" "$test_root/rejected-f" \
    >"$test_root/rejected-f.log" 2>&1; then
    echo "test-matched-wgs: legacy -F was accepted in matched mode" >&2
    exit 1
fi
if "$dwgsim_bin" -N 10 --somatic-rate 0.01 \
    "$test_root/reference.fa" "$test_root/rejected-somatic" \
    >"$test_root/rejected-somatic.log" 2>&1; then
    echo "test-matched-wgs: matched-only option was accepted without --matched" >&2
    exit 1
fi
if "$dwgsim_bin" --matched --tumor-only -N 10 \
    "$test_root/reference.fa" "$test_root/rejected-two-modes" \
    >"$test_root/rejected-two-modes.log" 2>&1; then
    echo "test-matched-wgs: matched and tumor-only modes were combined" >&2
    exit 1
fi
if "$dwgsim_bin" --tumor-only --normal-pairs 10 --tumor-pairs 10 \
    "$test_root/reference.fa" "$test_root/rejected-normal-pairs" \
    >"$test_root/rejected-normal-pairs.log" 2>&1; then
    echo "test-matched-wgs: tumor-only mode accepted normal pairs" >&2
    exit 1
fi
if "$dwgsim_bin" --tumor-only -N 10 -F 0.25 \
    "$test_root/reference.fa" "$test_root/rejected-tumor-f" \
    >"$test_root/rejected-tumor-f.log" 2>&1; then
    echo "test-matched-wgs: legacy -F was accepted in tumor-only mode" >&2
    exit 1
fi
cp "$test_root/reference.fa" "$test_root/unindexed.fa"
if "$dwgsim_bin" --matched -N 10 "$test_root/unindexed.fa" \
    "$test_root/rejected-index" >"$test_root/rejected-index.log" 2>&1; then
    echo "test-matched-wgs: unindexed matched reference was not rejected" >&2
    exit 1
fi
if "$dwgsim_bin" --tumor-only -N 10 "$test_root/unindexed.fa" \
    "$test_root/rejected-tumor-index" \
    >"$test_root/rejected-tumor-index.log" 2>&1; then
    echo "test-matched-wgs: unindexed tumor-only reference was not rejected" >&2
    exit 1
fi
if find "$test_root" -name '*.partial.*' -print -quit | grep -q .; then
    echo "test-matched-wgs: staging files remain after tests" >&2
    exit 1
fi

echo "test-matched-wgs: all checks passed"
