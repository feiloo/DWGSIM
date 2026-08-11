# Matched normal/tumor and tumor-only simulation

## Quick start

Matched mode generates two independent paired-end libraries from one indexed
reference:

```sh
samtools faidx reference.fa
./dwgsim --matched -N 1000000 -z 13 \
  -1 150 -2 150 \
  -r 0.001 \
  --somatic-rate 0.00001 \
  --tumor-vaf 0.25 \
  reference.fa matched
```

`-N` is the number of pairs in each library. To use different depths:

```sh
./dwgsim --matched \
  --normal-pairs 400000000 \
  --tumor-pairs 800000000 \
  -z 13 -r 0.001 --somatic-rate 0.00001 --tumor-vaf 0.25 \
  reference.fa matched
```

Tumor-only mode uses the same biological model but emits only the paired tumor
library:

```sh
./dwgsim --tumor-only -N 800000000 \
  -z 13 -r 0.001 --somatic-rate 0.00001 --tumor-vaf 0.25 \
  reference.fa tumor
```

In this mode `-N` is the tumor pair count; `--tumor-pairs` is an
equivalent sample-specific spelling. `--normal-pairs` is rejected, and
`--matched` and `--tumor-only` are mutually exclusive.

Both modes imply no random reads and BWA per-end output. They use all online
logical CPUs by default; `-t INT` selects another worker count. BGZF
level 4 is the default size/runtime balance and `-l 1` is the
faster/larger profile.

## Output contract

One successful `--matched` run publishes:

| File | Contents |
| --- | --- |
| `<prefix>.normal.bwa.read1.fastq.gz` | normal R1 |
| `<prefix>.normal.bwa.read2.fastq.gz` | normal R2 |
| `<prefix>.tumor.bwa.read1.fastq.gz` | tumor R1 |
| `<prefix>.tumor.bwa.read2.fastq.gz` | tumor R2 |
| `<prefix>.germline.vcf` | phased variants shared by NORMAL and TUMOR |
| `<prefix>.somatic.vcf` | tumor-only phased somatic variants |
| `<prefix>.matched.complete` | completion and byte-count manifest |

A successful `--tumor-only` run publishes:

| File | Contents |
| --- | --- |
| `<prefix>.tumor.bwa.read1.fastq.gz` | tumor R1 |
| `<prefix>.tumor.bwa.read2.fastq.gz` | tumor R2 |
| `<prefix>.germline.vcf` | phased germline baseline |
| `<prefix>.somatic.vcf` | tumor-only phased somatic variants |
| `<prefix>.tumor-only.complete` | completion and byte-count manifest |

It does not create empty or placeholder normal FASTQs. The emitted FASTQs are
staged, synchronized, and closed before publication. Both VCFs are also
staged. The applicable completion manifest is renamed last; consumers should
require `<prefix>.matched.complete` or
`<prefix>.tumor-only.complete` for the requested mode.

Tumor-only VCFs deliberately retain both `NORMAL` and `TUMOR` genotype
columns. `NORMAL` describes the germline baseline and makes the somatic
difference explicit even when no normal library was requested. For the same
seed and simulation options, these truth VCFs and the two tumor FASTQs have
the same bytes as the truth VCFs and tumor side of a matched run.

Normal and tumor libraries have independent deterministic fragment locations,
errors, and qualities. Within either sample, R1 and R2 are created as one
fragment: their normalized identifiers match, and both mates use the same
haplotype and tumor-clone membership.

Variant-mode names begin with `normal_` or `tumor_`. The remaining
coordinate, strand, sequencing-error, biological-substitution, biological-
indel, pair-ordinal, and `/1`/`/2` fields retain the normal DWGSIM
layout.

## Biological model

### Germline

`-r FLOAT` is the per-reference-base germline event rate. Events are
SNVs unless selected as indels by `-R`. For each germline event:

- one third are homozygous alternate (`1|1`);
- two thirds are heterozygous and are assigned deterministically to haplotype
  one (`1|0`) or two (`0|1`); and
- NORMAL and TUMOR receive the same phased genotype.

Each fragment selects one haplotype uniformly, and both mates are extracted
from it.

### Somatic

`--somatic-rate FLOAT` is the per-reference-base tumor-only event
rate. Every accepted somatic event is heterozygous and assigned to one phased
haplotype. The somatic truth VCF reports `NORMAL=0|0`, the phased
TUMOR genotype, the `SOMATIC` flag, and `EXPECTED_VAF`.

`--tumor-vaf V` accepts `0 <= V <= 0.5`. For each tumor fragment,
the simulator selects the somatic-carrying clone with probability `2V` and
then selects either haplotype uniformly. A heterozygous somatic allele is
therefore observed with expected probability:

```text
P(clone) * P(event haplotype) = (2V) * (1/2) = V
```

All somatic events share that per-fragment clone state. This is a simple
diploid mixture model. It does not model locus-specific cancer cell fractions,
multiple subclones, normal contamination as a separate parameter, loss of
heterozygosity, or copy-number changes. Those require an extended genome/
clone model rather than reinterpretation of `--tumor-vaf`.

Sequencing errors are applied after germline and somatic alleles. Consequently
the observed FASTQ allele fraction can differ from the biological expectation
because of finite sampling, errors, indels, ambiguous reference bases, and
placement coverage.

### SNVs and indels

`-R` sets the indel fraction, `-X` the geometric extension
probability, and `-I` the minimum indel length. Insertions and deletions
are selected equally within the indel fraction. Generated events are
limited to 1,048,576 bases.

Events overlapping another accepted germline or somatic event are rejected.
Events requiring ambiguous reference bases are also rejected. Deletions use
the preceding reference base as the VCF anchor; insertions use their reference
position as the anchor. VCF coordinates are one-based and alleles use the
standard anchored representation.

The rates are event rates, not alternate-allele fractions. Rejection at
ambiguous bases, contig boundaries, and existing events means the final count
can be slightly below `rate * eligible_reference_length`.

## Parallel implementation

The optimized path avoids the legacy pair of full mutated-genome arrays:

1. FASTA-index records define immutable contigs in canonical order.
2. Contigs are loaded with independent indexed reads.
3. Contig workers jump geometrically between candidate event positions, so
   preparation is proportional to event count rather than a serial genome
   scan for ordinary rates.
4. Accepted germline and somatic events are stored once in one sparse sorted
   vector per contig.
5. Fixed 8,192-pair tasks independently generate the requested logical
   samples: normal plus tumor, or tumor only.
   Read extraction performs one binary search and then walks nearby events.
6. Each worker compresses four task-local BGZF streams for matched mode or two
   for tumor-only mode.
7. One ordered appender per stream consumes the bounded result ring; a result
   is freed only after every requested R1/R2 stream has consumed it.

The reference and sparse variants are shared read-only. Memory therefore
scales with the encoded reference, event count, active task buffers, and
bounded reorder window—not with the number of workers times genome size.

### Full-reference implementation check

The optimized implementation was measured on the complete 3,298,430,636-base
GRCh38.p14 assembly with eight logical CPUs, BGZF level 4, one million 2x150
pairs per sample, germline rate 0.001, somatic rate 0.00001, and tumor VAF
0.25. The timed run included reference loading, sparse variant preparation,
both truth VCFs, read/error/quality generation, compression, and file writes:

| Metric | Result |
| --- | ---: |
| Accepted germline events | 3,133,920 |
| Accepted somatic events | 31,519 |
| Combined sample-pairs/s | 329,992.02 |
| Individual reads/s | 659,984.04 |
| Simulated bases/s | 99.00 million |
| Internal elapsed time | 6.061 s |
| Process elapsed time | 6.27 s |
| CPU utilization | 575% |
| Peak RSS | 3,469,080 KiB (3.31 GiB) |
| Four FASTQs | 575,862,807 bytes |
| Germline plus somatic VCFs | 269,447,090 bytes |

All four FASTQs passed gzip integrity validation. This is evidence for the
development host and its warm-cache/storage state, not a universal throughput
guarantee.

## Determinism

For the same executable/build platform, reference and index bytes, seed,
non-thread simulation options, and BGZF level:

- every complete FASTQ file (four matched or two tumor-only) has identical
  bytes at every worker count;
- both truth VCFs have identical bytes at every worker count; and
- repeated runs reproduce those bytes.

The tumor RNG identity is independent of whether the normal output streams
are active. Consequently, an equivalent matched and tumor-only run has
byte-identical tumor FASTQs and truth VCFs.

Thread count changes scheduling only. Variant generation, genotypes, fragment
locations, haplotypes, clone selection, errors, qualities, task boundaries,
record order, BGZF block boundaries, and truth order use stable identities.
The completion manifest records the requested thread count, so the manifest
itself is intentionally metadata rather than a cross-thread hash target.

Hash identity is scoped to the same build platform and zlib implementation;
cross-architecture/compressor identity is not currently promised.

## Current constraints

Matched and tumor-only modes currently require:

- a readable FASTA plus `<reference>.fai`;
- paired Illumina nucleotide reads with outer distance;
- untargeted WGS;
- BWA per-end FASTQ output; and
- explicit pair counts through `-N`, both sample-specific options for
  matched mode, or `--tumor-pairs` for tumor-only mode.

They reject coverage-driven `-C`, `-x` BED regions, supplied
`-m`/`-b`/`-v` mutation files, random reads, haploid mode,
amplicon mode, inner-distance pairs, BFAST output, SOLiD, and Ion Torrent.
Rejection is deliberate: neither mode silently falls back to the legacy
`-F` approximation. BED/WES, supplied truth, CNAs, purity, and multi-clone
tumors are future extensions.

## Validation

Run:

```sh
make test-matched-wgs
```

The test creates 160 fixed tasks and checks:

- byte identity for matched and tumor-only output at 1, 8, repeated 8, 128
  requested workers, and the automatic all-online-CPU default;
- exact counts, 150-base sequences/qualities, synchronized mates, gzip and
  BGZF decoding, and canonical EOF blocks for every emitted FASTQ;
- deterministic phased VCFs containing SNVs, insertions, and deletions;
- shared germline and tumor-only somatic genotypes;
- independent normal/tumor pair counts and empty-sample task chunks;
- exact equality between tumor-only files and the tumor/truth files from an
  equivalent matched run, with no normal files published;
- zero somatic alternate observations in an error-free normal SNV run and a
  measured tumor allele fraction within a statistical tolerance of the
  configured VAF; and
- explicit rejection of incompatible mode combinations, normal counts in
  tumor-only mode, BED input, and unindexed variant-mode input.

The implementation has also passed AddressSanitizer/UndefinedBehaviorSanitizer
for the complete test and ThreadSanitizer at 8 workers and a 160-task,
128-requested-worker stress case.
