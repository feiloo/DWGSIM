# Simulating Reads

<!---toc start-->
  * [Overview](#overview)
  * [Targeted reads and BED files](#targeted-reads-and-bed-files)
  * [Matched normal/tumor and tumor-only WGS](#matched-normaltumor-and-tumor-only-wgs)
  * [Error rates explained](#error-rates-explained)
  * [Read names explained](#read-names-explained)
  * [Mate pair or paired end modes](#mate-pair-or-paired-end-modes)
  * [Output mutations file](#output-mutations-file)
  * [Output FASTQ files](#output-fastq-files)

<!---toc end-->

## Overview

The `dwgsim` tool simulates reads from a reference genome FASTA.
It can be used to evaluate both mapping and variant calling.
See the `dwgsim_eval` tool for [Evaluating Mappings](04_Evaluating_Mappings.md).

Use `dwgsim -h` to see the full set of command line options (usage).

An example command is as follows:

```console
dwgsim -N 10000 -1 150 -2 150 -y 0 phix.fasta output
```

This will simulate 10,000 read pairs (`-N 10000`), which are paired-end 2x150 bp (`-1 150` for R1 and `-2 150` for R2), with no random reads (`-y 0`),
from the given genome FASTA (`phix.fasta`), producing output with prefix `output`.

The default is paired-end 2x150 bp. R1 and R2 are emitted as synchronized
`output.bwa.read1.fastq.gz` and `output.bwa.read2.fastq.gz` files; use `-o 1`
when those two per-end files are the only FASTQ outputs needed.

The following output will be created:

| Name | Description |
| --- | --- |
| output.bfast.fastq.gz | Interleaved FASTQ containing both read one and read two |
| output.bwa.read1.fastq.gz | FASTQ containing only read one |
| output.bwa.read2.fastq.gz | FASTQ containing only read two |
| output.mutations.vcf | VCF containing simulated mutations |
| output.mutations.txt | TXT in a custom format containing simulated mutations (see [below](#output-mutations-file)) |


Notes:

- matched and tumor-only modes limit generated insertion/deletion events to
  1,048,576 bases
- The `-H` mode will simulate a haploid genome, whereas the default is to simulate a diploid genome. 

## Matched normal/tumor and tumor-only WGS

Use `--matched` for a true matched-sample model instead of the legacy
`-F` haplotype-sampling bias:

```console
samtools faidx reference.fa
dwgsim --matched -N 1000000 -z 13 -r 0.001 \
  --somatic-rate 0.00001 --tumor-vaf 0.25 reference.fa matched
```

The command emits four paired-end BGZF FASTQs, phased germline and somatic
truth VCFs, and `matched.matched.complete`. `-N` applies to each
sample; `--normal-pairs` and `--tumor-pairs` allow unequal
libraries. Matched read names add `normal_` or `tumor_` before
the contig field. The remaining fields and `/1`/`/2` pairing rules are
unchanged.

This optimized profile currently requires indexed, untargeted, paired
Illumina WGS with outer distance and BWA per-end output. It does not accept
`-C`, BED regions, supplied mutation files, random reads, haploid mode,
amplicon mode, or BFAST output. Unsupported combinations fail explicitly.
See [Matched normal/tumor simulation](08_Matched_Normal_Tumor.md) for the
biological model and exact output contract.

Use `--tumor-only` to emit only the tumor R1/R2 FASTQs while retaining the
same generated germline and somatic truth:

```console
dwgsim --tumor-only -N 1000000 -z 13 -r 0.001 \
  --somatic-rate 0.00001 --tumor-vaf 0.25 reference.fa tumor
```

The outputs are `tumor.tumor.bwa.read1.fastq.gz`,
`tumor.tumor.bwa.read2.fastq.gz`, the two truth VCFs, and
`tumor.tumor-only.complete`. No normal FASTQs are created. The VCFs keep
their `NORMAL` and `TUMOR` columns to describe the germline baseline and
somatic difference even though only tumor reads were requested.

## Targeted reads and BED files

Use `-x FILE` to restrict reference-derived reads to target intervals:

```console
dwgsim -x targets.bed -N 1000000 -1 150 -2 150 -d 350 -s 50 -y 0 reference.fa targeted
```

The regions file must follow this contract:

- BED3 or wider text; the first three fields are contig, start, and end. Extra fields are ignored.
- Coordinates are zero-based and half-open: `[start,end)`.
- Contig names must exactly match the FASTA names.
- Records must follow FASTA contig order, then nondecreasing start coordinate within each contig.
- Every interval must satisfy `0 <= start < end <= contig_length`.
- Blank lines, lines beginning with `#`, and UCSC `track` and `browser` lines are accepted.
- Input must be uncompressed. Decompress `.bed.gz` files before passing them to `-x`.

Malformed coordinates, empty files, unknown contigs, invalid bounds, and ordering errors are rejected with the BED line number. Overlapping or directly adjacent intervals are merged in memory.

DWGSIM requires the entire outer paired-end fragment—not merely one read or any overlap—to fit within one merged interval. Raw exon/CDS intervals should therefore be padded and merged by at least the intended fragment flank. This is targeted placement, not a capture model: it does not reproduce bait efficiency, GC bias, uneven depth, duplicate formation, or off-target capture. Set `-y 0` when all reads must come from the target intervals.

`-N` distributes an exact read-pair count over the target bases; `-C` computes mean coverage from their total length. Regions BED mode cannot be combined with amplicon mode (`-a`).

The reproducible GRCh38 RefSeq CDS and blacklist-complement examples used by the smoke tests are described in [Human reference smoke tests](06_Human_Reference_Smoke_Tests.md).

## Error rates explained 

The `-e` and `-E` options accept a uniform error rate (i.e. `-e 0.01` for 1%), or a uniformly increasing/decreasing error rate (i.e. `-e 0.01-0.1` for an error rate of 1% at the start of the read increasing to 10% at the end of the read).

## Read names explained

Read names are of the form:

```
 @<#1>_<#2>_<#3>_<#4>_<#5>_<#6>_<#7>_<#8>:<#9>:<#10>_<#11>:<#12>:<#13>_<#14>
```

| Field | Description |
| --- | --- |
| 1 | contig name (chromsome name) |
| 2 | start read 1 (one-based) |
| 3 | start read 2 (one-based) |
| 4 | strand read 1 (0 - forward, 1 - reverse) |
| 5 | strand read 2 (0 - forward, 1 - reverse) |
| 6 | random read 1 (0 - from the mutated reference, 1 - random) |
| 7 | random read 2 (0 - from the mutated reference, 1 - random) |
| 8 | number of sequencing errors read 1 (color errors for colorspace) |
| 9 | number of SNPs read 1 |
| 10 | number of indels read 1 |
| 11 | number of sequencing errors read 2 (color errors for colorspace) |
| 12 | number of SNPs read 2 |
| 13 | number of indels read 2 |
| 14 | read number (unique within a given contig/chromsome) |

Read 1 and read 2 correspond to the first and second reads from a paired-end/mate-pair read respectively.

## Mate pair or paired end modes

This utility can generate mate pair or paired end reads using the `-S` option.
By default, Illumina (nucleotide) data are paired end, and SOLiD (color space) data are mate pair.
For clarity, lets call the first end sequence E1 and the second end E2.

Paired end reads have the following orientation:

```
 5' E1 -----> ....             3'
 3'           .... <------- E2 5'
```

Above, the start co-ordinate of E1 is less than E2, with E1 and E2 reported on opposite strands.

Mate pair reads have following orientation

```
 5' E2 -----> .... E1 -------> 3'
 3'           ....             5'
```

Above, the start co-ordinate of E1 is greater than E2, with E1 and E2 reported on the same strand.

So for SOLiD mate pair reads, the R3 tag (E2) is listed before the F3 tag (E1).
For SOLiD paired end reads, the F3 tag (E1) is listed before the F5 tag (E2).

## Output mutations file

The locations of introduced mutations are given in a `<prefix>.mutations.txt` text file.
There are file columns:

1. the chromosome/contig name
2. the one-based position
3. the original reference base
4. the new reference base(s)
5. the variant strand(s)

SNPs are represented on one line, and in the case of heterozygous mutations, the new reference base is an [IUPAC](http://www.bioinformatics.org/sms/iupac.html) code.

```
 contig4   4   T   K   1
```

The above shows a heterozygous mutation at position 4 of contig4 on the first strand, mutating the T base to a heterozygous K (G or T) SNP.

Insertions are represented on one line, where the reference base is missing (indicated by a '-' in the third column).

```
 contig5   13   -   TAC   3
```

The above shows a homozygous insertion of TAC prior to position 13.

Each base of a deletion is represented on one line, where the new reference base is missing and represented by a '-'.

```
 contig6   22   A   -   2
```

The above shows a heterozygous deletion of T at position 22 on the second strand.
Multi-base deletions are show on consecutive lines.

```
 contig6   22   A   -   2
 contig6   23   C   -   2
```

The above shows a two base homozygous deletion of positions 22 and 23 on the second strand.

## Output FASTQ files

Three FASTQ files are produced, for use with BFAST (interleaved FASTQ) and BWA (one FASTQ per read end).

FASTQ output uses gzip-compatible BGZF compression at level 4 by default. Ordinary gzip readers continue to work; BGZF-aware tools can also process the independent compressed blocks. Use `-l INT` to select level 1-9; level 1 favors speed while level 4 is the default size/runtime balance.

Use `-t INT` to set the generation/compression worker count. The default is the number of online logical CPUs detected at startup. Indexed mutation-free BWA WGS and indexed generated-variant `--matched`/`--tumor-only` WGS use the deterministic parallel generator. Fixed tasks, task-local BGZF chunks, and ordered appenders make all complete FASTQs byte-identical across worker counts and repeated runs for identical non-thread inputs. Mutation-free mode publishes `<prefix>.dwgsim.complete`; matched mode publishes `<prefix>.matched.complete` after four FASTQs and two truth VCFs are present; tumor-only mode publishes `<prefix>.tumor-only.complete` after two tumor FASTQs and both truth VCFs are present.

Other modes retain the legacy simulator. In that path, `-t` is distributed across BGZF helpers while read simulation remains serial. BED-restricted WGS/WES, single-sample mutation generation, mutation input, random reads, BFAST/combined output, inner-distance pairs, and non-Illumina data currently use this fallback.

The implementation and validation strategy are documented in [BGZF FASTQ output](05_BGZF_FASTQ_Output.md) and [Deterministic parallel generation](07_Deterministic_Parallel_Generation.md).

The FASTQ for BFAST is formatted so that the multi-end reads (paired end or mate pair) occur consecutively in the FASTQ (interleaved), with the read that is 5' of the other listed first.
For paired end reads, this means that E1 is always listed before E2.
For mate pair reads, this means that E2 is always listed before E1.

The FASTQs for BWA are split into two files, the first file for one end, the second file for the other, with the read that is 5' of the other in the first file.
For paired end reads, this means that E1 is in the first file and E2 is in the second file.
For mate pair reads, this means that E2 is in the first file and E1 is in the second file.
