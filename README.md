[![Build Status](https://github.com/nh13/DWGSIM/actions/workflows/testing.yml/badge.svg)](https://github.com/nh13/DWGSIM/actions/workflows/testing.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/dwgsim/README.html)
[![License](http://img.shields.io/badge/license-GPLv2-blue.svg)](https://github.com/nh13/dwgsim/blob/main/LICENSE)
[![Language](http://img.shields.io/badge/language-C-brightgreen.svg)](https://en.wikipedia.org/wiki/C_(programming_language))


Welcome to DWGSIM.
----

Documentation can be found in the [docs folder](docs/01_Introduction.md)

## Building

```sh
make
make help
```

`make help` lists the build, test, download, and configurable full-reference targets.

## BGZF FASTQ compression

DWGSIM writes FASTQ files as gzip-compatible BGZF using the implementation bundled with this repository. Compression level 4 is the default size/runtime balance; select levels 1-9 with `-l` (`-l 1` is the fast profile). No external `bgzip` executable is required at runtime.

Use `-t` to set the generation/compression worker count:

```sh
./dwgsim -t 4 reference.fa output
```

By default, DWGSIM uses all online logical CPUs. There are two deterministic
optimized WGS profiles:

- mutation-free, random-read-free, reads-only BWA pairs
  (`-r 0 -y 0 -M 1 -o 1`); and
- generated germline plus somatic variants in matched normal/tumor
  (`--matched`) or tumor-only (`--tumor-only`) mode.

Both load indexed contigs in parallel, generate fixed 8,192-pair tasks with
schedule-independent RNG, compress task-local BGZF chunks, and commit streams
in canonical order. The mutation-free and tumor-only profiles commit R1/R2;
matched mode commits normal R1/R2 and tumor R1/R2 together. Complete FASTQs
are byte-identical across worker counts and repeated runs for the same build,
seed, reference/index, non-thread options, and compression level. Successful
runs publish a completion manifest only after every required file is present.

BED-restricted WGS/WES, supplied mutation files, random reads,
BFAST/combined output, inner-distance pairs, and non-Illumina modes still use
the compatibility simulator. Build a missing FASTA index with
`samtools faidx reference.fa`. See [the BGZF implementation
notes](docs/05_BGZF_FASTQ_Output.md), [deterministic parallel generation
design/status](docs/07_Deterministic_Parallel_Generation.md), and
[matched normal/tumor simulation](docs/08_Matched_Normal_Tumor.md).

## Matched and tumor-only WGS

`--matched` generates independent paired libraries that share one phased
germline truth set while only the tumor can carry somatic events:

```sh
samtools faidx reference.fa
./dwgsim --matched -N 1000000 -z 13 \
  -r 0.001 --somatic-rate 0.00001 --tumor-vaf 0.25 \
  reference.fa matched
```

`-N` is the pair count for each sample. Different library sizes can be
requested with `--normal-pairs` and `--tumor-pairs`. Output is:

- `matched.normal.bwa.read1.fastq.gz` and
  `matched.normal.bwa.read2.fastq.gz`;
- `matched.tumor.bwa.read1.fastq.gz` and
  `matched.tumor.bwa.read2.fastq.gz`;
- `matched.germline.vcf` and `matched.somatic.vcf`; and
- `matched.matched.complete`.

Both mates of a fragment use the same haplotype and tumor-clone state.
`--tumor-vaf` is the expected diploid somatic allele fraction (0 to
0.5); it does not model copy number, arbitrary tumor purity, or multiple
clones. The exact model, VCF representation, limitations, and validation are
documented in [Matched normal/tumor simulation](docs/08_Matched_Normal_Tumor.md).

To generate only the paired tumor library with the same germline/somatic
model, use:

```sh
./dwgsim --tumor-only -N 1000000 -z 13 \
  -r 0.001 --somatic-rate 0.00001 --tumor-vaf 0.25 \
  reference.fa tumor
```

This publishes `tumor.tumor.bwa.read1.fastq.gz`,
`tumor.tumor.bwa.read2.fastq.gz`, both truth VCFs, and
`tumor.tumor-only.complete`; it does not create normal FASTQs. Use
`--tumor-pairs` instead of `-N` when a sample-specific count is clearer.
The truth VCFs retain both `NORMAL` and `TUMOR` genotype columns so the
germline baseline and tumor-only scope remain explicit. With identical
simulation options, tumor-only FASTQs and truth VCFs are byte-identical to
the tumor side and truth files of a matched run.

## Performance benchmark

Run the self-contained reads-only benchmark with:

```sh
make benchmark
```

By default it generates 250,000 paired-end 2x150 bp reads from the bundled small reference and writes its report to `build/benchmark/benchmark.txt`. The timed section includes DWGSIM startup, simulation, gzip compression, and both FASTQ writes; post-run FASTQ integrity and record-count checks are excluded.

For a representative full-reference WGS measurement and a nominal 100x projection, use:

```sh
make benchmark-wgs
```

This uses the complete pinned GRCh38.p14 assembly, disables mutations and random off-reference reads, and emits only the paired BWA FASTQs to avoid duplicate interleaved output. The default sample is five million paired-end 2x150 bp pairs, large enough to distinguish sustained generation from fixed reference work. Paired ends are written as synchronized `*.bwa.read1.fastq.gz` and `*.bwa.read2.fastq.gz` files. After the sample, it times a one-pair full-reference run and subtracts that fixed scan/setup cost from throughput before projecting production runtime. The report includes both raw and startup-adjusted rates, exact reference length, estimated pairs, duration, and compressed output size for `BENCHMARK_ESTIMATE_COVERAGE` (100 by default).

The report includes read pairs/s, individual reads/s, bases/s, compressed-output throughput, elapsed and CPU time, CPU utilization, peak resident memory, compression level, and compressed output size. The workload and reference are configurable, for example:

```sh
make benchmark BENCHMARK_READ_PAIRS=1000000
make benchmark BENCHMARK_REFERENCE=reference/GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna
make benchmark BENCHMARK_THREADS=4
make benchmark BENCHMARK_COMPRESSION_LEVEL=1
make benchmark-wgs WGS_BENCHMARK_READ_PAIRS=5000000
```

On the eight-core development host, the implemented path generated five million GRCh38.p14 2x150 pairs at 484,070.87 startup-adjusted pairs/s (968,141.74 reads/s, 145.22 Mbases/s) with BGZF level 4. Peak RSS was 3,250 MiB, CPU utilization was 686%, and the resulting 100x projection was 37m52s and 293.05 GiB. The report is retained at `build/benchmark-parallel-wgs/benchmark.txt`; results depend strongly on CPU and storage.

Use `make help` to list all benchmark variables. Results are most comparable when run on an otherwise idle machine with the same compiler flags, reference, read lengths, and storage.

## Full human-reference smoke test

The full-reference workflow is pinned to the official NCBI RefSeq GRCh38.p14 assembly (`GCF_000001405.40`). It includes three reads-only smoke profiles:

- WGS over the complete assembly;
- WES-like reads over the complete published GIAB RefSeq CDS set, padded for fragment placement;
- WGS over the complement of the released ENCODE GRCh38 blacklist (`ENCFF356LFX`).

The reference and vendored BED sources are checksum-verified. Preparing the BEDs only uses shell tools and the bundled `samtools faidx`; GATK, bedtools, and bcftools are not required.

```sh
make download
make prepare-human-regions
make test-human-reference
```

Each profile can also be run independently:

```sh
make test-human-wgs
make test-human-wes
make test-human-wgs-filtered
```

These are small smoke simulations (100 pairs by default), not 30x production datasets. They scan the full assembly, disable mutations and random off-target reads, verify all BGZF FASTQs, and check exact record counts. The BED downloads, transformations, checksums, limitations, and publication citations are documented in [Human reference smoke tests](docs/06_Human_Reference_Smoke_Tests.md).
