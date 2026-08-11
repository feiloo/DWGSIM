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

## Performance benchmark

Run the self-contained reads-only benchmark with:

```sh
make benchmark
```

By default it generates 250,000 pairs of 100 bp reads from the bundled small reference and writes its report to `build/benchmark/benchmark.txt`. The timed section includes DWGSIM startup, simulation, gzip compression, and both FASTQ writes; post-run FASTQ integrity and record-count checks are excluded.

The report includes read pairs/s, individual reads/s, bases/s, compressed-output throughput, elapsed and CPU time, CPU utilization, peak resident memory, and compressed output size. The workload and reference are configurable, for example:

```sh
make benchmark BENCHMARK_READ_PAIRS=1000000
make benchmark BENCHMARK_REFERENCE=reference/GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna
```

Use `make help` to list all benchmark variables. Results are most comparable when run on an otherwise idle machine with the same compiler flags, reference, read lengths, and storage.

## Full human-reference smoke test

The full-reference workflow is pinned to the official NCBI RefSeq GRCh38.p14 assembly (`GCF_000001405.40`) and verifies the download against NCBI's published MD5 checksum.

The download needs `curl`, `md5sum`, `gzip`, and roughly 5 GB of free disk space while unpacking. No additional genome-analysis toolkit is required.

```sh
make download
make test-human-reference
```

The smoke test runs DWGSIM against the complete assembly in reads-only mode with mutations disabled. It verifies all three compressed FASTQ outputs and checks that each contains the expected number of records. Use `make help` to see variables for the output location and read-pair count.
