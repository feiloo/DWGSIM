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

## Full human-reference validation

The full-reference workflow is pinned to the official NCBI RefSeq GRCh38.p14 assembly (`GCF_000001405.40`) and verifies the download against NCBI's published MD5 checksum.

The download needs `curl`, `md5sum`, `gzip`, and roughly 5 GB of free disk space while unpacking. The validation needs GATK 4 available as `gatk` (or supplied with `GATK=/path/to/gatk`).

```sh
make download
make test-human-reference
```

The test builds the FASTA index, creates a GATK sequence dictionary, simulates reads across the complete assembly, checks the generated FASTQ archives, and runs GATK `ValidateVariants` against the generated VCF. Use `make help` to see variables for the output location, read-pair count, mutation rate, and GATK executable.
