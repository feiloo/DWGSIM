# Human reference smoke tests

## Commands

Download and checksum-verify the NCBI GRCh38.p14 reference and the pinned BED sources:

```console
make download
```

Build reference-ordered, RefSeq-name inclusion BEDs:

```console
make prepare-human-regions
```

Run all three profiles, or one profile:

```console
make test-human-reference
make test-human-wgs
make test-human-wes
make test-human-wgs-filtered
```

The defaults generate 100 pairs per profile. Override this with `HUMAN_SMOKE_READ_PAIRS`; this remains a smoke test, not a full-depth WGS or WES dataset.

## Profiles

| Target | BED | Purpose |
| --- | --- | --- |
| `test-human-wgs` | none | Exercise untargeted reads over the complete 705-sequence GRCh38.p14 RefSeq assembly. |
| `test-human-wes` | padded GIAB RefSeq CDS | Load and validate the complete published coding-target set across chromosomes 1–22, X, and Y. |
| `test-human-wgs-filtered` | complement of ENCODE `ENCFF356LFX` | Exercise a complete-assembly inclusion BED with known problematic hg38 regions removed. |

All profiles use reads-only mode, mutations disabled (`-r 0 -M 1`), and random reads disabled (`-y 0`). They verify each BGZF stream with `gzip --test`, check exact FASTQ line counts, and reject any random-read headers.

## Pinned sources and citations

The untouched source files are committed under `testdata/regions/upstream/`. Exact URLs, MD5 and SHA-256 hashes, dataset identifiers, limitations, and full citations are listed in [the region-data README](../testdata/regions/README.md).

In summary:

- WES-like targets use GIAB genome stratifications v3.6 `GRCh38_refseq_cds.bed.gz` from the NCBI GIAB release tree. Cite Olson et al., *Cell Genomics* (2022), [doi:10.1016/j.xgen.2022.100129](https://doi.org/10.1016/j.xgen.2022.100129), and the NIST dataset [doi:10.18434/mds2-2499](https://doi.org/10.18434/mds2-2499).
- Filtered WGS uses the released ENCODE GRCh38 blacklist file [ENCFF356LFX](https://www.encodeproject.org/files/ENCFF356LFX/), from annotation [ENCSR636HFF](https://www.encodeproject.org/annotations/ENCSR636HFF/). Cite Amemiya, Kundaje & Boyle, *Scientific Reports* (2019), [doi:10.1038/s41598-019-45839-z](https://doi.org/10.1038/s41598-019-45839-z).

## Reproducible transformations

`scripts/prepare_human_regions.sh` performs the following without GATK or bedtools:

1. Index the pinned NCBI FASTA with the bundled `samtools faidx` to obtain exact contig order and lengths.
2. Translate source `chr*` names to NCBI RefSeq accessions using the pinned GRCh38.p14 assembly report.
3. Pad every GIAB CDS interval by `HUMAN_WES_PADDING` (250 bp by default), clip it to the chromosome, and merge overlapping or touching intervals.
4. Subtract the 910 intervals in ENCODE file `ENCFF356LFX` from the matching primary chromosomes while retaining every other FASTA contig in full.
5. Validate BED syntax, bounds, FASTA contig order, interval order, and nonempty output.

At the default padding, the WES BED has 153,551 merged intervals covering 121,496,105 bases. The filtered-WGS BED has 1,615 intervals covering 3,226,860,351 bases. The latter equals the 3,298,430,636-base FASTA minus 71,570,285 blacklisted bases.

The generated files are:

```text
reference/GRCh38.p14/regions/grch38-p14-refseq-cds-padded.bed
reference/GRCh38.p14/regions/grch38-p14-without-encode-blacklist.bed
```

DWGSIM requires the complete paired-end fragment to be contained by one target interval. Padding is therefore required for practical exon simulation. These profiles test target handling; they do not model capture efficiency, off-target capture, GC bias, or a specific commercial WES panel.
