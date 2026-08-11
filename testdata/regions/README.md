# GRCh38 region sources

The files in `upstream/` are unmodified, pinned source data used by the full-reference smoke tests. `make download-human-regions` verifies them and restores a missing file from the exact URL below.

## RefSeq coding regions for the WES-like profile

- Local file: `upstream/GRCh38_refseq_cds.bed.gz`
- Dataset: Genome in a Bottle (GIAB) genome stratifications v3.6, `Functional/GRCh38_refseq_cds.bed.gz`
- Exact URL: <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/Functional/GRCh38_refseq_cds.bed.gz>
- MD5: `c27674ff0559893df02833e4e7b0e7ca`
- SHA-256: `825dad3d0c0d01d8a3555e4b3393eebb13456d1d68ddee5784b31ecf36feb577`
- Dataset DOI: [NIST Genome Stratifications, 10.18434/mds2-2499](https://doi.org/10.18434/mds2-2499)
- Recommended citation: Olson N. et al. “precisionFDA Truth Challenge V2: Calling variants from short- and long-reads in difficult-to-map regions.” *Cell Genomics* 2(5), 100129 (2022). [doi:10.1016/j.xgen.2022.100129](https://doi.org/10.1016/j.xgen.2022.100129).

This is a complete RefSeq CDS stratification, not a commercial exome-capture bait design. GIAB derived it from the GRCh38.p13 RefSeq annotation; primary chromosome coordinates are unchanged in GRCh38.p14.

## Problematic regions for the filtered-WGS profile

- Local file: `upstream/ENCFF356LFX.bed.gz`
- Dataset: released ENCODE GRCh38 blacklist file [ENCFF356LFX](https://www.encodeproject.org/files/ENCFF356LFX/), annotation [ENCSR636HFF](https://www.encodeproject.org/annotations/ENCSR636HFF/)
- Exact URL: <https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz>
- MD5: `393688b4f06c9ce26165d47433dd8c37`
- SHA-256: `a9d086ce90ca67f933b29adfbff56ea768bbce2ec4b0d51e2d405e5a1d61bb56`
- Citation: Amemiya H.M., Kundaje A. & Boyle A.P. “The ENCODE Blacklist: Identification of Problematic Regions of the Genome.” *Scientific Reports* 9, 9354 (2019). [doi:10.1038/s41598-019-45839-z](https://doi.org/10.1038/s41598-019-45839-z).

The source has 910 intervals covering 71,570,285 bases. It was designed for anomalous signal in functional-genomics assays; this repository uses it as a compact, citable example of a WGS exclusion workflow, not as a universal sequencing or variant-calling mask.

## Assembly name mapping

- Local file: `upstream/GCF_000001405.40_GRCh38.p14_assembly_report.txt`
- Exact NCBI URL: <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt>
- MD5: `21f3ac4aa8245a99eb874082051b9dde`
- SHA-256: `64318ddff470b69b261a667d813210044f60d4ce654253a547db80ff73638d38`

The preparation script uses this report to translate UCSC `chr*` names to the RefSeq accessions in the NCBI FASTA. See [the human smoke-test documentation](../../docs/06_Human_Reference_Smoke_Tests.md) for the exact transformations.
