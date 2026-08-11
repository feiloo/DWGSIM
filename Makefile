.DEFAULT_GOAL := all

PACKAGE_VERSION=0.1.17-dev
CC=			gcc
CFLAGS=		-g -Wall -O3 #-m64 #-arch ppc
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -DPACKAGE_VERSION=\"$(PACKAGE_VERSION)\"
DWGSIM_AOBJS = src/dwgsim_opt.o src/mut.o src/contigs.o src/regions_bed.o \
			   src/mut_txt.o src/mut_bed.o src/mut_vcf.o src/mut_input.o src/fastq_writer.o \
			   src/parallel_wgs.o src/dwgsim.o
DWGSIM_EVAL_AOBJS = src/dwgsim_eval.o \
					samtools/knetfile.o \
					samtools/bgzf.o samtools/kstring.o samtools/bam_aux.o samtools/bam.o samtools/bam_import.o samtools/sam.o samtools/bam_index.o \
					samtools/bam_pileup.o samtools/bam_lpileup.o samtools/bam_md.o samtools/razf.o samtools/faidx.o samtools/bedidx.o \
					samtools/bam_sort.o samtools/sam_header.o samtools/bam_reheader.o samtools/kprobaln.o samtools/bam_cat.o
DWGSIM_MUT_TO_VCF_AOBJS = src/dwgsim_mut_to_vcf.o
DWGSIM_PILEUP_EVAL_AOBJS = src/dwgsim_pileup_eval.o

PROG=		dwgsim dwgsim_eval dwgsim_mut_to_vcf dwgsim_pileup_eval
INCLUDES=	-I.
SUBDIRS=	samtools . 
CLEAN_SUBDIRS=	samtools src
LIBPATH=

CURL ?= curl
MD5SUM ?= md5sum
AVAILABLE_CPU_THREADS ?= $(shell getconf _NPROCESSORS_ONLN 2>/dev/null || printf '1')
REFERENCE_ROOT ?= reference
HUMAN_REFERENCE_RELEASE ?= GCF_000001405.40_GRCh38.p14
HUMAN_REFERENCE_BASENAME ?= $(HUMAN_REFERENCE_RELEASE)_genomic
HUMAN_REFERENCE_DIR ?= $(REFERENCE_ROOT)/GRCh38.p14
HUMAN_REFERENCE_ARCHIVE ?= $(HUMAN_REFERENCE_DIR)/$(HUMAN_REFERENCE_BASENAME).fna.gz
HUMAN_REFERENCE_FASTA ?= $(HUMAN_REFERENCE_DIR)/$(HUMAN_REFERENCE_BASENAME).fna
HUMAN_REFERENCE_URL ?= https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/$(HUMAN_REFERENCE_BASENAME).fna.gz
HUMAN_REFERENCE_MD5 ?= c30471567037b2b2389d43c908c653e1
HUMAN_REGIONS_SOURCE_DIR ?= testdata/regions/upstream
HUMAN_WES_SOURCE ?= $(HUMAN_REGIONS_SOURCE_DIR)/GRCh38_refseq_cds.bed.gz
HUMAN_WES_SOURCE_URL ?= https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/Functional/GRCh38_refseq_cds.bed.gz
HUMAN_WES_SOURCE_MD5 ?= c27674ff0559893df02833e4e7b0e7ca
HUMAN_BLACKLIST_SOURCE ?= $(HUMAN_REGIONS_SOURCE_DIR)/ENCFF356LFX.bed.gz
HUMAN_BLACKLIST_SOURCE_URL ?= https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
HUMAN_BLACKLIST_SOURCE_MD5 ?= 393688b4f06c9ce26165d47433dd8c37
HUMAN_ASSEMBLY_REPORT ?= $(HUMAN_REGIONS_SOURCE_DIR)/GCF_000001405.40_GRCh38.p14_assembly_report.txt
HUMAN_ASSEMBLY_REPORT_URL ?= https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt
HUMAN_ASSEMBLY_REPORT_MD5 ?= 21f3ac4aa8245a99eb874082051b9dde
HUMAN_REGIONS_DIR ?= $(HUMAN_REFERENCE_DIR)/regions
HUMAN_WES_BED ?= $(HUMAN_REGIONS_DIR)/grch38-p14-refseq-cds-padded.bed
HUMAN_WGS_FILTERED_BED ?= $(HUMAN_REGIONS_DIR)/grch38-p14-without-encode-blacklist.bed
HUMAN_WES_PADDING ?= 250
HUMAN_SMOKE_DIR ?= build/human-reference-smoke
HUMAN_SMOKE_READ_PAIRS ?= 100
HUMAN_SMOKE_SEED ?= 13
HUMAN_SMOKE_THREADS ?= $(AVAILABLE_CPU_THREADS)
DWGSIM_BIN ?= ./dwgsim
TIME_BIN ?= /usr/bin/time
BENCHMARK_REFERENCE ?= samtools/examples/ex1.fa
BENCHMARK_DIR ?= build/benchmark
BENCHMARK_READ_PAIRS ?= 250000
BENCHMARK_READ_LENGTH_1 ?= 150
BENCHMARK_READ_LENGTH_2 ?= 150
BENCHMARK_SEED ?= 13
BENCHMARK_THREADS ?= $(AVAILABLE_CPU_THREADS)
BENCHMARK_COMPRESSION_LEVEL ?= 4
BENCHMARK_ESTIMATE_COVERAGE ?= 100
BENCHMARK_MEASURE_STARTUP ?= 0
WGS_BENCHMARK_READ_PAIRS ?= 5000000
WGS_BENCHMARK_DIR ?= build/benchmark-wgs

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all-recur lib-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)" $$target || exit 1; \
			cd $$wdir; \
		done;

src/dwgsim_mut_to_vcf.o: src/dwgsim_mut_to_vcf.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

src/dwgsim_pileup_eval.o: src/dwgsim_pileup_eval.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

src/fastq_writer.o: src/fastq_writer.c src/fastq_writer.h samtools/bgzf.h
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

src/parallel_wgs.o: src/parallel_wgs.c src/parallel_wgs.h src/dwgsim_opt.h samtools/bgzf.h
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

.PHONY:all lib clean cleanlocal test test-unit test-integration test-bgzf test-parallel-wgs test-matched-wgs test-bed clean-tests help
.PHONY:download download-human-reference download-human-regions prepare-human-regions samtools-program
.PHONY: test-human-reference test-human-wgs test-human-wes test-human-wgs-filtered
.PHONY:benchmark benchmark-wgs
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

help:
	@printf 'Usage: make <target> [VARIABLE=value]\n\n'
	@printf 'Build and test targets:\n'
	@printf '  %-26s %s\n' 'all' 'Build all DWGSIM executables (default).'
	@printf '  %-26s %s\n' 'test' 'Run unit, integration, BGZF, parallel/variant-WGS, and BED tests.'
	@printf '  %-26s %s\n' 'benchmark' 'Measure reads-only simulation throughput and resource use.'
	@printf '  %-26s %s\n' 'benchmark-wgs' 'Benchmark full GRCh38 WGS and estimate 100x generation.'
	@printf '  %-26s %s\n' 'test-bgzf' 'Test BGZF compatibility, modes, and thread determinism.'
	@printf '  %-26s %s\n' 'test-parallel-wgs' 'Test deterministic paired WGS across worker counts.'
	@printf '  %-26s %s\n' 'test-matched-wgs' 'Test deterministic matched and tumor-only WGS and truth.'
	@printf '  %-26s %s\n' 'test-bed' 'Test BED boundaries, headers, and validation errors.'
	@printf '  %-26s %s\n' 'download' 'Download/verify GRCh38.p14 and pinned human BED sources.'
	@printf '  %-26s %s\n' 'prepare-human-regions' 'Build RefSeq-name WES and blacklist-complement BEDs.'
	@printf '  %-26s %s\n' 'test-human-reference' 'Run all full-reference WGS/WES smoke profiles.'
	@printf '  %-26s %s\n' 'test-human-wgs' 'Smoke-test untargeted WGS on full GRCh38.p14.'
	@printf '  %-26s %s\n' 'test-human-wes' 'Smoke-test padded GIAB RefSeq CDS targets.'
	@printf '  %-26s %s\n' 'test-human-wgs-filtered' 'Smoke-test WGS excluding ENCODE ENCFF356LFX.'
	@printf '  %-26s %s\n' 'clean' 'Remove compiled programs, objects, and regular test artifacts.'
	@printf '  %-26s %s\n' 'help' 'Show this help.'
	@printf '\nBenchmark settings:\n'
	@printf '  %-30s %s\n' 'BENCHMARK_REFERENCE=...' 'Input FASTA (default: $(BENCHMARK_REFERENCE)).'
	@printf '  %-30s %s\n' 'BENCHMARK_READ_PAIRS=...' 'Read pairs to generate (default: $(BENCHMARK_READ_PAIRS)).'
	@printf '  %-30s %s\n' 'BENCHMARK_READ_LENGTH_1=...' 'First-read length (default: $(BENCHMARK_READ_LENGTH_1)).'
	@printf '  %-30s %s\n' 'BENCHMARK_READ_LENGTH_2=...' 'Second-read length (default: $(BENCHMARK_READ_LENGTH_2)).'
	@printf '  %-30s %s\n' 'BENCHMARK_SEED=...' 'Random seed (default: $(BENCHMARK_SEED)).'
	@printf '  %-30s %s\n' 'BENCHMARK_THREADS=...' 'DWGSIM thread budget (default: $(BENCHMARK_THREADS)).'
	@printf '  %-30s %s\n' 'BENCHMARK_COMPRESSION_LEVEL=...' 'BGZF level, 1-9 (default: $(BENCHMARK_COMPRESSION_LEVEL)).'
	@printf '  %-30s %s\n' 'BENCHMARK_ESTIMATE_COVERAGE=...' 'Coverage projection (default: $(BENCHMARK_ESTIMATE_COVERAGE)x).'
	@printf '  %-30s %s\n' 'BENCHMARK_MEASURE_STARTUP=0|1' 'Subtract a one-pair fixed-cost run (default: $(BENCHMARK_MEASURE_STARTUP)).'
	@printf '  %-30s %s\n' 'BENCHMARK_DIR=...' 'Output and report directory (default: $(BENCHMARK_DIR)).'
	@printf '  %-30s %s\n' 'WGS_BENCHMARK_READ_PAIRS=...' 'Full-GRCh38 sample pairs (default: $(WGS_BENCHMARK_READ_PAIRS)).'
	@printf '  %-30s %s\n' 'WGS_BENCHMARK_DIR=...' 'Full-GRCh38 report directory (default: $(WGS_BENCHMARK_DIR)).'
	@printf '\nHuman-reference settings:\n'
	@printf '  %-26s %s\n' 'HUMAN_REFERENCE_DIR=...' 'Reference destination (default: $(HUMAN_REFERENCE_DIR)).'
	@printf '  %-26s %s\n' 'HUMAN_SMOKE_READ_PAIRS=...' 'Read pairs to simulate (default: $(HUMAN_SMOKE_READ_PAIRS)).'
	@printf '  %-26s %s\n' 'HUMAN_SMOKE_THREADS=...' 'Thread budget for each smoke profile (default: $(HUMAN_SMOKE_THREADS)).'
	@printf '  %-26s %s\n' 'HUMAN_SMOKE_DIR=...' 'Smoke-test output directory (default: $(HUMAN_SMOKE_DIR)).'
	@printf '  %-26s %s\n' 'HUMAN_WES_PADDING=...' 'Padding added around RefSeq CDS targets (default: $(HUMAN_WES_PADDING)).'

download: download-human-reference download-human-regions

download-human-reference:
	CURL="$(CURL)" MD5SUM="$(MD5SUM)" /bin/bash scripts/download_human_reference.sh \
		"$(HUMAN_REFERENCE_URL)" "$(HUMAN_REFERENCE_MD5)" \
		"$(HUMAN_REFERENCE_ARCHIVE)" "$(HUMAN_REFERENCE_FASTA)"

download-human-regions:
	CURL="$(CURL)" MD5SUM="$(MD5SUM)" /bin/bash scripts/download_human_regions.sh \
		"$(HUMAN_WES_SOURCE_URL)" "$(HUMAN_WES_SOURCE_MD5)" "$(HUMAN_WES_SOURCE)" \
		"$(HUMAN_BLACKLIST_SOURCE_URL)" "$(HUMAN_BLACKLIST_SOURCE_MD5)" "$(HUMAN_BLACKLIST_SOURCE)" \
		"$(HUMAN_ASSEMBLY_REPORT_URL)" "$(HUMAN_ASSEMBLY_REPORT_MD5)" "$(HUMAN_ASSEMBLY_REPORT)"

samtools-program:
	$(MAKE) -C samtools samtools

prepare-human-regions: samtools-program download-human-reference download-human-regions
	HUMAN_WES_PADDING="$(HUMAN_WES_PADDING)" \
		/bin/bash scripts/prepare_human_regions.sh \
		"$(HUMAN_REFERENCE_FASTA)" "./samtools/samtools" \
		"$(HUMAN_ASSEMBLY_REPORT)" "$(HUMAN_WES_SOURCE)" "$(HUMAN_BLACKLIST_SOURCE)" \
		"$(HUMAN_WES_BED)" "$(HUMAN_WGS_FILTERED_BED)"

test-human-wgs: all
	$(MAKE) --no-print-directory download-human-reference
	DWGSIM_BIN="$(DWGSIM_BIN)" \
		HUMAN_SMOKE_READ_PAIRS="$(HUMAN_SMOKE_READ_PAIRS)" \
		HUMAN_SMOKE_SEED="$(HUMAN_SMOKE_SEED)" \
		HUMAN_SMOKE_THREADS="$(HUMAN_SMOKE_THREADS)" \
		/bin/bash scripts/smoke_test_human_reference.sh \
		"$(HUMAN_REFERENCE_FASTA)" "$(HUMAN_SMOKE_DIR)" wgs

test-human-wes: all
	$(MAKE) --no-print-directory prepare-human-regions
	DWGSIM_BIN="$(DWGSIM_BIN)" \
		HUMAN_SMOKE_READ_PAIRS="$(HUMAN_SMOKE_READ_PAIRS)" \
		HUMAN_SMOKE_SEED="$(HUMAN_SMOKE_SEED)" \
		HUMAN_SMOKE_THREADS="$(HUMAN_SMOKE_THREADS)" \
		/bin/bash scripts/smoke_test_human_reference.sh \
		"$(HUMAN_REFERENCE_FASTA)" "$(HUMAN_SMOKE_DIR)" wes "$(HUMAN_WES_BED)"

test-human-wgs-filtered: all
	$(MAKE) --no-print-directory prepare-human-regions
	DWGSIM_BIN="$(DWGSIM_BIN)" \
		HUMAN_SMOKE_READ_PAIRS="$(HUMAN_SMOKE_READ_PAIRS)" \
		HUMAN_SMOKE_SEED="$(HUMAN_SMOKE_SEED)" \
		HUMAN_SMOKE_THREADS="$(HUMAN_SMOKE_THREADS)" \
		/bin/bash scripts/smoke_test_human_reference.sh \
		"$(HUMAN_REFERENCE_FASTA)" "$(HUMAN_SMOKE_DIR)" wgs-filtered "$(HUMAN_WGS_FILTERED_BED)"

test-human-reference: all
	$(MAKE) --no-print-directory prepare-human-regions
	DWGSIM_BIN="$(DWGSIM_BIN)" HUMAN_SMOKE_READ_PAIRS="$(HUMAN_SMOKE_READ_PAIRS)" \
		HUMAN_SMOKE_SEED="$(HUMAN_SMOKE_SEED)" HUMAN_SMOKE_THREADS="$(HUMAN_SMOKE_THREADS)" \
		/bin/bash scripts/smoke_test_human_reference.sh \
		"$(HUMAN_REFERENCE_FASTA)" "$(HUMAN_SMOKE_DIR)" wgs
	DWGSIM_BIN="$(DWGSIM_BIN)" HUMAN_SMOKE_READ_PAIRS="$(HUMAN_SMOKE_READ_PAIRS)" \
		HUMAN_SMOKE_SEED="$(HUMAN_SMOKE_SEED)" HUMAN_SMOKE_THREADS="$(HUMAN_SMOKE_THREADS)" \
		/bin/bash scripts/smoke_test_human_reference.sh \
		"$(HUMAN_REFERENCE_FASTA)" "$(HUMAN_SMOKE_DIR)" wes "$(HUMAN_WES_BED)"
	DWGSIM_BIN="$(DWGSIM_BIN)" HUMAN_SMOKE_READ_PAIRS="$(HUMAN_SMOKE_READ_PAIRS)" \
		HUMAN_SMOKE_SEED="$(HUMAN_SMOKE_SEED)" HUMAN_SMOKE_THREADS="$(HUMAN_SMOKE_THREADS)" \
		/bin/bash scripts/smoke_test_human_reference.sh \
		"$(HUMAN_REFERENCE_FASTA)" "$(HUMAN_SMOKE_DIR)" wgs-filtered "$(HUMAN_WGS_FILTERED_BED)"

benchmark: dwgsim
	DWGSIM_BIN="$(DWGSIM_BIN)" \
		TIME_BIN="$(TIME_BIN)" \
		BENCHMARK_READ_PAIRS="$(BENCHMARK_READ_PAIRS)" \
		BENCHMARK_READ_LENGTH_1="$(BENCHMARK_READ_LENGTH_1)" \
		BENCHMARK_READ_LENGTH_2="$(BENCHMARK_READ_LENGTH_2)" \
		BENCHMARK_SEED="$(BENCHMARK_SEED)" \
		BENCHMARK_THREADS="$(BENCHMARK_THREADS)" \
		BENCHMARK_COMPRESSION_LEVEL="$(BENCHMARK_COMPRESSION_LEVEL)" \
		BENCHMARK_ESTIMATE_COVERAGE="$(BENCHMARK_ESTIMATE_COVERAGE)" \
		BENCHMARK_MEASURE_STARTUP="$(BENCHMARK_MEASURE_STARTUP)" \
		/bin/bash scripts/benchmark.sh \
		"$(BENCHMARK_REFERENCE)" "$(BENCHMARK_DIR)"

benchmark-wgs: dwgsim download-human-reference
	$(MAKE) --no-print-directory benchmark \
		BENCHMARK_REFERENCE="$(HUMAN_REFERENCE_FASTA)" \
		BENCHMARK_READ_PAIRS="$(WGS_BENCHMARK_READ_PAIRS)" \
		BENCHMARK_DIR="$(WGS_BENCHMARK_DIR)" \
		BENCHMARK_ESTIMATE_COVERAGE="$(BENCHMARK_ESTIMATE_COVERAGE)" \
		BENCHMARK_MEASURE_STARTUP=1

dwgsim:lib-recur $(DWGSIM_AOBJS)
	$(CC) $(CFLAGS) -o $@ $(DWGSIM_AOBJS) samtools/libbam.a -lm -lz -lpthread

dwgsim_eval:lib-recur $(DWGSIM_EVAL_AOBJS)
	$(CC) $(CFLAGS) -o $@ $(DWGSIM_EVAL_AOBJS) -Lsamtools -lm -lz -lpthread

dwgsim_mut_to_vcf: $(DWGSIM_MUT_TO_VCF_AOBJS)
	$(CC) $(CFLAGS) -o $@ $(DWGSIM_MUT_TO_VCF_AOBJS)

dwgsim_pileup_eval: $(DWGSIM_PILEUP_EVAL_AOBJS)
	$(CC) $(CFLAGS) -o $@ $(DWGSIM_PILEUP_EVAL_AOBJS)

cleanlocal:
		rm -vfr gmon.out *.o a.out *.exe *.dSYM razip bgzip $(PROG) *~ *.a *.so.* *.so *.dylib; \
		wdir=`pwd`; \
		list='$(CLEAN_SUBDIRS)'; for subdir in $$list; do \
			if [ -d $$subdir ]; then \
				cd $$subdir; \
				pwd; \
				rm -vfr gmon.out *.o a.out *.exe *.dSYM razip bgzip $(PROG) *~ *.a *.so.* *.so *.dylib; \
				cd $$wdir; \
			fi; \
		done;

clean:cleanlocal-recur clean-tests

dist:clean
	if [ -f dwgsim-${PACKAGE_VERSION}.tar.gz ]; then \
        rm -rv dwgsim-${PACKAGE_VERSION}.tar.gz; \
	fi; \
	if [ -f dwgsim-${PACKAGE_VERSION}.tar ]; then \
        rm -rv dwgsim-${PACKAGE_VERSION}.tar; \
	fi; \
	if [ -d dwgsim-${PACKAGE_VERSION} ]; then \
        rm -rv dwgsim-${PACKAGE_VERSION}; \
	fi; \
    mkdir dwgsim-${PACKAGE_VERSION}; \
	cp -r INSTALL LICENSE Makefile README scripts src dwgsim-${PACKAGE_VERSION}/.; \
	tar -vcf dwgsim-${PACKAGE_VERSION}.tar dwgsim-${PACKAGE_VERSION}; \
	gzip -9 dwgsim-${PACKAGE_VERSION}.tar; \
	rm -rv dwgsim-${PACKAGE_VERSION};

# Run all tests
test: test-unit test-integration test-bgzf test-parallel-wgs test-matched-wgs test-bed

# Integration tests
test-integration: dwgsim
	if [ -d tmp ]; then rm -r tmp; fi
	/bin/bash testdata/test.sh

# BGZF output tests
test-bgzf: dwgsim
	$(MAKE) -C samtools bgzip
	DWGSIM_BIN="$(DWGSIM_BIN)" BGZIP_BIN="./samtools/bgzip" \
		/bin/bash tests/test_bgzf_output.sh "$(BENCHMARK_REFERENCE)"

# Fixed-task generation, pairing, and byte determinism tests
test-parallel-wgs: dwgsim samtools-program
	$(MAKE) -C samtools bgzip
	DWGSIM_BIN="$(DWGSIM_BIN)" SAMTOOLS_BIN="./samtools/samtools" BGZIP_BIN="./samtools/bgzip" \
		/bin/bash tests/test_parallel_wgs.sh "$(BENCHMARK_REFERENCE)"

# Sparse germline/somatic matched and tumor-only determinism tests
test-matched-wgs: dwgsim samtools-program
	$(MAKE) -C samtools bgzip
	DWGSIM_BIN="$(DWGSIM_BIN)" SAMTOOLS_BIN="./samtools/samtools" BGZIP_BIN="./samtools/bgzip" \
		/bin/bash tests/test_matched_wgs.sh "$(BENCHMARK_REFERENCE)"

# Regions BED parser and placement tests
test-bed: dwgsim
	DWGSIM_BIN="$(DWGSIM_BIN)" /bin/bash tests/test_regions_bed.sh

# Unit test target
TEST_OBJS = tests/test_main.o src/fastq_writer.o src/regions_bed.o src/contigs.o
TEST_PROG = tests/run_tests

tests/test_main.o: tests/test_main.c tests/test_framework.h src/fastq_writer.h src/regions_bed.h src/contigs.h
	$(CC) -c $(CFLAGS) $(DFLAGS) -I. -Isrc tests/test_main.c -o $@

$(TEST_PROG): lib-recur $(TEST_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_OBJS) samtools/libbam.a -lm -lz -lpthread

test-unit: $(TEST_PROG)
	./$(TEST_PROG)

clean-tests:
	rm -f tests/*.o $(TEST_PROG)
