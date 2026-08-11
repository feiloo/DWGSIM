.DEFAULT_GOAL := all

PACKAGE_VERSION=0.1.17-dev
CC=			gcc
CFLAGS=		-g -Wall -O3 #-m64 #-arch ppc
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -DPACKAGE_VERSION=\"$(PACKAGE_VERSION)\"
DWGSIM_AOBJS = src/dwgsim_opt.o src/mut.o src/contigs.o src/regions_bed.o \
			   src/mut_txt.o src/mut_bed.o src/mut_vcf.o src/mut_input.o src/fastq_writer.o src/dwgsim.o
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
REFERENCE_ROOT ?= reference
HUMAN_REFERENCE_RELEASE ?= GCF_000001405.40_GRCh38.p14
HUMAN_REFERENCE_BASENAME ?= $(HUMAN_REFERENCE_RELEASE)_genomic
HUMAN_REFERENCE_DIR ?= $(REFERENCE_ROOT)/GRCh38.p14
HUMAN_REFERENCE_ARCHIVE ?= $(HUMAN_REFERENCE_DIR)/$(HUMAN_REFERENCE_BASENAME).fna.gz
HUMAN_REFERENCE_FASTA ?= $(HUMAN_REFERENCE_DIR)/$(HUMAN_REFERENCE_BASENAME).fna
HUMAN_REFERENCE_URL ?= https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/$(HUMAN_REFERENCE_BASENAME).fna.gz
HUMAN_REFERENCE_MD5 ?= c30471567037b2b2389d43c908c653e1
HUMAN_SMOKE_DIR ?= build/human-reference-smoke
HUMAN_SMOKE_READ_PAIRS ?= 100
HUMAN_SMOKE_SEED ?= 13
DWGSIM_BIN ?= ./dwgsim
TIME_BIN ?= /usr/bin/time
BENCHMARK_REFERENCE ?= samtools/examples/ex1.fa
BENCHMARK_DIR ?= build/benchmark
BENCHMARK_READ_PAIRS ?= 250000
BENCHMARK_READ_LENGTH_1 ?= 100
BENCHMARK_READ_LENGTH_2 ?= 100
BENCHMARK_SEED ?= 13
BENCHMARK_THREADS ?= 1

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

all:$(PROG)

.PHONY:all lib clean cleanlocal test test-unit test-integration test-bgzf clean-tests help
.PHONY:download download-human-reference test-human-reference
.PHONY:benchmark
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

help:
	@printf 'Usage: make <target> [VARIABLE=value]\n\n'
	@printf 'Build and test targets:\n'
	@printf '  %-26s %s\n' 'all' 'Build all DWGSIM executables (default).'
	@printf '  %-26s %s\n' 'test' 'Run unit, integration, and BGZF output tests.'
	@printf '  %-26s %s\n' 'benchmark' 'Measure reads-only simulation throughput and resource use.'
	@printf '  %-26s %s\n' 'test-bgzf' 'Test BGZF compatibility, modes, and thread determinism.'
	@printf '  %-26s %s\n' 'download' 'Download and verify the full NCBI GRCh38.p14 reference.'
	@printf '  %-26s %s\n' 'test-human-reference' 'Run the reads-only DWGSIM smoke test on full GRCh38.p14.'
	@printf '  %-26s %s\n' 'clean' 'Remove compiled programs, objects, and regular test artifacts.'
	@printf '  %-26s %s\n' 'help' 'Show this help.'
	@printf '\nBenchmark settings:\n'
	@printf '  %-30s %s\n' 'BENCHMARK_REFERENCE=...' 'Input FASTA (default: $(BENCHMARK_REFERENCE)).'
	@printf '  %-30s %s\n' 'BENCHMARK_READ_PAIRS=...' 'Read pairs to generate (default: $(BENCHMARK_READ_PAIRS)).'
	@printf '  %-30s %s\n' 'BENCHMARK_READ_LENGTH_1=...' 'First-read length (default: $(BENCHMARK_READ_LENGTH_1)).'
	@printf '  %-30s %s\n' 'BENCHMARK_READ_LENGTH_2=...' 'Second-read length (default: $(BENCHMARK_READ_LENGTH_2)).'
	@printf '  %-30s %s\n' 'BENCHMARK_SEED=...' 'Random seed (default: $(BENCHMARK_SEED)).'
	@printf '  %-30s %s\n' 'BENCHMARK_THREADS=...' 'DWGSIM thread budget (default: $(BENCHMARK_THREADS)).'
	@printf '  %-30s %s\n' 'BENCHMARK_DIR=...' 'Output and report directory (default: $(BENCHMARK_DIR)).'
	@printf '\nHuman-reference settings:\n'
	@printf '  %-26s %s\n' 'HUMAN_REFERENCE_DIR=...' 'Reference destination (default: $(HUMAN_REFERENCE_DIR)).'
	@printf '  %-26s %s\n' 'HUMAN_SMOKE_READ_PAIRS=...' 'Read pairs to simulate (default: $(HUMAN_SMOKE_READ_PAIRS)).'
	@printf '  %-26s %s\n' 'HUMAN_SMOKE_DIR=...' 'Smoke-test output directory (default: $(HUMAN_SMOKE_DIR)).'

download: download-human-reference

download-human-reference:
	CURL="$(CURL)" MD5SUM="$(MD5SUM)" /bin/bash scripts/download_human_reference.sh \
		"$(HUMAN_REFERENCE_URL)" "$(HUMAN_REFERENCE_MD5)" \
		"$(HUMAN_REFERENCE_ARCHIVE)" "$(HUMAN_REFERENCE_FASTA)"

test-human-reference: all
	$(MAKE) --no-print-directory download-human-reference
	DWGSIM_BIN="$(DWGSIM_BIN)" \
		HUMAN_SMOKE_READ_PAIRS="$(HUMAN_SMOKE_READ_PAIRS)" \
		HUMAN_SMOKE_SEED="$(HUMAN_SMOKE_SEED)" \
		/bin/bash scripts/smoke_test_human_reference.sh \
		"$(HUMAN_REFERENCE_FASTA)" "$(HUMAN_SMOKE_DIR)"

benchmark: dwgsim
	DWGSIM_BIN="$(DWGSIM_BIN)" \
		TIME_BIN="$(TIME_BIN)" \
		BENCHMARK_READ_PAIRS="$(BENCHMARK_READ_PAIRS)" \
		BENCHMARK_READ_LENGTH_1="$(BENCHMARK_READ_LENGTH_1)" \
		BENCHMARK_READ_LENGTH_2="$(BENCHMARK_READ_LENGTH_2)" \
		BENCHMARK_SEED="$(BENCHMARK_SEED)" \
		BENCHMARK_THREADS="$(BENCHMARK_THREADS)" \
		/bin/bash scripts/benchmark.sh \
		"$(BENCHMARK_REFERENCE)" "$(BENCHMARK_DIR)"

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
test: test-unit test-integration test-bgzf

# Integration tests
test-integration: dwgsim
	if [ -d tmp ]; then rm -r tmp; fi
	/bin/bash testdata/test.sh

# BGZF output tests
test-bgzf: dwgsim
	$(MAKE) -C samtools bgzip
	DWGSIM_BIN="$(DWGSIM_BIN)" BGZIP_BIN="./samtools/bgzip" \
		/bin/bash tests/test_bgzf_output.sh "$(BENCHMARK_REFERENCE)"

# Unit test target
TEST_OBJS = tests/test_main.o src/fastq_writer.o
TEST_PROG = tests/run_tests

tests/test_main.o: tests/test_main.c tests/test_framework.h src/fastq_writer.h
	$(CC) -c $(CFLAGS) $(DFLAGS) -I. -Isrc tests/test_main.c -o $@

$(TEST_PROG): lib-recur $(TEST_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_OBJS) samtools/libbam.a -lm -lz -lpthread

test-unit: $(TEST_PROG)
	./$(TEST_PROG)

clean-tests:
	rm -f tests/*.o $(TEST_PROG)
