.DEFAULT_GOAL := all

PACKAGE_VERSION="0.1.17-dev"
CC=			gcc
CFLAGS=		-g -Wall -O3 #-m64 #-arch ppc
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -DPACKAGE_VERSION=\\\"${PACKAGE_VERSION}\\\"
DWGSIM_AOBJS = src/dwgsim_opt.o src/mut.o src/contigs.o src/regions_bed.o \
			   src/mut_txt.o src/mut_bed.o src/mut_vcf.o src/mut_input.o src/dwgsim.o
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

GATK ?= gatk
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
HUMAN_TEST_DIR ?= build/human-reference-test
HUMAN_TEST_READ_PAIRS ?= 1000
HUMAN_TEST_MUTATION_RATE ?= 0.000001
HUMAN_TEST_SEED ?= 13
DWGSIM_BIN ?= ./dwgsim
SAMTOOLS_BIN ?= samtools/samtools

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

all:$(PROG)

.PHONY:all lib clean cleanlocal test test-unit test-integration clean-tests help
.PHONY:download download-human-reference test-human-reference
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

help:
	@printf 'Usage: make <target> [VARIABLE=value]\n\n'
	@printf 'Build and test targets:\n'
	@printf '  %-26s %s\n' 'all' 'Build all DWGSIM executables (default).'
	@printf '  %-26s %s\n' 'test' 'Run unit and integration tests.'
	@printf '  %-26s %s\n' 'download' 'Download and verify the full NCBI GRCh38.p14 reference.'
	@printf '  %-26s %s\n' 'test-human-reference' 'Run DWGSIM on GRCh38.p14 and validate its VCF with GATK.'
	@printf '  %-26s %s\n' 'clean' 'Remove compiled programs, objects, and regular test artifacts.'
	@printf '  %-26s %s\n' 'help' 'Show this help.'
	@printf '\nHuman-reference settings:\n'
	@printf '  %-26s %s\n' 'GATK=gatk' 'GATK launcher or absolute executable path.'
	@printf '  %-26s %s\n' 'HUMAN_REFERENCE_DIR=...' 'Reference destination (default: $(HUMAN_REFERENCE_DIR)).'
	@printf '  %-26s %s\n' 'HUMAN_TEST_READ_PAIRS=...' 'Read pairs to simulate (default: $(HUMAN_TEST_READ_PAIRS)).'
	@printf '  %-26s %s\n' 'HUMAN_TEST_MUTATION_RATE=...' 'Mutation rate (default: $(HUMAN_TEST_MUTATION_RATE)).'
	@printf '  %-26s %s\n' 'HUMAN_TEST_DIR=...' 'Test output directory (default: $(HUMAN_TEST_DIR)).'

download: download-human-reference

download-human-reference:
	CURL="$(CURL)" MD5SUM="$(MD5SUM)" /bin/bash scripts/download_human_reference.sh \
		"$(HUMAN_REFERENCE_URL)" "$(HUMAN_REFERENCE_MD5)" \
		"$(HUMAN_REFERENCE_ARCHIVE)" "$(HUMAN_REFERENCE_FASTA)"

test-human-reference: all
	@if ! command -v "$(GATK)" >/dev/null 2>&1; then \
		echo "GATK was not found: $(GATK)" >&2; \
		echo "Install GATK or run make test-human-reference GATK=/path/to/gatk" >&2; \
		exit 1; \
	fi
	$(MAKE) --no-print-directory download-human-reference
	GATK="$(GATK)" DWGSIM_BIN="$(DWGSIM_BIN)" SAMTOOLS_BIN="$(SAMTOOLS_BIN)" \
		HUMAN_TEST_READ_PAIRS="$(HUMAN_TEST_READ_PAIRS)" \
		HUMAN_TEST_MUTATION_RATE="$(HUMAN_TEST_MUTATION_RATE)" \
		HUMAN_TEST_SEED="$(HUMAN_TEST_SEED)" \
		/bin/bash scripts/test_human_reference.sh \
		"$(HUMAN_REFERENCE_FASTA)" "$(HUMAN_TEST_DIR)"

dwgsim:lib-recur $(DWGSIM_AOBJS)
	$(CC) $(CFLAGS) -o $@ $(DWGSIM_AOBJS) -lm -lz -lpthread

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

# Run all tests (unit + integration)
test: test-unit test-integration

# Integration tests
test-integration:
	if [ -d tmp ]; then rm -r tmp; fi
	/bin/bash testdata/test.sh

# Unit test target
TEST_OBJS = tests/test_main.o
TEST_PROG = tests/run_tests

tests/test_main.o: tests/test_main.c tests/test_framework.h
	$(CC) -c $(CFLAGS) $(DFLAGS) -I. -Isrc tests/test_main.c -o $@

$(TEST_PROG): $(TEST_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_OBJS) -lm

test-unit: $(TEST_PROG)
	./$(TEST_PROG)

clean-tests:
	rm -f tests/*.o $(TEST_PROG)
