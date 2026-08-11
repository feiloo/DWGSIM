/*
 * DWGSIM Unit Test Runner
 *
 * This file runs all unit tests for DWGSIM.
 * Individual test files are included below.
 */

#include "test_framework.h"
#include "fastq_writer.h"
#include "contigs.h"
#include "regions_bed.h"

/* Include test files here as they are added */
/* #include "test_position.c" */
/* #include "test_memory.c" */
/* #include "test_parsing.c" */

/*
 * Placeholder tests to verify the test framework works
 */
TEST(test_framework_assert) {
    ASSERT(1 == 1);
    ASSERT(0 == 0);
}

TEST(test_framework_assert_eq) {
    ASSERT_EQ(42, 42);
    ASSERT_EQ(-1, -1);
    ASSERT_EQ(0, 0);
}

TEST(test_framework_assert_ne) {
    ASSERT_NE(1, 2);
    ASSERT_NE(-1, 1);
}

TEST(test_framework_assert_str_eq) {
    ASSERT_STR_EQ("hello", "hello");
    ASSERT_STR_EQ("", "");
}

TEST(test_framework_assert_float_eq) {
    ASSERT_FLOAT_EQ(3.14, 3.14, 0.001);
    ASSERT_FLOAT_EQ(0.1 + 0.2, 0.3, 0.0001);
}

TEST(test_framework_assert_comparisons) {
    ASSERT_LT(1, 2);
    ASSERT_LE(1, 1);
    ASSERT_LE(1, 2);
    ASSERT_GT(2, 1);
    ASSERT_GE(2, 2);
    ASSERT_GE(3, 2);
}

TEST(test_framework_assert_null) {
    int *ptr = NULL;
    int val = 42;
    ASSERT_NULL(ptr);
    ASSERT_NOT_NULL(&val);
}

TEST(fastq_writer_thread_distribution) {
    ASSERT_EQ(1, fastq_writer_threads_for_stream(1, 2, 0));
    ASSERT_EQ(1, fastq_writer_threads_for_stream(1, 2, 1));

    ASSERT_EQ(3, fastq_writer_threads_for_stream(4, 2, 0));
    ASSERT_EQ(2, fastq_writer_threads_for_stream(4, 2, 1));

    ASSERT_EQ(2, fastq_writer_threads_for_stream(4, 3, 0));
    ASSERT_EQ(2, fastq_writer_threads_for_stream(4, 3, 1));
    ASSERT_EQ(2, fastq_writer_threads_for_stream(4, 3, 2));

    ASSERT_EQ(0, fastq_writer_threads_for_stream(0, 2, 0));
    ASSERT_EQ(0, fastq_writer_threads_for_stream(1, 0, 0));
    ASSERT_EQ(0, fastq_writer_threads_for_stream(1, 2, 2));
}

TEST(regions_bed_half_open_boundaries) {
    contigs_t *contigs = contigs_init();
    regions_bed_txt *regions;
    FILE *bed = tmpfile();
    uint32_t position = UINT32_MAX;

    ASSERT_NOT_NULL(contigs);
    ASSERT_NOT_NULL(bed);
    contigs_add(contigs, "chr1", 100);
    contigs_add(contigs, "chr2", 100);

    fputs("# BED comments are accepted\n", bed);
    fputs("track name=targets\n", bed);
    fputs("browser position chr1:1-100\n", bed);
    fputs("chr1\t0\t10\tfirst\n", bed);
    fputs("chr1\t10\t20\tadjacent\n", bed);
    fputs("chr1\t30\t40\n", bed);
    fputs("chr2\t5\t15\n", bed);
    rewind(bed);

    regions = regions_bed_init(bed, contigs);
    ASSERT_NOT_NULL(regions);
    ASSERT_EQ(3, regions->n);
    ASSERT_EQ(0, regions->start[0]);
    ASSERT_EQ(20, regions->end[0]);
    ASSERT_EQ(30, regions->start[1]);
    ASSERT_EQ(40, regions->end[1]);

    ASSERT_EQ(1, regions_bed_query(regions, 0, 0, 20));
    ASSERT_EQ(0, regions_bed_query(regions, 0, 0, 21));
    ASSERT_EQ(0, regions_bed_query(regions, 0, 20, 30));
    ASSERT_EQ(1, regions_bed_query(regions, 0, 30, 40));

    ASSERT_EQ(1, regions_bed_map_offset(regions, 0, 0, &position));
    ASSERT_EQ(0, position);
    ASSERT_EQ(1, regions_bed_map_offset(regions, 0, 19, &position));
    ASSERT_EQ(19, position);
    ASSERT_EQ(1, regions_bed_map_offset(regions, 0, 20, &position));
    ASSERT_EQ(30, position);
    ASSERT_EQ(1, regions_bed_map_offset(regions, 0, 29, &position));
    ASSERT_EQ(39, position);
    ASSERT_EQ(0, regions_bed_map_offset(regions, 0, 30, &position));
    ASSERT_EQ(1, regions_bed_map_offset(regions, 1, 0, &position));
    ASSERT_EQ(5, position);

    regions_bed_destroy(regions);
    fclose(bed);
    contigs_destroy(contigs);
}

int main(void) {
    printf("DWGSIM Unit Tests\n");
    printf("========================================\n");

    TEST_SUITE("Test Framework Verification");
    RUN_TEST(test_framework_assert);
    RUN_TEST(test_framework_assert_eq);
    RUN_TEST(test_framework_assert_ne);
    RUN_TEST(test_framework_assert_str_eq);
    RUN_TEST(test_framework_assert_float_eq);
    RUN_TEST(test_framework_assert_comparisons);
    RUN_TEST(test_framework_assert_null);

    TEST_SUITE("BGZF FASTQ writer");
    RUN_TEST(fastq_writer_thread_distribution);

    TEST_SUITE("Regions BED");
    RUN_TEST(regions_bed_half_open_boundaries);

    /* Add test suites here as they are created */
    /* TEST_SUITE("Position Calculation Tests"); */
    /* Include tests from test_position.c */

    TEST_SUMMARY();

    return test_failures > 0 ? 1 : 0;
}
