#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <fcntl.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>

#include "parallel_wgs.h"
#include "samtools/bgzf.h"

#define PARALLEL_FORMAT_VERSION 1
#define PAIRS_PER_TASK 8192U
#define MAX_FRAGMENT_ATTEMPTS 10000U
#define REFERENCE_READ_BUFFER (1024U * 1024U)
#define QUAL_MAX 40

#define DOMAIN_FRAGMENT UINT64_C(0x465241474d454e54)
#define DOMAIN_STRAND   UINT64_C(0x535452414e442020)
#define DOMAIN_ERROR    UINT64_C(0x4552524f52202020)
#define DOMAIN_QUALITY  UINT64_C(0x5155414c49545920)

typedef struct {
    unsigned char *data;
    size_t length;
    size_t capacity;
} byte_buffer_t;

typedef struct {
    char *name;
    uint64_t length;
    uint64_t offset;
    uint32_t line_bases;
    uint32_t line_width;
    uint8_t *bases;
    uint64_t pair_count;
    int eligible;
} reference_contig_t;

typedef struct {
    reference_contig_t *contigs;
    size_t count;
    uint64_t total_length;
} reference_manifest_t;

typedef struct {
    uint64_t id;
    uint32_t contig_index;
    uint64_t first_pair;
    uint32_t pair_count;
} pair_task_t;

typedef struct {
    pair_task_t *tasks;
    size_t count;
} task_manifest_t;

typedef struct {
    uint64_t state;
    double saved_normal;
    int has_saved_normal;
} local_rng_t;

typedef struct {
    uint64_t task_id;
    uint32_t pair_count;
    byte_buffer_t compressed[2];
} paired_task_result_t;

typedef struct {
    paired_task_result_t *result;
    unsigned consumed;
} result_slot_t;

typedef struct {
    const dwgsim_opt_t *opt;
    const reference_manifest_t *reference;
    const task_manifest_t *tasks;
    result_slot_t *slots;
    size_t window;
    uint64_t next_claim;
    uint64_t released;
    int failed;
    char error_message[256];
    pthread_mutex_t mutex;
    pthread_cond_t changed;
} pipeline_t;

typedef struct {
    pipeline_t *pipeline;
    int mate;
    FILE *output;
    int status;
} appender_arg_t;

static const unsigned char bgzf_header[18] = {
    31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0, 0, 0
};

static const unsigned char bgzf_eof[28] = {
    31, 139, 8, 4, 0, 0, 0, 0, 0, 255, 6, 0, 66, 67, 2, 0,
    27, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

static int
buffer_reserve(byte_buffer_t *buffer, size_t additional)
{
    size_t required;
    size_t capacity;
    unsigned char *data;

    if(additional > SIZE_MAX - buffer->length) return -1;
    required = buffer->length + additional;
    if(required <= buffer->capacity) return 0;

    capacity = buffer->capacity ? buffer->capacity : 4096;
    while(capacity < required) {
        if(capacity > SIZE_MAX / 2) {
            capacity = required;
            break;
        }
        capacity *= 2;
    }

    data = realloc(buffer->data, capacity);
    if(NULL == data) return -1;
    buffer->data = data;
    buffer->capacity = capacity;
    return 0;
}

static int
buffer_append(byte_buffer_t *buffer, const void *data, size_t length)
{
    if(0 == length) return 0;
    if(NULL == data || 0 != buffer_reserve(buffer, length)) return -1;
    memcpy(buffer->data + buffer->length, data, length);
    buffer->length += length;
    return 0;
}

static int
buffer_appendf(byte_buffer_t *buffer, const char *format, ...)
{
    va_list args;
    va_list copy;
    int needed;

    va_start(args, format);
    va_copy(copy, args);
    needed = vsnprintf(NULL, 0, format, args);
    va_end(args);
    if(needed < 0 || 0 != buffer_reserve(buffer, (size_t)needed + 1)) {
        va_end(copy);
        return -1;
    }
    if(vsnprintf((char *)buffer->data + buffer->length,
                 (size_t)needed + 1, format, copy) != needed) {
        va_end(copy);
        return -1;
    }
    va_end(copy);
    buffer->length += (size_t)needed;
    return 0;
}

static void
buffer_destroy(byte_buffer_t *buffer)
{
    free(buffer->data);
    memset(buffer, 0, sizeof(*buffer));
}

static uint64_t
mix64(uint64_t value)
{
    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return value;
}

static void
rng_init(local_rng_t *rng, uint64_t seed, uint64_t domain,
         uint64_t contig, uint64_t pair, uint64_t extra)
{
    uint64_t key = mix64(seed ^ UINT64_C(0x44574753494d0001));
    key = mix64(key ^ domain);
    key = mix64(key ^ contig);
    key = mix64(key ^ pair);
    rng->state = mix64(key ^ extra);
    rng->saved_normal = 0.0;
    rng->has_saved_normal = 0;
}

static uint64_t
rng_next(local_rng_t *rng)
{
    rng->state += UINT64_C(0x9e3779b97f4a7c15);
    return mix64(rng->state);
}

static double
rng_uniform(local_rng_t *rng)
{
    return (double)(rng_next(rng) >> 11) * 0x1.0p-53;
}

static uint64_t
rng_bounded(local_rng_t *rng, uint64_t bound)
{
    uint64_t threshold;
    uint64_t value;

    if(bound <= 1) return 0;
    threshold = (uint64_t)(-bound) % bound;
    do {
        value = rng_next(rng);
    } while(value < threshold);
    return value % bound;
}

static double
rng_normal(local_rng_t *rng)
{
    double factor;
    double radius;
    double value1;
    double value2;

    if(rng->has_saved_normal) {
        rng->has_saved_normal = 0;
        return rng->saved_normal;
    }

    do {
        value1 = 2.0 * rng_uniform(rng) - 1.0;
        value2 = 2.0 * rng_uniform(rng) - 1.0;
        radius = value1 * value1 + value2 * value2;
    } while(1.0 <= radius || 0.0 == radius);

    factor = sqrt(-2.0 * log(radius) / radius);
    rng->saved_normal = value1 * factor;
    rng->has_saved_normal = 1;
    return value2 * factor;
}

static char *
path_with_suffix(const char *path, const char *suffix)
{
    size_t path_length;
    size_t suffix_length;
    char *result;

    if(NULL == path || NULL == suffix) return NULL;
    path_length = strlen(path);
    suffix_length = strlen(suffix);
    if(path_length > SIZE_MAX - suffix_length - 1) return NULL;
    result = malloc(path_length + suffix_length + 1);
    if(NULL == result) return NULL;
    memcpy(result, path, path_length);
    memcpy(result + path_length, suffix, suffix_length + 1);
    return result;
}

static int
parse_uint64_field(const char *text, uint64_t *value)
{
    char *end = NULL;
    unsigned long long parsed;

    if(NULL == text || '\0' == *text || '-' == *text) return -1;
    errno = 0;
    parsed = strtoull(text, &end, 10);
    if(ERANGE == errno || end == text ||
       ('\0' != *end && '\n' != *end && '\r' != *end)) return -1;
    *value = (uint64_t)parsed;
    return 0;
}

static int
manifest_append_contig(reference_manifest_t *manifest,
                       const reference_contig_t *contig,
                       size_t *capacity)
{
    reference_contig_t *contigs;
    size_t next_capacity;

    if(manifest->count == *capacity) {
        next_capacity = *capacity ? *capacity * 2 : 64;
        if(next_capacity < *capacity ||
           next_capacity > SIZE_MAX / sizeof(*contigs)) return -1;
        contigs = realloc(manifest->contigs,
                          next_capacity * sizeof(*contigs));
        if(NULL == contigs) return -1;
        manifest->contigs = contigs;
        *capacity = next_capacity;
    }
    manifest->contigs[manifest->count++] = *contig;
    return 0;
}

static int
read_reference_manifest(const char *reference_path,
                        reference_manifest_t *manifest)
{
    char *fai_path = path_with_suffix(reference_path, ".fai");
    char *line = NULL;
    size_t line_capacity = 0;
    size_t contig_capacity = 0;
    ssize_t line_length;
    FILE *fai = NULL;
    int status = -1;

    if(NULL == fai_path) goto cleanup;
    fai = fopen(fai_path, "r");
    if(NULL == fai) {
        fprintf(stderr, "[dwgsim_parallel] cannot open FASTA index %s: %s\n",
                fai_path, strerror(errno));
        goto cleanup;
    }

    while(0 <= (line_length = getline(&line, &line_capacity, fai))) {
        reference_contig_t contig;
        char *save = NULL;
        char *fields[5];
        uint64_t line_bases;
        uint64_t line_width;
        int field;

        if(0 == line_length) continue;
        memset(&contig, 0, sizeof(contig));
        for(field = 0; field < 5; ++field) {
            fields[field] = strtok_r(0 == field ? line : NULL,
                                     "\t\n\r", &save);
            if(NULL == fields[field]) break;
        }
        if(5 != field || '\0' == fields[0][0] ||
           0 != parse_uint64_field(fields[1], &contig.length) ||
           0 != parse_uint64_field(fields[2], &contig.offset) ||
           0 != parse_uint64_field(fields[3], &line_bases) ||
           0 != parse_uint64_field(fields[4], &line_width) ||
           0 == contig.length || 0 == line_bases ||
           UINT32_MAX < line_bases || UINT32_MAX < line_width ||
           line_width < line_bases) {
            fprintf(stderr, "[dwgsim_parallel] malformed FASTA index line\n");
            goto cleanup;
        }
        contig.name = strdup(fields[0]);
        if(NULL == contig.name) goto cleanup;
        contig.line_bases = (uint32_t)line_bases;
        contig.line_width = (uint32_t)line_width;
        if(UINT64_MAX - manifest->total_length < contig.length ||
           0 != manifest_append_contig(manifest, &contig,
                                       &contig_capacity)) {
            free(contig.name);
            goto cleanup;
        }
        manifest->total_length += contig.length;
    }
    if(ferror(fai) || 0 == manifest->count) goto cleanup;
    status = 0;

cleanup:
    if(NULL != fai) fclose(fai);
    free(line);
    free(fai_path);
    return status;
}

static void
destroy_reference_manifest(reference_manifest_t *manifest)
{
    size_t index;

    for(index = 0; index < manifest->count; ++index) {
        free(manifest->contigs[index].name);
        free(manifest->contigs[index].bases);
    }
    free(manifest->contigs);
    memset(manifest, 0, sizeof(*manifest));
}

static uint8_t
encode_base(unsigned char base)
{
    switch(base) {
      case 'A': case 'a': return 0;
      case 'C': case 'c': return 1;
      case 'G': case 'g': return 2;
      case 'T': case 't': return 3;
      default: return 4;
    }
}

static int
load_reference_contig(int reference_fd, reference_contig_t *contig)
{
    unsigned char *input = NULL;
    uint64_t last_base;
    uint64_t physical_length;
    uint64_t physical_position = 0;
    uint64_t sequence_position = 0;
    uint64_t non_acgt = 0;
    double non_acgt_fraction;
    int status = -1;

    if(!contig->eligible) return 0;
    if(contig->length > SIZE_MAX) return -1;
    last_base = contig->length - 1;
    if(last_base / contig->line_bases >
       (UINT64_MAX - (last_base % contig->line_bases) - 1) /
       contig->line_width) return -1;
    physical_length = (last_base / contig->line_bases) *
                      contig->line_width +
                      (last_base % contig->line_bases) + 1;

    contig->bases = malloc((size_t)contig->length);
    input = malloc(REFERENCE_READ_BUFFER);
    if(NULL == contig->bases || NULL == input) goto cleanup;

    while(physical_position < physical_length) {
        uint64_t remaining = physical_length - physical_position;
        size_t requested = remaining < REFERENCE_READ_BUFFER ?
                           (size_t)remaining : REFERENCE_READ_BUFFER;
        ssize_t count;
        size_t index;

        do {
            count = pread(reference_fd, input, requested,
                          (off_t)(contig->offset + physical_position));
        } while(count < 0 && EINTR == errno);
        if(count <= 0) goto cleanup;
        for(index = 0; index < (size_t)count &&
                       sequence_position < contig->length; ++index) {
            uint8_t encoded;
            if('\n' == input[index] || '\r' == input[index]) continue;
            encoded = encode_base(input[index]);
            contig->bases[sequence_position++] = encoded;
            if(4 == encoded) non_acgt++;
        }
        physical_position += (uint64_t)count;
    }

    if(sequence_position != contig->length) goto cleanup;
    non_acgt_fraction = (double)non_acgt / (double)contig->length;
    if(0.95 < non_acgt_fraction) {
        fprintf(stderr,
                "[dwgsim_parallel] skipping %s: %.2f%% non-ACGT bases\n",
                contig->name, 100.0 * non_acgt_fraction);
        free(contig->bases);
        contig->bases = NULL;
        contig->eligible = 0;
    }
    status = 0;

cleanup:
    if(0 != status) {
        free(contig->bases);
        contig->bases = NULL;
    }
    free(input);
    return status;
}

typedef struct {
    reference_manifest_t *manifest;
    int reference_fd;
    volatile size_t next_contig;
    volatile int failed;
} reference_loader_t;

static void *
reference_loader_worker(void *data)
{
    reference_loader_t *loader = data;

    for(;;) {
        size_t index = __sync_fetch_and_add(&loader->next_contig, 1);
        if(loader->manifest->count <= index) break;
        if(__sync_fetch_and_add(&loader->failed, 0)) continue;
        if(0 != load_reference_contig(loader->reference_fd,
                                     &loader->manifest->contigs[index])) {
            fprintf(stderr, "[dwgsim_parallel] failed to load contig %s\n",
                    loader->manifest->contigs[index].name);
            __sync_lock_test_and_set(&loader->failed, 1);
        }
    }
    return NULL;
}

static int
prepare_reference(const dwgsim_opt_t *opt, const char *reference_path,
                  reference_manifest_t *manifest)
{
    reference_loader_t loader;
    pthread_t *threads = NULL;
    size_t worker_count;
    size_t created = 0;
    size_t index;
    int reference_fd;
    double minimum_length = opt->dist + 3.0 * opt->std_dev;
    int32_t maximum_read = opt->length[0] > opt->length[1] ?
                           opt->length[0] : opt->length[1];

    for(index = 0; index < manifest->count; ++index) {
        reference_contig_t *contig = &manifest->contigs[index];
        contig->eligible = contig->length <= UINT32_MAX &&
                           contig->length >= (uint64_t)maximum_read &&
                           (double)contig->length >= minimum_length;
    }

    reference_fd = open(reference_path, O_RDONLY);
    if(reference_fd < 0) {
        fprintf(stderr, "[dwgsim_parallel] cannot open %s: %s\n",
                reference_path, strerror(errno));
        return -1;
    }

    memset(&loader, 0, sizeof(loader));
    loader.manifest = manifest;
    loader.reference_fd = reference_fd;
    worker_count = (size_t)opt->compression_threads;
    if(manifest->count < worker_count) worker_count = manifest->count;
    if(0 == worker_count) worker_count = 1;

    if(1 == worker_count) {
        reference_loader_worker(&loader);
    }
    else {
        threads = calloc(worker_count, sizeof(*threads));
        if(NULL == threads) {
            __sync_lock_test_and_set(&loader.failed, 1);
        }
        else {
            while(created < worker_count) {
                if(0 != pthread_create(&threads[created], NULL,
                                       reference_loader_worker, &loader)) {
                    fprintf(stderr,
                            "[dwgsim_parallel] failed to create reference worker\n");
                    __sync_lock_test_and_set(&loader.failed, 1);
                    break;
                }
                created++;
            }
        }
        for(index = 0; index < created; ++index) {
            pthread_join(threads[index], NULL);
        }
    }

    free(threads);
    close(reference_fd);
    return __sync_fetch_and_add(&loader.failed, 0) ? -1 : 0;
}

typedef struct {
    size_t contig_index;
    long double remainder;
} allocation_remainder_t;

static int
compare_remainders(const void *left, const void *right)
{
    const allocation_remainder_t *a = left;
    const allocation_remainder_t *b = right;

    if(a->remainder > b->remainder) return -1;
    if(a->remainder < b->remainder) return 1;
    if(a->contig_index < b->contig_index) return -1;
    if(a->contig_index > b->contig_index) return 1;
    return 0;
}

static int
allocate_pairs(reference_manifest_t *manifest, uint64_t requested_pairs)
{
    allocation_remainder_t *remainders = NULL;
    uint64_t eligible_length = 0;
    uint64_t assigned = 0;
    size_t eligible_count = 0;
    size_t index;

    for(index = 0; index < manifest->count; ++index) {
        reference_contig_t *contig = &manifest->contigs[index];
        if(!contig->eligible) continue;
        if(UINT64_MAX - eligible_length < contig->length) return -1;
        eligible_length += contig->length;
        eligible_count++;
    }
    if(0 == eligible_length || 0 == eligible_count) {
        fprintf(stderr, "[dwgsim_parallel] reference has no eligible contigs\n");
        return -1;
    }

    remainders = calloc(eligible_count, sizeof(*remainders));
    if(NULL == remainders) return -1;
    eligible_count = 0;
    for(index = 0; index < manifest->count; ++index) {
        reference_contig_t *contig = &manifest->contigs[index];
        long double exact;
        uint64_t integral;

        if(!contig->eligible) continue;
        exact = (long double)requested_pairs * contig->length /
                eligible_length;
        integral = (uint64_t)floorl(exact);
        contig->pair_count = integral;
        assigned += integral;
        remainders[eligible_count].contig_index = index;
        remainders[eligible_count].remainder = exact - integral;
        eligible_count++;
    }
    if(assigned > requested_pairs ||
       requested_pairs - assigned > eligible_count) {
        free(remainders);
        return -1;
    }
    qsort(remainders, eligible_count, sizeof(*remainders),
          compare_remainders);
    for(index = 0; index < (size_t)(requested_pairs - assigned); ++index) {
        manifest->contigs[remainders[index].contig_index].pair_count++;
    }
    free(remainders);

    for(index = 0; index < manifest->count; ++index) {
        reference_contig_t *contig = &manifest->contigs[index];
        if(0 == contig->pair_count) {
            free(contig->bases);
            contig->bases = NULL;
        }
    }
    return 0;
}

static int
build_tasks(const reference_manifest_t *reference,
            task_manifest_t *manifest)
{
    size_t task_count = 0;
    size_t task_index = 0;
    size_t contig_index;

    for(contig_index = 0; contig_index < reference->count; ++contig_index) {
        uint64_t pairs = reference->contigs[contig_index].pair_count;
        uint64_t count = (pairs + PAIRS_PER_TASK - 1) / PAIRS_PER_TASK;
        if(count > SIZE_MAX - task_count) return -1;
        task_count += (size_t)count;
    }
    if(0 == task_count || task_count > SIZE_MAX / sizeof(*manifest->tasks)) {
        return -1;
    }
    manifest->tasks = calloc(task_count, sizeof(*manifest->tasks));
    if(NULL == manifest->tasks) return -1;
    manifest->count = task_count;

    for(contig_index = 0; contig_index < reference->count; ++contig_index) {
        uint64_t first_pair;
        uint64_t pairs = reference->contigs[contig_index].pair_count;
        for(first_pair = 0; first_pair < pairs;
            first_pair += PAIRS_PER_TASK) {
            pair_task_t *task = &manifest->tasks[task_index];
            uint64_t remaining = pairs - first_pair;
            task->id = task_index;
            task->contig_index = (uint32_t)contig_index;
            task->first_pair = first_pair;
            task->pair_count = remaining < PAIRS_PER_TASK ?
                               (uint32_t)remaining : PAIRS_PER_TASK;
            task_index++;
        }
    }
    return task_index == task_count ? 0 : -1;
}

static void
destroy_tasks(task_manifest_t *manifest)
{
    free(manifest->tasks);
    memset(manifest, 0, sizeof(*manifest));
}

static void
extract_read(const reference_contig_t *contig, uint64_t start,
             int reverse, uint8_t *read, int32_t length)
{
    int32_t index;

    if(!reverse) {
        memcpy(read, contig->bases + start, (size_t)length);
        return;
    }
    for(index = 0; index < length; ++index) {
        uint8_t base = contig->bases[start + (uint64_t)(length - index - 1)];
        read[index] = base < 4 ? (uint8_t)(3 - base) : 4;
    }
}

static int
count_n_bases(const uint8_t *read, int32_t length)
{
    int32_t index;
    int count = 0;
    for(index = 0; index < length; ++index) {
        if(4 == read[index]) count++;
    }
    return count;
}

static int
apply_errors(const dwgsim_opt_t *opt, uint32_t contig_index,
             uint64_t pair_index, int mate, int strand,
             uint8_t *read, int32_t length)
{
    local_rng_t rng;
    int errors = 0;
    int32_t step = strand ? -1 : 1;
    int32_t index = strand ? length - 1 : 0;

    rng_init(&rng, (uint32_t)opt->seed, DOMAIN_ERROR, contig_index,
             pair_index, (uint64_t)mate);
    while(0 <= index && index < length) {
        uint8_t base = read[index];
        double error_rate = opt->e[mate].start + opt->e[mate].by * index;
        if(base < 4 && rng_uniform(&rng) < error_rate) {
            base = (uint8_t)((base + (uint8_t)(rng_uniform(&rng) * 3.0 + 1.0)) & 3);
            read[index] = base;
            errors++;
        }
        index += step;
    }
    return errors;
}

static void
generate_qualities(const dwgsim_opt_t *opt, uint32_t contig_index,
                   uint64_t pair_index, int mate, char *quality,
                   int32_t length)
{
    local_rng_t rng;
    int32_t index;

    if(NULL != opt->fixed_quality) {
        memset(quality, opt->fixed_quality[0], (size_t)length);
        return;
    }

    rng_init(&rng, (uint32_t)opt->seed, DOMAIN_QUALITY, contig_index,
             pair_index, (uint64_t)mate);
    for(index = 0; index < length; ++index) {
        double error_rate = opt->e[mate].start + opt->e[mate].by * index;
        int score;
        if(0.0 < error_rate) {
            score = (int)(-10.0 * log(error_rate) / log(10.0) + 0.499);
        }
        else {
            score = QUAL_MAX;
        }
        if(0.0 < opt->quality_std) {
            score += (int)(rng_normal(&rng) * opt->quality_std + 0.5);
        }
        if(score < 0) score = 0;
        if(QUAL_MAX < score) score = QUAL_MAX;
        quality[index] = (char)('!' + score);
    }
}

static int
append_fastq_record(byte_buffer_t *output, const dwgsim_opt_t *opt,
                    const reference_contig_t *contig, uint64_t pair_index,
                    int mate, uint64_t coordinate1, uint64_t coordinate2,
                    int strand1, int strand2, int errors1, int errors2,
                    const uint8_t *read, const char *quality,
                    int32_t length)
{
    static const char bases[] = "ACGTN";
    int32_t index;

    if(0 != buffer_appendf(output,
            "@%s%s%s_%" PRIu64 "_%" PRIu64 "_%d_%d_0_0_%d:0:0_%d:0:0_%" PRIx64 "/%d\n",
            NULL == opt->read_prefix ? "" : opt->read_prefix,
            NULL == opt->read_prefix ? "" : "_",
            contig->name, coordinate1 + 1, coordinate2 + 1,
            strand1, strand2, errors1, errors2, pair_index, mate + 1)) {
        return -1;
    }
    if(0 != buffer_reserve(output, (size_t)length * 2 + 5)) return -1;
    for(index = 0; index < length; ++index) {
        output->data[output->length++] =
            (unsigned char)bases[read[index]];
    }
    memcpy(output->data + output->length, "\n+\n", 3);
    output->length += 3;
    memcpy(output->data + output->length, quality, (size_t)length);
    output->length += (size_t)length;
    output->data[output->length++] = '\n';
    return 0;
}

static int
generate_pair(const dwgsim_opt_t *opt, const reference_contig_t *contig,
              uint32_t contig_index, uint64_t pair_index,
              uint8_t *read1, uint8_t *read2,
              char *quality1, char *quality2,
              byte_buffer_t *output1, byte_buffer_t *output2)
{
    uint64_t start[2] = {0, 0};
    int strand[2];
    uint32_t attempt;
    int errors[2];
    local_rng_t strand_rng;

    rng_init(&strand_rng, (uint32_t)opt->seed, DOMAIN_STRAND,
             contig_index, pair_index, 0);
    switch(opt->read_one_strand) {
      case 1: strand[0] = 0; break;
      case 2: strand[0] = 1; break;
      default: strand[0] = rng_uniform(&strand_rng) < 0.5 ? 1 : 0; break;
    }
    strand[1] = 1 - strand[0];

    for(attempt = 0; attempt < MAX_FRAGMENT_ATTEMPTS; ++attempt) {
        local_rng_t fragment_rng;
        double sampled_distance;
        uint64_t distance;
        uint64_t position;
        uint64_t minimum_distance = (uint64_t)opt->length[0] +
                                    (uint64_t)opt->length[1];
        uint64_t range;

        rng_init(&fragment_rng, (uint32_t)opt->seed, DOMAIN_FRAGMENT,
                 contig_index, pair_index, attempt);
        sampled_distance = rng_normal(&fragment_rng) * opt->std_dev +
                           opt->dist;
        if(sampled_distance <= (double)minimum_distance) {
            distance = minimum_distance;
        }
        else if(sampled_distance >= (double)contig->length) {
            distance = contig->length;
        }
        else {
            distance = (uint64_t)(sampled_distance + 0.5);
        }
        if(contig->length < distance) distance = contig->length;
        range = contig->length - distance + 1;
        position = rng_bounded(&fragment_rng, range);

        if(0 == strand[0]) {
            start[0] = position;
            start[1] = position + distance - (uint64_t)opt->length[1];
        }
        else {
            start[0] = position + distance - (uint64_t)opt->length[0];
            start[1] = position;
        }
        extract_read(contig, start[0], strand[0], read1, opt->length[0]);
        extract_read(contig, start[1], strand[1], read2, opt->length[1]);
        if(count_n_bases(read1, opt->length[0]) <= opt->max_n &&
           count_n_bases(read2, opt->length[1]) <= opt->max_n) {
            break;
        }
    }
    if(MAX_FRAGMENT_ATTEMPTS == attempt) {
        fprintf(stderr,
                "[dwgsim_parallel] could not place pair %s:%" PRIu64
                " after %u attempts\n",
                contig->name, pair_index, MAX_FRAGMENT_ATTEMPTS);
        return -1;
    }

    errors[0] = apply_errors(opt, contig_index, pair_index, 0, strand[0],
                             read1, opt->length[0]);
    errors[1] = apply_errors(opt, contig_index, pair_index, 1, strand[1],
                             read2, opt->length[1]);
    generate_qualities(opt, contig_index, pair_index, 0, quality1,
                       opt->length[0]);
    generate_qualities(opt, contig_index, pair_index, 1, quality2,
                       opt->length[1]);

    if(0 != append_fastq_record(output1, opt, contig, pair_index, 0,
                                start[0], start[1], strand[0], strand[1],
                                errors[0], errors[1], read1, quality1,
                                opt->length[0]) ||
       0 != append_fastq_record(output2, opt, contig, pair_index, 1,
                                start[0], start[1], strand[0], strand[1],
                                errors[0], errors[1], read2, quality2,
                                opt->length[1])) {
        return -1;
    }
    return 0;
}

static void
pack_uint16(unsigned char *output, uint16_t value)
{
    output[0] = (unsigned char)value;
    output[1] = (unsigned char)(value >> 8);
}

static void
pack_uint32(unsigned char *output, uint32_t value)
{
    output[0] = (unsigned char)value;
    output[1] = (unsigned char)(value >> 8);
    output[2] = (unsigned char)(value >> 16);
    output[3] = (unsigned char)(value >> 24);
}

static int
compress_bgzf_block(const unsigned char *input, size_t input_length,
                    int compression_level, unsigned char *output,
                    size_t *output_length)
{
    z_stream stream;
    uint32_t crc;
    size_t compressed_length;
    int result;

    memset(&stream, 0, sizeof(stream));
    stream.next_in = (Bytef *)input;
    stream.avail_in = (uInt)input_length;
    stream.next_out = output + sizeof(bgzf_header);
    stream.avail_out = BGZF_MAX_BLOCK_SIZE - sizeof(bgzf_header) - 8;
    result = deflateInit2(&stream, compression_level, Z_DEFLATED, -15, 8,
                          Z_DEFAULT_STRATEGY);
    if(Z_OK != result) return -1;
    result = deflate(&stream, Z_FINISH);
    {
        int end_result = deflateEnd(&stream);
        if(Z_STREAM_END != result || Z_OK != end_result) return -1;
    }

    compressed_length = (size_t)stream.total_out + sizeof(bgzf_header) + 8;
    if(BGZF_MAX_BLOCK_SIZE < compressed_length) return -1;
    memcpy(output, bgzf_header, sizeof(bgzf_header));
    pack_uint16(output + 16, (uint16_t)(compressed_length - 1));
    crc = (uint32_t)crc32(crc32(0L, Z_NULL, 0), input,
                          (uInt)input_length);
    pack_uint32(output + compressed_length - 8, crc);
    pack_uint32(output + compressed_length - 4, (uint32_t)input_length);
    *output_length = compressed_length;
    return 0;
}

static int
compress_bgzf_buffer(const byte_buffer_t *input, int compression_level,
                     byte_buffer_t *output)
{
    unsigned char block[BGZF_MAX_BLOCK_SIZE];
    size_t offset = 0;

    while(offset < input->length) {
        size_t block_input = input->length - offset;
        size_t block_output = 0;
        if(BGZF_BLOCK_SIZE < block_input) block_input = BGZF_BLOCK_SIZE;
        if(0 != compress_bgzf_block(input->data + offset, block_input,
                                    compression_level, block,
                                    &block_output) ||
           0 != buffer_append(output, block, block_output)) {
            return -1;
        }
        offset += block_input;
    }
    return 0;
}

static void
destroy_task_result(paired_task_result_t *result)
{
    if(NULL == result) return;
    buffer_destroy(&result->compressed[0]);
    buffer_destroy(&result->compressed[1]);
    free(result);
}

static paired_task_result_t *
generate_task(const pipeline_t *pipeline, const pair_task_t *task)
{
    const reference_contig_t *contig =
        &pipeline->reference->contigs[task->contig_index];
    paired_task_result_t *result = NULL;
    byte_buffer_t raw[2] = {{0}, {0}};
    uint8_t *read1 = NULL;
    uint8_t *read2 = NULL;
    char *quality1 = NULL;
    char *quality2 = NULL;
    uint32_t index;

    result = calloc(1, sizeof(*result));
    read1 = malloc((size_t)pipeline->opt->length[0]);
    read2 = malloc((size_t)pipeline->opt->length[1]);
    quality1 = malloc((size_t)pipeline->opt->length[0]);
    quality2 = malloc((size_t)pipeline->opt->length[1]);
    if(NULL == result || NULL == read1 || NULL == read2 ||
       NULL == quality1 || NULL == quality2) goto error;

    result->task_id = task->id;
    result->pair_count = task->pair_count;
    for(index = 0; index < task->pair_count; ++index) {
        uint64_t pair_index = task->first_pair + index;
        if(0 != generate_pair(pipeline->opt, contig, task->contig_index,
                              pair_index, read1, read2, quality1, quality2,
                              &raw[0], &raw[1])) goto error;
    }

    if(0 != compress_bgzf_buffer(&raw[0], pipeline->opt->compression_level,
                                 &result->compressed[0]) ||
       0 != compress_bgzf_buffer(&raw[1], pipeline->opt->compression_level,
                                 &result->compressed[1])) goto error;

    buffer_destroy(&raw[0]);
    buffer_destroy(&raw[1]);
    free(read1);
    free(read2);
    free(quality1);
    free(quality2);
    return result;

error:
    buffer_destroy(&raw[0]);
    buffer_destroy(&raw[1]);
    free(read1);
    free(read2);
    free(quality1);
    free(quality2);
    destroy_task_result(result);
    return NULL;
}

static void
pipeline_fail_locked(pipeline_t *pipeline, const char *message)
{
    if(!pipeline->failed) {
        pipeline->failed = 1;
        snprintf(pipeline->error_message, sizeof(pipeline->error_message),
                 "%s", message);
    }
    pthread_cond_broadcast(&pipeline->changed);
}

static void
pipeline_fail(pipeline_t *pipeline, const char *message)
{
    pthread_mutex_lock(&pipeline->mutex);
    pipeline_fail_locked(pipeline, message);
    pthread_mutex_unlock(&pipeline->mutex);
}

static void *
generation_worker(void *data)
{
    pipeline_t *pipeline = data;

    for(;;) {
        uint64_t task_id;
        paired_task_result_t *result;
        result_slot_t *slot;

        pthread_mutex_lock(&pipeline->mutex);
        while(!pipeline->failed &&
              pipeline->next_claim < pipeline->tasks->count &&
              pipeline->next_claim >= pipeline->released + pipeline->window) {
            pthread_cond_wait(&pipeline->changed, &pipeline->mutex);
        }
        if(pipeline->failed ||
           pipeline->tasks->count <= pipeline->next_claim) {
            pthread_mutex_unlock(&pipeline->mutex);
            break;
        }
        task_id = pipeline->next_claim++;
        pthread_mutex_unlock(&pipeline->mutex);

        result = generate_task(pipeline, &pipeline->tasks->tasks[task_id]);
        if(NULL == result) {
            pipeline_fail(pipeline, "read generation or BGZF compression failed");
            break;
        }

        pthread_mutex_lock(&pipeline->mutex);
        if(pipeline->failed) {
            pthread_mutex_unlock(&pipeline->mutex);
            destroy_task_result(result);
            break;
        }
        slot = &pipeline->slots[task_id % pipeline->window];
        if(NULL != slot->result) {
            pipeline_fail_locked(pipeline, "internal result-ring collision");
            pthread_mutex_unlock(&pipeline->mutex);
            destroy_task_result(result);
            break;
        }
        slot->result = result;
        slot->consumed = 0;
        pthread_cond_broadcast(&pipeline->changed);
        pthread_mutex_unlock(&pipeline->mutex);
    }
    return NULL;
}

static int
write_bytes(FILE *output, const void *data, size_t length)
{
    return length == fwrite(data, 1, length, output) ? 0 : -1;
}

static void *
ordered_appender(void *data)
{
    appender_arg_t *arg = data;
    pipeline_t *pipeline = arg->pipeline;
    uint64_t expected;

    arg->status = -1;
    for(expected = 0; expected < pipeline->tasks->count; ++expected) {
        paired_task_result_t *result;
        result_slot_t *slot;

        pthread_mutex_lock(&pipeline->mutex);
        slot = &pipeline->slots[expected % pipeline->window];
        while(!pipeline->failed &&
              (NULL == slot->result || slot->result->task_id != expected)) {
            pthread_cond_wait(&pipeline->changed, &pipeline->mutex);
        }
        if(pipeline->failed) {
            pthread_mutex_unlock(&pipeline->mutex);
            goto close_output;
        }
        result = slot->result;
        pthread_mutex_unlock(&pipeline->mutex);

        if(0 != write_bytes(arg->output, result->compressed[arg->mate].data,
                            result->compressed[arg->mate].length)) {
            pipeline_fail(pipeline, "failed to append compressed FASTQ task");
            goto close_output;
        }

        pthread_mutex_lock(&pipeline->mutex);
        slot->consumed |= 1U << arg->mate;
        if(3U == slot->consumed) {
            if(pipeline->released != expected) {
                pipeline_fail_locked(pipeline,
                                     "internal paired-commit ordering failure");
                pthread_mutex_unlock(&pipeline->mutex);
                goto close_output;
            }
            destroy_task_result(slot->result);
            slot->result = NULL;
            slot->consumed = 0;
            pipeline->released++;
        }
        pthread_cond_broadcast(&pipeline->changed);
        pthread_mutex_unlock(&pipeline->mutex);
    }

    if(0 != write_bytes(arg->output, bgzf_eof, sizeof(bgzf_eof)) ||
       0 != fflush(arg->output) || 0 != fsync(fileno(arg->output))) {
        pipeline_fail(pipeline, "failed to finish BGZF FASTQ output");
        goto close_output;
    }
    arg->status = 0;

close_output:
    if(0 != fclose(arg->output)) {
        arg->status = -1;
        pipeline_fail(pipeline, "failed to close BGZF FASTQ output");
    }
    arg->output = NULL;
    return NULL;
}

static int
run_pipeline(const dwgsim_opt_t *opt,
             const reference_manifest_t *reference,
             const task_manifest_t *tasks,
             FILE *output1, FILE *output2)
{
    pipeline_t pipeline;
    appender_arg_t appender[2];
    pthread_t appender_threads[2];
    pthread_t *workers = NULL;
    size_t worker_count = (size_t)opt->compression_threads;
    size_t workers_created = 0;
    int appenders_created = 0;
    size_t index;
    int status = -1;

    memset(&pipeline, 0, sizeof(pipeline));
    pipeline.opt = opt;
    pipeline.reference = reference;
    pipeline.tasks = tasks;
    pipeline.window = worker_count * 2;
    if(pipeline.window < 8) pipeline.window = 8;
    if(tasks->count < pipeline.window) pipeline.window = tasks->count;
    pipeline.slots = calloc(pipeline.window, sizeof(*pipeline.slots));
    if(NULL == pipeline.slots ||
       0 != pthread_mutex_init(&pipeline.mutex, NULL) ||
       0 != pthread_cond_init(&pipeline.changed, NULL)) {
        free(pipeline.slots);
        fclose(output1);
        fclose(output2);
        return -1;
    }

    if(tasks->count < worker_count) worker_count = tasks->count;
    if(0 == worker_count) worker_count = 1;
    workers = calloc(worker_count, sizeof(*workers));
    if(NULL == workers) {
        pipeline_fail(&pipeline, "could not allocate generation threads");
    }

    memset(appender, 0, sizeof(appender));
    appender[0].pipeline = &pipeline;
    appender[0].mate = 0;
    appender[0].output = output1;
    appender[1].pipeline = &pipeline;
    appender[1].mate = 1;
    appender[1].output = output2;

    if(NULL != workers) {
        for(index = 0; index < 2; ++index) {
            if(0 != pthread_create(&appender_threads[index], NULL,
                                   ordered_appender, &appender[index])) {
                pipeline_fail(&pipeline,
                              "could not create FASTQ appender thread");
                break;
            }
            appenders_created++;
        }
    }
    while(NULL != workers && appenders_created == 2 &&
          workers_created < worker_count) {
        if(0 != pthread_create(&workers[workers_created], NULL,
                               generation_worker, &pipeline)) {
            pipeline_fail(&pipeline, "could not create generation worker");
            break;
        }
        workers_created++;
    }

    for(index = 0; index < workers_created; ++index) {
        pthread_join(workers[index], NULL);
    }
    for(index = 0; index < (size_t)appenders_created; ++index) {
        pthread_join(appender_threads[index], NULL);
    }
    for(index = (size_t)appenders_created; index < 2; ++index) {
        if(NULL != appender[index].output) fclose(appender[index].output);
        appender[index].output = NULL;
    }

    if(!pipeline.failed && 0 == appender[0].status &&
       0 == appender[1].status && pipeline.released == tasks->count) {
        status = 0;
    }
    else {
        fprintf(stderr, "[dwgsim_parallel] pipeline failed: %s\n",
                pipeline.error_message[0] ? pipeline.error_message :
                "incomplete paired output");
    }

    for(index = 0; index < pipeline.window; ++index) {
        destroy_task_result(pipeline.slots[index].result);
    }
    free(workers);
    free(pipeline.slots);
    pthread_cond_destroy(&pipeline.changed);
    pthread_mutex_destroy(&pipeline.mutex);
    return status;
}

static int
publish_outputs(const char *staging1, const char *staging2,
                const char *output1, const char *output2,
                const char *manifest_path, const char *manifest_staging,
                const dwgsim_opt_t *opt, size_t task_count)
{
    struct stat stat1;
    struct stat stat2;
    FILE *manifest;

    if(0 != rename(staging1, output1)) {
        fprintf(stderr, "[dwgsim_parallel] cannot publish %s: %s\n",
                output1, strerror(errno));
        return -1;
    }
    if(0 != rename(staging2, output2)) {
        fprintf(stderr, "[dwgsim_parallel] cannot publish %s: %s\n",
                output2, strerror(errno));
        return -1;
    }
    if(0 != stat(output1, &stat1) || 0 != stat(output2, &stat2)) return -1;

    manifest = fopen(manifest_staging, "w");
    if(NULL == manifest) return -1;
    if(0 > fprintf(manifest,
                   "format=dwgsim-deterministic-wgs-v%d\n"
                   "seed=%d\nread_pairs=%" PRId64 "\nthreads=%d\n"
                   "compression_level=%d\ntasks=%zu\n"
                   "read1_bytes=%" PRIu64 "\nread2_bytes=%" PRIu64 "\n",
                   PARALLEL_FORMAT_VERSION, opt->seed, opt->N,
                   opt->compression_threads, opt->compression_level,
                   task_count, (uint64_t)stat1.st_size,
                   (uint64_t)stat2.st_size) ||
       0 != fflush(manifest) || 0 != fsync(fileno(manifest))) {
        fclose(manifest);
        unlink(manifest_staging);
        return -1;
    }
    if(0 != fclose(manifest)) {
        unlink(manifest_staging);
        return -1;
    }
    if(0 != rename(manifest_staging, manifest_path)) {
        unlink(manifest_staging);
        return -1;
    }
    return 0;
}

int
dwgsim_parallel_wgs_supported(const dwgsim_opt_t *opt,
                              const char *reference_path)
{
    char *fai_path;
    int supported;

    if(NULL == opt || NULL == reference_path ||
       0 == strcmp(reference_path, "-") ||
       OUTPUT_TYPE_READS != opt->output_type ||
       READS_OUTPUT_TYPE_BWA != opt->reads_output_type ||
       0 != opt->data_type || 0 >= opt->length[1] || 0 >= opt->N ||
       0.0 != opt->mut_rate || NULL != opt->fn_muts_input ||
       0.0 != opt->rand_read || NULL != opt->fn_regions_bed ||
       opt->amplicons || opt->is_inner ||
       (0 != opt->strandedness && 2 != opt->strandedness)) {
        return 0;
    }
    fai_path = path_with_suffix(reference_path, ".fai");
    if(NULL == fai_path) return 0;
    supported = 0 == access(reference_path, R_OK) &&
                0 == access(fai_path, R_OK);
    free(fai_path);
    return supported;
}

int
dwgsim_parallel_wgs_run(const dwgsim_opt_t *opt,
                        const char *reference_path,
                        const char *output_prefix)
{
    reference_manifest_t reference = {0};
    task_manifest_t tasks = {0};
    char *output1 = NULL;
    char *output2 = NULL;
    char *manifest_path = NULL;
    char *staging1 = NULL;
    char *staging2 = NULL;
    char *manifest_staging = NULL;
    char staging_suffix[64];
    FILE *file1 = NULL;
    FILE *file2 = NULL;
    struct timespec start;
    struct timespec finish;
    double elapsed;
    int status = 1;

    clock_gettime(CLOCK_MONOTONIC, &start);
    output1 = path_with_suffix(output_prefix, ".bwa.read1.fastq.gz");
    output2 = path_with_suffix(output_prefix, ".bwa.read2.fastq.gz");
    manifest_path = path_with_suffix(output_prefix, ".dwgsim.complete");
    snprintf(staging_suffix, sizeof(staging_suffix), ".partial.%ld",
             (long)getpid());
    if(NULL != output1) staging1 = path_with_suffix(output1, staging_suffix);
    if(NULL != output2) staging2 = path_with_suffix(output2, staging_suffix);
    if(NULL != manifest_path) {
        manifest_staging = path_with_suffix(manifest_path, staging_suffix);
    }
    if(NULL == output1 || NULL == output2 || NULL == manifest_path ||
       NULL == staging1 || NULL == staging2 || NULL == manifest_staging) {
        fprintf(stderr, "[dwgsim_parallel] output path allocation failed\n");
        goto cleanup;
    }
    unlink(manifest_path);

    fprintf(stderr,
            "[dwgsim_parallel] deterministic WGS v%d, %d workers, BGZF level %d\n",
            PARALLEL_FORMAT_VERSION, opt->compression_threads,
            opt->compression_level);
    if(0 != read_reference_manifest(reference_path, &reference)) goto cleanup;
    fprintf(stderr, "[dwgsim_parallel] %zu sequences, total length: %" PRIu64
                    "\n", reference.count, reference.total_length);
    if(0 != prepare_reference(opt, reference_path, &reference) ||
       0 != allocate_pairs(&reference, (uint64_t)opt->N) ||
       0 != build_tasks(&reference, &tasks)) goto cleanup;

    fprintf(stderr,
            "[dwgsim_parallel] planned %" PRId64
            " pairs in %zu fixed tasks\n", opt->N, tasks.count);
    file1 = fopen(staging1, "wb");
    file2 = fopen(staging2, "wb");
    if(NULL == file1 || NULL == file2) {
        fprintf(stderr, "[dwgsim_parallel] cannot create staging FASTQs: %s\n",
                strerror(errno));
        if(NULL != file1) fclose(file1);
        if(NULL != file2) fclose(file2);
        file1 = file2 = NULL;
        goto cleanup;
    }
    if(0 != run_pipeline(opt, &reference, &tasks, file1, file2)) {
        file1 = file2 = NULL;
        goto cleanup;
    }
    file1 = file2 = NULL;

    if(0 != publish_outputs(staging1, staging2, output1, output2,
                            manifest_path, manifest_staging, opt,
                            tasks.count)) {
        fprintf(stderr, "[dwgsim_parallel] failed to publish paired outputs\n");
        goto cleanup;
    }

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (double)(finish.tv_sec - start.tv_sec) +
              (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    fprintf(stderr,
            "[dwgsim_parallel] complete: %" PRId64
            " pairs in %.3f s (%.2f pairs/s)\n",
            opt->N, elapsed,
            0.0 < elapsed ? (double)opt->N / elapsed : 0.0);
    status = 0;

cleanup:
    if(NULL != file1) fclose(file1);
    if(NULL != file2) fclose(file2);
    if(0 != status) {
        if(NULL != staging1) unlink(staging1);
        if(NULL != staging2) unlink(staging2);
        if(NULL != manifest_staging) unlink(manifest_staging);
    }
    destroy_tasks(&tasks);
    destroy_reference_manifest(&reference);
    free(output1);
    free(output2);
    free(manifest_path);
    free(staging1);
    free(staging2);
    free(manifest_staging);
    return status;
}
