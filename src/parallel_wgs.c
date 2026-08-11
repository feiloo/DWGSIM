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
#define MATCHED_FORMAT_VERSION 1
#define PAIRS_PER_TASK 8192U
#define MAX_FRAGMENT_ATTEMPTS 10000U
#define MAX_VARIANT_LENGTH 1048576U
#define REFERENCE_READ_BUFFER (1024U * 1024U)
#define QUAL_MAX 40
#define MAX_SAMPLES 2
#define MAX_STREAMS (MAX_SAMPLES * 2)

#define SAMPLE_NORMAL 0
#define SAMPLE_TUMOR 1

#define VARIANT_SNV 0
#define VARIANT_INSERTION 1
#define VARIANT_DELETION 2

#define VARIANT_GERMLINE 1
#define VARIANT_SOMATIC 2

#define DOMAIN_FRAGMENT UINT64_C(0x465241474d454e54)
#define DOMAIN_STRAND   UINT64_C(0x535452414e442020)
#define DOMAIN_ERROR    UINT64_C(0x4552524f52202020)
#define DOMAIN_QUALITY  UINT64_C(0x5155414c49545920)
#define DOMAIN_HAPLOTYPE UINT64_C(0x4841504c4f545950)
#define DOMAIN_GERMLINE UINT64_C(0x4745524d4c494e45)
#define DOMAIN_SOMATIC  UINT64_C(0x534f4d4154494320)

typedef struct {
    uint32_t position;
    uint32_t length;
    uint8_t type;
    uint8_t haplotype_mask;
    uint8_t scope;
    uint8_t alternate_base;
    uint8_t *inserted_bases;
} variant_event_t;

typedef struct {
    variant_event_t *events;
    size_t count;
    size_t capacity;
} variant_set_t;

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
    uint64_t pair_count[MAX_SAMPLES];
    variant_set_t variants;
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
    uint64_t first_pair[MAX_SAMPLES];
    uint32_t pair_count[MAX_SAMPLES];
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
    uint32_t pair_count[MAX_SAMPLES];
    byte_buffer_t compressed[MAX_STREAMS];
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
    unsigned sample_count;
    unsigned stream_count;
    uint64_t next_claim;
    uint64_t released;
    int failed;
    char error_message[256];
    pthread_mutex_t mutex;
    pthread_cond_t changed;
} pipeline_t;

typedef struct {
    pipeline_t *pipeline;
    unsigned stream;
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

static uint64_t
sample_random_extra(unsigned sample, uint64_t extra)
{
    if(SAMPLE_NORMAL == sample) return extra;
    return mix64(extra ^ UINT64_C(0x54554d4f5253414d));
}

static void
variant_set_destroy(variant_set_t *set)
{
    size_t index;

    if(NULL == set) return;
    for(index = 0; index < set->count; ++index) {
        free(set->events[index].inserted_bases);
    }
    free(set->events);
    memset(set, 0, sizeof(*set));
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
        variant_set_destroy(&manifest->contigs[index].variants);
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

static uint64_t
variant_end(const variant_event_t *event)
{
    return (uint64_t)event->position +
           (VARIANT_DELETION == event->type ? event->length : 1U);
}

static int
variant_set_append(variant_set_t *set, const variant_event_t *event)
{
    variant_event_t *events;
    size_t capacity;

    if(set->count == set->capacity) {
        capacity = set->capacity ? set->capacity * 2 : 1024;
        if(capacity < set->capacity ||
           capacity > SIZE_MAX / sizeof(*events)) return -1;
        events = realloc(set->events, capacity * sizeof(*events));
        if(NULL == events) return -1;
        set->events = events;
        set->capacity = capacity;
    }
    set->events[set->count++] = *event;
    return 0;
}

static size_t
variant_lower_bound(const variant_set_t *set, uint32_t position)
{
    size_t low = 0;
    size_t high = set->count;

    while(low < high) {
        size_t middle = low + (high - low) / 2;
        if(set->events[middle].position < position) low = middle + 1;
        else high = middle;
    }
    return low;
}

static int
variant_conflicts(const variant_set_t *set, const variant_event_t *candidate)
{
    size_t index = variant_lower_bound(set, candidate->position);
    uint64_t end = variant_end(candidate);

    if(0 < index &&
       variant_end(&set->events[index - 1]) > candidate->position) return 1;
    if(index < set->count && set->events[index].position < end) return 1;
    return 0;
}

static int
reference_interval_is_acgt(const reference_contig_t *contig,
                           uint64_t position, uint64_t length)
{
    uint64_t index;

    if(position > contig->length || length > contig->length - position) {
        return 0;
    }
    for(index = 0; index < length; ++index) {
        if(4 <= contig->bases[position + index]) return 0;
    }
    return 1;
}

static uint32_t
draw_indel_length(const dwgsim_opt_t *opt, local_rng_t *rng)
{
    uint64_t length = (uint32_t)opt->indel_min;

    if(0.0 < opt->indel_extend && length < MAX_VARIANT_LENGTH) {
        uint64_t extra;
        if(1.0 <= opt->indel_extend) {
            extra = MAX_VARIANT_LENGTH - length;
        }
        else {
            double uniform = rng_uniform(rng);
            double sampled;
            if(uniform <= 0.0) uniform = 0x1.0p-53;
            sampled = floor(log(uniform) / log(opt->indel_extend));
            extra = sampled >= (double)(MAX_VARIANT_LENGTH - length) ?
                    MAX_VARIANT_LENGTH - length : (uint64_t)sampled;
        }
        length += extra;
    }
    return (uint32_t)length;
}

static int
generate_variant_scope(const dwgsim_opt_t *opt,
                       reference_contig_t *contig,
                       uint32_t contig_index, uint8_t scope,
                       const variant_set_t *avoid, variant_set_t *output)
{
    double rate = VARIANT_GERMLINE == scope ?
                  opt->mut_rate : opt->somatic_rate;
    uint64_t cursor = 0;
    local_rng_t rng;

    if(rate <= 0.0 || !contig->eligible) return 0;
    rng_init(&rng, (uint32_t)opt->seed,
             VARIANT_GERMLINE == scope ? DOMAIN_GERMLINE : DOMAIN_SOMATIC,
             contig_index, 0, 0);

    while(cursor < contig->length) {
        variant_event_t event;
        uint64_t position;
        uint64_t remaining = contig->length - cursor;

        if(rate < 1.0) {
            double uniform = rng_uniform(&rng);
            double gap = floor(log1p(-uniform) / log1p(-rate));
            if(!isfinite(gap) || gap >= (double)remaining) break;
            position = cursor + (uint64_t)gap;
        }
        else {
            position = cursor;
        }
        cursor = position + 1;

        memset(&event, 0, sizeof(event));
        event.position = (uint32_t)position;
        event.length = 1;
        event.scope = scope;
        if(rng_uniform(&rng) < opt->indel_frac) {
            event.type = rng_uniform(&rng) < 0.5 ?
                         VARIANT_DELETION : VARIANT_INSERTION;
            event.length = draw_indel_length(opt, &rng);
        }
        else {
            event.type = VARIANT_SNV;
        }

        if(VARIANT_DELETION == event.type) {
            if(0 == position || event.length > contig->length - position ||
               !reference_interval_is_acgt(contig, position - 1,
                                            (uint64_t)event.length + 1)) {
                continue;
            }
        }
        else if(!reference_interval_is_acgt(contig, position, 1)) {
            continue;
        }
        if((NULL != avoid && variant_conflicts(avoid, &event)) ||
           variant_conflicts(output, &event)) {
            continue;
        }

        if(VARIANT_GERMLINE == scope && rng_uniform(&rng) < (1.0 / 3.0)) {
            event.haplotype_mask = 3;
        }
        else {
            event.haplotype_mask = rng_uniform(&rng) < 0.5 ? 1 : 2;
        }
        if(VARIANT_SNV == event.type) {
            uint8_t reference_base = contig->bases[position];
            event.alternate_base = (uint8_t)((reference_base + 1U +
                rng_bounded(&rng, 3)) & 3U);
        }
        else if(VARIANT_INSERTION == event.type) {
            uint32_t index;
            event.inserted_bases = malloc(event.length);
            if(NULL == event.inserted_bases) return -1;
            for(index = 0; index < event.length; ++index) {
                event.inserted_bases[index] = (uint8_t)rng_bounded(&rng, 4);
            }
        }
        if(0 != variant_set_append(output, &event)) {
            free(event.inserted_bases);
            return -1;
        }
        if(VARIANT_DELETION == event.type) cursor = variant_end(&event);
    }
    return 0;
}

static int
merge_variant_sets(reference_contig_t *contig,
                   variant_set_t *germline, variant_set_t *somatic)
{
    variant_set_t merged = {0};
    size_t germline_index = 0;
    size_t somatic_index = 0;
    size_t output_index = 0;

    if(germline->count > SIZE_MAX - somatic->count) return -1;
    merged.count = germline->count + somatic->count;
    merged.capacity = merged.count;
    if(0 < merged.count) {
        merged.events = calloc(merged.count, sizeof(*merged.events));
        if(NULL == merged.events) return -1;
    }
    while(germline_index < germline->count ||
          somatic_index < somatic->count) {
        if(somatic_index == somatic->count ||
           (germline_index < germline->count &&
            germline->events[germline_index].position <
            somatic->events[somatic_index].position)) {
            merged.events[output_index++] = germline->events[germline_index++];
        }
        else {
            merged.events[output_index++] = somatic->events[somatic_index++];
        }
    }
    free(germline->events);
    free(somatic->events);
    memset(germline, 0, sizeof(*germline));
    memset(somatic, 0, sizeof(*somatic));
    contig->variants = merged;
    return 0;
}

static int
generate_contig_variants(const dwgsim_opt_t *opt,
                         reference_contig_t *contig,
                         uint32_t contig_index)
{
    variant_set_t germline = {0};
    variant_set_t somatic = {0};
    int status = -1;

    if(0 != generate_variant_scope(opt, contig, contig_index,
                                   VARIANT_GERMLINE, NULL, &germline) ||
       0 != generate_variant_scope(opt, contig, contig_index,
                                   VARIANT_SOMATIC, &germline, &somatic) ||
       0 != merge_variant_sets(contig, &germline, &somatic)) {
        goto cleanup;
    }
    status = 0;

cleanup:
    variant_set_destroy(&germline);
    variant_set_destroy(&somatic);
    return status;
}

typedef struct {
    const dwgsim_opt_t *opt;
    reference_manifest_t *manifest;
    volatile size_t next_contig;
    volatile int failed;
} variant_loader_t;

static void *
variant_loader_worker(void *data)
{
    variant_loader_t *loader = data;

    for(;;) {
        size_t index = __sync_fetch_and_add(&loader->next_contig, 1);
        if(loader->manifest->count <= index) break;
        if(__sync_fetch_and_add(&loader->failed, 0)) continue;
        if(0 != generate_contig_variants(loader->opt,
                                         &loader->manifest->contigs[index],
                                         (uint32_t)index)) {
            fprintf(stderr, "[dwgsim_parallel] failed to prepare variants for %s\n",
                    loader->manifest->contigs[index].name);
            __sync_lock_test_and_set(&loader->failed, 1);
        }
    }
    return NULL;
}

static int
prepare_variants(const dwgsim_opt_t *opt, reference_manifest_t *manifest)
{
    variant_loader_t loader;
    pthread_t *threads = NULL;
    size_t worker_count = (size_t)opt->compression_threads;
    size_t created = 0;
    size_t index;

    memset(&loader, 0, sizeof(loader));
    loader.opt = opt;
    loader.manifest = manifest;
    if(manifest->count < worker_count) worker_count = manifest->count;
    if(0 == worker_count) worker_count = 1;
    threads = calloc(worker_count, sizeof(*threads));
    if(NULL == threads) return -1;
    while(created < worker_count) {
        if(0 != pthread_create(&threads[created], NULL,
                               variant_loader_worker, &loader)) {
            __sync_lock_test_and_set(&loader.failed, 1);
            break;
        }
        created++;
    }
    for(index = 0; index < created; ++index) {
        pthread_join(threads[index], NULL);
    }
    free(threads);
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
allocate_pairs(reference_manifest_t *manifest, uint64_t requested_pairs,
               unsigned sample)
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
        contig->pair_count[sample] = integral;
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
        manifest->contigs[remainders[index].contig_index].pair_count[sample]++;
    }
    free(remainders);

    return 0;
}

static void
discard_unused_contigs(reference_manifest_t *manifest, unsigned sample_count)
{
    size_t contig_index;

    for(contig_index = 0; contig_index < manifest->count; ++contig_index) {
        reference_contig_t *contig = &manifest->contigs[contig_index];
        unsigned sample;
        int used = 0;
        for(sample = 0; sample < sample_count; ++sample) {
            if(0 < contig->pair_count[sample]) used = 1;
        }
        if(!used) {
            free(contig->bases);
            contig->bases = NULL;
        }
    }
}

static int
build_tasks(const reference_manifest_t *reference,
            task_manifest_t *manifest, unsigned sample_count)
{
    size_t task_count = 0;
    size_t task_index = 0;
    size_t contig_index;

    for(contig_index = 0; contig_index < reference->count; ++contig_index) {
        uint64_t count = 0;
        unsigned sample;
        for(sample = 0; sample < sample_count; ++sample) {
            uint64_t pairs = reference->contigs[contig_index].pair_count[sample];
            uint64_t sample_tasks =
                (pairs + PAIRS_PER_TASK - 1) / PAIRS_PER_TASK;
            if(count < sample_tasks) count = sample_tasks;
        }
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
        uint64_t task_ordinal;
        uint64_t contig_tasks = 0;
        unsigned sample;
        for(sample = 0; sample < sample_count; ++sample) {
            uint64_t pairs = reference->contigs[contig_index].pair_count[sample];
            uint64_t sample_tasks =
                (pairs + PAIRS_PER_TASK - 1) / PAIRS_PER_TASK;
            if(contig_tasks < sample_tasks) contig_tasks = sample_tasks;
        }
        for(task_ordinal = 0; task_ordinal < contig_tasks; ++task_ordinal) {
            pair_task_t *task = &manifest->tasks[task_index];
            task->id = task_index;
            task->contig_index = (uint32_t)contig_index;
            for(sample = 0; sample < sample_count; ++sample) {
                uint64_t pairs = reference->contigs[contig_index].pair_count[sample];
                uint64_t first_pair = task_ordinal * PAIRS_PER_TASK;
                uint64_t remaining = first_pair < pairs ? pairs - first_pair : 0;
                task->first_pair[sample] = first_pair;
                task->pair_count[sample] = remaining < PAIRS_PER_TASK ?
                    (uint32_t)remaining : PAIRS_PER_TASK;
            }
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
variant_is_active(const variant_event_t *event, unsigned haplotype,
                  int somatic_clone)
{
    if(0 == (event->haplotype_mask & (1U << haplotype))) return 0;
    return VARIANT_GERMLINE == event->scope || somatic_clone;
}

static int
extract_variant_read_forward(const reference_contig_t *contig,
                             uint64_t start, unsigned haplotype,
                             int somatic_clone, uint8_t *read,
                             int32_t length, int *substitutions,
                             int *indels)
{
    const variant_set_t *set = &contig->variants;
    size_t event_index = variant_lower_bound(set, (uint32_t)start);
    uint64_t position = start;
    int32_t output_index = 0;

    if(0 < event_index &&
       VARIANT_DELETION == set->events[event_index - 1].type &&
       variant_end(&set->events[event_index - 1]) > position) {
        event_index--;
    }
    while(output_index < length && position < contig->length) {
        const variant_event_t *event = NULL;
        int active;

        while(event_index < set->count &&
              variant_end(&set->events[event_index]) <= position) {
            event_index++;
        }
        if(event_index < set->count) event = &set->events[event_index];
        if(NULL == event || position < event->position) {
            read[output_index++] = contig->bases[position++];
            continue;
        }
        active = variant_is_active(event, haplotype, somatic_clone);
        if(VARIANT_DELETION == event->type &&
           position < variant_end(event)) {
            if(active) {
                position = variant_end(event);
                (*indels)++;
                event_index++;
            }
            else {
                read[output_index++] = contig->bases[position++];
            }
            continue;
        }
        if(position != event->position) {
            event_index++;
            continue;
        }
        if(VARIANT_SNV == event->type) {
            read[output_index++] = active ? event->alternate_base :
                                           contig->bases[position];
            if(active) (*substitutions)++;
        }
        else {
            uint32_t insertion_index;
            read[output_index++] = contig->bases[position];
            if(active) {
                (*indels)++;
                for(insertion_index = 0;
                    insertion_index < event->length &&
                    output_index < length; ++insertion_index) {
                    read[output_index++] = event->inserted_bases[insertion_index];
                }
            }
        }
        position++;
        event_index++;
    }
    return output_index == length ? 0 : -1;
}

static int
extract_variant_read_reverse(const reference_contig_t *contig,
                             uint64_t start, unsigned haplotype,
                             int somatic_clone, uint8_t *read,
                             int32_t length, int *substitutions,
                             int *indels)
{
    const variant_set_t *set = &contig->variants;
    int64_t position = (int64_t)(start + (uint64_t)length - 1);
    size_t event_cursor = variant_lower_bound(set, (uint32_t)(position + 1));
    int32_t output_index = 0;

    while(output_index < length && 0 <= position) {
        const variant_event_t *event = NULL;
        int active;

        while(0 < event_cursor &&
              (uint64_t)set->events[event_cursor - 1].position >
              (uint64_t)position) {
            event_cursor--;
        }
        if(0 < event_cursor) event = &set->events[event_cursor - 1];
        if(NULL == event ||
           (VARIANT_DELETION != event->type &&
            event->position < (uint64_t)position) ||
           (VARIANT_DELETION == event->type &&
            variant_end(event) <= (uint64_t)position)) {
            uint8_t base = contig->bases[position--];
            read[output_index++] = base < 4 ? (uint8_t)(3 - base) : 4;
            continue;
        }
        active = variant_is_active(event, haplotype, somatic_clone);
        if(VARIANT_DELETION == event->type &&
           event->position <= (uint64_t)position &&
           (uint64_t)position < variant_end(event)) {
            if(active) {
                position = (int64_t)event->position - 1;
                (*indels)++;
                event_cursor--;
            }
            else {
                uint8_t base = contig->bases[position--];
                read[output_index++] = (uint8_t)(3 - base);
                if(position < (int64_t)event->position) event_cursor--;
            }
            continue;
        }
        if(event->position != (uint64_t)position) {
            event_cursor--;
            continue;
        }
        if(VARIANT_SNV == event->type) {
            uint8_t base = active ? event->alternate_base :
                                    contig->bases[position];
            read[output_index++] = (uint8_t)(3 - base);
            if(active) (*substitutions)++;
        }
        else {
            if(active) {
                uint32_t insertion_index = event->length;
                (*indels)++;
                while(0 < insertion_index && output_index < length) {
                    uint8_t base = event->inserted_bases[--insertion_index];
                    read[output_index++] = (uint8_t)(3 - base);
                }
            }
            if(output_index < length) {
                uint8_t base = contig->bases[position];
                read[output_index++] = (uint8_t)(3 - base);
            }
        }
        position--;
        event_cursor--;
    }
    return output_index == length ? 0 : -1;
}

static int
extract_biological_read(const reference_contig_t *contig, uint64_t start,
                        int reverse, unsigned haplotype, int somatic_clone,
                        uint8_t *read, int32_t length, int *substitutions,
                        int *indels)
{
    *substitutions = 0;
    *indels = 0;
    if(0 == contig->variants.count) {
        extract_read(contig, start, reverse, read, length);
        return 0;
    }
    if(reverse) {
        return extract_variant_read_reverse(contig, start, haplotype,
                                            somatic_clone, read, length,
                                            substitutions, indels);
    }
    return extract_variant_read_forward(contig, start, haplotype,
                                        somatic_clone, read, length,
                                        substitutions, indels);
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
             uint64_t pair_index, unsigned sample, int mate, int strand,
             uint8_t *read, int32_t length)
{
    local_rng_t rng;
    int errors = 0;
    int32_t step = strand ? -1 : 1;
    int32_t index = strand ? length - 1 : 0;

    rng_init(&rng, (uint32_t)opt->seed, DOMAIN_ERROR, contig_index,
             pair_index, sample_random_extra(sample, (uint64_t)mate));
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
                   uint64_t pair_index, unsigned sample, int mate, char *quality,
                   int32_t length)
{
    local_rng_t rng;
    int32_t index;

    if(NULL != opt->fixed_quality) {
        memset(quality, opt->fixed_quality[0], (size_t)length);
        return;
    }

    rng_init(&rng, (uint32_t)opt->seed, DOMAIN_QUALITY, contig_index,
             pair_index, sample_random_extra(sample, (uint64_t)mate));
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
                    unsigned sample, int mate,
                    uint64_t coordinate1, uint64_t coordinate2,
                    int strand1, int strand2, int errors1, int errors2,
                    int substitutions1, int substitutions2,
                    int indels1, int indels2,
                    const uint8_t *read, const char *quality,
                    int32_t length)
{
    static const char bases[] = "ACGTN";
    int32_t index;

    if(0 != buffer_appendf(output,
            "@%s%s%s%s_%" PRIu64 "_%" PRIu64 "_%d_%d_0_0_%d:%d:%d_%d:%d:%d_%" PRIx64 "/%d\n",
            NULL == opt->read_prefix ? "" : opt->read_prefix,
            NULL == opt->read_prefix ? "" : "_",
            opt->matched ? (SAMPLE_NORMAL == sample ? "normal_" : "tumor_") : "",
            contig->name, coordinate1 + 1, coordinate2 + 1,
            strand1, strand2,
            errors1, substitutions1, indels1,
            errors2, substitutions2, indels2,
            pair_index, mate + 1)) {
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
              uint32_t contig_index, uint64_t pair_index, unsigned sample,
              uint8_t *read1, uint8_t *read2,
              char *quality1, char *quality2,
              byte_buffer_t *output1, byte_buffer_t *output2)
{
    uint64_t start[2] = {0, 0};
    int strand[2];
    uint32_t attempt;
    int errors[2];
    int substitutions[2] = {0, 0};
    int indels[2] = {0, 0};
    unsigned haplotype;
    int somatic_clone;
    local_rng_t strand_rng;
    local_rng_t haplotype_rng;

    rng_init(&haplotype_rng, (uint32_t)opt->seed, DOMAIN_HAPLOTYPE,
             contig_index, pair_index, sample_random_extra(sample, 0));
    haplotype = rng_uniform(&haplotype_rng) < 0.5 ? 0U : 1U;
    somatic_clone = SAMPLE_TUMOR == sample &&
                    rng_uniform(&haplotype_rng) < 2.0 * opt->tumor_vaf;

    rng_init(&strand_rng, (uint32_t)opt->seed, DOMAIN_STRAND,
             contig_index, pair_index, sample_random_extra(sample, 0));
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
                 contig_index, pair_index,
                 sample_random_extra(sample, attempt));
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
        if(0 == extract_biological_read(contig, start[0], strand[0],
                                        haplotype, somatic_clone, read1,
                                        opt->length[0], &substitutions[0],
                                        &indels[0]) &&
           0 == extract_biological_read(contig, start[1], strand[1],
                                        haplotype, somatic_clone, read2,
                                        opt->length[1], &substitutions[1],
                                        &indels[1]) &&
           count_n_bases(read1, opt->length[0]) <= opt->max_n &&
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

    errors[0] = apply_errors(opt, contig_index, pair_index, sample, 0, strand[0],
                             read1, opt->length[0]);
    errors[1] = apply_errors(opt, contig_index, pair_index, sample, 1, strand[1],
                             read2, opt->length[1]);
    generate_qualities(opt, contig_index, pair_index, sample, 0, quality1,
                       opt->length[0]);
    generate_qualities(opt, contig_index, pair_index, sample, 1, quality2,
                       opt->length[1]);

    if(0 != append_fastq_record(output1, opt, contig, pair_index, sample, 0,
                                start[0], start[1], strand[0], strand[1],
                                errors[0], errors[1], substitutions[0],
                                substitutions[1], indels[0], indels[1],
                                read1, quality1,
                                opt->length[0]) ||
       0 != append_fastq_record(output2, opt, contig, pair_index, sample, 1,
                                start[0], start[1], strand[0], strand[1],
                                errors[0], errors[1], substitutions[0],
                                substitutions[1], indels[0], indels[1],
                                read2, quality2,
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
    unsigned stream;

    if(NULL == result) return;
    for(stream = 0; stream < MAX_STREAMS; ++stream) {
        buffer_destroy(&result->compressed[stream]);
    }
    free(result);
}

static paired_task_result_t *
generate_task(const pipeline_t *pipeline, const pair_task_t *task)
{
    const reference_contig_t *contig =
        &pipeline->reference->contigs[task->contig_index];
    paired_task_result_t *result = NULL;
    byte_buffer_t raw[MAX_STREAMS] = {{0}};
    uint8_t *read1 = NULL;
    uint8_t *read2 = NULL;
    char *quality1 = NULL;
    char *quality2 = NULL;
    unsigned sample;
    unsigned stream;

    result = calloc(1, sizeof(*result));
    read1 = malloc((size_t)pipeline->opt->length[0]);
    read2 = malloc((size_t)pipeline->opt->length[1]);
    quality1 = malloc((size_t)pipeline->opt->length[0]);
    quality2 = malloc((size_t)pipeline->opt->length[1]);
    if(NULL == result || NULL == read1 || NULL == read2 ||
       NULL == quality1 || NULL == quality2) goto error;

    result->task_id = task->id;
    for(sample = 0; sample < pipeline->sample_count; ++sample) {
        uint32_t index;
        result->pair_count[sample] = task->pair_count[sample];
        for(index = 0; index < task->pair_count[sample]; ++index) {
            uint64_t pair_index = task->first_pair[sample] + index;
            if(0 != generate_pair(pipeline->opt, contig, task->contig_index,
                                  pair_index, sample, read1, read2,
                                  quality1, quality2,
                                  &raw[sample * 2],
                                  &raw[sample * 2 + 1])) goto error;
        }
    }

    for(stream = 0; stream < pipeline->stream_count; ++stream) {
        if(0 != compress_bgzf_buffer(&raw[stream],
                                     pipeline->opt->compression_level,
                                     &result->compressed[stream])) goto error;
    }

    for(stream = 0; stream < MAX_STREAMS; ++stream) {
        buffer_destroy(&raw[stream]);
    }
    free(read1);
    free(read2);
    free(quality1);
    free(quality2);
    return result;

error:
    for(stream = 0; stream < MAX_STREAMS; ++stream) {
        buffer_destroy(&raw[stream]);
    }
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
    if(0 == length) return 0;
    if(NULL == data) return -1;
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

        if(0 != write_bytes(arg->output, result->compressed[arg->stream].data,
                            result->compressed[arg->stream].length)) {
            pipeline_fail(pipeline, "failed to append compressed FASTQ task");
            goto close_output;
        }

        pthread_mutex_lock(&pipeline->mutex);
        slot->consumed |= 1U << arg->stream;
        if((1U << pipeline->stream_count) - 1U == slot->consumed) {
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
             FILE **outputs, unsigned sample_count)
{
    pipeline_t pipeline;
    appender_arg_t appender[MAX_STREAMS];
    pthread_t appender_threads[MAX_STREAMS];
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
    pipeline.sample_count = sample_count;
    pipeline.stream_count = sample_count * 2;
    pipeline.window = worker_count * 2;
    if(pipeline.window < 8) pipeline.window = 8;
    if(tasks->count < pipeline.window) pipeline.window = tasks->count;
    pipeline.slots = calloc(pipeline.window, sizeof(*pipeline.slots));
    if(NULL == pipeline.slots ||
       0 != pthread_mutex_init(&pipeline.mutex, NULL) ||
       0 != pthread_cond_init(&pipeline.changed, NULL)) {
        free(pipeline.slots);
        for(index = 0; index < pipeline.stream_count; ++index) {
            if(NULL != outputs[index]) fclose(outputs[index]);
        }
        return -1;
    }

    if(tasks->count < worker_count) worker_count = tasks->count;
    if(0 == worker_count) worker_count = 1;
    workers = calloc(worker_count, sizeof(*workers));
    if(NULL == workers) {
        pipeline_fail(&pipeline, "could not allocate generation threads");
    }

    memset(appender, 0, sizeof(appender));
    for(index = 0; index < pipeline.stream_count; ++index) {
        appender[index].pipeline = &pipeline;
        appender[index].stream = (unsigned)index;
        appender[index].output = outputs[index];
    }

    if(NULL != workers) {
        for(index = 0; index < pipeline.stream_count; ++index) {
            if(0 != pthread_create(&appender_threads[index], NULL,
                                   ordered_appender, &appender[index])) {
                pipeline_fail(&pipeline,
                              "could not create FASTQ appender thread");
                break;
            }
            appenders_created++;
        }
    }
    while(NULL != workers &&
          appenders_created == (int)pipeline.stream_count &&
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
    for(index = (size_t)appenders_created;
        index < pipeline.stream_count; ++index) {
        if(NULL != appender[index].output) fclose(appender[index].output);
        appender[index].output = NULL;
    }

    if(!pipeline.failed && pipeline.released == tasks->count) {
        status = 0;
        for(index = 0; index < pipeline.stream_count; ++index) {
            if(0 != appender[index].status) status = -1;
        }
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

static char
decoded_base(uint8_t base)
{
    static const char bases[] = "ACGTN";
    return bases[base < 5 ? base : 4];
}

static const char *
variant_type_name(uint8_t type)
{
    switch(type) {
      case VARIANT_SNV: return "SNV";
      case VARIANT_INSERTION: return "INS";
      default: return "DEL";
    }
}

static const char *
phased_genotype(uint8_t haplotype_mask)
{
    switch(haplotype_mask) {
      case 1: return "1|0";
      case 2: return "0|1";
      default: return "1|1";
    }
}

static int
write_vcf_header(FILE *output, const reference_manifest_t *reference,
                 const dwgsim_opt_t *opt, uint8_t scope)
{
    size_t contig_index;

    if(0 > fprintf(output,
            "##fileformat=VCFv4.2\n"
            "##source=DWGSIM-deterministic-matched-v%d\n"
            "##seed=%d\n"
            "##germline_event_rate=%.17g\n"
            "##somatic_event_rate=%.17g\n"
            "##tumor_expected_vaf=%.17g\n"
            "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant event type\">\n"
            "##INFO=<ID=SCOPE,Number=1,Type=String,Description=\"Biological variant scope\">\n"
            "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Tumor-only somatic event\">\n"
            "##INFO=<ID=EXPECTED_VAF,Number=1,Type=Float,Description=\"Expected tumor alternate allele fraction\">\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotype\">\n",
            MATCHED_FORMAT_VERSION, opt->seed, opt->mut_rate,
            opt->somatic_rate, opt->tumor_vaf)) return -1;
    for(contig_index = 0; contig_index < reference->count; ++contig_index) {
        const reference_contig_t *contig = &reference->contigs[contig_index];
        if(0 > fprintf(output, "##contig=<ID=%s,length=%" PRIu64 ">\n",
                       contig->name, contig->length)) return -1;
    }
    if(0 > fprintf(output,
            "##truth_scope=%s\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n",
            VARIANT_GERMLINE == scope ? "GERMLINE" : "SOMATIC")) return -1;
    return 0;
}

static int
write_vcf_event(FILE *output, const reference_contig_t *contig,
                const variant_event_t *event, uint64_t ordinal,
                const dwgsim_opt_t *opt)
{
    uint64_t vcf_position = VARIANT_DELETION == event->type ?
                            event->position : (uint64_t)event->position + 1;
    uint32_t index;
    char anchor;

    if(VARIANT_DELETION == event->type) {
        anchor = decoded_base(contig->bases[event->position - 1]);
    }
    else {
        anchor = decoded_base(contig->bases[event->position]);
    }
    if(0 > fprintf(output, "%s\t%" PRIu64 "\t%s.%" PRIu64 "\t",
                   contig->name, vcf_position,
                   VARIANT_GERMLINE == event->scope ? "germline" : "somatic",
                   ordinal)) return -1;
    if(EOF == fputc(anchor, output)) return -1;
    if(VARIANT_DELETION == event->type) {
        for(index = 0; index < event->length; ++index) {
            if(EOF == fputc(decoded_base(contig->bases[event->position + index]),
                            output)) return -1;
        }
    }
    if(EOF == fputc('\t', output)) return -1;
    if(VARIANT_SNV == event->type) {
        if(EOF == fputc(decoded_base(event->alternate_base), output)) return -1;
    }
    else if(VARIANT_INSERTION == event->type) {
        if(EOF == fputc(anchor, output)) return -1;
        for(index = 0; index < event->length; ++index) {
            if(EOF == fputc(decoded_base(event->inserted_bases[index]), output)) {
                return -1;
            }
        }
    }
    else if(EOF == fputc(anchor, output)) {
        return -1;
    }
    if(VARIANT_GERMLINE == event->scope) {
        if(0 > fprintf(output,
                "\t.\tPASS\tTYPE=%s;SCOPE=GERMLINE\tGT\t%s\t%s\n",
                variant_type_name(event->type),
                phased_genotype(event->haplotype_mask),
                phased_genotype(event->haplotype_mask))) return -1;
    }
    else if(0 > fprintf(output,
             "\t.\tPASS\tSOMATIC;TYPE=%s;SCOPE=SOMATIC;EXPECTED_VAF=%.17g"
             "\tGT\t0|0\t%s\n",
             variant_type_name(event->type), opt->tumor_vaf,
             phased_genotype(event->haplotype_mask))) {
        return -1;
    }
    return 0;
}

static int
close_synced_file(FILE **output)
{
    int status = 0;

    if(NULL == *output) return 0;
    if(0 != fflush(*output) || 0 != fsync(fileno(*output))) status = -1;
    if(0 != fclose(*output)) status = -1;
    *output = NULL;
    return status;
}

static int
write_truth_vcfs(const reference_manifest_t *reference,
                 const dwgsim_opt_t *opt,
                 const char *germline_path, const char *somatic_path,
                 uint64_t *germline_count, uint64_t *somatic_count)
{
    FILE *germline = NULL;
    FILE *somatic = NULL;
    size_t contig_index;
    int status = -1;

    *germline_count = 0;
    *somatic_count = 0;
    germline = fopen(germline_path, "w");
    somatic = fopen(somatic_path, "w");
    if(NULL == germline || NULL == somatic ||
       0 != write_vcf_header(germline, reference, opt, VARIANT_GERMLINE) ||
       0 != write_vcf_header(somatic, reference, opt, VARIANT_SOMATIC)) {
        goto cleanup;
    }
    for(contig_index = 0; contig_index < reference->count; ++contig_index) {
        const reference_contig_t *contig = &reference->contigs[contig_index];
        size_t event_index;
        for(event_index = 0; event_index < contig->variants.count;
            ++event_index) {
            const variant_event_t *event = &contig->variants.events[event_index];
            FILE *output;
            uint64_t ordinal;
            if(VARIANT_GERMLINE == event->scope) {
                output = germline;
                ordinal = ++*germline_count;
            }
            else {
                output = somatic;
                ordinal = ++*somatic_count;
            }
            if(0 != write_vcf_event(output, contig, event, ordinal, opt)) {
                goto cleanup;
            }
        }
    }
    if(0 != close_synced_file(&germline) ||
       0 != close_synced_file(&somatic)) goto cleanup;
    status = 0;

cleanup:
    if(NULL != germline) fclose(germline);
    if(NULL != somatic) fclose(somatic);
    if(0 != status) {
        unlink(germline_path);
        unlink(somatic_path);
    }
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

static int
publish_matched_outputs(char **staging, char **outputs,
                        char **truth_staging, char **truth_outputs,
                        const char *manifest_path,
                        const char *manifest_staging,
                        const dwgsim_opt_t *opt, size_t task_count,
                        uint64_t germline_count, uint64_t somatic_count)
{
    static const char *stream_names[MAX_STREAMS] = {
        "normal_read1", "normal_read2", "tumor_read1", "tumor_read2"
    };
    struct stat stream_stats[MAX_STREAMS];
    struct stat truth_stats[2];
    FILE *manifest = NULL;
    unsigned renamed_streams = 0;
    unsigned renamed_truth = 0;
    unsigned index;
    int status = -1;

    for(index = 0; index < MAX_STREAMS; ++index) {
        if(0 != rename(staging[index], outputs[index])) {
            fprintf(stderr, "[dwgsim_matched] cannot publish %s: %s\n",
                    outputs[index], strerror(errno));
            goto cleanup;
        }
        renamed_streams++;
    }
    for(index = 0; index < 2; ++index) {
        if(0 != rename(truth_staging[index], truth_outputs[index])) {
            fprintf(stderr, "[dwgsim_matched] cannot publish %s: %s\n",
                    truth_outputs[index], strerror(errno));
            goto cleanup;
        }
        renamed_truth++;
    }
    for(index = 0; index < MAX_STREAMS; ++index) {
        if(0 != stat(outputs[index], &stream_stats[index])) goto cleanup;
    }
    for(index = 0; index < 2; ++index) {
        if(0 != stat(truth_outputs[index], &truth_stats[index])) goto cleanup;
    }

    manifest = fopen(manifest_staging, "w");
    if(NULL == manifest ||
       0 > fprintf(manifest,
                   "format=dwgsim-deterministic-matched-v%d\n"
                   "seed=%d\nnormal_read_pairs=%" PRId64 "\n"
                   "tumor_read_pairs=%" PRId64 "\nthreads=%d\n"
                   "compression_level=%d\ntasks=%zu\n"
                   "germline_event_rate=%.17g\n"
                   "somatic_event_rate=%.17g\n"
                   "tumor_expected_vaf=%.17g\n"
                   "germline_variants=%" PRIu64 "\n"
                   "somatic_variants=%" PRIu64 "\n",
                   MATCHED_FORMAT_VERSION, opt->seed,
                   opt->normal_pairs, opt->tumor_pairs,
                   opt->compression_threads, opt->compression_level,
                   task_count, opt->mut_rate, opt->somatic_rate,
                   opt->tumor_vaf, germline_count, somatic_count)) {
        goto cleanup;
    }
    for(index = 0; index < MAX_STREAMS; ++index) {
        if(0 > fprintf(manifest, "%s_bytes=%" PRIu64 "\n",
                       stream_names[index],
                       (uint64_t)stream_stats[index].st_size)) goto cleanup;
    }
    if(0 > fprintf(manifest,
                   "germline_vcf_bytes=%" PRIu64 "\n"
                   "somatic_vcf_bytes=%" PRIu64 "\n",
                   (uint64_t)truth_stats[0].st_size,
                   (uint64_t)truth_stats[1].st_size) ||
       0 != fflush(manifest) || 0 != fsync(fileno(manifest))) {
        goto cleanup;
    }
    if(0 != fclose(manifest)) {
        manifest = NULL;
        goto cleanup;
    }
    manifest = NULL;
    if(0 != rename(manifest_staging, manifest_path)) goto cleanup;
    status = 0;

cleanup:
    if(NULL != manifest) fclose(manifest);
    if(0 != status) {
        unlink(manifest_staging);
        for(index = 0; index < renamed_streams; ++index) {
            unlink(outputs[index]);
        }
        for(index = 0; index < renamed_truth; ++index) {
            unlink(truth_outputs[index]);
        }
    }
    return status;
}

int
dwgsim_parallel_wgs_supported(const dwgsim_opt_t *opt,
                              const char *reference_path)
{
    char *fai_path;
    int supported;

    if(NULL == opt || NULL == reference_path ||
       0 == strcmp(reference_path, "-")) {
        return 0;
    }
    if(opt->matched) {
        if(OUTPUT_TYPE_ALL != opt->output_type ||
           READS_OUTPUT_TYPE_BWA != opt->reads_output_type ||
           0 != opt->data_type || 0 >= opt->length[1] ||
           opt->normal_pairs <= 0 || opt->tumor_pairs <= 0 ||
           NULL != opt->fn_muts_input || 0.0 != opt->rand_read ||
           NULL != opt->fn_regions_bed || opt->amplicons ||
           opt->is_inner || opt->is_hap ||
           (0 != opt->strandedness && 2 != opt->strandedness)) {
            return 0;
        }
    }
    else if(OUTPUT_TYPE_READS != opt->output_type ||
            READS_OUTPUT_TYPE_BWA != opt->reads_output_type ||
            0 != opt->data_type || 0 >= opt->length[1] ||
            0 >= opt->N || 0.0 != opt->mut_rate ||
            NULL != opt->fn_muts_input || 0.0 != opt->rand_read ||
            NULL != opt->fn_regions_bed || opt->amplicons ||
            opt->is_inner ||
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

static int
run_matched_wgs(const dwgsim_opt_t *opt,
                const char *reference_path,
                const char *output_prefix)
{
    static const char *stream_suffixes[MAX_STREAMS] = {
        ".normal.bwa.read1.fastq.gz",
        ".normal.bwa.read2.fastq.gz",
        ".tumor.bwa.read1.fastq.gz",
        ".tumor.bwa.read2.fastq.gz"
    };
    static const char *truth_suffixes[2] = {
        ".germline.vcf", ".somatic.vcf"
    };
    reference_manifest_t reference = {0};
    task_manifest_t tasks = {0};
    char *outputs[MAX_STREAMS] = {0};
    char *staging[MAX_STREAMS] = {0};
    char *truth_outputs[2] = {0};
    char *truth_staging[2] = {0};
    FILE *files[MAX_STREAMS] = {0};
    char *manifest_path = NULL;
    char *manifest_staging = NULL;
    char staging_suffix[64];
    uint64_t germline_count = 0;
    uint64_t somatic_count = 0;
    struct timespec start;
    struct timespec finish;
    uint64_t total_pairs;
    double elapsed;
    unsigned index;
    int status = 1;

    clock_gettime(CLOCK_MONOTONIC, &start);
    snprintf(staging_suffix, sizeof(staging_suffix), ".partial.%ld",
             (long)getpid());
    for(index = 0; index < MAX_STREAMS; ++index) {
        outputs[index] = path_with_suffix(output_prefix, stream_suffixes[index]);
        if(NULL != outputs[index]) {
            staging[index] = path_with_suffix(outputs[index], staging_suffix);
        }
        if(NULL == outputs[index] || NULL == staging[index]) goto cleanup;
    }
    for(index = 0; index < 2; ++index) {
        truth_outputs[index] =
            path_with_suffix(output_prefix, truth_suffixes[index]);
        if(NULL != truth_outputs[index]) {
            truth_staging[index] =
                path_with_suffix(truth_outputs[index], staging_suffix);
        }
        if(NULL == truth_outputs[index] || NULL == truth_staging[index]) {
            goto cleanup;
        }
    }
    manifest_path = path_with_suffix(output_prefix, ".matched.complete");
    if(NULL != manifest_path) {
        manifest_staging = path_with_suffix(manifest_path, staging_suffix);
    }
    if(NULL == manifest_path || NULL == manifest_staging) goto cleanup;
    unlink(manifest_path);

    fprintf(stderr,
            "[dwgsim_matched] deterministic matched v%d, %d workers, "
            "BGZF level %d\n",
            MATCHED_FORMAT_VERSION, opt->compression_threads,
            opt->compression_level);
    if(0 != read_reference_manifest(reference_path, &reference)) goto cleanup;
    fprintf(stderr, "[dwgsim_matched] %zu sequences, total length: %" PRIu64
                    "\n", reference.count, reference.total_length);
    if(0 != prepare_reference(opt, reference_path, &reference) ||
       0 != prepare_variants(opt, &reference) ||
       0 != write_truth_vcfs(&reference, opt,
                             truth_staging[0], truth_staging[1],
                             &germline_count, &somatic_count) ||
       0 != allocate_pairs(&reference, (uint64_t)opt->normal_pairs,
                           SAMPLE_NORMAL) ||
       0 != allocate_pairs(&reference, (uint64_t)opt->tumor_pairs,
                           SAMPLE_TUMOR) ||
       0 != build_tasks(&reference, &tasks, MAX_SAMPLES)) {
        goto cleanup;
    }
    discard_unused_contigs(&reference, MAX_SAMPLES);

    fprintf(stderr,
            "[dwgsim_matched] planned %" PRId64 " normal and %" PRId64
            " tumor pairs in %zu fixed tasks; %" PRIu64
            " germline and %" PRIu64 " somatic variants\n",
            opt->normal_pairs, opt->tumor_pairs, tasks.count,
            germline_count, somatic_count);
    for(index = 0; index < MAX_STREAMS; ++index) {
        files[index] = fopen(staging[index], "wb");
        if(NULL == files[index]) {
            fprintf(stderr,
                    "[dwgsim_matched] cannot create staging FASTQ: %s\n",
                    strerror(errno));
            goto cleanup;
        }
    }
    if(0 != run_pipeline(opt, &reference, &tasks, files, MAX_SAMPLES)) {
        for(index = 0; index < MAX_STREAMS; ++index) files[index] = NULL;
        goto cleanup;
    }
    for(index = 0; index < MAX_STREAMS; ++index) files[index] = NULL;

    if(0 != publish_matched_outputs(staging, outputs,
                                    truth_staging, truth_outputs,
                                    manifest_path, manifest_staging,
                                    opt, tasks.count,
                                    germline_count, somatic_count)) {
        fprintf(stderr, "[dwgsim_matched] failed to publish matched outputs\n");
        goto cleanup;
    }

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (double)(finish.tv_sec - start.tv_sec) +
              (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    total_pairs = (uint64_t)opt->normal_pairs +
                  (uint64_t)opt->tumor_pairs;
    fprintf(stderr,
            "[dwgsim_matched] complete: %" PRIu64
            " total pairs in %.3f s (%.2f pairs/s)\n",
            total_pairs, elapsed,
            0.0 < elapsed ? (double)total_pairs / elapsed : 0.0);
    status = 0;

cleanup:
    for(index = 0; index < MAX_STREAMS; ++index) {
        if(NULL != files[index]) fclose(files[index]);
    }
    if(0 != status) {
        for(index = 0; index < MAX_STREAMS; ++index) {
            if(NULL != staging[index]) unlink(staging[index]);
        }
        for(index = 0; index < 2; ++index) {
            if(NULL != truth_staging[index]) unlink(truth_staging[index]);
        }
        if(NULL != manifest_staging) unlink(manifest_staging);
    }
    destroy_tasks(&tasks);
    destroy_reference_manifest(&reference);
    for(index = 0; index < MAX_STREAMS; ++index) {
        free(outputs[index]);
        free(staging[index]);
    }
    for(index = 0; index < 2; ++index) {
        free(truth_outputs[index]);
        free(truth_staging[index]);
    }
    free(manifest_path);
    free(manifest_staging);
    return status;
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

    if(opt->matched) {
        return run_matched_wgs(opt, reference_path, output_prefix);
    }

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
       0 != allocate_pairs(&reference, (uint64_t)opt->N, SAMPLE_NORMAL) ||
       0 != build_tasks(&reference, &tasks, 1)) goto cleanup;
    discard_unused_contigs(&reference, 1);

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
    {
        FILE *files[2] = {file1, file2};
        if(0 != run_pipeline(opt, &reference, &tasks, files, 1)) {
            file1 = file2 = NULL;
            goto cleanup;
        }
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
