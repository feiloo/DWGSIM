#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define BAM_FUNMAP 0x4
#define BAM_FREAD1 0x40

typedef struct {
  char **lines;
  size_t length;
  size_t capacity;
} read_group_t;

typedef struct {
  uint64_t groups;
  uint64_t found_correct[2];
  uint64_t is_correct[2];
  uint64_t found_pair_correct;
  uint64_t is_pair_correct;
} totals_t;

typedef struct {
  uint64_t pos[2];
  char *chromosome;
} origin_t;

static void
usage(FILE *stream)
{
  fprintf(stream,
          "Usage: dwgsim_pileup_eval [-bcp] [-g GAP] [-n READ_PAIRS] "
          "[in.sam ...]\n"
          "\t-p\t\tprint incorrect mapped alignments\n"
          "\t-c\t\tcolor-space alignments\n"
          "\t-g\tINT\tmaximum position difference [5]\n"
          "\t-n\tINT\tnumber of raw input paired-end reads (compatibility option)\n"
          "\t-b\t\talignments are from BWA (for color-space data)\n"
          "\t-h\t\tprint this help message\n");
}

static int
parse_uint64(const char *text, uint64_t *value)
{
  char *end = NULL;
  unsigned long long parsed;

  if('\0' == text[0] || '-' == text[0]) return 0;
  errno = 0;
  parsed = strtoull(text, &end, 10);
  if(0 != errno || '\0' != *end) return 0;
  *value = (uint64_t)parsed;
  return 1;
}

static int
parse_int64(const char *text, int64_t *value)
{
  char *end = NULL;
  long long parsed;

  if('\0' == text[0]) return 0;
  errno = 0;
  parsed = strtoll(text, &end, 10);
  if(0 != errno || '\0' != *end) return 0;
  *value = (int64_t)parsed;
  return 1;
}

static int
split_sam_line(char *line, char ***fields_out, size_t *count_out)
{
  char **fields = NULL;
  size_t count = 0, capacity = 16, i;

  fields = malloc(capacity * sizeof(*fields));
  if(NULL == fields) return 0;
  fields[count++] = line;
  for(i = 0; '\0' != line[i]; ++i) {
      if('\t' == line[i]) {
          line[i] = '\0';
          if(count == capacity) {
              char **resized;
              if(SIZE_MAX / 2 < capacity) {
                  free(fields);
                  return 0;
              }
              capacity *= 2;
              resized = realloc(fields, capacity * sizeof(*fields));
              if(NULL == resized) {
                  free(fields);
                  return 0;
              }
              fields = resized;
          }
          fields[count++] = line + i + 1;
      }
      else if('\r' == line[i] || '\n' == line[i]) {
          line[i] = '\0';
          break;
      }
  }
  *fields_out = fields;
  *count_out = count;
  return 1;
}

static int
parse_origin(const char *qname, origin_t *origin)
{
  static const char separators[] = "_::_::_______";
  char *copy, *tokens[14], *save = NULL, *token;
  uint64_t ignored;
  size_t i, separator_index = 0, token_count = 0;

  memset(origin, 0, sizeof(*origin));
  copy = strdup(qname);
  if(NULL == copy) return -1;

  for(i = strlen(copy); 0 < i && separator_index < sizeof(separators) - 1; --i) {
      if(copy[i - 1] == separators[separator_index]) {
          copy[i - 1] = ' ';
          separator_index++;
      }
  }
  if(separator_index != sizeof(separators) - 1) {
      free(copy);
      return 0;
  }

  token = strtok_r(copy, " ", &save);
  while(NULL != token && token_count < 14) {
      tokens[token_count++] = token;
      token = strtok_r(NULL, " ", &save);
  }
  if(14 != token_count || NULL != token) {
      free(copy);
      return 0;
  }
  if(!parse_uint64(tokens[1], &origin->pos[0]) ||
     !parse_uint64(tokens[2], &origin->pos[1])) {
      free(copy);
      return 0;
  }
  for(i = 3; i <= 12; ++i) {
      if(!parse_uint64(tokens[i], &ignored)) {
          free(copy);
          return 0;
      }
  }
  origin->chromosome = strdup(tokens[0]);
  free(copy);
  return NULL == origin->chromosome ? -1 : 1;
}

static int
within_gap(uint64_t expected, uint64_t observed, uint64_t gap)
{
  if(expected < observed) return observed - expected <= gap;
  return expected - observed <= gap;
}

static int
read_group_push(read_group_t *group, const char *line)
{
  char **resized;

  if(group->length == group->capacity) {
      size_t new_capacity = 0 == group->capacity ? 8 : group->capacity * 2;
      if(new_capacity < group->capacity) return 0;
      resized = realloc(group->lines, new_capacity * sizeof(*group->lines));
      if(NULL == resized) return 0;
      group->lines = resized;
      group->capacity = new_capacity;
  }
  group->lines[group->length] = strdup(line);
  if(NULL == group->lines[group->length]) return 0;
  group->length++;
  return 1;
}

static void
read_group_clear(read_group_t *group)
{
  size_t i;
  for(i = 0; i < group->length; ++i) free(group->lines[i]);
  group->length = 0;
}

static void
read_group_destroy(read_group_t *group)
{
  read_group_clear(group);
  free(group->lines);
  group->lines = NULL;
  group->capacity = 0;
}

static int
process_group(const read_group_t *group, totals_t *totals, uint64_t gap,
              int color_space, int bwa, int print_incorrect)
{
  int has_max_score[2] = {0, 0};
  int has_correct_score[2] = {0, 0};
  int64_t max_score[2] = {INT64_MIN, INT64_MIN};
  int64_t best_correct_score[2] = {INT64_MIN, INT64_MIN};
  int found_correct[2] = {0, 0};
  int best_is_correct[2] = {0, 0};
  size_t i;

  for(i = 0; i < group->length; ++i) {
      char *copy = strdup(group->lines[i]);
      char **fields = NULL;
      size_t field_count = 0, j;
      uint64_t flag, left, mapq;
      int64_t alignment_score;
      int end, is_alignment_correct = 0;
      int has_alignment_score = 0;
      origin_t origin;
      int origin_status;

      if(NULL == copy || !split_sam_line(copy, &fields, &field_count)) {
          free(copy);
          free(fields);
          fprintf(stderr, "dwgsim_pileup_eval: out of memory\n");
          return 0;
      }
      if(field_count < 11 || !parse_uint64(fields[1], &flag) ||
         !parse_uint64(fields[3], &left) || !parse_uint64(fields[4], &mapq)) {
          fprintf(stderr, "dwgsim_pileup_eval: malformed SAM record\n");
          free(fields);
          free(copy);
          return 0;
      }
      if((flag & BAM_FUNMAP) || 0 == strcmp(fields[2], "*")) {
          free(fields);
          free(copy);
          continue;
      }

      origin_status = parse_origin(fields[0], &origin);
      if(0 > origin_status) {
          fprintf(stderr, "dwgsim_pileup_eval: out of memory\n");
          free(fields);
          free(copy);
          return 0;
      }
      if(0 == origin_status) {
          fprintf(stderr,
                  "[dwgsim_pileup_eval] read '%s' was not generated by dwgsim?\n",
                  fields[0]);
          free(fields);
          free(copy);
          continue;
      }

      end = (flag & BAM_FREAD1) ? 0 : 1;
      if(color_space && bwa) {
          uint64_t temp = origin.pos[0];
          origin.pos[0] = origin.pos[1];
          origin.pos[1] = temp;
      }
      for(j = 11; j < field_count; ++j) {
          if(0 == strncmp(fields[j], "AS:i:", 5) &&
             parse_int64(fields[j] + 5, &alignment_score)) {
              has_alignment_score = 1;
              if(!has_max_score[end] || max_score[end] < alignment_score) {
                  has_max_score[end] = 1;
                  max_score[end] = alignment_score;
              }
              break;
          }
      }

      if(0 == strcmp(origin.chromosome, fields[2]) &&
         within_gap(origin.pos[end], left, gap)) {
          found_correct[end] = 1;
          if(has_alignment_score &&
             (!has_correct_score[end] || best_correct_score[end] < alignment_score)) {
              has_correct_score[end] = 1;
              best_correct_score[end] = alignment_score;
          }
          is_alignment_correct = 1;
      }
      if(print_incorrect && !is_alignment_correct && 0 < mapq) {
          size_t raw_length = strlen(group->lines[i]);
          fputs(group->lines[i], stderr);
          if(0 == raw_length || '\n' != group->lines[i][raw_length - 1]) {
              fputc('\n', stderr);
          }
      }

      free(origin.chromosome);
      free(fields);
      free(copy);
  }

  totals->groups++;
  for(i = 0; i < 2; ++i) {
      best_is_correct[i] = found_correct[i] && has_max_score[i] &&
                           has_correct_score[i] &&
                           best_correct_score[i] == max_score[i];
      totals->found_correct[i] += (uint64_t)found_correct[i];
      totals->is_correct[i] += (uint64_t)best_is_correct[i];
  }
  if(found_correct[0] && found_correct[1]) totals->found_pair_correct++;
  if(best_is_correct[0] && best_is_correct[1]) {
      totals->is_pair_correct++;
  }
  return 1;
}

static int
line_has_sam_fields(const char *line)
{
  size_t tabs = 0;
  for(; '\0' != *line && '\r' != *line && '\n' != *line; ++line) {
      if('\t' == *line) tabs++;
  }
  return 10 <= tabs;
}

static char *
copy_qname(const char *line)
{
  const char *tab = strchr(line, '\t');
  size_t length;
  char *qname;

  if(NULL == tab) return NULL;
  length = (size_t)(tab - line);
  qname = malloc(length + 1);
  if(NULL == qname) return NULL;
  memcpy(qname, line, length);
  qname[length] = '\0';
  return qname;
}

static int
consume_stream(FILE *input, const char *name, read_group_t *group,
               char **previous_qname, totals_t *totals, uint64_t gap,
               int color_space, int bwa, int print_incorrect,
               uint64_t *line_count)
{
  char *line = NULL;
  size_t capacity = 0;

  while(0 <= getline(&line, &capacity, input)) {
      char *qname;
      (*line_count)++;
      if(0 == *line_count % 10000) fprintf(stderr, "\r%" PRIu64, *line_count);
      if('@' == line[0] || !line_has_sam_fields(line)) continue;

      qname = copy_qname(line);
      if(NULL == qname) {
          fprintf(stderr, "dwgsim_pileup_eval: out of memory\n");
          free(line);
          return 0;
      }
      if(NULL != *previous_qname && 0 != strcmp(*previous_qname, qname)) {
          if(!process_group(group, totals, gap, color_space, bwa,
                            print_incorrect)) {
              free(qname);
              free(line);
              return 0;
          }
          read_group_clear(group);
      }
      if(NULL == *previous_qname || 0 != strcmp(*previous_qname, qname)) {
          free(*previous_qname);
          *previous_qname = qname;
      }
      else {
          free(qname);
      }
      if(!read_group_push(group, line)) {
          fprintf(stderr, "dwgsim_pileup_eval: out of memory\n");
          free(line);
          return 0;
      }
  }
  free(line);
  if(ferror(input)) {
      fprintf(stderr, "dwgsim_pileup_eval: could not read '%s': %s\n",
              name, strerror(errno));
      return 0;
  }
  return 1;
}

int
main(int argc, char *argv[])
{
  int option, bwa = 0, color_space = 0, print_incorrect = 0;
  uint64_t gap = 5, raw_input_pairs = 0, line_count = 0;
  read_group_t group = {0};
  totals_t totals = {0};
  char *previous_qname = NULL;
  int status = EXIT_FAILURE, i;

  while(-1 != (option = getopt(argc, argv, "bcpg:n:h"))) {
      switch(option) {
        case 'b': bwa = 1; break;
        case 'c': color_space = 1; break;
        case 'p': print_incorrect = 1; break;
        case 'g':
          if(!parse_uint64(optarg, &gap)) {
              fprintf(stderr, "dwgsim_pileup_eval: invalid gap '%s'\n", optarg);
              goto cleanup;
          }
          break;
        case 'n':
          if(!parse_uint64(optarg, &raw_input_pairs)) {
              fprintf(stderr, "dwgsim_pileup_eval: invalid read count '%s'\n", optarg);
              goto cleanup;
          }
          break;
        case 'h': usage(stdout); status = EXIT_SUCCESS; goto cleanup;
        default: usage(stderr); goto cleanup;
      }
  }
  (void)raw_input_pairs;

  if(optind == argc && isatty(STDIN_FILENO)) {
      usage(stderr);
      goto cleanup;
  }

  fprintf(stderr, "Analyzing...\nCurrently on:\n0");
  if(optind == argc) {
      if(!consume_stream(stdin, "standard input", &group, &previous_qname,
                         &totals, gap, color_space, bwa, print_incorrect,
                         &line_count)) {
          goto cleanup;
      }
  }
  else {
      for(i = optind; i < argc; ++i) {
          FILE *input;
          if(0 == strcmp(argv[i], "-")) {
              input = stdin;
          }
          else {
              input = fopen(argv[i], "r");
              if(NULL == input) {
                  fprintf(stderr, "dwgsim_pileup_eval: could not open '%s': %s\n",
                          argv[i], strerror(errno));
                  goto cleanup;
              }
          }
          if(!consume_stream(input, argv[i], &group, &previous_qname,
                             &totals, gap, color_space, bwa, print_incorrect,
                             &line_count)) {
              if(stdin != input) fclose(input);
              goto cleanup;
          }
          if(stdin != input && 0 != fclose(input)) {
              fprintf(stderr, "dwgsim_pileup_eval: could not close '%s': %s\n",
                      argv[i], strerror(errno));
              goto cleanup;
          }
      }
  }
  if(0 < group.length &&
     !process_group(&group, &totals, gap, color_space, bwa, print_incorrect)) {
      goto cleanup;
  }

  fprintf(stderr, "\r%" PRIu64 "\n", line_count);
  fprintf(stderr,
          "%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64
          "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
          totals.groups, totals.found_correct[0], totals.found_correct[1],
          totals.is_correct[0], totals.is_correct[1],
          totals.found_pair_correct, totals.is_pair_correct);
  fprintf(stderr, "Analysis complete.\n");
  status = EXIT_SUCCESS;

cleanup:
  free(previous_qname);
  read_group_destroy(&group);
  return status;
}
