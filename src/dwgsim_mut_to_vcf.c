#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  char *name;
  uint64_t pos;
  char *ref;
  uint64_t strand;
  size_t length;
} deletion_t;

static void
usage(FILE *stream)
{
  fprintf(stream, "Usage: dwgsim_mut_to_vcf <in.mutations.txt>\n");
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
split_mutation_line(char *line, char *fields[5])
{
  size_t i, n = 1;

  line[strcspn(line, "\r\n")] = '\0';
  fields[0] = line;
  for(i = 0; '\0' != line[i]; ++i) {
      if('\t' == line[i]) {
          if(5 <= n) return 0;
          line[i] = '\0';
          fields[n++] = line + i + 1;
      }
  }
  if(5 != n) return 0;
  for(i = 0; i < 5; ++i) {
      if('\0' == fields[i][0]) return 0;
  }
  return 1;
}

static const char *
converted_alt(const char *ref, const char *alt)
{
  static const struct {
    const char *key;
    const char *alt;
  } conversions[] = {
    {"AR", "G"}, {"GR", "A"},
    {"CY", "T"}, {"TY", "C"},
    {"CS", "G"}, {"GS", "C"},
    {"AW", "T"}, {"TW", "A"},
    {"GK", "T"}, {"TK", "G"},
    {"AM", "C"}, {"CM", "A"}
  };
  char key[3];
  size_t i;

  if(1 != strlen(ref) || 1 != strlen(alt)) return NULL;
  key[0] = ref[0];
  key[1] = alt[0];
  key[2] = '\0';
  for(i = 0; i < sizeof(conversions) / sizeof(conversions[0]); ++i) {
      if(0 == strcmp(key, conversions[i].key)) return conversions[i].alt;
  }
  return NULL;
}

static void
clear_deletion(deletion_t *deletion)
{
  free(deletion->name);
  free(deletion->ref);
  deletion->name = NULL;
  deletion->ref = NULL;
  deletion->pos = 0;
  deletion->strand = 0;
  deletion->length = 0;
}

static int
flush_deletion(deletion_t *deletion)
{
  if(0 == deletion->length) return 1;
  if(0 > fprintf(stdout, "%s\t%" PRIu64 "\t.\t%s\t.\t100\tPASS\tpl=%" PRIu64 "\n",
                 deletion->name, deletion->pos, deletion->ref,
                 deletion->strand)) {
      return 0;
  }
  clear_deletion(deletion);
  return 1;
}

static int
start_deletion(deletion_t *deletion, const char *name, uint64_t pos,
               const char *ref, uint64_t strand)
{
  deletion->name = strdup(name);
  deletion->ref = strdup(ref);
  if(NULL == deletion->name || NULL == deletion->ref) {
      clear_deletion(deletion);
      return 0;
  }
  deletion->pos = pos;
  deletion->strand = strand;
  deletion->length = 1;
  return 1;
}

static int
append_deletion(deletion_t *deletion, const char *ref)
{
  size_t old_length = strlen(deletion->ref);
  size_t added_length = strlen(ref);
  char *new_ref;

  if(SIZE_MAX - old_length <= added_length) return 0;
  new_ref = realloc(deletion->ref, old_length + added_length + 1);
  if(NULL == new_ref) return 0;
  memcpy(new_ref + old_length, ref, added_length + 1);
  deletion->ref = new_ref;
  deletion->length++;
  return 1;
}

int
main(int argc, char *argv[])
{
  FILE *input;
  char *line = NULL;
  size_t line_capacity = 0;
  ssize_t line_length;
  uint64_t line_number = 0;
  deletion_t deletion = {0};
  int status = EXIT_FAILURE;

  if(2 != argc) {
      usage(stderr);
      return EXIT_FAILURE;
  }
  input = fopen(argv[1], "r");
  if(NULL == input) {
      fprintf(stderr, "dwgsim_mut_to_vcf: could not open '%s': %s\n",
              argv[1], strerror(errno));
      return EXIT_FAILURE;
  }

  printf("##fileformat=VCFv4.1\n");
  printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample\n");

  while(0 <= (line_length = getline(&line, &line_capacity, input))) {
      char *fields[5];
      const char *name, *ref, *alt, *output_alt;
      uint64_t pos, strand;
      int continues_deletion;

      (void)line_length;
      line_number++;
      if(!split_mutation_line(line, fields) ||
         !parse_uint64(fields[1], &pos) ||
         !parse_uint64(fields[4], &strand)) {
          fprintf(stderr,
                  "dwgsim_mut_to_vcf: input line %" PRIu64
                  " is not in the proper format\n",
                  line_number);
          goto cleanup;
      }
      name = fields[0];
      ref = fields[2];
      alt = fields[3];

      if(0 == strcmp(alt, "-")) {
          continues_deletion = 0 < deletion.length &&
                               0 == strcmp(deletion.name, name) &&
                               deletion.pos + deletion.length == pos &&
                               deletion.strand == strand;
          if(continues_deletion) {
              if(!append_deletion(&deletion, ref)) {
                  fprintf(stderr, "dwgsim_mut_to_vcf: out of memory\n");
                  goto cleanup;
              }
          }
          else {
              if(!flush_deletion(&deletion)) goto write_error;
              if(!start_deletion(&deletion, name, pos, ref, strand)) {
                  fprintf(stderr, "dwgsim_mut_to_vcf: out of memory\n");
                  goto cleanup;
              }
          }
          continue;
      }

      if(!flush_deletion(&deletion)) goto write_error;
      if(0 == strcmp(ref, "-")) {
          if(0 > printf("%s\t%" PRIu64 "\t.\t.\t%s\t100\tPASS\tpl=%" PRIu64 "\n",
                        name, pos, alt, strand)) {
              goto write_error;
          }
          continue;
      }

      output_alt = alt;
      if(3 != strand) {
          output_alt = converted_alt(ref, alt);
          if(NULL == output_alt) {
              fprintf(stderr,
                      "dwgsim_mut_to_vcf: could not convert %s%s on line %" PRIu64 "\n",
                      ref, alt, line_number);
              goto cleanup;
          }
      }
      if(0 > printf("%s\t%" PRIu64 "\t.\t%s\t%s\t100\tPASS\tpl=%" PRIu64 "\n",
                    name, pos, ref, output_alt, strand)) {
          goto write_error;
      }
  }

  if(ferror(input)) {
      fprintf(stderr, "dwgsim_mut_to_vcf: could not read '%s': %s\n",
              argv[1], strerror(errno));
      goto cleanup;
  }
  if(!flush_deletion(&deletion)) goto write_error;
  if(0 != fflush(stdout)) goto write_error;
  status = EXIT_SUCCESS;
  goto cleanup;

write_error:
  fprintf(stderr, "dwgsim_mut_to_vcf: could not write output: %s\n",
          strerror(errno));

cleanup:
  clear_deletion(&deletion);
  free(line);
  fclose(input);
  return status;
}
