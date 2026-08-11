/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   */

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>
#include "contigs.h"
#include "regions_bed.h"

static char *
regions_bed_next_field(char **cursor)
{
  char *field, *p = *cursor;

  while('\0' != *p && isspace((unsigned char)*p)) p++;
  if('\0' == *p) {
      *cursor = p;
      return NULL;
  }

  field = p;
  while('\0' != *p && !isspace((unsigned char)*p)) p++;
  if('\0' != *p) *p++ = '\0';
  *cursor = p;
  return field;
}

static int32_t
regions_bed_is_header(const char *line)
{
  const char *p = line;

  while('\0' != *p && isspace((unsigned char)*p)) p++;
  if('\0' == *p || '#' == *p) return 1;
  if(0 == strncmp(p, "track", 5) &&
     ('\0' == p[5] || isspace((unsigned char)p[5]))) return 1;
  if(0 == strncmp(p, "browser", 7) &&
     ('\0' == p[7] || isspace((unsigned char)p[7]))) return 1;
  return 0;
}

static int32_t
regions_bed_parse_coordinate(const char *text, uint32_t *value)
{
  char *end = NULL;
  unsigned long long parsed;

  if(NULL == text || '\0' == text[0] || '-' == text[0] || '+' == text[0]) {
      return 0;
  }
  errno = 0;
  parsed = strtoull(text, &end, 10);
  if(0 != errno || end == text || '\0' != *end || UINT32_MAX < parsed) {
      return 0;
  }
  *value = (uint32_t)parsed;
  return 1;
}

static int32_t
regions_bed_find_contig(contigs_t *c, const char *name)
{
  int32_t i;

  for(i=0;i<c->n;i++) {
      if(0 == strcmp(name, c->contigs[i].name)) return i;
  }
  return -1;
}

static void
regions_bed_grow(regions_bed_txt *r)
{
  uint32_t *temp_contig, *temp_start, *temp_end;

  r->mem <<= 1;
  temp_contig = realloc(r->contig, r->mem * sizeof(uint32_t));
  temp_start = realloc(r->start, r->mem * sizeof(uint32_t));
  temp_end = realloc(r->end, r->mem * sizeof(uint32_t));
  if(NULL == temp_contig || NULL == temp_start || NULL == temp_end) {
      fprintf(stderr, "Error: memory allocation failed in regions_bed_init\n");
      exit(1);
  }
  r->contig = temp_contig;
  r->start = temp_start;
  r->end = temp_end;
}

regions_bed_txt *regions_bed_init(FILE *fp, contigs_t *c)
{
  regions_bed_txt *r = NULL;
  char *line = NULL;
  size_t line_capacity = 0;
  ssize_t line_length;
  uint64_t line_number = 0;
  uint32_t start, end, prev_start = 0, prev_end = 0;
  int32_t contig, prev_contig = -1;

  r = calloc(1, sizeof(regions_bed_txt));
  if(NULL == r) {
      fprintf(stderr, "Error: memory allocation failed in regions_bed_init\n");
      exit(1);
  }
  r->n = 0;
  r->mem = 4;
  r->contig = malloc(r->mem * sizeof(uint32_t));
  r->start = malloc(r->mem * sizeof(uint32_t));
  r->end = malloc(r->mem * sizeof(uint32_t));
  if(NULL == r->contig || NULL == r->start || NULL == r->end) {
      fprintf(stderr, "Error: memory allocation failed in regions_bed_init\n");
      exit(1);
  }

  while(0 <= (line_length = getline(&line, &line_capacity, fp))) {
      char *cursor, *name, *start_text, *end_text;
      (void)line_length;
      line_number++;
      if(regions_bed_is_header(line)) continue;

      cursor = line;
      name = regions_bed_next_field(&cursor);
      start_text = regions_bed_next_field(&cursor);
      end_text = regions_bed_next_field(&cursor);
      if(NULL == name || NULL == start_text || NULL == end_text) {
          fprintf(stderr, "Error: BED line %llu has fewer than 3 fields\n",
                  (unsigned long long)line_number);
          exit(1);
      }
      if(0 == regions_bed_parse_coordinate(start_text, &start) ||
         0 == regions_bed_parse_coordinate(end_text, &end)) {
          fprintf(stderr, "Error: BED line %llu has a non-integer or out-of-range coordinate\n",
                  (unsigned long long)line_number);
          exit(1);
      }
      if(start >= end) {
          fprintf(stderr, "Error: BED line %llu must satisfy start < end [%s,%u,%u]\n",
                  (unsigned long long)line_number, name, start, end);
          exit(1);
      }

      contig = (prev_contig >= 0 &&
                0 == strcmp(name, c->contigs[prev_contig].name))
                   ? prev_contig : regions_bed_find_contig(c, name);
      if(contig < 0) {
          fprintf(stderr, "Error: BED line %llu contig not found in reference [%s]\n",
                  (unsigned long long)line_number, name);
          exit(1);
      }
      if(contig < prev_contig) {
          fprintf(stderr, "Error: BED line %llu is not in reference contig order [%s]\n",
                  (unsigned long long)line_number, name);
          exit(1);
      }
      if(contig == prev_contig && start < prev_start) {
          fprintf(stderr, "Error: BED line %llu is not sorted by start [%s,%u,%u]\n",
                  (unsigned long long)line_number, name, start, end);
          exit(1);
      }
      if(end > (uint32_t)c->contigs[contig].len) {
          fprintf(stderr, "Error: BED line %llu end is outside reference contig [%s,%u,%u; length=%d]\n",
                  (unsigned long long)line_number, name, start, end,
                  c->contigs[contig].len);
          exit(1);
      }

      if(prev_contig == contig && start <= prev_end) {
          if(prev_end < end) {
              r->end[r->n-1] = end;
              prev_end = end;
          }
      }
      else {
          prev_end = end;
          if(r->mem <= r->n) regions_bed_grow(r);
          r->contig[r->n] = contig;
          r->start[r->n] = start;
          r->end[r->n] = end;
          r->n++;
      }
      prev_contig = contig;
      prev_start = start;
  }
  if(ferror(fp)) {
      fprintf(stderr, "Error: failed while reading regions BED after line %llu\n",
              (unsigned long long)line_number);
      exit(1);
  }
  free(line);
  if(0 == r->n) {
      fprintf(stderr, "Error: regions BED contains no intervals\n");
      exit(1);
  }
  return r;
}

void regions_bed_destroy(regions_bed_txt *r)
{
  free(r->contig);
  free(r->start);
  free(r->end);
  free(r);
}

int32_t regions_bed_query(regions_bed_txt *r, uint32_t contig, uint32_t start, uint32_t end)
{
  int32_t low, high, mid;
  if(NULL == r) return 1;

  low = 0;
  high = r->n-1;

  while(low <= high) {
      mid = low + (high - low) / 2;
      if(contig < r->contig[mid] ||
         (contig == r->contig[mid] && start < r->start[mid])) {
          high = mid - 1;
      }
      else if(r->contig[mid] < contig ||
              (r->contig[mid] == contig && r->end[mid] < end)) {
          low = mid + 1;
      }
      else if(r->contig[mid] == contig && r->start[mid] <= start && end <= r->end[mid]) {
          return 1;
      }
      else {
          break;
      }
  }
  return 0;
}

int32_t regions_bed_map_offset(regions_bed_txt *r, uint32_t contig,
                               uint64_t offset, uint32_t *position)
{
  uint32_t i;

  if(NULL == r || NULL == position) return 0;
  for(i=0;i<r->n;i++) {
      uint64_t length;
      if(r->contig[i] < contig) continue;
      if(contig < r->contig[i]) break;
      length = (uint64_t)r->end[i] - r->start[i];
      if(offset < length) {
          *position = r->start[i] + (uint32_t)offset;
          return 1;
      }
      offset -= length;
  }
  return 0;
}
