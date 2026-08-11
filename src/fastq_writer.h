#ifndef FASTQ_WRITER_H
#define FASTQ_WRITER_H

#include <stddef.h>

#include "samtools/bgzf.h"

typedef struct {
    BGZF *fp;
    unsigned char buffer[BGZF_BLOCK_SIZE];
    size_t length;
    int error;
} fastq_writer_t;

fastq_writer_t *fastq_writer_open(const char *path, int compression_threads,
                                  int compression_level);
int fastq_writer_flush(fastq_writer_t *writer);
int fastq_writer_write(fastq_writer_t *writer, const void *data, size_t length);
int fastq_writer_printf(fastq_writer_t *writer, const char *format, ...)
    __attribute__((format(printf, 2, 3)));
int fastq_writer_close(fastq_writer_t *writer);

int fastq_writer_threads_for_stream(int total_threads, int stream_count,
                                    int stream_index);

static inline int
fastq_writer_putc(fastq_writer_t *writer, int value)
{
    if(NULL == writer || writer->error) return -1;
    if(sizeof(writer->buffer) == writer->length &&
       0 != fastq_writer_flush(writer)) {
        return -1;
    }
    writer->buffer[writer->length++] = (unsigned char)value;
    return (unsigned char)value;
}

#endif
