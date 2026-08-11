#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fastq_writer.h"

#define FASTQ_BGZF_MODE "w1"
#define FASTQ_BGZF_SUB_BLOCKS 64

fastq_writer_t *
fastq_writer_open(const char *path, int compression_threads)
{
    fastq_writer_t *writer;

    if(NULL == path || compression_threads < 1) return NULL;

    writer = calloc(1, sizeof(*writer));
    if(NULL == writer) return NULL;

    writer->fp = bgzf_open(path, FASTQ_BGZF_MODE);
    if(NULL == writer->fp) {
        free(writer);
        return NULL;
    }

    if(1 < compression_threads &&
       0 != bgzf_mt(writer->fp, compression_threads,
                    FASTQ_BGZF_SUB_BLOCKS)) {
        bgzf_close(writer->fp);
        free(writer);
        return NULL;
    }

    return writer;
}

int
fastq_writer_flush(fastq_writer_t *writer)
{
    ssize_t written;

    if(NULL == writer || writer->error) return -1;
    if(0 == writer->length) return 0;

    written = bgzf_write(writer->fp, writer->buffer, writer->length);
    if(written != (ssize_t)writer->length) {
        writer->error = 1;
        return -1;
    }

    writer->length = 0;
    return 0;
}

int
fastq_writer_write(fastq_writer_t *writer, const void *data, size_t length)
{
    const unsigned char *input = data;

    if(NULL == writer || (NULL == data && 0 < length) || writer->error) {
        return -1;
    }

    while(0 < length) {
        size_t available = sizeof(writer->buffer) - writer->length;
        size_t to_copy;

        if(0 == available) {
            if(0 != fastq_writer_flush(writer)) return -1;
            available = sizeof(writer->buffer);
        }

        to_copy = length < available ? length : available;
        memcpy(writer->buffer + writer->length, input, to_copy);
        writer->length += to_copy;
        input += to_copy;
        length -= to_copy;
    }

    return 0;
}

int
fastq_writer_printf(fastq_writer_t *writer, const char *format, ...)
{
    va_list args;
    va_list args_copy;
    size_t available;
    int needed;

    if(NULL == writer || NULL == format || writer->error) return -1;
    if(sizeof(writer->buffer) == writer->length &&
       0 != fastq_writer_flush(writer)) {
        return -1;
    }

    available = sizeof(writer->buffer) - writer->length;
    va_start(args, format);
    va_copy(args_copy, args);
    needed = vsnprintf((char *)writer->buffer + writer->length,
                       available, format, args);
    va_end(args);

    if(needed < 0) {
        va_end(args_copy);
        writer->error = 1;
        return -1;
    }

    if((size_t)needed < available) {
        writer->length += needed;
        va_end(args_copy);
        return needed;
    }

    if(0 != fastq_writer_flush(writer)) {
        va_end(args_copy);
        return -1;
    }

    if((size_t)needed < sizeof(writer->buffer)) {
        int second_result = vsnprintf((char *)writer->buffer,
                                      sizeof(writer->buffer), format,
                                      args_copy);
        va_end(args_copy);
        if(second_result != needed) {
            writer->error = 1;
            return -1;
        }
        writer->length = needed;
    }
    else {
        char *dynamic_buffer = malloc((size_t)needed + 1);
        int second_result;

        if(NULL == dynamic_buffer) {
            va_end(args_copy);
            writer->error = 1;
            return -1;
        }
        second_result = vsnprintf(dynamic_buffer, (size_t)needed + 1,
                                  format, args_copy);
        va_end(args_copy);
        if(second_result != needed ||
           0 != fastq_writer_write(writer, dynamic_buffer,
                                   (size_t)needed)) {
            free(dynamic_buffer);
            writer->error = 1;
            return -1;
        }
        free(dynamic_buffer);
    }

    return needed;
}

int
fastq_writer_close(fastq_writer_t *writer)
{
    int status = 0;

    if(NULL == writer) return 0;
    if(0 != fastq_writer_flush(writer)) status = -1;
    if(0 != bgzf_close(writer->fp)) status = -1;
    free(writer);
    return status;
}

int
fastq_writer_threads_for_stream(int total_threads, int stream_count,
                                int stream_index)
{
    int helper_threads;

    if(total_threads < 1 || stream_count < 1 || stream_index < 0 ||
       stream_count <= stream_index) {
        return 0;
    }

    helper_threads = total_threads - 1;
    return 1 + helper_threads / stream_count +
           (stream_index < helper_threads % stream_count);
}
