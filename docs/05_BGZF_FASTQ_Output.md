# BGZF FASTQ output plan

## Goals

- Replace zlib's single-stream gzip FASTQ writer with the BGZF implementation already bundled under `samtools/`.
- Keep output readable by ordinary gzip/FASTQ tools while allowing parallel block compression.
- Use compression level 1, which is the best measured speed/size tradeoff for DWGSIM output.
- Add `-t INT` as a total DWGSIM thread budget without multiplying that value for every output file.
- Preserve read names, sequences, qualities, output filenames, and all non-FASTQ output formats.
- Avoid a new runtime dependency on an external `bgzip` executable or modern HTSlib.

## Design

`src/fastq_writer.[ch]` provides the small output abstraction that BGZF itself does not provide. It owns a `BGZF *`, buffers small writes, supports formatted output and individual characters, reports write/close errors, and opens files with BGZF compression level 1. Buffering keeps the existing per-base generation loops inexpensive while passing block-sized writes to BGZF.

The FASTQ handles in `dwgsim_opt_t` become `fastq_writer_t *`. Existing `gzprintf`, `gzputc`, `gzopen`, and `gzclose` calls are replaced by the corresponding writer operations. The `dwgsim` executable links the bundled `libbam.a`, which already contains `bgzf.o` and its required support objects.

`-t` counts the main DWGSIM thread plus BGZF helper threads. Its default is 1. For `T` requested threads and `S` active FASTQ streams, the `T - 1` helpers are distributed across streams as evenly as possible:

```text
helpers = T - 1
stream_threads[i] = 1 + helpers / S + (i < helpers % S)
```

The leading 1 represents the shared caller for a stream; only the remainder creates helper threads. This keeps the process-wide helper count at `T - 1`. With paired BWA output, for example, `-t 4` allocates two helpers to one output and one helper to the other instead of creating eight compression workers.

## Compatibility and tradeoffs

BGZF is a blocked form of gzip, so existing `.fastq.gz` readers continue to work. The compressed bytes and file sizes change, and level 1 usually produces larger files than zlib's former default level 6. BGZF's independent blocks can also be slightly larger than an equivalent continuous gzip stream.

Only FASTQ compression is parallel. Read simulation and mutation generation remain single-threaded. Small jobs may see little benefit because thread setup and BGZF block batching have fixed overhead.

## Validation

1. Build all executables with the existing compiler warnings enabled.
2. Run unit and integration tests.
3. Exercise BWA-only, BFAST-only, and combined FASTQ output modes.
4. Check every output with both ordinary `gzip --test` and BGZF-aware validation.
5. Confirm seeded single- and multi-thread runs decompress to identical FASTQ content.
6. Benchmark `-t 1` and a representative multi-thread setting, reporting reads/s, CPU use, memory, and compressed size.
