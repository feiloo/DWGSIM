# Deterministic parallel generation target

> Status: design target, not the current implementation.
>
> The current implementation parallelizes BGZF compression while read
> simulation and reference traversal remain serial. This document defines the
> target architecture for parallel read generation and deterministic output.

## Goals

- Scale paired-end read generation beyond 100 CPU cores.
- Keep peak resident memory below 500 GiB, with an operational target below
  150 GiB for a 128-worker GRCh38 run.
- Produce exactly two synchronized paired-end FASTQs in BWA output mode.
- Produce byte-identical BGZF FASTQs across repeated runs and thread counts
  when all non-thread inputs are identical.
- Parallelize reference traversal, mutation preparation, read generation,
  quality/error generation, FASTQ formatting, and compression.
- Write directly to the final paired FASTQs without per-shard FASTQ files or
  a whole-output merge pass.
- Keep gzip and BGZF compatibility.
- Make compression level explicit and reproducible.
- Preserve the read-pair contents of the new deterministic single-thread
  baseline in every parallel execution.

The first optimized milestone targets reads-only, BWA-only, paired-end 2x150
bp generation. The same task and RNG model must subsequently support WGS,
BED-restricted WGS, WES-like BED input, mutation input, and mutation output.

## Non-goals

- Reproducing the exact reads or compressed hashes produced by the legacy
  global `drand48()` implementation.
- Preserving the legacy global read order. A new canonical order may be
  established once, but it must then remain stable.
- Guaranteeing identical compressed bytes across different compressor
  versions, compiler floating-point modes, architectures, or format versions.
  Cross-host reproducibility requires pinning those components.
- Using thread completion order as output order.
- Writing uncompressed FASTQ and compressing it in a second pass.

## Determinism contract

For the same executable and format version, reference bytes and index, BED or
mutation inputs, random seed, simulation options, and compression level, the
following must be true:

```text
sha256(R1, -t 1)   == sha256(R1, -t 8)
sha256(R1, -t 8)   == sha256(R1, -t 128)
sha256(R2, -t 1)   == sha256(R2, -t 8)
sha256(R2, -t 8)   == sha256(R2, -t 128)
sha256(repeat one) == sha256(repeat two)
```

Thread count is a scheduling parameter only. It must not affect:

- reference or interval ordering;
- pair-count allocation;
- task boundaries or task identifiers;
- random-number streams;
- read names, sequences, qualities, or mutation annotations;
- R1/R2 record order;
- BGZF block boundaries;
- compression level or compression input;
- final EOF blocks.

R1 and R2 are treated as one logical pair throughout generation. At every
record position, their base identifiers must match and differ only by the
expected `/1` and `/2` suffixes.

The new canonical files may differ from legacy DWGSIM once. After the new
single-thread baseline is accepted, every worker count must reproduce its
hashes exactly.

## Baseline and decision inputs

The full-reference benchmark on the development host used the complete
3,298,430,636-base GRCh38.p14 assembly, five million 2x150 bp pairs, eight
logical CPUs, reads-only BWA output, and BGZF level 1:

| Metric | Measured value |
| --- | ---: |
| Startup-adjusted pairs/s | 102,699.76 |
| Startup-adjusted bases/s | 30.81 million |
| Full-reference fixed work | 67.15 s |
| Estimated 100x runtime | 2 h 59 m 33 s |
| Estimated 100x BGZF output | 314.80 GiB |
| Peak RSS | 4,076.77 MiB |
| Average CPU including startup | 119% |

The low CPU utilization shows that the serial generator, rather than the
number of available BGZF helpers, is the present critical path.

A controlled one-million-pair chromosome 22 benchmark measured the current
zlib BGZF levels:

| Level | Eight-thread wall time | Size | Change from level 1 |
| ---: | ---: | ---: | ---: |
| 1 | 11.00 s | 307.1 MB | baseline |
| 4 | 11.75 s | 285.8 MB | +6.8% time, -6.9% size |
| 6 | 18.48 s | 276.8 MB | +68.0% time, -9.9% size |
| 9 | 24.94 s | 274.8 MB | +126.7% time, -10.5% size |

Level 4 is the provisional size-balanced default for the optimized design;
level 1 remains the fast profile. The final default must be confirmed after
generation and compression share all cores, because level 4 was 22% slower
than level 1 in a single-core-per-shard measurement.

## Architecture

```text
FASTA index + normalized regions + options + seed
                         |
                         v
               deterministic planner
                         |
       +-----------------+------------------+
       |                                    |
       v                                    v
parallel reference/mutation preparation   fixed pair tasks
       |                                    |
       +-----------------+------------------+
                         |
                         v
                generation workers
          {task id, R1 records, R2 records}
                         |
                         v
             task-local BGZF compression
          {task id, R1 chunks, R2 chunks}
                         |
                         v
                 bounded result ring
                   /             \
                  v               v
         ordered R1 appender   ordered R2 appender
                  |               |
                  v               v
             final R1 BGZF    final R2 BGZF
```

Execution order is unconstrained. Output order is determined solely by the
planner and task identifier.

## Deterministic planning

### Reference and region manifest

The planner reads the FASTA index and creates an immutable manifest in FASTA
contig order. Workers access reference data with indexed random reads,
`pread()`, or an equivalent thread-safe mechanism; they do not share a
mutable sequential `FILE *`.

BED input is normalized under the existing contract:

- zero-based, half-open intervals;
- reference contig order;
- nondecreasing coordinates;
- merged overlap and adjacency;
- exact contig-name matching;
- complete paired fragment containment.

Eligible bases and fragment-placement ranges are calculated before any read
task is dispatched.

### Pair-count allocation

For `-N`, allocation must sum exactly to the requested pair count. For
coverage-driven generation, the target pair count is calculated before
allocation. A stable largest-remainder allocation is used:

1. Calculate every contig or region's exact fractional share.
2. Assign each share's integer floor.
3. Distribute remaining pairs by descending fractional remainder.
4. Break equal remainders by canonical contig or region ordinal.

Allocation cannot depend on traversal order, worker count, or task completion.

### Pair tasks and canonical order

The provisional task size is 8,192 read pairs. Each task contains:

```text
format version
task id
contig ordinal
first pair ordinal
pair count
simulation option fingerprint
```

Task identifiers are assigned in canonical contig order and increasing pair
ordinal. The task size and ordering are format decisions: after golden hashes
are established, changing either requires a deliberate format-version change.

Large contigs naturally produce many tasks. Small contigs may share planning
batches, but their pair identities and canonical ordering remain independent
of how workers are grouped.

## Random-number model

The current global `drand48()` stream and static Box-Muller cache cannot be
used concurrently or partitioned reliably because the number of random draws
depends on rejection loops, errors, and indels.

The target uses a counter-based or independently splittable PRNG with
domain-separated keys. Random decisions are keyed by stable identities such
as:

```text
(root seed, format version, "mutation", contig, position/event)
(root seed, format version, "fragment", contig, pair ordinal)
(root seed, format version, "mate", contig, pair ordinal, mate)
(root seed, format version, "error", contig, pair ordinal, mate, cycle)
(root seed, format version, "quality", contig, pair ordinal, mate, cycle)
```

Retries consume only the relevant local domain. A rejected fragment placement
must not shift mutation, mate, error, or quality randomness for another pair.

Mutation preparation and read generation use separate streams. No mutable RNG
state is shared between workers. The same task therefore produces the same
paired records on every worker and under every scheduling order.

Any replacement for the existing normal-quality generator must first be
validated against the intended distribution. Performance work must not
silently change the accepted single-thread read sequences or qualities.

## Parallel reference and mutation preparation

Reference preparation is a dependency graph rather than one serial traversal:

1. Load indexed contigs independently.
2. Build or apply each contig's diploid mutation representation using its
   deterministic mutation stream.
3. Publish the immutable contig state.
4. Release all pair tasks for that contig.

The scheduler prepares contigs ahead of the generation window. Large contigs
are prioritized early to avoid a long tail. Immutable contig state is shared
by every pair worker using it and reference-counted until the last task
completes.

Reads-only `-r 0` runs should bypass unnecessary mutated-reference copies.
Mutation VCF/TXT output retains its own canonical contig/coordinate writer;
read-order relaxation does not permit unsorted variant output.

## Generation and paired records

A worker receives one fixed task and a read-only contig state. All mutable
scratch data is worker-local. It generates each read pair as one unit:

```c
struct paired_task_result {
    uint64_t task_id;
    uint32_t pair_count;
    buffer_t read1;
    buffer_t read2;
};
```

The two buffers must:

- contain the same number of records;
- use the same pair ordering;
- have matching normalized read identifiers;
- contain whole FASTQ records;
- contain no records belonging to another task.

FASTQ records should be assembled into contiguous buffers rather than emitted
one character at a time. This enables block writes and isolates generation
from filesystem latency.

## Deterministic BGZF chunks

Each task's R1 and R2 byte streams are compressed independently into BGZF
chunks using:

- a fixed compression implementation and version;
- an explicit compression level;
- deterministic block packing;
- no timestamps or scheduling metadata;
- no per-task BGZF EOF marker.

Task boundaries are compression boundaries. A final partial BGZF block is
allowed at the end of a task, making compressed task results independent and
appendable. The modest compression-ratio cost is bounded by choosing tasks
large enough to contain many full 64 KiB blocks.

Compression can run on the generation worker or a shared task pool. In either
case, compression input and block layout are fixed by task content, not by the
number of compression workers.

One canonical BGZF EOF block is appended to each final FASTQ after every task
has been committed.

## Ordered commit without a global generation lock

Completed compressed results enter a bounded ring indexed by task ID. The
appenders maintain `next_task_to_write` and consume only that slot.

R1 and R2 may use separate appender threads, but they consume the same task
sequence. A result is freed only after both appenders have consumed it.
Successful files therefore have matching pair order even if physical writes
finish at different times.

The scheduling window may scale with worker count because it does not alter
task identities or output order. Initial guidance is 4-8 outstanding tasks per
worker. Backpressure stops dispatch when the window or memory budget is full.
The scheduler prioritizes the lowest missing task to limit head-of-line stalls.

Files are written in a staging area. Both FASTQs are closed, synchronized, and
validated before publication. Because POSIX cannot atomically rename two
independent files, successful publication also writes one small completion
manifest atomically after both final names exist. Consumers that require
transactional behavior treat outputs without that manifest as incomplete.
Task failures occur before publication, and rename or crash recovery removes
or repairs an incomplete pair on the next invocation.

## Scaling beyond 100 cores

The task graph contains far more tasks than workers for a 100x human genome.
Scaling therefore depends on avoiding centralized compute and allocation:

- use worker-local scratch buffers and compressor state;
- use indexed result-ring slots rather than one global completed-task list;
- batch queue operations;
- keep reference and mutation state immutable;
- use NUMA-local worker queues and memory placement;
- steal locally before crossing NUMA nodes;
- perform only ordered appends in the writer threads;
- use large writes or `writev()` to minimize system-call overhead.

At the measured level-4 single-core rate of about 54,500 pairs/s, an ideal 100
cores would generate 1,099,476,879 pairs in approximately 3.4 minutes.
Practical compute time is expected to be 4-7 minutes after coordination and
NUMA overhead.

At that speed, output storage becomes the likely bottleneck. Approximately
293 GiB of level-4 output implies these lower bounds:

| Sustained output bandwidth | Output-only lower bound |
| ---: | ---: |
| 500 MiB/s | about 10 minutes |
| 1 GiB/s | about 5 minutes |
| 2 GiB/s | about 2.5 minutes |

Meaningful scaling beyond 100 cores therefore requires fast NVMe or a parallel
filesystem. Output is split naturally between R1 and R2 appenders, but their
logical task order remains identical.

## Memory budget

Reference and mutation representations are shared per contig; they are never
copied per read worker. The scheduler can retain the complete prepared genome
or evict completed contigs without affecting deterministic output.

A conservative 128-worker target is:

| Component | Target |
| --- | ---: |
| Reference and prepared diploid state | 55-100 GiB |
| Active uncompressed task buffers | less than 2 GiB |
| Compressed result/reorder window | 1-5 GiB |
| Compressor and worker scratch | less than 2 GiB |
| Index, BED, mutation, and task metadata | less than 10 GiB |
| Allocator, queues, and safety margin | 20-40 GiB |
| Expected peak | 90-150 GiB |
| Hard acceptance limit | less than 500 GiB |

Memory limits affect only cache residency and scheduling. They must never
change task boundaries, RNG keys, canonical order, or output hashes.

## Compression policy

Compression level becomes an explicit simulation option and a reported
benchmark field. Proposed profiles are:

- level 1: fastest output;
- level 4: provisional default size/runtime balance;
- levels 5-9: opt-in only.

The default must be chosen from an end-to-end parallel benchmark, not from
compression in isolation. Higher levels cannot solve the FASTQ size problem
alone: level 9 saved only about 10.5% versus level 1 while more than doubling
wall time in the controlled benchmark.

Quality binning or fixed qualities can reduce output much more, but they change
the simulated quality model and remain explicit scientific options rather than
compression defaults.

## Validation and acceptance gates

### Correctness

- Generate the same workload with `-t 1`, `-t 2`, `-t 8`, `-t 64`,
  `-t 128`, and automatic threads where hardware permits.
- Repeat each representative worker count.
- Compare SHA-256 hashes of complete R1 and R2 BGZF files.
- Validate every output with ordinary gzip and a BGZF-aware reader.
- Verify exact FASTQ record counts and 150-base sequence/quality lengths.
- Verify every normalized R1 identifier matches R2 at the same record index.
- Validate WGS, BED-restricted WGS, and WES-like BED workloads.
- Validate mutation-free, generated-mutation, and supplied-mutation modes.
- Inject reversed scheduling, random task delays, and one deliberately slow
  lowest-ID task; hashes must remain unchanged.
- Confirm only one canonical BGZF EOF block is required and accepted.
- Confirm a failed task never publishes one completed mate file without the
  other.
- Interrupt each publication step and confirm that the completion manifest
  never identifies a missing, partial, or mismatched mate file as complete.

### Performance

- Report pair/s, read/s, base/s, CPU, peak RSS, output bytes, compression
  level, and storage throughput.
- On the existing eight-core development host, target at least 350,000
  default-quality pairs/s for paired BWA output at the selected default
  compression level.
- Measure compute scaling separately from filesystem scaling with a bounded
  sink benchmark.
- Demonstrate useful scaling at 64 and 128 workers without unbounded queue
  growth.
- Keep expected 128-worker RSS below 150 GiB and enforce the 500 GiB hard
  limit.

### Determinism scope

The test matrix records the compiler, architecture, libc/libm, compressor
implementation and version, compression level, format version, and reference
checksums. Cross-host hash guarantees are added only after these dependencies
are pinned or replaced with platform-independent implementations.

## Decisions to close before freezing golden hashes

- Select and version the counter-based PRNG and domain-key encoding.
- Confirm the provisional 8,192-pair task size with compression-ratio,
  straggler, and reorder-memory measurements.
- Select and pin the BGZF compressor implementation.
- Confirm level 4 or choose another default from the completed parallel
  pipeline benchmark.
- Choose the staging-directory and completion-manifest naming contract.
- Decide whether cross-host hashes are part of the first release contract.

Changing any choice that affects pair identities, canonical order, FASTQ
formatting, task compression boundaries, or compressor bytes requires a new
format version and a deliberately regenerated golden baseline.

## Implementation stages

1. **Planner and baseline:** implement stable pair allocation, task IDs, and
   domain-separated RNG in single-thread mode; accept new golden hashes.
2. **Paired task workers:** generate fixed R1/R2 task buffers in parallel while
   retaining the current ordered compressor.
3. **Task-local BGZF:** compress deterministic task chunks independently and
   append them through the result ring.
4. **Parallel reference preparation:** replace sequential FASTA traversal with
   indexed contig tasks and immutable shared mutation state.
5. **Memory and NUMA controls:** add bounded dispatch, contig eviction, local
   queues, and cross-node work stealing.
6. **Compression selection:** expose the level, confirm the default, and test
   hashes for every supported level.
7. **Scale validation:** run the complete determinism, memory, and 64/128-core
   performance gates.

Each stage must preserve the accepted single-thread paired-read contents and
must pass hash determinism before the next stage begins.
