# Runtime profile and allocation pass (2026-07-21)

## Scope and reproducibility

The baseline is dev commit `9a036db`. Profiling and timing used LongBench H69
PacBio name-collated BAMs, four threads, `--coverage-model auto`,
`--seq-tech pac-bio`, `--endpoint-weight 0.3`, and `--em-accel none`. Disabling
acceleration deliberately exposes the cost of repeated fixed-point evaluations.
Raw profiler data is retained outside the repository at
`oarfish-evaluation-data/profiling-20260721/h69-250k.data`; repeated benchmark
tables are in `runtime-opt-baseline-9a036db-20260721` and
`runtime-opt-pass1-20260721` under the same evaluation-data directory.

## Baseline profile

On H69 250k, the median/mean run is about 5.0 seconds. Instrumented phase means
from three baseline runs were:

| Phase | Seconds | Share of 5.04 s wall time |
|---|---:|---:|
| Alignment parsing/filtering | 0.794 | 15.8% |
| Coverage modeling, including anchor inference | 2.007 | 39.8% |
| Final EM | 1.597 | 31.7% |
| Digest, output, and other overhead | 0.642 | 12.7% |

The sampling profile attributes roughly 53% of CPU samples to the inlined EM
fixed-point work. Other visible costs include zlib inflation (8.4%), Rayon
scheduling/reduction (7.1%), memory clearing (4.2%), memory copies (2.5%), BAM
decoding (1.4%), and the physical-endpoint likelihood path (about 1.4% including
children). Thus the two inference passes, not endpoint fitting itself, are the
main optimization target.

## Retained changes

1. Parallel EM now writes each fixed-point result directly into the driver's
   destination buffer. The old path accumulated into a second transcript-sized
   vector and copied it after every evaluation (2,000 copies for the two
   unaccelerated PacBio auto inference passes in this benchmark).
2. Hybrid/adaptive coverage combination reuses one scratch vector across reads
   rather than allocating and freeing a vector for every read. PacBio auto calls
   this path twice, so H69 250k avoids about 400,000 short-lived allocations.
3. Physical-endpoint likelihoods are appended and normalized directly in their
   final packed vector rather than first allocating a temporary vector per read.
4. Logistic coverage normalization reserves its known final packed capacity.

The optimized and baseline H69 250k quantification files are byte-identical.
All 34 release tests pass.

Five-run `perf stat` measurements show the retained pass reduces aggregate CPU
task time from 14.110 to 13.482 seconds (-4.5%), cycles from 67.72 to 65.46
billion (-3.3%), and cache misses from 824.3 to 806.3 million (-2.2%). Wall time
was scheduling-sensitive on the shared machine (5.30 +/- 0.09 seconds baseline
versus 5.46 +/- 0.15 optimized in that ordered run), so no wall-clock claim is
made from that noisy comparison. The smaller H69 50k three-run phase means did
move in the expected direction: coverage 0.687 to 0.670 seconds (-2.5%) and EM
0.565 to 0.555 seconds (-1.7%), with essentially unchanged total wall time and
RSS.

## Rejected experiment

The raw-read producer and mapper batch vectors currently retain thread-local
capacity but clone batches when sending them through channels. Moving the
buffers instead removed those copies and reduced one 50k SIRV raw-read run from
7.62 to 7.47 seconds, but peak RSS increased from 106 to 120 MiB as replacement
buffers churned through the allocator. Quantification differed only by
`4.55e-13` reads (parallel floating-point order). This tradeoff was reverted.
A future version should use an explicit buffer-return pool if this path is
revisited, preserving both ownership transfer and thread-local capacity.

## Next targets

The dominant remaining memory traffic is clearing every dense per-shard
transcript accumulator and reducing all shards across all transcripts on every
EM evaluation. A sparse touched-transcript scheme is promising only when the
active set is small; it needs workload-sensitive selection because dense human
transcriptomes can make bookkeeping slower than `memset`. BAM decompression is
the next independent phase worth testing with parallel block decompression or a
faster input representation.
