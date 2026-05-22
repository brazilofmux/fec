# Reed-Solomon Decoder Benchmarks

This file captures `rstest -b` results across hardware platforms so future
sessions on other boxes can sanity-check their numbers against a known baseline.

## How to reproduce

```
cd rs/librs && make
cd ../rstest && make && ./rstest -b
```

`rstest -b` iterates `tt ∈ {1..8, 16, 24, 32, 40, 48, 56, 64}` and, for each
`tt`, runs 100k iterations of each decoder variant on the same codeword.
Before timing it also runs an internal cross-decoder correctness probe — if
Auto / General (LFSR) / scalar Direct / Direct NEON / Direct AVX2 / Direct AVX-512
disagree on the corrected output, the run aborts with
`ERROR: decoder disagreement at tt=N`.

The `General (LFSR)` column is the historical baseline: the LFSR/Horner-form
decoder (`RS_DECODER_GENERAL`) that Auto routed to for any unspecialized `tt`
before the Direct family landed. Same-run comparison against it is the speedup
vs what production used to ship.

All times are µs/decode (lower is better). Two columns per decoder:
`w/err` = codeword with `tt/2` injected errors, `clean` = codeword with no errors.
The `[Auto phases]` line breaks down where time is spent in the Auto decoder:
`synd` (syndrome calculation — the SIMD-vectorized part), `BM` (Berlekamp-Massey),
`Chien` (Chien search), `Forney` (Forney's algorithm).

## Recording a new run

Append a new section below with:
- Date and commit hash (`git rev-parse --short HEAD`)
- CPU model and the SIMD flags actually used by the build
- Compiler and version
- Whatever cpufreq/thermal context is relevant (a VM has none — say so)
- The raw `rstest -b` output, fenced

Format the table as Markdown if you want it pretty, but the fenced raw output
is the source of truth — keep that intact for future grep-based comparisons.

---

## x86_64 — Cascade Lake (Xeon 8259CL) — 2026-05-17

- **Commit:** `1d974a8` (Fix `_mm512_broadcast_i128` → `_mm512_broadcast_i32x4`)
- **CPU:** Intel Xeon Platinum 8259CL @ 2.50 GHz (turbo to ~3.10 GHz observed)
- **SIMD features:** AVX2, AVX-512F, AVX-512BW, AVX-512DQ, AVX-512CD, AVX-512VL
- **Compiler:** g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
- **Build flags:** `-O3 -fomit-frame-pointer` (+ `-mavx2` on the AVX2 TU,
  `-mavx512f -mavx512bw` on the AVX-512 TU)
- **Host:** Linux VM, no cpufreq governor file (cloud guest); single-thread run

### Headline

| tt  | scalar Direct (w/err) | NEON (w/err, scalar fallback) | AVX2 (w/err) | AVX-512 (w/err) | AVX2 vs scalar |
|----:|----------------------:|------------------------------:|-------------:|----------------:|---------------:|
|   1 |                  2.49 |                          2.27 |         0.90 |            0.98 |          2.8× |
|   8 |                 12.52 |                          7.55 |         1.18 |            1.36 |         10.6× |
|  16 |                 25.47 |                         14.57 |         1.79 |            2.25 |         14.2× |
|  32 |                 52.16 |                         28.86 |         4.02 |            4.36 |         13.0× |
|  64 |                104.24 |                         62.44 |        11.40 |           12.17 |          9.1× |

AVX2 beats AVX-512 by ~5–15 % on this CPU at every `tt`. Cascade Lake clocks
down on AVX-512 workloads and the per-iteration `vbroadcasti32x4` adds cost the
AVX2 path doesn't pay; net result is AVX2 wins on Skylake-X / Cascade Lake even
though AVX-512 has 2× the lanes. Expect AVX-512 to overtake on Ice Lake and
later (lighter license penalty + faster `vbroadcasti32x4`).

The NEON decoder runs its portable scalar fallback on x86 (the NEON intrinsics
compile out), so its "NEON" column here is roughly 2× faster than scalar Direct
because of incidental codegen differences, not vector ISA.

### Raw output

```
RS Test Configuration:
  Mode: Benchmark
  Number of codewords: 10000
  Verbosity level: Normal
  Random seed: 1093
  Verify correction: Yes
//                Auto             Direct (scalar)    Direct (NEON)        Direct (AVX2)       Direct (AVX512)
// tt  Encode  w/err   clean      w/err   clean      w/err   clean      w/err   clean      w/err   clean
//  1    0.81     1.10    1.09       2.49    2.59       2.27    2.29       0.90    0.85       0.98    0.98
//     [Auto phases @ tt=1]  synd=1.07  BM=0.19  Chien=0.33  Forney=0.10  (total=2.02)  errs=0
//  2    0.69     1.02    1.08       3.94    3.92       2.85    2.87       0.87    0.86       1.07    1.02
//     [Auto phases @ tt=2]  synd=2.77  BM=0.81  Chien=2.32  Forney=0.58  (total=6.88)  errs=1
//  3    1.51     1.12    1.11       5.26    5.60       3.53    3.53       0.89    0.89       1.09    1.09
//     [Auto phases @ tt=3]  synd=1.91  BM=0.80  Chien=2.03  Forney=0.73  (total=5.80)  errs=1
//  4    1.86     1.11    1.16       6.66    6.87       4.30    4.31       0.99    0.94       1.17    1.14
//     [Auto phases @ tt=4]  synd=1.44  BM=1.05  Chien=2.06  Forney=0.59  (total=5.47)  errs=2
//  5    1.59     1.19    1.18       8.26    8.07       5.26    5.32       0.97    0.97       1.25    1.19
//     [Auto phases @ tt=5]  synd=1.78  BM=1.24  Chien=2.67  Forney=0.69  (total=6.71)  errs=2
//  6    1.99     1.23    1.32       9.50    9.64       6.00    5.76       1.01    1.06       1.21    1.27
//     [Auto phases @ tt=6]  synd=2.08  BM=1.61  Chien=2.86  Forney=0.87  (total=7.73)  errs=3
//  7    1.52     1.32    1.36      13.61   13.88       6.79    6.70       1.07    1.10       1.29    1.27
//     [Auto phases @ tt=7]  synd=1.75  BM=2.04  Chien=2.68  Forney=0.99  (total=7.77)  errs=3
//  8    1.78     1.42    1.54      12.52   12.74       7.55    7.63       1.18    1.21       1.36    1.43
//     [Auto phases @ tt=8]  synd=2.23  BM=2.11  Chien=3.40  Forney=1.17  (total=9.22)  errs=4
// 16    1.70     2.05    2.11      25.47   24.89      14.57   14.44       1.79    1.82       2.25    2.16
//     [Auto phases @ tt=16]  synd=2.41  BM=4.61  Chien=5.88  Forney=2.28  (total=15.50)  errs=8
// 24    1.73     3.04    3.07      39.00   39.74      21.40   21.25       3.13    3.09       3.20    3.07
//     [Auto phases @ tt=24]  synd=2.36  BM=10.97  Chien=25.82  Forney=6.54  (total=46.12)  errs=12
// 32    1.45     4.28    4.38      52.16   51.93      28.86   29.17       4.02    4.12       4.36    4.31
//     [Auto phases @ tt=32]  synd=1.73  BM=11.94  Chien=11.28  Forney=6.43  (total=31.72)  errs=16
// 40    1.39     6.31    6.46      65.25   64.26      36.17   36.52       5.72    5.90       6.35    6.40
//     [Auto phases @ tt=40]  synd=3.10  BM=17.34  Chien=14.06  Forney=10.00  (total=44.86)  errs=20
// 48    1.28     8.22    8.04      77.34   76.50      43.74   44.12       7.39    7.20       8.65    8.45
//     [Auto phases @ tt=48]  synd=3.54  BM=23.28  Chien=16.51  Forney=12.81  (total=56.47)  errs=24
// 56    1.25    10.17    9.97      89.31   90.71      55.79   54.70       9.17    9.71      10.23   10.48
//     [Auto phases @ tt=56]  synd=3.50  BM=29.30  Chien=19.27  Forney=16.66  (total=69.07)  errs=28
// 64    1.22    12.07   11.99     104.24  103.39      62.44   61.98      11.40   11.17      12.17   12.12
//     [Auto phases @ tt=64]  synd=4.61  BM=53.64  Chien=23.01  Forney=22.80  (total=104.42)  errs=32
```

### What to look for in a re-run on the same family of hardware

- **Cross-decoder correctness probe:** absence of any `ERROR: decoder disagreement`
  line is the primary signal. A regression here is a correctness bug.
- **Headline speedups:** at `tt=64`, AVX2 should land in the 10–15 µs band on
  any post-Skylake-X x86 part. A 50+ µs result means SIMD isn't actually being
  used — check `__builtin_cpu_supports("avx2")` / `avx512bw`, the makefile flags,
  and that the runtime check in `RS_DECODER_DIRECT_AVX2::cpu_has_avx2()` returns
  true.
- **Phase breakdown:** `synd` time should be small relative to BM+Chien+Forney
  at high `tt`. If `synd` blows up disproportionately, the SIMD path is broken
  or the scalar fallback is being taken.
- **Auto vs forced AVX-512:** Auto should track AVX-512 closely (≤ 0.2 µs
  overhead from the codec/vtable indirection). A large gap means
  `create_best_decoder` is selecting the wrong path.

---

## aarch64 — Apple M5 Max (MacBook Pro) — 2026-05-22

- **Commit:** `b32ff9d` (Add General (LFSR) column to rstest -b benchmark)
- **CPU:** Apple M5 Max, 18 cores (macOS reports `Mac17,6`)
- **SIMD features:** NEON (always-on on aarch64; no compile flag needed)
- **Compiler:** Apple clang 21.0.0 (clang-2100.1.1.101), target `arm64-apple-darwin25.4.0`
- **Build flags:** `-O3 -fomit-frame-pointer` (librs) / `-O3 -finline-functions -Wall` (rstest);
  no `-mavx2` / `-mavx512*` on aarch64 (Makefile gates those on `uname -m == x86_64`)
- **Host:** macOS 26.4.1 on bare metal, default scheduler, single-thread run

### Headline

Auto routes to `RS_DECODER_T1` at `tt=1` (specialized) and to `RS_DECODER_DIRECT_NEON`
for every other `tt`, so the Auto and NEON columns coincide for `tt ≥ 2`. The
AVX2 / AVX-512 columns are the portable scalar fallback on this host (the
intrinsics compile out), included for cross-decoder correctness only.

Same-machine speedup of Auto over the historical General (LFSR) baseline:

| tt  | General (LFSR, w/err) | Auto (w/err) | Direct NEON (w/err) | Auto vs General |
|----:|----------------------:|-------------:|--------------------:|----------------:|
|   1 |                  0.90 |         0.47 |                0.22 |           1.9×  |
|   8 |                  1.43 |         0.33 |                0.33 |           4.3×  |
|  16 |                  2.58 |         0.64 |                0.64 |           4.0×  |
|  32 |                  5.17 |         1.40 |                1.39 |           3.7×  |
|  64 |                 12.45 |         3.91 |                3.90 |           3.2×  |

At `tt=64` the Auto-phase line shows `synd=1.67` out of `total=34.92` — syndromes
are <5 % of total decode time once NEON has had its way with them. Further wins
at high `tt` would have to come from BM / Chien / Forney, not syndrome calc.

The encode column bounces (`tt=1: 0.42`, `tt=3: 1.10`, `tt=4: 0.56`, `tt=5: 1.07`)
because only `tt ∈ {1, 2, 4, 8, 14, 16, 32, 64}` have hand-written
`RS_FLIPPED_ENCODER_T*` specializations; everything else falls through to
`RS_FLIPPED_ENCODER_GENERAL` at roughly 2× the cost. No SIMD on the encode side.

### Raw output

```
RS Test Configuration:
  Mode: Benchmark
  Number of codewords: 10000
  Verbosity level: Normal
  Random seed: 1093
  Verify correction: Yes
//                Auto           General (LFSR)     Direct (scalar)    Direct (NEON)        Direct (AVX2)       Direct (AVX512)
// tt  Encode  w/err   clean      w/err   clean      w/err   clean      w/err   clean      w/err   clean      w/err   clean
//  1    0.42     0.47    0.48       0.90    0.90       0.62    0.63       0.22    0.22       0.63    0.62       0.63    0.62
//     [Auto phases @ tt=1]  synd=0.67  BM=0.04  Chien=0.04  Forney=0.04  (total=0.88)  errs=0
//  2    0.45     0.23    0.23       0.97    0.97       1.04    1.04       0.23    0.23       0.82    0.78       0.86    0.84
//     [Auto phases @ tt=2]  synd=0.58  BM=0.29  Chien=0.88  Forney=0.21  (total=2.08)  errs=1
//  3    1.10     0.25    0.26       1.11    1.11       1.41    1.40       0.26    0.26       1.00    1.01       1.05    1.05
//     [Auto phases @ tt=3]  synd=0.58  BM=0.29  Chien=0.88  Forney=0.21  (total=2.00)  errs=1
//  4    0.56     0.25    0.25       1.22    1.23       1.78    1.77       0.25    0.25       1.20    1.21       1.22    1.31
//     [Auto phases @ tt=4]  synd=0.71  BM=0.46  Chien=1.00  Forney=0.17  (total=2.42)  errs=2
//  5    1.07     0.27    0.26       1.33    1.33       2.19    2.19       0.26    0.27       1.42    1.43       1.44    1.43
//     [Auto phases @ tt=5]  synd=0.83  BM=0.50  Chien=0.96  Forney=0.21  (total=2.50)  errs=2
//  6    1.07     0.29    0.29       1.41    1.42       2.62    2.62       0.29    0.29       1.66    1.66       1.68    1.68
//     [Auto phases @ tt=6]  synd=0.67  BM=0.33  Chien=0.96  Forney=0.21  (total=2.29)  errs=3
//  7    1.06     0.32    0.32       1.43    1.44       3.08    3.08       0.32    0.32       1.93    1.94       1.94    1.95
//     [Auto phases @ tt=7]  synd=0.58  BM=0.88  Chien=0.96  Forney=0.38  (total=2.83)  errs=3
//  8    0.61     0.33    0.34       1.43    1.45       3.44    3.46       0.33    0.33       2.18    2.19       2.18    2.18
//     [Auto phases @ tt=8]  synd=0.58  BM=1.00  Chien=1.12  Forney=0.25  (total=3.08)  errs=4
// 16    0.63     0.64    0.65       2.58    2.59       7.64    7.62       0.64    0.65       4.26    4.28       4.28    4.26
//     [Auto phases @ tt=16]  synd=0.96  BM=2.21  Chien=1.83  Forney=0.75  (total=5.88)  errs=8
// 24    0.98     1.04    1.05       3.93    3.99      10.41   10.54       1.03    1.05       6.55    6.46       6.54    6.56
//     [Auto phases @ tt=24]  synd=0.92  BM=4.29  Chien=3.00  Forney=1.21  (total=9.46)  errs=12
// 32    0.57     1.40    1.41       5.17    5.18      13.74   13.74       1.39    1.38       8.61    8.62       8.70    8.71
//     [Auto phases @ tt=32]  synd=1.08  BM=6.42  Chien=3.46  Forney=1.83  (total=12.88)  errs=16
// 40    0.84     1.92    1.91       6.62    6.63      17.04   17.10       1.90    1.90      10.92   11.03      10.92   10.93
//     [Auto phases @ tt=40]  synd=1.17  BM=9.12  Chien=3.92  Forney=2.83  (total=17.12)  errs=20
// 48    0.77     2.50    2.49       8.98    8.99      20.25   20.26       2.50    2.49      13.24   13.14      13.48   13.47
//     [Auto phases @ tt=48]  synd=1.33  BM=11.71  Chien=5.21  Forney=3.50  (total=21.88)  errs=24
// 56    0.70     3.12    3.12       9.72    9.72      23.78   23.89       3.14    3.13      15.72   15.79      15.83   16.33
//     [Auto phases @ tt=56]  synd=2.62  BM=16.21  Chien=5.50  Forney=4.42  (total=28.79)  errs=28
// 64    0.52     3.91    3.93      12.45   12.47      27.34   27.54       3.90    3.90      18.21   18.24      18.31   18.46
//     [Auto phases @ tt=64]  synd=1.67  BM=20.29  Chien=6.79  Forney=6.04  (total=34.92)  errs=32
```

### Note on the older x86_64 entry above

The Cascade Lake entry from 2026-05-17 predates the `General (LFSR)` column.
Its raw output has 10 columns instead of 12; that's a known format gap, not a
data problem. When that box (or another x86 host) re-runs, the new layout will
overwrite this asymmetry naturally.
