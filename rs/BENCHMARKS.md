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
Auto / scalar Direct / Direct NEON / Direct AVX2 / Direct AVX-512 disagree on
the corrected output, the run aborts with `ERROR: decoder disagreement at tt=N`.

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

## aarch64 — _placeholder_

When the M-series / Graviton box re-runs, add a section here mirroring the
format above. The dedup commit (`721a55e`) reported partial data
(`synd` NEON: 1.34 µs vs scalar 13.49 µs at tt=32; 3.95 µs vs 27.4 µs at tt=64)
but not a full headline / phase table, so leaving this blank rather than
inventing rows.
