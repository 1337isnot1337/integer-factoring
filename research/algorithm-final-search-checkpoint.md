# Completed entrywise, finite-frame, mixed-place, and hot-code searches

This checkpoint preserves the completed searches that accompany the
constructive notes. It does not assert a global optimum, an exclusion of
unsearched tensors, or a new factoring exponent. The strongest exact q96
bound is stated separately in [the rational interval certificate](exact-natural-pressure-intervals.md),
with [independent integer-power controls](hot96-rational-pressure.md).

The portable archive is
[`algorithm-final-searches.json.gz`](../evidence/algorithm-final-searches.json.gz).
It contains the original JSON bytes of each report and each external input.
Historical absolute paths remain provenance strings; the replay reads their
contents from the archive. The compact result is
[`algorithm-final-searches-checkpoint.json`](../evidence/algorithm-final-searches-checkpoint.json).
The [checkpoint tool](../tools/checkpoint_algorithm_final_searches.py) performs
the exact and numerical replays described below.

## 1. Exact entrywise and finite-frame targets

The [entrywise interpolation search](entrywise-reflection-jet-interpolation.md)
completed its entire list of 5,760 targets. There were 1,542 feasible targets.
The archive retains 111 distinct complete tensors selected by degree, jet
order, protection, and numerical cost. The checker solves every retained
target again from the original reflection tensor and requires the entire
coefficient tensor to agree with its saved hash.

Three [regular-frame finite-jet searches](finite-jet-smith-fibers.md) each
completed all 234 declared targets. With infinity relaxation windows 0, 4,
and 6, respectively, they found 0, 72, and 90 feasible targets. These are
finite template counts, not universal infeasibility statements.

The last two searches retain 80 and 88 primitive tensors. For every retained
tensor, the checker restores its saved scalar factor with the endpoint
normalization, recovers the full original tensor hash, verifies every
prescribed Hasse jet and regular-frame row target, and recomputes the two
complete integer profiles. It then repeats the full scalar gcd division
and the complete profile at the finite place. This contributes 8,064
integer channels across 168 independent retained witnesses. The primitive
tensor is never assigned the pre-division degree or frame certificate.

## 2. Full-cost mixed-place refinements

The archive includes complete seed inputs and disjoint-place specifications
for the dense affine, quadratic-place, Bell-jet, finite-frame, truncated,
and stabilized hot16 refinements. Balanced and asymmetric costs are recorded
separately. Asymmetric support must be below 1 to pass its transfer gate;
the balanced threshold is \(1/\sqrt2\).

For each archived explicit mixed-place run, the checker reconstructs its
best recorded seed, repeats scalar normalization, validates the chosen
partition, and evaluates its saved optimizer parameters through the full
local-profile interface. It requires every complete profile to have the
correct number of channels and binds each saved conjugate label to the
validated place group. The common physical types, degree, and allocation
costs remain included. These are floating-point cost replays, not exact
inequality certificates.

The main completed lists contain 1,200 dense-affine starts, 108 quadratic
starts, 175 balanced and 175 asymmetric Bell-jet starts, 190 balanced and
190 asymmetric finite-frame starts, and 240 truncated-q4 starts. A separate
900-second affine asymmetric run attempted 270 starts and completed 269
optimizer sequences. The stabilized hot16 pilot attempted all 12 declared
starts; each reached its 90-second local cap, so **zero full optimizer
sequences completed**. Its best retained value remains the natural control
\(\rho\approx0.7147366313\).

The earlier 600-second affine-word run and its complete finite-rank input
are retained as provenance: 1,001 starts were attempted and 1,000 optimizer
sequences completed. This archive does not reconstruct every affine word
or repeat its optimization. The later explicit-seed replays above and the
separate exact affine witnesses have their own stated checks.

The coherent-boundary refinement supplied by the independent audit lane is
also preserved. It contains six declared trials, with three completed
optimizer sequences and three capped sequences. For this supplied report,
the checkpoint validates the source/input hashes and accounting only; it
does not repeat optimization or its winner calculations. Its exact
construction and finite q8 certificate are in
[the coherent-boundary note](coherent-hot-boundary-crossings.md).

## 3. Adaptive catalogues and interrupted-worker accounting

The [adaptive fixed-type worker](../tools/search_adaptive_fixed_type_gauges.py)
completed all 68 declared q4 starts. The q16
pilot used its 900-second budget: 39 starts were attempted and 18 optimizer
sequences completed. The checker independently reconstructs all 111 saved
local catalogue profiles and their type-mean shifts. Their sampled curve
maxima remain numerical diagnostics; they do not replace an exact finite
block certificate. The [finite typed wrapper](finite-typed-adaptive-certificates.md)
states the sufficient exact gate.

The q16 coordinate-coupled and square-zero searches were interrupted during
the session. Their original prefixes are preserved unchanged, and their
continuations use declared new random seeds. The checker binds each
continuation to the exact prefix hash and checks that its saved elapsed
time covers the remaining declared budget. The recorded combined times
are approximately 3,606.998 and 3,601.457 seconds, respectively. They are
separate prefixes and continuations, not restarts falsely counted as
additional completed original work. The archive retains all timeout and
unresolved-query counters. The q16 truncation worker completed its declared
60-screen and refinement lists before using the entire wall-clock cap.

## 4. Hot-code follow-up

The final hot-code records include the full
q32, q48, q64, and q96 charge profiles; the 30-minute
q24 coordinate-downset walk; and the 20-minute weighted compound-cone
search. The compact result records the actual completed sizes and trial
counts. No partially computed large charge block is retained as a profile.

The completed energy-prefix profiles give these numerical diagnostics:

| Alphabet size | Complete channels | Natural \(\rho\) |
| --- | ---: | ---: |
| 32 | 1,024 | 0.7121199592 |
| 48 | 2,304 | 0.7111215197 |
| 64 | 4,096 | 0.7105903584 |
| 96 | 9,216 | 0.7100058749 |

The [q96 rational certificate](hot96-rational-pressure.md) proves
\(\rho<49/69\) independently of the displayed floating evaluation. The
checkpoint binds its full profile to the archived constructor output and
repeats both stored strict integer comparisons.

The downset winner is reconstructed through the exact charge-block engine.
Every completed recorded code is checked for deletion closure and the full
channel count. The weighted winner's rational parameters satisfy the
integrality cone by an integer comparison, and its complete profile and
common Gibbs-type cost are recomputed. Integer cone membership certifies
admissibility; its numerical pressure does not itself certify a strict
factoring improvement.

The downset walk completed 2,521 distinct proposals in 1,800.222 seconds;
every complete-profile call resolved. The weighted search completed 3,041
full-profile queries and 360 search iterations in 1,200.110 seconds, with
zero unresolved calls and zero failed linear programs. Both retained the
initial q24 value \(\rho\approx0.7129641788\). These are completed bounded
searches, not proofs that other codes or admissible gauges cannot improve it.

Large retained charge profiles receive fresh pressure and dimension checks
in this archive replay. The expensive native elimination is not repeated
by this checkpoint. The q64 profile is additionally required to equal the
profile in its separately replayed rational certificate. This distinction
keeps the checks reproducible without labeling a saved-profile comparison
as an independent large-matrix reconstruction.

## 5. Reproduction and scope

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python tools/checkpoint_algorithm_final_searches.py
```

The default command reads only committed files. `--build` is an archival
operation for the original session paths and is not required for replay.
`--include-hot` adds the final completed hot workers when building an archive.
The archive contains the exact source hashes used by its original workers;
the checker requires those sources to agree before replay.

The result separates exact reconstruction counts, integer profile channels,
floating profile channels, and optimization accounting. A finite search
can miss a better chamber, physical basis, code, type, or place allocation.
The retained tensors and complete profiles remain constructive inputs for
such follow-up work.
