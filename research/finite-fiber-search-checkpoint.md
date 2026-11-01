# Completed finite-fiber basis and iteration searches

## Result and coverage

Two thirty-minute workers completed their declared budgets on the
primitive finite-fiber crossings. Neither produced a factoring-threshold
candidate. Their full records are retained as compressed archives,
and every retained winner has been freshly reconstructed and reprofiled.

| Worker | Recorded coverage | Elapsed seconds |
|---|---|---:|
| Constant physical bases and full gauges | 877 proposals, 30 duplicate bases, 847 retained trial sequences | 1,800.106 |
| F4 Smith-fiber iteration after deflation | 27,784 proposals, 3,455 distinct solved target vectors, 6,910 complete profile calls | 1,800.068 |

The physical-basis worker completed 836 trial sequences; eleven reached
a local or global time limit. A completed sequence is not a proof of
optimizer convergence or global optimality. The archive retains each
optimizer outcome separately.

The basis worker recorded 811,415 objective calls, including 262,937
unresolved native calls that received an unfavorable search penalty.
Those calls supplied no partial-profile bound and do not exclude the
corresponding gauges. All four retained winners resolve their complete
profiles in the fresh replay.

The iteration worker recorded 6,059 distinct profile-class records.
It kept forty complete polynomial witnesses. During the search it ran
236 full frame replays, including witnesses subsequently replaced in
the retained list. The new checkpoint independently replays the forty
retained witnesses; it does not claim a fresh coefficient reconstruction
of every discarded sample.

## Retained diagnostic values

These are floating complete-profile pressure evaluations. The two
source indices refer to
[`finite-fiber-basis-input.json`](../evidence/finite-fiber-basis-input.json)
and its
[exactly deflated counterpart](../evidence/deflated-finite-fiber-input.json).

| Source | Best basis balanced value | Best positive supporting-line value | Best recorded deflated iteration value |
|---:|---:|---:|---:|
| 0 | 0.7381234955 | 1.0020967840 | 0.7381196578 |
| 1 | 0.7356448409 | 1.0021226035 | 0.7374853070 |

The balanced and supporting-line thresholds are `1/sqrt(2)` and `1`,
respectively. None of these values passes its threshold. Finite
optimizer outputs are not global lower bounds for either family.

## Exact replay and archived sources

The [checkpoint tool](../tools/checkpoint_finite_fiber_searches.py)
checks both completed budgets, all frozen source and input hashes,
proposal and duplicate accounting, and profile-class counts. It then
reconstructs every winning constant physical basis and all forty
retained full-entry fiber witnesses. Exact scalar division is repeated
before profiling a primitive witness; old frame bounds are not assigned
to the divided tensor without re-evaluation.

Fresh controls include 21,504 strict endpoint entries, 1,408 complete
profile channels, forty full frame-witness checks, and four reconstructed
basis winners. Every saved tensor hash agrees. The 1,280 integral fiber
channels agree exactly; the 128 basis-winner channels agree at their
saved floating heights. Field arithmetic and endpoint checks are exact.
Real-weight chamber selection and final transcendental pressures remain
numerical evaluations.

The complete byte-preserving archives are
[physical bases](../evidence/finite-fiber-physical-basis-search.json.gz)
and [deflated iteration](../evidence/deflated-f4-smith-fiber-search.json.gz).
Their 3,789,662 and 10,064,914 original bytes compress to 758,686 and
820,818 bytes. The
[compact checkpoint](../evidence/finite-fiber-search-checkpoint.json)
binds both compressed and uncompressed hashes and stores all fresh
winner profiles.

```sh
OPENBLAS_NUM_THREADS=1 python tools/checkpoint_finite_fiber_searches.py
```

The command reads the committed gzip archives directly. The frozen
search sources are
[`search_explicit_physical_bases.py`](../tools/search_explicit_physical_bases.py)
and [`search_iterated_f4_smith_fibers.py`](../tools/search_iterated_f4_smith_fibers.py).
The construction and profile guarantees are in
[f4-iterated-smith-fibers.md](f4-iterated-smith-fibers.md) and
[finite-jet-smith-fibers.md](finite-jet-smith-fibers.md).
