# Weighted closed-place refinement of paired restrictions

## Result

The complete catalogue of 200
[paired physical restrictions](finite-fiber-paired-restrictions.md)
was screened across 1,978 distinct nonempty assignments of its ordinary
singular conjugate places. A subsequent full weighted refinement selected
24 seeds and completed all 402 planned local optimization trials. It
used the [unequal paired-alphabet pressure](paired-alphabet-mixed-place-crossings.md)
and the [closed-place oracle](closed-place-weighted-oracles.md), including
irreducible factors of degrees up to seven.

The best recorded bound among these selected partitions with additional
places was `rho = 0.8515910806692913`. It does not improve the best
previously recorded restriction, approximately `0.7416677`, or the
reflection construction. No improved factoring exponent is claimed.

Twelve selected seeds did improve their own earlier recorded two-end
values. For example, a two-by-two restriction with several extension-field
places moved from approximately `0.996408` to `0.881191`. These are useful
checks that the new places can contribute. They do not establish global
optimality of either optimizer.

## 1. Scope of the completed search

The unweighted stage uses the exact complete determinant catalogues as
input. Conjugates of the same irreducible polynomial have identical
unweighted profiles. Therefore it enumerates how many conjugates are
assigned to the row cut, removing only permutations of identical
conjugates. There are 1,978 valid partitions with both sides nonempty.
Constant degree-zero controls are handled directly. The numerical
capacity grids select candidates; they are not outward interval proofs.

The weighted stage alternates seeds with promising unweighted adaptive
values and seeds with promising earlier weighted values. It retains 24
distinct tensors with additional ordinary singular places. For each it
tries a catalogue-based partition, selected single-factor extensions of
the two natural ends, and an opposite-factor partition when available.
The resulting 134 seed-partition pairs each receive three starts.

Every optimization uses one shared pair of physical types, all physical
letter heights, the common outer scale, and the place allocations. The
finite coefficient degree is charged in full. Selected conjugates are
counted separately in the Hermite scale, while restriction-of-scalars
channel multiplicity is divided out exactly as in the oracle theorem.

The run completed all 402 trials in `331.264724612` seconds, below its
1,200-second time ceiling. There were zero unresolved objective calls.
All 24 two-end warm-start controls agreed with their previous native
values before additional places were introduced.

This does not exhaust weighted gauges, physical bases, all possible
place partitions, or weighted effects at ordinarily regular places.
A positive constant-frequency allocation to an additional place can
make a bound worse; the comparison therefore retains the earlier
two-end result separately. Adaptive flags can ignore a place below the
cost at which its jets become useful.

## 2. Representative retained records

The following are floating diagnostic pressures. The source indices
refer to the two deflated primitive tensors in the restriction catalogue.

| Source | `X` labels | `Y` labels | Degree | Earlier two-end value | Selected extra-place value |
|---:|---|---|---:|---:|---:|
| 0 | 0, 2, 3 | 0, 2, 3 | 10 | 0.7606981113 | 0.8515910807 |
| 1 | 0, 1 | 0, 1, 3 | 9 | 0.7603367752 | 0.8540361316 |
| 0 | 0, 2 | 0, 2, 3 | 9 | 0.7603857771 | 0.8696311992 |
| 0 | 1, 3 | 0, 2 | 7 | 0.9964078298 | 0.8811905193 |

The first row is the best among this run's additional-place partitions.
The fourth illustrates improvement on an otherwise weak tensor.

## 3. Exact replay and retained evidence

The checkpoint binds the complete raw archive to every input and source
digest, verifies the complete trial-index set, and checks each retained
winner against the corresponding trial. The raw record is compressed
without changing its bytes: 1,895,784 uncompressed bytes become a
315,480-byte gzip archive.

For all 24 retained winners, an independent exact bordered-minor oracle
recomputes the complete local profiles at rational physical heights and
scales. The type probabilities and group allocations are rationalized
and normalized exactly. Every local scale is then derived from the one
common outer scale and its allocation, rather than rounded independently.

The replay checks 112 complete profiles with 559 channels, 784 strict
endpoint entries, and 224 saved allocation identities. The original
native profiles replay with zero discrepancy. The largest discrepancy
between exact rational orders and the native oracle at those same
rational parameters is below `7.2e-15`. Reconstructed original pressures
agree within `5.6e-16`; rationalization changes a displayed pressure by
at most `4.7e-14`.

The algebraic orders and rational mean shifts are exact. Exponentials,
entropies, and the final displayed pressures remain floating evaluations.
No threshold inequality is inferred from those decimals.

The sources are
[`screen_restricted_closed_partitions.py`](../tools/screen_restricted_closed_partitions.py),
[`search_restricted_closed_pressure.py`](../tools/search_restricted_closed_pressure.py),
and
[`checkpoint_restricted_closed_search.py`](../tools/checkpoint_restricted_closed_search.py).
The retained data are the
[adaptive screen](../evidence/restricted-closed-adaptive-screen.json),
[complete compressed search](../evidence/restricted-closed-weighted-pressure-search.json.gz),
[compact summary](../evidence/restricted-closed-weighted-pressure-summary.json),
and [exact winner replay](../evidence/restricted-closed-weighted-pressure-replay.json).

```sh
python3 tools/checkpoint_restricted_closed_search.py \
  evidence/restricted-closed-weighted-pressure-search.json \
  --archive evidence/restricted-closed-weighted-pressure-search.json.gz \
  --summary evidence/restricted-closed-weighted-pressure-summary.json \
  --output evidence/restricted-closed-weighted-pressure-replay.json
```

The input JSON in this command is the byte-for-byte decompression of the
committed gzip archive.
