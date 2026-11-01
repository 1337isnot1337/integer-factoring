# Affine combinations across filtered sixteen-letter directions

## Result and scope

A completed forty-minute search tested 849 distinct affine combinations
of admissible rational-frame differences. It combined different
filtered physical directions and different pole orders. Every retained
tensor had degree 96 and the strict crossing endpoint. Every requested
natural-profile computation succeeded.

The best retained natural pressure was `0.8128910933978555`, above the
original seed's `0.7166536597281867`. This search therefore produced no
improvement. It is a bounded search at one physical-height shape, with
the common Gibbs type and one scale optimized. It does not exclude
other physical gauges, additional finite places, or other affine words.

The complete run is preserved in the
[compressed archive](../evidence/filtered-reflection16-affine-combinations-search.json.gz).
The [replay record](../evidence/filtered-reflection16-affine-combinations-replay.json)
binds that archive to the exact source versions and reconstructs every
retained polynomial tensor.

## 1. Exact component family

Start with the sixteen-letter reflection tensor `T` for interval size
eight. Its degree is 96. Use the five filtered rank-one directions
from the [one-direction checkpoint](filtered-reflection16-affine-search.md).
The requested pole orders were

```text
1, 2, 3, 4, 6, 8, 12, 16.
```

For a given direction and pole, let `T(c)` be the complete coupled
frame with nonzero scalar `c` in `F4`. Retain a pole only after all
three scalar frames are polynomial. The complete four-leg expansion
has the form

```math
T(c)-T=cU_1+c^2U_2+c^3U_3.                       \tag{1}
```

The fourth power of `c` has been combined with its first power,
using `c^4=c`. The three-by-three matrix evaluating `c,c^2,c^3`
at the three nonzero field elements is invertible. Consequently all
three components in (1) are polynomial whenever the three complete
scalar frames are polynomial. This checks the higher-pole cancellation
condition before an individual component enters the catalogue.

Remove zero components and identify scalar multiples. Normalize the
first nonzero coefficient to one. The resulting catalogue has 52
distinct nonzero components `Delta_i`. Each has degree less than 96,
zero endpoint at one, and zero paired traces.

For any finite coefficient word, define

```math
S=T+\sum_i c_i\Delta_i,\qquad c_i\in\mathbb F_4.   \tag{2}
```

Every such tensor has the original leading coefficient and strict
endpoint. The paired-trace identity remains
`(1+w^64) I` in both directions. Its common scalar content is one:
a common factor would divide `1+w^64=(1+w)^64`, while the endpoint
at one is nonzero. The construction does not remove degree or scalar
content after forming a word.

## 2. The completed run

The [search program](../tools/search_filtered_reflection16_affine_combinations.py)
started at `2026-10-31T23:00:37.972447+00:00` and finished at
`2026-10-31T23:40:38.817577+00:00`. Its actual elapsed time was
`2400.8451536269995` seconds, against a declared 2,400-second budget
and a separate 12,000-proposal ceiling. The time budget ended the run.

| Quantity | Completed count |
|---|---:|
| Exact normalized components | 52 |
| Proposed affine words | 856 |
| Duplicate tensors | 7 |
| Distinct retained tensors | 849 |
| Unresolved profile requests | 0 |
| Distinct profiles, including the original control | 675 |

Fresh words used two to five distinct components. The search also
mutated a retained pool of up to 24 distinct profile classes. Every
word, exact tensor hash, profile assignment, and proposal index is
recorded. The pseudorandom seed was `312401`.

Both complete weighted cuts were recomputed. In particular, the
infinity spectrum was not inherited from a component or assumed
unchanged under addition. Each profile records all 512 channels.
The physical energies were the subset sums of `(28,20,12,4)`.
The common type was obtained from the Gibbs pressure equation, and
the shared scale was optimized in the declared interval.

The best retained word was

```text
[[0, 1], [1, 3]]
```

in the recorded component catalogue. Its complete tensor hash is

```text
dedc96dd29bd7d7aa6f4e15e86847fc90a67becd3876a5265b4da2b44785e7f9
```

Its natural pressure was `0.8128910933978555`. The original control
was `0.7166536597281867`. These are numerical diagnostic pressures;
they are not rational crossing certificates.

## 3. Exact reconstruction and independent controls

The [checker](../tools/check_filtered_reflection16_affine_combinations.py)
checks all nine source hashes recorded by the search. It reconstructs
the catalogue and all 849 retained tensors, including every endpoint,
degree, scalar monomial content, and complete polynomial hash. This
covers 55,640,064 polynomial tensor entries.

For the winner, the checker independently inverts the three-scalar
matrix in (1) by exhaustive field arithmetic and reconstructs its
components from the three complete scalar frames. That replay does
not use the search's affine-component formula. It compares 131,072
complete polynomial entries.

Finally, the checker recomputes all 512 natural channels for the
original tensor and the best retained tensor, then re-evaluates their
reported pressures at the recorded scales. The other retained
numerical profiles remain the source-bound completed search record;
the checker does not claim to optimize every gauge or to replay every
numerical minimization.

```sh
python3 tools/check_filtered_reflection16_affine_combinations.py \
  --output evidence/filtered-reflection16-affine-combinations-replay.json
```

The archive preserves the exact completed runtime and all retained
trials. This is separate from the earlier one-direction search: it
permits sums across different directions and pole orders, while using
the same explicit polynomial and endpoint gates.
