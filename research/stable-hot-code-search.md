# Completed searches on stabilized reflection codes

## Result

The [stabilized hot-code construction](stable-hot-reflection-codes.md)
was screened on all eleven last-energy choices for physical alphabet
sizes 8, 12, 16, 20, and 24. Every tensor and scalar-content division was
constructed exactly. The best natural-gauge numerical value was
`rho = 0.7129641787845314`, on the twenty-four-letter code of degree 112.

A separate full physical-gauge and supporting-line search refined the
first eight-letter and sixteen-letter examples. Its best support value
was `1.0098488477041911`, still above the required value one. Neither
search supplies an improved factoring exponent or a global exclusion.
The [separate rational certificate](stable-hot-code-rational-pressure.md)
proves the constructive upper bound `rho<5/7` for the twenty-four-letter
example without relying on numerical optimization.

## 1. Complete natural-energy boundary screen

For each declared size, order interval-sixteen physical labels by
natural energy. Include every lower-energy label and enumerate every
subset of the final energy level that gives the required size. Thus
ties are retained as separate exact tensors rather than resolved by
one arbitrary convention.

| Alphabet size | Boundary choices | Best numerical `rho` | Primitive degree of that record |
|---:|---:|---:|---:|
| 8 | 2 | 0.7200216153415486 | 64 |
| 12 | 2 | 0.7164001436660765 | 80 |
| 16 | 3 | 0.7147357780605276 | 96 |
| 20 | 3 | 0.7136590278856072 | 104 |
| 24 | 1 | 0.7129641787845314 | 112 |

The scalar search uses the complete two natural cut profiles and
optimizes the common paired physical type through its Gibbs equation.
It includes every channel, the reciprocal infinity orders, and both
physical means. It does not use only a finite-point rank or an
unrestored profile.

All eleven declared cases completed, with 188 resolved pressure calls
and zero unresolved calls. The run began at
`2026-11-01T00:30:54.316411+00:00`, ended at
`2026-11-01T00:34:39.105955+00:00`, and recorded
`224.789258026` elapsed seconds. Its declared ceiling was 900 seconds.

The degree-88 sixteen-letter example in the construction note is the
first boundary choice. The best degree-96 boundary choice has a
slightly lower natural pressure; both exact tensors are retained by
their recipes and hashes.

## 2. Full physical-gauge refinement

The refinement started from the first eight-letter and sixteen-letter
codes, with natural values `0.7200277848498949` and
`0.7147366313237278`. For each it used three initial supporting-line
angles, `0.25`, `pi/4`, and `1.32`. Every physical letter height, the
shared scale, the supporting-line angle, and the common Gibbs type
were allowed to change.

Each optimization had a 120-second ceiling. The three sixteen-letter
trials reached that ceiling and retained their best resolved values.
The three eight-letter optimizations completed. Six unresolved
objective calls were penalized and were not interpreted as rank gains.

| Starting tensor | Best retained support value |
|---|---:|
| Sixteen letters, degree 88 | 1.0098488477041911 |
| Eight letters, degree 64 | 1.0163097282391127 |

The run began at `2026-11-01T00:30:04.510176+00:00`, ended at
`2026-11-01T00:36:19.240105+00:00`, and recorded
`374.729903989` elapsed seconds. This checkpoint contains no full-gauge
optimization of the twenty-letter or twenty-four-letter tensors.

## 3. Independent replay and retained records

The checkpoint reconstructs all eleven exact tensors and scalar
contents, covering 1,058,048 complete polynomial entries. It recomputes
the 22 complete natural integer profiles through separately assembled
matrix interfaces, retaining all 5,920 channels. Every one of the 188
recorded natural pressure calls is then replayed from those integer
profiles. The largest discrepancy is below `8.5e-14`.

Both retained full-gauge winners are also recomputed through the matrix
interface, covering four complete profiles and 640 channels. Their
physical means and shared type entropies are independently restored.
The largest order discrepancy is below `3.6e-15`; the support values
agree within `2.3e-16`. These arbitrary weighted orders and final costs
remain numerical controls. The natural integer profiles and polynomial
coefficient reconstructions are exact.

The raw natural search, full-gauge search, and original full-gauge input
are preserved byte for byte in deterministic gzip archives. Input and
source hashes are checked before archival. The replay also verifies the
entire declared boundary-choice list and all six gauge-trial keys.

Sources:

- [Natural code search](../tools/search_implicit_reflection_codes.py)
- [Existing full-gauge worker](../tools/refine_finite_fiber_restrictions.py)
- [Archive and replay checker](../tools/checkpoint_implicit_reflection_search.py)

Evidence:

- [Natural search archive](../evidence/implicit-reflection-code-search-natural.json.gz)
- [Full-gauge search archive](../evidence/implicit-reflection-code-search-gauge.json.gz)
- [Original full-gauge input](../evidence/implicit-reflection-code-search-gauge-input.json.gz)
- [Compact summary](../evidence/implicit-reflection-code-search-summary.json)
- [Complete replay profiles](../evidence/implicit-reflection-code-search-replay.json)
