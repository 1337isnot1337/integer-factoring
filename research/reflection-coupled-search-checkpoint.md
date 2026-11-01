# Completed coordinate-frame search and affine crossing paths

## Results and scope

The four-letter coordinate-frame search completed its full finite
candidate list. It did not improve the original reflection candidate.
The evidence retains every candidate, every optimizer outcome, and the
best resolved parameters from each trial. Independent reconstruction
and rational local-profile checks passed.

Affine combinations of the same frames give a larger exact family.
A concrete combination below changes a complete cut determinant and
creates two additional singular interpolation points without increasing
the polynomial degree. Its full contraction cost is still being tested.
Neither result establishes a factoring exponent above `1/2`.

## 1. The completed finite search

Start with the degree-24 reflection tensor on four physical letters.
In Bell coordinates use `G=c E_(i,j)`, with `i!=j` and `c` one of
the three nonzero elements of `F4`. The complete Laurent test accepts
102 negative-power frames. The
[dense rank-one catalogue](dense-rank-one-coupled-frames.md) independently
covers every positive pole order of this coordinate subfamily.

The search tests each accepted frame and the identity in both the
original and Bell physical bases. Thus it contains

```math
2(102+1)=206
```

candidate cases. All 206 cases completed their six optimization loops:
four balanced starts and two asymmetric starts. There were 1,236
returned optimizer calls, of which 925 reported successful convergence.
The remaining calls retained their best resolved complete-profile
values; they are not treated as convergence certificates. No candidate
case was omitted by the time limit. The recorded search duration is
954.942408119 seconds, ending when the finite list was exhausted.

The asymmetric weights have the form `(cos(theta),sin(theta))`, with
`0.03<=theta<=pi/2-0.03`. All physical-symbol heights and one common
pair of types are included. A balanced value below `1/sqrt(2)` or an
asymmetric supporting value below one would require an independent
strict certificate before using the
[rectangular factoring transfer](asymmetric-modular-transfer.md).

The retained numerical results are:

| Candidate scope | Balanced value | Asymmetric supporting value |
|---|---:|---:|
| Entire tested list, including the original seed | 0.7310809359734765 | 1.0014328543559963 |
| Nonidentity negative frames only | 0.7546022288178558 | 1.0019637253557476 |

These are achieved numerical values, not lower bounds over all gauges,
types, or supporting angles. The nonidentity balanced winner uses
`(i,j,c,k)=(1,2,1,2)` in original coordinates. The corresponding
asymmetric winner uses `(1,2,2,3)` in original coordinates.

The [replay](../tools/check_reflection4_coupled_search.py) reconstructs
all candidates through the separate arbitrary-matrix constructor,
checks their complete coefficient hashes and 52,736 endpoint entries,
and verifies every finite-catalogue key. It freshly evaluates four
retained winners and checks eight complete profiles after rationalizing
their height-to-scale ratios. Those exact local controls validate the
profile engine near the numerical witnesses; they do not certify a
global optimizer minimum.

Files:

- [Complete search record](../evidence/reflection4-coupled-frame-search.json).
- [Independent replay record](../evidence/reflection4-coupled-search-check.json).
- [Search source](../tools/search_reflection4_coupled_frames.py).

The separate [dense pressure record](../evidence/dense-rank-one-pressure-search.json)
screens all 7,176 admissible rank-one tensors, obtaining 153 distinct
natural complete profiles, then performs 64 selected gauge/type
refinements. It completes that requested work in 425.377722623 recorded
seconds. Its best refined balanced value is approximately 0.7546022267,
and its best supporting value is approximately 1.0019637816. Neither
improves the original seed. The
[archive checker](../tools/check_dense_rank_one_pressure.py) binds every
screen and refinement to the published exact catalogue and freshly
replays both winners. It preserves the historical input digest: the
later catalogue checks a larger pole range but has exactly the same
generator and candidate records.

## 2. Affine closure preserves the strict path

Let `T` and `T_1,...,T_s` be polynomial crossing paths over one field.
Assume that each has the same strict endpoint at one, the same degree
`D`, and the same coefficient tensor at degree `D`. For any constants
`a_1,...,a_s`, define

```math
T_a=T+\sum_{i=1}^s a_i(T_i-T).                    \tag{1}
```

Every difference is polynomial, vanishes at one, and has degree less
than `D`. Therefore (1) has the original endpoint, degree, and leading
coefficient. This proof is valid over any field.

The accepted negative-power physical frames satisfy these hypotheses:
they become identity at infinity and preserve the crossing at one.
Thus they supply an effective affine family with no extra degree.
The constructor adds every complete polynomial coefficient before any
rank calculation. It does not combine separate rank bounds as though
ranks were linear.

In particular, individual physical frames preserve the two cost-cut
determinants, but an affine combination need not do so. Complete
weighted profiles, common types, and all interpolation costs must be
computed again. Even a weighted profile preserved by every individual
frame need not remain unchanged under their sum.

## 3. An exact determinant-changing witness

Write `F4=F2(omega)`, with `omega^2+omega+1=0`. For the base tensor
`T`, let `T_c` be the coupled frame with Bell generator `c E_(2,1)`
and pole order one, restored to the original physical coordinates.
Take

```math
T_* = T+\omega(T_\omega-T)
         +\omega(T_{\omega^2}-T).                \tag{2}
```

In the frozen 102-frame implementation catalogue this is the word
`[[56,2],[58,2]]`, where the integer field encoding is
`omega=2`, `omega^2=3`. The complete resulting tensor has degree 24,
zero monomial content, and the exact crossing endpoint.

At each of the two points `z=omega,omega^2`, the complete unweighted
Smith profiles of (2) are:

| Cut | Local orders |
|---|---|
| `(x,y)` versus `(u,v)` | sixteen zeros |
| `(x,u)` versus `(y,v)` | eleven zeros, four twos, and one four |

The original tensor has rank 16 in both cuts at both points. The
new second cut has rank 11 and determinant valuation 12 at each point.
Therefore its determinant has changed. It cannot be obtained merely
by another pair of inverse-transpose physical frames, since every
such change preserves both complete cut determinants as rational
functions.

The [exact checker](../tools/check_affine_coupled_frames.py) independently
reconstructs (2) using the general square-zero matrix implementation,
compares all 256 coefficients and leading entries, and verifies the
endpoint. It translates the polynomials to each finite point using
an independent `F4` implementation, computes the ranks, and obtains
all four profiles by fraction-free elimination. The
[evidence](../evidence/affine-coupled-frames.json) binds the full tensor
and source hashes.

The two extra singular points can be used by a mixed-place construction,
subject to its full weighted rank and degree accounting. Their presence
alone is not a contraction certificate.

## 4. Reproduction and ongoing work

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/check_reflection4_coupled_search.py \
  --output /tmp/reflection4-coupled-search-check.json
OPENBLAS_NUM_THREADS=1 python3 tools/check_affine_coupled_frames.py \
  --output /tmp/affine-coupled-frames-check.json
```

The separate [affine search](../tools/search_affine_coupled_reflection_frames.py)
screens exact combinations and then optimizes complete two-place
balanced and asymmetric costs. New finite singularities also motivate
mixed-place profiles. These later searches are separate from the
completed 206-case record and are not claimed complete here.
