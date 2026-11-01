# Exact paired restrictions of the primitive finite-fiber cores

All 200 coordinate restrictions to paired alphabets of size two or three
have been constructed and screened. The catalogue includes unequal
paired alphabets. Every restricted endpoint and scalar-content division
is exact. No tested restriction passed the stronger crossing threshold.

The smallest recorded balanced pressure is `0.7416676987094379`. A
separate supporting-line refinement on eight selected restrictions
remained above one. These are bounded numerical search results, not
global lower bounds on the restricted families.

## 1. Construction and complete scope

The two parents are the degree-ten and degree-eleven cores in
[`deflated-finite-fiber-input.json`](../evidence/deflated-finite-fiber-input.json).
For each parent choose subsets `S_x,S_y` of its four physical letters,
each of size two or three. Restrict the `x,v` legs to `S_x`, and the
`y,u` legs to `S_y`.

The endpoint remains

```math
T(1)_{xyuv}=\delta_{xv}\delta_{yu}
```

because each paired inclusion and projection compose to the identity.
The restricted alphabets have dimensions
`(|S_x|,|S_y|,|S_y|,|S_x|)`. There are

```math
2\left(\binom42+\binom43\right)^2=200
```

cases: 72 with dimensions `2,2,2,2`, 32 with dimensions `3,3,3,3`,
and 96 with unequal paired sizes.

Each tensor is divided by its complete common polynomial factor `g`.
The endpoint implies `g(1)!=0`; multiplying the quotient by `g(1)`
restores the exact endpoint. This normalization is included before
any local profile or pressure calculation.

The hot restriction `S_x=S_y={0,2}` is identical for both parents.
After removing its scalar `w^2`, it has degree four. Its complete
profiles, although smaller, do not provide a better recorded pressure.

## 2. New ordinary singular places

The full parent cores have only zero and infinity as ordinary singular
places, as proved in the
[parent catalogue](finite-fiber-singular-catalogues.md). Restrictions
can introduce new determinant factors: 115 of the 200 cases have at
least one nonzero finite singular place.

The [complete restricted catalogue](../evidence/finite-fiber-restriction-catalogue.json)
contains the two determinant polynomials, their complete factorizations
over `F4`, both local Smith profiles at every finite closed point, and
both homogeneous infinity profiles for every case. The nonzero factors
have residue degrees from one through seven.

The [catalogue program](../tools/catalogue_finite_fiber_restrictions.py)
first factors over `F16`, descends conjugate factor pairs by the fourth
power Frobenius, verifies each descended polynomial irreducible over
`F4`, and reconstructs each determinant exactly. Local profiles use
the original polynomial matrices and those factors directly.

For example, the first parent's restriction to
`S_x=S_y={0,1,2}` has degree nine and two new row-cut places,
`w=omega` and `w=omega^2`. Both have profile `[0 x 8,2]`; their column
profiles are all zero. These are additional search candidates, with
all other channels retained in their costs.

This catalogue concerns ordinary local singularities. Physical weights
may still exploit typed support patterns at regular points. Neither
the catalogue nor the completed two-end search excludes that possibility.

## 3. Recorded pressure searches

The balanced screen ran from `2026-11-01T00:06:19.034156Z` to
`2026-11-01T00:11:48.391074Z`, with saved elapsed time
`329.357161567` seconds. All 200 candidate tensors were processed.
There were 131 completed per-candidate optimizer runs and 69 which
reached their declared 2.5-second cap. Every candidate retained at
least one valid complete-profile evaluation.

The source is
[`search_finite_fiber_restrictions.py`](../tools/search_finite_fiber_restrictions.py),
and the [full record](../evidence/finite-fiber-restriction-search.json)
retains every tensor, subset, removed scalar, parameter vector, common
type, outcome, and source hash. It records 98,523 objective calls,
including 25,920 unresolved calls. An unresolved call receives a
penalty and is never interpreted as a rank saving.

The best recorded balanced values across the two parents are:

| Paired alphabet sizes | Best recorded rho |
| --- | ---: |
| 2 by 2 | 0.7602278138449967 |
| 2 by 3 | 0.7551307613607177 |
| 3 by 2 | 0.7551383967925278 |
| 3 by 3 | 0.7416676987094379 |

A frozen [selection of eight candidates](../evidence/finite-fiber-restriction-selected-input.json)
was then refined with unequal positive cut weights. The
[supporting-line worker](../tools/refine_finite_fiber_restrictions.py)
completed all 24 declared starts in `46.094739478` seconds; one start
reached its six-second cap. Its smallest recorded support value was
`1.0006629847738142`, above the required strict value one.
The [full asymmetric record](../evidence/finite-fiber-restriction-asymmetric-search.json)
preserves all heights, scales, angles, paired types, and outcomes.

These searches do not include every new closed place. The
[closed-place weighted oracle](closed-place-weighted-oracles.md)
provides the exact interface for that next search, including the
unequal paired dimensions and higher-degree residue fields.

## 4. Independent replay

The [checker](../tools/check_finite_fiber_restrictions.py) reconstructs
all restrictions by direct coordinate indexing, independently of the
restriction routine. It verifies the scalar reconstruction and all
7,200 endpoint entries.

For every cut determinant, it evaluates the matrix at more distinct
`F256` points than its degree bound. The resulting 20,454 ordinary
field determinants independently prove all 400 polynomial identities.
Twenty-two quotient-ring ranks cover representatives of every residue
degree from one through seven.

The checker then treats every winning height and scale as an exact
decimal rational, computes all 32 complete local profiles independently,
normalizes the common paired types, and restores their means. It
recomputes the balanced and asymmetric costs numerically, including
the homogeneous infinity correction. These final logarithmic costs
are floating values, not interval certificates.

The [replay evidence](../evidence/finite-fiber-restriction-replay.json)
binds every input, search, catalogue, selected-prefix case, and source.
The selected-prefix digest is preserved, while each selected tensor
and parameter vector is also matched to its immutable record in the
completed search.

```sh
python3 tools/check_finite_fiber_restrictions.py \
  --output evidence/finite-fiber-restriction-replay.json
```
