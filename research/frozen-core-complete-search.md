# Complete finite searches of the polynomial frozen-row family

## Result and scope

The [polynomial frozen-core family](polynomial-frozen-core-crossings.md)
has twelve coefficients in the four-element field and one endpoint
constraint, `g_0+g_10=1`. Thus it contains exactly
`4^11=4,194,304` parameter rows at each value of the field parameter.
The [compact compiler](frozen-core-row-profiles.md) makes both active
Smith matrices linear in these twelve coefficients.

The search in `tools/search_frozen_core_complete_pressure.py` covers
**every endpoint-valid row at each recorded gauge and rational scale**.
It uses exact affine counts, complete polynomial matrices, and outward
integer interval arithmetic. Rows with a potentially favorable lower
bound receive their complete integer-weighted Smith profiles. No row
is discarded using a floating heuristic.

This is a finite comparison with the Gaussian member at the **same**
gauges. It does not cover untested gauges or scales. It is not an
obstruction for the entire continuous search family, and it does not
change the factoring exponent.

## 1. A short coefficient flag determines the first order

Let `A(g,z)` be a matrix of Laurent polynomials, linear in a vector
`g` in a field of dimension `d`. Include the fixed integral row and
column gauges in the Laurent exponents. Every coefficient of every
entry is a linear form in `g`.

Scan these forms in nondecreasing exponent order. Retain a form only
when it is independent of the preceding retained forms. This produces
at most `d` forms, with nondecreasing weights:

\[
(\ell_0,w_0),\ldots,(\ell_{r-1},w_{r-1}),\qquad r\le d.
\]

For any nonzero specialized matrix, let `i` be the first retained
form with `ell_i(g)` nonzero. Then

\[
\min_{a,b}\operatorname{ord}A_{ab}(g,z)=w_i.
\]

Indeed, all earlier coefficient forms lie in the span of earlier
retained forms, so they vanish on `g`. The retained form itself is an
actual coefficient at weight `w_i`. Forms of equal weight can be
processed in any fixed order. A vanishing coefficient vector produces
the zero matrix and has no first nonzero form.

In our family the endpoint identity excludes a zero specialized cut.
Two coefficient flags therefore partition every valid parameter row
into at most `d^2=144` cells. Cell `(i,j)` imposes the earlier forms
from both flags as zero and requires both selected forms to be nonzero.

## 2. Exact affine counting and parameter fibers

Intersect each cell with the endpoint affine space. If `V` is the
space satisfying the earlier zero constraints, and `ell,h` are the
selected forms, the exact cell size is

\[
|V|-|V\cap\ker\ell|-|V\cap\ker h|
   +|V\cap\ker\ell\cap\ker h|.
\]

Every term is either zero for an inconsistent system or a power of
four obtained by finite-field elimination. The implementation checks
that the disjoint cell counts sum to `4^11`.

Distinct parameter rows can define the same pair of matrices. Let
`K` be the common kernel of all their coefficient forms. We first
verify that the endpoint functional lies in the span of those forms.
Consequently, every fiber of the coefficient map stays within one
endpoint level. If the joint coefficient rank is `r`, each valid
fiber has size `4^(12-r)`.

An echelon basis selects pivot columns of the coefficient map. Setting
all other parameter coordinates to zero gives a section: the map is
injective on these pivot coordinates, and its image has dimension `r`.
We enumerate the section and multiply every count by the common fiber
size. At field parameter zero, the tested matrices have rank nine,
so each specialized tensor represents 64 original parameter rows.
At field parameter one the joint rank is ten, giving fibers of size 16.

## 3. Certified lower bounds reject whole cells

Let the two exact minimum entry orders in a cell be `m_0,m_inf`.
They are also the first Smith orders. For a scale `s>1`, define the
active profile partition functions

\[
Z_\nu(g;s)=\sum_i s^{-\lambda_{\nu,i}(g)}.
\]

Every cell member satisfies

\[
Z_0(g;s)Z_\infty(g;s)\ge s^{-m_0-m_\infty}.
\]

If this lower bound is at least the Gaussian active product at the
same gauges and scale, the entire cell cannot improve that product.
The Gaussian spectators and physical type partition functions are
identical for all rows at these gauges, so this comparison also
preserves the direction of the full pressure comparison.

There is a stronger row-specific bound. Normalize each cut matrix
by its least entry order and reduce modulo the local parameter. If
the resulting ranks are `r_0,r_inf`, exactly these numbers of Smith
orders equal the respective minimum orders. Hence

\[
Z_0(g;s)Z_\infty(g;s)
\ge r_0r_\infty s^{-m_0-m_\infty}.
\]

These ranks use exact arithmetic in the four-element field. Every
row receiving a complete profile also checks both the minimum order
and its multiplicity against the independent leading matrix.

## 4. Complete remaining profiles and exact comparisons

All coefficient matrices come from the complete polynomial compiler.
The search rejects a retained-series matrix before using it. Fixed
rational gauges are accepted only when their weighted matrix orders
are integral; another fixed ramification can handle other rational
choices.

For every row not certified by a lower bound, the native matrix engine
computes the complete integer Smith orders. An unresolved precision
request is recorded as unresolved, never as an exclusion. Time and
query caps likewise leave explicit unresolved counts. A run is marked
complete only if this count is zero.

The scale is rational. Every partition comparison first uses 192-bit
outward dyadic intervals. A strict separation decides the sign. If
the intervals overlap, the code evaluates both products as exact
rational numbers. Identical profile lists are recognized directly.
Thus roundoff cannot create either a false improvement or a false
exclusion in a completed run.

## 5. Reproduction

For example:

```sh
python3 tools/search_frozen_core_complete_pressure.py \
  --M 4 --beta 1 --grid=-2,-1,-1/2,0,1/2,1,2 \
  --scale 40813173/25000000 --output frozen-M4.json
python3 tools/search_frozen_core_complete_pressure.py \
  --M 32 --beta 0 --grid=-1,0,1 \
  --scale 58275883/50000000 --output frozen-M32-zero.json
python3 tools/search_frozen_core_complete_pressure.py \
  --M 64 --beta 1 --grid=-1,0,1 \
  --scale 13703179/12500000 --output frozen-M64-one.json
```

The field-orbit theorem in the constructor identifies every nonzero
field parameter with parameter one after a corresponding row change.
Therefore the parameter-zero and parameter-one searches represent the
two field-parameter orbits. The conclusion remains limited to the
stated gauges and scales in both orbits.

## 6. Completed finite results

The [complete evidence](../evidence/frozen-core-complete-search.json)
records the following 86 gauge/scale/orbit cases. Every case covers all
4,194,304 valid parameter rows. No case found a product below the
same-gauge Gaussian product; every unresolved count is zero.

| Interval `M` | Field parameter | Each gauge ranges over | Gauge pairs | Exact rational scale |
| ---: | ---: | --- | ---: | --- |
| 4 | 0 | `0` | 1 | `40813173/25000000` |
| 4 | 1 | `-2,-1,-1/2,0,1/2,1,2` | 49 | `40813173/25000000` |
| 32 | 0 and 1 | `-1,0,1` | 9 per orbit | `58275883/50000000` |
| 64 | 0 and 1 | `-1,0,1` | 9 per orbit | `13703179/12500000` |

The [checker](../tools/check_frozen_core_complete_search.py) independently
compares 480 flagged minimum orders with complete coefficient scans,
checks 100 rational comparison signs, and repeats both complete M4
zero-gauge orbit searches. Its
[control record](../evidence/frozen-core-complete-search-controls.json)
keeps representative profile queries distinct from the original
parameter rows covered by their coefficient fibers.
