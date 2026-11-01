# Complete singular-place catalogues of two finite-fiber cores

The two primitive seeds retained by the finite-frame search have no
nonzero finite singular places, over any field extension. Their complete
cut determinants are

| Deflated degree | Row determinant | Column determinant |
| --- | --- | --- |
| 10 | `w^112` | `1` |
| 11 | `w^120` | `1` |

The tensors are the two source-bound cases in
[`deflated-finite-fiber-input.json`](../evidence/deflated-finite-fiber-input.json).
They come from the degree-twenty and degree-twenty-two primitive
[finite-frame candidates](finite-jet-smith-fibers.md), respectively, by
the exact substitution `T(w)=R(w^2)`.

This calculation closes a specific search gap: neither core has an
unexamined singularity at a quadratic or higher-degree closed point.
Every invertible constant physical basis change multiplies a cut
determinant by a nonzero constant, so it preserves this conclusion.
It does not constrain subsequent polynomial modifications of the tensor.
The catalogue is a statement about ordinary local singularities. Physical
weights can still exploit a typed support pattern at a regular place;
determinant nonvanishing does not exclude every weighted finite-place
construction.

## 1. Complete local orders

The following orders use homogeneous reversal `t^d R(1/t)` at infinity,
where `d` is the displayed deflated degree.

| Degree | Cut and place | Complete local Smith orders |
| --- | --- | --- |
| 10 | Row, zero | `0,2,4,4,4,6,6,6,8,8,8,10,10,10,12,14` |
| 10 | Row, infinity | `0,0,1,1,2,2,3,3,3,3,4,4,5,5,6,6` |
| 10 | Column, zero | sixteen zeros |
| 10 | Column, infinity | `0,2,4,5,6,7,9,9,11,11,13,14,15,16,18,20` |
| 11 | Row, zero | `0,2,4,4,4,6,6,8,8,8,8,10,12,12,12,16` |
| 11 | Row, infinity | `0,1,1,2,2,3,3,3,4,4,4,5,5,6,6,7` |
| 11 | Column, zero | sixteen zeros |
| 11 | Column, infinity | `0,3,4,6,7,7,8,11,11,13,13,15,17,19,19,23` |

All other finite places have sixteen unit channels in both cuts. For
the original degree-twenty and degree-twenty-two seeds, the displayed
zero and infinity orders are doubled. The determinant exponents are
therefore 224 and 240, respectively, in the row cut; both column
determinants remain one.

## 2. Exact polynomial certificates

The [checker](../tools/check_finite_fiber_closed_places.py) first computes
each determinant and each complete profile by fraction-free local
elimination. It then verifies the determinant identities independently.

A sixteen-by-sixteen polynomial matrix with entry degree at most `d`
has determinant degree at most `16d`. The checker constructs `F256`,
embeds `F4` by an explicit root of `X^2+X+1`, and evaluates the matrix
at `16d+1` distinct field points. Ordinary field elimination agrees
with the claimed monomial at every point. Since `16d+1<=177<256`,
these equalities prove the polynomial identity. This check does not
use the polynomial elimination or its pivot sequence.

For each of the eight displayed local matrices, the checker also
forms its multiplication matrix modulo `t^L`. Its rank is compared
with

```math
\sum_{i=1}^{16}\max(L-\nu_i,0)
```

for every `L` through one beyond the largest asserted order. These
ranks determine the entire Smith-order multiset, and provide a second
calculation of every local profile. The evidence records all ranks,
field data, tensor hashes, and source hashes.

The determinant identities imply the claimed complete catalogue:
a polynomial matrix is singular at a closed point exactly when its
determinant vanishes there. The monomial has only the zero place,
and the constant determinant has none. Homogeneous infinity is
handled by its separately computed profile.

## 3. Scope and replay

These are complete algebraic controls for two explicit candidates.
They do not improve their optimized crossing pressure, and do not
establish a factoring exponent above one half. The ongoing physical
basis search may use the existing zero and infinity places with
confidence that no extra finite place has been omitted.

```sh
python3 tools/check_finite_fiber_closed_places.py \
  --output evidence/finite-fiber-closed-places.json
```

The [replay record](../evidence/finite-fiber-closed-places.json) binds the
calculation to the complete deflated input artifact.
