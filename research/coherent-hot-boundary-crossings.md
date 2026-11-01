# Coherent energy-boundary codes and an eight-letter certificate

## Result

An explicit paired compression of a stabilized reflection tensor gives
an eight-letter strict polynomial crossing of degree 56 in `w`, with
the uniform balanced bound

```math
\rho<\frac{18}{25}=0.72.                              \tag{1}
```

All tensor coefficients and both complete natural profiles are exact.
The proof of (1) is a strict integer comparison. The accompanying
finite controls improve the recorded coordinate-boundary choices;
they do not establish optimality among eight-letter crossings. Bound
(1) remains above `1/sqrt(2)` and does not improve the current factoring
exponent.

The code construction allows coherent linear combinations at a fixed
physical energy. It extends the coordinate restrictions used in
[stable-hot-reflection-codes.md](stable-hot-reflection-codes.md).

## 1. Paired compression at a common energy

Fix a finite hot energy level `E`. Retain every physical state with
energy below `E`, and let `W` be the full coordinate space of states
at energy `E`. Choose any two matrices over `F4`,

```math
K:W\longrightarrow F_4^r,\qquad J:F_4^r\longrightarrow W,
\qquad KJ=I_r.                                       \tag{2}
```

Extend them by the identity on all lower-energy states. On a paired
crossing apply the four maps `(K_X,K_Y,J_Y^T,J_X^T)` to `(x,y,u,v)`.
The literal crossing endpoint is preserved by (2), coefficient by
coefficient. Both maps commute with the natural diagonal weights,
because each boundary state has the same energy.

The natural compound reduction remains exact. To see this, write the
regular natural matrix as `U E_B V`, as in
[reflection-natural-cut-duality.md](reflection-natural-cut-duality.md).
Every nonidentity entry of `U` or `V` deletes occupied bits. Any
deletion from a boundary state lands strictly below energy `E`, where
every coordinate state is retained. A Gaussian deletion coefficient
depends only on the removed bits; a leg with no deletion is acted on
by the identity. Equivalently, these factors are polynomials in the
individual bit-deletion operators. Hence the row space of `K` is
preserved by the relevant right action of `U`, and the image of `J`
is preserved by the left action of `V`. Their product-code versions
give

```math
K_R U=U_C K_R,\qquad VJ_C=J_C V_C,
\qquad K_R U E_B VJ_C=U_C(K_R E_BJ_C)V_C.              \tag{3}
```

The induced matrices `U_C,V_C` are triangular by physical energy with
diagonal one. They are polynomial-unimodular, so (3) preserves all
Smith orders. This proof permits arbitrary state coefficients in
(2); it does not require coordinate boundary maps.

At infinity the natural duality exchanges `y,u`. Consequently its
compound matrix uses the dual boundary data
`(K_X,J_X; J_Y^T,K_Y^T)`. These maps again satisfy (2). The zero and
infinity compound profiles can differ; they coincide automatically
for self-dual data `K=J^T`. Their opposite scalar monomial shifts
still cancel in the complete two-place pressure.

If the boundary basis mixes states of different occupation numbers,
the compressed matrix need not split into the old occupation blocks.
The complete profile must then be recomputed without that split.

## 2. The explicit eight-letter crossing

Use interval `M=16`, hot energies `E(a)=8a+4`, and retain the seven
coordinate states

```text
0, 1, 2, 3, 4, 5, 8.
```

Their energies are `0,4,12,16,20,24,28`. The energy-32 boundary has
the two states `6,9`. On this boundary use

```math
K=(1\ \ 1),\qquad J=\binom{0}{1}.                    \tag{4}
```

The same maps are used for the two paired alphabets. Thus the final
physical energy list is

```text
0, 4, 12, 16, 20, 24, 28, 32.
```

All maps and the resulting tensor have binary coefficients. The
complete tensor has scalar content `w^104`; after division its actual
degree is 56 and its SHA-256 digest is

```text
17254d71e2fece406c199fba1988424f3619911d5cc93f01c900d43a5c8bdbf9
```

The vacuum coefficient before division is `w^128`, so the common
scalar gcd is a monomial. The checker computes its exact exponent
from every coefficient. Its two natural primitive profiles are
`-24+sigma_infinity` and `24+sigma_0`, with 64 orders each. Both
nonnegative lists appear in the evidence and are independently
recovered from the literal tensor and the two different compound
compressions in (3).

## 3. Strict rational pressure and all-width transfer

Set `z=141/142`, `A=36`, `B=25`, `R=A/B`, and
`s=A log(1/z)`. Define

```math
P_\infty=\sum_i z^{A\sigma_{\infty,i}},\qquad
P_0=\sum_i z^{A\sigma_{0,i}},\qquad
Q=\sum_x z^{B E(x)}.
```

The checker proves the exact inequality

```math
(P_\infty P_0)^B<Q^{2A}.                              \tag{5}
```

After denominator clearing, both positive integers have 720,752
bits. Their difference is strictly positive. No numerical minimum
or logarithm decides (5).

Choose the positive rational common type `p(x)=z^(BE(x))/Q`.
Its entropy satisfies `H=log Q+(s/R) E_p`. The opposite scalar shifts
cancel, while the full physical mean contributes `2s E_p`. Hence

```math
K_{\rm full}=\log(P_\infty P_0)+2sE_p<2R H,
\qquad K_{\rm full}/(4H)<R/2=18/25.                  \tag{6}
```

The exact type counts and strict integer gap make a sufficient finite
typed block effective. The coefficient prefactors and interpolation
piece count are polynomial in its length, and the strict gap absorbs
them. Applying
[finite-block-crossing-amplification.md](finite-block-crossing-amplification.md)
then gives (1) for every external label width, including all query
and factorization costs specified there.

## 4. Reproduction and controls

```sh
python tools/check_coherent_hot_boundary.py \
  --output evidence/coherent-hot8-pressure.json
```

The [checker](../tools/check_coherent_hot_boundary.py) compares 12,288
complete coefficients from the accelerated compression with the
independent four-leg retraction implementation. Two native interfaces
agree on all 128 natural channels of the primitive eight-letter seed.
Literal Pfaffian tensors and compressed compound matrices agree on
328 channels in total, including a boundary that mixes one- and
three-particle states. Its complete matrix is used without an invalid
occupation split.

An unequal-energy countercontrol remains a strict crossing but fails
the unchanged-natural-gauge shortcut. The fixed-energy hypothesis is
therefore enforced in the theorem. The
[evidence](../evidence/coherent-hot8-pressure.json) contains the full
primitive tensor, every profile, the exact type counts, integer hashes,
and source hashes.

The constructor [coherent_hot_boundary.py](../tools/coherent_hot_boundary.py)
uses object-array XOR only to accelerate exact `F4` polynomial arithmetic.
The separate [finite search](../tools/search_coherent_hot_boundaries.py)
enumerates projective right vectors and every left vector pairing to
one. It records complete profiles; every new numerical winner is also
reconstructed from the original literal reflection tensor before it
is retained. These finite numerical searches are not optimality proofs.

## 5. Completed finite boundary grid

The [full archive](../evidence/coherent-hot-boundary-search.json.gz)
contains every rank-one idempotent over `F4` on the indicated boundary.
Normalize the first nonzero coordinate of the right vector to one,
then enumerate every left vector pairing to one. A boundary of
dimension `b` has `((4^b-1)/3)4^(b-1)` such idempotents. This gives
the following exact finite coverage.

| Physical letters | Boundary dimension | Idempotents | Best recorded natural pressure | Winner degree |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 2 | 20 | 0.7199753178555 | 56 |
| 10 | 2 | 20 | 0.7177834990515 | 72 |
| 12 | 2 | 20 | 0.7164001436661 | 80 |
| 16 | 3 | 336 | 0.7147049456406 | 96 |

The last case of the larger worker was absent from its saved prefix.
That prefix remains marked incomplete. One explicit continuation
computed the missing case; the checker verifies that their disjoint
union exhausts the grid. It does not label the prefix as an
uninterrupted completed run. All three original JSON byte streams,
including their timing and source metadata, are preserved in the
single compressed archive.

```sh
python tools/checkpoint_coherent_hot_boundaries.py
```

The [archive checker](../tools/checkpoint_coherent_hot_boundaries.py)
verifies all 396 distinct idempotents and recomputes every saved
pressure at its recorded scale. It independently reconstructs the
four full winning tensors from the literal reflection source, with
100,368 coefficient comparisons and 1,128 complete profile channels.
The [replay](../evidence/coherent-hot-boundary-search-replay.json)
separates these exact checks from numerical pressure evaluations.
It does not claim fresh literal-tensor reconstruction of every
discarded idempotent, or global optimality of any scale minimization.
