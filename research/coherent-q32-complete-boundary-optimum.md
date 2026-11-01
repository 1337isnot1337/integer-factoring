# The complete binary height-16 boundary optimum

## Exact finite-class theorem

The 32-letter coherent map with boundary kernel `7` and covector `15`
minimizes the natural product pressure at every `0 < r < 1` among all
496 binary rank-four retractions on the five-dimensional height-16
boundary. Its transposed projector, with kernel `15` and covector `7`,
has the same product polynomial. Every other map has strictly larger
product pressure throughout that interval.

This extends the previous 121-map occupation-preserving theorem to the
entire fixed boundary class, including all 375 mixed-occupation maps.
The exact certificate remains `rho < 1781/2501`; no factoring exponent
above one half follows. The theorem concerns the declared binary class
and natural gauges. It does not cover extension fields, changes below
the boundary, other ranks, or general physical letter gauges.

## 1. Why the class has exactly 496 projectors

Retain the 28 coordinate states of height less than 16. The five states
at height 16 are the two-particle words `24,36,66,129` and the
four-particle word `15`. On their five-dimensional span, choose an
arbitrary rank-four idempotent over `F_2`.

Every such idempotent has the unique form

\[
\Pi=I-w\ell,\qquad w\ne0,\qquad \ell(w)=1.
\]

There are 31 nonzero vectors `w` and 16 functionals taking value one on
each `w`, giving `31*16=496` projectors. Selecting a basis of the image
and its paired retraction gives `P,L` with `LP=I` and `PL=Pi`.
Changing that basis multiplies every physical leg by an invertible
constant matrix. It preserves both natural Smith profiles, and all four
boundary labels still have height 16. Thus counting projectors covers
every basis choice within this class.

Every boundary derivative has smaller height and lies in the fully
retained lower part. Both map spaces are therefore deletion closed.
Mixed occupation changes the complete block split from occupation to
parity, as proved in the [mixed-boundary construction](coherent-hot-mixed-boundary.md).

## 2. Complete profile accounting

The [finished catalogue](../evidence/coherent-q32-mixed-boundary.json)
contains all 496 exact map recipes and both complete 1,024-order profiles.
It reuses the already replayed 121 occupation-preserving profiles and
calculates the other 375 maps with the exact parity-block constructor.
Every new cut has blocks of sizes 514 and 510. All computations finish
with zero unresolved maps.

The literal endpoint separately proves each block's generic rank. Series
elimination must resolve every promised pivot, so a partial list cannot
be recorded as a complete profile. The catalogue retains the precision
used for every newly calculated block.

The [proof checker](../tools/check_coherent_q32_mixed_catalogue.py) rebuilds
every map, checks the exact retraction and deletion conditions, verifies
the complete profile lengths, and checks the source hashes of the
finished native computation. It uses no numerical pressure values in
the optimum proof.

## 3. Integer polynomial proof at every scale

For a candidate map, form its complete counting polynomials `Z_0(r)`
and `Z_infinity(r)`. Let `W_0,W_infinity` be those of the displayed
winner. Every physical height multiset is identical, so it is enough
to prove

\[
Z_0(r)Z_\infty(r)-W_0(r)W_\infty(r)\ge0.
\]

For each nonzero difference, the checker removes its exact power of `r`
and its exact power of `1-r`. The remaining integer polynomial has a
nonnegative coefficient expansion of one of two kinds:

| Certificate | Number of maps |
| --- | ---: |
| Nonnegative ordinary power coefficients | 326 |
| Nonnegative unnormalized Bernstein coefficients | 168 |
| Identical product polynomial | 2 |

For a Bernstein certificate of degree `d`, the checker proves the exact
identity

\[
Q(r)=\sum_{i=0}^d b_i r^i(1-r)^{d-i},\qquad b_i\ge0,
\]

by expanding the right side independently. Each nonzero certificate
has a positive coefficient, so it is strictly positive on `(0,1)`.
No root approximation or finite sampling of `r` enters this argument.

The [proof receipt](../evidence/coherent-q32-mixed-boundary-optimum.json)
retains every integer quotient and every nonnegative basis certificate.
Its two equality cases are exactly `(7,15)` and `(15,7)`.

## 4. Replay and limits

The default proof reads completed native profiles and checks all integer
identities. Recomputing every full native profile is an explicit optional
operation:

```sh
python3 tools/check_coherent_q32_mixed_catalogue.py \
  --output /tmp/coherent-q32-mixed-boundary-optimum.json
python3 tools/check_coherent_q32_mixed_catalogue.py --recompute-profiles \
  --output /tmp/coherent-q32-mixed-boundary-full-replay.json
```

The checked data cover 992 complete profiles, 1,015,808 Smith orders, and
253,952 physical deletion conditions. The finite catalogue now closes
the mixed-occupation version of this particular boundary search. A new
gain must change some part of that class, its coefficient field, or its
profile and type optimization.
