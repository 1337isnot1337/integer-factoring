# Moving physical maps and protected-core obstructions

## Construction and research outcome

Polynomial physical embeddings and retractions give a new effective
crossing family in which both represented spaces can vary with the
crossing parameter. This escapes the fixed-embedding hypothesis of the
[joint rational row model](joint-hot-rational-retraction-module.md).
It does not, by itself, escape a protected physical subtensor.

For the first eight-letter family below, an ordinary proof and an exact
certificate establish the uniform natural-pressure floor

```math
\rho\ge\frac57>\frac1{\sqrt2}.
```

The theorem allows arbitrary polynomial degree and arbitrary changes
outside the specified seven-letter subtensor. It is stronger than a
failure of the finite screen, but concerns the stated common-scale,
natural physical gauges. Other gauges, changes to the lower subtensor,
and different coefficient-rank arguments remain open. The complete
factoring score stays **1/2**; no new Lean work is included.

The subsequent [five-letter certificate](#five-protected-letters-already-suffice)
strengthens the protected-subtensor scope. A separate finite-prefix
certificate covers 16,384 completed frames at every common natural scale
and physical type; it does not exhaust the 65,536-case search.

This work starts from repository commit `941cc1d` and builds on the
[coherent boundary construction](coherent-hot-boundary-crossings.md),
the [protected-submatrix argument](coherent-boundary-protected-floor.md),
and the original reflection constructor cited there. Substantial OpenAI
assistant support was used for the mathematics, implementation, and checks.

## Polynomial paired maps

Over any characteristic-two field, choose polynomials a(w), b(w), and set

```math
L(w)=(1,a(w)),\qquad
P(w)=(1+a(w)b(w),b(w))^{\mathsf T}.
```

Then LP=1 identically. Extend both maps by the identity on retained
coordinate states and apply L,L,P-transposed,P-transposed to the four
physical legs of a polynomial crossing. At w=1 the endpoint contracts to
`(LP)_(x,v)(LP)_(y,u)`, so it is exactly the new crossing. There is no
evaluation of an undefined internal rational expression.

If the ambient degree is D, the resulting degree is at most
`D+4 deg(a)+2 deg(b)`, interpreting the degree of a zero parameter as zero
for this upper bound. The maps and tensor are explicit over one fixed
finite field. Tensor powers of a fixed choice still have linear degree
and polynomially many coefficient pieces. A useful complete pressure or
rank theorem is still necessary before claiming an improved evaluator.

The [constructor](../tools/moving_boundary_crossings.py) acts in the
original w variable on the literal polynomial tensor. It does not use
the constant-map compound-profile shortcut for varying maps.

## Complete finite affine screen

Use the ambient nine-letter reflection crossing at interval 16, on hot
occupation words `0,1,2,3,4,5,8,6,9`. The last two words have equal height
eight. Retain the first seven coordinates and compress this last pair to
one letter. The inherited physical heights are

```text
0, 1, 3, 4, 5, 6, 7, 8.
```

All **256** choices of affine a,b over F4 were constructed. Both complete
64-order natural Smith lists were resolved for every tensor. There are
256 distinct tensors and 23 pairs of profiles. The largest resulting
raw degree is 164. The best retained numerical pressure is approximately
`0.719975317855538`, attained by constant a=b=1. These fits are numerical
search evidence, without a global optimization guarantee or an exact
pressure certificate. They do not improve the prior constant construction.

The [complete compressed archive](../evidence/moving-boundary-f4-affine.json.gz)
records every parameter pair, tensor digest, both profiles, source hashes,
and fit result. A quadratic screen retained **2,112 of 4,096** planned
pairs before the uniform obstruction below justified stopping it. Its
[partial archive](../evidence/moving-boundary-f4-quadratic-partial.json.gz)
explicitly records that it was not exhausted. No finite negative theorem
is inferred from that partial run.

## A uniform obstruction independent of the maps

Let C be the seven coordinate letters of heights below eight. Consider
any polynomial eight-letter crossing T(w) that agrees with the original
reflection tensor whenever all four physical letters lie in C. Changes
elsewhere are unrestricted. In particular, the theorem covers every
polynomial pair a,b above, without a degree bound.

Use the natural weighted cut matrices: at zero the row weight is
`2(h(x)+h(y))`; at infinity the weight is `2(h(u)+h(v))`. Their common
physical charge is the inherited height h. Both matrices contain the same
fixed 49-by-49 weighted coordinate submatrix. Its two profiles, in w,
are `128+2 nu` and `-128+2 nu`, where the complete nu histogram is

```text
order:         0 1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22 24 28
multiplicity:  1 1 1 1 1 1 1 2 2  2  2  3  3  3  4  3  3  4  4  2  2  1  2
```

Write Z(r)=sum r^nu with these multiplicities, and
`Q(r)=1+r+r^3+r^4+r^5+r^6+r^7+r^8`.
The [Smith counting-function monotonicity proof](coherent-boundary-protected-floor.md#2-smith-counting-functions-increase-under-taking-a-larger-matrix)
applies to both fixed submatrices. It also applies to Laurent matrices
after a common scalar shift of the matrix and submatrix. Thus, at any
scale s>0, the product of the two complete channel sums is at least
`Z(exp(-2s))^2`; the opposite 128 shifts cancel. This step requires no
deletion invariance of the changed tensor or its maps.

The exact all-scale inequality needed for the floor is

```math
Z(z^5)^7>Q(z^7)^{10}\qquad(0<z<1).                 \tag{1}
```

For completeness, put `F(z)=7 log Z(z^5)-10 log Q(z^7)`. Its derivative
has the sign of the integer polynomial

```math
G(z)=Z'(z^5)Q(z^7)-2z^2 Q'(z^7)Z(z^5).
```

G has degree 191 and G(0)=1. The coefficients of
`(1+u)^191 G(u/(1+u))` have exactly two sign variations. Descartes' rule
therefore gives at most two roots in (0,1), counted with multiplicity.
For the exact fractions

```text
l = 90000111/100000000,
u = 90000112/100000000,
```

the checker proves G(l)<0<G(u). Together with G(0)>0 this forces two
roots, with the sole interior local minimum of F in [l,u]. Positive
coefficients make Z and Q increasing, and an exact integer comparison
after denominator clearing proves `Z(l^5)^7>Q(u^7)^10`. This proves F>0
at that minimum. Since F(0)=0 and F initially increases, (1) follows
throughout the interval. All decisions use integers; no floating root
finder enters the certificate.

Substitute `z=exp(-tau/5)` in (1). It gives

```math
\log Z(e^{-\tau})>
2(5/7)\log Q(e^{-\tau/(5/7)}).
```

Combine this with the two protected channel sums and the Gibbs entropy
inequality for any paired physical type. The resulting natural pressure
is at least 5/7, including types with zero probabilities and independently
chosen types on the two paired alphabets. Finally `2*5^2-7^2=1>0`, so
this floor exceeds the required threshold exactly.

The proof protects a literal subtensor; varying only the boundary with
higher degrees cannot remove it. The constructive next search must
change lower states or change the pressure method.

## Verification

The [contraction checker](../tools/check_moving_boundary_crossings.py)
compares 65,536 constant-specialization entries with the existing
constructor, 12,288 full moving-map polynomial entries with independent
coefficientwise F4 convolution, and 500 arithmetic/splitting controls.
It separately checks all 12,288 endpoint entries of the moving examples
and resolves both complete profiles of each.

The [floor checker](../tools/check_moving_boundary_core_floor.py) reconstructs
the protected core and compares both 49-order profiles with independent
binary-polynomial Smith elimination and the native F4 implementation.
Its [certificate](../evidence/moving-boundary-core-floor.json) includes
the 98 order comparisons, every derivative coefficient and transformed
sign, the exact root box, and integer comparison hashes. The
[contraction receipt](../evidence/moving-boundary-controls.json) records
the separate construction checks and source bindings.

```sh
make verify-moving-boundary
```

The longer finite screen can be replayed independently:

```sh
python3 tools/search_moving_boundary_crossings.py \
  --field 4 --degree 1 --output /tmp/moving-boundary-f4-affine.json
```

The mathematical obstruction is an ordinary proof supported by exact
finite certificates. It is not a stronger factoring bound, a Lean theorem,
or an independent mathematical review.

The [exact-source verification receipt](../evidence/moving-map-verification/receipt.json)
records a completed default `make verify` at
`410e3e2cda3b416dffe358be67fc8e1413ee66a9`: 341 Python command blocks,
66 unit tests, five crossing tests, and links in 973 Markdown files. The
frozen worktree remained clean. The
[complete compressed log](../evidence/moving-map-verification/full-410e3e2.log.gz)
is retained. The same receipt records all four focused moving-map checks
at the later frame source `1dba7bad8628b7ce3b3448da76bb79361c316cb5`, with
their [focused log](../evidence/moving-map-verification/focused-1dba7ba.log).
Later integrations and the separate full frame search are outside the
default-gate source pin. No Lean replay is claimed by these ordinary checks.

## Moving lower states and complete interior-pair controls

The boundary obstruction motivates a different search: move any pair of
the nine physical states, including the vacuum and unequal-height pairs.
The [interior driver](../tools/search_moving_interior_crossings.py) relabels
the literal ambient tensor before applying the polynomial maps. It
computes both complete weighted cuts directly. There is no assumption that
the maps preserve height, occupation number, or a Gaussian factorization.

A complete binary affine screen covers all 36 physical pairs and all 16
a,b choices per pair. For unequal-height pairs the assigned physical
weight is tested at both original heights and their midpoint. The unique
equal-height pair needs only one weight. This gives exactly **1,696**
full two-profile evaluations. The
[complete archive](../evidence/moving-interior-binary-affine.json.gz) records
the declared parameters and every complete profile. Its best numerical
pressure remains approximately `0.7199753178555385`, attained by the old
constant map on the height-eight pair. No all-family theorem follows from
these finite fits.

For low-entropy Gibbs roots, the implementation first cancels exact scalar
shifts of the two spectra and the physical weights, then evaluates the
remaining positive sums with `log1p`. This avoids losing small positive
terms to floating cancellation. The original unnormalized numerical
routine failed on a pilot launch; that launch is not counted as the
completed screen. The successful screen was rerun with the normalized
routine. Its [independent controls](../evidence/moving-interior-controls.json)
check 19,683 literal relabeling entries, 12,288 complete polynomial entries,
and root invariance under opposite cut shifts of 100,000 and a common
physical shift of 300. These numerical root checks are separate from the
exact tensor and profile checks.

## Codimension-one frames with no protected lower coordinate core

For q-vectors a(w),b(w), define polynomial physical maps

```math
L=[I_q,a],\qquad
P=\begin{bmatrix}I_q+a b^{\mathsf T}\\b^{\mathsf T}\end{bmatrix}.
```

Characteristic two gives `LP=I_q`. Applying L,L,P-transposed,P-transposed
therefore supplies the exact q-letter endpoint as before. Unlike a single
boundary cell, these maps can change all retained physical states. The
seven-letter protected-corner theorem applies only if the resulting
tensor actually retains that corner.

The [frame constructor](../tools/moving_physical_frames.py) applies these
maps by rank-one updates, without materializing dense polynomial maps.
On a left leg its ith output is `input_i+a_i input_last`. On a right leg
its jth output is

```math
input_j+b_j\left(input_{last}+\sum_i a_i input_i\right).
```

These identities are purely polynomial. No matrix rank estimate is
inferred from the low-rank form of a physical map.

### Reciprocal Laurent flows

The maps also admit monomial denominators. Write
`a_i=w^(s_i) A_i(w)` and `b_j=w^(t_j) B_j(w)`, where A_i,B_j are
polynomials. Choose

```math
m_L=\max(0,\max_{A_i\ne0}(-s_i)),\qquad
m_P=\max(0,\max_{B_j\ne0}(m_L-t_j)).
```

Then `w^m_L L` and `w^m_P P` are polynomial maps, and their product is
`w^(m_L+m_P) I`. Applying the four cleared maps multiplies the Laurent
crossing by exactly `w^(2(m_L+m_P))`. At w=1 the endpoint is unchanged.
This common monomial only reindexes the finitely many coefficient
tensors; it does not change their ranks. In the two-place counting
functions it adds a common shift to one spectrum and the opposite shift
to the other, so their product is unchanged. The assigned physical
weights and their full entropy cost are still included.

Every fixed set of Laurent shifts therefore yields an effective
polynomial crossing after clearing. This argument removes only a common
monomial. It does not discard a nonmonomial denominator or an uncharged
physical cost.

The [dense checker](../tools/check_moving_physical_frames.py) independently
constructs the cleared maps, verifies their exact scaled split identity,
and contracts all four legs by coefficientwise F4 convolution. It checks
24 dense frames with positive and negative shifts on three ambient
letters, comprising 384 complete entries. Three reflection specializations
compare all 12,288 entries with the separate boundary constructor. The
[receipt](../evidence/moving-physical-frame-controls.json) distinguishes
these dense small controls from full eight-letter endpoint checks.

### Declared pilot cells and next search

The [pilot archive](../evidence/moving-physical-frame-pilots.json.gz) contains
nine complete 16-case cells. In each cell, A and B range independently over
the masks `0,1,128,255`. Uniform flows use `a=w^c A,b=w^(-c) B` for
`c=-2,-1,0,1,2`. Height flows instead use

```math
a_i=A_iw^{c(8-h_i)},\qquad b_i=B_iw^{-c(8-h_i)},
\qquad c=-2,-1,1,2,
```

on physical heights `0,1,3,4,5,6,7,8`. All 144 candidate tensors have
complete cut profiles. Their best numerical pressure remains the
constant boundary control; these small cells are not an exclusion of
the full frame family.

The [full driver](../tools/search_moving_physical_frames.py) exposes both
flow families and an affine polynomial family. Without restricted mask
lists it covers all 65,536 binary pairs in a deterministic recorded order.
It writes each completed case to a separate JSON-lines record and saves
progress without overwriting another run. Source hashes, both profiles,
tensor digests, and the clearing power are retained. A live partial run
must not be described as an exhausted search.

```sh
python3 tools/search_moving_physical_frames.py --height-flow 1 \
  --output /tmp/moving-physical-frames.json \
  --records /tmp/moving-physical-frames.jsonl
```

A numerical pressure below `1/sqrt(2)` would trigger an exact candidate
audit and uniform-transfer accounting. None is established by the
completed controls above. The next positive task is a realizing frame,
tensor, or coefficient decomposition with a complete strict bound.

## Only the endpoint split is necessary

The all-parameter identity LP=I is a sufficient construction rule, not a
requirement of the crossing transfer. One may instead use

```math
L=[I,a],\qquad P=\begin{bmatrix}I+r b^{\mathsf T}\\b^{\mathsf T}\end{bmatrix},
\qquad LP=I+(a+r)b^{\mathsf T},
```

with independent polynomial vectors a,r, provided `a(1)=r(1)`.
The endpoint is still exactly the crossing after all four physical maps.
The [constructor](../tools/moving_physical_frames.py) supports this case
and continues to check every complete endpoint entry.

An explicit binary family is

```math
a=e_7+wA,\qquad r=w(A+e_7),\qquad b=B.
```

It has `LP=I+(1+w)e_7B^T`. If `B_7=1`, its determinant is w, so it
loses one rank at zero while remaining the identity at one. This
intermediate degeneration was disallowed by the earlier all-parameter
split constraint. It is not itself a favorable crossing-rank estimate.

The [endpoint-only pilot](../tools/search_endpoint_physical_frames.py)
covers all 100 pairs in the declared ten-by-ten mask grid. Its
[complete archive](../evidence/endpoint-moving-frame-pilot.json.gz) retains
both complete profiles for every tensor. The best numerical pressure is
approximately `0.7199780715909008`, so this pilot gives no improved bound.
Independent dense controls include 24 separate endpoint-only frames and
384 complete polynomial entries. Four full reflection controls compare
all 512 natural Smith orders with independent binary-polynomial
elimination, including an endpoint-only dense frame.

## Five protected letters already suffice

The protected-core obstruction extends from seven letters to just the
five original hot words `0,1,2,3,4`, of heights `0,1,3,4,5`.
Allow arbitrary polynomial changes to every other entry of an
eight-letter tensor, keeping the physical heights
`0,1,3,4,5,6,7,8` and the crossing endpoint. If the tensor retains this
five-letter reflection subtensor, its natural common-scale pressure is
still at least **5/7**, over every characteristic-two field, at any
polynomial degree and for all physical types. This includes independent
types on the two paired alphabets. Other physical gauges and better
coefficient-rank methods remain outside the theorem.

The common core has the following 25 complete orders in the deflated
variable:

```text
0,1,3,4,5,6,7,8,8,9,9,10,10,11,11,12,12,12,13,13,14,15,15,16,20.
```

Both cuts are independently checked against native F4 profiles. With Z
the counting polynomial of this list and Q the unchanged eight-height
polynomial, the proof again establishes `Z(z^5)^7>Q(z^7)^10` for
every `0<z<1`. Its derivative-sign polynomial has degree 151 and exactly
two Descartes sign variations after the same interval transformation.
The sole interior minimum is enclosed by
`91056849/10^8` and `91056850/10^8`; the signs and the strict lower
power comparison are integer calculations. The
[checker](../tools/check_moving_five_core_floor.py) and
[certificate](../evidence/moving-five-core-floor.json) retain all 50
independent order comparisons and the complete all-scale data.

The submatrix and Gibbs entropy arguments are identical to the
seven-letter proof above. The smaller protected set broadens the
excluded class: moving only higher states cannot suffice by this
natural-pressure route. A positive construction must change the stated
low subtensor or use information outside that pressure criterion.

## An exact quarter-prefix exclusion

The source-frozen height-flow-one run has an immutable
[16,384-case archive](../evidence/moving-frame-height1-prefix16384.jsonl.gz)
and [source/order manifest](../evidence/moving-frame-height1-prefix16384.json).
It covers the first 64 declared A masks and all 256 B masks at source
`1dba7ba`. This is a quarter of the planned finite search; the remaining
cases are not claimed covered by this snapshot.

For every retained frame, let P(r) be the product of its two complete
natural profile counting functions, using the original w orders. Let
C(r)=Z_7(r^2)^2 for the previously certified seven-letter core.
The [prefix checker](../tools/check_moving_frame_prefix.py) proves
`P(r)>C(r)` for all `0<r<1` by exact integer coefficients.

Specifically, write `P-C=sum_k d_k r^k`, allowing negative integer
exponents, and let `S_k=sum_(j<=k)d_j`. Every S_k is nonnegative in
every one of the 16,384 cases, and the final cumulative excess is
`64^2-49^2=1695`. For smallest/largest exponents m,M,

```math
P(r)-C(r)=(1-r)\sum_{k=m}^{M-1}S_k r^k+S_Mr^M>0.
```

The inherited core theorem therefore gives natural pressure at least
5/7 for **every retained frame, every common positive scale, and every
physical type**. No numerical fit is used in this finite exclusion.
The [replay receipt](../evidence/moving-frame-height1-prefix16384-floor.json)
binds the complete integer profiles, their exact source and parameter
order, and the core certificate. It does not exclude other flows,
nonbinary coefficients, other physical gauges, or the unfinished cases.

## Bounded full-gauge refinements

The [refinement driver](../tools/refine_moving_physical_frames.py) tested
12 declared height-flow frames with three full-gauge starts each, at an
80-iteration budget. Its [archive](../evidence/moving-frame-gauge-refinements.json.gz)
retains all 36 attempts. One optimizer reported convergence and 35 reached
their iteration limit. There were no arithmetic failures; the 7,826
objective evaluations did not produce a pressure below the old constant
control. The best retained numerical value is approximately
`0.7199753178555406`.

Those iteration-limited fits are not optimization exclusions or global
optimality certificates. They remain distinct from the exact natural
quarter-prefix obstruction and the all-degree protected-subtensor proofs.

## Asymmetric supporting-line screen of the frozen prefix

The balanced product obstruction does not exclude asymmetric supporting
weights. The [existing supporting-line theorem](asymmetric-support-line-crossings.md)
and [complete rectangular transfer](asymmetric-modular-transfer.md) give a
different sufficient criterion: positive weights a,b with `a^2+b^2=1`
and a certified support value R<1 would yield a stronger exponent.
The coefficient Chernoff slopes must be `t/a` and `t/b`, so the coefficient
index cancels. Independently choosing those slopes is invalid.

For the retained natural profiles, set `E_i=4h_i`. Both crossing-pair
physical charges are `t E_i`. After eliminating both physical types, the
rate-one test is

```math
a\log\sum_j e^{-t\nu_{\infty,j}/a}
+b\log\sum_j e^{-t\nu_{0,j}/b}
-2\log\sum_i e^{-tE_i}<0.                         \tag{2}
```

The [data-only screen](../tools/screen_asymmetric_moving_profiles.py)
uses every complete profile in the 16,384-case frozen prefix. It removes
opposite scalar shifts exactly before evaluating the logarithms.
Six unsymmetric native-oracle comparisons and two balanced comparisons
verify its normalization: at a=b=1/sqrt(2) and t=s/sqrt(2), its support
value equals sqrt(2) times the ordinary balanced rho. The maximum
discrepancy in those eight controls is `2.23e-15`.

The [complete screen archive](../evidence/asymmetric-moving-frame-prefix-screen.json.gz)
records **11,010,048 numerical grid evaluations**, with 13 declared
rational half-angle choices, the balanced angle, and 48 tilts between
0.002 and 1.5. None gives a negative rate-one gap. Twenty retained cases
were then refined in angle and tilt, within the stated bounds; all
twenty optimizers reported convergence. The smallest retained support
value is approximately **1.0004187269519174**, above the required one.

This is a finite numerical search, not an all-angle or all-scale theorem.
It does not extend the exact quarter-prefix exclusion to arbitrary
supporting weights, and it does not establish a new factoring exponent.

```sh
python3 tools/screen_asymmetric_moving_profiles.py \
  --output /tmp/asymmetric-moving-frame-prefix-screen.json
```
