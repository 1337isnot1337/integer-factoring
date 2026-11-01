# Moving physical maps and the seven-letter core obstruction

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
