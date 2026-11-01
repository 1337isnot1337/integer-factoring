# A uniform natural-pressure floor for height-25 boundary retractions

## Result and precise scope

Fix the 92 coordinate occupation words of height below 25. Add four
paired labels from the eleven-dimensional height-25 space, using any
height-preserving embedding and retraction with LP=I. The two physical
sides may use different maps. Over every characteristic-two field,
the natural zero/infinity typed Smith-pressure cost satisfies

```math
\boxed{\rho\ge\frac{17}{24}>\frac1{\sqrt2}.}             \tag{1}
```

The boundary maps may be nonbinary, nonorthogonal, and may mix the
three- and five-particle words at height 25. This is not a finite
catalogue statement. It includes the
[improved F4 hyperoval code](f4-hyperoval-crossing-bound.md).

The conclusion concerns this natural two-place pressure method with
the lower 92-label core retained. It does not exclude changes to that
core, more general physical weights, other valuation places, or better
coefficient recombinations. In particular, it is not a lower bound on
integer factoring. The complete factoring exponent remains 1/2.

## 1. A common polynomial submatrix

Let

```math
h(w)=\sum_{a\in w}(2a+1),\qquad
\mathcal C=\{w:h(w)<25\}.
```

There are 92 such words, all using the first twelve hot modes. Every
candidate in the stated class retains their coordinate embeddings
and retractions. In both regular natural-cut compound matrices, the
rows and columns indexed by pairs from C are therefore exactly the
same 8,464-by-8,464 coordinate-core matrix.

For this core the embedding equals its retraction. The two compound
formulas are literally identical, so its two cut profiles coincide.
Write their common counting polynomial as

```math
Z(r)=\sum_{i=1}^{8464}r^{\nu_i}.
```

The [complete core receipt](../evidence/hot-height25-common-core.json)
records all orders. Its coefficients are binary, so extension of the
coefficient field does not change any minor valuation or Smith order.
The calculation resolves every endpoint-proved pivot, including orders
above the initial series window; no unresolved tail is discarded.

All four added physical labels have height 25. Thus every candidate
has the same physical height polynomial

```math
Q(r)=\sum_{w\in\mathcal C}r^{h(w)}+4r^{25}.            \tag{2}
```

Height homogeneity gives the coherent natural-cut factorization even
when a height-25 vector mixes occupations. Deletion closure is automatic:
every deletion of a height-25 word has smaller height and is already
in C. This is the [mixed-occupation interface](coherent-hot-mixed-boundary.md),
with the same linear identities over characteristic-two extensions.

## 2. Submatrix monotonicity of the full channel sum

For a polynomial matrix M with nonzero Smith orders nu_i, its truncated
multiplication operator modulo t to the power ell has rank

```math
R_M(\ell)=\sum_i\max(\ell-\nu_i,0).
```

If N is a polynomial submatrix of M, its truncated operator is a
submatrix too. Hence R_N(ell) is at most R_M(ell). For 0<r<1,

```math
\sum_i r^{\nu_i}
 =(1-r)^2\sum_{\ell\ge1}R_M(\ell)r^{\ell-1}.           \tag{3}
```

This follows by summing the geometric series separately for each
Smith order. The series converges because the rank grows at most
linearly in ell. Combining rank monotonicity and (3) proves that a
submatrix's channel sum is at most the full channel sum.

Apply this to the fixed core in each candidate cut. If its two sums
are Z_0 and Z_infinity, then

```math
Z_0(r)\ge Z(r),\qquad Z_\infty(r)\ge Z(r).             \tag{4}
```

The opposite scalar monomial shifts in the natural reflection cuts
cancel. The retained vacuum entry also prevents an additional common
nonmonomial scalar factor. Neither a partial profile nor a conjectured
interlacing of individual Smith orders is used in (4).

## 3. One exact inequality for every scale

The certificate proves the global polynomial inequality

```math
Z(z^{17})^{24}>Q(z^{24})^{34}\qquad(0<z<1).           \tag{5}
```

Here is a short verification argument that avoids expanding either
large power. Set

```math
F(z)=24\log Z(z^{17})-34\log Q(z^{24}).
```

Its derivative has the sign of the integer polynomial

```math
G(z)=Z'(z^{17})Q(z^{24})-2z^7Q'(z^{24})Z(z^{17}).      \tag{6}
```

The degree is 2,827 and G(0)=1. The checker computes every coefficient
of (6). For d=2827, the coefficients of

```math
(1+u)^d G\!\left(\frac{u}{1+u}\right)
```

are exactly

```math
b_k=\sum_{i=0}^k g_i\binom{d-i}{k-i}.
```

Their nonzero signs have exactly two variations. Descartes' rule
therefore gives at most two roots of G in (0,1), counted with
multiplicity. Exact integer evaluation gives opposite signs at

```math
l=\frac{98447622}{10^8},\qquad
u=\frac{98447623}{10^8},\qquad G(l)<0<G(u).
```

Together with G(0)>0, these signs force exactly two roots. The second
lies in [l,u], and F first increases, then decreases, then increases.
Its only interior local minimum is in that small interval.

Positive coefficients make Z and Q increasing. Outward integer power
bounds at 192 bits establish

```math
Z(l^{17})^{24}>Q(u^{24})^{34}.                         \tag{7}
```

The checker verifies (7) directly by an integer power comparison, and
also by independent outward logarithm bounds. The latter give lower
numerator

```text
1121601824999796586546104499855224182287520348184971837386
```

in units of 2 to the power minus 192 for F on the entire interval.
Thus its local minimum is positive. Since F(0)=0 and its initial
increase is strict, F(z)>0 throughout (0,1). This proves (5).
The endpoint fractions and all signs are checked exactly; the proof
does not depend on a numerical root finder or optimizer.

## 4. Entropy gives the pressure floor for every physical type

Let s>0 and substitute z=exp(-s/17) in (5). With rho_0=17/24,

```math
\log Z(e^{-s})>2\rho_0\log Q(e^{-s/\rho_0}).           \tag{8}
```

For any paired physical distributions p_X,p_Y with positive total
entropy, the natural pressure cost has numerator

```math
\log Z_0(e^{-s})+\log Z_\infty(e^{-s})
 +2s\bigl(\langle h,p_X\rangle+\langle h,p_Y\rangle\bigr)
```

and denominator 2(H(p_X)+H(p_Y)). The Gibbs entropy inequality applied
to (2) gives

```math
H(p)\le\log Q(e^{-s/\rho_0})+\frac{s}{\rho_0}\langle h,p\rangle.
```

Use this on both sides, then combine (4) and (8). The numerator is
strictly larger than rho_0 times the denominator. Taking infima gives
(1). This includes types with some zero probabilities; no interior
Gibbs-optimum assumption is needed. Finally, 2 times 17 squared is
578, strictly larger than 24 squared, 576, so the floor exceeds the
post-half threshold exactly.

## 5. A limited weighted corollary

The same conclusion holds for additional weights supported only on
the four boundary labels, with all lower-core weights zero, when each
cut satisfies the [deletion-operator integrality cone](weighted-coherent-hot-f4.md).
Every nonzero positive-height boundary vector has a nonzero deletion.
Pair that deletion edge with the retained mode-zero edge from word 1
to the vacuum on the other leg. Its kernel valuation is zero because
a bitwise-and zero is zero. The corresponding cone inequality forces
each boundary weight to be nonnegative.

The weighted core is consequently unchanged, and all additional
physical means are nonnegative. The argument above still applies.
This corollary does not permit arbitrary interior weights or weights
outside the sufficient cone.

## 6. Reproduction and research consequence

Run `python3 tools/check_hot_boundary_core_barrier.py`. The
[exact certificate](../evidence/hot-boundary-core-barrier.json) retains
the derivative coefficients, transformed sign sequence, rational
minimum box, integer comparison fingerprints and source bindings.
The separate `tools/compute_hot_common_core.py` reconstructs the full
core profile; its recorded run took about 516 seconds.

Boundary-only variations can still yield small improvements such as
the F4 hyperoval result. They cannot pass the post-half threshold by
this natural-pressure route. The substantive next targets are changes
to lower-height structure, different polynomial templates, or methods
outside the scope stated above.

Substantial OpenAI assistant support was used for the argument,
computation, verification, and text.
