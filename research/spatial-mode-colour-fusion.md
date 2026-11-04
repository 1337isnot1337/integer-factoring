# Fuse the leading spatial modes as well as the colour pairs

## Result and stopped fixed core

There is a uniform colour-and-spatial quartic construction whose two
fixed cut profiles begin, for every `q=2^s>=4`, with

```math
0^1,\qquad 6^{2(q-1)},\qquad 12^{5(q-1)}.                 \tag{1}
```

This changes the **whole** order-twelve multiplicity, rather than just
one selected column-sum restriction. At `q=8` it falls from 91 in the
[original fixed-curve tensor](fixed-curve-quartic-floor.md) to 35. The
endpoint, cold ranks, fixed parent wiring and homogeneous degree cost
are retained. The construction changes the colour-permutation symmetry
used by [colour alternation](colour-alternation-joint-cold-curve.md).

Returning layers nevertheless exclude this **fixed q8 core** from the
common-scale, unweighted paired-type target. Exact low jets and a
determinant-degree tail give a pressure floor `3/4>1/sqrt(2)`. Thus
neither (1) nor a local rank drop is a complete factoring improvement.
The complete proven exponent remains **1/2**.

The [prediction and stopping condition](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19363914)
were registered before computation; the
[checkpoint](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19363955)
records the successful prediction and returning layers. No larger
alphabet, curve, coefficient or full-profile sweep was launched.
The q8 certificate is not a theorem excluding every growing alphabet
in this new family, weighted gauges, or unequal scales.

## 1. One fixed bivariate tensor

Retain the PR97 Gaussian pair kernel, after the fixed substitution
`(X,Y,Z)=(x^2,y^2,xy)`, and the fixed pairing of z-neighbour vertices
into four parent superletters. Write this pair kernel as `K~`.
Each pair coefficient has degree sixteen. Its leading binary spatial
matrices at the two cuts are

```math
Q_x=[x^6y^{10}]\widetilde K,\qquad
Q_y=[x^{10}y^6]\widetilde K.
```

The cut submatrices have rank two. For example, in the x-cut order
`rows=(0,2,4,6)`, `columns=(7,5,3,1)`,

```
Qx = [1 1 0 0]
     [1 1 0 0]
     [0 0 1 1]
     [0 0 1 1].
```

For a fixed four-vertex support and assigned distinct nonzero colours
with XOR zero, define `alpha_x` as follows. It is zero unless the support
has two vertices in each x half. Otherwise, sort each half's two colours
in increasing numerical order. Match first colour to first colour and
second to second, and multiply the corresponding two entries of `Q_x`.
Define `alpha_y` in the same way at the fixed y cut. Add the creation
gate with coefficient

```math
C=\alpha_x x^{12}y^{20}
  +(\alpha_x+\alpha_y)x^{16}y^{16}
  +\alpha_y x^{20}y^{12}.                              \tag{2}
```

Include this gate for every such assigned support; no coefficient
parameters are searched. Numerical colour order is part of this one
uniform construction, so GL(s,2) colour transitivity is not assumed.
The two leading prescriptions are parts of one polynomial tensor,
without changing wiring or selecting a cut for an individual summand.

All creation gates commute; overlapping creations have product zero.
Every `C` has degree 32, valuation at least twelve at both axes, and
`C(1,1)=0`. The resulting full tensor therefore has:

- the same invertible antipodal endpoint, normalized by the same fixed
  internal z swap on the parent legs;
- cold rank one at each axis, since all added quartics vanish there;
- homogeneous degree **eight per occupied old vertex**, including all
  products of disjoint gates.

The complete local table has `q^8` entries, with maximum degree 64.
Pair/quartic expansion constructs it in `q^O(1)` operations at this fixed
eight-vertex size. No growing table or degree is treated as free.
The original q8 divided-square identity is not assumed for (2).

## 2. Proof of the all-alphabet leading prefix

Remove the creation gates supported wholly on one cut side. Each is
polynomially invertible, so this preserves Smith orders. Every remaining
pair has valuation at least six, and every remaining quartet at least
twelve. A term of occupied count `w` has valuation at least `3w`.
Through order twelve only the vacuum, one pair, or a four-occupied
entry can contribute. The vacuum is a disjoint unit; the one-particle
orders are `(6,6,18,18)` for each nonzero colour.

Fix a nonzero total colour charge `c`. Its unordered distinct nonzero
colour pairs `P={a,b}`, `a XOR b=c`, number `n=(q-2)/2`. Each pair has
twelve possible ordered placements on a four-vertex face. Canonical
colour order identifies these twelve states for every `P`.

Let `B` be the common twelve-by-twelve spatial matrix formed by matching
the two canonical colours with `Q_x`. It has rank four: it is the tensor
product of two rank-two spatial-mode pairings, and the distinct-vertex
placement maps still span all four ordered mode pairs. Each of the two
mode classes has two vertices, so equal-mode placement is also possible.

For `P=Q` the Gaussian leading block is `B`; no distinct-colour quartic
can change it. For `P!=Q`, equal nonzero XOR implies the pairs are
disjoint, and `alpha_x` sets the entire off-diagonal leading block to
the **same** `B`. Thus the full charge-c block at order twelve is

```math
J_n\otimes B,
```

of rank four for every `n>=1`, rather than four per pair. This argument
does not assume `n` odd. Charge zero has the unchanged same-colour
two-particle minors: each of the `q-1` colours contributes one rank-one
second-compound channel of the rank-two spatial kernel. Distinct-colour
quartics cannot connect these repeated-colour states.

A quartet split one/three has `alpha_x=0`, so it creates no additional
order-twelve entries between one- and three-occupied states. Larger
occupations start at eighteen. Consequently the whole order-twelve
multiplicity is `(q-1)+4(q-1)=5(q-1)`, proving (1). The y-cut proof
is identical, with `alpha_y`; the intervening monomial in (2) starts
at order sixteen and cannot affect this argument.

## 3. The fixed q8 return and complete pressure invoice

The [checker](../tools/check_spatial_mode_colour_fusion.py) computes all
eight charge blocks, each with 512 rows including repeated colours and
all occupations. No colour-orbit reduction is used. Both cuts have the
exact prefix through eighteen

```math
0^1,\quad6^{14},\quad12^{35},\quad15^{84},\quad
16^{70},\quad17^{14},\quad18^{208}.                     \tag{3}
```

Independent column ranks at truncation nineteen check the streamed row
ranks. Literal quartet-support/Pfaffian enumeration checks 213,280
entries per cut against the matching formula, using all seven nonzero
zero-XOR colour quartets and their extensions by an equal-colour pair.
Two quartics start at order 24 and cannot enter these jets; this is a
valuation argument, not the old divided-square shortcut.

There are `n=4096` finite Smith orders because the endpoint is invertible.
The homogeneous full determinant degree is still `56n`, by the charged
degree calculation in the original fixed-curve note. The 426 known
orders in (3) sum to 6866. The remaining 3670 orders have mean at most
`22251/367<61`. Convexity and monotonicity therefore give

```math
Z(t)\ge1+14t^6+35t^{12}+84t^{15}+70t^{16}
        +14t^{17}+208t^{18}+3670t^{61}=:L(t).
```

The difference `L(t)-(1+7t^(32/3))^3`, after `t=u^3` and removal of
the positive factor `u^18`, has all 301 degree-300 Bernstein coefficients
strictly positive. The checker uses exact rationals and reconstructs
every monomial coefficient, including zeros. Hence
`Z(t)>(1+7t^(32/3))^3` throughout `0<t<1` at each cut.

For a common paired-letter type `nu`, with entropy `Hnu>0` and mean
old-letter occupancy `Ew`, the charged two-cut envelope is

```math
\frac{(\log Z_x+\log Z_y)/2-16Ew\log t}{2Hnu}.
```

The paired-letter Gibbs partition is `(1+7t^(32/3))^2`. It proves a
floor at least `3/4`, for arbitrary correlated paired types. This is
an obstruction to this common-scale profile estimate, not an actual
coefficient-rank or factoring lower bound. Expanding the fixed q8
profile cannot repair this certified failure.

## Verification, credit and next structural question

`make verify-spatial-mode-fusion` reproduces the
[source-pinned receipt](../evidence/spatial-mode-colour-fusion.json).
These are author-checked ordinary arguments and finite exact controls;
no new Lean or independent mathematical review is claimed.

The [full default receipt](../evidence/spatial-full-verification.json) and
[complete compressed stdout](../evidence/spatial-full-verification.stdout.txt.xz)
record exit-zero `make verify` at clean frozen
`a5189cd125c5d651ac8235b2d4fe19654dd459e4`, tree
`c9e4bbbde40cce049c0484763e6ba61ef3b6da2e`, in 1406.23 seconds.
A fresh marker probe verified private `/tmp`, read-only host root and
writable frozen checkout; network was unshared. Tracked state remained
clean. PR110/113 and the later sampler-budget and archival prose are
outside this source pin. The latter prose retains focused documentation
and whitespace checks.

The pair construction retains the original homogeneous-Gaussian and
worker-4 colour-class attribution. The creation-operator model comes
from the primary worker's PR97. The fixed curve, paired-degree invoice,
Toeplitz/column checks, literal-control scaffolding and convex-tail/
Bernstein method retain Sol worker 2's PR111 credit. Sol worker 1
developed (2), its whole-leading-block proof and this fixed-core test
with substantial OpenAI Codex assistance. No novelty claim is made for
tensor products, canonical ordering, Smith truncations or Bernstein bases.

The structural gain is the uniform linear-in-q leading layer, not a
complete bound. A further proposal must control the returning layers
or change the transfer with all costs included; simply varying more
coefficients or enlarging this fixed-core replay is not the next step.
