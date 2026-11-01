# An energy-refined Newton ceiling for arbitrary rational gauges

The zero-cut Smith profile of the retained five-letter crossing is not
fixed under gauges singular at zero. However, every rational joint-energy
gauge obeys an explicit stronger-than-usual Newton ceiling. The argument
applies over every characteristic-two field, without degree, polynomiality,
nilpotence, pole-budget or strict-endpoint assumptions on the gauge.

This identifies a genuine gap between the current zero profile and an
invariant ceiling. It supplies no gauge attaining that ceiling and no
factoring improvement. The complete score remains `1/2`.

## 1. A commuting diagonal reveals more invariants

Use the joint operator `A(t)` and energy decomposition from
[the single-evaluation algebra note](joint-evaluation-algebra.md).
Write

```math
R(t)=D_H(t)A(t),\qquad D_H(t)=\operatorname{diag}(t^{H_r}).
```

The zero orders used in the existing profile calculation are exactly
the local Smith orders of `R`. An energy gauge is a matrix preserving
each equal-energy block, so it commutes with `D_H`. For any invertible
rational such matrix `U`,

```math
R_U=D_H U A U^{-1}=U R U^{-1}.
```

Ordinary characteristic coefficients already constrain these Smith
orders. More information survives if a different scalar is placed on
each energy block. In the retained block order choose

```text
(14,1,6,14,6,15,5,8,3,12,13)
```

in the fixed `F16` encoding. Let `D` be this constant block-scalar
diagonal matrix. It commutes with every energy gauge, hence

```math
D R_U=U(D R)U^{-1},\qquad
P(z,t)=\det(I+z D R_U)=\det(I+z D R).
```

The constant matrix `D` is invertible, so multiplying by it changes no
Smith orders. Its coefficients can be adjoined to any original
characteristic-two field; constant field extension changes no local
minor valuations or ranks. The result therefore holds over every such
original field, not only those already containing `F16`.

## 2. Exact characteristic polynomial and prefix ceilings

Write `P(z,t)=sum_(k=0)^25 p_k(t) z^k`. Exact polynomial arithmetic gives
the following `t`-valuations, including `p_0=1`:

```text
0,8,20,28,40,56,68,88,100,120,132,152,168,
188,204,224,244,264,284,304,328,348,372,392,416,444.
```

If the increasing Smith orders of `R_U` are `lambda_1,...,lambda_25`,
then every `k` by `k` minor has valuation at least their first-`k` sum.
The coefficient `p_k` is a sum of principal minors of `D R_U`. Thus

```math
\sum_{i=1}^k\lambda_i\le\operatorname{ord}_t p_k,
\qquad \sum_{i=1}^{25}\lambda_i=444.
```

The Smith-prefix polygon is convex. It therefore lies below the lower
convex hull of the listed points. Indeed, on each hull segment,
convexity bounds a prefix polygon by the chord between its endpoint
values, which in turn lies below the hull chord. The hull vertices are

```text
(0,0),(1,8),(3,28),(4,40),(6,68),(10,132),
(14,204),(19,304),(23,392),(24,416),(25,444).
```

Its 25 slopes are

```text
8,10,10,12,14,14,16,16,16,16,18,18,18,18,
20,20,20,20,20,22,22,22,22,24,28.
```

These ceilings are stronger than the ordinary untwisted characteristic
polygon. For example, the prefix at four is at most 40 rather than 42;
the prefix at 24 is at most 416 rather than 418. The current zero profile
is still strictly less balanced at several other prefixes. In particular,
this invariant does not prove global optimality of the current profile.

## 3. Meaning for the search and verification

For every positive scale, convexity of the exponential makes the hull's
slopes an optimistic lower bound on the zero partition sum. Using those
slopes together with the *unchanged reference infinity profile* in the
current numerical pressure routine gives approximately `rho=0.693317562`.
That is an optimistic comparison, below `1/sqrt(2)`, rather than an
attained or certified crossing parameter. No argument here preserves the
infinity profile while changing the zero profile, or constructs a gauge
realizing the ceiling. Those are the substantive remaining tasks.

The checker constructs all 26 characteristic coefficient polynomials.
It uses the two support components of sizes 13 and 12, a literal subset
determinant expansion in characteristic two, and four independent binary
coefficient planes for `F16`. Its variable spacing exceeds the entire
possible `t` degree, so no monomial aliasing occurs.

Two separate controls check the computation: 256 direct `F16` determinant
evaluations use ordinary Gaussian elimination, and the top coefficient
agrees with the full determinant from a fraction-free Smith calculation.
The finite evaluations are controls; the symbolic coefficient computation
supports the uniform statement.

```sh
make verify-joint-graded-newton
```

Sources:
[`check_joint_graded_newton.py`](../tools/check_joint_graded_newton.py) and
[`joint-graded-newton.json`](../evidence/joint-graded-newton.json).
This is ordinary mathematics with exact algebraic replay, without a new
Lean proof or independent mathematical review. Substantial OpenAI
assistant support.
