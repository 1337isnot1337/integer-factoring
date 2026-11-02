# Complete norm search near a known smooth prime power

## Result and exact limit

A rejected smooth leading norm makes the entire smaller-norm search cheap
when the target has a known smooth representation close to `sqrt(D)`.
This is an elementary geometry-of-numbers deduction for the actual
prime-power targets, rather than an assumed coordinate-localization theorem.
It does **not** prove that the search finds a relation.

Let `Q=(A,B,C)` be a primitive reduced positive form of discriminant
`-D`, and suppose it represents a `Y`-smooth integer `m` coprime to `D`.
If its leading coefficient is not both `Y`-smooth and coprime to `D`, then

```math
A\ge D/(4m).                                                   \tag{1}
```

All integer vectors with `0<Q(x,y)<m` can then be enumerated in

```math
O(1+m/\sqrt D)\,\operatorname{poly}(\log D+\log m)
                                                               \tag{2}
```

bit operations, before smoothness testing. Testing and decoding every
primitive eligible vector costs an additional factor polynomial in `Y`
and the written bit lengths. Rejected vectors and row setup are included.

For `Q=red(g_p^k)`, the primitive norm `m=p^k` is explicitly available.
A smooth primitive norm **strictly below** `m` yields a weighted dominant
relation. Thus (2) removes the coordinate-search obligation from a
near-height strict-descent hypothesis for these targets. The strict-descent
existence statement remains open. A counterfamily below shows that even
the first power above the reduced-height bound can have no smaller smooth
norm. No complete factoring improvement follows; the score stays **1/2**.

The [scope and stopping test](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19359593)
were registered before the fixed controls. This changes the arbitrary-target
assumptions in [smooth-norm localization](smooth-norm-localization.md): here
the target already has a known smooth norm, and the searched norm is close
to `sqrt(D)`, rather than of order `D`. The
[target-dependent norm-selector qualifications](norm-selector-coverage.md)
and the [principal-word gap](leading-power-window-graphs.md#1-the-graph-and-the-lattice)
remain applicable. Neither implies descent existence.

## 1. A failed leading test forces balance

Write the known representation as `m=Q(x_0,y_0)`. If `y_0=0`, then
`A` divides `m`, so it is smooth and coprime to `D`, contrary to the
failed leading test. Thus `y_0!=0`. Completing the square gives

```math
4Am=(2Ax_0+By_0)^2+Dy_0^2\ge D.
```

This proves (1). The balance deduction does not require the known
representation to be primitive; the prime-power application supplies a
primitive one anyway. The leading test must precede the claim: for an
arbitrary reduced form, very small `A` can give exponentially many points
on `y=0` under a cap comparable to `sqrt(D)`.

## 2. Exact enumeration and complete cost

For `0<Q(x,y)<m`, use integer bounds

```math
|y|\le\left\lfloor\sqrt{4A(m-1)/D}\right\rfloor,
\quad |2Ax+By|\le\left\lfloor\sqrt{4A(m-1)-Dy^2}\right\rfloor.
                                                               \tag{3}
```

At each integer `y`, the second inequality is one explicitly computable
integer interval for `x`. Integer square roots and ceiling/floor division
give every point, with no floating-point boundary decisions. Simultaneous
negation allows `y>=0`, retaining `x>0` when `y=0`.

For a conservative signed count, there are at most
`4sqrt(Am/D)+1` rows, each with at most `2sqrt(m/A)+1` points.
Put `M=max(1,m/sqrt(D))`. Since reducedness gives `A<=sqrt(D/3)` and
(1) gives `sqrt(m/A)<=2m/sqrt(D)`, their product is at most

```math
8m/\sqrt D+4\sqrt{Am/D}+2\sqrt{m/A}+1\le17M.             \tag{4}
```

The number of row setups is also `O(M)`. All endpoints, coordinates and
norms have `O(log D+log m)` bits. This proves (2), including empty rows.
Gcd tests filter primitive vectors and norms coprime to `D`. Trial division
through `Y` recognizes smoothness and supplies the factorization; the
[existing norm-word decoder](smooth-norm-relations.md) obtains orientations
from the represented vector, without a factoring oracle or balanced split.
Setup, all rejected tests and decoding therefore have total cost
`(1+m/sqrt(D))*poly(Y,log D+log m)`.

The [enumerator](../tools/near_height_smooth_descent.py) implements (3).
For a supplied unramified split-prime root it also constructs the oriented
`p^k` ideal by the existing Hensel lift and tracks the unimodular reduction
matrix to recover its known primitive vector. This costs polynomial time
in `log m`, `log D` and the factor-base setup cost. Only this near-height
use materializes `p^k`; it has `O(b)` bits in the application below.

## 3. What would suffice for a complete bound

Fix `alpha>1/2`, `b=bits(D)` and
`Y=ceil(exp((log b)^(1+alpha)))`, eventually above the pinned
[full generating cutoff](polylog-prime-form-generators.md).
Choose an explicitly computable threshold

```math
L_D\ge\sqrt D,\qquad
L_D/\sqrt D=\exp(O(b/(\log b)^\alpha)),
```

and, for every unramified split generator `p<=Y`, let `k_p` be the
least integer with `p^(k_p)>L_D`. Then
`m_p=p^(k_p)<=p L_D<=Y L_D`, and `log m_p=O(b)`.

The missing arithmetic input is precise:

> For every sufficiently large relevant discriminant and every such
> target `red(g_p^(k_p))`, either the leading norm is eligible, or there
> is a primitive eligible smooth represented norm `n<m_p`.

This is **unproved**. Existence of a norm `n<=m_p` is automatic from
`p^(k_p)` and is insufficient. Global smooth-class density, density on
covered subgroups, and the known boundary norm do not establish the strict
inequality. Unlike the older short-box interface, no additional coordinate
localization assumption is needed if this strict-descent input holds.

For an eligible leading norm use `(1,0)`; otherwise exhaust (3).
The decoded word `v_p` satisfies `prod q^|v_(p,q)|=n<m_p`, so
`k_p e_p-v_p` is a genuine relation, with weighted diagonal margin at
least `log(m_p/n)>0`. Distinct target primes therefore give a nonsingular
square relation matrix by the existing logarithm-weighted dominance proof.
The factor-base size is at most `Y`, each full search has cost
`(1+Y L_D/sqrt(D))*poly(Y,b)`, and the sum over all targets is
`exp(O(b/(log b)^alpha))`. All matrix entries have polynomial bit length.
The [determinant-free extractor](determinant-free-ambiguous-extraction.md)
has polynomial cost in this encoded matrix size. Full generation, under
the pinned analytic contract, supplies the required ambiguous directions;
primality, perfect-power handling and recursive splitting have polynomial
overhead. This is a conditional complete-cost transfer, not a new bound.

## 4. A known boundary norm need not descend

Let `m=2^k`, and choose an odd prime `a` in `(m/4,m/2)`. Bertrand's
postulate supplies such an `a` for all sufficiently large `k`. Set

```math
\Delta=1-4ma,\qquad D=4ma-1,\qquad Q=(a,-1,m).
```

This is a primitive strict reduced form; `2` is unramified and split.
The oriented norm-`m` form `(m,1,a)` reduces to `Q`, so
`Q=red(g_2^k)` for the norm-two orientation with root one. Also
`m/2<=floor(sqrt(D/3))<m`: this is exactly the first power used by
the original leading-power selector. The known vector `(0,1)` has norm
`m` and is primitive.

Nevertheless every norm below `m` has `y=0`. For `|y|=1`,
`a x^2-xy>=0`, giving `Q(x,y)>=m`. For `|y|>=2`, completing the
square gives `Q(x,y)>=(m-1/(4a))y^2>m`. At `y=0` every positive
norm is divisible by `a`, and the only primitive one is `a` itself.
Consequently there is **no** `Y`-smooth norm below `m` when `Y<a`.
In particular any cutoff `Y=2^o(b)` is eventually covered by this
order-discriminant counterfamily.

This refutes an unconditional first-power descent statement for all
quadratic orders, even though the entire search has constant size.
It is not a factoring lower bound: the discriminants have not been
restricted to inputs surviving primality and small-factor preprocessing.
It does not exclude descent at `m/sqrt(D)` growing subexponentially,
later prime powers, other target choices or useful partial presentations.
No larger parameter sweep is justified merely by this counterfamily.

## 5. Fixed evidence and attribution

```sh
make verify-near-height-smooth-descent
```

The [checker](../tools/check_near_height_smooth_descent.py) uses six
predeclared prime-power targets from the existing 60-, 80- and 92-bit
finite inputs, at the first eligible exponent and one subsequent power.
It compares the exact intervals with an independent direct signed-box
exhaustion, verifies known primitive norm witnesses, checks (1) and (4)
after leading rejection, and replays every decoded word with polynomial
ideal arithmetic. The [receipt](../evidence/near-height-smooth-descent.json)
also contains four explicit first-power failures, using checked rough
coefficients; their primality is unnecessary for the finite controls.
The 60-bit prime-seven target at exponent twelve has a rough leading
coefficient and two smaller smooth primitive norms. These are finite
examples, not a uniform density or complexity theorem.

The square-completion and lattice-point estimates are classical elementary
facts. The existing supplied-norm decoder and polynomial ideal arithmetic
retain their source attribution; the principal-word gap was already
recorded by orda-aska. Sol worker 1 developed this target-specific
enumeration deduction and counterfamily with substantial OpenAI assistant
support. This is an author-checked ordinary argument and finite replay,
without independent mathematical review or new Lean formalization.
