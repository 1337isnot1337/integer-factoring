# Exact polynomial synthesis with an all-length rank bound

This method constructs paths in the [Bell local space](bell-diagonal-crossings.md)
by exact affine constraints, then certifies a complete rank-width upper bound.
The search order is heuristic; every retained path and rational bound is
checked independently. This lane has not improved the factoring exponent.
It led to the explicit [folded-cube family](folded-bell-crossings.md), whose
all-length coefficient-span limit is proved separately.

## 1. Uniform bound from scalar polynomial orders

Take local physical dimension `d=2^b`, and a polynomial vector `s(z)`
over `F_3` of fixed degree at most `h`, with `s(1)=1`. Its opposite-cut
vector is `Ts(z)`. Every local eigenpolynomial is nonzero, since its
value at one is one. Write `l_j` for the orders at zero of `s_j`, and
`u_j` for the degrees of `(Ts)_j`.

For a rational `x>=1` put

```math
A(x)=\sum_j x^{-l_j},\qquad B(x)=\sum_j x^{u_j}.
```

At tensor length `t`, a word contributes to coefficient `i` only if
its summed row orders are at most `i`, and its summed opposite degrees
are at least `i`. Counting those words gives the exact upper bounds

```math
r_i\le x^i A(x)^t,\qquad
c_i\le x^{-i}B(x)^t,\qquad
r_i c_i\le(A(x)B(x))^t.
\tag{1}
```

In the first sum replace each indicator for total order `L<=i` by
`x^(i-L)`; in the second replace the indicator for total degree `U>=i`
by `x^(U-i)`. Both resulting sums factor into local products. The cut
directions may be exchanged. This bound does not assume every
coefficient inside a polynomial's support interval is nonzero.

A rational pressure `Q=A(x)B(x)` proves leading contraction
`rho=log2(Q)/(4b)`, uniformly for every `t`. The path has `ht+1`
pieces and effective factors by the parent construction. Arbitrary
label widths use its restriction. A successful improvement still
needs `rho<1/sqrt(2)`.

Receipts certify `rho<=a/1000` by `Q^1000<=2^(4ba)`. The improvement
flag tests `Q^577<2^(4b*408)`; `408/577<1/sqrt(2)` makes this a
rigorous sufficient criterion. Floating optimization selects a scale
only; the final pressure and inequalities use rational arithmetic.

## 2. Exact affine synthesis

Write `s(z)=sum_(i=0)^h a_i z^i`. To prescribe lower row orders `l_j`
and upper opposite degrees `u_j`, impose

```math
(a_i)_j=0\ (i<l_j),\qquad
(Ta_i)_j=0\ (i>u_j),\qquad
\sum_i(a_i)_j=1.
\tag{2}
```

These are affine equations over `F_3`. Feasible targets return a full
polynomial path. A Gaussian basis stores ternary rows as two disjoint
bit sets. Independent dense elimination checks feasibility, solutions,
and preservation of the old basis after rejected constraints.

The driver starts from tensor products of one-bit paths, relaxes order
constraints, and greedily rebuilds them. An adaptive schedule uses the
best prior polynomial as its next starting point. This is a heuristic
target order, not exhaustive polynomial optimization.

## 3. Retained schedules and their limits

| Physical dimension | Permitted degree | Schedule | Certified rho upper bound |
| --- | ---: | --- | ---: |
| 4 | 4 | 1,024 adaptive affine trials | 0.876 |
| 8 | 3 | 2,048 adaptive affine trials | 0.873 |

The [dimension-four receipt](../evidence/bell-synthesis-d4h4.json) has
31,425 exact constraint proposals and 22 distinct returned profile pairs.
The [dimension-eight receipt](../evidence/bell-synthesis-d8h3.json) has
199,594 proposals and 39 distinct profile pairs. The initial product
pressure gave rho `log2(12)/4≈0.896241`. These improvements within this
route remain above the threshold. The later folded-cube formula is
stronger than these particular witnesses and has a proved span obstruction.

```sh
python3 tools/search_bell_synthesis.py --bits 2 --degree 4 --trials 1024 --tensor-initial --adaptive --output /tmp/bell-synthesis-d4h4.json
python3 tools/search_bell_synthesis.py --bits 3 --degree 3 --trials 2048 --tensor-initial --adaptive --output /tmp/bell-synthesis-d8h3.json
python3 tools/check_bell_profiles.py evidence/bell-synthesis-d4h4.json evidence/bell-synthesis-d8h3.json
```

Readback checks every coefficient, transformed coefficient, endpoint,
rational pressure and integer width inequality. Literal length-two
tensor-word ranks independently satisfy (1); 2,448 dense affine controls
cover feasible and inconsistent systems. No independent proof review or
Lean formalization is claimed. Substantial OpenAI Codex assistance was used.

The next target is a structurally different uniform family, including
possible code quotients, rather than extrapolating these finite schedules.
