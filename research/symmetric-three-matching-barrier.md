# A sharp rate barrier for symmetric three-matching pieces

## Result and scope

Allowing all three coordinate pairings does not rescue a construction
that remains invariant under every permutation of its binary coordinates.
This holds for arbitrary linear recombinations and arbitrarily many
pieces in that invariant space.

Work over any field of characteristic two. On one-bit labels define

```math
E=[x=v][y=u],\qquad R=[x=y][v=u],\qquad
S=[x=u][v=y],\qquad F=E+R+S.
```

For `r+s+f=t`, let `T_(r,s,f)` be the sum of all **distinct tensor
words** with `r` copies of `R`, `s` copies of `S`, and `f` copies
of `F`. Each word occurs once; no factorial division or multinomial
coefficient is inserted. Put

```math
\mathcal V_t=\operatorname{span}\{T_{r,s,f}:r+s+f=t\}.
\tag{1}
```

These type tensors are linearly independent. The three local tensors
are independent, so their ordered tensor words are independent; the
types are sums of disjoint sets of those words. Thus their coefficients
are unambiguous even in characteristic two.

**Theorem.** Suppose any finite family `B_j in V_t` satisfies

```math
\sum_jB_j=dE^{\otimes t},\qquad d\ne0.
```

Writing `r(B),c(B)` for the two crossing flattening ranks,

```math
\max_j r(B_j)c(B_j)
\ge
\frac{B(t,\lfloor t/2\rfloor)B(t,\lceil t/2\rceil)}
     {16(t+1)^6},\qquad
B(t,i)=\sum_{a=0}^{i}\binom{2t}{a}.
\tag{2}
```

Consequently

```math
\max_j\log_2(r(B_j)c(B_j))
\ge4tH_2(1/4)-O(\log(t+2)).
\tag{3}
```

This leading coefficient is sharp: the elementary coefficient family
in [binary-coefficient-entropy.md](binary-coefficient-entropy.md)
belongs to `V_t` and attains the matching upper bound.
The exact leading rate for this entire class is therefore

```math
\rho=H_2(1/4)=0.8112781244\ldots>1/\sqrt2.
```

There is **no piece-count hypothesis** in (2). The theorem covers
coefficient families of any common polynomial curve
`(a(z)R+b(z)S+c(z)F)^(tensor t)`, and arbitrary linear combinations
of their coefficients, whenever their sum is the crossing.
It does not cover coordinate-dependent curves, arbitrary label subsets,
fixed-weight restrictions, or general crossing tensors.

The inclusion-rank input is Wilson's theorem, with its primary source
and exact formula stated in the companion binary-coefficient note.
The two filtrations and the walk-count comparison below give the new
extension from that one family to (1).

## 1. Local nilpotents and two rank-preserving changes of scale

In the row flattening, `E=I_4`. The matrices `R` and `T=I_4+S`
have rank one, square zero, and complementary two-dimensional supports.
In a fixed basis they are `N direct-sum 0_2` and
`0_2 direct-sum N`, where `N` is one nilpotent Jordan block.
In particular

```math
S=I+T,\qquad F=R+T.
\tag{4}
```

For nonzero `B in V_t`, let `s_*` be the largest `S` count with
a nonzero type coefficient. Set `D=t-s_*`. After (4), the least
total nilpotent degree is exactly `D`. Indeed the top `S` slice
becomes a nonzero combination of types `I^(s_*)R^kT^(D-k)`.
If its original coefficients are `a_r`, the new ones are

```math
b_k=\sum_{r\le k}\binom kr a_r.
\tag{5}
```

This Pascal transform is unitriangular, so it cannot annihilate a
nonzero slice. Its first nonzero coefficient has index equal to the
least `r` with `a_r!=0`; call it `k`.

Here is a precise rank argument for discarding higher degrees.
Each local two-dimensional support is a copy of
`K[x]/(x^2)`, with `N` acting by multiplication by `x`.
Conjugate by the diagonal scaling of degree `a` by an indeterminate
`z^a`. Multiplication of nilpotent degree `j` acquires factor `z^j`.
Divide by `z^D` and set `z=0`. The resulting matrix has rank no
larger than the original: over `K(z)` the scaling is invertible,
and specialization cannot increase matrix rank.

Now scale only the variables on the `R` supports. In the degree-`D`
matrix, divide by the least occurring power and specialize again.
This keeps exactly the type with `k` copies of `R` and `D-k`
copies of `T`. Both operations give rank lower bounds for the
original matrix. No evaluation at a special parameter is used to
assert a generic-rank upper bound.

## 2. Reduce the surviving type to inclusion matrices

Choose the `R` support at `a` coordinates and the `T` support at
the others. There are `binom(t,a)` such blocks. On that block,
the surviving operator is multiplication by

```math
e_k(x_1,\ldots,x_a)\,
e_{D-k}(y_1,\ldots,y_{t-a})
```

in the square-zero algebra on all `t` variables. Let `u(v,j)`
be the rank of multiplication by `e_j` on the algebra on `v`
variables. Ranks multiply for tensor products and add for direct sums.
Therefore

```math
r(B)\ge L(t,D,k):=
\sum_{a=0}^{t}\binom ta u(a,k)u(t-a,D-k).
\tag{6}
```

Define the binomial tail

```math
U(v,j)=\sum_{2b\le v-j}\binom vb,
```

with value zero when `j>v`. Wilson's formula, as applied in the
companion note, gives

```math
\frac{U(v,j)}{v+1}\le u(v,j)\le2U(v,j).
\tag{7}
```

Put

```math
W(t,D,k)=\sum_{a=0}^{t}\binom ta U(a,k)U(t-a,D-k).
```

Equations (6)--(7) imply

```math
L(t,D,k)\ge\frac{W(t,D,k)}{(t+1)^2},\qquad
L(t,D,0)\le4W(t,D,0).
\tag{8}
```

## 3. Splitting the nilpotent degree cannot lower the leading rate

The number `W(t,D,k)` counts length-`t` words in the four symbols
`R+`, `R-`, `T+`, `T-`, whose two net counts satisfy

```math
R_+-R_-\ge k,\qquad T_+-T_-\ge D-k.
```

Replace each symbol by a pair of signs:

| Symbol | First sign | Second sign |
|---|---:|---:|
| `R+` | `+1` | `+1` |
| `R-` | `-1` | `-1` |
| `T+` | `+1` | `-1` |
| `T-` | `-1` | `+1` |

This is a bijection with pairs of independent sign words. Let their
sums be `u,v`. The two inequalities become

```math
u\ge D,\qquad 2k-u\le v\le u-2D+2k.
\tag{9}
```

For fixed `u`, the allowed interval for `v` has radius `u-D`
and center `2k-D`. Its endpoints have the same parity as `t`.
The number of sign words with sum `v` is
`binom(t,(t+v)/2)`, a symmetric unimodal sequence.
Among intervals of this fixed radius whose centers range from `-D`
to `D` in steps of two, an endpoint center has the smallest mass.
One can check each shift directly: moving a nonnegative center farther
from zero removes a left endpoint and adds an endpoint of no greater
binomial weight.

It follows for every `0<=k<=D` that

```math
W(t,D,k)\ge W(t,D,0).
\tag{10}
```

The exact binary coefficient formula gives
`L(t,D,0)=r_(t-D)`, where `r_i` is the rank in the companion note.
Combining (6), (8), (10), and its binomial-tail lower bound proves

```math
\boxed{\quad
r(B)\ge\frac{r_{s_*}}{4(t+1)^2}
       \ge\frac{B(t,s_*)}{4(t+1)^3}.\quad}
\tag{11}
```

In the column flattening, `R` and `S` exchange roles. If `r_*`
is the largest `R` count in `B`, the same proof gives

```math
c(B)\ge\frac{B(t,r_*)}{4(t+1)^3}.
\tag{12}
```

## 4. A balanced type forces the crossing bound

Since `E=R+S+F`, every type coefficient of `E^(tensor t)` is
one. In a decomposition from the theorem, some summand has a nonzero
coefficient on

```math
T_{\lfloor t/2\rfloor,\,\lceil t/2\rceil,\,0}.
```

For that summand, `r_*>=floor(t/2)` and
`s_*>=ceil(t/2)`. Binomial tails are monotone. Multiplying
(11)--(12) proves (2), and the uniform entropy estimate proves (3).
No monotonicity assumption about the exact characteristic-two ranks
`r_i` is needed.

The argument is unchanged after extending the characteristic-two field,
including to a Laurent-series field. A controlled algebraic degeneration
that stays inside this symmetric span therefore has the same generic
leading obstruction. Its exact coefficient tensors also remain in the
span. A degeneration outside the span is a different question.
For an approximate identity
`sum_j B_j(z)=z^kappa E^(tensor t)+O(z^(kappa+1))`, the balanced
type coefficient is a nonzero Laurent series. Some summand has that
coefficient nonzero, so the same proof bounds its generic ranks even
though the generic sum need not be exactly a scalar multiple of the
crossing.

## 5. Verification and stopped searches

Run `python3 tools/check_symmetric_three_matching.py`. The checker
compares the filtration lower bounds with independently assembled full
tensors, tests the walk comparison by exact integer counts, and checks
the rank-to-tail inequalities. It includes a control showing why the
symmetry assumption matters: an individual nonsymmetric tensor word
can have a much smaller product of its two ranks.

The theorem was independently reviewed during this research session.
It made further polynomial-curve searches inside this span unnecessary.
Those searches were stopped rather than extrapolated to a better rate.
Their partial output does not certify exhaustive optimality at larger
parameters. The proof and checks used substantial OpenAI assistant
support; no new Lean formalization was run.
