# Polynomial Bell-basis decompositions of fixed degree fail the width target

## Result

Work in the row Bell basis of the `d=8` (`q=64`) crossing over `F_3`.
Let

```math
F(z)=\sum_{a=0}^{D} z^a A_a
```

be any polynomial tensor with `F(1)=cP_1`, where `P_1` is the Bell
crossing tensor and `c!=0`. In the row basis, `P_1` has coefficient `1`
on every word. Hence for each of the `q^t` words, the word coefficient
polynomial satisfies `f_\nu(1)=c!=0`, so it is nonzero and contributes
to at least one `z`-degree.

The row rank of `A_a` is the number of words whose coefficient at
degree `a` is nonzero. Therefore the sum of the row ranks over
`a=0,...,D` is at least `q^t`, and some coefficient has row rank at
least

```math
\left\lceil {q^t\over D+1}\right\rceil .
```

The target at length `t` is `q^{rho t}`. Since

```math
{q^t/(D+1)\over q^{rho t}}
={q^{(1-rho)t}\over D+1}\longrightarrow\infty,
```

for every **fixed** total degree `D` there is a finite `t` at which the
row cut alone exceeds the target. Thus no fixed-degree polynomial
Bell-basis decomposition can satisfy the two-cut width inequality at
every length.

At `d=8,t=2`, `q^t=4096`; the target is
`64^{2rho}<64^{sqrt(2)}\approx357`. The bound is
`ceil(4096/(D+1))`, which is `373` at `D=10` and `342` at `D=11`.
Hence every polynomial Bell decomposition of total degree at most ten
fails at `t=2`.

A Bell tensor-power path with local degree `D_path` has total degree
`tD_path`; the earlier fixed-degree path barrier is the special case
`D=tD_path`. The linear and monomial quadratic censuses below are
sharper finite special cases.

## Degree table

At `t=2`:

| total degree `D` | coefficients `D+1` | row bound `ceil(4096/(D+1))` |
| ---: | ---: | ---: |
| 1 | 2 | 2048 |
| 2 | 3 | 1366 |
| 3 | 4 | 1024 |
| 4 | 5 | 820 |
| 5 | 6 | 683 |
| 6 | 7 | 586 |
| 7 | 8 | 512 |
| 8 | 9 | 456 |
| 9 | 10 | 410 |
| 10 | 11 | 373 |
| 11 | 12 | 342 |

The bound crosses the target between total degree ten and eleven; this
argument does not exclude degree at least eleven.

The same calculation applies at other physical dimensions. For `d=32`
(`q=1024`) at `t=2`, there are `q^2=1,048,576` words and the target is
`1024^{sqrt(2)}≈18,064`. A cubic path (total degree three) has four
coefficients, so some rank is at least `ceil(1,048,576/4)=262,144`; a
quartic path has five coefficients and bound `209,716`. Both are far
above the target. Thus the dimension-32 cubic/quartic paths left open by
the outer-rank endpoint pressure test are excluded at actual rank by the
same word-count argument.

For a path of local degree `D_path`, the total degree is `2D_path`.
Thus local degree at most five is excluded at `t=2` (`D=10`), while
local degree six has total degree `12`; its first failure is `t=3`,
where `ceil(64^3/(18+1))=13798>64^{3/sqrt(2)}\approx6784`.

## Linear and monomial refinements

For the linear family `s_j(z)=1+(z-1)v_j`, write the value counts as
`(a,b,c)`, `a+b+c=64`. The length-two profile is
`linear_profile_ternary((a,b,c),2)`. Enumerating all `2,145` count
triples gives a minimum of the maximum profile entry

```math
\min_{a+b+c=64}\ \max_j P_2(j)=1848,
```

attained at `(22,42,0)` and `(42,22,0)`. The transform counts are not
needed for this lower bound: the row cut alone has rank at least
`1848`, and the column cut has the same minimum over its own count
triple. The earlier explicit path has rank exponents `0.904` to `0.998`,
all above `1/sqrt(2)`; its finite rank fractions below one were the
wrong criterion.

For monomial local polynomials `z^{d_j}` with `d_j in {0,1,2}` and counts
`n_0,n_1,n_2` summing to `64`, the length-two row profile is the number
of **ordered** words with `d_i+d_j=d`:

```math
(n_0^2,\ 2n_0 n_1,\ n_1^2+2n_0 n_2,\ 2n_1 n_2,\ n_2^2).
```

The ordered form is essential: `linear_profile_ternary` multiplies
literal words, not unordered pairs. Enumerating all `2,145` triples
gives minimum maximum `1173`, attained at `(13,17,34)` and
`(34,17,13)`. This is another special case of the universal bound.

## Relation to the earlier candidate

An explicit path was previously displayed because its **rank fraction**
`R/q^t` was below one. That was the wrong criterion: the target is
`R<=q^{rho t}`, so the allowed fraction at `t=2` is
`64^{-(1-rho)2}`, about `0.086` at `rho=1/sqrt(2)`, not `0.66`. The
candidate's ranks are far above the target, and the censuses above show
that no linear or monomial path does better.

## Cost consequence

Let `m` be the number of nonzero coefficient tensors. The same word
count gives `max_a rank(A_a) >= q^t/m`. To meet the width target
`q^{rho t}`, a family must therefore have

```math
m >= q^{(1-rho)t}=2^{(1-rho)lambda},
```

where `lambda=t log_2 q` is the label width. This is exponential in
`lambda` for every fixed `rho<1`. The Bell polynomial synthesis
enumerates its `D+1` coefficient tensors (or `Dt+1` path pieces);
if it meets the width target, that enumeration is exponential in
`lambda`, while the target complete bound is only
`2^{O(b/(log b)^{1/2})}`. Thus the polynomial synthesis route cannot
improve the complete bound under its published piece-cost contract.

A compact evaluator that produces the many coefficients without
enumerating them would be a different representation; none is supplied.
This is a method-specific cost obstruction, not a factoring lower
bound.

## What remains open

The escape from this barrier requires either unbounded total degree or a
non-polynomial family, and the cost consequence shows that a direct
polynomial synthesis of the necessary exponential coefficient list
cannot improve the bound. A compact non-enumerative evaluator, a
genuinely different Bell-basis representation, or a non-path
construction remains open.

## Reproduction

```sh
make verify-bell-d8-explicit-family
```

The checker verifies the total-degree table through `D=11`, the
path-local degree table through `D=6`, the linear minimum `1848`, the
monomial ordered-word minimum `1173`, and the explicit path's profiles.
It does not enumerate unbounded-degree or non-polynomial families.
