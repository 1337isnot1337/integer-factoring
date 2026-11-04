# Fixed-degree Bell-diagonal paths fail the width target

## Result

At physical dimension `d=8` (`q=64`) over `F_3`, consider any polynomial
path `s_j(z)` with `s_j(1)=1` for every local eigenvalue `j`, and let
`D` be the maximum local degree. At tensor length `t`, the tensor
`P_{s(z)}^{tensor t}` has at most `tD+1` coefficient tensors, indexed
by the powers `z^0,...,z^{tD}`.

In the row Bell basis, the rank of a coefficient is the number of words
whose coefficient polynomial is nonzero, with the multinomial
multiplicity of the word type. There are `q^t` words. For every word,
the product of its local polynomials is nonzero because it equals `1`
at `z=1`. Thus the word contributes to at least one coefficient, and the
sum of the row ranks over all `z`-degrees is at least `q^t`. Therefore
some coefficient has row rank at least

```math
\left\lceil {q^t\over tD+1}\right\rceil .
```

The target at `rho<1/sqrt(2)` allows at most `q^{rho t}`. Since

```math
{q^t/(tD+1)\over q^{rho t}}
={q^{(1-rho)t}\over tD+1}\longrightarrow\infty,
```

the average bound exceeds the target for some `t`, for every **fixed**
`D`. Hence no Bell-diagonal tensor-power path of bounded local degree
can satisfy the two-cut width inequality at every length, and none gives
a uniform crossing family.

This is a row-cut argument and does not use the Bell transform. At
`t=2` it already excludes every path of degree at most five; a degree-six
path fails at `t=3`. The linear and monomial quadratic censuses below
are sharper finite special cases.

## Explicit failures

At `t=2`:

| max degree `D` | coefficients `2D+1` | row bound `ceil(4096/(2D+1))` |
| ---: | ---: | ---: |
| 0 | 1 | 4096 |
| 1 | 3 | 1366 |
| 2 | 5 | 820 |
| 3 | 7 | 586 |
| 4 | 9 | 456 |
| 5 | 11 | 373 |
| 6 | 13 | 316 |

The target is `64^{2rho}<64^{sqrt(2)}\approx357`, so `D<=5` fails at
`t=2`. At `D=6`, the first failure is `t=3`:
`ceil(64^3/(18+1))=13798>64^{3/sqrt(2)}\approx6784`.

The general threshold is only an existence statement: for each fixed
`D`, some finite `t` fails. It does not claim a uniform bound on that
`t`.

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

## What remains open

The escape from the fixed-degree barrier would require the local degree
`D_t` to grow with `t`, at least fast enough that
`q^{(1-rho)t}/(tD_t+1)` stays bounded. That is a different
representation: the path must interpolate many more tensor coefficients,
and the decomposition cost of those coefficients must still be
controlled. Increasing a fixed degree, or searching a larger fixed
family, cannot work. A different Bell-diagonal representation or a
non-path construction is needed.

## Reproduction

```sh
make verify-bell-d8-explicit-family
```

The checker verifies the degree table, the linear minimum `1848`, the
monomial ordered-word minimum `1173`, the explicit path's profiles, and
the first fixed-degree failures through `D=12`. It does not enumerate
paths of unbounded degree.
