# Bell-diagonal d=8 paths of degree at most five fail at t=2

## Result

At physical dimension `d=8` (`q=64`) over `F_3`, consider any polynomial
path `s_j(z)` with `s_j(1)=1` for every local eigenvalue `j`, and let
`D` be the maximum local degree. The length-two tensor
`P_{s(z)}^{tensor 2}` has at most `2D+1` coefficient tensors, indexed by
the powers `z^0,...,z^{2D}`.

In the row Bell basis, the rank of a coefficient is the number of words
whose coefficient polynomial is nonzero, with the multinomial
multiplicity of the word type. There are `q^2=4096` words. For every
word, the product
`s_{j_1}(z)s_{j_2}(z)` is nonzero because it equals `1` at `z=1`.
Thus the word contributes to at least one coefficient, and the sum of
the row ranks over all `z`-degrees is at least `4096`. Therefore some
coefficient has row rank at least

```math
\left\lceil {4096\over 2D+1}\right\rceil .
```

For `D<=5` the smallest such bound is `ceil(4096/11)=373`, while the
target at `t=2` and `rho<1/sqrt(2)` allows at most
`64^{2rho}<64^{sqrt(2)}\approx 357`. Hence **no Bell-diagonal
tensor-power path of local degree at most five can meet the two-cut
width target at `t=2`**, and none gives a uniform crossing family.

This is a row-cut bound and does not use the Bell transform. The linear
and monomial quadratic censuses below are sharper special cases; the
universal bound is what closes all low-degree paths, including
non-monomial quadratic ones.

## Degree table

| max degree `D` | coefficients `2D+1` | row bound `ceil(4096/(2D+1))` |
| ---: | ---: | ---: |
| 0 | 1 | 4096 |
| 1 | 3 | 1366 |
| 2 | 5 | 820 |
| 3 | 7 | 586 |
| 4 | 9 | 456 |
| 5 | 11 | 373 |
| 6 | 13 | 316 |

The bound crosses the target between `D=5` and `D=6`; this argument
does not exclude degree at least six.

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

Only paths of local degree at least six survive the degree table. At
`D=6` the average bound is `316`, below the target, so the counting
argument no longer excludes them. A successful construction must have
both cut ranks at most `q^{rho t}` for every `t`; at `t=2` its profile
rank must be below `358`, more than three times smaller than the
monomial minimum. The next search should target genuinely
higher-degree local eigenpolynomials, or a different Bell-diagonal
representation; refining linear, monomial, or non-monomial quadratic
paths cannot work.

## Reproduction

```sh
make verify-bell-d8-explicit-family
```

The checker verifies the universal degree table through `D=5`, the
linear minimum `1848`, the monomial ordered-word minimum `1173`, and
the explicit path's profiles. It does not enumerate nonlinear paths of
degree at least six.
