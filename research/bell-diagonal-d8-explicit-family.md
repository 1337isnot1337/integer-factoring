# Linear and monomial Bell-diagonal paths at d=8 fail the t=2 width target

## Result

At physical dimension `d=8` (`q=64`) over `F_3`, every **linear** path

```math
s_j(z)=1+(z-1)v_j
```

with a ternary vector `v` of length `64` has a length-two profile rank
of at least `1848`. The allowed rank at `rho<1/sqrt(2)` and `t=2` is
`q^{2rho}<64^{sqrt(2)}<1848`. Therefore no linear ternary path can meet
the two-cut width target at `t=2`, and hence no such path gives a
uniform crossing family.

The same ordered-word census also excludes **monomial** local polynomials
`z^{d_j}` with `d_j in {0,1,2}`: the minimum possible maximum row rank is
`1173`, also above the target. This closes the linear and monomial
quadratic path families, not the whole Bell-diagonal `d=8` construction.
Non-monomial quadratic paths and other characteristics remain open.

## Count census

For a path with value counts `(a,b,c)`, `a+b+c=64`, the length-two rank
profile is `linear_profile_ternary((a,b,c),2)`. Enumerating all 2,145
count triples gives a minimum of the maximum profile entry

```math
\min_{a+b+c=64}\ \max_j P_2(j)=1848,
```

attained at `(22,42,0)` and `(42,22,0)`.
The transform counts are not needed for this lower bound: the row cut
alone has rank at least `1848`, and the column cut has the same
minimum over its own count triple. Hence every linear ternary path has
`max(R_0,R_infinity)\ge1848`.

The target at `t=2` is

```math
64^{2/\sqrt2}=64^{\sqrt2}\approx 357.0<1848.
```

The inequality is strict and the census is finite.

For monomial local polynomials `z^{d_j}`, with `n_0,n_1,n_2` counts
summing to `64`, the length-two coefficient rank at degree `d` is the
number of **ordered** words with `d_i+d_j=d`:

```math
(n_0^2,\ 2n_0 n_1,\ n_1^2+2n_0 n_2,\ 2n_1 n_2,\ n_2^2).
```

The ordered form is essential: `linear_profile_ternary` multiplies
literal words, not unordered pairs. Enumerating all `2,145` triples
gives a minimum of the maximum entry `1173`, attained at
`(n_0,n_1,n_2)=(13,17,34)` and `(34,17,13)`. This is also above
`64^{\sqrt2}`; every monomial quadratic path is excluded at `t=2` on
the row cut.

## Relation to the earlier candidate

An explicit path was previously displayed because its **rank fraction**
`R/q^t` was below one. That was the wrong criterion: the target is
`R<=q^{rho t}`, so the allowed fraction at `t=2` is
`64^{-(1-rho)2}`, about `0.086` at `rho=1/sqrt(2)`, not `0.66`. The
candidate's ranks are in fact far above the target, and the census
shows that no linear path does better than `1848`.

## What remains open

The remaining quadratic paths have at least one non-monomial local
polynomial. A successful construction must have both cut ranks at most
`q^{rho t}` for every `t`, so in particular at `t=2` its profile rank
must be below `358`, more than three times smaller than the monomial
minimum. The next search should target non-monomial local
eigenpolynomials, not refine linear or monomial paths.

## Reproduction

```sh
make verify-bell-d8-explicit-family
```

The checker enumerates all count triples, verifies the linear minimum
`1848` and the monomial ordered-word minimum `1173`, checks the strict
inequalities against `64^{sqrt2}`, and verifies the explicit path's
profiles. It does not enumerate nonlinear paths.
