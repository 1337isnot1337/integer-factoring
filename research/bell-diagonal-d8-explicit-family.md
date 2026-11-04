# Linear Bell-diagonal paths at d=8 fail the t=2 width target

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

This closes the simple linear path family, not the whole Bell-diagonal
`d=8` construction. Nonlinear paths, position-dependent families, and
other characteristics remain open.

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

## Relation to the earlier candidate

An explicit path was previously displayed because its **rank fraction**
`R/q^t` was below one. That was the wrong criterion: the target is
`R<=q^{rho t}`, so the allowed fraction at `t=2` is
`64^{-(1-rho)2}`, about `0.086` at `rho=1/sqrt(2)`, not `0.66`. The
candidate's ranks are in fact far above the target, and the census
shows that no linear path does better than `1848`.

## What remains open

The Bell-diagonal `d=8` note leaves open nonlinear or
position-dependent families; this result does not exclude them. A
successful construction must have both cut ranks at most
`q^{rho t}` for every `t`, so in particular at `t=2` its profile rank
must be below `358`, more than five times smaller than the linear
minimum. The next search should target genuinely nonlinear local
eigenpolynomials; refining linear paths cannot work.

## Reproduction

```sh
make verify-bell-d8-explicit-family
```

The checker enumerates all count triples, verifies the minimum `1848`,
verifies the explicit path and its profiles, and checks the strict
inequality against `64^{sqrt2}`. It does not enumerate nonlinear paths.
