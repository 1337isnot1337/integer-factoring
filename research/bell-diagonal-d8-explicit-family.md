# A d=8 Bell-diagonal finite profile is not a width crossing

## Correction and result

An earlier version of this note compared the **rank fraction**
`R/q^t` with `1/sqrt(2)`. That is not the width criterion. The
`pieces-to-factoring` target requires the rank itself to be at most
`q^{rho t}` with `rho<1/sqrt(2)`, equivalently a rank fraction at most
`q^{-(1-rho)t}`, which is exponentially small. A fraction below `1` is
not enough.

For the explicit linear path

```math
s_j(z)=1+(z-1)v_j
```

over `F_3` at `d=8` (`q=64`), with `v` the ternary vector in the
[receipt](../evidence/bell-d8-explicit-family.json), the exact profile
ranks for lengths `1..100` give rank fractions between about `0.43` and
`0.67`. The corresponding **rank exponents**

```math
e_t=\frac{\log_64 R_t}{t}
```

are `0.904` at `t=1`, `0.950` at `t=2`, `0.984` at `t=10`, `0.995` at
`t=31`, and `0.998` at `t=100`. They are all above `1/sqrt(2)`; the
path's ranks grow almost like the full `q^t`, not like `q^{0.707t}`.

The path therefore does **not** approach the required width. The
Bell-diagonal `d=8` constructive target remains open; this finite
profile is a negative screen, not a candidate.

## Why the fraction is misleading

The local total at length `t` is `q^t=64^t`. A crossing piece of rank
`R` must satisfy `R <= 64^{rho t}`; at `rho=0.707`, the allowed fraction
is `64^{-0.293t}`: about `0.086` at `t=2`, `0.0013` at `t=10`, and
`1.6e-30` at `t=100`. The observed fractions `0.43..0.67` are many
orders of magnitude above that. The determinant floor
`gamma_64=0.656961...` is a lower bound on the **exponent** of a forced
rank, not an upper bound on a construction; it does not make a rank
fraction near `0.66` favorable.

## Reproduction

```sh
make verify-bell-d8-explicit-family
```

The checker verifies the transform and counts, the `T^2=I` identity,
the finite profiles, the rank fractions, and the rank exponents for
lengths `1..100`. It asserts that every exponent is above
`1/sqrt(2)` and records the negative result. It does not prove a
uniform all-length statement; it proves only that this explicit path is
not a crossing.
