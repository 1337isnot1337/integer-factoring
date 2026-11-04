# An explicit d=8 Bell-diagonal path with finite profile below 1/sqrt(2)

## Result

The Bell-diagonal construction has an explicit linear path at physical
dimension `d=8` (local dimension `q=64`) over `F_3` whose two-cut rank
profiles are below `1/sqrt(2)` for all tested tensor lengths `1` through
`100`. The maximum rank ratio is

```math
\frac{43}{64}=0.671875<0.7071067812...
```

The path is

```math
s_j(z)=1+(z-1)v_j,
```

with `v` the ternary vector in the
[receipt](../evidence/bell-d8-explicit-family.json). Its value counts are
`(31,21,12)` and the counts of `Tv` are `(22,30,12)`. At `z=1` every
local factor is `1`, so the path specializes to the exact scale-one
crossing; the Bell transform identity `T^2=I` holds over `F_3`.

## Profile table

The maximum rank ratios for lengths `1` through `12` are

```text
t=1  43/64              = 0.671875
t=2  2712/4096          = 0.662109375
t=3  164475/262144      = 0.627422...
t=4  9891936/16777216   = 0.589605...
t=5  633279582/1073741824 = 0.589787...
t=6  37479881439/68719476736 = 0.545404...
t=7  2433895866099/4398046511104 = 0.553404...
t=8  150591989206440/281474976710656 = 0.535010...
t=9  9519800180198400/18014398509481984 = 0.528455...
t=10 602643695601065634/1152921504606846976 = 0.522710...
t=11 37842957437614534656/73786976294838206464 = 0.512868...
t=12 2431681430557255489800/4722366482869645213696 = 0.514929...
```

The maximum is at `t=1`; for `t=31..100` the maximum is
`0.510285...` at `t=31`, and at `t=100` it is `0.429518...`. The
sequence is not monotone but stays far below the `t=1` value. The
determinant lower floor at `q=64` is `0.656961...`, so this path is
slightly above the floor but strictly below the target.

## What this does and does not establish

It establishes an explicit finite profile with both cut ranks below
`1/sqrt(2)` at every tested length. The pieces are the monomial
coefficients of the path; there are `t+1` of them, and each has both
cut ranks bounded by the displayed maximum. A uniform all-length bound
of the same form would give a two-cut width inequality with exponent at
most `0.672+o(1)`, below the Bell-diagonal barrier at `d=2,4` and below
`1/sqrt(2)`.

It does **not** establish the all-length bound, the full
`pieces-to-factoring` transfer, or a factoring exponent. The finite
table is evidence for a promising construction, not a proof; the
all-length profile bound and the transfer costs remain open. The next
concrete task is to prove that this path's maximum profile ratio stays
below a fixed rational `<1/sqrt(2)` for every `t`, or find a path with a
provable uniform bound.

## Reproduction

```sh
make verify-bell-d8-explicit-family
```

The checker recomputes the transform, the value and transform counts,
the `T^2=I` identity, and the profiles for lengths `1..100`, and verifies
the `43/64` maximum against the receipt. It does not prove the
all-length statement.
