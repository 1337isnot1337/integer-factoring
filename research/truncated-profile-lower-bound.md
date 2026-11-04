# A resolved truncated profile is a lower bound on the complete pressure

## Lemma

Fix a coherent two-cut map with nonnegative orders and a height-zero
letter. Let `O_0` and `O_infinity` be its complete order multisets, and
suppose a resolved truncation retains subsets
`O'_0 subset O_0` and `O'_infinity subset O_infinity`. For a scale
`t>0`, define

```math
C(t)=log sum_{o in O_0} e^{-t o}+log sum_{o in O_infinity} e^{-t o},
qquad
C'(t)=log sum_{o in O'_0} e^{-t o}+log sum_{o in O'_infinity} e^{-t o},
```

and

```math
H_t(rho)=log sum_h e^{-h/rho}.
```

The pressure equation is `C(t)=4 rho H_t(rho)` (respectively with `C'`).
Because all orders are nonnegative, `C'(t) <= C(t)`. The function
`F_t(rho)=4 rho H_t(rho)` is strictly increasing for `rho>0`: its
derivative is

```math
4 H_t(rho)+4 rho H_t'(rho)
=4 H_t(rho)+4/rho * E_t[h],
```

where `E_t[h] >= 0`, and `H_t(rho) >= 0` because the height-zero letter
contributes `e^0=1`. Hence the root `rho'(t)` is at most `rho(t)` for
every `t`. Taking the minimum over `t` preserves the inequality:

```math
min_t rho'(t) <= min_t rho(t).
```

So the numerical pressure computed from any resolved subset of orders,
with the same heights and scale search, is a **lower bound** on the
complete pressure.

## Consequence

A truncated profile above `1/sqrt(2)` **excludes** the complete map from
crossing: the complete pressure is at least the truncated value. A
truncated profile below `1/sqrt(2)` does **not** prove a crossing; the
unresolved orders can raise the complete pressure above the threshold.
This is why the q94 mode-11 pilot is only an exclusion, not a candidate.

## Application to the q94 mode-11 drop

The q94 map of
[the pilot note](f8-hyperoval-q94-mode11-pilot.md) has a height-zero
letter (`min(heights)=0`) and nonnegative recorded orders. Its cutoff-47
pilot pressure is `0.709928941871497`, above `1/sqrt(2)` by about
`0.002822`. By the lemma, the complete q94 pressure is at least that
value, so the q94 drop cannot cross. A full q94 profile would compute a
larger number, not a crossing; it is unnecessary for the crossing
question.

## Scope

This is an analytic monotonicity statement about a fixed map and a
resolved subset of its orders. It does not compare different maps, does
not supply upper bounds for unresolved tails, and does not change the
factoring exponent. The numerical pilots cited remain finite checks.

## Reproduction

```sh
make verify-truncated-profile-lower-bound
```

The checker validates the q94 map, checks the height-zero letter and
the monotonicity inequality on a grid, and replays the recorded pilot
pressure and its defining equation.
