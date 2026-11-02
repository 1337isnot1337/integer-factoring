# The height-38 protected-core barrier does not extend to height 39

## Exact result

The [height-38 Jensen-tail floor](jensen-tail-protected-core-floor.md) is
`29/41 > 1/sqrt(2)`. At height 39 the same necessary condition no longer
excludes the subthreshold value

```math
\rho_0=\frac{70}{99}=0.707070707\ldots<1/\sqrt2,
\qquad 2\cdot70^2=9800<99^2=9801.
```

The comparison is exact. Let

```math
\widetilde Z(r)=\sum_{\text{known}}r^{\nu_i}+3202\,r^{82},
\qquad
Q_H(t)=\left[\prod_{\substack{1\le j\le H\\j\text{ odd}}}(1+t^j)\right]_{\deg\le H}.
```

At `z=899/900`, `a=70`, `b=99`, the test value
`Z_tilde(z^a)^b - Q_H(z^b)^(2a)` is positive for `H=37,38` and negative
for `H=39,40`. Positive means the Jensen-tail necessary condition excludes
`rho_0`; negative means it does not. The capacities are
`Q_37(1)=372`, `Q_38(1)=409`, `Q_39(1)=450`, `Q_40(1)=496`.

Thus the height-38 conclusion is not inherited by height 39. This is a
statement about the lower-bound argument, not a construction: no split
with `rho=70/99` is known, and the true height-39 floor may be lower. What
fails is the claim that every support-through-38 protected split is above
`1/sqrt2`; at height 39 the same test admits a value below it.

## What this changes

The best explicit two-cut constructions, including the complete q96
profiles, use support height 25. The height-38 floor applies to them:
retaining the height-below-24 core and supporting only through height 38
forces balanced natural `rho >= 29/41 > 1/sqrt2`. Refining a height-25
retraction changes the boundary vectors and the last decimals of the
pressure, but the height-38 floor still forces balanced natural
`rho >= 29/41`; such a retraction cannot cross `1/sqrt2`. The observed
plateau near `0.7100` is therefore structural, not a search-resolution
issue.

An escape within the protected-core family needs support at least 39, or
a change to the protected low block itself. A height-39 map would still
have to satisfy the literal core, constant height-homogeneous maps, the
common balanced cut scale, and the full uniform transfer. The exact test
here only says that the height-38 barrier no longer rules it out.

## Scope

The assumptions are those of the parent note: literal protection of the
height-below-24 coordinate core, constant height-homogeneous paired maps,
`LP=I`, independent paired types, all hot modes allowed, and a common
balanced natural cut scale. The comparison is a single exact rational
point for a necessary condition. It is not a global optimum, a split, an
exact pressure certificate, or a score change. The complete factoring
exponent remains `1/2`.

## Reproduction

```sh
python3 tools/check_protected_core_height39.py \
  --output evidence/protected-core-height39-transition.json
make verify-protected-core-height39-transition
```

The [receipt](../evidence/protected-core-height39-transition.json) records
the capacities, the exact sign pattern, side hashes, the core-receipt
hash, and source pins. The checker recomputes the integer comparison and
does not treat the unresolved tail orders as known.
