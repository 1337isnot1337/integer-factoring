# A uniform cubic pressure floor from the endpoint alone

## Claim and scope

For every odd-characteristic field, let a degree-three polynomial path
in the [Bell local space](bell-diagonal-crossings.md) have the exact
crossing endpoint. Suppose its constant coefficient has **row rank one**
and its leading coefficient has **column rank one**. At every physical
dimension `d=2^k>=8`, its order/degree pressure parameter is strictly
greater than `5/7>1/sqrt(2)`.

This excludes the displayed sufficient criterion for those cubic targets,
without requiring self-reciprocity, prescribed middle coefficients, or
a finite search. It does not lower-bound actual tensor-power coefficient
ranks or arbitrary coefficient recombinations. Quartics, different outer
ranks, and sharper actual-rank estimates remain open. The supported
factoring exponent stays `1/2`.

The argument uses the exact Bell cut transform and the complete
[polynomial pressure contract](bell-polynomial-synthesis.md). It is ordinary
mathematics with substantial OpenAI Codex assistance; independent review
and any Lean formalization are not claimed.

## 1. Endpoint forcing in both cuts

Write `q=d^2` and

```math
s(z)=a_0+a_1z+a_2z^2+a_3z^3,\qquad s(1)=\boldsymbol1,
\qquad b_i=Ta_i.
```

Here `T` is the regular signed Hadamard Bell transform, `T^2=I` and
`T1=1`. The rank hypotheses say that `a_0` and `b_3` each have exactly
one nonzero coordinate. Each column of `T` has entries `+1/d,-1/d`,
with respective integer counts `(q+d)/2,(q-d)/2`.

Since `a_3=Tb_3`, its entries are `+alpha,-alpha` for a nonzero scalar
`alpha`. If `alpha` is neither `1` nor `-1`, every entry differs from
one. Otherwise at least `(q-d)/2` entries differ from one. Discarding
the single coordinate supporting `a_0`, at least

```math
M=(q-d)/2-1
```

remaining coordinates have `a_0=0` and `a_3!=1`. The endpoint then forces
`a_1+a_2=1-a_3!=0`. Each of these coordinates has row order at most two.
The support of `a_0` itself has order zero.

In the column cut, `b_0=Ta_0` has the same two-valued property.
Outside the singleton support of `b_3`, at least `M` coordinates have
`b_0!=1`; hence `b_1+b_2!=0`. Each has column degree at least one,
while the singleton support of `b_3` has degree three.

All statements hold in every odd characteristic, including extensions.
The two signs are distinct; no numerical norm or positivity of field
elements is used. Only integer dimensions are counted.

Exchanging the two cut directions does not bypass the argument: `b_0`
and `a_3` have no zero coordinates, so every opposite-cut order is zero
and every opposite upper degree is three. That orientation has pressure
`q^2 x^3>=q^2`, hence parameter at least one.

## 2. Every scale has too much pressure

For `x>=1`, let `A(x)` sum `x^(-row_order)` over physical eigenvalues,
and let `B(x)` sum `x^(column_degree)`. The preceding disjoint sets give

```math
A(x)\ge1+M/x^2,\qquad
B(x)\ge x^3+Mx,
\qquad Q(x)=A(x)B(x)\ge x^3(1+M/x^2)^2.
\tag{1}
```

The minimum over all positive `x` of the right side is

```math
16(M/3)^{3/2}.
\tag{2}
```

For example, minimize its positive square root
`x^(3/2)+M*x^(-1/2)`; its derivative vanishes at `x^2=M/3`.
Using all positive scales only weakens the lower bound, so (2) is valid
for the permitted interval `x>=1` as well.

For `d>=8`,

```math
64M-27d^2=(d-8)(5d+8)\ge0.
```

Thus (1)--(2) imply

```math
Q(x)\ge\frac{27}{32}d^3>d^{20/7}.
\tag{3}
```

The strict comparison is exact: at the smallest allowed `d=8` it reduces
to `3^21>2^32`, which follows from `3^7>2^11`. The ratio only increases
with `d`. Therefore

```math
\rho_{\rm pressure}(x)=\frac{\log_2 Q(x)}{4\log_2d}
>\frac57>\frac1{\sqrt2}.
```

At `d=8`, the explicit relaxed floor is `Q>=432` and
`rho>=log2(432)/12≈0.729574`. This is a floor for the pressure upper-bound
method, not a lower bound on the actual evaluation time or actual ranks.

## 3. Why quartics remain a different target

The same endpoint reasoning for degree `h>2` gives

```math
Q(x)\ge x^h(1+M/x^{h-1})^2.
```

Its optimized leading contraction floor as `d` grows is
`h/(2(h-1))`: `3/4` for cubics but only `2/3` for quartics. The quartic
value is below the critical threshold, so this argument does not close
that target. It supplies no quartic construction and no stronger bound.

The restarted synthesis therefore keeps quartic targets and cubics that
escape the stated outer-rank conditions in scope. Improving cubic middle
coefficients while retaining these two outer ranks cannot succeed under
this pressure criterion.

## 4. Reproduction

```sh
python3 tools/check_bell_cubic_floor.py --output evidence/bell-cubic-outer-rank-floor.json
```

The [exact replay](../evidence/bell-cubic-outer-rank-floor.json) constructs
216 literal endpoint-valid cubic paths in characteristics `3,5,7`, checks
both forced support counts and 864 rational pressure comparisons, and
checks the dimension inequality through `d=2^32`. The infinite-dimensional
and all-field assertions are the proof above, not a finite extrapolation.
