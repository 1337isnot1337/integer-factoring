# Exact valuation-separated jets at integer centers

## Constructive theorem

The [global packet circuits](interval-packet-circuits.md) make modular
input polynomials inexpensive, but their large powers of two can hide
all useful coefficients at modest raw precision. This note computes
that local content separately, keeping its exponent in binary.

Let `F=P_I` or `F=Q_(N,J)` with `N>=1` from the
[interval orbit-gcd interface](interval-orbit-gcd.md). Given an integer
center `x`, precision `p>=1`, and jet degree `d>=0`, there is a
deterministic algorithm that computes integers `k,E,z` and the jet

```math
U(Z)=\frac{F(x+2^kZ)}{2^E Z^z}
\pmod{(2^p,Z^{d+1})},\qquad U(0)\text{ odd}.        \tag{1}
```

The quotient before reduction is an integer polynomial. Here `z` is
the exact order of the zero at the center, and `E` is the exact
2-adic content after removing that zero. For these products `z` is
zero or one. The exponent `E` is stored in binary and can be as large
as the numerical interval length. If all input integers have `O(n)`
bits, then `k=O(n)`, `E` has `O(n+log n)` bits, and the complete bit
cost is polynomial in `n,p,d+1`.

The algorithm evaluates a specified center. It does not search all
centers or extract a global common factor. The strongest complete
factoring score remains `1/2`.

## 1. A scaled affine progression

Both substitutions have the form

```math
\prod_{t=L}^{U-1}\bigl(a+bt+2^k(c+et)Z\bigr),
```

with `(a,b,c,e)=(x,-1,1,0)` for `P`, and
`(N,-x,0,-1)` for `Q`. Choose

```math
k=\max\bigl(1,\operatorname{bitlength}
  (\max(|a+bL|,|a+b(U-1)|))\bigr).
```

Every nonzero constant factor has valuation strictly below `k`,
because an affine function's maximum absolute value on an interval
is attained at an endpoint.

If a constant factor is zero at `t_0`, extract

```math
2^k(c+et_0)Z
=2^{k+v_2(c+et_0)}Z\cdot\operatorname{oddpart}(c+et_0).
```

The remaining indices form at most two intervals. The coefficient
`c+et_0` is nonzero for the present `P,Q` products. Handle the empty
product as one. Thus it suffices to treat an interval with no constant
zeros.

Remove `g=min(v_2(a),v_2(b))` from every factor, using the nonzero
argument when one of `a,b` is zero. Replace `a,b` by their quotients
and `k` by `k-g`, and add `g*(U-L)` to `E`. At least one of `a,b`
is now odd, and the new `k` is still greater than all remaining
constant valuations.

## 2. Only `O(n)` valuation classes are needed

If `b` is even, `a` is odd, so all factors already have odd constant
term. Their base and progression step are jets

```math
A(Z)=a+2^kcZ,\qquad B(Z)=b+2^keZ.
```

The base is a unit and every coefficient of `B` is even.

If `b` is odd, valuation `v` occurs exactly at

```math
t\equiv r_v=(2^v-a)b^{-1}\pmod{2^{v+1}}.
```

The possible `v` are bounded by the bit length of the largest
nonzero affine constant, hence there are `O(n)` of them. Set
`t=r_v+2^(v+1)w`; the indices `w` form an ordinary interval obtained
by ceiling division. Divide this bucket's factors by `2^v`. They
become

```math
A_v(Z)+wB_v(Z),
```

where

```math
A_v(Z)=\frac{a+br_v}{2^v}+2^{k-v}(c+er_v)Z,
\qquad
B_v(Z)=2b+2^{k+1}eZ.
```

The base constant is odd and every coefficient of the step is even.
Add `v` times the bucket's cardinality to `E`. The algorithm checks
that the bucket cardinalities sum to the interval length, without
visiting their individual elements.

## 3. Unit progression products have short exact jets

For a unit jet `A`, an even jet `B`, and `w` in an integer interval
of length `ell`, the finite product identity gives

```math
\prod_w(A+wB)
=A^\ell\sum_{j=0}^{\min(\ell,p-1)}e_j(w)(BA^{-1})^j
\pmod{2^p}.                                     \tag{2}
```

Higher terms vanish because `B` is even coefficientwise. The
[unit-jet algorithm](interval-unit-jets.md) computes the necessary
`e_j` from exact integer power sums and Newton identities. It pays
only polynomial bit cost in the binary endpoints and `p`. Jet
inversion and multiplication cost polynomial time in `p,d+1`.

Multiply (2) over the `O(n)` buckets and the extracted zero factor.
Every resulting unit has odd constant term, so their product does
too. This proves both the exact content statement and the full
algorithmic bound in (1). At no point is a ring of raw precision
`E+p` allocated; only `E` itself is stored as a binary integer.

## 4. Remaining global obligation

This supplies effective normalized local scalars at integer centers,
including centers where one product has a genuine zero. It avoids
paying the numerical size of the enormous local content. It does
not identify a favorable center: the explicit scale can describe a
very small 2-adic neighborhood, and finding the actual common
integer roots remains the orbit-elimination problem.

In particular, treating these local unit jets as a global factor
certificate would repeat the local-divisibility issue in the earlier
unit-point construction. A global algorithm must retain polynomial
degree information or another exact certificate while using these
normalized evaluations.

## 5. Evidence and replay

```sh
python3 tools/check_normalized_interval_jets.py --output /tmp/normalized-jets.json
cmp /tmp/normalized-jets.json evidence/normalized-interval-jets.json
make verify-arithmetic
```

The [implementation](../tools/normalized_interval_jets.py) returns the
binary content, zero order, scale, and unit coefficients. The
[checker](../tools/check_normalized_interval_jets.py) compares all
three with full dense integer polynomials in 1,153 finite controls,
including a true factor at `N=143,x=11`. It then verifies six
long-range split identities with the exact scale adjustments.

Two factorial controls independently check the valuation by Legendre's
formula and the unit by a fixed-precision residue-period table. Two
more controls check `Q`'s valuation using cumulative divisibility
counts instead of exact-valuation buckets. The
[receipt](../evidence/normalized-interval-jets.json) includes an
interval of length `2^200+12` with a 201-bit content exponent; its
unit jet uses only 12 bits of precision. These are exact local data,
not a factoring timing claim.
