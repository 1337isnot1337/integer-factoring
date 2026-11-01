# Rational finite adaptive bounds without channel convolution

## Result and scope

The [finite typed adaptive certificate](finite-typed-adaptive-certificates.md)
can use exact rational coefficient majorants in place of the full channel
convolution. Every reported jet bound remains an integer upper bound,
including the first positive jet after a free prefix. A finite catalogue
may select both a physical gauge and a rational coefficient parameter
independently for each jet.

This is a computational alternative to the exact coefficient sum. It
cannot improve that sum itself. Its purpose is to evaluate large finite
blocks and the exact amplification gate without expanding all channel
count coefficients. No threshold-crossing block or new factoring
exponent is supplied here.

## 1. Positive coefficient identity

Use a certified local profile with integer ramification `L`, integral
orders `nu_i`, exact mean shift `A`, and tensor length `m` with `mA`
integral. Put

```math
b=\min_i\nu_i,\qquad
P(x)=\sum_i x^{\nu_i-b},\qquad
Q=Lj+m(A-b).                                         \tag{1}
```

All exponents of `P` are nonnegative integers. The exact channel bound
from the weighted catalogue theorem is

```math
c_g(j)=[x^Q]\frac{P(x)^m}{(1-x)^2}                    \tag{2}
```

when `Q>=0`, and is zero when `Q<0`. Indeed, the coefficient of `x^k`
in `(1-x)^-2` is `k+1`. Each channel word therefore contributes exactly
`max{Q-sum(nu_i-b)+1,0}`.

For any rational `0<t<1`, nonnegativity of every coefficient gives

```math
c_g(j)\le t^{-Q}\frac{P(t)^m}{(1-t)^2}.               \tag{3}
```

Thus, for the exact typed cut dimension `d`, a valid integer jet bound is

```math
B_{g,t}(j)=
\begin{cases}
0,&Q<0,\\
\min\left\{d,
\left\lceil t^{-Q}P(t)^m/(1-t)^2\right\rceil\right\},&Q\ge0.
\end{cases}                                          \tag{4}
```

The zero test in (4) is exact. In particular, at `Q=0` the expression is
strictly positive. It never replaces the first positive jet by the limit
from the free side of a capacity curve.

## 2. Monotonicity and implementation

For fixed gauge and `t`, increasing `j` by one multiplies the uncapped
rational expression by `t^-L>1`. Its ceiling and dimension cap preserve
monotonicity. Taking a minimum over any fixed finite collection of
gauges and rational parameters also preserves monotonicity. The same
greedy initial-jet merge and opposite-flag reconstruction therefore
apply to (4).

The implementation computes `P(t)^m/(1-t)^2` exactly once. It starts at
the first nonnegative `Q`, and advances the numerator and denominator
by the fixed integer powers of the numerator and denominator of `t`.
Once the type-dimension cap is reached, all later bounds equal that cap.
It uses no floating coefficient or rank estimate.

For a fixed local profile and rational parameter, the integers involved
have `O(m)` bits and there are `O(m)` jets. The arithmetic and storage
are polynomial in `m`. No channel coefficient array is constructed.
The final finite rank table can be passed directly to the sharp
sum-of-square-roots amplification gate.

The physical bases, types, place distinctness, and complete local
profiles have exactly the same independent proof obligations as in the
finite typed wrapper. A numerical local optimizer is not a substitute
for those inputs.

## 3. Exact controls

[`rational_typed_adaptive_bounds.py`](../tools/rational_typed_adaptive_bounds.py)
provides the local bound, catalogue merge, and full finite certificate.
The checker verifies 1,824 comparisons with exact channel convolutions
and 1,824 translated-gauge identities. It includes negative orders,
three ramification values, small and large exact dimension caps, and
four rejected coefficient parameters.

Its free-prefix control has two channels of order one, tensor length
three, and `t=1/2`. The resulting bounds begin
`0,0,0,32,64,128`. The first positive jet is retained.

Four actual ternary tensor controls, including two-by-three and
three-by-two paired alphabets, recompute 24 complete local profiles.
They reconstruct all opposite-flag pieces, verify 52 complete cut
ranks, and check 1,960 endpoint entries. Every bound dominates its
corresponding actual rank.

Finally, an exact length-600 uniform ternary block produces all 1,201
piece rank pairs. Each external type class has a 942-bit cardinality,
and the largest rank-product upper bound has 3,159 bits. The complete
integer calculation took about 0.33 seconds in the recorded run. Its
threshold test fails, as expected; the record retains its exact rank
table digest rather than a large repeated decimal table.

```sh
python3 tools/check_rational_typed_adaptive_bounds.py \
  --output evidence/rational-typed-adaptive-bound-controls.json
```

The [checker](../tools/check_rational_typed_adaptive_bounds.py) and
[replay record](../evidence/rational-typed-adaptive-bound-controls.json)
bind every source used in these controls.
