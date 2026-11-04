# Bounds on unresolved orders extend the protected-core barrier

## Results

For constant height-homogeneous paired splits that retain every coordinate
word of height below 24, over any characteristic-two field:

- Support through height **37** forces balanced natural `rho >= 17/24`.
- Support through height **38** forces balanced natural `rho >= 29/41`.

Both constants exceed `1/sqrt(2)`. All hot modes are allowed. The two
paired alphabets and their physical types may differ. The high map spaces
need not be deletion closed. The assumptions are literal protection of
the low coordinate block, `LP=I`, constant height-preserving maps, and the
common balanced natural cut scale.

This extends the [81-state height-36 result](height36-protected-core-floor.md).
It does not exclude changes to the protected low block, higher supports,
nonconstant maps, other physical gauges, or unequal cut scales and
supporting angles. The height-38 conclusion does not extend to height 39:
the same test admits `70/99 < 1/sqrt(2)` there; see the
[height-39 transition](protected-core-height39-transition.md). The
factoring exponent remains `1/2`.

## 1. Recover a lower bound from the unknown tail

The protected core has 81 words and 6,561 channels. Its exact retained
profile gives 3,359 orders through 64, whose sum is 173,663. The remaining
3,202 orders are unknown. No individual values are assigned to them.

Every entry of the binary one-particle hot matrix has degree at most the
sum of its row and column mode heights. A selected compound entry has
that same bound with the heights summed over its occupied modes.
For the core channel matrix, each row is a pair of physical words.
Bounding each determinant permutation by its row and column height sums
therefore gives

```math
\deg\det M\le 2\sum_{x,y}(h_x+h_y)
 =4\cdot81\sum_x h_x=433836,
\qquad \sum_xh_x=1339.
```

The endpoint proves that the determinant is nonzero. The sum of its local
Smith orders is its order of vanishing and is at most its degree.
Thus the unknown orders have sum at most

```math
433836-173663=260173.
```

For `0<r<1`, convexity of `r^x` gives an aggregate Jensen bound

```math
\sum_{\rm unknown}r^{\nu_i}
\ge3202\,r^{260173/3202}
\ge3202\,r^{82}.
```

Define the lower polynomial

```math
\widetilde Z(r)=\sum_{\rm known}r^{\nu_i}+3202r^{82}.
```

The term at 82 is **not a claim that any unknown order equals or is
bounded individually by 82**. Only their total degree budget and the
resulting partition lower bound are used. This distinction is essential.

## 2. Protect literal tensors rather than the high map spaces

Use the literal-submatrix argument from the
[general height-25 protected-state theorem](coherent-protected-height25.md).
Only the fixed low core needs deletion closure. Constant height-preserving
maps commute with the natural diagonal normalizations, and they act as
identities on that core. Both larger regular cuts therefore contain the
same fixed literal restrictions, even when the high map spaces are not
closed under deletion.

The protected vacuum entry also ensures that removable common scalar
content is a monomial, up to a field unit. Its opposite local shifts
cancel in the two-cut product. The binary core profiles are unchanged
by characteristic-two field extension. The protected-submatrix inequality
and the Jensen bound give

```math
Z_0(r)Z_\infty(r)\ge\widetilde Z(r)^2.
```

This inherits the literal-tensor and scalar-shift argument from the cited
work, rather than assuming that an arbitrary high split has the same
profile as its compressed compound matrix.

## 3. All-mode physical entropy budgets

For a height cap `H`, every allowable word uses odd mode heights at most
`H`. The maximal physical height polynomial is

```math
Q_H(t)=\left[\prod_{\substack{1\le j\le H\\j\text{ odd}}}(1+t^j)\right]_{\deg\le H}.
```

The two capacities are `Q_37(1)=372` and `Q_38(1)=409`.
The rank condition `LP=I` bounds each physical height-cell dimension by
its coordinate capacity, so both physical partition sums are at most
`Q_H(t)`. An ambient reflection interval of size 64 suffices for these
finite height ranges.

For independently chosen paired types, their two Gibbs inequalities give

```math
H_X+H_Y+\frac{\log r}{\rho_0}(E_X+E_Y)
\le2\log Q_H(r^{1/\rho_0}).
```

Combining this with the channel lower bound proves the stated pressure
floor whenever

```math
b\log\widetilde Z(z^a)-2a\log Q_H(z^b)>0,
\qquad \rho_0=a/b,
\qquad 0<z<1.
```

## 4. Exact global certificates

For `(H,a,b)=(37,17,24)` and `(38,29,41)`, the checker constructs the
integer derivative-sign polynomial

```math
G(z)=\sum_{i,j}(i-2j)c_i d_j z^{ai+bj-a},
```

where `c_i,d_j` are the coefficients of `Z_tilde,Q_H`.
Its degrees are 2,265 and 3,907. Each transformed polynomial
`(1+x)^degree G(x/(1+x))` has exactly two coefficient sign variations.
Exact rational endpoint evaluations isolate a positive-to-negative root
and a negative-to-positive root, exhausting the roots in `(0,1)`.

The latter is the only interior minimum. Outward 192-bit logarithm
intervals prove a positive gap throughout its isolating interval in
both cases. The limiting gap at zero is zero, and endpoint positivity
is checked by `6561^b > Q_H(1)^(2a)`. The integer tests
`2*17^2 > 24^2` and `2*29^2 > 41^2` establish that both floors exceed the
factoring threshold.

Run:

```sh
python3 tools/check_jensen_tail_core_floor.py
```

The [receipt](../evidence/jensen-tail-core-floors.json) contains all degree
budgets, root intervals, sign hashes and positive gap numerators.
The check uses the source-pinned
[low-core data](../evidence/height24-core-low-orders.json); it does not
reconstruct the native core matrices. Those can be rebuilt with
`python3 tools/check_height36_core_floor.py --replay-core`.
No missing Smith orders are treated as known, and no Lean claim is made.
