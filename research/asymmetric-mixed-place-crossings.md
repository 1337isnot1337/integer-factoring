# Supporting lines from opposite flags at several places

## Result and scope

The [opposite Hasse-jet flags](mixed-place-hermite-crossings.md) support
an asymmetric bound on the two cut widths. Each cut can use several
places, independent physical gauges, and its own reconstruction
allocations. Every place must use the same pair of physical types.

This extends the [two-place supporting-line criterion](asymmetric-support-line-crossings.md).
The formula includes the full polynomial degree exactly once. It also
counts every conjugate root as a separate reconstruction place. It is
a sufficient criterion; no new seed in this note satisfies the strict
threshold needed to improve the factoring exponent.

## 1. Complete local costs

Let a fixed polynomial crossing tensor have degree `D`, equal physical
alphabet size `q`, and strict endpoint `[x=v][y=u]` at one. Choose
disjoint nonempty row and column place sets, excluding one. All places
are rational over one fixed coefficient extension. At infinity use
the homogeneous section `tau^D K(1/tau)`.

For a cut and a place `j`, write

```math
\lambda_{j,i}(s,h),\qquad 1\le i\le q^2,
```

for its complete weighted Smith orders. An entry with ordinary local
valuation `v` and physical symbols `x,y,u,v'` has weighted valuation

```math
s\,v+h_x(x)+h_y(y)+h_u(u)+h_v(v').
```

Thus `s>0` is the local logarithmic scale, and `h` contains the already
scaled physical heights. The orders are positively homogeneous:

```math
\lambda_{j,i}(c s,c h)=c\lambda_{j,i}(s,h)\qquad(c>0).       \tag{1}
```

For rational `h/s`, these are ordinary Smith orders after fixed
ramification, multiplied by `s`. Negative orders are permitted.
Define the complete local logarithmic channel sum

```math
\mathcal L_j(s,h)=\log\sum_{i=1}^{q^2}
                   \exp(-\lambda_{j,i}(s,h)).               \tag{2}
```

Use rational paired types `pi_X,pi_Y` on the `x,v` and `y,u` words.
Put

```math
H=H(\pi_X)+H(\pi_Y)>0,\qquad
\mu=(\pi_X,\pi_Y,\pi_Y,\pi_X).                            \tag{3}
```

The type shift at a place is exactly `mu·h`. It cannot be optimized
with different types at different places.

## 2. The mixed supporting-line formula

Fix positive `a,b` with `a^2+b^2=1`, and one intrinsic scale `t>0`.
On each cut, group any conjugate places that have identical complete
profiles under the chosen physical gauge. A group `j` contains `e_j`
distinct places. Give it positive allocation `alpha_j`, with
allocations summing to one separately on each cut. Each root receives
allocation `alpha_j/e_j`.

For the row cut set `a_R=a`; for the column cut set `a_C=b`.
Choose an independent branch-weighted physical height vector `b_j`
at every group. Its actual height in (2) is `b_j/a_cut`. Define

```math
F_j=a_{\rm cut}\,
\mathcal L_j\!\left(\frac{t e_j}{a_{\rm cut}\alpha_j},
                    \frac{b_j}{a_{\rm cut}}\right)
       +\mu\cdot b_j.                                     \tag{4}
```

Then the opposite-flag construction gives square crossing pieces
whose normalized cut widths satisfy

```math
\boxed{
ar+bc\le
\frac{D t+\max_{j\in R}F_j+\max_{j\in C}F_j}{H}+o(1).
}                                                          \tag{5}
```

The width normalization is the one in the supporting-line note:
`r=log(rank_R)/(2 log Q)` and `c=log(rank_C)/(2 log Q)` for the final
square physical alphabet of size `Q`.

To prove (5), let the tensor-power degree be `Dm`. A piece with row
flag prefix `k+1` has a representation by initial jets at every row
place. At a root in group `j`, the largest needed jet index is
at most `(alpha_j/e_j)k+O(1)`. Use local scale
`t e_j/(a alpha_j)` and local heights `b_j/a`. The complete
typed Smith bound gives

```math
a\log\operatorname{rank}_R A_{m,k}
 \le t k+m\max_{j\in R}F_j+O(\log(m+1)).                  \tag{6}
```

The opposite column flag has complementary prefix `Dm-k+1`, so

```math
b\log\operatorname{rank}_C A_{m,k}
 \le t(Dm-k)+m\max_{j\in C}F_j+O(\log(m+1)).               \tag{7}
```

The two representations define the same piece. Adding (6) and (7)
cancels the index and produces `Dt` once. Pairing the restricted path
with its physical pair swap gives equal alphabets, retaining the two
piece indices separately. Its alphabet logarithm is `mH+O(log(m+1))`.
The doubled rank bound cancels the factor two in the normalized
widths, giving (5).

The flags and coefficient reconstruction are effective over the fixed
extension, including in positive characteristic. The prior Hasse-flag
proof supplies the stable basis construction and polynomial overhead.
Admissible tensor powers clear all rational type, allocation, and
height-ramification denominators.

A strict rationally certified right side below one can be combined
with the [complete asymmetric transfer](asymmetric-modular-transfer.md).
That conclusion still requires the exact path, complete local
profiles, common types, and effective all-width factors. A floating
search value alone supplies none of those missing certificates.

## 3. Normalization and angle derivative

At `a=b=1/sqrt(2)`, divide the balanced mixed-place intrinsic scale
and all already scaled heights by `sqrt(2)`. The right side of (5)
is exactly `sqrt(2)` times the balanced parameter `rho`.

For one group, write `a_cut=c(theta)`. By (1), the first term in
(4) is the perspective

```math
c\log\sum_i\exp\left[-\lambda_{j,i}
              \left(\frac{t e_j}{\alpha_j},b_j\right)/c\right].
```

The physical mean `mu·b_j` is independent of `c`. Within a fixed
attaining-minor region, its angle derivative is

```math
c'(\theta)\left(\log Z+\sum_i p_i\lambda_i/c\right),       \tag{8}
```

where `p_i=exp(-lambda_i/c)/Z`. The expression in parentheses is
the channel entropy. A selected attaining region supplies the
corresponding derivative away from ties. The search uses a smooth
upper approximation to each maximum; it records the unsmoothed value
separately.

The [implementation](../tools/asymmetric_mixed_place_pressure.py) uses
`b_j=t B_j` to keep normalized native valuation windows controlled.
It optimizes common types, local gauges, allocations, the intrinsic
scale, and the supporting angle jointly. The [explicit-seed worker](../tools/search_asymmetric_explicit_places.py)
accepts complete polynomial tensors and saves every completed trial
atomically. Its outputs remain numerical search evidence.

## 4. Independent controls

Run

```sh
python3 tools/check_asymmetric_mixed_place_pressure.py \
  --output evidence/asymmetric-mixed-place-pressure-controls.json
```

The [checker](../tools/check_asymmetric_mixed_place_pressure.py) uses
the original reflection seed and two different affine deformations.
Its [evidence](../evidence/asymmetric-mixed-place-pressure-controls.json)
contains nine balanced-normalization identities, nine comparisons
against the independent native two-place Gibbs implementation, and
348 finite-difference derivative coordinates. The controls include
conjugate groups and both height parameterizations.

These checks validate the objective and search interface. They do not
establish a global optimum, a rational pressure certificate, or an
improved factoring score.
