# q96 unequal cut scales: a boundary-layer screen

## Result and scope

The screened complete 96-letter F8 profiles stay above the asymmetric
supporting-line threshold under unequal cut scales. A deterministic
screen of six completed profiles finds the best retained interior
support value at

\[
R=1.0037692354525074
\]

on the refined `789` retraction. The corresponding balanced value is
`1.0040998136`, so the scale asymmetry removes only about `3.3e-4`
from the support. Every screened interior value is strictly above one.

The screen also exposes the previously reported near-one values. As
one supporting weight tends to zero, the support decreases toward one
from above. This is the degenerate one-cut boundary layer, not a
fixed-aspect crossing. The published moving-frame best
`1.0004187269519174` occurs at supporting coefficient
`cos(1.5697963267948964) = 0.000778`, in the same layer.

No global optimum, exact pressure certificate, or stronger factoring
score is claimed. The complete exponent remains `1/2`.

## 1. Support equation

For a complete two-cut profile write

\[
U(s)=\sum_i e^{-s\kappa_{\infty,i}},\qquad
V(s)=\sum_i e^{-s\kappa_{0,i}},
\]

with integer raw orders and physical heights `h_j`. For positive
`a,b` with `a^2+b^2=1` and intrinsic scale `t>0`, set

\[
\Phi(t,a,b)=a\log U(t/a)+b\log V(t/b).
\]

The two paired legs on each side carry the same height list, so the
Gibbs physical charges are `G_X=G_Y=2th`. The common type is therefore
also the optimal independent paired type. The support value `R(t,a,b)`
is the positive root of

\[
\Phi(t,a,b)=R\left(\log\sum_j e^{-2th_j/R}
                    +\log\sum_j e^{-2th_j/R}\right).       \tag{1}
\]

At `a=b=1/sqrt(2)` the value `R` is exactly `sqrt(2)` times the natural
balanced pressure `rho` of the
[supporting-line normalization](asymmetric-support-line-crossings.md).
The screen checks this against the certified or recorded pressure of
every displayed profile.

## 2. Complete profiles and calibration

The six inputs are the completed q96 profiles already in the
repository:

| profile | recorded/certified `rho` | balanced `R/sqrt(2)` | interior `R` | angle | `(a,b)` | scale |
|---|---:|---:|---:|---:|---:|---:|
| `f8_hyperoval_full` | 0.710005835 | 0.7100058347 | 1.0037757401 | 1.4417958661 | (0.1286430, 0.9916910) | 0.1684814 |
| `coherent_four_boundary` | 0.7100058653 | 0.7100058653 | 1.0037792126 | 1.4414662072 | (0.1289699, 0.9916485) | 0.1686107 |
| `coherent_coordinate` | — | 0.7100058749 | 1.0037809117 | 1.4413289147 | (0.1291060, 0.9916308) | 0.1686636 |
| `coherent_hot_f4_hyperoval` | 0.7100058606 | 0.7100058606 | 1.0037791260 | 1.4414761581 | (0.1289600, 0.9916498) | 0.1686067 |
| `f8_retraction_678` | 0.7100057890 | 0.7100057890 | 1.0037699802 | 1.4423150079 | (0.1281281, 0.9917576) | 0.1682764 |
| `f8_retraction_789` | 0.7100057872 | 0.7100057872 | **1.0037692355** | 1.4423973626 | (0.1280465, 0.9917682) | 0.1682433 |

The balanced calibration error is `3.1e-10` for the first profile and
at most `1.5e-14` for the other five. The first recorded value is a
certified strict upper bound, so the numerical balanced minimum may lie
slightly below it. The two values agree to the displayed precision.

Swapping the two cuts maps the retained angle to its mirror and leaves
the support unchanged to numerical precision. The table keeps the
representative with the smaller coefficient at `a`.

## 3. The one-cut boundary layer

The screen keeps `a,b >= 0.05` for the interior comparison. At smaller
angles the value decreases toward one. For the `789` profile:

| angle | support |
|---:|---:|
| 0.03 | 1.004634583 |
| 0.01 | 1.002980157 |
| 0.003 | 1.000914015 |
| 0.001 | 1.000306548 |
| 0.0003 | 1.000092161 |
| 0.0001 | 1.000030739 |

The excess `R-1` is approximately `0.307` times the angle over these
points. It approaches zero only as one cut is switched off and the
intrinsic scale also tends to zero. The
[geometric recurrence](asymmetric-support-line-crossings.md) needs a
fixed aspect ratio whose two sides both grow. A sequence with `a -> 0`
has orientation denominator tending to one and cannot supply a uniform
`R<1` gap for any `alpha>1/2`. The near-one values are therefore a
boundary artifact, not a proposed crossing.

The same reading applies to the best retained moving-frame support
`1.0004187269519174`: the published record has supporting coefficient
about `7.8e-4` and the constant masks `a=b=0`, so it sits in the
degenerate layer. This note does not reprove that record's global
scope; it records the common normalization and the boundary behavior.

## 4. Reproduction and limits

```sh
python3 tools/screen_q96_asymmetric_scales.py \
  --output evidence/q96-asymmetric-scale-screen.json
make verify-q96-asymmetric-scale
```

The checker replays the full deterministic grid and the boundary
sequence, verifies the profile hashes, the balanced calibration, and
the strict positivity of every retained interior support. The
[evidence](../evidence/q96-asymmetric-scale-screen.json) records every
profile hash, angle, scale, coefficient pair, and boundary point.

The screen is finite. It uses complete profiles but only the two-place
supporting-line criterion with a common type. It does not exclude
other gauges, mixed-place allocations, other profiles, or a different
coefficient construction. Numerical refinement is not a global
optimality proof, and no exact rational pressure certificate is
supplied. The factoring exponent remains `1/2`.
