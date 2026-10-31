# Uniform exclusion of balanced gauges above the logarithmic scale

## Result and scope

The [uniform entry certificates](rectangular-uniform-ray-obstruction.md)
exclude more than fixed asymptotic rays. The gauge vector may vary
arbitrarily with the interval length and the profile scale.

**Theorem.** Fix `c>1/sqrt(2)`. For all sufficiently large power-of-two
`M`, every positive profile scale `ell=log s` satisfying

\[
\tau=M\ell\ge c\log M                                      \tag{1}
\]

and every real 24-component gauge vector `theta` satisfying

\[
\theta_0+\theta_1=\theta_2+\theta_3=0                       \tag{2}
\]

have strictly positive critical pressure for the mixed rectangular
seed with `alpha=1`, `beta=omega`, its original quartic numerator,
and its partial output repair.

The theorem covers the full 25-parameter family and its smaller
19-parameter restriction. The core weights need not be bounded,
converge, or follow an affine ray. There is no upper bound on `ell`.
The eventual threshold depends on `c` and the fixed seed, but not on
the gauge vector or profile scale.

The result remains restricted to balanced aggregate bulk charges.
It does not cover a different repair, the full four-by-four physical
alphabet, small isolated interval lengths, or scales below (1).
It does not improve the factoring exponent. Its purpose is to remove
a large asymptotic search region from this specific construction.

## 1. A single uniform entry-offset constant

The 28 algebraic witnesses from the predecessor give constants `B>=0`
and `M_0` such that, for every power-of-two `M>=M_0`, every witness
with slope label `a_i` has unweighted entry valuation at most

\[
a_iM+B.                                                   \tag{3}
\]

One common `B` suffices because there are finitely many witnesses.
Their proof supplies eventual affine valuations; an explicit value
of `B` is unnecessary here.

Let `lambda_infinity,lambda_0` be the smallest weighted active Smith
orders multiplied by `ell`. Equivalently, they are the minimum
weighted entry valuations in logarithmic units. Regular homogeneous
coefficient maps preserve the physical weight labels, so (3) gives

\[
\lambda_e\le a_i\tau+f_i\cdot\theta+B\ell                 \tag{4}
\]

for every associated feature vector `f_i`, uniformly over all real
weights. Put

\[
A=\theta/\tau,\qquad
\widehat\gamma_e=(\lambda_e-B\ell)/\tau.                   \tag{5}
\]

These variables satisfy the exact 28 leading-model inequalities and
the balanced-charge equalities. No limit of `A` is being taken.

## 2. The gap duals work away from the zero-cost face

Use the physical core energies with their `O(ell)` endpoint offsets
removed:

\[
\begin{aligned}
E_X^0&=(0,-2\tau+v_l+x_l,\ 2\tau+v_t+x_t,\ v_b+x_b),\\
E_Y^0&=(0,-2\tau+y_l+u_l,\ 2\tau+y_t+u_t),
\end{aligned}                                             \tag{6}
\]

where each core weight is summed over the two places. Let `m_X^0`
and `m_Y^0` be their minima, and define

\[
T=m_X^0+m_Y^0-\lambda_\infty-\lambda_0+2B\ell.             \tag{7}
\]

The first 12 exact dual certificates imply

\[
T\ge0.                                                    \tag{8}
\]

If either physical minimum in (6) is non-low, their table gives

\[
T\ge2\tau.                                                \tag{9}
\]

The five gap certificates have an especially simple additional
property: their multiplier for the `C=0` equality is exactly `-1`.
Keeping that term instead of setting it to zero yields, whenever
both selected minima in (6) are low,

\[
E_X^0(i)-E_X^0(\mathrm{low})\ge g_i\tau-T,\qquad
E_Y^0(j)-E_Y^0(\mathrm{low})\ge g_j\tau-T,                  \tag{10}
\]

with `g=(2,4,2)` for `X` empty/top/both and `g=(2,4)` for `Y`
empty/top. This is an exact consequence of the same dual identities,
without a constraint `T=0`.

The actual low energy adds `2ell` to (6); the top energy subtracts
`2ell`; the empty and both energies have zero offset. Consequently
any actual gap above low is at least its gap in (10) minus `4ell`.

## 3. Nonpositive pressure would force quantitative freezing

Write `R=sqrt(2)` for the critical paired-type pressure multiplier.
The critical pressure decomposes as

\[
F=S+\log W_\infty+\log W_0
  -R\log Z_{X,\mathrm{core}}-R\log Z_{Y,\mathrm{core}},      \tag{11}
\]

where `S` contains the scalar Gaussian spectators and the physical
interior modes. Under (2), `S` is independent of all core gauges.
The elementary bounds

\[
\log W_e\ge-\lambda_e,\qquad
\log Z_{X,\mathrm{core}}\le-m_X/R+\log4,\qquad
\log Z_{Y,\mathrm{core}}\le-m_Y/R+\log3                    \tag{12}
\]

use the actual core minima `m_X,m_Y`.

Compare `S` with the complete Gaussian interval of length `M`.
There are 20 active Gaussian factors, ten at each place. The sum of
its negative active one-particle orders is exactly `-4M`:

```text
infinity negative orders: 1-2M, 3-2M, -3, -1
zero negative orders:     none
```

Thus its active logarithmic profile partition is at most
`4tau+20 log2`. Each physical Gaussian core includes the low state
of energy `-2tau+2ell`, so the two physical core terms contribute at
least `4tau-4ell`. The known nonnegative Gaussian interval pressure
therefore gives

\[
S\ge-20\log2-4\ell.                                      \tag{13}
\]

Combining (7), (11)--(13), and the endpoint-energy offsets yields

\[
F\ge T-(2B+8)\ell-K_0,
\qquad K_0=20\log2+R\log12.                               \tag{14}
\]

Suppose for contradiction that `F<=0`. Then

\[
T\le K_0+(2B+8)\ell.                                     \tag{15}
\]

For all sufficiently large `M`, uniformly under (1), the right side
is less than `2tau`: divide the coefficient of `ell` by `M`, and use
`tau>=c log M`. Equation (9) therefore forces both minima in (6)
to be low.

Apply (10) and include the actual endpoint offsets. Every non-low
actual energy then exceeds the low energy by at least

\[
2\tau-K_0-(2B+12)\ell.                                   \tag{16}
\]

This is the needed freezing estimate. It holds for every candidate
with nonpositive pressure, even when its gauges vary without bound.

## 4. A Gaussian lower bound valid at every profile scale

Put `L=M-2`, which is even. The
[all-low background theorem](rectangular-low-background.md) gives
an exact Gaussian interval of length `L` after conditioning on all
four low extreme states. The profile partitions dominate those of
the restricted submatrices. Using (16) only to estimate the physical
partition correction gives

\[
F\ge F_{G,L}-C_0\exp\!\left(-\frac{(2M-K_1)\ell}{R}\right),
\quad K_1=2B+12,\quad C_0=5R e^{K_0/R}.                  \tag{17}
\]

It remains to lower-bound `F_{G,L}` uniformly over the unbounded
range of scales in (1).

Let `h(x)=log(1+exp(-x))`. The even infinite-product continuation
at criticality is

\[
F_{G,\infty}(\ell)
 =4\sum_{j\ge1}h((2j-1)\ell)
  -4\sqrt2\sum_{j\ge1}h(\sqrt2(2j-1)\ell).               \tag{18}
\]

The [interval floor theorem](schur-one-species-interval.md) makes
(18) strictly positive for every `ell>0`. Its
[small-scale expansion](quantitative-interval-frontier.md) is
`ell/6+O(exp(-k/ell))`; at large scale its first term gives
`4 exp(-ell)(1+o(1))`. Positivity and continuity on compact scale
intervals therefore give a constant `c_0>0` such that

\[
F_{G,\infty}(\ell)\ge c_0\min(\ell,e^{-\ell})
\quad\text{for every }\ell>0.                            \tag{19}
\]

The positive omitted tail in the finite pressure is

\[
P=4\sqrt2\sum_{j>L/2}h(\sqrt2(2j-1)\ell),
\]

and the negative omitted tail is

\[
N=4\sum_{j>L}h((2j-1)\ell).
\]

The inequalities `x/(1+x)<=log(1+x)<=x`, followed by geometric
summation, give

\[
\frac PN\ge
\sqrt2\frac{1-e^{-2\ell}}{1-e^{-2\sqrt2\ell}}
\frac{\exp(((2-\sqrt2)L+1-\sqrt2)\ell)}
 {1+e^{-\sqrt2(L+1)\ell}}
\ge
\frac{\exp(((2-\sqrt2)L+1-\sqrt2)\ell)}
 {1+e^{-\sqrt2(L+1)\ell}}.                               \tag{20}
\]

For the second inequality, use
`1-exp(-sqrt(2)x)<=sqrt(2)(1-exp(-x))`. Under (1),
`L ell` tends to infinity uniformly. The final expression in (20)
therefore eventually exceeds one, uniformly over all permitted
`ell`. Thus

\[
F_{G,L}\ge F_{G,\infty}(\ell)
          \ge c_0\min(\ell,e^{-\ell}).                    \tag{21}
\]

For `ell<=1`, the error in (17), divided by `ell`, is at most

\[
C_0\frac M\tau
 \exp\!\left(-\left(\sqrt2-\frac{K_1}{\sqrt2M}\right)\tau\right).
\tag{22}
\]

It tends to zero uniformly under (1), because `sqrt(2)c>1`.
For `ell>=1`, divide instead by `exp(-ell)`; the ratio is at most

\[
C_0\exp\!\left(-\left(\frac{2M-K_1}{R}-1\right)\ell\right),
\tag{23}
\]

which also tends to zero uniformly. Equations (17) and (21)--(23)
contradict `F<=0`. This proves the theorem.

## 5. Exact certificate checks

The [dual-gap checker](../tools/check_rectangular_balanced_dual_gaps.py)
replays all supported rational duals and verifies the facts needed
above: the unique zero-cost minimum, the `-1` off-face multipliers,
and the 20-factor Gaussian comparison. It uses no numerical solver.

```sh
python tools/check_rectangular_balanced_dual_gaps.py \
  --output evidence/rectangular-balanced-dual-gaps.json
```

The algebraic entry bounds themselves are checked by the predecessor's
[uniform coefficient checker](../tools/check_rectangular_uniform_entry_bounds.py).
The explicit inequalities in Sections 3--4 supply the remaining
uniform analytic argument. No finite profile fit is used to extend
the result to changing gauge vectors.
