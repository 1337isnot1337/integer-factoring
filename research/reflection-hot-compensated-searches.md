# Search checkpoint for cold-compensated reflection hot cores

## Result and scope

The completed searches below use the exact
[compensated negative-height construction](reflection-hot-compensated-gauges.md).
They retain the full Gaussian spectators, the added cold-bit type
energies, the original polynomial degree, and both complete active cuts.
None improved the Gaussian pressure in its recorded search scope.

Exact field coefficients and complete Smith profiles are distinct from
the floating scale and height optimizations. These records do not prove
an all-gauge or all-scale obstruction. The separate
[physical-leg bulk cone](reflection-hot-leg-bulk-gauges.md) has its own
search records and is not part of the counts below.

## 1. Continuous compensated-height searches

The driver is
[`search_reflection_hot_compensated_gauges.py`](../tools/search_reflection_hot_compensated_gauges.py).
For one selected hot bit it uses nonnegative bulk heights and eight
nonnegative hot-state slacks, with
\(h_a(1)=\mathrm{slack}_a-2b_a\). Thus negative total hot heights
are permitted while every removed hot-cold Gaussian factor remains
regular. Bulk heights lie in \([0,8]\), slacks in \([0,24]\), and
the logarithmic scale lies between one quarter and three times the
natural Gaussian optimum. Each coefficient row uses four starts and
a nominal 100-evaluation limit per start; line-search bookkeeping can
slightly exceed that limit, and every actual count is recorded.

| Interval | Complete coefficient scope | Bulk parameters | Best pressure | Best non-Gaussian pressure |
|---:|---|---|---:|---:|
| 16 | All 234 saved natural-profile representatives | One bulk per cut | 0.7106044045201707 | 0.7580701412490927 |
| 32 | Nine fixed rows and sixteen seeded random rows | One bulk per boundary per cut | 0.7083157559380195 | 0.7618422781092767 |
| 64 | The same 25 rows | One bulk per boundary per cut | 0.7075072113149573 | 0.7455948452332484 |

The Gaussian attained the best value in each row of the table. The
234 representatives exhaust the natural-profile classes, not the
coefficient tensors at arbitrary gauges: two tensors that agree at
the natural gauge can separate after heights change. The other two
runs are explicitly finite sampled coefficient sets.

Saved complete records:

- [M16, 234 representatives](../evidence/reflection-hot-compensated-search-M16.json).
- [M32, 25 rows](../evidence/reflection-hot-compensated-search-M32.json).
- [M64, 25 rows](../evidence/reflection-hot-compensated-search-M64.json).

## 2. A complete coefficient grid on seven fixed negative-height rays

The driver
[`search_reflection_hot_compensated_rays.py`](../tools/search_reflection_hot_compensated_rays.py)
reuses cached exact coefficient components and enumerates all
\(4^8=65,536\) coefficient rows in the established one-hot-bit
\(\mathbb F_4\) family at \(M=64\).

Each ray specifies a two-by-two bulk array: its rows are the zero and
infinity cuts, and its columns are the left and right boundary bulks.
Every hot relative height is exactly minus twice its boundary bulk.
The seven arrays are

\[
\begin{pmatrix}1&0\\0&1\end{pmatrix},\quad
\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
\begin{pmatrix}1&1\\0&0\end{pmatrix},\quad
\begin{pmatrix}0&0\\1&1\end{pmatrix},\quad
\begin{pmatrix}1&1\\1&1\end{pmatrix},\quad
\begin{pmatrix}2&0\\0&2\end{pmatrix},\quad
\begin{pmatrix}2&2\\2&2\end{pmatrix}.
\tag{1}
\]

This is 917,504 complete cut profiles: two cuts for every coefficient
row and every ray. The profile-pair class counts on the seven rays
are respectively

\[
541,\ 541,\ 496,\ 471,\ 986,\ 557,\ 971.
\tag{2}
\]

Each complete integer profile was evaluated at three logarithmic scales,
equal to one half, one, and two times the natural Gaussian optimum:

\[
0.02378746726020672,\quad
0.04757493452041344,\quad
0.09514986904082688.
\tag{3}
\]

The minimum relative logarithmic cost against the Gaussian at the
same heights was zero; no sampled value was negative. The minimum
sampled baseline-equation excess was
\(0.15027753642807617\), and the minimum critical-equation excess
was \(0.1597889057715136\). These are floating evaluations of exact
profiles, not interval-certified all-scale inequalities.

Forty-two independent component-compiler controls compare the cached
linear construction with the full compensated model. The
[complete record](../evidence/reflection-hot-compensated-rays-M64.json)
contains all scope parameters, the class counts, source hashes, and
a digest of the complete profile stream.

## 3. A nonlinear two-bit OR-state height family

For two selected hot bits, the driver
[`search_reflection_hot_or_heights.py`](../tools/search_reflection_hot_or_heights.py)
uses the same height on each physical leg at a given cut:

\[
h_\nu(s)=u_\nu\,\mathbf1_{s\ne0}-2b_\nu|s|,
\qquad u_\nu,b_\nu\ge0.
\tag{4}
\]

The OR-state slack in (4) is nonlinear in the two occupations. It is
monotone, so the compensated regularity condition holds. The finite
grid takes \(b_0,b_\infty\in\{0,1,2,4\}\) and
\((u_0,u_\infty)=r(b_0,b_\infty)\) for
\(r\in\{0,1,2\}\). Removing duplicates leaves 46 height pairs.

At \(M=64\), all 46 pairs were evaluated for the Gaussian and
three explicit exterior defects from the saved two-hot-bit catalogue.
This gives 184 complete two-cut height evaluations. Each profile was
calculated at an exactly representable scale and converted back to
integer orders; subsequent scale optimization reused those complete
orders. The scale search range was one quarter through four times the
natural Gaussian optimum.

Every tensor's best recorded height pair had zero bulk and zero slack.
The best pressures were

| Tensor | Best pressure in the recorded scale interval |
|---|---:|
| Gaussian | 0.7075072113149573 |
| First defect | 0.7201975818156903 |
| Second defect | 0.7213132255225648 |
| Third defect | 0.7204731192972906 |

The three defect scale optima were at the upper edge of the permitted
scale interval, so these values are not claims of their unrestricted
scale minima. The complete factors, every height-grid value, and both
best profiles are retained in the
[saved record](../evidence/reflection-hot-or-heights-M64.json).
An independent original \(M=8\) nonlinear control verifies both
weighted cuts and all 65,536 endpoint coefficients before this search.

## Interpretation

These searches broaden the earlier positive-height experiments to
negative hot heights with all compensating cold costs included. Their
negative results are scoped to the listed coefficient sets, height
ranges, and scales. They leave other exterior polynomials, other state
height functions, the separate-leg bulk cone, and gauges outside these
regularity cones open.
