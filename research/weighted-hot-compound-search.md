# Weighted hot-code profile and search interfaces

## Complete block oracle

The [weighted integrality theorem](weighted-hot-compound-profiles.md)
allows additional full-letter diagonal weights while retaining the
occupation-block computation. The weights can be non-additive functions
of a physical state. The theorem supplies a sufficient linear cone;
it does not permit arbitrary weights without checking that cone.

[`weighted_hot_code_oracle.py`](../tools/weighted_hot_code_oracle.py)
stores the complete selected polynomial matrix of every occupation
block. Each native matrix query must resolve every endpoint-proved
channel. The implementation increases the series precision after an
unresolved query; it never returns a partial profile as complete.
The present weighted interface supports block order at most 1,024.

The method `log_profile_with_gradient(heights,scale)` takes one already
scaled height per physical letter on each of the four legs. It returns
the logarithm of the full channel sum, its height derivatives, its scale
derivative, and every scaled order. Compound polynomials use `t=w^2`,
so their native valuation scale is twice the supplied original scale.

The method `cone(scale)` returns the sufficient scaled inequalities
for a stabilized ambient realization. Every finite hot code admits
such a realization. Smaller intervals with within-leg reflection edges
need the additional same-leg constraints from the theorem.

## Mirrored gauges and optimized paired types

One may use an extra zero-place gauge `g=(gx,gy,gu,gv)` and its mirrored
infinity gauge `(gx,gu,gy,gv)`. The two compound channel sums coincide;
write their common logarithm as `log Z`. Their scalar shifts cancel,
their additional physical means double, and their natural physical
mean remains `s(EX+EY)`.

Define the physical energy vectors

```math
c_X(x)=g_x(x)+g_v(x)+\tfrac{s}{2}E(x),\qquad
c_Y(y)=g_y(y)+g_u(y)+\tfrac{s}{2}E(y).
```

After dividing both sides of the balanced pressure ratio by two, the
cost at common paired types is

```math
\rho(p_X,p_Y)
 =\frac{\log Z+\langle p_X,c_X\rangle+
                       \langle p_Y,c_Y\rangle}{H_X+H_Y}.       \tag{1}
```

At its positive interior type optimum, `rho` is the root of

```math
\log Z=\rho\left(
 \log\sum_x e^{-c_X(x)/\rho}
 +\log\sum_y e^{-c_Y(y)/\rho}\right),                         \tag{2}
```

and the two types are the corresponding normalized Gibbs weights.
[`hot_compound_pressure.py`](../tools/hot_compound_pressure.py) uses
this equation. It rejects a missing positive root or degenerate
entropy, rather than assigning a favorable numerical value.

On a smooth Smith chamber, the height derivative is the sum of the
channel-sum derivative and the four paired type marginals, divided by
`HX+HY`. The scale derivative also adds
`(mean EX+mean EY)/2` before dividing by that entropy. These are envelope
derivatives of (1); they do not differentiate the optimizing types a
second time.

## Controls and bounded search

```sh
OPENBLAS_NUM_THREADS=1 python tools/check_weighted_hot_code_oracle.py \
  --output evidence/weighted-hot-code-oracle-controls.json
```

The [checker](../tools/check_weighted_hot_code_oracle.py) binds all 1,240
channels in the independently checked literal reflection controls to
the new block interface, with zero order discrepancy. It checks 2,480
cone inequalities, four complete common-type pressure values against
the original full-tensor implementation, and 102 derivative coordinates.
The largest finite-difference derivative error was below `6.1e-10`.
The [record](../evidence/weighted-hot-code-oracle-controls.json) binds
both the input controls and every implementation source.

The worker
[`search_weighted_hot_cone.py`](../tools/search_weighted_hot_cone.py)
uses a constrained linear proposal and a complete-profile line search.
Before a profile call it rounds the scale and all heights to a common
denominator `10^10`. Any violated integer cone constraint is repaired
by adding a nonnegative occupation weight. Each Gaussian edge adds
two occupied bits, so the repair increases each constraint by twice
that integer increment. The saved minimum cone numerator is therefore
an exact integer membership check.

The resulting channel orders and Gibbs pressures are still evaluated
numerically. They require independent rational order replay and an
exact pressure comparison before a threshold claim. Completed calls,
source hashes, failures, restarts, and the declared time budget are
recorded separately. The joint gauge, scale, and type search makes no
global optimality claim.
