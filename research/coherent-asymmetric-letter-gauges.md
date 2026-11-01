# Independent letter gauges at the two crossing cuts

The [shared-gauge construction](coherent-shared-letter-gauges.md) extends
to two independent letter-potential vectors at a common profile scale.
Both physical costs must be charged. The implementation and exact
rational-profile controls are retained; the bounded 32-letter numerical
search did not find an improvement.

## 1. A sufficient two-place cone

Keep the coherent physical maps, heights `h_j`, and induced deletion
edges of the shared construction. Fix `u>0`. Choose two potential vectors
`f^(0), f^(infinity)`, each satisfying

\[
 (f_j-f_k)+(f_l-f_m)+2u(a\mathbin{\&}b)\ge0
\]

for every pair of induced edges `(a,j,k),(b,l,m)`. Normalize each vacuum
potential to zero. The two vectors need not agree.

Each vector belongs to the previously proved shared regularity cone.
It therefore gives regular, invertible descended Gaussian factors at its
assigned cut. This is a sufficient domain for independent gauges; no
claim is made that independently optimizing the two larger cut-specific
regularity regions would give the same domain.

For every nonvacuum label, the deletion-path argument gives
`u h_j+2 f_j^(c)>0` at each place. Averaging these inequalities shows
that the physical energy used below is positive.

## 2. Physical cost and derivatives

Let `lambda_(c,i)` be the complete weighted Smith orders using the
potential vector for cut `c`. Each cut has exactly `q^2` channels.
Define

\[
 e_j=u h_j+f_j^{(0)}+f_j^{(\infty)}.
\]

The two cut gauges add the total cost
`4 sum_j p_j(f_j^(0)+f_j^(infinity))`, alongside the natural height cost
`4u sum_j p_j h_j`. The balanced pressure equation is consequently

\[
 \log\sum_i e^{-\lambda_{0,i}}+
 \log\sum_i e^{-\lambda_{\infty,i}}
 =4\rho\log\sum_j e^{-e_j/\rho}.
\]

Putting both potential vectors equal recovers the shared-gauge equation
with `e_j=u h_j+2f_j`. Holding the physical energies fixed while changing
one cut's gauge would omit its cost and is not the implemented objective.

Write `L_c=log sum_i exp(-lambda_(c,i))` and let
`p_j=exp(-e_j/rho)/sum_k exp(-e_k/rho)`. In one fixed pivot region,
implicit differentiation gives

\[
 \frac{\partial\rho}{\partial f_j^{(c)}}=
 \frac{\partial L_c/\partial f_j^{(c)}+4p_j}{4H(p)},
 \qquad
 \frac{\partial\rho}{\partial u}=
 \frac{\partial L_0/\partial u+\partial L_\infty/\partial u
       +4\sum_jp_jh_j}{4H(p)}.
\]

At a pivot boundary these are derivatives of the selected region, so a
local optimizer has no global guarantee. At equal potentials the two
potential derivatives add to the derivative of the shared objective.

The [oracle](../tools/coherent_asymmetric_gauges.py) reuses the complete
compound blocks and precision checks. It changes only the two gauge
vectors, their combined physical cost, and the corresponding derivatives.

## 3. Independent controls

The [checker](../tools/check_coherent_asymmetric_gauges.py) uses a
six-letter inverse system, a nonorthogonal eight-letter map, and a
nonlinear eight-letter quotient. It chooses two distinct rational
potential vectors at `u=1/2`, including negative potentials.

For each place separately, the checker reconstructs literal ambient
Pfaffians, ramifies by 32, and applies the integer letter shifts.
Independent binary-polynomial Smith elimination agrees with all six
complete weighted profiles, containing 328 orders in total. Three
equal-potential comparisons recover the shared pressure, physical type,
and combined derivative. Eighteen central finite differences at interior
points check the independent-potential gradients to error below `10^-6`.
The literal profiles are exact; the derivative checks are numerical.

The [control receipt](../evidence/coherent-asymmetric-gauge-controls.json)
retains both potential vectors, all integral profiles, and source hashes.

## 4. Bounded 32-letter search

The [search driver](../tools/search_coherent_asymmetric_gauges.py) was run
on the published binary boundary map with kernel 7 and covector 15.
Its 1,200-second budget permitted 317 objective evaluations across three
starts. Two optimizer attempts completed and the third reached the
time budget. The best retained numerical value was the initial natural
gauge, approximately `0.7121150127573352`.

The [complete search record](../evidence/coherent-asymmetric-q32.json)
contains the maps, parameters, full profiles, optimizer outcomes, and
timing. This finite failure is not a global optimality result. No floating
value in that record is presented as an exact pressure certificate, and
the complete factoring exponent remains unchanged.

```sh
python3 tools/check_coherent_asymmetric_gauges.py \
  --output /tmp/coherent-asymmetric-gauge-controls.json
```
