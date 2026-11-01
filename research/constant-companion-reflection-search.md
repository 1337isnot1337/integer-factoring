# Constant companions of the complete dense rank-one catalogue

## Result and scope

The [complete four-letter rank-one catalogue](dense-rank-one-coupled-frames.md)
contains 7,176 polynomial coupled frames with positive pole orders.
Adding each of the three nonzero `F4` constants to the Laurent parameter
gives 21,528 further explicit polynomial tensors. All their coefficients,
strict endpoints, exact degree 24, and distinct tensor hashes have been
reconstructed.

The completed numerical search found no better crossing rate. Its best
balanced value was approximately `0.731104093726`, compared with the
original four-letter value `0.731080935974`. Its best asymmetric support
value was approximately `1.001436002791`, above the required threshold
one. These are reproducible search observations, not global optimality
certificates.

## 1. The companion identity

Let `G^2=0` in Bell coordinates and let `f=w^-k`. The coupled physical
frame uses `I+fG` on the first and third legs and its inverse transpose
on the second and fourth legs. For any field constant `c`,

```math
I+(f+c)G=(I+cG)(I+fG),\qquad (I+cG)^{-1}=I+cG.
```

Thus the companion is a constant physical change of basis applied to
the accepted Laurent frame. It is polynomial whenever the original
frame is polynomial. A constant invertible map preserves exact degree,
and the paired dual action preserves the strict crossing endpoint.

The leading coefficient may change under this constant frame. Complete
weighted local profiles must therefore be recomputed. Endpoint
preservation does not imply preservation of a fixed physical type or
its weighted pressure.

If the complete four-leg expansion is

```math
T(f)=\sum_{r=0}^{4} f^r T_r,
```

then the translated expansion is computed exactly by

```math
T(f+c)=
\sum_{j=0}^{4}f^j
\sum_{r=j}^{4}\binom rj c^{r-j}T_r.              \tag{1}
```

The implementation uses the binary parity of each binomial
coefficient. It retains all four-leg interactions and does not require
intermediate single-leg actions to be polynomial.

## 2. The completed search

The worker started at `2026-10-31 22:43:24.872813 UTC` and finished at
`22:50:25.515404 UTC`, after `420.640597907` measured seconds. It stopped
because its declared finite screen and refinement limits were reached.
Its 1,800-second time limit was not exhausted.

| Quantity | Completed value |
|---|---:|
| Exact companion tensors | 21,528 |
| Distinct tensor hashes | 21,528 |
| Distinct natural two-cut profile vectors | 155 |
| Full gauge and common-type refinements | 80 |
| Optimizer loops completed | 400 |
| Loops reporting optimizer convergence | 250 |
| Resolved objective calls | 171,763 |
| Unresolved objective calls | 3,264 |

Every refinement completed three balanced loops and two asymmetric
loops. Calls for which the local oracle was unresolved received no
pressure claim. Optimizer convergence flags are retained separately
from the best resolved evaluations.

The natural-profile screen optimized a common scalar and Gibbs types
for each complete profile. It obtained a best value
`0.8391297016169927`. Full physical gauge and common-type optimization
then found the values reported above. Both winners use catalogue
generator 4035, pole four, and companion constant one.

The best balanced winner replayed as `0.7311040937261724`, differing
from its recorded value by about `1.33e-13`. The asymmetric winner
replayed exactly at the stored floating precision. Neither replay is
an exact inequality certificate for the continuous optimum.

## 3. Exact binding and independent controls

The [archive](../evidence/dense-coupled-companion-pressure.json) records
the source hashes, generator indices, poles, constants, complete
natural profile vectors, numerical parameters, optimizer outcomes,
and actual start and finish times.

The [checker](../tools/check_dense_coupled_companions.py) binds every
record to the published complete rank-one catalogue, verifies every
recorded search-source hash, and reconstructs all 21,528 complete
tensors. It also binds all 80 refinement records to their exact input
tensors.

An independent exponent-dictionary implementation directly applies
the four matrices with entries `I+(w^-k+c)G`, using separate finite-field
arithmetic. It checks 5,376 complete tensor entries across seven
generators and all three constants, including an order-eight frame.
This checks the binomial translation in (1) against the actual coupled
action. Fresh complete-profile winner evaluations supply the numerical
replay.

The [replay evidence](../evidence/dense-coupled-companion-replay.json)
records all these checks. The source catalogue digest is preserved
exactly; this search used the published all-pole catalogue directly.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
python3 tools/check_dense_coupled_companions.py \
  --output /tmp/dense-coupled-companion-replay.json
```

The construction enlarges the tested physical-basis family. This finite
search does not establish a lower bound for arbitrary constant bases,
other Laurent coefficient functions, or affine combinations of
different generators.
