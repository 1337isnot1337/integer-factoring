# An all-scale floor for the homogeneous middle-fiber family

Every one of the 120 homogeneous restrictions in
[the middle-fiber experiment](reflection-homogeneous-restrictions.md)
has common-type Smith-pressure rate at least `5/7`. In particular, none
of these restrictions reaches `1/sqrt(2)` at any scale. This statement
uses the exact saved profiles of that complete finite family. It does
not exclude other homogeneous subspaces or different physical gauges.

The positive comparison in the earlier note remains valid: one
nonorthogonal map improves a coordinate deletion at the same scale.
The present calculation shows that this improvement is insufficient
for the desired asymptotic threshold.

## 1. The polynomial to certify

Let `a_i` and `b_i` be the restored Smith orders at the two endpoints.
Let `H_j` be the fifteen retained physical heights. Write

$$Z(r)=\sum_i r^{a_i},\qquad
  I(r)=\sum_i r^{b_i},\qquad
  A(x)=\sum_j x^{H_j}.$$

The previous construction proves that the pressure is
`B(s)=Z(1/s)I(1/s)` and the total type energies are `4H_j`.
Consequently, a trial rate `rho` uses the comparison

$$B(s)\le
\left(\sum_j s^{-H_j/\rho}\right)^{4\rho}. \tag{1}$$

At `rho=5/7`, put `x=s^(-1/5)`. The opposite inequality is equivalent to

$$P(x)=\bigl[Z(x^5)I(x^5)\bigr]^7-A(x^7)^{20}\ge0. \tag{2}$$

The two endpoint polynomials must both be retained in (2). Some of the
nonorthogonal maps have different profiles at the two endpoints.

For each fixed finite scale, the right side of (1) increases with
`rho`: the derivative of its logarithm is four times the entropy of
the corresponding positive distribution. Thus strict positivity in
(2) excludes every trial rate at most `5/7`. Taking an infimum over
scales still gives a rate floor of `5/7`.

## 2. A finite coefficient certificate

For an integer polynomial `P(x)=sum_(j=0)^D p_j x^j`, form the two
cumulative sums

$$a_m=\sum_{j=0}^m p_j,\qquad
  b_m=\sum_{i=0}^m a_i\qquad(0\le m\le D). \tag{3}$$

The numbers `b_m` are the coefficients of `P(x)/(1-x)^2` through
degree `D`. Beyond that degree they are exactly

$$b_{D+k}=b_D+k a_D\qquad(k\ge0). \tag{4}$$

Therefore the finite conditions

$$b_m\ge0\quad(0\le m\le D),\qquad a_D>0 \tag{5}$$

certify `P(x)>0` for every `0<x<1`. Equation (4) supplies the entire
infinite tail. At `x=1`, positivity follows separately from
`P(1)=a_D>0`.

Every one of the 24 distinct profile pairs satisfies (5). The largest
degree in (2) is 3150. Each endpoint has 225 Smith orders and the
physical alphabet has fifteen letters, so the final first sum is
always

$$a_D=225^{14}-15^{20}=15^{28}-15^{20}>0. \tag{6}$$

These are exact integer inequalities. No numerical optimization,
floating-point root bound, or sample of scales enters the certificate.

## 3. Replay and provenance

The checker
[`check_homogeneous_restriction_floor.py`](../tools/check_homogeneous_restriction_floor.py)
reads the previously constructed complete profile evidence. It checks
that all 120 projector pairs occur, reconstructs their embeddings,
retractions and physical heights, and verifies the source engine and
checker hashes. It then constructs (2) using ordinary integer
polynomial convolution and checks all coefficients in (5).

The result is saved in
[`homogeneous-restriction-pressure-floor.json`](../evidence/homogeneous-restriction-pressure-floor.json).
It contains a certificate summary for each distinct pair, the complete
source-evidence hash, and both source-code hashes. This replay does not
repeat the earlier four-minute Smith-profile construction; that
construction and its independent Toeplitz controls remain the source
of the input profiles.

```sh
python3 tools/check_homogeneous_restriction_floor.py \
  --output evidence/homogeneous-restriction-pressure-floor.json
```

The coefficient replay takes less than a second in the recorded
environment. The floor is a statement about this finite family and its
inherited type-pressure calculation, not a lower bound on every
factoring algorithm or every possible coefficient-rank method.
