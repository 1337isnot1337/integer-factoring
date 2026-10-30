# Computed Euler and charge-filtered Ward bounds

The [Euler identity](noncommutative-crossing-ward-euler.md) and
[one-sided charge gauge](noncommutative-crossing-ward-charge-filtration.md)
give smaller, explicit factorizations for the coefficients of radical-code
crossings. The computations below use those factorizations with the
[stable rectangular evaluator](ward-code-rank-bounds.md).

The new bounds improve the recorded coefficient-cost estimates. They do
**not** prove a factoring exponent above `1/2`. Failure of an upper bound
to pass the criterion is not a lower bound on the actual coefficient ranks.

## Cost convention

For a crossing with `d` logical labels, write

\[
 S=\sum_{a,b}\sqrt{r_{a,b}c_{a,b}},\qquad
 \rho=\frac{\log_2 S}{2\log_2 d}.
\]

The seed criterion is `S < d^sqrt(2)`, or equivalently
`rho < 1/sqrt(2) = 0.70710678...`. The programs compute integer upper
bounds by summing the ceilings of the square roots of certified rank
products. They take the minimum of each new cut bound and the previous
graded Ward bound. Every displayed `rho` is the logarithm of that integer
upper bound, not an exact rank-cost measurement.

The one-piece crossing itself always gives the independent upper bound
`d^2`. A coefficient-cost bound greater than `d^2` can therefore be
discarded when choosing a decomposition.

## One global form

For a global symplectic form, the incidence operator becomes a parameter
derivative. Homogeneity turns each derivative coefficient into a scalar
parity. Both cuts then factor through only two grade rectangles, without
a state-channel multiplier.

The labels below are the full homogeneous radical codes; `P` is the
number of physical mode pairs and `L` is the exterior grade.

| `P` | `L` | Labels `d` | Previous coefficient-cost bound | Euler coefficient-cost bound | New `rho` |
|---:|---:|---:|---:|---:|---:|
| 16 | 10 | 8,811,817 | 9,303,652,473,389 | 2,683,995,769,279 | 0.89479227 |
| 24 | 14 | 57,671,029,749 | 9,582,897,190,018,480,730 | 2,594,735,082,693,044,913 | 0.85559731 |
| 32 | 18 | 402,083,973,823,997 | 9,051,057,704,148,292,919,789,805 | 2,392,370,413,548,346,825,207,815 | 0.83464461 |

The largest explicit matrices used in these three calculations have
dimensions 32, 48, and 80. Integer direct-sum multiplicities account for
the much larger logical spaces. In the last case the old bound had
`rho = 0.85442886`; the new bound is still above the required threshold.

The complete profiles, integer costs, and source hashes are in
[ward-differential-code-bounds.jsonl](../evidence/ward-differential-code-bounds.jsonl).
The small `P=4,L=4,d=26` control is included there too. Its coefficient-cost
bound drops from 3,399 to 1,696, while the independent one-piece bound
is already 676.

Reproduce the four recorded cases with:

```sh
python tools/ward_differential_bounds.py \
  --output /tmp/ward-differential-code-bounds.jsonl
```

## Several forms without charge-channel multiplicities

For `k` forms, the dual-incidence identity moves the complete Ward
operator to one side of either cut. It works with arbitrary form weights
in characteristic two. The resulting estimate uses at most
`binom(k+2,2)` grade rectangles, with no factor for the number of charge
values, selected creation operators, or physical sites. For binary forms,
the independent charge-gauge proof shows that only the maximal permitted
grade shift survives; each constant interface entry is one `k` by `k`
minor.

The recorded examples use `k` disjoint indicator forms, each on `b`
physical pairs, with `P=b*k`. Their output maps preserve the admissible
syndromes. Since each block size is even, the total form Gram matrix
vanishes, and the joint vanishing order is at least `2*k`.

| `b` | `k` | `L` | Labels `d` | Previous generic coefficient-cost bound | New coefficient-cost bound |
|---:|---:|---:|---:|---:|---:|
| 4 | 2 | 8 | 678 | 7,254,199 | 5,230,575 |
| 4 | 3 | 10 | 2,031 | 113,969,022 | 85,783,992 |
| 8 | 2 | 12 | 3,623,919 | 76,519,828,278,838 | 28,948,403,307,646 |
| 8 | 3 | 18 | 9,409,503,908 | 918,232,013,534,284,889,596 | 212,607,881,453,648,885,777 |

Here “previous generic” means the earlier Ward bound without a
disjoint-block-specific improvement to its channel count. All four new
coefficient bounds are still greater than `d^2`. Thus these examples
verify the general improvement to the factorization, but do not improve
the best available decomposition cost for their codes.

The complete profiles and hashes are in
[ward-channel-free-full-bounds.jsonl](../evidence/ward-channel-free-full-bounds.jsonl).
Reproduce them with:

```sh
python tools/ward_channel_free_bounds.py \
  --output /tmp/ward-channel-free-full-bounds.jsonl
```

## Independent exact controls

[The determinant checker](../tools/check_ward_charge_filtration.py) does
not infer the new theorem from the large numerical bounds. It checks
the algebraic components and compares the resulting estimates with
literal small tensors:

- The symmetric determinant translation identity is checked as a formal
  polynomial through matrix size six. Repeated variables remain genuine
  polynomial powers; the computation does not pass to a Boolean quotient.
- 800 deterministic exact cases check the rank-one finite-difference
  identity against its bordered determinant, including repeated columns.
- 830 cases check the operator-coefficient augmentation pairing bounds.
- Every centered jet of the full and orientation-half codes at
  `(P,L)=(3,3)` and `(4,4)` is compared with both new bounds.
- All recorded centered jets of the overlapping two-form constructions
  on 6 and 34 labels satisfy the new bound and the applicable joint
  vanishing-order restriction.

The separate [Euler checker](../tools/check_noncommutative_crossing_ward_euler.py)
checks complete polynomial tensor identities, arbitrary grade-preserving
output maps in small controls, and 1,096 typed Gram/minor identities.
The [local gauge checker](../tools/check_pfaffian_charge_gauges.py) checks
all 256 entries of each of the one-sided and two-sided polynomial
identities, together with 1,201 augmentation pairing ranks.
The [dual-incidence checker](../tools/check_ward_dual_incidence.py) adds
512 first-order identities, 256 incidence/derivative identities, 2,560
multisite product controls, 1,320 maximal-shift finite differences, and
2,240 combined-minor identities.

Run the inexpensive proof controls with:

```sh
python tools/check_noncommutative_crossing_ward_euler.py
python tools/check_pfaffian_charge_gauges.py
python tools/check_ward_charge_filtration.py
python tools/check_ward_dual_incidence.py
```

The corresponding records are
[Euler evidence](../evidence/noncommutative-crossing-ward-euler.json),
[gauge evidence](../evidence/pfaffian-charge-gauges.json), and
[charge-filtration evidence](../evidence/ward-charge-filtration.json).
The sharper descent is also checked in
[dual-incidence evidence](../evidence/ward-dual-incidence.json).

The remaining issue is the rank after imposing the actual label maps.
The full-rectangle bounds retain directions that a particular code can
kill. The new identities remove a channel-count loss, but the recorded
calculations do not yet control those surviving correlations sharply
enough to meet the seed criterion.
