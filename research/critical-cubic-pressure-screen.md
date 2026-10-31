# Reproducible critical-cubic pressure screens

## Completed two-mode screen

The [critical cubic construction](full-alphabet-monomial-residuals.md) supplies
four rank-restoring pairs, three nonzero normalized coefficient products over
`F4`, and partial or complete inverse repairs. The first screen covered all
24 resulting two-mode seeds, in both complementary cut assignments.

The best floating result was `rho = 0.795392218267972`, for pair D,
coefficient product 3, the partial repair, and cut 1. This is above the
required `1/sqrt(2)`. It is a local numerical result; neither a global
optimum nor an exact pressure certificate is claimed. Several seeds have
poor local optima, so their recorded values are not lower bounds.

The [complete screen record](../evidence/critical-cubic-pressure-screen.json)
retains every constructor specification, exact seed hash, local denominator,
final weights, paired Gibbs types, both Smith lists, and restart outcomes.
The coefficient constructor remains exact. The weighted pivot choices and
pressure optimization are floating point.

Reproduce the first screen with

```sh
python3 tools/screen_full_alphabet_monomial_pressure.py \
  --modes 2 --restarts 4 --maxfun 200 --seconds 900 \
  --candidate-seconds 45 --output /tmp/cubic-full-M2.json
```

Timing and local optimization can vary between machines. The saved candidate
weights permit separate replay; a future strict improvement must still be
checked with rational weights and exact pressure arithmetic.

## Search interface

The [driver](../tools/screen_full_alphabet_monomial_pressure.py) accepts
nonconsecutive modes, either complementary cut assignment, coefficient
products, both repairs, and independent Laurent powers for the two cubics.
It restores the canceled local denominator exactly once, in the profile
oracle, and checks that its direct pressure evaluation agrees with its
paired Gibbs root. Native unresolved profiles are counted as failures.
They are never accepted as low-cost candidates.

Three- and four-mode direct matrices grow rapidly. The bounded active-space
coordinate constructors are preferable for large intervals; their retained
coefficient precision must be respected when reusing exported arrays.

## A complementary polynomial perturbation interface

Starting from an already polynomial exterior tensor `C(w)`, the constructor
[ramified_polynomial_crossings.py](../tools/ramified_polynomial_crossings.py)
forms

\[
C(v^d)\left(1+\sum_j a_jv^{k_j}(1+v)^{h_j}H_j\right)
\]

over characteristic two, with positive integer `d,h_j`, nonnegative `k_j`,
and arbitrary nonempty exterior monomials `H_j`. At `v=1`, every correction
vanishes and `C(v^d)=C(1)`. Thus the unweighted polynomial endpoint is
unchanged entry by entry, and any typed crossing endpoint is inherited.
This operation is not justified by applying it to the numerator of a
rational tensor whose poles have not already been canceled.

The [screening driver](../tools/screen_ramified_polynomial_crossings.py)
starts with the polynomial two-mode frozen Gaussian, scans balanced even
monomials or opposite-charge cubic pairs at explicit Gaussian gauges, and
optimizes the best candidates over all physical-letter weights. It records
ramification, coefficient hashes, and source hashes. These are additional
search families, without a score claim. The full fixed-seed degree is kept;
ramification is not declared free input-dependent overhead.
