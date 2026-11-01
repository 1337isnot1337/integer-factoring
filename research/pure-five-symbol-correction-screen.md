# Certified finite screen of pure five-symbol corrections

## Result and scope

The completed search tests 192,238 endpoint-zero pure-tensor
correction queries around the retained degree-18 five-symbol crossing.
Every tested correction has a complete profile product at least as
large as the reference product at every positive scale, at both of the
two displayed integer gauges. This comparison is exact. It follows
from certified complete Smith orders and integer prefix inequalities,
without a numerical scale search for the candidates.

The finite result covers all 167,238 Frobenius-orbit representatives
in a specified two-dimensional physical subspace and another 25,000
seeded queries over the full physical alphabet. It is not an exclusion
of arbitrary pure corrections, other scalar polynomials, other
physical gauges, or sums of pure terms. It produces no new factoring
exponent.

The [historical search record](../evidence/pure-five-symbol-correction-search.json)
used the original structured accelerator. Its scope explicitly leaves
unreplayed profiles uncertified. The separate
[complete stream certificate](../evidence/pure-five-symbol-correction-stream-certificate.json)
recomputes every profile using the tracked-error engine. Thus the
certificate strengthens the evidence without changing the historical
record or its source hashes.

## 1. Literal endpoint-preserving family

Work over `F4 = F2[omega]/(omega^2+omega+1)`. Let `T0(z)` be the exact
625-entry reference from the
[recentered crossing construction](f4-recentered-crossing-lattices.md),
with degree at most 18. Its label order has natural energies

\[
(0,10,6,2,8).
\]

For four nonzero constant vectors `a,b,c,d` in `F4^5`, a nonzero
scalar `lambda`, and a polynomial `f` of degree at most 18 with
`f(1)=0`, define

\[
T_{xyuv}(z)=T_{0,xyuv}(z)+\lambda f(z)a_xb_yc_ud_v.
\]

Every tensor entry is retained literally, and

\[
T_{xyuv}(1)=\delta_{xv}\delta_{yu}.
\]

At the infinity cut `(X,U)|(Y,V)`, the correction has row factor
`a tensor c` and column factor `b tensor d`. At the zero cut
`(X,Y)|(U,V)`, the factors are `a tensor b` and `c tensor d`.
It is rank one at both cuts over the rational-function field.

The two tested gauges are the natural energies above and the two-zero
energies `(0,12,6,0,8)`. At infinity the row and column heights are
`E(u)/2` and `E(v)/2`; at zero the row height is `(E(x)+E(y))/2`
and the column height is zero. The full physical energy is charged in
the paired Gibbs equation. The two-zero gauge has zero energy on
physical labels zero and three.

## 2. Exact fixed-frame reduction

The [frame adapter](../tools/pure_five_symbol_frames.py) constructs the
complete local reference matrices. At infinity it replaces `T(z)` by
`z^18 T(1/z)` and restores the scalar shift of 18 in the reported
orders. No infinity shift is discarded.

For each place and gauge, exact polynomial row and column operations
give

\[
L A_0 R=D\pmod {z^{128}},
\]

where `D` is diagonal and `L(0),R(0)` are invertible. The certificate
checks all 2,500 entries of these four identities and all eight
constant-frame ranks. A pure correction becomes

\[
D+(Lu)(v^T R)\pmod {z^{128}}.
\]

The scalar polynomial, physical monomial heights, and global field
coefficient are included in these vectors. F4 scalar polynomials are
represented by two binary coefficient planes; the validator requires
both planes to have even endpoint parity.

The [certified structured engine](certified-structured-smith.md)
tracks absolute precision through all divisions and Schur updates.
Every accepted pivot must be below the tracked error bound and must
be a certified global minimum in the residual matrix. The stream
replay used work precision 512 throughout; no higher-precision retry
was needed. Its minimum pivot error margin was 488.

The largest accepted local Smith order was 55, leaving a margin of
73 below the frame precision 128. Two matrices congruent modulo
`z^128`, whose complete exposed invariant orders are all below 128,
have the same invariant orders. One proof uses the equal ranks of
their maps on every finite quotient up to that precision. This
transfers the certified structured profiles to the complete original
polynomial matrices. It does not assume that the truncated frames
are identities as untruncated polynomial matrices.

## 3. Integer prefix comparison at every scale

Let the complete local orders be `mu` and `nu`. Sort their pairwise
sums into an increasing vector

\[
A=\operatorname{sort}\{\mu_i+\nu_j:1\leq i,j\leq25\}.
\]

Let `B` be the analogous 625-entry vector for the reference. Suppose

\[
\sum_{i=1}^k A_i\leq\sum_{i=1}^k B_i
\qquad(1\leq k\leq625). \tag{1}
\]

Then the candidate's profile product is no smaller at every positive
scale. The direction of the comparison matters: smaller low-order
prefixes increase a decreasing convex exponential sum.

To prove this directly, set

\[
H_A(t)=\sum_i(t-A_i)_+
=\max_{0\leq k\leq625}
\left(kt-\sum_{i=1}^k A_i\right).
\]

Equation (1) gives `H_A(t) >= H_B(t)` for every real `t`. For any
positive `s`, the exact integral identity

\[
e^{-sa}=s^2\int_a^\infty(t-a)e^{-st}\,dt
\]

therefore implies

\[
\sum_i e^{-sA_i}\geq\sum_i e^{-sB_i}.
\]

The sums factor as

\[
\left(\sum_i e^{-s\mu_i}\right)
\left(\sum_j e^{-s\nu_j}\right).
\]

This argument allows negative infinity orders and does not require
equal total sums. The
[implementation](../tools/convolution_prefix_dominance.py) checks (1)
using integer arithmetic. Physical energies are unchanged by each
correction. Increasing the profile-product logarithm can only increase
the paired Gibbs root, so the comparison also excludes a better ratio
at the same scale and energies.

This is a comparison with the reference at every scale. The decimal
reference minima stored in the search record remain numerical values;
the prefix certificate does not by itself certify their global
minimum to a stated number of decimal places.

## 4. Completed finite schedule

The structured part uses every projective F4 vector in the span of
physical labels zero and three:

\[
e_0,\quad e_3,\quad e_0+e_3,\quad
e_0+\omega e_3,\quad e_0+\omega^2e_3.
\]

There are `5^4 * 3 = 1875` normalized pure tensors with a nonzero
global field coefficient. Frobenius fixes `3^4 = 81` of them, so there
are `(1875+81)/2 = 978` conjugacy orbits. Frobenius preserves all
Smith valuations and fixes the binary reference. The driver tests
every orbit with every binary endpoint-zero binomial

\[
f(z)=z^a+z^b,\qquad0\leq a<b\leq18.
\]

This gives `978 * binomial(19,2) = 167238` queries. It is exhaustive
for this stated finite family, including all its field coefficients.

The remaining 25,000 queries use random projective vectors over the
whole five-dimensional alphabet. The seed is `2026103217`. Each leg
has a sampled nonempty support and first nonzero coefficient one.
Three quarters of scalar draws are binary binomials. The remaining
draws use genuine F4 trinomials

\[
f(z)=(z^a+z^c)+\omega(z^b+z^c),\qquad a<b<c.
\]

Both coefficient planes vanish at one. The random schedule is a
finite list of queries; no claim of distinctness or exhaustive
coverage is made for it.

| Gauge | Certified queries | Distinct profile pairs | Queries satisfying (1) |
|---|---:|---:|---:|
| Natural | 192,238 | 21,493 | 192,238 |
| Two-zero | 192,238 | 26,042 | 192,238 |

Thus every tested candidate is excluded by the exact relative
all-scale comparison. There were no candidates requiring a numerical
pressure search after this test.

## 5. Complete historical replay

The [search driver](../tools/search_pure_five_symbol_corrections.py)
completed in about 538 seconds. During that run, it also made 820
comparisons with full original matrix profiles and checked 128,125
literal endpoint entries. Those controls were useful for the search;
they were not used as a substitute for full stream certification.

The [independent stream checker](../tools/check_pure_five_symbol_stream.py)
regenerates the parameter schedule from its frozen source and seed.
It recomputes all 768,952 complete matrix profiles through the
tracked-error engine and the checked frame transfer. It then rebuilds
the rolling parameter/profile digest. The digest exactly matches the
historical search:

```text
bc3649c16e5c1861fa869569a966c03b518bcd164c8e7a1ac336c1f2d2713688
```

It also certifies 96 additional profiles of 24 literal tensors using
direct complete polynomial matrices and independently constructed
regular frames. Those tensors check 15,000 literal endpoint entries.
The full replay completed in about 401 seconds.

Historical source hashes resolve either to the matching worktree or
to a matching read-only Git blob. The checker does not rewrite old
hashes when a driver changes. The certificate separately records its
own certified-engine source hashes.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/search_pure_five_symbol_corrections.py --random-queries 25000 --seconds 1800 --seed 2026103217 --output /tmp/pure-five-symbol-correction-search.json
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/check_pure_five_symbol_stream.py evidence/pure-five-symbol-correction-search.json --output /tmp/pure-five-symbol-correction-stream-certificate.json
```

Sums of several pure terms and independent physical gauges remain
outside this finite comparison. The correlated-head construction
tests the former possibility separately.
