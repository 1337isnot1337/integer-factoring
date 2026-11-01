# Exact channel-support feasibility by bounded binary carries

## Result and scope

The admissibility condition in the
[exact matching-type growth theorem](exact-channel-type-growth.md) is
decidable by a finite, elementary search. The search either supplies the
required positive integer edge counts or supplies a finite closed set of
carry states proving that no such counts exist. No integer-programming
oracle or unstated bound on a smallest solution is required.

Together with the entropy theorem's rational lower certificates and its
reciprocal-weight upper certificates, this makes the exact exponential
matching-type count rate uniformly computable from finite labeled graphs
and prescribed allowed endpoints. The algorithms can be expensive. This
is an effectiveness statement about a channel-count bound, not an improved
factoring exponent or a lower bound on tensor ranks.

## 1. A finite carry graph for a nonnegative integer system

Let `A` be an integer `r` by `m` matrix and let `b` be an integer vector.
Consider the system

```math
Ax=b,\qquad x\in\mathbb Z_{\ge0}^{m}.
\tag{1}
```

For each equation put

```math
M_i=\max\left\{|b_i|,\sum_{j=1}^{m}|A_{ij}|\right\},
\qquad
\mathcal B=\prod_{i=1}^{r}\{-M_i,-M_i+1,\ldots,M_i\}.
\tag{2}
```

The carry graph has vertex set `B`, initial vertex `-b`, and terminal
vertex zero. For each binary digit vector `u in {0,1}^m`, include the
directed edge

```math
c\longrightarrow c'=(c+Au)/2
\tag{3}
```

exactly when all coordinates of the numerator are even. Every such
successor remains in `B`, since

```math
|c'_i|\le\frac{|c_i|+\sum_j|A_{ij}|}{2}\le M_i.
```

**Theorem.** System (1) has a solution if and only if zero is reachable
from `-b` in this finite graph.

**Proof.** A path of length `L` with edge labels `u_0,...,u_(L-1)` satisfies

```math
2^L c_L=-b+A\sum_{j=0}^{L-1}2^j u_j.
\tag{4}
```

If `c_L=0`, the nonnegative integer vector
`x=sum_j 2^j u_j` solves (1).

Conversely, write any nonnegative solution in binary as
`x=sum_(j<L)2^j u_j`. The carries in (4) are integers at every step:
`-b+A sum_(j<s)2^j u_j` equals minus `A` times the remaining, divisible
part of `x`. Thus every required edge is permitted. The invariant bound
above keeps the entire path in `B`, and the final carry is zero. This
proves the equivalence.

Writing

```math
S=|\mathcal B|=\prod_i(2M_i+1),
```

a shortest path has at most `S-1` edges unless the initial state is
already zero, in which case the zero solution suffices. Consequently a
feasible system has a solution with every coordinate less than
`2^(S-1)` when `S>1`. The path, rather than this coarse bound, is what the
algorithm constructs. Breadth-first search inspects at most `S 2^m`
candidate transitions, each using finitely many integer operations on
the displayed input and carry bounds.

## 2. Finite certificates for both outcomes

A positive certificate consists of the solution vector itself, or of
the path and its binary edge labels. Equation (1) verifies the former
directly; equation (4) verifies the latter.

A negative certificate is any finite subset `C` of `B` with these three
properties:

1. `-b` belongs to `C`;
2. zero does not belong to `C`;
3. every permitted transition (3) from a member of `C` has its successor
   in `C`.

Every reachable carry then stays in `C`, so the theorem proves
infeasibility. Exhaustive reachability returns exactly such a set when
it fails to reach zero. Checking this certificate needs only the finite
integer transition calculation. It does not require trusting a claimed
search timeout or a numerical infeasibility result.

## 3. Application to positive Euler counts on exact supports

Fix row and column edge supports and their endpoints. First check the
finite graph requirements from the entropy theorem: both supports are
nonempty, weakly connected after unused vertices are removed, and
contain their prescribed endpoints.

Let `a,b` denote the edge-count vectors. Their equations are:

```math
\operatorname{out}(a)-\operatorname{in}(a)=1_{s_R}-1_{t_R},
\quad
\operatorname{out}(b)-\operatorname{in}(b)=1_{s_C}-1_{t_C},
\quad
\operatorname{type}(a)-\operatorname{type}(b)=0.
\tag{5}
```

Every selected edge must occur at least once. Write
`(a,b)=1+x`, where `x` has one coordinate per selected edge. If the
matrix and right side of (5) are `D` and `e`, respectively, positivity
becomes precisely

```math
Dx=e-D\boldsymbol1,\qquad x\ge0\text{ integral}.
\tag{6}
```

Apply the carry theorem to (6). A solution gives positive Euler counts
with equal piece totals. The connected Euler-trail lemma in the entropy
theorem turns these into two actual same-type walks. A closed carry set
rules out this support pair exactly.

For `E` total selected edges and `r` displayed equations, every matrix
entry in (6) has absolute value at most one, and every carry bound is
at most `E+1`. Thus `(2E+3)^r` is an explicit finite upper bound on the
number of carry states. Redundant balance equations may be retained.
The bound is intentionally coarse; no efficient general feasibility
claim follows.

### A mixed-endpoint parity obstruction

Give every edge the same piece label. The row graph has the two-edge
cycle `0 -> 2 -> 0` and the final edge `0 -> 1`, with start zero and end
one. Positive Euler counts are necessarily `(q,q,1)` with `q>=1`, so
every row walk has odd length `2q+1`.

The column graph is the two-edge cycle `0 -> 1 -> 0`, with both endpoints
zero. Its positive counts are `(s,s)`, with even length `2s`. Therefore
the support pair is infeasible, although the stationary unit-flow
polytopes have matching piece totals. Even the positive real version
of the base equations is feasible, for example `q=1,s=3/2`.

The carry checker produces and verifies a finite negative certificate
for exactly these graphs. It also checks the two-label bridge example
from the entropy theorem, where a positive integer base exists although
the limiting stationary flow on the bridge is zero.

## 4. The exponential count rate is a computable real

There are finitely many exact support pairs and allowed endpoint
choices. Sections 1--3 decide which pairs are admissible and retain a
positive integer base for each surviving pair.

For one admissible support pair, the entropy theorem supplies two
families of exact bounds on its rate `B_sigma`.

* A stationary pair of nonnegative integer flows `u,v` of common total
  `D>0` and equal piece totals gives the lower bound
  `(K(u)K(v))^(1/(2D))`, with the integer-product ratios `K` defined in
  that theorem. Rational stationary flows are dense in its compact
  rational polytope, so these lower bounds approach `B_sigma` whenever
  that polytope is nonempty. Zero is always a valid lower bound.
* Positive rational piece weights `z`, positive rational vectors
  `v_R,v_C`, and positive rationals `beta_R,beta_C` satisfying
  `T_R(z)v_R<=beta_R v_R` and
  `T_C(z^(-1))v_C<=beta_C v_C` give the upper bound
  `sqrt(beta_R beta_C)`. Section 8 of the entropy theorem proves that
  these rational certificates approach `B_sigma` from above, including
  the case `B_sigma=0` when the stationary matching polytope is empty.

Every certificate check uses integer or rational arithmetic. For
example, a rational proposed lower bound `l>=0` is checked by raising
it to `2D` and comparing with `K(u)K(v)` after clearing denominators.
A rational upper bound `U` is checked from
`U^2>=beta_R beta_C`.

Enumerate both certificate families in parallel. Retain the largest
certified rational lower bound and the smallest certified rational
upper bound. Their gap eventually falls below any prescribed positive
rational tolerance. Taking the maximum over the finitely many
admissible supports computes the overall rate. If there are no
admissible supports, the rate is zero directly.

This procedure does not need to decide equality with a proposed
threshold. Nor does it claim a polynomial running time or a uniform
simple convergence rate. It gives a terminating algorithm for each
requested precision using only finite data and exact certificates.

## Reproduction

Run from the repository root:

```sh
python3 tools/check_channel_support_feasibility.py --output /tmp/channel-support-feasibility.json
```

The checker compares the carry decision with an independent extended-gcd
and interval calculation for all 637 signed two-variable equations in
its stated range. It also compares bounded-digit reachability with direct
enumeration in 200 small multiequation systems, validates every positive
witness and every negative closed-set certificate, and records the two
channel-support examples above. The finite checks support the separate
all-size proof; they are not used to assume a solution-size bound.
