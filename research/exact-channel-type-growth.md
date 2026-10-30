# Exact exponential growth of matching channel types

## Result and scope

The [weighted channel bound](twisted-channel-walk-transfer.md) can lose
information about connectivity and exact type compatibility. The
[component-path refinement](component-channel-walk-transfer.md) removes
some connectivity losses. This note gives an exact variational formula
for the exponential growth of the underlying matching-type count.

The additional condition is concrete: a pair of edge supports must admit
one pair of actual walks with the same piece totals. Once such a pair
exists, every pair of stationary integer flows with matching totals can
be added repeatedly. An elementary count of Euler trails then realizes
the predicted entropy up to polynomial factors.

Thus stationary flows give lower bounds on **actual matching-type
counts**, provided a connected integer base is also supplied. Without
that base they only bound the weighted-supervector method, as the
[flow-dual controls](channel-flow-dual-certificates.md) demonstrate.

The theorem concerns a sufficient rank-count bound. A lower bound on
that count is not a lower bound on the tensor ranks that it upper-bounds.
No improved factoring exponent follows from this note.

## 1. Fixed labeled graphs and exact supports

Let `G_R,G_C` be finite directed multigraphs. Each edge has one of the
same `m` piece labels. Parallel edges remain distinct. First fix an
initial and a terminal vertex in each graph. Let `R_nu(k)` and `C_nu(k)`
count walks of `k` edges with piece-count vector `nu`, and put

```math
\mathcal C_k=\sum_{|\nu|=k}\sqrt{R_\nu(k)C_\nu(k)}.
\tag{1}
```

The same result holds for finitely many initial and terminal choices
with fixed nonnegative weights. Discard zero weights and take the
maximum of the resulting exponential rates. In particular, terminal
bilinear ranks in the [product-space automata](scalarized-channel-automata.md)
are permitted. A vertex-labeled channel word is converted to an
edge-labeled walk by adding an initial edge carrying its first label.

For an edge vector `a`, write `out_v(a)` and `in_v(a)` for its outgoing
and incoming totals, and `type(a)` for its piece totals. An **admissible
support pair** consists of nonempty edge sets `F_R,F_C`, the prescribed
endpoints, and positive integer vectors `a^0,b^0` on those supports such
that:

- each support is weakly connected after unused vertices are removed,
  and contains its initial and terminal vertices;
- `out(a^0)-in(a^0)=1_s-1_t` in the row graph, and the analogous
  equation holds in the column graph;
- `type(a^0)=type(b^0)`.

All edges of the selected supports have positive base count. The two
base lengths are equal because every edge has exactly one piece label.
Connectivity is essential: balance alone allows disconnected cycles
that no one walk can visit.

These conditions are equivalent to the existence of two walks with
the same type and exactly these edge supports. Necessity is immediate.
Sufficiency is the elementary Euler-trail argument in Section 3.
There are only finitely many support pairs and endpoint choices.
An integer base is a finite positive certificate of admissibility;
discarding a proposed pair requires an actual infeasibility argument.

## 2. The stationary entropy formula

For an admissible support pair `sigma`, define the compact polytope

```math
\mathcal P_\sigma=\left\{(x,y):
\begin{array}{l}
x,y\ge0,\quad |x|=|y|=1,\\
\operatorname{out}(x)=\operatorname{in}(x),\quad
\operatorname{out}(y)=\operatorname{in}(y),\\
\operatorname{type}(x)=\operatorname{type}(y)
\end{array}\right\}.
\tag{2}
```

Flows use the selected supports. They may be zero on some support
edges, even though the base vectors are positive on every such edge.
This permits a walk to spend arbitrarily long in several components
while traversing a connecting edge only a fixed number of times.

Using natural logarithms and `0 log 0=0`, set

```math
\mathsf h(x)=\sum_v\operatorname{out}_v(x)
                         \log\operatorname{out}_v(x)
                  -\sum_e x_e\log x_e,
\qquad
\Lambda_\sigma=
\max_{(x,y)\in\mathcal P_\sigma}
             \frac{\mathsf h(x)+\mathsf h(y)}2.
\tag{3}
```

Empty polytopes contribute exponential base zero. If there are no
nonempty polytopes, the maximum below is defined to be zero.

**Theorem.** For the fixed finite graphs and endpoint data,

```math
\boxed{\limsup_{k\to\infty}\mathcal C_k^{1/k}
   =\max_{\sigma\text{ admissible}}\exp(\Lambda_\sigma).}
\tag{4}
```

This is a limsup statement. Exact congruences can make the count zero
on infinitely many lengths. Replacing limsup by a limit would be false.
The theorem is an exact finite-dimensional description of the counting
rate; it does not assert that its entropy maximization is inexpensive.

## 3. Counting one connected edge-multiplicity class

Fix a support with `V` used vertices and `E` used edges, endpoints
`s,t`, and positive integer counts `a` satisfying Euler balance. Let
`k=|a|`, and let `W(a)` be the number of edge sequences with precisely
these counts. Put

```math
M(a)=\frac{\prod_v\operatorname{out}_v(a)!}{\prod_e a_e!},
\qquad
\mathsf H(a)=\sum_v\operatorname{out}_v(a)
                         \log\operatorname{out}_v(a)
                  -\sum_e a_e\log a_e.
```

An outgoing queue at each vertex specifies an ordering of its outgoing
edge multiset. There are exactly `M(a)` queue choices. Starting at `s`
and consuming the next edge whenever a vertex is visited determines
at most one walk. Therefore `W(a)<=M(a)`.

Every used vertex has a directed path to `t`. One way to see this is
to contract strongly connected components. A sink component not
containing `t` has no outgoing flow, and Euler balance forces it to
have no incoming flow either (a sink containing `s` but not `t` would
have positive net outgoing flow, which is already impossible).
This contradicts weak connectivity unless it is the sole component.
The same argument applies when `s=t`, where a nontrivial condensation
would have an impossible sink.

Choose a directed tree of support edges toward `t`. At every vertex
other than `t`, reserve one copy of its tree edge as the last outgoing
edge in the queue. The fraction of queue choices meeting this condition
at vertex `v` is `a_(tree edge)/out_v(a)>=1/(k+1)`.

Every resulting queue system gives a complete Euler trail. Balance
prevents it from stopping at a vertex other than `t`. If it stopped
at `t` before using every edge, the unused edge multiset would be
balanced. At each vertex with an unused edge, the reserved last tree
edge would remain unused. Its target would then also have unused
outgoing flow. Following the tree would show that `t` still has an
unused outgoing edge, contradicting termination. Consequently

```math
(k+1)^{-V}M(a)\le W(a)\le M(a).
\tag{5}
```

This also proves the sufficiency of connected positive Euler counts.
It avoids any assumption about distinguishable repeated copies of
an edge: the queues use edge categories with their stated multiplicities.

The multinomial probability bound gives `M(a)<=exp(H(a))`. The
elementary integral estimates for `log(n!)` give
`n!>=(n/e)^n` and `n!<=en(n/e)^n` for `n>=1`.
Since `en<=(k+1)^2` when `1<=n<=k`, they also give

```math
(k+1)^{-2E}\exp(\mathsf H(a))\le M(a)
                                      \le\exp(\mathsf H(a)).
```

Together with (5), this proves the uniform estimate

```math
\boxed{(k+1)^{-(2E+V)}\exp(\mathsf H(a))
                 \le W(a)\le\exp(\mathsf H(a)).}
\tag{6}
```

All constants depend only on the fixed finite support. Also
`H(a)=k h(a/k)`, because the terms `k log k` cancel.

## 4. Upper bound in the theorem

Partition walks by their exact edge support and endpoints. For each
piece type, the square root of the product of the two finite sums
is at most the sum of the square roots for all support pairs. It is
also at least the contribution of any one pair. The total growth
rate is therefore the maximum over the admissible support pairs.

For one support pair, there are at most
`(k+1)^(E_R+E_C)` possible pairs of edge-count vectors. Thus (6)
bounds its matching-type sum by a fixed polynomial in `k` times

```math
\max_{a,b}\exp\bigl((\mathsf H(a)+\mathsf H(b))/2\bigr),
\tag{7}
```

where the maximum is over actual same-type counts of length `k`.
To justify this even when several flows have the same type, use
`sqrt((sum_a u_a)(sum_b v_b))<=sum_(a,b) sqrt(u_a v_b)`.

Along any unbounded sequence of feasible lengths, a convergent
subsequence of `(a/k,b/k)` exists in the two unit simplices. Its
limit is stationary, since the endpoint imbalances divided by `k`
tend to zero. Its two piece totals still agree. Hence the limit
belongs to (2). The entropy functions are continuous at zero,
so (7) has exponential rate at most `exp(Lambda_sigma)`.
If (2) is empty, no unbounded sequence of feasible lengths exists.
This proves the upper bound in (4).

## 5. Matching lower bound from an integer base

Fix an admissible support pair and its positive base vectors
`a^0,b^0`, of common length `k_0`. The polytope (2) is rational,
and its rational points are dense: express any point as a convex
combination of rational vertices and approximate the coefficients
by nonnegative rational coefficients summing to one.

Choose a rational point `(x,y)` whose objective in (3) is as close
to the maximum as desired. Clear a common denominator `D` and put
`u=Dx,v=Dy`. Then `u,v` are stationary nonnegative integer vectors
of total `D`, and their piece totals agree. For every integer `n>=0`,

```math
a^{(n)}=a^0+nu,\qquad b^{(n)}=b^0+nv
\tag{8}
```

remain positive on every support edge, satisfy the endpoint balances,
and have the same piece totals and length `k_n=k_0+nD`.
Their supports stay connected, even if the stationary flows vanish
on the connecting edges. Section 3 therefore supplies actual walks.

The common type contributes at least
`sqrt(W(a^(n)) W(b^(n)))` to (1). Equation (6), followed by
`a^(n)/k_n -> x` and `b^(n)/k_n -> y`, gives a limsup base at least
`exp((h(x)+h(y))/2)`. Taking rational points arbitrarily close to
the maximum proves the lower bound in (4).

## 6. An exact integer lower certificate

For a stationary integer flow `u` of total `D`, define

```math
K(u)=\frac{\prod_v\operatorname{out}_v(u)^{\operatorname{out}_v(u)}}
             {\prod_e u_e^{u_e}},\qquad 0^0=1.
```

Supply an admissible integer base and stationary flows `u,v` on
those supports, with equal total `D` and equal piece totals. The
preceding proof gives the exact certificate

```math
\limsup_k\mathcal C_k^{1/k}\ge
                         (K(u)K(v))^{1/(2D)}.
\tag{9}
```

The right side is compared with proposed rational or algebraic
thresholds by integer powers. This is the additional connectivity
and arithmetic hypothesis needed to turn a flow-dual method bound
into an actual count-growth bound.

For example, take two parallel loops with the same label in the
row graph and three parallel loops in the column graph. Use base
and stationary counts `(3,3)` and `(2,2,2)`, with `D=6`.
Their entropy ratios are `2^6` and `3^6`. Thus (9) gives base
`sqrt(6)`, exactly matching `C_k=sqrt(6)^k`.

A connecting edge may have positive base count but zero stationary
count. In the row graph, put two `A` loops at the left vertex, two
`B` loops at the right, and one `B` bridge from left to right; use
those two endpoints. In the column graph, put two `A` edges from
zero to one and two `B` edges from one to zero; start at one and
finish at zero. Base counts `(1,1,1,1,1)` and `(1,1,1,2)` have
common length five. Stationary increments `(1,1,1,1,0)` and
`(1,1,1,1)` have common length four and entropy product `2^8`.
At every odd length `k`, the matching type has row count `2^(k-1)`
and column count `2^k`. Thus the exact limsup base is two. The bridge
is used once, while repeated stationary flows occur on the loops.
For this example the optimal base is an infimum of positive-vector
certificates: a vector cannot absorb a positive bridge while keeping
both recurrent diagonal blocks at their spectral radius exactly.

## 7. Arithmetic controls and practical use

Give every edge the same piece label. Let both graphs be directed
two-cycles. In the row graph, require start and finish at vertex zero;
only even lengths occur. In the column graph, require start zero and
finish one; only odd lengths occur. There is no admissible support
pair with matching counts, and `C_k=0` for every length. Nevertheless
both individual spectral radii are one. Component-path refinement
alone cannot detect this exact congruence obstruction.

If both graphs instead require a closed walk, the count is one at
even lengths and zero at odd lengths. Its limsup base is one, while
its ordinary limit does not exist. This verifies why the theorem
uses a limsup and why an arithmetic progression suffices in its proof.

For an upper-bound application, the original matching-type sum may
be covered by all admissible support pairs. Each pair can receive its
own reciprocal weights and rational positive-vector certificate,
exactly as in the component-path refinement. The number of pairs is
fixed by the finite seed. A proposed exclusion of a support pair must
be justified by its endpoint, connectivity, or integer type equations;
finite sampling alone is insufficient.

## 8. Optimality of independent reciprocal weights on exact supports

The entropy formula also identifies precisely what a sufficiently
refined spectral certificate can attain. For a selected support `F`,
let `T_F(z)` be its adjacency matrix, summing `z_i` for every edge
with piece label `i`. Then, for every admissible support pair,

```math
\boxed{\exp(\Lambda_\sigma)=
 \inf_{z_i>0}\sqrt{\operatorname{spr}T_{F_R}(z)
                         \operatorname{spr}T_{F_C}(z^{-1})}.}
\tag{10}
```

When the matching stationary polytope is empty, the right side is
zero. If either support is acyclic, its spectral radius is zero
for every weighting, and the assertion follows directly.

Here is a proof for the remaining case. For a graph with a cycle,
its stationary unit-flow polytope `P_F` is nonempty. The pressure
identity is

```math
\log\operatorname{spr}T_F(e^u)
 =\max_{x\in P_F}\bigl(\mathsf h(x)
                               +u\cdot\operatorname{type}(x)\bigr).
\tag{11}
```

For the upper inequality, a positive supervector `T h<=beta h`
and the weighted arithmetic-geometric mean inequality give
`h(x)+u.type(x)<=log beta`; the vertex-vector factors cancel by
stationarity. Let `beta` decrease to the spectral radius. Positive
supervectors for every strictly larger `beta` are supplied by
`(beta I-T)^(-1) 1`.

For the reverse inequality, choose a strongly connected component
whose spectral radius `rho` is maximal. Let `l,r` be its positive
left and right Perron vectors, normalized by `l.r=1`. For each edge
`e:s->t` in that component, put
`x_e=l_s exp(u_(label e)) r_t/rho`. This is a stationary unit flow.
Substitution in its entropy gives equality in (11). Extend it by
zero outside the component. More generally, a stationary flow
decomposes across recurrent components; entropy is homogeneous on
disjoint components, so mixtures cannot exceed the largest pressure.

Using (11), the logarithm of the right side of (10) is the infimum
over `u` of

```math
\max_{x\in P_{F_R},\,y\in P_{F_C}}
 \left\{\frac{\mathsf h(x)+\mathsf h(y)}2
       +\frac u2\cdot(\operatorname{type}(x)-\operatorname{type}(y))
 \right\}.                                      \tag{12}
```

The entropy objective is concave and continuous on a compact
rational polytope. Concavity follows vertex by vertex from the
perspective of the ordinary probability entropy: the contribution
at a vertex is its outgoing mass times the entropy of its conditional
outgoing-edge distribution. Its value as a function of the type difference
has a closed convex hypograph. If type difference zero is feasible,
separate the point `(0,Lambda_sigma+epsilon)` from that hypograph.
The separating functional has a strictly positive vertical
coefficient: the hypograph contains `(0,Lambda_sigma)`, so a
horizontal functional could not strictly separate these two points.
Dividing by its vertical coefficient supplies a finite multiplier
in (12) with value at most `Lambda_sigma+epsilon`. The reverse
inequality follows simply by restricting the maximum to type
difference zero. Letting `epsilon` tend to zero proves (10), without
an interior-point assumption or attainment of the infimum.

If zero is not feasible, the two compact convex type sets are
strictly separable. Scaling the separating multiplier sends (12)
to minus infinity; the bounded entropy term does not prevent this.
This proves the zero-base case as well.

Consequently the exact-support version of the reciprocal-weight
method loses **no exponential factor in the matching-type count**.
Each admissible support pair may use its own weights. Given any
strict upper base larger than (4), finitely many choices suffice.
Continuity permits positive rational piece weights. For a rational
transition matrix, choose rational `beta` strictly above its
spectral radius and use the positive rational vector
`h=(beta I-T)^(-1) 1`; it satisfies `T h=beta h-1`.
Thus the strict upper rate is witnessed by finite rational
positive-vector inequalities, after the support exclusions have
been justified. This does not remove the need to check exact
integer admissibility, and does not make the underlying tensor
rank bound an equality.

## Verification

The checker `tools/check_exact_channel_growth.py` enumerates small
edge-count classes, verifies the exact queue and entropy inequalities,
checks integer base/flow certificates, and includes disconnected and
parity controls. The all-size equality (4) follows from the proof,
not from these bounded computations.
