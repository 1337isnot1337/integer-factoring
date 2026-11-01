# Faster baselines for the recursive clique evaluator

## Result and scope

The clique evaluator can use fast matrix multiplication in two concrete
places. A three-way contraction improves its initial exponential constant,
and materializing the grouped tables avoids paying their construction cost
for every query made by the recursive child. These changes do not improve
the power `1/2` of the logarithm in the present square-grid recursion.

The construction works over every commutative unital ring. With
`omega0=log_2(7)`, a `K`-vertex clique with `lambda` bits on each edge can
be evaluated in

```math
2^{\omega_0\lceil K/3\rceil^2\lambda}
\operatorname{poly}(K,\lambda)                         \tag{1}
```

ring operations and input-table queries. It uses

```math
O\!\left(2^{2\lceil K/3\rceil^2\lambda}\right)          \tag{2}
```

live ring elements, in addition to polynomial-length loop indices and the
workspace of a primitive ring operation. Thus the limiting coefficient in
(1), normalized by `K^2 lambda`, is `log_2(7)/9=0.311928...`, in place of
the enumeration coefficient `1/2`. An unequal grouping gives a slightly
smaller constant, with its rectangular dimensions charged explicitly.

There is also a separate low-memory improvement. For `K` a power of
two, balanced recursive conditioning uses exactly
`K 2^((K^2-1)lambda/3)` input-table queries and `O(K^2 lambda)` bits
of assignment and control storage. It gives initial coefficient `1/3`
without dense matrices. Section 1.1 proves this result; it can be used
with the original streamed recursion.

This is an ordinary application of Strassen's algorithm, with all tensor
dimensions specified below. It is not a new matrix-multiplication theorem.
The dense matrices require exponential storage. This note does not replace
the separate polynomial-space evaluator or assert a better factoring time
exponent.

The source recursion is Appendix A of
[the pinned `anthropics/formal-math` paper](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/paper.pdf).
Its arithmetic-circuit interpretation over rings is proved in
[scale-one modular transfer](scale-one-modular-transfer.md).

## 1. Three boundary blocks give one matrix product

Let `Q=2^lambda`. A clique assigns a label in `[Q]` to every unordered
edge. Each vertex has an arbitrary table on its `K-1` incident labels.
The desired value is the sum, over all edge assignments, of the product
of the vertex-table values.

Partition the vertices into `A,B,C`, each of size at most
`a=ceil(K/3)`. Write `i,j,k` for the complete label assignments on the
edge blocks `AB,BC,CA`, respectively. Their cardinalities are

```math
D_{AB}=Q^{|A||B|},\qquad D_{BC}=Q^{|B||C|},\qquad
D_{CA}=Q^{|C||A|}.                                    \tag{3}
```

For fixed boundary labels, sum over the internal edges in each group.
This gives matrices `T_A[k,i]`, `T_B[i,j]`, and `T_C[j,k]`. Every original
edge is either internal to exactly one group or occurs in exactly one of
the three boundary blocks. Finite distributivity therefore gives

```math
Z=\sum_{i,j,k}T_A[k,i]T_B[i,j]T_C[j,k]
 =\operatorname{tr}(T_AT_BT_C).                        \tag{4}
```

Zero-pad the matrices to square order
`D=2^(a^2 lambda)`. Padding introduces no nonzero term into (4).
Compute `T_A T_B` and then take its entrywise pairing with the transpose
of `T_C`. Only one matrix product is needed.

Building all entries of the table for a group of size `s<=a` takes at
most

```math
Q^{s(K-s)+\binom{s}{2}}\operatorname{poly}(K,\lambda)
\le Q^{(5a^2-a)/2}\operatorname{poly}(K,\lambda)        \tag{5}
```

ring operations and original table queries. The first exponent counts
both the boundary assignments and the internal assignments; there is no
unaccounted sum inside a matrix entry. Since `omega0>5/2`, (5) is bounded
by the cost of multiplying the padded matrices.

Strassen's seven-product identity uses addition, subtraction, and
multiplication only. Recursive use at order `D=2^d` takes `7^d` scalar
multiplications and `6(7^d-4^d)` scalar additions or subtractions.
It is valid over rings with zero divisors. See
[Strassen's original paper](https://doi.org/10.1007/BF02165411).
The final pairing costs `O(D^2)` operations. This proves (1).

A depth-first implementation retains a constant number of matrices of
order `D`, together with one active recursive multiplication. Its storage
recurrence is `M(D)<=M(D/2)+O(D^2)`, so `M(D)=O(D^2)`. The construction
of an entry in (5) needs only an accumulator and polynomial-length
indices. This proves (2). When one ring value has a fixed `b`-bit
representation, multiply the ring-storage bound by `b` and separately
charge the workspace and bit cost of ring operations.

### 1.1. Recursive conditioning keeps linear assignment space

For `K=2^m`, split the vertices into equal halves `A,B`. Fix the labels
on their `K^2/4` crossing edges. The remaining sum is a product of two
conditional clique values, one on each half. Evaluate those two values
sequentially by applying the same split again. At a singleton, query its
original vertex table with all incident labels already assigned.

If `Q_K` counts original table queries, this gives the exact recurrence

```math
Q_1=1,\qquad Q_K=2^{K^2\lambda/4}\,2Q_{K/2},
\qquad Q_K=K\,2^{(K^2-1)\lambda/3}.                   \tag{5a}
```

The exponent is the geometric sum
`lambda(K^2/4+K^2/16+...+1)`. Finite distributivity proves the recurrence;
no cancellation, positivity, or matrix multiplication is used. Loop
indexing and additions and multiplications add only a polynomial factor.

The coefficient `1/3` is optimal within this direct binary-conditioning
recurrence. If `a+b=K` and `b>=a`, induction on the larger child gives

```math
ab+\frac{b^2-1}{3}-\frac{K^2-1}{3}
=\frac{a(b-a)}3\ge0.
```

Thus an unequal split cannot improve the exponent obtained by summing
all crossing assignments and recursively recomputing both children.
This statement concerns that recurrence, not arbitrary clique algorithms.

Use contiguous binary intervals for the groups. Keep the edge labels
in one packed `binom(K,2)lambda`-bit string, and compute an edge's offset
from its two vertex indices. A stack frame keeps one crossing-assignment
counter of length `s^2 lambda/4` for its current `s`-vertex interval,
constant many ring accumulators, and `O(log K)` index bits. These counter
lengths sum to `O(K^2 lambda)` along a live path. There are `O(log K)`
ring elements and `O(log^2 K)` additional index bits. Thus, over a fixed
finite ring, the evaluator uses `O(K^2 lambda)` work bits, apart from its
read-only input-query interface.

One can retain the source's on-demand grouped-table evaluation and set
its initial coefficient to `(K_0^2-1)/(3K_0^2)`, or simply `1/3`.
The subsequent recurrence still includes its query-inlining term
`3/K_j`. Its leading decay is still `Theta(rho^j)`. This low-memory
variant does not use the materialization step of Section 3. Its packed
assignments also fit the geometric volume accounting of
[linear-space half-rate factoring](linear-space-half-rate.md).

### 1.2. An unequal rectangular split

Let the three sizes be `xK,zK,zK`, where `x=1-2z` and `1/3<=z<1/2`;
assume these sizes are integers. The two boundary dimensions touching
the smaller group are
`r=2^(xz K^2 lambda)`, and the remaining dimension is
`s=2^(z^2 K^2 lambda)`. Therefore

```math
T_A:r\times r,\qquad T_B:r\times s,\qquad T_C:s\times r.
```

Since `s/r` is an integer power of two, split the long dimension into
blocks of length `r`. Form one pair `B_j,C_j` of square blocks at a
time and add `tr(T_A B_j C_j)`. Strassen takes
`(s/r)r^omega0` operations. The total table-formation exponent is at
most `(z-z^2/2)K^2 lambda`, by (5). Hence the leading time coefficient is

```math
\gamma(z)=\max\{z-z^2/2,\
 z^2+(\omega_0-1)z(1-2z)\},                          \tag{5b}
```

with only `O(r^2)` live ring elements. In particular, `z=3/8` gives
groups `K/4,3K/8,3K/8` and

```math
\gamma(3/8)=\frac{6\log_2 7+3}{64}=0.310064\ldots,
\qquad \log_2(r^2)/(K^2\lambda)=3/16.                \tag{5c}
```

The formation coefficient is `39/128`, which is smaller than (5c):
the required inequality is `log_2 7>11/4`, or `7^4>2^11`.
The source's `K0` is divisible by eight, so no rounding is needed.
For the materialized recursion below, one may replace the initial
time coefficient by (5c) and the initial memory coefficient by `3/16`.
The same dominance arguments apply.

The two terms of (5b) meet at
`z=(omega0-2)/(2omega0-7/2)`, where the coefficient is
`0.308902...`. Fixed rational splits can approach this value; such
rounding affects only the initial constant. No asymptotic matrix
multiplication conjecture or logarithmic-power gain is used here.

## 2. A grouped table entry is a weighted two-variable network

The source's grouping lemma collects `q=g/K` nodes into one child
vertex. Each collected node has at most `2^P` possible settings. Once a
child boundary assignment `Theta` is fixed, its table value has the form

```math
\widetilde W(\Theta)
=\sum_{z_1,\ldots,z_q}
  \prod_a W_a^{\Theta}(z_a)
  \prod_{a<b}\Psi_{ab}(z_a,z_b).                       \tag{6}
```

The factors connecting a node to an outside group have been included in
`W_a^Theta`. The internal pair factors remain `Psi_ab`. Padded boundary
bits are checked first. This is exactly the source's grouping formula;
its proof does not require the factors to be nonnegative.

Split these `q` nodes into three clusters of size at most `b=ceil(q/3)`.
Each cluster has at most `2^(bP)` complete settings. For each cluster,
multiply its one-node factors and all pair factors internal to that
cluster. Form a matrix between every two clusters from their cross-pair
factors, placing one cluster weight in each of the three matrices.
Their trace product is (6). Constructing each matrix takes at most
`2^(2bP) poly(q,P)` factor queries and ring operations. Consequently,

```math
\text{one grouped entry costs }
2^{\omega_0\lceil q/3\rceil P}\operatorname{poly}(q,P) \tag{7}
```

operations, in addition to the cost of a primitive factor query. The
matrix storage is `O(2^(2 ceil(q/3)P))` ring elements.

The estimate counts arbitrary pair functions. It assumes no special
algebraic form for a grouped input table.

## 3. Build all child tables once

At level `j+1` of the source recursion, put

```math
K=K_j,\quad g=g_j,\quad N=g^4\lambda,\quad P=g^3\lambda,
\quad q=g/K,\quad \eta_j=\rho+3\zeta_j.
```

The chosen orientation and grouping give a child width `L` with

```math
K^2L\le\eta_jN+\frac c4g^4\log_2(2\lambda)+K^2.        \tag{8}
```

Here `c` is the effective crossing-family constant. The equal group
size `q=g/K` is the conclusion of case (ii) of the source grouping
lemma, not an additional balancing assumption. The implementation is in
[Group/CaseII.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Group/CaseII.lean).

For each crossing-piece pattern, materialize all `K` child tables. Each
table has at most `2^((K-1)L)` entries. Equation (7), with the crossing
factor query cost `2^(c lambda)`, bounds this construction by

```math
K\,2^{(K-1)L+\omega_0\lceil q/3\rceil P+c\lambda}
\operatorname{poly}(g,\lambda).                        \tag{9}
```

The recursive child now reads stored entries. Its cost is **added** to
(9). It is not multiplied by (9). Different piece patterns are still
processed one at a time, and their values are added as in the source.

To state the resulting bound, use the same schedule

```math
K_0=2^{x_0},\quad g_j=2^{3j+x_0}K_j,\quad
K_{j+1}=g_j^2,\quad\zeta_j=2^{-j-5},
```

where `x0>=18`, `2^x0>=c+5`, and `1/2<=rho<1`. Define

```math
\widehat e_0=\omega_0\frac{\lceil K_0/3\rceil^2}{K_0^2},
\qquad
\widehat e_{j+1}=\eta_j\widehat e_j.                   \tag{10}
```

**Proposition.** For some fixed constant `C`, the modified evaluator uses
at most

```math
2^{\widehat e_jK_j^2\lambda}
\bigl(K_j(\lambda+1)\bigr)^{CK_j^2}                   \tag{11}
```

ring operations and input-table queries. Structural computation and
effective piece generation have the same fixed polynomial overhead as
in the source circuit argument.

To verify the leading exponent, ignore only factors already covered by
the displayed polynomial in (11). From (8), the child contributes
`eta_j widehat e_j N`. The contribution from (9) is at most `b_j N`,
where

```math
b_j=\frac{\eta_j+\omega_0/3}{K_j}
    +\frac{2\omega_0}{3g_j}+\frac c{g_j^4}.            \tag{12}
```

The ceiling inequality used here is `ceil(q/3)<=q/3+2/3`.
The logarithmic term and padding term in (8), the number of piece
patterns, the grouping search, and the child polynomial overhead all
fit `bigl(g^2(lambda+1)bigr)^(C g^4)`, by the same uniform induction as
the source. In particular, the additional table-entry factor in (9)
adds only `O(c g^4 log(2lambda)/K)` to that polynomial exponent.

Thus the leading recurrence is the maximum of
`eta_j widehat e_j` and `b_j`. Under the stated schedule, `b_j<4/K_j`.
Also `widehat e_0>1/4`, `widehat e_j>=widehat e_0 rho^j`, and
`K_j>=K_0^(2^j)>=2^(j+5)`. Hence

```math
b_j<4/K_j\le2^{-j-3}
 <\widehat e_0\rho^{j+1}\le\eta_j\widehat e_j.
```

The recursive term dominates, proving (10)--(11). The bound
`widehat e_j<1`, needed when absorbing the width slack, follows from
`sum_j 3zeta_j=3/16` and

```math
\widehat e_j\le\widehat e_0\rho^j
\exp\!\left(\frac3{16\rho}\right).
```

This is an arithmetic-operation and table-query improvement. An entry
read from a materialized child table is charged as one query in (11).
A finite-tape bit-complexity use must also account for that storage-access
pattern, as well as the representation and arithmetic costs of the ring.
Equation (11) alone does not assert its leading constant for bit time.

## 4. Space of the modified recursion

Let `M_j(lambda)` count live ring elements, excluding externally supplied
input-table storage. For one pattern at the next level, retain its `K`
child tables, then run the child. The local matrices can be released
before the child starts. Therefore the memory bound is

```math
M_{j+1}(\lambda)
\le K2^{(K-1)L}
   +\max\{M_j(L),\ O(2^{2\lceil q/3\rceil P})\}
   +\operatorname{poly}(g,\lambda).                  \tag{13}
```

Piece-generation workspace is additionally bounded by its `2^(c lambda)`
time bound. Pattern iteration is sequential, so there is no factor equal
to the number of patterns in (13). Parent tables remain available while
the child runs, as explicitly counted by the first term.

For completeness, the same induction gives a bound of the form (11)
on memory with leading coefficient

```math
s_0=2\frac{\lceil K_0/3\rceil^2}{K_0^2},\qquad
s_{j+1}=\eta_js_j.                                   \tag{14}
```

Indeed the two extra leading coefficients in (13) are at most
`eta_j/K_j` and `2/(3K_j)+4/(3g_j)`, while coefficient-generation space
contributes `c/g_j^4`. They are dominated by `eta_j s_j` under the
stated large initial group count. Finite sums of memory bounds and the
width slack are absorbed in the polynomial factor. This proves an
explicit exponential-space bound; it does not prove polynomial space.

## 5. What the batching changes, and what it leaves open

Equations (10) and (14) have exact products

```math
\widehat e_j=\widehat e_0\prod_{i<j}(\rho+3\zeta_i),
\qquad s_j=s_0\prod_{i<j}(\rho+3\zeta_i).
```

Since `sum zeta_i` is finite, both coefficients are `Theta(rho^j)`.
The unchanged schedule has `log K_j=Theta(2^j)`. Consequently, the
logarithmic saving remains `alpha=log_2(1/rho)`, which is `1/2` when
`rho=1/sqrt(2)`.

Section 3 uses the original logarithmic width slack in (8). For a
family supplied with the weaker moment-based slack, the
[volume-cutoff transfer](subexponential-piece-transfer.md) adds its
bounded product of width errors. That extension retains the same
logarithmic power; it should not be confused with the exact product
formula (10) under the stronger interface.

The reason is specific: the child still has volume `K^2L` bounded by
`(rho+o(1))N`. Matrix multiplication accelerates the initial scalar
evaluator and the construction of its input tables. It does not change
the width argument passed to the scalar recursive call. Even eliminating
the local construction and pattern-enumeration costs from this displayed
bound would leave that recursive term.

This is a limitation of the retained recurrence and interface. It is
not a lower bound for clique evaluation or factoring, and does not
exclude an algorithm that evaluates a structured family of child calls
together. A proposed tensor-valued child interface must account for its
boundary dimensions and show how its factors remain local. Treating all
boundary entries as free outputs would omit precisely the dimensions
in (3), (5), and (9).

There is a smaller elementary query lower bound worth separating from
this recurrence discussion. Fix every non-star edge of a clique to zero
by indicator weights, and leave only the table at the star center
arbitrary. The clique value is the sum of all `2^((K-1)L)` entries of
that table. Any exact deterministic entry-query algorithm over a
nontrivial ring must inspect every entry in the worst case. An
unqueried entry can be changed by one while keeping the entire observed
transcript fixed. This justifies a possible table-size cost; its exponent
is much smaller than the recursive exponent and supplies no half-rate
barrier.

## Reproduction

```sh
python tools/check_fast_matrix_clique.py \
  --output evidence/fast-matrix-clique.json
```

The checker compares the three-boundary construction with direct
edge-assignment sums, and compares (6)--(7) with direct setting sums.
It also checks the packed recursive-conditioning evaluator and the
streamed rectangular-block trace against those independent sums.
It tests finite fields and rings with zero divisors. It also checks the
seven-product arithmetic counts and finite exact inequalities used to
separate the table-construction term from the child term. These finite
checks support the implementations; the uniform bounds above have
independent proofs.
