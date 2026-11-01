# Exact rank-sum searches and dual label restrictions

## Results and scope

The sufficient finite-seed criterion is

\[
 E_d=\sum_j B_j,\qquad
 S=\sum_j\sqrt{r(B_j)c(B_j)}<d^{\sqrt2}.                           \tag{1}
\]

The amplification theorem in
[finite-seed-amplification.md](finite-seed-amplification.md) turns
such a certificate over a fixed field into a uniform contraction
strictly better than \(1/\sqrt2\). A small maximum rank product of
the original pieces is insufficient.

This note gives two constructive operations for seeking (1):
four separate dual label restrictions, and exact optimization over a
finite coefficient span. The implementation encloses the optimal
rank sum with integer bounds and returns an explicit decomposition.
It uses no floating-point rank or shortest-path calculation.

Applying it to the 24-label
[split-permutation seed](split-permutation-seed.md) gives an exact
optimum \(S=325\) in its uniform coefficient span. This does not
meet (1). The additional finite searches below are also negative,
with their scopes stated explicitly. The full 141-dimensional
coefficient span has not been exhaustively searched. A subsequent
[single-cycle block argument](permutation-span-cycle-obstruction.md)
excludes the half threshold in that entire original span analytically.

## 1. Four dual label maps preserve a crossing

Use the tensor convention
\[
 E_D(a,b,u,y)=[a=b][u=y].
\]
The row cut is \((a,y)\mid(b,u)\), and the column cut is
\((a,u)\mid(b,y)\); the crossing matrix is the identity in both.

Let \(P_A,P_B,P_U,P_Y\) be \(D\)-by-\(d\) matrices over the
coefficient field, satisfying

\[
 P_A^{\mathsf T}P_B=I_d,\qquad
 P_U^{\mathsf T}P_Y=I_d.                                         \tag{2}
\]

Apply these four linear maps to the corresponding external label
coordinates of every tensor in a crossing decomposition. Summing
the transformed tensors gives \(E_d\), because the two equality
factors contract to the two identities in (2).

More explicitly, the row matrix of a tensor \(T\) becomes

\[
 (P_A\otimes P_Y)^{\mathsf T}\,
       T_{\rm row}\,
 (P_B\otimes P_U),                                               \tag{3}
\]

and its column matrix becomes

\[
 (P_A\otimes P_U)^{\mathsf T}\,
       T_{\rm col}\,
 (P_B\otimes P_Y).                                               \tag{4}
\]

Thus neither flattening rank can increase. The new label dimension is
smaller, so (1) must still be checked with that new dimension.
These statements are independent of how the seed was constructed.

For a fixed finite seed, the four matrix contractions and exact rank
factorizations have finite computable cost. They are therefore
compatible with the effective finite-seed amplification theorem.

### Two concrete choices over \(\mathbb F_2\)

Suppose the old labels are split into disjoint nonempty groups.
The vector that sums the coordinates in a group has dot product one
with itself exactly when the group size is odd. Sums over disjoint
odd groups are therefore orthonormal. Taking all four maps equal to
their embedding preserves crossing.

There is also a construction for even groups. Let \(P\) sum each
group and let \(Q\) select one representative from each group.
Then \(P^{\mathsf T}Q=I\), regardless of the group sizes. Use \(P,Q\)
as either ordered dual pair in (2). This gives asymmetric label
restrictions; it does not require a self-dual label basis.

For example, sums over the left or right cosets of a three-cycle
subgroup in \(S_4\) give eight orthonormal labels from the original
24 permutation labels.

## 2. Exact optimization in a coefficient span

Let \(V\) be an \(m\)-dimensional \(\mathbb F_2\)-space of tensors
containing \(E_d\), and choose a basis. Identify \(V\) with
\(\mathbb F_2^m\). For every nonzero coordinate vector \(v\), compute
both actual flattening ranks and assign cost

\[
                         \kappa(v)=\sqrt{r(v)c(v)}.
\]

Consider the graph on \(\mathbb F_2^m\) having an edge from \(x\)
to \(x+v\) of cost \(\kappa(v)\). A path from zero to the
coordinate vector of \(E_d\) is exactly a decomposition of \(E_d\)
into tensors in \(V\), with the same rank sum. Conversely every
decomposition gives such a path. Thus the shortest path is the
optimal rank sum among all decompositions using \(V\).

An optimal decomposition can contain no nonempty zero-sum subfamily:
removing one lowers its positive cost. In particular, an optimum
uses at most \(m\) linearly independent nonzero tensors. The graph
formulation nevertheless avoids having to guess that subfamily.

### Integer enclosures

Fix a positive integer scale \(M\). For each \(v\ne0\), compute

\[
 \ell_v=\left\lfloor M\sqrt{r(v)c(v)}\right\rfloor,\qquad
 u_v=\left\lceil M\sqrt{r(v)c(v)}\right\rceil.
\]

Integer square roots give these values exactly. Run integer Dijkstra
twice, once with each cost family. If the resulting distances are
\(L,U\), then

\[
                       L/M\le S_{\min}\le U/M.                   \tag{5}
\]

The upper search supplies a concrete tensor decomposition.
The lower search considers every tensor in the entire declared
span, so its lower bound applies to arbitrary linear recombinations
inside that span.

There is no numerical comparison with an irrational threshold.
For a positive rational \(a/b<\sqrt2\), the integer inequality

\[
                         U^b<M^b d^a                            \tag{6}
\]

certifies (1). To certify a negative result within the span, one may
choose \(a/b>\sqrt2\) and verify

\[
                         L^b\ge M^b d^a.                        \tag{7}
\]

The implementation includes both exact power tests. Its declared
exhaustive-search limit is checked before enumerating a span.

For repeated searches inside a larger known coefficient space, the
rank cache can use exact coordinates in that larger space. This is
an injective linear encoding, not a probabilistic hash. It avoids
storing a large packed matrix separately for every cached rank.

## 3. The universal three-pairing baseline

Work over \(\mathbb F_2\) in this section.

Write

\[
 R=[a=y][b=u],\qquad
 C=[a=u][b=y],\qquad
 H=E_d+R+C.
\]

The rank pairs of \(R,C\) are \((1,d^2)\) and \((d^2,1)\).
In either cut, \(H\) is the sum of the identity, the swap, and
the outer product of the diagonal vector with itself.
Each unordered off-diagonal pair gives a two-by-two all-ones
block of rank one; the diagonal coordinates give one further
all-ones block. Therefore

\[
 r(H)=c(H)=\binom d2+1.
\]

This gives the general decomposition \(E_d=R+C+H\), with rank sum

\[
 S_{\rm base}=2d+\binom d2+1=\frac{(d+1)(d+2)}2.                  \tag{8}
\]

This baseline is included to identify what the finite search found.
It is not presented as a new crossing identity or a half-exponent
improvement.

## 4. Exact 24-label computations

For the split-permutation family on \(S_4\), denote its coefficient
tensors by \(T_{S,T}\), \(S,T\subseteq[4]\).
There are 252 nonzero pieces. Their full binary linear span has
dimension 141.

Set all \(P\)- and \(Q\)-variables equal to \(z\), and call the
resulting coefficient span \(V_{\rm unif}\). It has dimension six.
Its constant and leading coefficients are \(R,C\), and its
coefficients sum to \(E_{24}\). Hence the baseline (8) lies in
this span. Exhaustive optimization gives

\[
 \boxed{S_{\min}(V_{\rm unif})=325.}                            \tag{9}
\]

An optimal decomposition has rank pairs
\((277,277),(576,1),(1,576)\). The upper and lower integer bounds
coincide, so (9) is exact.

The same search was performed on the following orthonormal label
restrictions of the uniform span:

| Label dimension | Restriction | Exact optimal rank sum |
| --- | --- | ---: |
| 7 | Coordinate prefix | 36 |
| 8 | Coordinate prefix | 45 |
| 10 | Coordinate prefix | 66 |
| 12 | Coordinate prefix | 91 |
| 16 | Coordinate prefix | 153 |
| 20 | Coordinate prefix | 231 |
| 8 | Left three-cycle coset sums | 45 |
| 8 | Right three-cycle coset sums | 45 |

Every listed optimum equals (8). For each case, the checker proves
\(S_{\min}\ge d^{3/2}>d^{\sqrt2}\) by an integer power comparison.
This excludes (1) within those particular spans.

### Adding original coefficient directions

Two further exhaustive finite checks go beyond the uniform curve:

1. Every space
   \(\operatorname{span}\{E_{24},T_{S,T},T_{S',T'}\}\).
   The minimum rank sum across all of them is exactly 325.
2. Every space
   \(V_{\rm unif}+\operatorname{span}\{T_{S,T}\}\).
   Every one of these spaces has exact optimum 325.

These searches use the physical \(S_4\times S_4\) symmetry, which
preserves both ranks and the crossing. For the first family, make
the first \(S,T\) initial intervals. The second pair is classified
by the number of its elements inside and outside those intervals,
independently on the two sides. This gives 1,225 parameter
representatives before duplicate tensor pairs are removed, and
766 distinct representative pairs afterwards. Thus every pair of
original coefficient pieces is covered. The second family has
25 representatives, indexed by \((|S|,|T|)\).

These finite calculations alone do not exclude decompositions that use
more independent directions of the full 141-dimensional span. The
subsequent single-cycle argument supplies that stronger exclusion.
Other label maps and other crossings remain outside that argument.

## 5. Implementation and reproduction

- [tools/crossing_rank_sum.py](../tools/crossing_rank_sum.py) provides
  exact coefficient-span optimization, integer certificates, and
  four dual label restrictions.
- [tools/check_seed_subspace_recombination.py](../tools/check_seed_subspace_recombination.py)
  checks the shortest-path method against exhaustive subset
  enumeration, checks the crossing identity under independent dual
  embeddings, and reproduces the finite results above.
- [evidence/seed-subspace-recombination.json](../evidence/seed-subspace-recombination.json)
  saves the lower and upper bounds and explicit rank-pair witnesses.

Run

\[
 \texttt{python tools/check\_seed\_subspace\_recombination.py
 --output evidence/seed-subspace-recombination.json}.
\]

The optional \(\texttt{--extended}\) flag also searches the uniform
span with every pair of original coefficient directions added. Those
larger searches are outside the default finite certificate summarized
above.
