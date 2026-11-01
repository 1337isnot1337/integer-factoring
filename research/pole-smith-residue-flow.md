# Exact pole-kernel Smith profiles from a convex flow on residues

## Constructive result

The [stabilization theorem](pole-smith-stabilization.md) separates the
true factor part of a Cauchy pole kernel from its nonresonant Smith
orders. This note gives an exact algorithm for all the orders when the
node sets are supplied explicitly. It uses a binary residue tree and
integer shortest paths; it never constructs the rational matrix.

Assume `N` and every node in the distinct positive sets `A,B` are odd.
Write `r=|A|`, `s=|B|`, `m=r+s`, `q=min(r,s)`, and choose `h>=1` with
`|N-ab|<2^h` throughout the rectangle. For `k>=h`, use

\[
 M_k(a,b)=\frac{2^kw}{N+2^kw-ab},\qquad w^2+w+1=0.
 \tag{1}
\]

The algorithm computes every Smith order of (1). Once the residues have
been constructed at depth `h`, its number of integer arithmetic and
comparison operations is

\[
 O(m+qh).
 \tag{2}
\]

The residue tree has at most `2m` vertices. If all input integers have
at most `b` bits, take `k=h=O(b)`. Residue generation and compression
cost `m poly(b)` bit operations. The flow integers have
`O(log(h+1)+log(m+1))` bits up to a fixed constant factor, so their bit
cost is also charged by (2) times a polynomial in that bit length.

The full bit cost is therefore

\[
 \boxed{m\,\mathrm{poly}(b).}
 \tag{3}
\]

This is polynomial in the lengths of the explicit node lists. Those
lengths can be exponential in the binary interval description. The
algorithm does not provide a compressed-interval evaluator or a score
above `1/2`.

The implementation is [pole_smith_flow.py](../tools/pole_smith_flow.py).
For odd `N`, restricting a factoring rectangle to its odd nodes retains
every factor pair, so the parity assumption does not discard factors.

## 1. A minor is a sum of squared residue imbalances

For equally sized subsets `R` of `A` and `S` of `B`, of size `j`, define
for `1<=ell<=k` and each residue `c modulo 2^ell`

\[
 r_{\ell,c}=\#\{a\in R:a\equiv c\pmod{2^\ell}\},\qquad
 s_{\ell,c}=\#\{b\in S:Nb^{-1}\equiv c\pmod{2^\ell}\}.
 \tag{4}
\]

The inverse exists because every column node is odd. Then

\[
 \boxed{\quad
 v_2\bigl(\det M_k[R,S]\bigr)
   =\frac12\sum_{\ell=1}^k\sum_c
       (r_{\ell,c}-s_{\ell,c})^2.
 \quad}
 \tag{5}
\]

### Proof

Put `z=N+2^k w`. The Cauchy determinant identity gives

\[
 v_2\bigl(\det M_k[R,S]\bigr)
 =jk+v_2(V(R))+v_2(V(S))
       -\sum_{a\in R,b\in S}v_2(z-ab).
 \tag{6}
\]

The factor `z` in the Cauchy numerator is a unit. Every cross difference
has valuation at most `k`: a true factor pair gives exactly `k`, and a
nonzero integer difference has valuation below `h`. Within-set
differences also have valuation at most `k`. For example, fixing any
`b` in the nonempty column set gives

\[
 |a-a'|\le |(a-a')b|
   =|(N-a'b)-(N-ab)|<2^{h+1},
\]

so a nonzero `a-a'` has valuation at most `h<=k`. The column argument is
identical. Empty minors have determinant one and require no argument.

For a nonzero difference of valuation at most `k`, its valuation equals
the number of levels through `k` at which its two nodes share a residue.
Moreover inversion by `b -> N/b` preserves all two-adic distances on odd
nodes, since `N` and both denominators are units. At level `ell<=k`, the
three contributions in (6) are therefore

\[
 \sum_c\binom{r_{\ell,c}}2,
 \qquad \sum_c\binom{s_{\ell,c}}2,
 \qquad \sum_c r_{\ell,c}s_{\ell,c}.
\]

Their combined contribution is

\[
 \frac12\sum_c(r_{\ell,c}-s_{\ell,c})^2-j.
\]

The `-j` terms cancel `jk` in (6), proving (5). The sum of the
imbalances is zero at every level. Its sum of squares is consequently
even, which also explains the integer half in (5).

## 2. The determinantal ideals are convex flow minima

Insert each row residue `a mod 2^k` and each column residue
`N b^(-1) mod 2^k` into a binary tree, reading bits from least to most
significant. A depth-`ell` vertex is one occupied residue class modulo
`2^ell`. Mark a leaf with its row and column multiplicities; distinct
input integers can have the same residue at the final depth.

Choose `j` rows and `j` columns. Send one unit from each chosen row leaf
toward the root, and one unit from the root to each chosen column leaf.
Cancel opposite flows on each edge. The signed flow from a child to its
parent is exactly the selected row count minus the selected column
count below that child. Therefore (5) says that the minor valuation is
one half the sum of the squared edge flows.

Conversely, add a source with an edge to each row leaf, of capacity
equal to its multiplicity, and a sink with the analogous incoming
column edges. Any integral source-to-sink flow of value `j` selects `j`
rows and `j` columns and has precisely the imbalances above. Thus, if
`Gamma_j` is the minimum valuation among size-`j` minors,

\[
 2\Gamma_j
   =\min_{\substack{\text{integral flow}\text{ of value }j}}
       \sum_e f_e^2.
 \tag{7}
\]

The Smith orders are `Gamma_j-Gamma_(j-1)` for `j=1,...,q`, with
`Gamma_0=0`. No minor is enumerated by the flow algorithm.

Suppress a maximal path with no branch and no marked internal vertex.
All its edge flows are equal. A path of length `ell_e` therefore becomes
one edge with cost `ell_e*f_e^2`. There are at most `m` distinct occupied
leaves and at most `m-1` branching vertices. Keeping a possible
degree-one root gives at most `2m` tree vertices. Compression preserves
every optimization problem (7) exactly. The code constructs this
compressed tree directly: sort the residues after reversing their
`h` bits, identify each group's common low-bit prefix, and split at its
first differing bit. It does not allocate the suppressed paths. Sorting,
prefix tests, and modular inverses have the input cost in (3).

## 3. Exact path augmentation computes all cardinalities

Orient every tree edge from child to parent. If its current signed flow
is `f`, the residual costs of increasing and decreasing it by one are

\[
 \ell_e(2f+1),\qquad \ell_e(-2f+1).
 \tag{8}
\]

Source and sink edges have zero cost, with their usual integral residual
capacities. Start with zero flow. For each successive cardinality, find
a cheapest residual source-to-sink path and augment one unit.

Here is the optimality invariant used by the implementation. Maintain
vertex potentials so that every residual edge has nonnegative reduced
cost. The zero potentials work initially. Dijkstra's exact integer
distances update the potentials. Every edge on the chosen shortest path
then has zero reduced cost. After augmentation its new reverse edge
also has zero reduced cost, while its next forward marginal cost
increases by `2*ell_e`. All other reduced costs remain nonnegative.
Source and sink capacity changes obey the same zero-cost reversal rule.
Thus there is no negative residual cycle after any augmentation.

For convex integral edge costs, this certifies optimality at the current
cardinality. To see this directly, compare with another integral flow
of the same cardinality. Its difference decomposes into residual cycles.
Convexity bounds the actual cost change below by the sum of the current
marginal costs along those cycles, and the potential terms telescope.
Every such sum is nonnegative. Hence the current flow minimizes (7).

The implementation retains this binary-heap algorithm as an independent
reference. It needs `O(qm log(m+1))` integer operations. The next section
removes the repeated scan of the tree.

## 4. Only two ancestor paths change after an augmentation

Every simple residual source-to-sink path starts at an available row
leaf, follows the unique tree path, and ends at an available column
leaf. It cannot revisit the source or pass through the sink before its
end. Maintain three values at each tree vertex:

- `R(v)`: the cheapest directed path from an available row leaf below
  `v` to `v`, together with its leaf;
- `B(v)`: the cheapest directed path from `v` to an available column
  leaf below it, together with its leaf;
- `P(v)`: the cheapest row-to-column path wholly within that subtree,
  together with both endpoint leaves.

At a marked leaf these values are zero whenever the required remaining
capacity is positive, and absent otherwise. For a child `u` of `v`, its
candidate contributions to `R(v),B(v)` are

\[
 R(u)+\ell_u(2f_u+1),\qquad
 B(u)+\ell_u(-2f_u+1).
 \tag{9}
\]

Take their minima. For `P(v)`, compare the two child-internal values and
the paths that use a row from one child and a column from the other.
There are at most two children, so this costs a constant number of
integer operations per vertex. Negative directed edge costs cause no
problem: the tree path between two endpoints is unique. The invariant
in the preceding section still certifies global flow optimality.

The pair in `P(root)` is the next augmentation. Decrease the two
endpoint capacities by one. Increase the signed flow on the row path
up to their lowest common ancestor, and decrease it on the column
path up to that ancestor. No other edge changes. Recompute the two
exclusive ancestor paths from bottom to top, followed by their shared
ancestor path to the root. This updates every affected dynamic-programming
value after its affected children.

There are at most `h` edges on a compressed root-to-leaf path. Thus
initialization refreshes at most `2m` vertices, and each augmentation
refreshes at most `2h+1` vertices. The code records this count and checks
the bound. Its total is `O(m+qh)`, proving (2). The corresponding
residual path cost is exactly twice the next Smith order, so the
algorithm writes the entire profile while it augments.

Each edge flow has absolute value at most `q`; its length is at most
`h`. Every directed path cost and every stored minimum has absolute
value `O(hq+h)`. The accumulated energy is at most `O(hq^2)`. Therefore
all flow arithmetic uses `O(log(h+1)+log(m+1))` bits. Combining this
with `h=O(b)`, `q<=m`, and the charged modular inverses and trie
construction proves (3).

### Independent interpretation as minimum-valuation Schur elimination

There is also a direct algebraic check of the path rule. Suppose the
previous pivots selected row set `U` and column set `V`. The exact
Cauchy Schur formula gives a remaining entry's valuation as

\[
 \begin{aligned}
 k-v_2(z-ab)
 &+\sum_{u\in U}v_2(a-u)-\sum_{v\in V}v_2(z-av)\\
 &+\sum_{v\in V}v_2(b-v)-\sum_{u\in U}v_2(z-ub).
 \end{aligned}
 \tag{10}
\]

The row correction in (10) is the sum of `ell_e*f_e` on the row's
ancestor path; the column correction is its negative on the column's
path. The common path cancels. Since every leaf has original depth
`h`, twice (10) is exactly the directed tree-path cost in (8), with
the computation performed at `k=h`.

Thus the selected path identifies a minimum-valuation entry of the
actual Schur complement, even though the algorithm has never formed
that entry. Minimum-valuation Gaussian pivoting uses integral row and
column operations and produces the Smith orders in increasing order.
This gives a second justification of the algorithm through the exact
Cauchy elimination identity, separate from the convex-flow optimality
argument.

## 5. Higher displacements and the remaining compression problem

Once the profile at `h` has been computed, the
[stabilization theorem](pole-smith-stabilization.md) gives every profile
at `k>=h`: leave all zero orders unchanged and add `k-h` to each positive
order. The implementation uses this rule before forming any modulus
`2^k`. Thus a very large binary-described displacement can be handled
with its written output cost instead of constructing a tree of depth
equal to that displacement.

This does not remove dependence on the number of supplied nodes.
Building the reciprocal residue tree still visits every column node.
For long intervals, a residue class's two capacities can be queried
cheaply from endpoints and a modular inverse, but the algorithm above
does not bound how many such classes must be opened. An efficient
implicit-tree optimizer would need a separate proof covering every
visited class and every shortest-path query. That is the remaining
global evaluation question; the explicit algorithm does not assume its
solution.

## 6. Exact checks

Run:

```
python3 tools/check_pole_smith_flow.py --output /tmp/pole-smith-flow.json
cmp /tmp/pole-smith-flow.json evidence/pole-smith-flow.json
```

The checker compares 2,000 complete flow profiles with independent
valuation-minimizing Gaussian elimination over `Q(w)`. For 200 of those
cases, it separately enumerates every pair of row and column subsets
and minimizes the literal residue energies in (5). Five hundred cases
also compare the dynamic-programming algorithm with the independent
binary-heap flow reference.

Six explicit controls use a 201-bit input and dimensions 16, 64, 256,
1,024, 16,384, and 65,536. Up to dimension 1,024, they check the complete
determinant content by literal Vandermonde and cross-difference
valuations. All six independently evaluate the same content by linear
congruence counts in the column progression, with one charged pass over
the explicit row list.

The largest case uses 262,142 tree vertices and 1,572,860 vertex
refreshes. Its suppressed paths would have produced 13,303,562 vertices;
they are counted from their lengths without being built. The receipt
records every full-profile hash, content, and update count. A separate
control requests displacement `k=10^80`, using the same fourteen-node
tree as its input-precision profile. These controls start from known
node lists and do not claim to find unknown 100-bit factors.
