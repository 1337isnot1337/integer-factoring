# Uniform orthogonal graph codes from a cubic differential

## Result and scope

Let the coefficient field have characteristic two. Put

\[
F=aD_0D_1D_5+bD_0D_2D_4+cD_1D_2D_3,
\qquad T=D_7+F,
\]

where `D_i` deletes occupation bit `i`. Every summand of `T` has
inherited height 15. The operator mixes particle numbers one and three.

There is an explicit, height-homogeneous graph basis of `ker T`.
Its coordinate retraction is coherent for every `a,b,c`. More strongly,

\[
1\notin\mathrm{span}_{\mathbb F_2}\{a,b,c\}
\tag{1}
\]

guarantees a coherent **orthogonal** retraction at every inherited-height
cutoff, with arbitrarily many additional free modes. This is an exact
criterion for all cutoffs to work: if (1) fails, the full finite core is
degenerate.

The construction works already over F4, for example with
`(a,b,c)=(omega,omega,omega)`. It also works over F16 with the binary
basis coefficients `(2,4,8)` used in the implementation. These are
concrete mixed-occupation candidates that change states below height 16.

The first completed 51-letter profiles have numerical pressure fit
`rho = 0.71143591915`. Their F4 and F16 capped profiles agree exactly.
They do not improve the existing pressure frontier or the factoring
exponent. Larger profiles are separate searches, not consequences of
the Gram theorem.

## 1. Explicit graph coordinates

Write the occupation algebra as a direct sum according to bit 7. For
a vector `u` containing no bit 7, define

\[
P(u)=u+x_7F(u).
\tag{2}
\]

All deletion operators commute and square to zero. Hence `F^2=0` in
characteristic two, and

\[
(D_7+F)(u+x_7F(u))=F(u)+F(u)+x_7F^2(u)=0.
\]

Conversely, write any vector as `u+x_7v`. The equation `T(u+x_7v)=0`
forces `v=F(u)`; its remaining component then vanishes because `F^2=0`.
Thus (2) parametrizes the whole kernel, not just a subspace of it.

Each monomial of `F` lowers inherited height by 15. Multiplication by
`x_7` raises it by 15. Therefore `P` preserves height, even though a
nontrivial term changes occupation number by two. Taking only graph
labels of height at most `H` gives the full kernel within that height
truncation.

Let `L` discard all physical coordinates containing bit 7. Then
`LP=I`. The graph space is stable under every deletion because `T`
commutes with each `D_i`. The row space of `L` is the coordinate
downset without bit 7, so it is deletion closed as well. These are the
coherent-map hypotheses from the existing
[graded quotient construction](coherent-graded-quotients.md).

The height inventory is independent of the coefficients:

\[
A_H(t)=\left[\prod_{j\ge0,\ j\ne7}(1+t^{2j+1})\right]_{\deg\le H}.
\tag{3}
\]

Only finitely many factors contribute to a fixed cutoff. Through
heights 20, 24, 25, and 26, this gives 51, 81, 91, and 101 labels.

## 2. The complete core Gram determinant

First omit the unused mode 6 and every mode above 7. The virtual core
is the 64-dimensional occupation space on modes 0 through 5. Arrange
physical coordinates according to the absence or presence of bit 7.
Then `P` is the vertical block matrix with upper block `I` and lower
block `F`, and its Gram matrix is

\[
G=I+F^TF.
\]

Let the three cubic supports be

\[
A=\{0,1,5\},\qquad B=\{0,2,4\},\qquad C=\{1,2,3\}.
\]

Each pair has intersection size one and union size five. Split the
virtual space by occupation number `k`. The map `F` lowers `k` by
three, so `G` respects this splitting.

For `k<=2`, the Gram block is the identity. At `k=3`, `F` has one
nonzero output row with coefficients `a,b,c`; at `k=6`, it has one
input column with the same three coefficients. The rank-one determinant
identity gives, in both cases,

\[
\det G_k=1+a^2+b^2+c^2=(1+a+b+c)^2.
\]

For `k=4`, no four-element input set contains two of `A,B,C`.
Every input column of `F_4` therefore has at most one nonzero entry.
The six-by-six matrix `F_4F_4^T` is diagonal, with entries

\[
c^2,\ b^2,\ a^2,\ a^2+b^2,\ a^2+c^2,\ b^2+c^2.
\]

Using `det(I+UV)=det(I+VU)` shows that `det G_4` is the product of
one plus these six entries.

For `k=5`, an output has two elements. Such an output cannot avoid
the union of two cubic supports, because that union already has five
of the six available elements. Hence each output row of `F_5` has
at most one nonzero entry. The six-by-six matrix `F_5^TF_5` is diagonal
with the same six entries, in a possibly different order. Its Gram
determinant equals the one for `k=4`.

Combining the four nontrivial occupation blocks gives

\[
\det G_{64}
=\prod_{\varnothing\ne S\subseteq\{a,b,c\}}
  \left(1+\sum_{s\in S}s\right)^4.
\tag{4}
\]

Here subsets refer to the three indexed coefficient positions, even
when two coefficients have equal values. Restoring the free mode 6
makes two orthogonal copies of this core. Thus the full 128-column
core determinant is

\[
\det G_{128}
=\prod_{\varnothing\ne S\subseteq\{a,b,c\}}
  \left(1+\sum_{s\in S}s\right)^8.
\tag{5}
\]

Equations (4) and (5) prove criterion (1).

## 3. Every height cutoff and every free tail

Every graph column is height homogeneous. Distinct height cells have
disjoint physical supports, so their Gram matrices are orthogonal
blocks. When the complete finite core Gram matrix is invertible,
every one of its height-cell blocks is invertible.

Additional modes above 7 are free coordinate factors. Different free
occupation words are orthogonal, and the whole graph space is the
tensor product of the finite core with this coordinate space. A height
cutoff keeps an orthogonal direct sum of complete core height cells,
one for each allowed free occupation word. Each retained cell remains
nondegenerate.

Consequently

\[
L_{\perp}=(P^TP)^{-1}P^T
\]

exists at every cutoff under (1). It preserves height. Its row space
equals the graph column space, so it is deletion closed. Both physical
maps therefore remain coherent.

The usual equality of zero and infinity profiles also applies to
these mixed-occupation orthogonal maps. Indeed, write `G=P^TP` and

\[
B=(P^T\otimes P^T)E(P\otimes P),
\]

where `E` is the common regular compound representative. After the
same constant reordering of their column labels, the two cuts are
respectively

\[
(G^{-1}\otimes G^{-1})B,
\qquad
(G^{-1}\otimes I)B(G^{-T}\otimes I).
\]

Both are equivalent to `B` by constant invertible matrices. Since
`G` preserves inherited height, these equivalences respect the regular
normalizations. No preservation of the exact particle number is needed.

## 4. Exact implementation and finite evidence

The [graph constructor](../tools/tensor_f16_cubic_graph_codes.py) writes
every sparse physical vector, checks `LP=I`, checks each deletion, and
forms the actual Gram inverses. The
[mixed native constructor](../tools/tensor_hot_compounds_f16_mixed.cpp)
contracts terms of equal occupation and assigns zero to pairs of
unequal occupation. Within a parity block, the monomial exponent is
one half the common row-plus-column height. It uses the existing
[exact capped-jet procedure](tensor-f16-capped-search.md).

The [independent checker](../tools/check_tensor_hot_compounds_f16_mixed.py)
compares 16,384 literal polynomial entries, verifies 500 independent
binary Toeplitz-rank identities, and checks 1,536 capped orders against
complete orders. Its mixed orthogonal examples directly compute both
cuts and verify their equality. It also compares (5) with actual
height-block Gram determinants for 131 coefficient triples, including
every triple over the embedded F4 and deterministic F16 controls.
The [retained evidence](../evidence/tensor-hot-compounds-f16-mixed-controls.json)
contains the complete core maps and source hashes.

The two completed height-20 searches use coefficients `(2,2,2)` over
F4 and `(2,4,8)` over F16. Each has 51 labels. Their 2,601 computed
capped orders are identical; the infinity profile follows from the
proved orthogonal identity. The numerical fit to the upper counting
polynomial is `0.711435919154`, and the fit after deleting the capped
tail is `0.711435896349`. These are numerical fits to rigorous
polynomial bounds, not certified global minima.

A [coordinate control](../evidence/tensor-coordinate-cubic-graph-height20-profiles.json)
uses the same 51 graph labels with `(a,b,c)=(0,0,0)`. Its numerical
upper-polynomial fit, `0.711431092142`, is slightly better than the
nonlinear examples. The new orthogonal family therefore supplies a
constructive search route, but these small pilots do not establish
a numerical advantage from the cubic terms.

Reproduce the finite checks and the F4 pilot with:

```sh
python3 tools/check_tensor_hot_compounds_f16_mixed.py
python3 tools/tensor_f16_cubic_graph_profiles.py \
  --coefficients 2 2 2 --height 20 \
  --output evidence/tensor-f4-cubic-graph-height20-profiles.json
```

The completed profile receipts are
[F4](../evidence/tensor-f4-cubic-graph-height20-profiles.json) and
[F16](../evidence/tensor-f16-cubic-graph-height20-profiles.json).
The driver also accepts higher cutoffs and explicit removal of top
basis vectors, subject to fresh closure and Gram checks. Such finite
subspace truncations do not follow automatically from the uniform
full-height-cell theorem.

No new factoring exponent or Lean formalization is claimed.
