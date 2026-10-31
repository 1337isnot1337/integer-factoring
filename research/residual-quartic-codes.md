# Polynomial codes for a singular quartic endpoint

This note constructs an explicit polynomial crossing after restricting a rational,
non-Gaussian tensor to a smaller physical alphabet. The construction has a typed
crossing endpoint and a quartic correction that survives the type restriction.
It is a candidate for pressure optimization, not a certificate of an improved
factoring exponent.

The rational quartic and its residual endpoint are developed in
[`localized-quartic-profiles.md`](localized-quartic-profiles.md). The general
strict-energy endpoint criterion is in
[`positive-energy-quartic-endpoints.md`](positive-energy-quartic-endpoints.md).
Here the energy inequality is an equality. A small physical code repairs the
resulting endpoint.

## 1. The Gaussian tensor and its singular perturbation

Work over \(\mathbb F_2\). Let \(M\ge2\) be a power of two and put
\(T=M-1\), \(c=1+w\), and \(z=w^2\). Physical labels are subsets of
\(\{0,\ldots,T\}\). Write

\[
 P_{ij}=[j\subseteq i],\qquad Q=\mathrm{diag}(c^{2i}),\qquad
 H_{ij}=w^{2(T-(i\mathbin{|}j))}.
\]

The inclusion relation is bitwise. The one-particle matrix, with rows \((V,Y)\)
and columns \((U,X)\), is

\[
 D(w)=\begin{pmatrix}cQP^T&w^{2M}P\\ I&cH\end{pmatrix}.
\]

Let \(G=\Gamma(D)\) be its Gaussian exterior tensor. Thus the coefficient of
\(x_Xy_Yu_Uv_V\) is the minor of \(D\) with rows \(V\cup Y\) and columns
\(U\cup X\). In characteristic two, the exterior generators commute and square
to zero. Consider the rational perturbation

\[
 \mathcal T(w)=G(w)\left(1+\frac{x_0y_Tu_0v_T}{c^2}\right). \tag{1}
\]

The column flattening of (1) has a generic rank loss, but (1) has a pole at the
desired endpoint \(w=1\). The rest of the construction resolves the pole and
repairs the endpoint simultaneously.

## 2. A two-letter code on the extreme modes

Order the four extreme-mode states as

\[
 (\varnothing,\{0\},\{T\},\{0,T\}).
\]

Keep the bulk modes \(1,\ldots,T-1\) unchanged. On each input leg \(x,u\), use
the polynomial embedding

\[
 B(w)e_0=|\varnothing\rangle,\qquad
 B(w)e_1=c|\{0\}\rangle+|\{0,T\}\rangle. \tag{2}
\]

On each output leg \(y,v\), use the constant projection

\[
 Ce_\varnothing=e_0,\qquad Ce_{\{T\}}=Ce_{\{0,T\}}=e_1,
 \qquad Ce_{\{0\}}=0. \tag{3}
\]

Equivalently, the two rows of \(C\) are
\(\langle\varnothing|\) and \(\langle\{T\}|+\langle\{0,T\}|\).
Apply these maps independently on all four physical legs of (1). Call the
result \(K_M(w)\). Its logical alphabet has size

\[
 q=2^{M-1}.
\]

The low \(M-2\) bits of a logical label encode bulk occupation. Its highest bit
encodes \(e_0\) or \(e_1\) in (2) and (3).

### Proposition 1: the coded tensor is polynomial

Every entry of \(K_M(w)\) belongs to \(\mathbb F_2[w]\).

**Proof.** Only the quartic part requires a denominator cancellation. A nonzero
quartic term forces the logical core bit to be one on all four legs. For each
input leg there are two choices in (2).

If the low-only state is chosen, its embedding coefficient supplies a factor
\(c\). If the both-occupied state is chosen, removing the low input variable
consumed by the quartic leaves a top input column. Both top output rows have
already been consumed by the quartic. The remaining \(u_T\) column is divisible
by \(c\): its only column-block identity entry lies in the removed \(y_T\) row,
and its other entries are \(c^{2i+1}\). The remaining \(x_T\) column is likewise
divisible by \(c\): its only Pascal entry lies in the removed \(v_T\) row,
and all its remaining entries equal \(c\).

Thus each of the two input legs supplies one factor of \(c\), independently of
which state is chosen. Every determinant term cancels the denominator \(c^2\)
in (1). This is a termwise cancellation, so no cancellation between different
physical code choices is needed. \(\square\)

## 3. A triangular endpoint gives a crossing on common types

For a logical label \(a\), define its integer height

\[
 h(a)=\sum_{i=1}^{T-1}i\,[i\text{ is occupied in }a]
       +T\,[\text{the logical core bit of }a\text{ is one}]. \tag{4}
\]

### Proposition 2: endpoint support

If \(K_M(1)_{xyuv}\ne0\), then

\[
 h(y)+h(v)-h(x)-h(u)\ge0. \tag{5}
\]

When equality holds in (5), the entry is exactly
\([x=v][y=u]\). Consequently, for every tensor power, restriction of all four
legs to one common logical word type has the crossing endpoint. The same is
true for any common subspace spanned by words with one fixed total height.

**Proof.** Introduce a temporary variable \(t\), set
\(w=1+t^T\), and multiply each physical tensor entry by
\(t^{\text{output index sum}-\text{input index sum}}\). This makes the
one-particle entries polynomial:

\[
 \begin{aligned}
 A_{ij}&=t^{T(2i+1)+i-j}P^T_{ij},\\
 (VX)_{ij}&=w^{2M}t^{i-j}P_{ij},\\
 (YU)_{ij}&=[i=j],\\
 B_{ij}&=t^{T+i-j}w^{2(T-(i\mathbin{|}j))}.
 \end{aligned} \tag{6}
\]

All exponents are nonnegative on the indicated support. The quartic in (1)
has index increase \(2T\), so its coefficient after this gauge is exactly one.
At \(t=0\), (6) is the crossing matrix together with only the two internal
edges \(u_Tv_0\) and \(x_Ty_0\). The complete extreme-mode endpoint is therefore

\[
 R=\Gamma(\text{crossing})
       (1+u_Tv_0)(1+x_Ty_0)(1+x_0y_Tu_0v_T). \tag{7}
\]

On this gauged tensor use the constant input columns
\((|\varnothing\rangle,|\{0\}\rangle+|\{0,T\}\rangle)\)
and the output rows (3). Direct expansion of the seven factors in (7)
gives the identity \(4\times4\) transfer on the two logical core bits. The bulk
endpoint is already a crossing. This sixteen-entry calculation is reproduced
exactly by the checker.

Let \(\widetilde K_M(t)\) denote this coded, gauged tensor. For each logical
entry, (2) and the constant gauged code give the exact relation

\[
 \widetilde K_M(t)_{xyuv}
 =t^{h(y)+h(v)-h(x)-h(u)}K_M(1+t^T)_{xyuv}. \tag{8}
\]

The left side is polynomial and has the crossing endpoint. Proposition 1
shows that the right-side tensor is polynomial before substitution. If the
height difference is negative, polynomiality of (8) forces \(K_M(1)_{xyuv}=0\).
If it is zero, (8) identifies the endpoint entry with the crossing. This proves
(5) and its equality statement.

In a product of entries between words of the same total height, the local
height differences sum to zero. Every nonzero local endpoint entry has
nonnegative difference, so every difference must be zero. The product is
therefore precisely the word crossing. \(\square\)

The endpoint need not be a crossing before type restriction. For \(M=4\),
there are exactly 64 additional local endpoint entries, each of positive height.
Rejecting this seed on that basis would discard a valid typed crossing.

## 4. The quartic survives the type restriction

For \(M=4\), the logical alphabet has eight letters and the polynomial degree
is 32. Its quartic contribution is nonzero in 186 entries. In particular it is
nonzero on each of the diagonal quadruples

\[
 (x,y,u,v)=(5,5,5,5),\ (6,6,6,6),\ (7,7,7,7).
\]

These entries have zero charge for every physical word type. Thus the
perturbation is not removed by the one-sided charge-pruning theorem.
This statement only establishes that the candidate is distinct after the
restriction; it does not establish a pressure improvement.

The superficially similar constant code with input columns
\((|\{0\}\rangle,|\varnothing\rangle+|\{0,T\}\rangle)\)
and output rows \((\langle\{0\}|,\langle\varnothing|)\) behaves differently.
Its correction has strictly one-sided occupation charge and disappears on a
common word type. That completion reproduces the truncated Gaussian family.

## 5. Finite classification of the residual endpoint

The untwisted operator (7), with rows \((y,v)\) and columns \((u,x)\), has
rank 12 on a sixteen-dimensional space over \(\mathbb F_2\).

There are 35 two-dimensional subspaces of \(\mathbb F_2^4\). Testing every
subspace with every left-inverse output map gives 560 same-code choices and
exactly thirteen crossings. The constant version of (2) and (3) is one of
them. This thirteen-code count is limited to the specified same-code,
left-inverse normalization.

For three-dimensional alphabets the checker makes a stronger, complete binary
test. It chooses the two input subspaces and the two output rowspaces
independently. There are fifteen possibilities for each, hence \(15^4=50{,}625\)
four-subspace choices. For each choice it forms the restricted \(9\times9\)
transfer. An exact realignment test asks whether this transfer is
\(A\otimes B\) with both \(3\times3\) factors invertible. Such a product can be
normalized to a crossing by changing the two output bases. Conversely, any
three-letter restriction to a crossing must pass this test, regardless of
which bases were used to enumerate the subspaces.

No choice passes. Therefore (7) has no three-letter crossing restriction by
four independent constant physical maps over \(\mathbb F_2\). This finite
statement does not exclude extension-field maps, polynomial maps, larger-block
codes, or approximate restrictions.

## Reproduction

```sh
python tools/check_residual_quartic_codes.py \
  --output evidence/residual-quartic-codes.json \
  --seed-output evidence/direct-residual-code6-n4.json
```

The checker verifies the residual classification, termwise polynomial
division, and all local endpoint entries for \(M=2,4\). An independent expansion
of the quadratic exterior factors checks every exported coefficient against
the determinant constructor. It also checks every two-copy fixed-mode type.
The compact eight-letter seed stores polynomial bit masks at array index
\(x+qy+q^2u+q^3v\), with endpoint \(w=1\), for direct use by pressure searches.
