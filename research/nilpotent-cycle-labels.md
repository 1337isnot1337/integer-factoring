# Self-dual nilpotent cycle labels

## Result and scope

This note constructs a uniform family of self-dual label subspaces for
the twenty-piece three-label crossing from the
[pinned PR #24 certificate](https://github.com/1337isnot1337/integer-factoring/tree/e96ee50560284301428d72421330d8773155d25b).
The construction uses linear combinations of word labels, so it is
different from a coordinate restriction of the tensor power.

The complete family has

\[
D_t=\frac{3^t-(-1)^t}{4}
\]

labels. Its Gram matrix has an explicit inverse, its coefficient
pieces are effective in \(2^{O(t)}\) time, and their number is polynomial
in \(t\). Nevertheless, at \(t=3m+1\) one literal piece satisfies

\[
r\ge \frac14\frac{(3m)!}{(m!)^3}8^m,
\qquad
c\ge \frac{(3m)!}{(m!)^3}8^m.                 \tag{1}
\]

Consequently the complete growing label family has joint-width
coefficient at least

\[
\rho\ge\frac{\log 6}{2\log 3}=0.81546\ldots.
\]

This is a bound for the full label family and for restrictions retaining
the displayed block. It does **not** establish the same bound for an
arbitrary binary prefix after normalization, or for arbitrary
recombinations of the coefficient pieces.

## 1. The local basis and differential

Use the exact matrices \(A_0,\ldots,A_{19}\) recorded in
[finite-seed-amplification.md](finite-seed-amplification.md), with row
indices \((x,y)\), column indices \((u,v)\), and crossing
\(E(x,y,u,v)=\delta_{xv}\delta_{yu}\).
Over \(\mathbb F_2\), change basis by

\[
S=(e_0+e_1,\ e_2,\ e_0),\qquad
A_i'=(S^T\otimes S^T)A_i(S^{-T}\otimes S^{-T}).       \tag{2}
\]

The transformed pieces still sum to \(E\). Call the new coordinate
letters \(v,w,q\), in that order. The active letters are \(v,q\);
\(w\) is passive. Define

\[
Nq=v,\qquad Nv=Nw=0,
\qquad d_t=\sum_{j=1}^t N_j.
\]

Here \(N_j\) acts on the \(j\)-th tensor factor. Since the field has
characteristic two and \(N^2=0\), we have \(d_t^2=0\).

Fix an active-position set \(K\) of odd cardinality, and let \(j\) be
its first position. On all positions outside \(K\), put \(w\). For each
word \(a\) in \(v,q\) on \(K\setminus\{j\}\), define

\[
Z_K(a)=v_j\otimes a+q_j\otimes d_{K\setminus\{j\}}a. \tag{3}
\]

The passive factors are suppressed in (3). These vectors are a basis
of \(\ker d_t\) within this active sector: a vector
\(v_j\otimes a+q_j\otimes b\) lies in the kernel precisely when
\(b=da\). Distinct active sectors have disjoint supports.
Take the direct sum over all odd \(|K|\). Its dimension is

\[
\sum_{k\text{ odd}}\binom tk2^{k-1}
=\frac{3^t-(-1)^t}{4}.
\]

## 2. An explicit nondegenerate pairing

Use the ordinary coordinate pairing in the new basis. In one sector,
put \(D=d_{K\setminus\{j\}}\). The Gram matrix of (3) is

\[
G_K=I+D^TD.                                         \tag{4}
\]

There are an even number of active positions in \(D\). On one active
position, \(NN^T+N^TN=I\); operators on different positions commute.
Thus

\[
DD^T+D^TD=(|K|-1)I=0.
\]

Together with \(D^2=(D^T)^2=0\), this gives
\((D^TD)^2=0\), and therefore

\[
G_K^2=I.                                            \tag{5}
\]

Let \(Z\) be the matrix of all columns (3), and let
\(G=Z^TZ\), the direct sum of the matrices (4). Contracting all four
legs of a tensor coefficient with \(Z\) gives a raw crossing sum
\(G_{xv}G_{yu}\). Replacing the columns on the \(u,v\) legs by \(ZG\)
normalizes it to \(\delta_{xv}\delta_{yu}\), because \(G^2=I\).
This normalization is an invertible change on two external legs, so
it preserves both flattening ranks of every complete piece.

## 3. Uniform construction costs

Each vector (3) has at most \(t\) nonzero word coefficients. The dual
vectors \(ZG\) are also polynomially sparse: (4) applies at most two
one-position differentials, followed by (3). All these coefficients
are explicit bits.

The coefficient pieces come from

\[
\left(\sum_{i=0}^{19}X_iA_i'\right)^{\otimes t}.
\]

There are at most \(\binom{t+19}{19}\) piece-count vectors. An entry
can be computed by expanding the polynomially many sparse word terms
and using a coefficient dynamic program on the twenty counters. Its
time is polynomial in \(t\), with a fixed exponent. Full tensor tables
and exact Gaussian rank factorizations cost \(2^{O(t)}\).

For binary length \(\lambda\), choose the least \(t\) with
\(D_t\ge2^\lambda\), and retain any matching subset of normalized
primal and dual labels. Then \(t=\Theta(\lambda)\), the crossing is
still exact, and the complete construction costs \(2^{O(\lambda)}\).
This cost assertion does not supply a favorable rank bound for those
binary restrictions.

## 4. A large block survives the differential

Set \(t=3m+1\). Within the cycle-label basis, retain the vectors

\[
Z(a)=v\otimes a+q\otimes da,
\]

where \(a\) is a word of length \(3m\) containing exactly \(m\)
passive letters \(w\). It has \(2m\) active letters, so the full vector
belongs to an odd active sector. The operator \(d\) preserves this
word space.

Let \(B\) be the coefficient tensor of type
\(4^m5^m17^m\) on the remaining \(3m\) positions, restricted to these
word labels. The exact core-block calculation in
[balanced-three-label-barrier.md](balanced-three-label-barrier.md)
gives, in both cuts,

\[
\operatorname{rank}B=R_m:=\frac{(3m)!}{(m!)^3}8^m.    \tag{6}
\]

The active change of basis (2) preserves these ranks and the prescribed
passive-position count.

Consider the full type \(16\,4^m5^m17^m\). On the first position,
all four cycle vectors use only active letters. None of pieces
\(4,5,17\) has an all-active entry. Piece \(16\) is the identity on
the active pair: its nonzero all-active entries are

\[
(v,v,v,v),\ (v,q,v,q),\ (q,v,q,v),\ (q,q,q,q).
\]

The first position must therefore use piece \(16\), and the contracted
tensor is exactly

\[
B'=(1+d_xd_u)(1+d_yd_v)B.                            \tag{7}
\]

The notation means substitution of \(da\) on the indicated external
legs, not multiplication of scalar entries.

In the column cut \((x,u)\mid(y,v)\), the two factors in (7) are
invertible row and column operators. Each squares to the identity.
Thus \(\operatorname{rank}_{\rm col}B'=R_m\).

In the row cut \((x,y)\mid(u,v)\), expanding (7) gives four terms,
each obtained by left and right matrix multiplication. The same
four-term operation is its own inverse. Applying rank subadditivity
to that inverse gives

\[
R_m=\operatorname{rank}_{\rm row}B
\le4\operatorname{rank}_{\rm row}B'.
\]

This proves (1). Since
\(\log_2R_m=3m\log_2 6-O(\log(m+2))\) and
\(\log_2D_{3m+1}=(3m+1)\log_2 3-O(1)\), the stated asymptotic
coefficient follows along this subsequence.

The restriction in this proof is a direct sum of entire active-sector
Gram blocks. Normalizing by \(G\) therefore preserves its ranks too.

## 5. Exact finite checks

[check_nilpotent_cycle_labels.py](../tools/check_nilpotent_cycle_labels.py)
checks the Gram involution through \(t=6\), the complete coefficient
profile at \(t=3\), and the core identity (7) at \(m=1\). The complete
three-position family has seven labels, 205 nonzero types, and maximum
rank product \(144=12\cdot12\). The displayed twelve-label block at
\(t=4\) has rank \(48\) in both cuts, before and after (7).

These finite ranks motivate the construction and verify the stated
identities. The all-size limitation uses the argument in Section 4,
not an extrapolation from those ranks.

```sh
python tools/check_nilpotent_cycle_labels.py \
  --output evidence/nilpotent-cycle-labels.json
```
