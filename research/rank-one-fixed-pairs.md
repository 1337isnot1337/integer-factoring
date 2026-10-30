# Rank-one fixed pairs of invertible rank-two transpose channels

## Claim and scope

Let `F` be any field, let `q >= 2`, and let

\[
U(X)=P_1X^TQ_1+P_2X^TQ_2
\]

be an invertible linear map on `Mat_q(F)`. Zero or redundant summands are
allowed. Suppose nonzero rank-one matrices `A,B` satisfy

\[
U(A)=A,\qquad \operatorname{Tr}(B U(X))=\operatorname{Tr}(BX)
\quad\text{for every }X,\qquad \operatorname{Tr}(BA)=0.
\]

Then

\[
\boxed{AB=BA=0.}
\]

No positivity, identity preservation, trace preservation, nilpotence, or
algebraic closure is assumed. This is a restriction on the two fixed matrices,
not a bound on the Jordan block size. Section 3 gives regular examples of every
binary physical size. These results do not improve the factoring exponent.

## 1. Proof of the fixed-pair restriction

Write `A=a b^T` and `B=c d^T`. A change of coordinates
`X=C Y D`, with invertible `C,D`, replaces the fixed pair by

\[
A'=C^{-1}AD^{-1},\qquad B'=DBC.
\]

It preserves the transpose-channel form, invertibility, and the assertions
`AB=0`, `BA=0`. It also preserves the two physical matrix ranks. Choose the
coordinates so that `A'=E_{00}`. The condition
`Tr(B'A')=0` says `(d^T a)(b^T c)=0` before this change of coordinates.

If exactly one of these two factors vanishes, rescaling and extending vectors
to bases reduce the pair to one of

\[
(A',B')=(E_{00},E_{01}),\qquad (E_{00},E_{10}).
\]

These reductions take place over the original field. If both factors vanish,
the conclusion already holds; such a pair can be reduced to
`(E_{00},E_{11})`. We therefore need only exclude the first two cases.

Suppose `A=E_{00}` and `B=E_{01}`. Define the columns and rows

\[
p_i=P_i e_0,\quad r_i=e_0^TQ_i,\quad
q_i=Q_i e_0,\quad s_i=e_1^TP_i\qquad(i=1,2).
\]

The two fixed-point conditions are exactly

\[
\sum_i p_i r_i=E_{00},\qquad
\sum_i q_i s_i=E_{10}.
\tag{1}
\]

For the second identity, use
`Tr(B P_i X^T Q_i)=Tr((Q_i B P_i)^T X)`.

The columns `p_1,p_2` must be independent. Otherwise their span, which contains
the range of `E_{00}`, is `F e_0`. This gives
`s_i e_0=e_1^Tp_i=0`, contradicting the second identity in (1) applied to `e_0`.
Similarly, `q_1,q_2` must be independent: dependence would put them in `F e_1`,
so `r_i e_0=e_0^Tq_i=0`, contradicting the first identity in (1).

Independence in (1) now implies

\[
r_i=\alpha_i e_0^T,\qquad s_i=\beta_i e_0^T.
\]

Taking the `(1,0)` entry of the first identity gives
`alpha_1 beta_1+alpha_2 beta_2=0`. On the `q`-dimensional subspace
`{e_0 v^T : v in F^q}`, the channel acts as

\[
U(e_0v^T)=M v e_0^T,\qquad M=\alpha_1P_1+\alpha_2P_2.
\]

Injectivity of `U` forces the square matrix `M` to be invertible. However,

\[
e_1^TM=(\alpha_1\beta_1+\alpha_2\beta_2)e_0^T=0,
\]

a contradiction. Conjugating the entire map by matrix transposition excludes
the pair `(E_{00},E_{10})` as well. Therefore `AB=BA=0`.

### Consequence for a regular nilpotent resolvent

If `S=I-U` is one nilpotent Jordan block of size `q^2`, its nonzero right and
left fixed matrices automatically satisfy `Tr(BA)=0`: the fixed vector belongs
to `im S`, and the fixed functional vanishes on `im S`. Consequently, if both
physical fixed matrices have rank one and `rank Phi(U) <= 2`, their products
must vanish in both orders.

The four independent physical changes of basis that preserve the crossing
endpoint reduce, on this operator description, to a separable similarity
`X -> C X D`, up to a nonzero scalar. Such changes preserve physical matrix
rank. They cannot turn a rank-one fixed matrix into the identity when `q>1`.
Thus identity-and-trace normalization excludes a real class of regular roots.

## 2. A nonpermutation regular witness

Over `F_2`, take `q=4` and

\[
U(X)=X^TR+A X^TB,
\]

with matrices encoded in binary row-major order by

\[
A=522,\qquad B=15319,\qquad R=50673.
\]

For clarity, these mean

\[
A=\begin{pmatrix}0&1&0&1\\0&0&0&0\\0&1&0&0\\0&0&0&0\end{pmatrix},
\quad B=\begin{pmatrix}1&1&1&0\\1&0&1&1\\1&1&0&1\\1&1&0&0\end{pmatrix},
\quad R=\begin{pmatrix}1&0&0&0\\1&1&1&1\\1&0&1&0\\0&0&1&1\end{pmatrix}.
\]

The checker verifies

\[
S^{15}(E_{11})=E_{00},\qquad S^{16}(E_{11})=0,
\quad S=I+U.
\]

These equations force a Jordan chain of length 16 in a 16-dimensional space,
so `S` is regular nilpotent. Its right fixed matrix is `E_{00}` and its left
fixed trace matrix is `E_{11}`. Its dyadic Choi ranks are

\[
\bigl(\operatorname{rank}\Phi(U),\operatorname{rank}\Phi(U^2),
\operatorname{rank}\Phi(U^4),\operatorname{rank}\Phi(U^8)\bigr)
=(2,12,12,16).
\]

The second rank is already too large for the existing first-coefficient gate.
This witness establishes the feasibility of the fixed-pair geometry, without
supplying a better local crossing profile.

## 3. Regular rank-two permutation roots for every binary size

Let `q=2^a`, and identify the matrix units with pairs `(i,j)` in
`(Z/qZ)^2`. Start with

\[
U_0(i,j)=(j,i+1).
\]

This map has `q/2` cycles, all of length `2q`. A cycle is determined by the
unordered pair `{d,q-1-d}` of coordinate differences. Its intersections with
output column zero are the states `(d,0)` and `(q-1-d,0)`.

Put `h=q/2`. Let `P` be the row permutation that swaps

\[
q-1-k\quad\text{with}\quad k+1,
\qquad k=0,\ldots,h-2.
\]

These transpositions are disjoint. Let `G` apply `P` only in output column
zero, and put `U=G U_0`. The swap for index `k` joins the original cycles
labelled `k` and `k+1`. Successive joins form a path through all original
cycles. Each join merges two distinct current cycles, so `U` is one cycle of
length `q^2`.

Let `Q` be the cyclic matrix for which `U_0(X)=X^TQ`. Then

\[
G(X)=X+(P-I)XE_{00},\qquad
U(X)=X^TQ+(P-I)X^TQE_{00}.
\]

Thus `rank Phi(U) <= 2`. It equals two for `q>=4` and one for `q=2`.
Since `q^2` is a power of two, `S=I+U` has one nilpotent Jordan block of size
`q^2`. Both fixed matrices are the all-ones matrix `J`, which has physical
rank one and satisfies `J^2=qJ=0` in characteristic two.

This family proves that the fixed-pair restriction is compatible with a
regular block of the largest possible size. It does not give an invariant
direct summand that would obstruct regularity. Its small exact controls have
second Choi rank `q^2`, so it does not supply the low early ranks needed for
an exponent improvement.

| q | Full cycle / nilpotence index | Choi ranks at powers 1,2,4,... |
|---|---:|---|
| 2 | 4 | 1, 4 |
| 4 | 16 | 2, 16, 16, 16 |
| 8 | 64 | 2, 64, 62, 48, 40, 64 |
| 16 | 256 | 2, 256, 254, 238, 188, 157, 169, 160 |

## 4. Exact verification

Run

```sh
python tools/check_rank_one_fixed_pairs.py --output /tmp/rank-one-fixed-pairs.json
```

The independent standard-library checker enumerates all 65,536 linear maps on
`Mat_2(F_2)`. Of the 20,160 invertible maps, 1,656 have Choi rank at most two.
It checks all 504 orthogonal rank-one fixed pairs among them and verifies both
zero products. The 612 nonorthogonal fixed pairs are counted separately and
are outside the theorem's hypotheses. It also checks the explicit regular
nonpermutation witness and the controlled-cycle family through `q=16`.

The all-field fixed-pair result and the all-binary-size cycle construction are
proved above. Finite enumeration is an independent control, not their proof.
