# Arbitrary identity selections one below a prime power

## Result and physical hypothesis

Let \(N=p^k\), \(n=N-1\), and let \(S\) be any set of
\(r\) distinct actual identity indices in \(\{0,\ldots,n-1\}\),
with \(r+1\le n\). The indices need not be consecutive or at the
ends of the interval. An occupation-preserving kernel satisfying the
existing invertible energy-filtered crossing endpoint condition has
a regular zero reduction to **\(2r+1\) active modes**.

The one-dummy [Weyl construction](fixed-deficit-identity-cores.md)
then uses a matrix of order \(p^{r+1}\), followed by exact division
of its known Gaussian digit factor. The original physical tensor,
endpoint repair, degree bound, infinity matrix, and paired type costs
still use only the \(r\) actual selected indices.

This extends the one-deficit reduction from end selections to every
selection satisfying the stated size condition. It is a constructive
profile theorem, not a factoring-exponent improvement.

## 1. How to choose the retained auxiliary modes

Append the formal column \(n\) and write

\[
J=S\cup\{n\},\qquad R=|J|=r+1.
\tag{1}
\]

Let \(a\) be the smallest index absent from \(J\). It exists
because \(R<N\). Thus \(0,\ldots,a-1\) belong to \(J\),
while every remaining member of \(J\) is greater than \(a\).
Define

\[
L=\{n-j:j\in J,\ j>a\},\qquad
H=\{n-a,\ldots,n-1\},\qquad T=L\cup H.
\tag{2}
\]

All indices in \(L\) are smaller than \(n-a\), so the sets in
(2) are disjoint. There are exactly \(R\) retained rows, and all
are actual physical auxiliary indices below \(n\).

For example, at \(p=3,n=26\):

| Actual selection \(S\) | Initial-prefix length \(a\) | Retained auxiliary indices \(T\) |
|---|---:|---|
| \((7,18)\) | 0 | \((0,8,19)\) |
| \((0,13)\) | 1 | \((0,13,25)\) |
| \((0,1)\) | 2 | \((0,24,25)\) |

## 2. Every required principal minor is weighted unimodular

Use \(\tau=z^{2P}\), \(c=1-\tau\), and

\[
f_j(i;c)=\sum_h\binom jh\binom ih c^h.
\tag{3}
\]

Let \(H_N=(f_j(i;c))_{0\le i,j<N}\), let \(R_N\) reverse
the parent indices, and put \(B=H_NR_N\). The parent prime-power
[reflected-basis theorem](prime-power-joint-identity-filters.md) gives

\[
\operatorname{ord}_\tau B_{ik}\ge\max(0,i-k),
\qquad B_{tt}(0)=(-1)^t.
\tag{4}
\]

For any index subset \(L\), the principal minor \(B_{LL}\)
is upper triangular modulo \(\tau\), with nonzero diagonal.
It is therefore a regular unit. Its conjugate by
\(D_L=\operatorname{diag}_{t\in L}(\tau^t)\) is also regular
and has the same unit determinant. Both inverses are regular, so

\[
\operatorname{ord}_\tau(B_{LL}^{-1})_{\ell t}
\ge\max(0,\ell-t).
\tag{5}
\]

The selected columns outside the initial prefix are, up to ordering,
the columns \(B_{*,L}\). Interpolating them on rows \(L\)
therefore gives \(r_i=B_{i,L}B_{LL}^{-1}\), with

\[
\operatorname{ord}_\tau(r_i)_t
\ge\min_{\ell\in L}\bigl(\max(0,i-\ell)+\max(0,\ell-t)\bigr)
\ge\max(0,i-t).
\tag{6}
\]

In particular, every coefficient of \(r_u\), \(u\in H\),
has positive valuation because \(u>\max L\). Empty initial
prefixes simply omit the high-row step below.

## 3. The initial prefix is absorbed by a regular Schur complement

Let \(I\) denote the initial columns \(0,\ldots,a-1\), and
write \(R_H=(r_u)_{u\in H}\). At \(\tau=0\),

\[
f_j(i;1)=\binom{i+j}{j}.
\tag{7}
\]

The initial minor \(I_H(0)\) has consecutive rows and initial
columns, hence determinant one over the integers. This remains a
unit in every characteristic. The complete retained minor is a
zero unit because its noninitial high block vanishes modulo
\(\tau\), while both diagonal blocks are units.

For a discarded actual row \(i\), the full interpolation satisfies

\[
B_{i,H}(I_H-R_HI_L)=I_i-r_iI_L,
\qquad
B_{i,L}=r_i-B_{i,H}R_H.
\tag{8}
\]

The Schur matrix in (8) is a zero unit, so \(B_{i,H}\) is
regular. Every discarded row lies below every high row. Thus the
correction to each low coefficient has order greater than \(i-t\)
whenever that lower bound is positive, and is regular otherwise.
Combining with (6) proves

\[
\boxed{\operatorname{ord}_\tau B_{i,t}\ge\max(0,i-t)}.
\tag{9}
\]

This is the complete weighted concentration certificate. The
polynomial adjugate construction, Gaussian companion, dummy pair,
and exact profile division are now precisely those of the
[fixed-deficit theorem](fixed-deficit-identity-cores.md), with \(d=1\).
No contiguity assumption on the noninitial selected indices remains.

## 4. A compact contraction for any selection

Put

\[
v=(\tau^j)_{j\in J},\quad s=\sum_{j\in J}\tau^j,
\quad d_0=1-\tau^N,\quad e=d_0-cs.
\tag{10}
\]

The missing-index interpretation is especially useful:

\[
e=c\sum_{j\notin J}\tau^j=\tau^a e_0,
\qquad e_0(0)=1.
\tag{11}
\]

The one-row Gram correction and Sherman--Morrison identity give

\[
\boxed{
C^{-1}=
\frac{e_0I+\tau^{N-a}c\,v\mathbf1^{\mathsf T}}{d_0 e_0}.
}
\tag{12}
\]

Both denominators in (12) are local zero units. Since
\(a<N\) and every \(j\ge0\), both \(C^{-1}\) and its
weighted conjugate are the identity modulo \(\tau\). Formula
(12) is checked against the separate full adjugate construction
in the implementation.

## 5. The non-scalar contraction has one local eigenchannel

Let \(m=\min J\). Form a square matrix \(U\) whose first
column is \(\tau^{-m}v\), and whose remaining columns are
\(e_j-e_m\), using the ordered augmented labels. Then

\[
\det U=\tau^{-m}s,
\tag{13}
\]

which is a zero unit. Both \(U\) and its conjugate by
\(A=\operatorname{diag}_{j\in J}(\tau^j)\) are regular
with unit determinant. In the first column the weighted monomial
orders cancel; the other columns have only a diagonal entry and
an entry in the minimum-index row. Thus both inverse changes are
regular as well.

The first column spans the rank-one correction in (12); every
other column has coordinate sum zero. Consequently

\[
\boxed{
U^{-1}CU=
\operatorname{diag}\left(
\frac{e_0}{e_0+\tau^{N-a}},
1-\tau^N,\ldots,1-\tau^N\right).
}
\tag{14}
\]

This isolates the non-scalar local contraction into one direction.
It does not grant a free global change of the actual occupation
kernel. The basis can mix actual and formal directions. In addition,
\(\det U(1)=R\) in the coefficient field, so it is singular at
the crossing endpoint whenever \(p\mid R\). In particular, the
two-actual-mode qutrit case has \(R=3\), and its local eigenbasis
is endpoint-singular. Even when this determinant is nonzero, physical
endpoint filtration must still be checked before using the basis
as part of a global construction.

## 6. Exact constructor and scope of verification

The constructor is
[`selected_one_deficit_geometry.py`](../tools/selected_one_deficit_geometry.py).
It returns the retained rows, complete polynomial form matrices,
weighted interpolation certificate, compact and adjugate contraction
parameters, and the local weighted eigenbasis.

The independent
[`checker`](../tools/check_selected_one_deficit_geometry.py) verifies
these identities over characteristics two, three, and five. Its
[`saved record`](../evidence/selected-one-deficit-geometry.json)
includes **every one of the 325 two-mode selections at \(n=26\)**,
as well as small complete selection families and selected triples.

```bash
python3 tools/check_selected_one_deficit_geometry.py
```

This package certifies the geometry and local algebra. It does not
claim that every new arbitrary selection has already received a full
original physical-tensor replay or an operator pressure search.
