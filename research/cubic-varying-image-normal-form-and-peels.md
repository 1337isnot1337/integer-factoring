# Structure and bounded optimization of the varying-image space

The [154-map varying-image construction](cubic-varying-image-affordable-escape.md) passes the affordable-span prerequisite and has a cost-36 witness over F2.
The results here give an exact normal form for its 177-dimensional space, exclude two natural matrix-peeling improvements over all characteristic-two fields, and record a bounded F2/F4 search with no cost below 36.
They do not prove that 36 is optimal in the full space.

## Exact normal form

Let \(F\) be a field of characteristic two. Let \(D_i=e_ie_i^T\), let \(A_{ij}=e_ie_j^T+e_je_i^T\), and let \(I=\sum_iD_i\).
In this note, \(A\otimes B\) denotes the tensor with column flattening \(\operatorname{vec}(A)\operatorname{vec}(B)^T\).
Let \(A_0\) be the baseline raw \(pr^2\) coefficient.

Define the 154-dimensional symmetric tensor space \(Q\) with basis
\[
 D_i\otimes D_j\quad(49\text{ directions}),
\]
\[
 A_{ij}\otimes D_i,\ A_{ij}\otimes D_j,
 \quad D_i\otimes A_{ij},\ D_j\otimes A_{ij}
 \quad(84\text{ directions}),
\]
\[
 A_{ij}\otimes A_{ij}\quad(21\text{ directions}).
\]
The full mixed coefficient space is exactly
\[
 \boxed{\quad
 \mathcal V=Q+\operatorname{Sym}_7(F)\otimes I+\operatorname{span}\{E,A_0\}.
 \quad}
\]
The second term adds precisely the 21 directions \(A_{ij}\otimes I\), and the final two directions are independent modulo the symmetric tensor space.
The checker constructs all 177 basis tensors, proves their independence, and reconstructs all 770 physical raw directions in both flattenings.
It also verifies that those physical directions span all 177 coordinates.
These are F2 matrix identities and therefore persist over every characteristic-two extension.

## First matrix peel: arbitrary alternating perturbations

Let
\[
 K=E+H+N,
\]
the balanced tensor whose two cut ranks are \((22,22)\).
Consider peeling off
\[
 U=(I+A)\otimes I,
\]
where \(A\) is any alternating seven-by-seven matrix: it is symmetric with zero diagonal in characteristic two.
This entire family belongs to the normal form above.

The row flattening of the remainder acts on matrices as
\[
 R(K+U):X\longmapsto X+AX^T+\operatorname{tr}(X)I.
\]
We prove its rank is at least 27 over every characteristic-two field.
Put \(r=\operatorname{rank}A\). Since \(A\) is alternating and its size is odd, \(r\leq6\).
The kernel's trace-zero subspace satisfies \(X=AX^T\).

Write \(A=PJP^T\), where \(P\) has full column rank \(r\) and \(J\) is invertible.
From \(X=AX^T\), first write \(X=PY\). Substitution then shows that the row space of \(Y\) lies in the row space of \(P^T\), so \(X=PZP^T\).
The equation becomes
\[
 Z=MZ^T,\qquad M=J(P^TP).
\]
For any fixed \(M\), the solution space of this equation has dimension at most \(r(r+1)/2\).
Indeed, its upper-triangular entries determine \(Z\): a solution with all those entries zero is strictly lower triangular; its first row is zero, the equation makes its first column zero, and induction removes every row and column.

Allowing a nonzero trace can add at most one kernel dimension. Hence
\[
 \dim\ker R(K+U)\leq r(r+1)/2+1\leq22,
 \qquad \operatorname{rank}R(K+U)\geq27.
\]
The peel has column rank one and row rank \(7\operatorname{rank}(I+A)\geq7\).
The remainder's column rank is at least \(22-1=21\).
Therefore this replacement of \(K\) costs at least
\[
 \sqrt7+\sqrt{27\cdot21}=10\sqrt7>22.
\]
Thus no member of this 21-parameter alternating family improves the balanced tensor's cost.

## Second matrix peel: a single alternating edge

Fix an edge matrix \(A=A_{ij}\) and any field element \(t\).
Consider
\[
 U_t=A\otimes(A+tI).
\]
The exact ranks are

| Parameter | Peel ranks | Remainder \(K+U_t\) ranks |
|---|---|---|
| \(t=0\) | \((4,1)\) | \((25,21)\) |
| \(t=1\) | \((12,1)\) | \((34,21)\) |
| \(t\ne0,1\) | \((14,1)\) | \((35,21)\) |

Here is an all-field derivation of the remainder row ranks.
After permuting labels, write the two edge coordinates first and let
\(J_2=\left(\begin{smallmatrix}0&1\\1&0\end{smallmatrix}\right)\).
For
\[
 X=\begin{pmatrix}P&Q\\R&S\end{pmatrix},\qquad\tau=\operatorname{tr}(X),
\]
the kernel equations imply \(\tau=0\), \(S=S^T\), \(R=Q^T\), and \(tJ_2Q=0\).
The remaining two-by-two equation is
\[
 P+P^T+J_2P^TJ_2+tJ_2P^T=0.
\]
Writing \(P=\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)\), its four entries are
\[
 d+tb,\quad c+td,\quad b+ta,\quad a+tc.
\]
They force \((1+t)^4b=0\).
Thus \(P=0\) unless \(t=1\), when \(P\) is an arbitrary scalar multiple of the all-ones two-by-two matrix.
The symmetric five-by-five matrix \(S\) has trace zero and dimension fourteen.
The block \(Q\) contributes ten free entries only when \(t=0\).
The kernel dimensions are therefore 24, 15, and 14, giving the displayed row ranks.

The column matrix of \(K\) is the sum of the outer squares of 22 independent vectors: the diagonal vector \(\operatorname{vec}(I)\) and the 21 alternating-edge vectors.
In this rank factorization, the update removes the selected edge's diagonal coefficient and may add an entry in the diagonal-vector column.
That column already has its own pivot, so the coefficient matrix has rank exactly 21 for every \(t\). The remainder column rank is therefore 21.
All three replacement costs exceed 22. For example, at zero the cost is \(2+5\sqrt{21}>22\), and in the other cases the remainder alone already costs more than 22.

The checker verifies the operator formulas coefficientwise and checks every edge at \(t=0,1\) and at a nonbinary F4 value as independent rank controls.
The field quantifier follows from the symbolic kernel calculation, not from those finite controls.

## Bounded regrouping search

The candidate list consists of all distinct physical raw directions from the 154 maps, the 21 alternating edge pieces, \(H,N,K\), and zero.
Deduplication leaves 538 tensors.
For each unordered pair, form
\[
 T=aT_i+bT_j.
\]
Compare \(T\) and its remainder against the four targets
\[
 E,\quad E+H,\quad E+N,\quad K,
\]
with fixed costs \(0,7,7,14\), respectively.
The search covers all nonzero coefficients over F2 and F4; pairs with zero supply every one-generator case.
F4 conjugation reduces its nine coefficient pairs to five representatives, because all source tensors and targets are binary.

| Field | Coefficient-pair cases | Exact residual comparisons | Cost below 36 |
|---|---:|---:|---:|
| F2 | 144,453 | 533,823 | 0 |
| F4 | 722,265 | 2,495,874 | 0 |

For rank products \(p,q\) and fixed cost \(c\), the comparison is exact:
\[
 c+\sqrt p+\sqrt q<36
\]
holds precisely when
\[
 D=(36-c)^2-p-q>0,\qquad4pq<D^2.
\]
Column elimination stops early only after the exact admissible product bound has been exceeded.
The archive records native source and input hashes. Python independently verifies 64 full rank controls with two pivot orders.

This covers the stated two-generator templates around the pairing decomposition. It does not search all vectors of the 177-dimensional space, arbitrary three-piece decompositions, or larger field extensions.

## Reproduction

- [Normal-form and peel checker](../tools/check_cubic_balanced_peels.py)
- [Structural certificate](../evidence/cubic-balanced-peels.json)
- [Regrouping explorer](../tools/explore_cubic_pairing_regroupings.py)
- [Native evaluator](../tools/search_cubic_pairing_regroupings.cpp)
- [Finite search archive](../evidence/cubic-pairing-regroupings.json)

```bash
python tools/check_cubic_balanced_peels.py --output /tmp/cubic-balanced-peels.json
c++ -O3 -std=c++20 -Wall -Wextra -Werror tools/search_cubic_pairing_regroupings.cpp -o /tmp/cubic-regroupings
python tools/explore_cubic_pairing_regroupings.py \
  --native /tmp/cubic-regroupings --output /tmp/cubic-pairing-regroupings.json
```

The existing cost-36 witness remains an upper bound. These results leave its optimality and the strict seed question open in the full varying-image space.
