# A strict mixed-cut budget excludes every local crossing through q=12

## Statement

The existing [geometric cut-budget theorem in PR #27](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/geometric-cut-rank-budget.md) already excludes exact finite seeds through 95 labels. This note gives a shorter independent argument for a smaller range; it does not extend that cutoff. A native finite-size exclusion does not by itself exclude asymptotic constructions obtained from tensor powers of a small polynomial source.

Let `F` be any field, `V=F^q`, and let `E` be the swap crossing

\[
E(x,y,u,v)=\delta_{xv}\delta_{yu}.
\]

For a tensor `T`, write `r(T)` for its rank across `(x,y)|(u,v)` and `c(T)`
for its rank across `(x,u)|(y,v)`. Consider any finite decomposition

\[
E=\sum_i T_i,\qquad r_i=r(T_i),\quad c_i=c(T_i).
\]

Zero pieces can be discarded. For every `q>=3`,

\[
\boxed{\ \sum_i\min(r_i,c_i)\ge q+1.\ }
\tag{1}
\]

Together with the two ordinary cut budgets, this implies

\[
\boxed{\ S:=\sum_i\sqrt{r_i c_i}
 \ge q-1+2\sqrt{q^2-q}.\ }
\tag{2}
\]

In particular, no such decomposition at physical size `2<=q<=12` can have
`S<q^sqrt(2)`, the local crossing cost needed to improve the current factoring
exponent through this transfer mechanism. This statement is universal over
fields and decompositions. It does not impose nilpotence, a polynomial
resolvent, a restricted coefficient family, or a characteristic assumption.
It does not exclude larger local crossings or other factoring methods.

The elementary refinement is the strict `q+1` budget. The size-two exception is real:
over `F_2`, write `E=(E-I)+I`; the residual row rank and correction column rank
are both one, so their sum is two.

## 1. A two-piece mixed-cut inequality

Suppose `E=A+B`, with `r(A)=r` and `c(B)=c`. A minimal column-cut
decomposition writes the row matrix of `B` as

\[
B=\sum_{j=1}^c H_j\otimes J_j,
\]

where the square matrices `H_j` are linearly independent. For any nonzero
`xi in V`, put `W_xi=span{H_j xi}`. Restrict the input to `xi tensor V` and
project the first output factor onto `V/W_xi`. The correction `B` vanishes,
whereas the swap induces

\[
\eta\longmapsto[\eta]\otimes\xi.
\]

This map has rank `q-dim W_xi`. Therefore

\[
r\ge q-\dim W_\xi.
\tag{3}
\]

This already proves `r+c>=q`. We now show `r+c>=q+1` when `q>=3`.

### 1.1 At least two column-cut terms

Assume `c>=2`. Matrix ranks and linear independence are unchanged by scalar
extension, so work over an algebraic closure. The span of two independent
square matrices contains a nonzero singular matrix. Indeed, this is immediate
if one is singular; otherwise an eigenvalue of `H_1^(-1) H_2` supplies a
singular combination `H_2-lambda H_1`, which is nonzero by independence.

Choose a nonzero kernel vector `xi` of that combination. Then the vectors
`H_j xi` satisfy a nontrivial linear relation, so `dim W_xi<=c-1`. Equation
(3) gives

\[
r+c\ge q+1.
\]

The inequality descends to the original field because the two ranks do not
change under extension.

### 1.2 One column-cut term

Assume `c=1`, so `B=H tensor J`, with both factors nonzero. If `H` is singular,
take `xi in ker H` in (3), giving `r>=q`. If `J` is singular, restrict to
`V tensor xi` instead; the same argument gives `r>=q`.

Suppose both factors are invertible, and denote the swap matrix by `P`.
Multiplying on the left and right by invertible Kronecker matrices gives

\[
(H^{-1}\otimes I)(P-H\otimes J)(I\otimes H)
 =P-I\otimes(JH).
\tag{4}
\]

Over an algebraic closure, triangularize `JH`. Conjugation by the same change
of basis in both tensor factors leaves `P` unchanged. We are reduced to
`P-I tensor R` with `R` upper triangular.

Select rows `(i,j)` with `i<j` and columns `(l,k)` with `k<l`, using the
matching order `(l,k)=(j,i)`. The swap submatrix is the identity. The
corresponding entry of `I tensor R` is `delta_(i,l) R_(j,k)`. A nonzero
Kronecker delta would force `k<i<j`, so that entry of the upper-triangular
matrix is zero. Thus the selected minor is an identity matrix of size
`binomial(q,2)`, over every characteristic. Consequently,

\[
r\ge\binom q2\ge q\qquad(q\ge3),
\]

and again `r+c>=q+1`.

### 1.3 Zero column-cut terms and grouping

If `c=0`, then `A=E` has row rank `q^2>=q+1`. This completes the two-piece
proof.

For an arbitrary decomposition, group into `A` the pieces with `r_i<=c_i`
and into `B` the remaining pieces. Rank subadditivity and the preceding
inequality imply

\[
\sum_i\min(r_i,c_i)
 =\sum_{r_i\le c_i}r_i+\sum_{r_i>c_i}c_i
 \ge r(A)+c(B)\ge q+1.
\]

This proves (1).

## 2. Converting the budgets to a cost bound

Put `n=q^2`, `x_i=min(r_i,c_i)`, `d_i=|r_i-c_i|`, and `m=sum x_i`.
Both cuts of `E` have rank `n`, so

\[
\sum_i r_i\ge n,\qquad \sum_i c_i\ge n.
\tag{5}
\]

For fixed `d>=0`, define

\[
f_x(d)=\sqrt{x^2+xd}-x
       =\frac d{\sqrt{1+d/x}+1}.
\]

It increases with `x>=1`, so `f_x(d)>=f_1(d)`. The function
`f_1(d)=sqrt(1+d)-1` is increasing and subadditive on nonnegative reals.
For example, its derivative decreases, so
`f_1(a+b)-f_1(a)<=f_1(b)-f_1(0)`.

When `m<=n`, (5) says that the total row-heavy surplus and the total
column-heavy surplus are each at least `n-m`. Applying subadditivity within
the two groups gives

\[
S=m+\sum_i f_{x_i}(d_i)
 \ge m+2\sqrt{n-m+1}-2.
\tag{6}
\]

The last expression increases with `m` on `m<=n`, since its derivative is
`1-1/sqrt(n-m+1)>=0`. Use `m>=q+1` from (1) to obtain (2). If `m>=n`, then
`S>=m>=n`, which is at least the right-hand side of (2) by
`2sqrt(n-q)<=n-q+1`. Thus no range of `m` is omitted.

### Optimality of this budget relaxation

For `q>=3`, the numerical rank pairs

\[
(q^2-q,1),\quad(1,q^2-q),\quad
q-1\text{ copies of }(1,1)
\]

have row and column sums `q^2`, mixed sum `q+1`, and cost exactly the
right-hand side of (2). Each pair separately satisfies the rank-one product
constraint, since `q^2-q=q(q-1)`. Therefore these three budgets, even together
with that product constraint, cannot yield a stronger universal numerical
cost floor. Their simultaneous realization as a decomposition of the
crossing is **not** asserted.

## 3. Exact cutoff certificates

At `q=12`, the lower bound is

\[
11+2\sqrt{132}=33.978250586\ldots,
\]

whereas `12^sqrt(2)=33.588665890...`. No floating-point comparison is needed:

\[
\sqrt2<\frac{17}{12},\qquad
\sqrt{132}>\frac{57}{5},\qquad
169^{12}>5^{12}\,12^{17}.
\]

These imply `11+2sqrt(132)>169/5>12^(17/12)>12^sqrt(2)`.
For `3<=q<=11`, the lower bound exceeds `3q-3`, and the finite integer
comparisons `(3q-3)^12>q^17` certify all the required inequalities. At `q=2`,
the original non-strict mixed budget yields `S>=2sqrt(3)>3>2^(17/12)`.

The resulting universal exclusion is therefore all physical sizes
`2<=q<=12`. At `q=13`, the new lower bound lies below `13^sqrt(2)`;
this elementary argument alone leaves that size open; the prior geometric theorem still excludes it.

## 4. Independent verification

```sh
python tools/check_strict_mixed_cut_gate.py --output /tmp/strict-mixed-cut-gate.json
```

The standard-library checker verifies every one-term opposite-cut correction
`H tensor J` with nonzero binary three-by-three factors: 261,121 pairs,
including all 28,224 pairs of invertible factors. It recomputes the residual
swap rank, checks the genuine size-two exception, verifies all cutoff
inequalities with integers, and checks the numerical relaxation witnesses.
The all-field strict budget is proved above; the finite check is independent
support for its delicate one-term branch.
