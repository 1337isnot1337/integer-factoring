# A swap/product rank bound and the finite rank-sum cutoff

## Results and scope

Let `F` be any field, let `V=F^d`, and let `S(v tensor w)=w tensor v`
be the swap on `V tensor V`. For arbitrary endomorphisms `P,Q`,

\[
\boxed{\quad
\operatorname{rank}(S-P\otimes Q)\ge {d\choose2}.
\quad}                                                        \tag{1}
\]

If either factor is singular, the stronger bound is

\[
\operatorname{rank}(S-P\otimes Q)\ge {d+1\choose2}.              \tag{2}
\]

Both bounds are sharp over every field. The proof also gives a
refinement in terms of the invertible Fitting part of `QP`.

These inequalities strengthen the universal restrictions on a finite
crossing seed. Suppose

\[
E_d(x,y,u,v)=[x=v][y=u]=\sum_i A_i(x,y,u,v),
\]

and let `r_i,c_i` be the ranks of `A_i` in the cuts
`(x,y)|(u,v)` and `(x,u)|(y,v)`. The proof and exact certificate below
show that, for every field and every `2<=d<=12`,

\[
\boxed{\qquad \sum_i\sqrt{r_ic_i}>d^{\sqrt2}.\qquad}             \tag{3}
\]

Consequently, a seed passing the strict sufficient criterion in the
[finite rank-sum transfer](finite-seed-amplification.md) needs at
least thirteen labels. This improves the earlier
[six-label exclusion](six-label-rank-sum-obstruction.md).
It does not limit arbitrary amplification methods, general crossing
families, or the complexity of integer factoring.

## 1. A triangular twisted swap

For any endomorphism `R` of an `r`-dimensional space,

\[
\dim\ker\bigl(I-(I\otimes R)S\bigr)\le\frac{r(r+1)}2.            \tag{4}
\]

Extend the field to an algebraic closure; this does not change rank.
Choose a basis in which
`R e_i=lambda_i e_i+sum_{k<i} R_ki e_k`. Filter the tensor square
by the weight `i+j` of `e_i tensor e_j`. The operator in (4)
preserves this increasing filtration. Its associated graded map is
a direct sum of one-dimensional blocks `1-lambda_i` and the
two-dimensional blocks

\[
\begin{pmatrix}1&-\lambda_j\\-\lambda_i&1\end{pmatrix}
\quad (i<j).                                                   \tag{5}
\]

Every block in (5) has rank at least one, including in
characteristic two. The rank of a filtered map is at least the rank
of its associated graded map: independent leading pivot images
remain independent before taking leading terms. There are
`binom(r,2)` blocks in (5), proving (4).

Equivalently, `S-I tensor R` has kernel dimension at most
`r(r+1)/2`. Multiplication by `S` turns it into
`I-(R tensor I)S`, and exchanging the tensor factors gives (4).
No diagonalizability or restriction on the eigenvalues is used.

## 2. The singular case and the Fitting refinement

Suppose `v` lies in the kernel of `S-P tensor Q`. Applying the
kernel equation twice gives

\[
v=(QP\otimes PQ)v.                                             \tag{6}
\]

Let

\[
U=\operatorname{im}(QP)^d,\qquad
W=\operatorname{im}(PQ)^d.
\]

The Fitting decompositions split each space into this invertible
part and a nilpotent part. On any tensor-product summand containing
a nilpotent part, `QP tensor PQ` is nilpotent, so its difference
from the identity is invertible. Equation (6) therefore implies
`v in U tensor W`.

The map `P` sends `U` into `W`, since
`P(QP)^d=(PQ)^d P`. It is injective there: `Pu=0` implies
`QP u=0`, whereas `QP` is invertible on `U`. Similarly, `Q`
injects `W` into `U`. Thus both restrictions are isomorphisms and

\[
\dim U=\dim W=:r.
\]

Identify `W` with `U` by the isomorphism `P:U -> W`.
On the domain `U tensor W` use the identification `I tensor P`,
and on the codomain `W tensor U` use `P tensor I`. Under these
identifications the equation `Sv=(P tensor Q)v` becomes

\[
Sz=(I\otimes R)z,\qquad R=(QP)|_U.
\]

Equation (4) bounds its solution space by `r(r+1)/2`. Hence

\[
\boxed{\quad
\operatorname{rank}(S-P\otimes Q)
 \ge d^2-\frac{r(r+1)}2.
\quad}                                                        \tag{7}
\]

Since `r<=d`, this proves (1). If either `P` or `Q` is singular,
then `r<=d-1`, proving (2).

For `P=Q=I`, the kernel of `S-I` has dimension `d(d+1)/2`:
there is one free coefficient per unordered pair, including the
diagonal. This proves sharpness of (1) in every characteristic.
Taking `P=Q=diag(I_(d-1),0)` gives exactly the symmetric tensors
on the first `d-1` coordinates as the kernel, proving sharpness
of (2).

## 3. Rank-one cuts force additional rank in the complement

Discard zero pieces and put

\[
N=d^2,\qquad D={d\choose2},\qquad
R_{\rm tot}=\sum_i r_i,\quad C_{\rm tot}=\sum_i c_i.
\]

Rank subadditivity gives `R_tot,C_tot>=N`. If a piece has column
rank one, it has the form `P(x,u)Q(y,v)`. Its row matrix is
`P tensor Q`, so its row rank is a product of two integers in
`{1,...,d}`. Its row complement is `S-P tensor Q`. Thus a
piece with rank pair `(a,1)` forces

\[
R_{\rm tot}\ge a+D_a,\qquad
D_a=N-\frac{h_d(a)(h_d(a)+1)}2,\qquad
h_d(a)=\max_{\substack{st=a\\1\le s,t\le d}}\min(s,t).           \tag{8}
\]

Indeed the stable dimension in (7) is at most
`min(rank P,rank Q)<=h_d(a)`. In particular `D_N=D`, while
`D_a>=D+d` for `a<N`; some products give a stronger bound.
The same argument with the cuts exchanged applies to a piece of
ranks `(1,b)`.

Let `a` be the largest row rank of a column-rank-one piece, or
zero if there is no such piece. Let `b` be the largest column
rank of a row-rank-one piece, or zero if there is no such piece.
Both belong to the finite set

\[
\mathcal P_d=\{0\}\cup\{uv:1\le u,v\le d\}.                     \tag{9}
\]

Define the valid total-rank budgets

\[
R=\begin{cases}N,&a=0,\\\max(N,a+D_a),&a>0,\end{cases}
\qquad
C=\begin{cases}N,&b=0,\\\max(N,b+D_b),&b>0.\end{cases}           \tag{10}
\]

No bound on the number of pieces is assumed.

## 4. Two complementary lower bounds for the cost

### A two-variable dual certificate

Take `u,v>=0` satisfying

\[
\begin{aligned}
Nu+2v&\le\sqrt{2N},&
2u+Nv&\le\sqrt{2N},\\
u+v&\le1 &&\text{if }a>0\text{ or }b>0,\\
au+v&\le\sqrt a &&\text{if }a>0,\\
u+bv&\le\sqrt b &&\text{if }b>0.
\end{aligned}                                                 \tag{11}
\]

Then every allowed piece obeys `u r_i+v c_i<=sqrt(r_i c_i)`,
and consequently

\[
\sum_i\sqrt{r_ic_i}\ge Ru+Cv.                                  \tag{12}
\]

Here are all the endpoint checks behind this implication. On the
rectangle `2<=r,c<=N`, the function
`u r+v c-sqrt(rc)` is convex in each variable separately, so its
maximum occurs at a corner. The first two inequalities in (11)
handle `(N,2)` and `(2,N)`. Their sum gives
`u+v<=2sqrt(2N)/(N+2)<=1`, which handles both remaining corners.
On an allowed rank-one interval `(r,1)`, `1<=r<=a`,
convexity reduces the check to `(1,1)` and `(a,1)`.
The other rank-one interval is identical.

The checker obtains a rational feasible pair by replacing each
square root in (11) by a rational lower bound and examining
intersections of the boundary lines and coordinate axes. Optimality
of the computed dual is unnecessary: any feasible pair proves (12).

### Actual maximum pieces and concavity

For any collection of positive rank pairs with total row rank at
least `A>=2` and total column rank at least `B>=2`,

\[
\sum_j\sqrt{s_jt_j}\ge\sqrt{A-1}+\sqrt{B-1}.                    \tag{13}
\]

For `k>=2` terms, decrease the totals to `max(A,k)` and
`max(B,k)`. Concavity in the row ranks permits all row surplus
above one to be placed on one coordinate; repeat for the column
ranks. Opposite surplus coordinates give the smaller value,

\[
\sqrt{\max(A,k)-k+1}+
\sqrt{\max(B,k)-k+1}+k-2.
\]

This expression is nondecreasing in `k>=2`, proving (13).
For one term, `sqrt(AB)>=sqrt(A-1)+sqrt(B-1)` follows by
squaring; the remaining difference is
`(sqrt((A-1)(B-1))-1)^2>=0`.

If `a>0`, a piece of ranks exactly `(a,1)` exists. Removing it
and applying (13) gives the additional bound

\[
\sum_i\sqrt{r_ic_i}
 \ge \sqrt a+\sqrt{R-a-1}+\sqrt{C-2},                           \tag{14}
\]

whenever `R-a>=2` and `C-1>=2`. There is a symmetric bound for
`b>0`. If `a,b>0` and `(a,b)!=(1,1)`, the pieces
`(a,1)` and `(1,b)` are distinct. Removing both gives

\[
\sum_i\sqrt{r_ic_i}
\ge\sqrt a+\sqrt b+\sqrt{R-a-2}+\sqrt{C-b-2},                   \tag{15}
\]

whenever both residual budgets are at least two. The checker tests
these conditions explicitly. It never counts a shared `(1,1)`
piece twice.

### When one or both cuts have no rank-one piece

If `a=b=0`, every piece has both ranks at least two. Apply (13)
to `r_i/2,c_i/2` and multiply the result by two:

\[
\sum_i\sqrt{r_ic_i}\ge2\sqrt{2(N-2)}.                          \tag{16}
\]

If `a=0,b>0`, every column rank is at least two. Remove an
actual piece `(1,b)`. The other pairs have row rank at least one
and column rank at least two. Apply (13) to their ranks
`r_i,c_i/2`, using budgets `N-1,(C-b)/2`. When `C-b>=4`,
this gives

\[
\sum_i\sqrt{r_ic_i}
\ge\sqrt b+\sqrt{C-b-2}+\sqrt{2N-4}.                           \tag{17}
\]

There is a symmetric bound for `a>0,b=0`. The profile `(a,b)=(0,1)`
is impossible, since its `(1,1)` piece would imply `a>=1`.

Finally, if `a=b=1`, every piece is either `(1,1)` or has
both ranks at least two. Let `j>=1` count the former pieces.
When `j<=N-4`, apply the scaled form of (13) to the other pieces
to obtain

\[
j+2\sqrt{2(N-j-2)}.
\]

This is nondecreasing for `1<=j<=N-4`. For `j>N-4`, the
remaining rank budgets are at most four. If there is one remaining
piece its cost is at least `N-j`; if there are at least two,
their cost is at least four; if there are none, `j>=N`.
These cases all give cost at least `N`. Therefore, also covering
`N=4` directly,

\[
\sum_i\sqrt{r_ic_i}\ge1+2\sqrt{2(N-3)}.                        \tag{18}
\]

## 5. Exact finite certificate

For each `2<=d<=12` and every unordered pair
`a,b in P_d`, the checker takes the maximum of a rational dual
bound (12) and the applicable bounds (14)--(18). The impossible
pair `(0,1)` is omitted. Every square
root `sqrt z` is replaced by

\[
\frac{\lfloor 10^7\sqrt z\rfloor}{10^7}
=\frac{\operatorname{isqrt}(10^{14}z)}{10^7}.                  \tag{19}
\]

All subsequent arithmetic uses exact integers and rational numbers.
For the smallest resulting rational bound `L_d`, it verifies

\[
L_d^{408}>d^{577},\qquad 577^2-2\cdot408^2=1.                  \tag{20}
\]

Since `577/408>sqrt2`, equations (19)--(20) prove (3).
The following decimal figures are only readable summaries of the
stored rational certificates:

| Labels `d` | Certified lower bound, rounded down | `d^(sqrt2)`, approximately |
|---:|---:|---:|
| 7 | 19.2166665 | 15.6728909 |
| 8 | 22.1131573 | 18.9305010 |
| 9 | 24.9957986 | 22.3615909 |
| 10 | 27.8683260 | 25.9545535 |
| 11 | 30.7332750 | 29.6996724 |
| 12 | 33.5924152 | 33.5886659 |

The certificate uses only the scalar maxima (9), two-variable
linear inequalities, and forced actual pieces. It does not
enumerate tensor entries or assume that every surviving rank
profile is realizable.

The same checker exhaustively tests (7) for all pairs of
two-by-two matrices over `F_2` and `F_3`, and tests singular,
invertible, and nilpotent examples in larger dimensions over
several finite fields. These are independent finite controls.
The all-field proof is Sections 1--2.

```sh
python tools/check_swap_product_rank_sum.py \
  --output evidence/swap-product-rank-sum.json
```
