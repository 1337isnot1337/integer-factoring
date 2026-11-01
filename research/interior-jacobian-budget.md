# Binary interior Jacobians extend the finite crossing cutoff to 96 labels

## Result and scope

Over every field of characteristic two, every exact decomposition of the
96-label crossing satisfies

\[
E_{96}=\sum_i T_i,
\qquad \sum_i\sqrt{r(T_i)c(T_i)}>96^{\sqrt2}.
\]

Here the row and column cuts are `(x,y)|(u,v)` and `(x,u)|(y,v)`. Combined
with the earlier all-field finite-seed exclusion through 95 labels, this
extends the characteristic-two cutoff to **96 labels**. The achieved
factoring exponent remains `1/2`. This is a lower bound for one finite
rank-sum seed criterion, not a factoring lower bound.

The physical size is that of the completed crossing decomposition. A
polynomial source on a smaller label space may yield a useful larger
crossing after taking tensor powers and grouping coefficients. This result
does not exclude such asymptotic source constructions.

The previous all-field result and the integer geometric compression
certificates are in PR27 at the pinned source revision
[`b49b9e1`](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/geometric-cut-rank-budget.md)
and its [subspace-compression proof](https://github.com/1337isnot1337/integer-factoring/blob/b49b9e19480c65d139ef15beecacb92064c1a960/research/geometric-subspace-compression.md).
This note adds ten binary full-rank Jacobians. Their minors are odd; their
integer determinants are not asserted to be units. The new conclusion is
therefore restricted to characteristic two.

## 1. The certificate being strengthened

For `m` matrices `P_i` on a `d`-dimensional space, seek subspaces `U,W` with

\[
\dim U=s,\qquad \operatorname{codim}W=q,
\qquad P_i U\subseteq W\quad(1\le i\le m).
\]

Set `p=d-s` and `ell=d-q`. At a standard pair of subspaces, write

\[
P_i=\begin{pmatrix}A_i&B_i\\0&D_i\end{pmatrix},
\qquad A_i\in\operatorname{Mat}_{\ell\times s},
\quad D_i\in\operatorname{Mat}_{q\times p}.
\]

The differential for moving the two subspaces is

\[
\mathcal J_{A,D}:\operatorname{Mat}_{p\times s}
\oplus\operatorname{Mat}_{q\times\ell}
\longrightarrow\bigoplus_{i=1}^m\operatorname{Mat}_{q\times s},
\qquad (X,Y)\longmapsto(D_iX-YA_i)_i.
\tag{1}
\]

A full row rank instance of (1) over `F_2` proves that every tuple over an
algebraically closed field of characteristic two admits the required
subspaces. For completeness, the moving graph equations are

\[
C_i+D_iX-YA_i-YB_iX=0.
\]

Solving for `C_i` gives an affine chart of the incidence variety. At
`X=Y=0`, the projection to matrix-tuple coordinates has full differential
rank exactly when (1) does. Its coordinate functions are then
algebraically independent: a relation of minimum degree would have
vanishing partial derivatives, and over the perfect coefficient field a
polynomial with all partial derivatives zero is a power of a lower-degree
polynomial. Thus the projection is dominant. The full incidence variety
is closed over a product of Grassmannians, which is projective, so its
image is closed. A dense closed image is the whole matrix-tuple space.
This is the same chart and properness argument as the preceding
subspace-compression theorem.

The resulting subspaces may require scalar extension. Matrix ranks do
not change under extension, so this is sufficient for rank inequalities
over the original field.

If `B=sum_i P_i tensor Q_i`, restriction to `U tensor V` followed by
projection of the first output factor to `V/W` annihilates `B`. The swap
crossing induces a map of rank `sq`. Consequently,

\[
\operatorname{rank}(E_d-B)\ge sq.
\tag{2}
\]

The same geometric compression survives a Laurent-series degeneration by
taking saturated subspace lattices, as in the preceding proof. This note
uses the exact-crossing consequence only; it does not update the recorded
formal-limit cutoff.

## 2. Ten explicit interior witnesses

The existing integer allocation certificates optimize over integers
`2<=A<=2m-2`, `B=2m-A`, with

\[
2s+Aq\le2d,\qquad Bs+2q\le2d.
\tag{3}
\]

The full dimension condition for (1) is less restrictive:

\[
msq\le s(d-s)+q(d-q).
\tag{4}
\]

Equation (4) alone is not used as an existence theorem. Each added pair
below has an explicit, independently replayed, full-rank binary
Jacobian. The evidence contains every `A_i` column and `D_i` row as a
hexadecimal binary mask. It also lists the omitted columns of a
nonsingular square minor; all remaining columns form that minor.

| Maps `m` | `s` | `q` | Jacobian rows | Jacobian columns | Old bound | New bound `sq` |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 34 | 43 | 4,386 | 4,387 | 1,444 | 1,462 |
| 5 | 25 | 30 | 3,750 | 3,755 | 729 | 750 |
| 7 | 12 | 38 | 3,192 | 3,212 | 451 | 456 |
| 8 | 10 | 38 | 3,040 | 3,064 | 378 | 380 |
| 9 | 8 | 41 | 2,952 | 2,959 | 320 | 328 |
| 10 | 7 | 41 | 2,870 | 2,878 | 280 | 287 |
| 12 | 6 | 38 | 2,736 | 2,744 | 225 | 228 |
| 14 | 4 | 47 | 2,632 | 2,671 | 184 | 188 |
| 17 | 3 | 50 | 2,550 | 2,579 | 147 | 150 |
| 18 | 3 | 47 | 2,538 | 2,582 | 138 | 141 |

No random success probability enters the theorem. Fixed seeds reproduce
the literal witnesses, and exact elimination checks their full ranks.
The checker also evaluates (1) directly on three independent inputs per
witness to verify the packed coordinate convention.

## 3. Replaying the exact weighted rank budget

Let `F(m)` be the previous exact complement table at `d=96`, strengthened
by the ten entries above. Thus `F(0)=d^2`, `F(1)=F(2)=d(d-1)/2`, and for
`m>=3` it is the maximum of the previous geometric bound and the new
certificates. Zero padding permits the monotone envelope. The table is
zero from `m=2d-1` onward.

For every subset `J` of crossing pieces, (2) implies

\[
\sum_{i\notin J}r_i\ge F\left(\sum_{i\in J}c_i\right),
\tag{5}
\]

and the same inequality holds after exchanging cuts. Define

\[
G(s;b)=\min_{b\le m\le\max(b,2d-1)}\{F(m)+s^2m\}.
\]

Choose a piece whose rank ratio or reciprocal is largest, exchange cuts
so its ranks are `(a,b)` with `a>=b`, and let `k=ceil(sqrt(a/b))`. Its
class satisfies

\[
a_0=(k-1)^2b+1\le a\le a_1=\min(k^2b,d^2),
\quad 2\le k\le d,
\quad 1\le b\le\left\lfloor\frac{d^2-1}{(k-1)^2}\right\rfloor.
\]

The equal-rank class is settled by `sum sqrt(r_i c_i)>=d^2`.
Set

\[
R=\max(d^2-a_1,F(b)),\qquad C=d^2-b.
\]

When `b=1`, also impose `R>=d^2-k(k+1)/2`, using the existing stable-dimension
swap/product theorem. The extra ordinary rank budget is

\[
e=\begin{cases}(k^2-3k+4)/2,&b=1,\\
\max(0,a_0+F(b)-d^2),&b\ge2.
\end{cases}
\]

For `1<=s<k`, put

\[
\mu_s=\frac{s^2k^2-1}{k^4-1},\qquad
\nu_s=\frac{k^2-s^2}{k^4-1},
\]

\[
U_s=\max\{G(s;b),s^2b+\mu_sR+\nu_sC\},\qquad
V_s=\max\{G(s;0),b+\nu_sR+\mu_sC\}.
\]

The preceding weighted-cut proof gives the entirely rational bound

\[
S\ge J(d,k,b)=\frac{2d^2+e}{2k-1}
 +\sum_{s=1}^{k-1}w_s(U_s+V_s),
\quad
w_1=\frac16,\quad
w_s=\frac2{(2s-1)(2s+1)}\ (s\ge2).
\tag{6}
\]

The new checker covers every one of the **15,014** scalar classes at
`d=96`. In 14,471 classes the nonnegative base term alone already exceeds
the running candidate minimum; the other 543 receive a full rational
evaluation of (6). The exact minimum is attained at `(k,b)=(52,3)`:

\[
\min J(96,k,b)=
\frac{405834227274447313387}{637288486874794755}
=636.813994968937\ldots.
\]

The strict comparison is certified by integers:

\[
405834227274447313387^{408}
>
96^{577}\,637288486874794755^{408},
\qquad 577^2-2\cdot408^2=1.
\]

Since `577/408>sqrt(2)`, this proves the stated target inequality.
For an independent limit control, the old table at the same class gives

\[
\frac{738176290889316273953}{1163304380803196775}
<96^{1393/985}<96^{\sqrt2},
\]

also checked by exact integers. The improvement therefore comes from
the added complement certificates.

## 4. Reproduction and limits

```sh
python tools/check_interior_jacobian_budget.py \
  --output evidence/interior-jacobian-budget.json
```

The standard-library replay takes about six seconds in the recorded
workspace. It reconstructs all literal binary matrices, verifies their
Jacobians, recomputes the strengthened complement table, checks every
scalar class, and proves the final power comparison without floating
point.

The new minors are characteristic-two certificates. The old all-field
cutoff remains 95. No stronger factoring exponent, realizable seed, or
improved formal-limit cutoff is asserted. A separate exploratory scan of
the dimension-feasible table at 97 leaves three below-target ratio
classes; dimension feasibility alone is not an additional theorem.
