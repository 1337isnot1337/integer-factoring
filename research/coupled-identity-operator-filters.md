# Operator kernels on coupled polynomial identity modes

## Result and scope

The bounded reduction for coupled identity-mode filters extends to an
arbitrary matrix of occupation transfers. A matrix kernel may transfer
occupation between selected modes, provided that it preserves total
occupation and equals the identity at the crossing endpoint. The zero
cut still has only \(2r\) active modes, and the infinity cut has an
explicit weighted \(p^r\)-by-\(p^r\) matrix.

This gives a strictly larger family than diagonal joint state filters.
The construction is polynomial after common monomial clearing, and
the degree and both profile shifts are accounted for below. It is a
search construction, with no claimed improvement to the factoring
exponent.

## 1. Definition and endpoint

Work in characteristic \(p\), with the truncated symmetric algebra in
which each physical variable has \(p\)-th power zero. Let
\(E(T)=\sum_{d=0}^{p-1}T^d/d!\). For the Pascal interval
\(0,\ldots,n-1\), write the one-particle matrix as

\[
\mathcal D(w)=\begin{pmatrix}A(w)&B(w)\\I&C(w)\end{pmatrix}.
\]

For selected identity directions \(j_0,\ldots,j_{r-1}\), put

\[
L_s=Y_{j_s}+\sum_iA_{i,j_s}V_i,\qquad
R_s=U_{j_s}+\sum_iC_{j_s,i}X_i,\qquad Q_s=L_sR_s.
\tag{1}
\]

The initial geometry uses \(j_s=s\). The mixed and terminal
geometries use the selections and regular reductions in
[Section 6 of the coupled reduction](coupled-identity-filter-reduction.md).
All formulas in this note apply to any selection for which the
simultaneous regular concentration of the individual forms is proved.

Let \(\mathcal A=\{0,\ldots,p-1\}^r\), and let
\(F(q)=(F_{ab}(q))_{a,b\in\mathcal A}\) have Laurent-polynomial
entries over the field. Impose

\[
F_{ab}=0\quad\text{if }|a|\ne|b|,\qquad F(1)=I.
\tag{2}
\]

Set \(q=-v\), choose a positive odd integer \(P\), and evaluate
the forms in (1) at \(w=v^P\). Define

\[
H_F=\sum_{a,b\in\mathcal A}F_{ab}(q)
       \frac{L_0^{a_0}\cdots L_{r-1}^{a_{r-1}}
             R_0^{b_0}\cdots R_{r-1}^{b_{r-1}}}
            {b_0!\cdots b_{r-1}!},
\qquad
T_F=G(v^P)\prod_{s=0}^{r-1}E(-Q_s)\,H_F.
\tag{3}
\]

The asymmetric factorial convention makes \(F\) itself the selected
coefficient matrix, up to constant column units. If \(F\) is
diagonal with entries \(\alpha_a\), (3) is exactly the previously
defined joint occupation filter.

At \(v=-1\), one has \(q=1\) and \(w=-1\). Equation (2)
gives

\[
H_F\big|_{q=1}=\prod_sE(Q_s),\qquad
T_F\big|_{v=-1}=G(-1).
\tag{4}
\]

The inverse relation \(E(Q_s)E(-Q_s)=1\) holds because
\(Q_s^p=0\), and every factorial below \(p\) is a field unit.
Thus the endpoint equality holds coefficient by coefficient, before
type restrictions. No endpoint repair is required.

### A larger endpoint class with a free physical repair

The identity condition in (2) can be relaxed. Put
\(E=F(1)\), and define the physical selected-mode energy
\(\epsilon(a)=\sum_sj_sa_s\). It suffices that \(E\) is
invertible, preserves total occupation, and satisfies

\[
E_{ab}=0\quad\text{whenever }\epsilon(a)<\epsilon(b).
\tag{4a}
\]

At the endpoint, every entry of \(A(-1)\) vanishes, so the
selected forms \(L_s\) are precisely the physical variables
\(Y_{j_s}\). The endpoint of (3) is therefore the original
Gaussian with \(E\) applied to its selected physical \(Y\)
occupation labels. Apply the constant label map \(E^{-1}\) to
the polynomial tensor. This restores (4) exactly and adds no
polynomial degree.

The repair preserves both computed profiles. At zero a selected
\(Y\) state has height
\(P(2\epsilon(a)+(1-J)|a|)\). Condition (4a) makes the
weighted versions of \(E\) and \(E^{-1}\) regular, with
regular inverses. At infinity that height is \(-nP|a|\), so
occupation preservation makes the weighted repair a constant unit.
Equal-energy blocks may be any invertible matrices. A diagonal
endpoint with all coordinates nonzero is a special case.

An arbitrary invertible endpoint matrix need not satisfy (4a).
Its constant physical repair can have poles in the zero weighted
basis and is not covered by this cost-free extension. The code
calls the proved extension `endpoint='filtered'` and explicitly
checks the energy filtration and invertibility.

## 2. Polynomial clearing and degree

Let \(h\ge0\) clear all negative powers in \(F\), and write
\(D_F=\max_{a,b}\deg(q^hF_{ab})\). Then \(q^hT_F\) is a
polynomial tensor. Since \(q=1\) at the endpoint, this scalar
clearing preserves (4). Its zero profile shifts by \(+2h\),
and its infinity profile shifts by \(-2h\), so the two
partition-function factors have cancelling scalar changes.

A direct degree bound is

\[
\deg_v(q^hT_F)
\le P\bigl(\deg_wG+2nr(p-1)\bigr)+D_F.
\tag{5}
\]

To verify it, first expand the product in (3). In each surviving
term, the exponent of each \(L_s\) and each \(R_s\) is at most
\(p-1\). Indeed, these linear forms have \(p\)-th power zero.
The Pascal formulas give
\(\deg L_s\le j_s+1\) and
\(\deg R_s\le2n-1-j_s\); their sum is at most \(2n\).
This proves (5), including the factors introduced by removing the
selected Gaussian modes. It is not valid to charge only the degree
of \(F\).

The construction retains the full finite alphabet and polynomial
degree. Reducing the matrices used to calculate profiles does not
remove those costs from the factoring construction.

## 3. Exact zero reduction

Use \(v=-z^2\) and the aggregate zero gauge \(J\). Define

\[
a_j=P(2j+1-J),\qquad b_i=P(2i+1+J).
\tag{6}
\]

For the initial geometry, the regular two-sided concentration from
[the coupled reduction](coupled-identity-filter-reduction.md)
acts on the two forms separately. In its notation, the retained
indices are \(t_\ell=n-r+\ell\), and the active forms are

\[
L_{s,0}=z^{a_s}Y_s+c\sum_{\ell=0}^{r-1}(WF)_{\ell s}V_\ell,
\qquad
R_{s,0}=U_s-z^{-a_s}\sum_{\ell=0}^{r-1}F_{\ell s}X_\ell,
\qquad c=1-z^{2P}.
\tag{7}
\]

The matrices \(F,W\) in (7) are the concentration matrices from
that note, not the occupation kernel \(F(q)\) in (2).
The active Gaussian is

\[
G_{0,\mathrm{act}}=
\prod_{s=0}^{r-1}E(z^{a_s}Y_sU_s)
\prod_{\ell=0}^{r-1}E(z^{b_{t_\ell}}X_\ell V_\ell).
\tag{8}
\]

Replace the forms and Gaussian in (3) by (7) and (8). The resulting
coefficient matrix, flattened as \((Y,X)\mid(U,V)\), gives the
complete active zero profile. The other modes are unchanged Gaussian
spectators. Its matrix order is \(p^{2r}\), independently of
\(n\).

This extension requires no new change of basis: the original proof
concentrates each \(L_s\) and each \(R_s\), and therefore
concentrates every monomial \(L^aR^b\), including \(a\ne b\).
In the zero normalization the individual physical forms carry
reciprocal common factors \(z^{PJ}\) and \(z^{-PJ}\) relative
to the displayed forms. Their contribution is
\(z^{PJ(|a|-|b|)}=1\), using the occupation condition in (2).
The same observation applies to the proved mixed and terminal
reductions. Their regular units must be retained in their active
Gaussian and forms, as in the implementation.

Total occupation preservation in (2) splits the active matrix into
charge blocks

\[
\gamma=|Y|-|X|=|U|-|V|.
\tag{9}
\]

For a fixed block, changing \(J=1\) to another gauge adds
\(P(1-J)\gamma\) to every Smith order. Thus one exact profile
calculation serves every aggregate gauge. Gauge changes still alter
the paired physical type partition and must be included there.

## 4. Exact infinity reduction

At infinity, the natural physical heights are

\[
h(V_i)=P(2i-n+2),\quad h(Y_j)=-nP,\quad
h(U_j)=P(2j-n+1),\quad h(X_i)=P(3n-1).
\tag{10}
\]

The elementary block elimination of \(\mathcal D\) isolates the
identity directions in (1). It remains a regular change in these
weighted coordinates. To check this, use

\[
\deg A_{ij}\le i+1,\qquad
\deg C_{ji}\le2n-1-\max(j,i).
\tag{11}
\]

After dividing \(L_s\) by \(z^{-nP}\), every other term has
nonnegative order. After dividing \(R_s\) by
\(z^{P(2j_s-n+1)}\), the same is true. These substitutions are
triangular with unit diagonal, so their inverses are regular as well.
They leave the Schur complement as the unselected Gaussian factor.

Consequently the active infinity matrix is precisely

\[
\mathcal K_\infty(z)_{ab}
 =z^{-nP|a|+\sum_sP(2j_s-n+1)b_s}
       F_{ab}(z^{-2}),
\tag{12}
\]

up to factorial units. Its exact Laurent Smith orders, together with
the unchanged Gaussian spectators, give the full infinity profile.
Because \(|a|=|b|\), the row weight can be moved to the columns:

\[
\mathcal K_\infty(z)_{ab}
 =F_{ab}(z^{-2})z^{\sum_s\kappa_s b_s},
\qquad \kappa_s=P(2j_s+1-2n).
\tag{13}
\]

For a diagonal kernel with \(F_{aa}=q^{-e_a}\), this recovers
the orders \(\sum_s\kappa_sa_s+2e_a\). Off-diagonal kernels
can change the Smith orders through cancellations between entries;
their entrywise minimum orders alone do not determine the profile.

## 5. A concrete additional degree of freedom

In characteristic two with \(r=2\), leave the vacuum and
full-occupation weights equal to one and use the one-occupation block

\[
\begin{pmatrix}1&1+q^{-1}\\1+q^{-1}&1\end{pmatrix}.
\tag{14}
\]

Its endpoint is the identity and its determinant is \(q^{-2}\).
The constant matrix at infinity is singular, giving a finite
infinity-order shift. The independent full-occupation coefficient
remains one, so this is not the Gaussian lift of its one-particle
block. One may also set that coefficient to \(q^{-1}\), still
preserving the endpoint.

Both variants are valid members of the construction. Their zero
profiles change as well; an infinity shift is not, by itself, a
pressure improvement. The initial numerical screens of these two
variants did not improve the fully optimized Gaussian pressure.

## 6. Exact implementation checks

The implementation is
[coupled_identity_operator_filters.py](../tools/coupled_identity_operator_filters.py).
It uses complete Laurent coefficient arrays. Binary matrices use the
native exact finite-field Smith engine with an adaptive retained
window; a missing promised pivot or unsupported complete source
raises an error. The full generic rank follows from the identity
endpoint and the invertible Gaussian endpoint, rather than from a
truncated numerical rank guess.

Run the independent original-tensor checks with

```sh
python3 tools/check_coupled_identity_operator_filters.py --extended \
  --output /tmp/coupled-identity-operator-filters-replay.json
```

The saved [evidence](../evidence/coupled-identity-operator-filters.json)
contains three exact diagonal-specialization comparisons and eighteen
full original cut-profile comparisons. These include three selected
binary modes, a qutrit transfer with unequal per-mode occupations,
ramification three, and a mixed binary example with original matrix
order 256. Three cases use nonidentity energy-filtered endpoints,
including a qutrit diagonal coefficient equal to two. Every replay checks
the coefficientwise repaired endpoint and
the degree bound (5). A coefficient at order 600 tests adaptation
beyond the initial native window, and an oversized source tests
explicit rejection rather than silent truncation.

## 7. An obstruction to a stronger auxiliary synthesis condition

One tempting sufficient condition for retaining the Gaussian zero
cost is to require the active matrix times the inverse diagonal
Gaussian matrix to be regular. That condition cannot coexist with
a strict improvement of the selected infinity ground order.

Let \(\mathbf t=(p-1,\ldots,p-1)\) be the full occupation
state. Set all active \(V,X\) variables to zero in (3). Then
the selected Gaussian factors cancel the factors \(E(-Q_s)\),
and the coefficient of \(Y^{\mathbf t}U^{\mathbf t}\) is,
up to a constant factorial unit,

\[
z^{\sum_sa_{j_s}(p-1)}F_{\mathbf t\mathbf t}(q).
\tag{15}
\]

The Gaussian diagonal entry in that column is the prefactor in
(15). Thus regularity after right division by the full Gaussian
diagonal forces \(\operatorname{ord}_{q=0}
F_{\mathbf t\mathbf t}\ge0\).

Every selected \(\kappa_s\) is negative, so the Gaussian
selected infinity ground is \(\sum_s\kappa_s(p-1)\).
The total-occupation block of \(\mathbf t\) is one dimensional.
Improving this ground order would require
\(\deg_qF_{\mathbf t\mathbf t}<0\), contradicting its
zero regularity and nonzero endpoint value. The argument also
applies after a filtered endpoint repair, whose full-occupation
block is a nonzero scalar.

This excludes the stated auxiliary column-divisibility condition.
It is not an obstruction for all operator kernels, nor a proof that
their zero Smith profiles must deteriorate by a particular amount.
Searches allowing explicitly charged zero losses remain outside
this argument.
