# A Weyl-algebra reduction for coupled identity operators

## Result and scope

The zero profile of an arbitrary occupation-preserving operator on
\(r\) selected identity directions reduces from a matrix of order
\(p^{2r}\) to a matrix of order \(p^r\), followed by explicit
auxiliary digit shifts. The reduction is exact in every prime
characteristic. It applies to the initial, terminal, and mixed
geometries from the
[coupled identity reduction](coupled-identity-filter-reduction.md).

The reduced matrix is a finite differential-operator transform of
the occupation kernel. For one diagonal mode, it recovers the
finite-difference weights \(\beta_k\). It also handles off-diagonal
occupation transfers and correlations among several modes. This is
an exact computational reduction, not a new factoring exponent.

## 1. The universal operator kernel

Let \(R\) be a discrete valuation ring of characteristic \(p\),
with fraction field \(K\) and parameter \(z\). A fixed Puiseux
extension may be used for rational heights. In each of the variable
families \(Y,X,U,V\), there are \(r\) variables, and each
variable has \(p\)-th power zero. Put

\[
\mathcal A=\{0,\ldots,p-1\}^r,\qquad
\mathbf t=(p-1,\ldots,p-1),\qquad
E(H)=\sum_{d=0}^{p-1}\frac{H^d}{d!}.
\tag{1}
\]

Let \(C\in\operatorname{GL}_r(R)\), and let
\(F=(F_{\mu\nu})_{\mu,\nu\in\mathcal A}\) be any
matrix over \(K\). Define the coefficient kernel

\[
K_F(Y,X;U,V)=
 \prod_{j=1}^rE(V_jX_j)
 \sum_{\mu,\nu\in\mathcal A}
 F_{\mu\nu}
 \frac{(Y+CV)^\mu(U-X)^\nu}{\nu!}.
\tag{2}
\]

Here vector powers mean products of the indicated component powers,
and \(\nu!=\prod_j\nu_j!\). Flatten (2) with row labels
\((Y,X)\) and column labels \((U,V)\), using ordinary monomial
coefficient bases. Denote the resulting matrix by \(\mathcal K_F\).

Choose heights \(h_1,\ldots,h_r\), and set

\[
A=\operatorname{diag}(z^{h_1},\ldots,z^{h_r}),\qquad
\mathcal D_A(Y^iX^j)=z^{h\cdot i-h\cdot j}Y^iX^j.
\tag{3}
\]

Assume the additional weighted regularity condition

\[
A^{-1}CA\in\operatorname{GL}_r(R).
\tag{4}
\]

Both parts of the assumption matter: a constant matrix need not be
regular after unequal variable heights have been applied.

## 2. The reduced differential operator

Let \(\partial_j\) be ordinary differentiation in the truncated
polynomial algebra
\(R[Y_1,\ldots,Y_r]/(Y_1^p,\ldots,Y_r^p)\).
Differentiation is well defined in characteristic \(p\), and
\([\partial_i,Y_j]=\delta_{ij}\) holds on the whole algebra,
including its top-degree basis vectors.

For a matrix unit \(E_{\mu\nu}\), define

\[
\boxed{
\Psi_C(E_{\mu\nu})
=(-1)^r\mu!\,
 \partial^{\mathbf t-\mu}
 \prod_{j=1}^r\bigl(-(C^{-1}Y)_j\bigr)^{p-1-\nu_j}.
}
\tag{5}
\]

The product on the right is a multiplication operator, applied
first; the derivatives are applied afterward. Extend (5) linearly
to \(\Psi_C(F)\). It is a matrix of order \(p^r\).

Let

\[
M_F=\operatorname{diag}_{i\in\mathcal A}(z^{h\cdot i})
          \Psi_C(F).
\tag{6}
\]

**Weighted Weyl reduction.** If the Smith orders of \(M_F\)
are \(\lambda_1,\ldots,\lambda_{p^r}\), including infinite
orders when appropriate, the Smith orders of
\(\mathcal D_A\mathcal K_F\) are exactly

\[
\boxed{
\left\{\lambda_\ell-h\cdot j:
       1\le\ell\le p^r,\ j\in\mathcal A\right\}.
}
\tag{7}
\]

In particular, when \(h=0\), every reduced order occurs exactly
\(p^r\) times. The theorem itself allows non-preserving
occupation kernels. For an actual crossing, endpoint admissibility
and all infinity row and column heights remain separate obligations.

## 3. The reference kernel is a regular unit

Write \(\mathcal K_* = \mathcal K_{E_{\mathbf t\mathbf t}}\).
Multiplication by any component of \(Y+CV\) or \(U-X\)
annihilates its generating polynomial, because the corresponding
linear form already has exponent \(p-1\). Consequently

\[
\mathcal K_* M_{U_j}^{\mathsf T}=M_{X_j}\mathcal K_*,
\qquad
\mathcal K_* M_{V_j}^{\mathsf T}
=-\sum_k(C^{-1})_{jk}M_{Y_k}\mathcal K_*.
\tag{8}
\]

The transposes occur because the variables in these identities
belong to the column side of a coefficient matrix.

The column labelled \(U^{\mathbf t}V^{\mathbf t}\) has a
constant row coefficient equal to

\[
\frac{\det(C)^{p-1}}{((p-1)!)^r},
\tag{9}
\]

which is a unit in \(R\). The identity
\(\prod_j(CV)_j^{p-1}=\det(C)^{p-1}V^{\mathbf t}\)
follows from elementary shears, permutations, and diagonal scalings;
each shear has no extra surviving term in top degree.

The column top vector is cyclic under the commuting lowering
operators \(M_U^{\mathsf T},M_V^{\mathsf T}\). Equations
(8) identify them with independent invertible linear combinations
of row multiplication operators. Its image is a polynomial with
unit constant term, so it is cyclic for the row truncated algebra
as well. This proves
\(\mathcal K_*\in\operatorname{GL}_{p^{2r}}(R)\).

## 4. An explicit derivative conjugation

Set

\[
\mathcal W=
\prod_{j=1}^r
 E\left(\partial_{Y_j}\sum_k C_{jk}\partial_{X_k}\right).
\tag{10}
\]

The factors commute, and changing every sign gives the inverse.
They are finite operators over \(R\). Direct commutation gives

\[
\mathcal W Y_j\mathcal W^{-1}
=Y_j+\sum_k C_{jk}\partial_{X_k},\qquad
\mathcal W\partial_{Y_j}\mathcal W^{-1}=\partial_{Y_j}.
\tag{11}
\]

There is no characteristic-\(p\) exponential remainder in (11).
For one factor with \(B=\sum_kC_{jk}\partial_{X_k}\), the
potential final term contains \(B^p=0\).

Differentiate the reference kernel in \(X_j\). The derivative
of its \(U-X\) factor is minus the \(U_j\) derivative, and
the derivative of \(\prod E(V_kX_k)\) contributes \(V_j\).
Using (8) gives

\[
\mathcal K_*\partial_{U_j}^{\mathsf T}\mathcal K_*^{-1}
=-\partial_{X_j}-(C^{-1}Y)_j.
\tag{12}
\]

Differentiating the reference kernel in \(Y\) and \(U\)
also gives the exact factorial identity

\[
\mathcal K_{E_{\mu\nu}}
=(-1)^r\mu!\,
 \partial_Y^{\mathbf t-\mu}\,
 \mathcal K_*\,
 (\partial_U^{\mathsf T})^{\mathbf t-\nu}.
\tag{13}
\]

This uses only \((p-1)!=-1\) in the prime field. Combining
(11)--(13) proves the matrix identity

\[
\boxed{
\mathcal K_F
=\mathcal W\bigl(I_X\otimes\Psi_C(F)\bigr)
 \mathcal W^{-1}\mathcal K_*.
}
\tag{14}
\]

The tensor-product notation in (14) means that the reduced operator
acts on the \(Y\) variables and leaves the \(X\) monomial
labels fixed; a harmless basis permutation gives the usual Kronecker
ordering.

## 5. Why the factorization respects the heights

A term \(C_{jk}\partial_{Y_j}\partial_{X_k}\) in (10)
changes the row height by \(-h_j+h_k\). Condition (4) says
exactly that its coefficient after this weighted change is regular.
It follows that

\[
\mathcal D_A\mathcal W^{\pm1}\mathcal D_A^{-1}
\in\operatorname{GL}_{p^{2r}}(R).
\tag{15}
\]

Multiply \(\mathcal D_A\mathcal K_F\) on the left by
\(\mathcal D_A\mathcal W^{-1}\mathcal D_A^{-1}\), and
on the right by \(\mathcal K_*^{-1}\mathcal W\). Both
operations preserve Smith orders. Equation (14) becomes

\[
\mathcal D_A(I_X\otimes\Psi_C(F)).
\tag{16}
\]

For a fixed auxiliary \(X\) occupation \(j\), its block is
\(z^{-h\cdot j}M_F\). This proves (7).

If \(F\) preserves total occupation, the number of derivatives
in (5) equals the multiplication degree. Hence \(\Psi_C(F)\)
preserves total \(Y\) occupation. One can calculate its Smith
orders in blocks of size at most
\([x^d]\,(1+x+\cdots+x^{p-1})^r\), and then use (7).
If a reduced block has occupation \(d\), the auxiliary state
\(j\) belongs to original charge \(d-|j|\).

## 6. Application to the coupled Pascal cores

The proved coupled reductions have an active Gaussian with entries
\(A=\operatorname{diag}(z^{h_j})\) on \(Y,U\), and a
diagonal matrix \(D\) on \(X,V\). Its diagonal entries have
orders \(b_i\), possibly multiplied by retained regular units.
Their individual forms are

\[
L=AY+cF_V^{\mathsf T}V,\qquad
R=U-A^{-1}F_X^{\mathsf T}X,
\tag{17}
\]

where \(F_X,F_V\) are regular invertible matrices and
\(c=1-z^{2P}\). In the mixed geometry, \(F_V\) and \(D\)
include the already-proved auxiliary unit clearing.

Make the regular changes
\(\widehat X=F_X^{\mathsf T}X\) and
\(\widehat V=cF_V^{\mathsf T}V\). The other Gaussian matrix
is then

\[
S=c^{-1}F_X^{-1}D F_V^{-\mathsf T}.
\tag{18}
\]

Temporarily write \(Y'=AY\), \(X'=A^{-1}\widehat X\),
retaining their row height matrix \(\mathcal D_A\) explicitly.
Removing the selected Gaussian directions leaves two regular
internal Gaussian factors and the mixed Gaussian
\(\Gamma(X'^{\mathsf T}(I+AS)\widehat V)\).
The internal row factor is \(\Gamma(Y^{\mathsf T}\widehat X)\),
which is regular in the original variables. The column factor is
regular as well.

Set

\[
Z=I+AS,\qquad C=Z^{-1},\qquad V'=Z\widehat V.
\tag{19}
\]

This gives exactly the universal kernel (2), with its explicit row
weight (3). Every step except that displayed row weight is regular.

For a selected identity index \(j\) and retained auxiliary
index \(i\), the physical heights satisfy

\[
h_j=P(2j+1-J),\qquad b_i=P(2i+1+J),\qquad
h_j+b_i=2P(j+i+1)>0.
\tag{20}
\]

Since the matrices in (18) are regular units, every entry of
\(AS\) and \(SA\) has positive valuation. Thus

\[
C=(I+AS)^{-1}\in\operatorname{GL}_r(R),\qquad
A^{-1}CA=(I+SA)^{-1}\in\operatorname{GL}_r(R).
\tag{21}
\]

This proves the weighted condition, for every aggregate gauge
\(J\), in all three implemented geometries. The unchanged
Gaussian spectators are then restored exactly as before. Infinity
profiles continue to use the original weighted occupation matrix;
the zero auxiliary units are not removed at infinity.

## 7. Recovery of the finite-difference theorem

For one mode and a diagonal occupation matrix
\(F_{dd}=\alpha_d\), put \(k=p-1-i\), where \(i\)
is the exponent of the reduced monomial \(Y^i\). Formula (5)
is diagonal, with eigenvalue

\[
\sum_{h=0}^k(-1)^h\binom kh c^{-h}\alpha_{p-1-h}
=c^{-k}\beta_k.
\tag{22}
\]

The factor \(c^{-k}\) is a unit. Equations (6)--(7) give
\(h(r-k)+\operatorname{ord}\beta_k\), for
\(0\le r,k<p\), exactly the earlier prime-state formula.
This independently checks the factorials, signs, and auxiliary
height shifts in the operator extension.

## 8. Complete-polynomial implementation

[coupled_identity_weyl_profiles.py](../tools/coupled_identity_weyl_profiles.py)
implements the reduction using polynomials in \(\tau=q^P\).
Let

\[
\Delta=c\det F_X\det F_V,
\qquad
N=\Delta I+A\operatorname{adj}(F_X)
                   D\operatorname{adj}(F_V)^{\mathsf T}.
\tag{23}
\]

Then \(C^{-1}=Z=N/\Delta\), and \(\Delta\) is a proved
zero unit. Multiplying the entire reduced operator by
\(\Delta^{r(p-1)}\) clears every local denominator in (5).
The code therefore builds complete polynomial matrices and sends
their full coefficients to the exact Smith engine. It does not
truncate a rational series or infer rank from a finite jet.

This auxiliary clearing changes neither the zero orders nor the
original physical polynomial tensor. It is a local representation
step, so its degree is not substituted for the degree of the
original crossing. The original endpoint and polynomial construction
remain those of the
[operator-filter theorem](coupled-identity-operator-filters.md).

The [checker](../tools/check_coupled_identity_weyl_profiles.py) and
[saved evidence](../evidence/coupled-identity-weyl-profiles.json)
verify 87 complete symbolic matrix identities over
\(\mathbb F_p[c,c^{-1}]\), for \(p=2,3,5,7\), and 87
diagonal-transform coefficient identities. The extended replay also
compares 16 old full-core profiles, containing 690 zero orders,
with their reduced profiles. It covers arbitrary polynomial
occupation transfers, nonidentity endpoint repairs, and all three
geometries.

An additional
[portable screen record](../evidence/coupled-identity-operator-screens.json)
contains 1,008 complete profiles calculated by the original
\(p^{2r}\) operator engine. These cover 64 initial and 64 mixed
two-mode kernels and 880 three-mode kernels at \(n=32\).
The checker reconstructs each finite catalogue and matches all
recorded charge profiles against the Weyl calculation. The associated
floating-point screens found no pressure improvement; their best
critical pressure was the Gaussian value, approximately
\(0.03306796250845\). This is finite search evidence, with no
claim that the entire operator family has been excluded.

Run the complete portable check with

```sh
python3 tools/check_coupled_identity_weyl_profiles.py --extended \
  --screen evidence/coupled-identity-operator-screens.json \
  --output /tmp/coupled-identity-weyl-profiles-replay.json
```

The finite-field calculations validate the stated formulas; the
all-prime and all-length claims follow from the proof above.
