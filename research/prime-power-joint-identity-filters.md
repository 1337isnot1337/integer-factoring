# Every identity direction separates at a prime-power interval length

## Result and scope

Let the coefficient field have characteristic \(p\), let the Pascal
interval length be \(n=p^e\), and let the ramification \(P\) be a positive
odd integer. A single regular change at the zero cut
separates **all** polynomial identity summands simultaneously. Thus any
set of \(r\) identity indices, including interior or irregularly spaced
indices, has exactly \(2r\) active zero modes. The selected indices need
not form an initial, terminal, or split interval.

For an arbitrary joint diagonal occupation filter, the complete zero
profile is a tensor-product finite-difference transform of its coefficient
polynomials. Its polynomial part is independent of the selected indices.
Only explicit linear heights depend on those indices. This permits one
coefficient calculation to be reused for every selected subset.

The result enlarges the constructive search family and removes a matrix
Smith calculation from this family. It does not assert an improved
factoring exponent. Endpoint repair and the full ordinary polynomial
degree are retained below.

## 1. A finite Pascal orthogonality identity

Write \(\tau=z^{2P}\), \(c=1-\tau\), and

\[
H_{ij}=f_j(i;c)=\sum_{h\ge0}\binom ih\binom jh c^h,
\qquad 0\le i,j<n.
\tag{1}
\]

All binomial coefficients are reduced in characteristic \(p\). Put
\(N=n-1\). Since \(n\) is a power of the characteristic,
\(\binom Nj=(-1)^j\). More generally, for \(0\le h\le j<n\),

\[
\binom{N-h}{j-h}=(-1)^{j-h}\binom jh.
\tag{2}
\]

For example, (2) follows by comparing coefficients below degree \(n\)
in \((1-X)^{N-h}=(1-X^n)(1-X)^{-h-1}\). Expanding the right side of
the following identity and applying (2) gives the finite generating
polynomial

\[
F_i(X):=\sum_{j=0}^{n-1}f_j(i;c)X^j
       =(1-X)^{N-i}(1-\tau X)^i.
\tag{3}
\]

Set \(A=(1-X)(1-Y)\) and
\(B=\tau^{-1}(1-\tau X)(1-\tau Y)\). Frobenius gives

\[
\begin{aligned}
\sum_{i=0}^{n-1}\tau^{-i}F_i(X)F_i(Y)
 &=\sum_{i=0}^{n-1}A^{N-i}B^i\\
 &=(A-B)^N\\
 &=c^N\tau^{-N}(1-\tau XY)^N\\
 &=c^N\tau^{-N}\sum_{j=0}^{n-1}\tau^j(XY)^j.
\end{aligned}
\tag{4}
\]

The second equality is the polynomial identity
\((A-B)\sum_iA^{N-i}B^i=A^n-B^n=(A-B)^n\).
The sign \((-1)^N\) is one in the coefficient field. Comparing the
coefficients of \(X^jY^k\) proves

\[
\boxed{
H^{\mathsf T}\mathrm{diag}_{0\le i<n}(\tau^{-i})H
 =\mathrm{diag}_{0\le j<n}
       \left(c^{n-1}\tau^{j+1-n}\right).
}
\tag{5}
\]

This identity is finite and exact. It involves no limiting series or
numerical approximation.

## 2. One regular change handles every selected subset

Let \(R\) be the permutation matrix that reverses the \(n\) indices,
and put \(B=HR\). Reversing the coefficients in (3) gives

\[
B_{ik}=f_{n-1-k}(i;c)
       =\tau^i f_k(i;-c/\tau).
\tag{6}
\]

Each term on the right has a factor \(\tau^{i-h}\) with
\(h\le\min(i,k)\). Consequently

\[
\mathrm{ord}_{\tau}B_{ik}\ge\max(0,i-k).
\tag{7}
\]

The ordinary Pascal factorization
\(H=\mathcal P\mathrm{diag}(c^h)\mathcal P^{\mathsf T}\), with
\(\mathcal P_{ih}=\binom ih\), shows

\[
\det B=\pm c^{n(n-1)/2}.
\tag{8}
\]

Thus \(B\) is regular with unit determinant at zero.

Use the standard aggregate zero charge \(J\), and write

\[
a_j=P(2j+1-J),\qquad b_i=P(2i+1+J),\qquad
D=\mathrm{diag}(z^{b_i}).
\tag{9}
\]

The original zero reduction from the
[coupled identity theorem](coupled-identity-filter-reduction.md) has
Gaussian \(X^{\mathsf T}DV\), together with the diagonal \(Y,U\)
Gaussian, and identity forms

\[
Q_j=\left(z^{a_j}Y_j+c(H^{\mathsf T}V)_j\right)
     \left(U_j-z^{-a_j}(H^{\mathsf T}X)_j\right).
\tag{10}
\]

Apply the Gaussian-preserving change

\[
X_{\rm old}=L^{\mathsf T}X_{\rm new},\qquad
V_{\rm old}=R_gV_{\rm new},\qquad
L=B^{-1},\quad R_g=D^{-1}BD.
\tag{11}
\]

The common scalar in \(D=z^{P(1+J)}\mathrm{diag}(\tau^i)\)
cancels in the conjugation. Inequality (7) makes every entry of
\(D^{-1}BD\) regular. Its determinant is the same unit as \(\det B\).
Therefore **all four matrices** \(L,L^{-1},R_g,R_g^{-1}\) are regular.
Their induced maps on the truncated symmetric algebras are regular and
invertible as well: the relations \(X_i^p=0\) are preserved under every
linear substitution in characteristic \(p\).

Now \(LH=R\). Equation (5), with the common scalar in \(D\) included,
also gives

\[
H^{\mathsf T}D^{-1}BD=c^{n-1}R.
\tag{12}
\]

Substituting (12) into (10), and using \(c^n=1-\tau^n\), proves the
simultaneous normal form

\[
\boxed{
Q_j=\left(z^{a_j}Y_j+(1-\tau^n)V_{n-1-j}\right)
     \left(U_j-z^{-a_j}X_{n-1-j}\right),
\qquad 0\le j<n.
}
\tag{13}
\]

The two Gaussian orders paired with \(Q_j\) are \(a_j\) and
\(b_{n-1-j}\), whose sum is \(2nP\). Distinct \(j\) use disjoint
variables. Any selected set \(S=\{j_1,\ldots,j_r\}\) therefore leaves
exactly these \(2r\) modes active, with every other mode a Gaussian
spectator. In the [Weyl reduction](identity-operator-weyl-reduction.md),
the associated parameter is simply

\[
C=(1-\tau^n)I_r.
\tag{14}
\]

## 3. The exact joint finite-difference profile

Put \(x=-v=z^2\), \(m=nP\), \(c_0=1-x^m\), and
\(\mathcal A=\{0,\ldots,p-1\}^r\). Let
\(\alpha_{\mathbf d}(x)\) be arbitrary Laurent polynomials over the
coefficient field, with every \(\alpha_{\mathbf d}(1)\ne0\).
For \(\mathbf k\in\mathcal A\), define

\[
\boxed{
\beta_{\mathbf k}(x)
=\sum_{0\le\mathbf h\le\mathbf k}
 (-1)^{|\mathbf k|-|\mathbf h|}
 \left(\prod_{\ell=1}^r\binom{k_\ell}{h_\ell}\right)
 c_0^{|\mathbf h|}
 \alpha_{\mathbf p-\mathbf1-\mathbf k+\mathbf h}(x).
}
\tag{15}
\]

This is the tensor product of the exact
[one-direction transform](prime-identity-filter-diagonalization.md).
It remains valid when the joint weights do not factor. Equivalently,
the diagonal entry of \(\Psi_C(\mathrm{diag}\alpha)\) on the
monomial indexed by \(\mathbf p-\mathbf1-\mathbf k\) is
\(c_0^{-|\mathbf k|}\beta_{\mathbf k}\). The factor
\(c_0^{-|\mathbf k|}\) is a zero unit and does not affect its order.

Let \(B_{\mathbf k}=\mathrm{ord}_x\beta_{\mathbf k}\). Every such
order is finite, because
\(\beta_{\mathbf k}(1)=(-1)^{|\mathbf k|}
\alpha_{\mathbf p-\mathbf1-\mathbf k}(1)\ne0\).
The exact active zero orders are

\[
\boxed{
\left\{2B_{\mathbf k}
 +\sum_{\ell=1}^r a_{j_\ell}(d_\ell-k_\ell):
       \mathbf d,\mathbf k\in\mathcal A\right\}.
}
\tag{16}
\]

The charge of this order is \(|\mathbf d|-|\mathbf k|\). At infinity
the active orders are

\[
\boxed{
\left\{P\sum_{\ell=1}^r(2j_\ell+1-2n)d_\ell
              -2\deg_x\alpha_{\mathbf d}:
       \mathbf d\in\mathcal A\right\}.
}
\tag{17}
\]

Equations (15)--(17), with the explicit unselected spectators, contain
both complete cuts. The tool
[`prime_power_identity_filters.py`](../tools/prime_power_identity_filters.py)
implements the formulas without series truncation or a Smith oracle.

For \(\alpha\equiv1\), (15) gives
\(\beta_{\mathbf k}=(-x^m)^{|\mathbf k|}\), so (16) recovers precisely
the selected Gaussian modes \(a_{j_\ell},2nP-a_{j_\ell}\).

## 4. Endpoint, field, degree, and search accounting

The physical polynomial filter is the one in the
[coupled construction](coupled-identity-filter-reduction.md), applied
to the original polynomial \(Q_{j_\ell}\). The change (11) is used
only to calculate the zero profile. It introduces no rational factor
into the physical seed.

At the endpoint, multiply each physical selected-Y occupation label
\(\mathbf d\) by the constant
\(\alpha_{\mathbf d}(1)^{-1}\). This restores the original crossing
coefficient by coefficient. The repair commutes with all diagonal
physical gauges, so it changes neither cut nor the ordinary degree.
A fixed extension of the coefficient field has the same properties.

One common Laurent monomial clearing adds the same integer to every
\(B_{\mathbf k}\) and every \(\deg\alpha_{\mathbf d}\); its two
profile shifts cancel in the product cost. A general rational
denominator must still be cleared once with its complete degree cost.
The established ordinary-degree bound remains

\[
P\deg G+2nPr(p-1)+\max_{\mathbf d}\deg\alpha_{\mathbf d}
\tag{18}
\]

after the common Laurent clearing. Fixed finite alphabets, ramification,
and this full degree remain charged in the crossing construction.

Finally, (15) is independent of \(S\). An ordered beta-valuation vector
and the ordered degree vector can be cached once and evaluated on every
selected subset. Continuous index or height optimization may help choose
candidates, but an admissible witness must use distinct actual indices
and must replay (16)--(18) with its explicit repaired endpoint.

## 5. Exact controls and a strict row improvement with its full cost

The [portable checker](../tools/check_prime_power_identity_filters.py)
and [saved record](../evidence/prime-power-identity-filters.json) contain:

- 505 complete polynomial Gram identities and 505 weighted-entry checks,
  over characteristics 2, 3, 5, and 7;
- 15 comparisons with the previous independent active-core engines,
  including F4 coefficients and a qutrit core;
- 14 complete original physical comparisons of **both** weighted cuts,
  including interior and nonconsecutive indices, three selected modes,
  all four binary modes at \(n=4\), and a qutrit matrix of order 729;
- 859,377 coefficientwise endpoint equalities, including a nontrivial
  qutrit endpoint repair.

The literal physical checks construct the original Pascal tensor before
any of the new changes of basis. They therefore test the physical gauges,
ordinary Laurent clearing, and endpoint, as well as the reduced formula.
The qutrit replay is enabled by `--include-qutrit`.

The polynomial lattice search also supplies a simple exact example of
why a stronger zero row target alone is insufficient. For binary digits,
choose index zero among the selected modes, set \(m=nP\), and take

\[
\alpha_{(0,\mathbf d)}=1+x^m+x^{m+1},\qquad
\alpha_{(1,\mathbf d)}=1.
\tag{19}
\]

Every endpoint value is one, and
\(B_{\mathbf k}=m|\mathbf k|+k_1\): all beta rows with first digit
one gain one full \(x\)-order. At \(n=32,P=1,S=(0,30,31),J=1\), the
eight reduced Gaussian row orders, grouped by occupation, change from

\[
(96,64,94,95,62,63,93,61)
\quad\hbox{to}\quad
(97,64,95,96,62,63,94,61).
\tag{20}
\]

However the four first-digit-zero weights in (19) have full degree
\(m+1\). The exact active cost ratios are

\[
R_0=\frac{1+s^{-b-2}}{1+s^{-b}},\qquad
R_\infty=\frac{s^{2m+2}+s^{2m-P}}{1+s^{2m-P}},
\qquad b=P(2n-1+J).
\tag{21}
\]

For every real \(J\) and \(s>1\), \(R_0\ge s^{-2}\). When
\(n\ge2\), one also has \(2m-P\ge P+2\), so
\(R_\infty\ge s^{P+2}\). Therefore
\(R_0R_\infty\ge s^P>1\). This feasible strict row improvement is
retained with its complete infinity cost; it is not a pressure witness.

## 6. Reusing one coefficient catalogue over all selected indices

The [fast pressure helper](../tools/joint_identity_filter_pressure.py)
factors out the unchanged selected-Y partition. The
[independent pressure controls](../tools/check_joint_identity_filter_pressure.py)
compare its values with complete physical partition sums, its compiled
beta valuations with literal sparse polynomials, and its index derivative
with finite differences. The
[saved pressure record](../evidence/joint-identity-filter-pressure.json)
contains 60 complete comparisons, with maximum floating discrepancy below
\(8\cdot10^{-15}\). These controls do not turn a floating sign into a
rigorous pressure certificate.

The [all-index driver](../tools/search_all_index_joint_filters.py) first
forms exact ordered beta/degree classes, then evaluates the same classes
over every stated subset. Its initial finite screens include:

| Coefficient family | Interval | Selected modes | Exact coefficient cases | Ordered profile classes | All selected subsets | Scale/charge points |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Complete stated F4 three-mode family | 32 | 3 | 369,603 | 3,850 | 4,960 | 15 |
| All binary two-level indicator tables, both signs | 16 | 4 | 65,536 | 65,535 | 1,820 | 4 |

The extra binary case is the repeated all-zero table. The three-mode
exponent family is closed under every permutation of the modes; its
nonzero colors are fully enumerated. The indicator family is also closed
under these permutations.

The records are
[the complete three-mode subset screen](../evidence/all-index-three-mode-f4-n32-screen.json)
and [the four-mode indicator screen](../evidence/all-index-four-mode-binary-n16-screen.json).
Neither screen found a same-gauge cost decrease or a negative critical
pressure beyond its stated \(10^{-9}\) detection tolerance. They cover
the listed coefficient families and finite scale/charge grids, not all
polynomial filters or all real scales.
