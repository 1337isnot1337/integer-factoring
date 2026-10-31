# Identity cores a fixed distance below a prime power

## Result

Let \(p\) be prime, \(N=p^k\), and \(n=N-d\). Select
\(a\) initial and \(b\) terminal identity directions of the
physical interval of length \(n\), and put \(r=a+b\). Assume

\[
d\ge0,\qquad a,b\ge0,\qquad r\ge1,\qquad r+d\le n.
\tag{1}
\]

For a physical crossing, require the operator kernel to satisfy the
existing invertible, energy-filtered endpoint condition from the
[operator construction](coupled-identity-operator-filters.md). Its
constant physical repair must be regular and invertible at both cuts.
Every such occupation-preserving operator on these \(r\) directions
has a regular zero reduction to **\(2r+d\) actual active modes**. Its
profile can be calculated through a Weyl matrix of order
**\(p^{r+d}\)** and exact division by \(d\) known Gaussian digit
factors. The physical crossing does not acquire the formal modes used
by that calculation.

The local algebraic identities also apply to arbitrary
occupation-preserving matrices, but such matrices can have singular
endpoints or infinite Smith orders. They are not automatically
admissible physical crossings.

Moreover, the resulting Weyl parameter \(C^{-1}\) differs from
\((1-\tau^N)^{-1}I\) by a matrix of rank at most \(d\) over the
rational function field. Thus a fixed deficit gives a bounded-rank
departure from the scalar prime-power geometry.

This is a constructive reduction theorem, not an improved factoring
bound. The accompanying code certifies the complete polynomial
geometry, contraction, and rank statement. It does not claim an
exhaustive operator or pressure search in the enlarged family.

## 1. Augmented feature columns and weighted interpolation

Use \(\tau=z^{2P}\), \(c=1-\tau\), with the existing odd
ramification \(P\), and

\[
f_j(i;c)=\sum_h\binom jh\binom ih c^h.
\tag{2}
\]

The actual selected indices are

\[
S=\{0,\ldots,a-1\}\cup\{n-b,\ldots,n-1\}.
\tag{3}
\]

Append the formal indices \(D=\{n,\ldots,N-1\}\). The augmented
terminal block is then exactly \(\{N-(b+d),\ldots,N-1\}\).
There are \(R=r+d\) augmented columns. Retain auxiliary rows

\[
T=\{0,\ldots,b+d-1\}\cup\{n-a,\ldots,n-1\}.
\tag{4}
\]

The two row blocks are disjoint by (1). Let \(E\) be the form
matrix on the actual rows \(0\le i<n\) and augmented columns
\(S\cup D\), and put \(F=E_T\).

The parent prime power gives

\[
f_{N-1-h}(i;c)=\tau^i f_h(i;-c/\tau).
\tag{5}
\]

Thus the augmented terminal columns, whose reversed degrees are
\(0,\ldots,b+d-1\), interpolate on the low retained rows with
coefficients \(\tau^{i-t}\) times constants. Their low square
minor is a zero unit. The initial columns have an invertible constant
minor on any consecutive block of \(a\) rows; its integer determinant
is one. Every augmented terminal entry on the high retained rows
vanishes modulo \(\tau\).

The proof of the [mixed-end interpolation theorem](split-identity-mode-reduction.md)
therefore applies unchanged with \(b+d\) low rows and \(a\) high
rows. It uses the parent reversal identity only; it does not require
the number of physical rows to equal \(N\). It proves that \(F\)
is a zero unit and that, for every discarded actual row,

\[
E_i=B_iF,\qquad
\operatorname{ord}_\tau(B_i)_t\ge\max(0,i-t).
\tag{6}
\]

For clarity, the high interpolation coefficients are obtained from a
regular Schur complement of the initial high minor. The low coefficients
are terminal interpolation minus a regular combination of higher-row
terminal interpolation. Each correction has order greater than the
required \(i-t\). This is also valid when one of the two actual
end blocks is empty.

## 2. The complete active data and weighted Weyl hypothesis

Put \(\delta=\det F\), \(N_i=E_i\operatorname{adj}(F)\).
As in the [one-deficit reduction](near-power-identity-operator-reduction.md),
the polynomial coefficient matrices are

\[
F_X=F,\qquad
(F_V)_{tj}=\delta F_{tj}
 +\sum_{i\notin T}\tau^{t-i}(N_i)_t E_{ij}.
\tag{7}
\]

They retain every coefficient. Formula (6) makes (7) polynomial.
The accompanying auxiliary Gaussian has diagonal
\(\delta\tau^{t+1}\) at \(J=1\).

Both matrices in (7) are zero units. Indeed,
\(F_V=\delta W F\), where

\[
W=I+D_TB^{\mathsf T}D_{T^c}^{-1}B,
\qquad D_i=\tau^{i+1}.
\tag{8}
\]

Its diagonal corrections have positive valuation and its entries below
the diagonal, ordered by retained index, have positive valuation.
Thus \(W\) is upper triangular with diagonal one modulo \(\tau\).
Equation (6) also makes the required Gaussian companion and its
inverse regular. Keeping only the \(r\) actual operator directions
leaves \(r\) identity pairs and \(R\) auxiliary pairs, or
\(2r+d\) active modes.

For the augmented square core, put
\(A=\operatorname{diag}_{j\in S\cup D}(\tau^j)\). Its Weyl
parameter has the form

\[
C^{-1}=I+A\mathcal S,\qquad
\mathcal S=c^{-1}F_X^{-1}
  \operatorname{diag}_{t\in T}(\delta\tau^{t+1})F_V^{-\mathsf T}.
\tag{9}
\]

Every entry of \(\mathcal S\) has positive valuation because
\(F_X,F_V,c,\delta\) are zero units and all \(t+1\ge1\).
All selected or formal indices are nonnegative. Therefore both
\(I+A\mathcal S\) and
\(A^{-1}C^{-1}A=I+\mathcal S A\) are the identity modulo
\(\tau\). Their inverses are regular. This proves both weighted
assumptions of the [Weyl theorem](identity-operator-weyl-reduction.md)
without a separate case analysis of the deficit.

## 3. Formal Gaussian factors are removed exactly

Append one independent local Gaussian \(E(\tau^jY_jU_j)\) for
each \(j\in D\), where

\[
E(H)=\sum_{h=0}^{p-1}\frac{H^h}{h!}.
\tag{10}
\]

Extend the actual occupation kernel by the identity on these \(d\)
formal digits. Its extra occupation sum is \(\prod_{j\in D}E(Q_j)\),
which cancels the removed factors \(\prod_{j\in D}E(-Q_j)\).
The equality \(E(Q)E(-Q)=1\) holds because \(Q^p=0\) in the
truncated algebra. The augmented tensor therefore equals the original
rectangular tensor times the \(d\) independent Gaussian pairs.

The square core has \(R\) identity and \(R\) auxiliary modes, so
the Weyl matrix has order \(p^R\). Its full reconstructed profile
has \(p^{2R}\) orders. At \(J=1\), its charge/order generating
polynomial satisfies

\[
\widehat{\mathcal Z}(u,x)=\mathcal Z(u,x)
 \prod_{j=n}^{N-1}\left(\sum_{h=0}^{p-1}u^h x^{2Pjh}\right).
\tag{11}
\]

The coefficients here are integer multiplicities, not field elements.
Divide one factor at a time in ascending charge. For a factor of
order \(v=2Pj\), use

\[
\mathcal Z_c=\widehat{\mathcal Z}_c
                  -\sum_{h=1}^{p-1}x^{hv}\mathcal Z_{c-h}.
\tag{12}
\]

Every recovered multiplicity is nonnegative and the final extra
charge blocks vanish. After all \(d\) divisions, exactly
\(p^{r+R}=p^{2r+d}\) orders remain. General \(J\) is restored
by adding \(P(1-J)\) times the charge to each recovered order.

The formal pairs never enter the global polynomial tensor. The
coefficientwise endpoint, physical alphabet, original degree bound,
and infinity matrix involve only \(r\) actual directions. The
extra \(d\) retained auxiliary modes are real zero-cut modes and
must remain in its partition until the usual spectator accounting.

## 4. The departure from the scalar parameter has rank at most d

Let \(U\) be the \(d\)-by-\(R\) matrix of the dropped parent
rows \(i=n,\ldots,N-1\), restricted to the augmented columns.
Let \(D_{\rm drop}=\operatorname{diag}(\tau^{i+1})\), and let
\(H=E^{\mathsf T}D^{-1}E\) be the actual selected Gram matrix.
The parent [Pascal Gram identity](prime-power-joint-identity-filters.md)
gives

\[
H=H_0-U^{\mathsf T}D_{\rm drop}^{-1}U,
\qquad H_0=\tau^{-N}c^{N-1}A.
\tag{13}
\]

The contraction in (9) is also \((cH)^{-1}\). Define
\(K_0=cH_0=\tau^{-N}(1-\tau^N)A\). Both \(H\) and
\(K_0\) are invertible over the rational function field. The inverse
identity gives

\[
(cH)^{-1}-K_0^{-1}
=(cH)^{-1}\,cU^{\mathsf T}D_{\rm drop}^{-1}U\,K_0^{-1}.
\tag{14}
\]

The right side has rank at most \(d\). Since
\(I+AK_0^{-1}=(1-\tau^N)^{-1}I\), equations (9) and (14)
prove

\[
\boxed{\operatorname{rank}_{\mathbb F_p(\tau)}
 \left(C^{-1}-(1-\tau^N)^{-1}I\right)\le d.}
\tag{15}
\]

This rank statement does not assert that the correction can be
discarded, or that its eigenvectors respect every occupation gauge.
It supplies a smaller rational representation for further exact
operator calculations. At \(d=0\) it recovers the scalar prime-power
parameter; at \(d=1\) it recovers the separately proved rank-one
formula.

## 5. Exact implementation and controls

[`deficit_identity_mode_geometry.py`](../tools/deficit_identity_mode_geometry.py)
constructs (7), its complete polynomial Weyl numerator and zero-unit
denominator, and the exact multiplicity division (12). Its
[`checker`](../tools/check_deficit_identity_mode_geometry.py) covers
characteristics two, three, and five; deficits zero through three;
and pure-terminal and mixed-end selections.

The saved
[`evidence`](../evidence/deficit-identity-mode-geometry.json) records
468 direct form identities, 174 parent Gram entries, 174 inverse
contraction identities, 246 vanishing minors certifying the stated
rank bounds in the checked cases, and nine exact dummy-factor
divisions. All polynomial comparisons are exact. These finite checks
support the formulas; the uniform statements follow from Sections
1--4 and the cited regular-concentration and Weyl theorems.

```bash
python3 tools/check_deficit_identity_mode_geometry.py
```

The full original physical-cut replay is currently supplied by the
separate binary one-deficit package. This general constructor does not
claim additional full-tensor replays or completed pressure searches.
