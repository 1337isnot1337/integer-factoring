# A sharper copied-wire obstruction from one truncated multiplication map

The [copied-wire theorem](shift-register-copy-wire-barrier.md) bounds
the bare two-place Smith criterion for every full-cycle reversible
shift register. This note improves its finite-size estimate from

\[
 B(s)\ge d^{m+E}/J^3
 \quad\hbox{to}\quad
 \boxed{B(s)\ge d^{m+E}/J},
 \qquad
 E=\frac m2-\left\lceil\frac{m}{2J}\right\rceil-J.
\]

The estimate holds for every integer \(1\le J<m\), every even
register length \(m\ge2\), every prime-power digit alphabet \(d\),
every balanced coordinate partition, and every scale \(s\ge1\).
The coefficient field has the same characteristic as \(d\).
No linearity assumption on the reversible feedback is needed.

In particular,

\[
 \rho_{\mathrm{bare}}\ge\frac{17}{24}>\frac1{\sqrt2}
 \qquad\text{for every even }m\ge380.
\]

The former cutoff was \(576\). This is a stronger obstruction within
the same declared bare-profile method, not an improved factoring
exponent. Physical word types, gauges, other physical bases, and other
polynomial generator transforms remain outside this conclusion.

## 1. One map controls every early coefficient combination

Let

\[
 H(w)=\sum_{k\ge0}w^k H_k
\]

be a rectangular matrix over a field. Let its finite Smith orders over
the power-series ring be \(\mu_i\). For an integer \(j\ge1\), reduce
modulo \(w^j\). Multiplication by \(H\) is a linear map

\[
 \mathcal T_j:
 (\mathbb F[w]/w^j)^b\longrightarrow
 (\mathbb F[w]/w^j)^a.
\]

As an ordinary matrix over \(\mathbb F\), this is the lower triangular
block Toeplitz matrix with block \(H_{r-c}\) at position \((r,c)\).
Smith unit factors remain invertible after reduction, so

\[
 \mathrm{rank}_{\mathbb F}\mathcal T_j
   =L_j:=\sum_i\max(0,j-\mu_i). \tag{1}
\]

For arbitrary scalars \(c_0,\ldots,c_{j-1}\), let \(C\) embed a
constant input vector as a truncated polynomial, and let \(P_c\) take
the indicated linear combination of the output coefficients. Then

\[
 \sum_{k=0}^{j-1}c_kH_k=P_c\mathcal T_j C.
\]

Consequently,

\[
 \boxed{\quad
 \mathrm{rank}\left(\sum_{k<j}c_kH_k\right)
 \le L_j\le jN_j,
 \qquad N_j=\#\{i:\mu_i<j\}.
 \quad} \tag{2}
\]

This counts the entire linear combination through one map. Summing
separate coefficient-rank bounds loses an unnecessary factor of \(j\).

## 2. A weighted channel bound

For the regular full-cycle crossing, put \(n=d^m=q^2\). The normalized
reversed row matrix is

\[
 H(w)=\Phi\bigl(U(I-wS)^{-1}\bigr),\qquad S=I-U.
\]

The row polynomial has degree \(n-1\), and its column zero orders are
\((0,n,\ldots,n)\). Thus its relevant bare product is exactly

\[
 B(s)=\sum_i g_{\mu_i}(s),\qquad
 g_\mu(s)=s^{n-1-\mu}+(n-1)s^{-1-\mu}. \tag{3}
\]

Suppose \(0\le\mu<j\) and \(j^2\le n\). We claim that

\[
 g_\mu(s)\ge \frac n j(j-\mu)
 \qquad(s\ge1). \tag{4}
\]

Write \(a=\mu+1\). Since \(a\le j\le\sqrt n<n\), the unrestricted
positive minimum, after setting \(x=s^n\), is

\[
 \min_{s>0}g_\mu(s)=
 n\left(\frac{n-1}{n-a}\right)^{1-a/n}a^{-a/n}
 \ge n\exp\left(-\frac{a\ln a}{n}\right).
\]

The elementary inequalities \(e^{-x}\ge1-x\) and
\(\ln(1+\mu)\le\mu\) give

\[
 \begin{aligned}
 g_\mu(s)
 &\ge n-a\ln a\\
 &\ge n-\mu(\mu+1)\\
 &\ge n-\mu j\\
 &\ge n\left(1-\frac\mu j\right).
 \end{aligned}
\]

The last inequality uses \(j^2\le n\). The \(\mu=0\) case gives
\(g_0(s)\ge n\) directly as well. This proves (4).

Summing only over the early channels in (3), then applying (1)--(2),
proves the reusable bound

\[
 B(s)\ge\frac n j L_j
       \ge\frac n j\,
          \mathrm{rank}\left(\sum_{k<j}c_kH_k\right). \tag{5}
\]

The weighted sum \(L_j\) is essential here. Replacing every early
channel by the same worst-case weight would lose a second factor.

## 3. Register powers and copied coordinates

The binomial identity in the coefficient field gives

\[
 U^j=U(I-S)^{j-1},\qquad
 \Phi(U^j)=\sum_{k<j}(-1)^k\binom{j-1}{k}H_k.
\]

Hence (5) applies with this explicit coefficient combination:

\[
 B(s)\ge\frac n j\,\mathrm{rank}\Phi(U^j). \tag{6}
\]

For completeness, the copied-wire argument selects such a power as
follows. Colour register coordinates by their two physical halves, and
let \(c_j\) count cyclic agreements with shift \(j\). Expanding the
sum of squares of every length-\((J+1)\) colour window shows that some
\(1\le j\le J\) has

\[
 c_j\ge\frac m2-\frac{m}{2J}.
\]

After \(j<m\) register steps, at least \(c_j-j\) same-colour digit
wires are copied literally. In the row reshuffling, they create
\(d^{c_j-j}\) disjoint nonzero diagonal blocks. Therefore

\[
 \mathrm{rank}\Phi(U^j)\ge d^E,
 \qquad E=\frac m2-\left\lceil\frac{m}{2J}\right\rceil-J.
\]

Every allowed \(J<m\) satisfies \(J^2\le(m-1)^2\le2^m\le d^m=n\).
The middle inequality follows by direct checks at \(m=2,3,4\) and
induction thereafter, since \(m^2\le2(m-1)^2\) for \(m\ge4\).
Thus the hypothesis of (6) is automatic. Substituting the selected
power gives

\[
 B(s)\ge\frac{nd^E}{j}\ge\frac{d^{m+E}}J. \tag{7}
\]

The opposite bare orientation remains at least \(q^4\), by the
leading-coefficient and determinant-degree argument already retained
in the copied-wire theorem. Thus it cannot defeat the stated lower
bound on the optimized bare rate.

## 4. A complete cutoff at 380 digits

First treat the infinite tail. For \(m\ge576\), let
\(r=\sqrt m\ge24\) and \(J=\lfloor r\rfloor\). Since
\(J\ge r-1\),

\[
 \left\lceil\frac{m}{2J}\right\rceil+J
 \le\frac{3r}{2}+2.
\]

Equation (7), with \(d\ge2\), gives

\[
 \log_d B(s)\ge
 \frac{3m}{2}-\frac{3r}{2}-2-\log_2r
 \ge\frac{3m}{2}-2r
 \ge\frac{17m}{12}.
\]

The middle inequality uses \(2+\log_2r\le r/2\) for \(r\ge24\).
At \(24\), use \(\log_2 24<5\); the difference then increases,
since \(\ln2\ge1/2\). The last inequality is \(r\ge24\).

For every even \(380\le m<576\), the checker retains an explicit
integer \(J\) and verifies

\[
 2^{12E-5m}\ge J^{12}. \tag{8}
\]

These are \(98\) exact integer comparisons. The chosen \(J\) is the
smallest minimizer of \(J+\lceil m/(2J)\rceil\). Its minimum equals
\(\lceil\sqrt{2m}\rceil\), because an integer sum \(L\) permits
\(J(L-J)\ge m/2\) precisely when \(\lfloor L^2/4\rfloor\ge m/2\).
No floating-point optimization is used.

At the first size, \(m=380\), use \(J=12,E=162\):

\[
 2^{44}=17592186044416
 \;>\;12^{12}=8916100448256.
\]

The older choice \(J=\lfloor\sqrt{m/2}\rfloor=13\) fails that
particular comparison, so the explicit witness matters. For every
larger digit alphabet, the nonnegative exponent in (8) makes the
binary comparison sufficient.

Finally,

\[
 \rho_{\mathrm{bare}}=\frac{\log_d B}{4\log_dq}
 =\frac{\log_d B}{2m}\ge\frac{17}{24},
 \qquad 2\cdot17^2>24^2.
\]

The checker also verifies that \(m=378\) has no witness for (8) among
all \(1\le J<m\). This identifies the last failure of this particular
integer certificate, not the exact boundary of all possible profile
obstructions. The asymptotic conclusion remains
\(\liminf\rho_{\mathrm{bare}}\ge3/4\).

## 5. Exact controls

```sh
python3 tools/check_shift_register_toeplitz_budget.py \
  --output /tmp/shift-register-toeplitz-budget.json
```

The checker verifies truncated multiplication ranks for rectangular
matrices with known Smith orders over three prime fields. It replays
the explicit constant-input and output-projection identity, complete
register coefficient combinations, rational channel inequalities,
and every finite cutoff witness. These controls accompany the
all-field and all-size proofs above. No new factoring saving or Lean
formalization is claimed.
