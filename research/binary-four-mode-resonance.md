# An exact four-mode resonance reduction

## Result and scope

At a power-of-two Pascal interval length, an occupation transfer between
two complementary, equal-energy two-particle states has an unusually
small exact profile. Only five scalar zero factors and one two-by-two
zero block change. At infinity, only the same two-by-two occupation
block changes. Every other selected state and every unselected Gaussian
mode remains explicit.

The reduction applies over any field of characteristic two, including
the fixed field \(\mathbb F_4\). It permits Laurent coefficient
matrices and an arbitrary invertible endpoint on the resonant block.
For the narrower family of polynomial matrices of degree at most
\(2nP\), neither cut can improve its Gaussian partition function at
any scale. That obstruction does not cover general Laurent matrices.
No improvement to the factoring exponent is asserted here.

## 1. The admissible resonant operator

Let \(n=2^e\), let \(P\) be a positive odd integer, and put

\[
q=z^2,\qquad m=nP,\qquad c=1-q^m.
\tag{1}
\]

Select four distinct identity indices \(j_1,j_2,j_3,j_4\) satisfying

\[
j_1+j_4=j_2+j_3.
\tag{2}
\]

Write \(A=\{1,4\}\) and \(B=\{2,3\}\) for the corresponding
binary occupation states. Thus \(A\) and \(B\) are complementary,
have occupation two, and have equal physical selected-mode energy.
Let \(K(q)\) be a two-by-two Laurent-polynomial matrix with
\(\det K(1)\ne0\). Define the full sixteen-state occupation kernel
\(F\) to be the identity except that its \(A,B\) block is \(K\).

Use \(F\) in the
[polynomial identity-operator construction](coupled-identity-operator-filters.md).
Its endpoint is admissible: the only nonidentity endpoint block
preserves both occupation and physical energy. Applying the constant
physical \(Y\)-label map \(K(1)^{-1}\) on this block restores the
crossing coefficient by coefficient. This repair is regular and
invertible at both weighted cuts and adds no degree.

All negative Laurent exponents are cleared by one common monomial.
The full ordinary degree is charged as in the operator construction;
the local units used below do not remove that cost. In particular,
for a clearing exponent \(h\) and
\(D_F=\max_{a,b}\deg(q^hF_{ab})\), its existing bound gives

\[
\deg_v(q^hT_F)
\le P\bigl(\deg_wG+8n\bigr)+D_F.
\tag{3}
\]

The degree bound is sufficient and need not be sharp for a particular
matrix. Common clearing shifts the two cuts by opposite scalar orders.

## 2. The exact reduced zero matrix

The [prime-power identity theorem](prime-power-joint-identity-filters.md)
separates all selected directions simultaneously. Its Weyl parameter
is \(C=cI_4\), independently of the locations of the four indices.
The [weighted Weyl reduction](identity-operator-weyl-reduction.md)
therefore gives, for binary occupation states \(a,b\),

\[
\Psi_C(E_{ab})
=c^{-(4-|b|)}\partial^{\mathbf1-a}Y^{\mathbf1-b}.
\tag{4}
\]

The multiplication operator acts before the derivatives. The variables
have square zero. In particular,

\[
\begin{aligned}
\Psi_C(E_{AA})&=c^{-2}\Pi_{\{s:s\subseteq A\}},&
\Psi_C(E_{BB})&=c^{-2}\Pi_{\{s:s\subseteq B\}},\\
\Psi_C(E_{AB})&=c^{-2}E_{AB},&
\Psi_C(E_{BA})&=c^{-2}E_{BA}.
\end{aligned}
\tag{5}
\]

For example, multiplication by \(Y_A\), followed by differentiation
in all variables of \(B\), is nonzero only on the input state \(B\).
Its output is exactly \(A\). A diagonal term instead projects onto
the subsets of its support. These observations prove (5) directly.

For the identity occupation kernel, the one-variable transformed
diagonal is \((q^m/c,1)\). Hence \(\Psi_C(I)\) on a state of
occupation \(d\) is \((q^m/c)^{4-d}\). Put

\[
\delta_A=K_{AA}-1,\qquad \delta_B=K_{BB}-1.
\tag{6}
\]

On each occupation-\(d\) block, multiply by the zero unit
\(c^{4-d}\). The resulting exact reduced blocks are:

| State or block | Reduced coefficient after the unit scaling |
|---|---|
| Vacuum | \(q^{4m}+c^2(\delta_A+\delta_B)\) |
| Singleton \(\{i\}\), \(i\in A\) | \(q^{3m}+c\delta_A\) |
| Singleton \(\{i\}\), \(i\in B\) | \(q^{3m}+c\delta_B\) |
| Resonant block \(A,B\) | \(R(q)=(q^{2m}-1)I_2+K(q)\) |
| Every other state \(s\) | \(q^{m(4-|s|)}\) |

There are two singleton states of each kind. Thus the table changes
five scalar entries and one two-by-two block in the sixteen-state
matrix. It retains all sixteen states.

For precise height accounting, define

\[
w(s)=P\sum_{i=1}^4j_i s_i,\qquad
a_i=P(2j_i+1-J).
\tag{7}
\]

Multiply each row indexed by \(s\) in the table by \(q^{w(s)}\).
Calculate its Smith orders with respect to \(q\), separately in
each occupation block; write one such order as \(\sigma\), and
write \(d\) for that block's occupation. The complete active zero
orders with respect to \(z\) are exactly

\[
\boxed{
2\sigma+P(1-J)d-\sum_{i=1}^4a_i\varepsilon_i,
\qquad \varepsilon\in\{0,1\}^4.
}
\tag{8}
\]

Each of the sixteen reduced orders supplies all sixteen auxiliary
shifts in (8), giving the full \(256\)-order active zero profile.
Equation (2) gives \(w(A)=w(B)\), so the resonant block requires
only the ordinary two-by-two Smith orders of \(R\), followed by
the common shift \(w(A)\). All unselected zero spectators are the
unchanged Gaussian modes from the prime-power theorem.

## 3. The exact infinity matrix

Put \(\kappa_i=P(2j_i+1-2n)\). The general infinity reduction
weights the column state \(s\) by
\(z^{\sum_i\kappa_i s_i}\). By (2), the two resonant weights
coincide:

\[
\kappa_A=\sum_{i\in A}\kappa_i
        =\sum_{i\in B}\kappa_i=\kappa_B.
\tag{9}
\]

Thus the changed infinity orders are precisely

\[
\boxed{
\kappa_A+\operatorname{Smith}_z K(z^{-2}).
}
\tag{10}
\]

The other fourteen selected occupation orders and every unselected
infinity spectator remain Gaussian. No series truncation is needed
to evaluate either changed two-by-two block: its first order is the
minimum entry valuation, and the sum of its orders is the valuation
of its nonzero determinant.

## 4. A polynomial-block obstruction

Suppose now that every entry of \(K\) is a polynomial of degree at
most \(2m\). At any fixed scale \(s>1\) and aggregate gauge \(J\),
the complete zero partition function is at least that of the same
Gaussian, and the same is true separately at infinity.

For the vacuum, the coefficient in the table has degree at most
\(4m\). It evaluates to one at \(q=1\), so it is nonzero and
its order at zero is at most \(4m\), the Gaussian order.
Each singleton coefficient similarly has degree at most \(3m\),
evaluates to one at the endpoint, and has order at most \(3m\).
All the other scalar orders are unchanged.

For the resonant block, every entry of \(R\) has degree at most
\(2m\), while \(\det R(1)=\det K(1)\ne0\). If its two
\(q\)-Smith orders are \(\lambda_1,\lambda_2\), then

\[
\lambda_1+\lambda_2
=\operatorname{ord}_{q=0}\det R
\le\deg\det R\le4m.
\tag{11}
\]

The equal row weights in (7), the common occupation, and each fixed
auxiliary state in (8) add a common shift to this pair. Convexity of
the exponential therefore gives

\[
s^{-2\lambda_1}+s^{-2\lambda_2}
\ge2s^{-(\lambda_1+\lambda_2)}
\ge2s^{-4m}.
\tag{12}
\]

The right side is the Gaussian pair with both \(q\)-orders equal
to \(2m\). This proves the zero comparison, including every
auxiliary shift and spectator.

At infinity, if \(\eta_1,\eta_2\) are the orders of
\(K(z^{-2})\), polynomiality gives

\[
\eta_1+\eta_2=-2\deg_q\det K\le0.
\tag{13}
\]

Consequently \(s^{-\eta_1}+s^{-\eta_2}\ge2\). The common
shift in (10) and all unchanged states preserve this inequality.
The physical type partition is the same as for the Gaussian, so
this family cannot improve the pressure at any fixed \(s,J\).
This is a comparison with the same charged Gaussian, not an
assumption about the optimal Gaussian gauge.

The argument includes constant and degree-one polynomial resonant
blocks. It does **not** include a general Laurent block such as
\(K=E+(1+q^{-1})N\). Clearing its negative exponents also multiplies
the other occupation blocks; one cannot apply the isolated-block
polynomial argument while leaving those blocks equal to the identity.
Higher-degree polynomial blocks and changes coupled across different
occupation levels also lie outside this obstruction.

## Verification status

The formulas follow from the proved prime-power and weighted Weyl
identities cited above. The exact implementation is
[`resonant_identity_operators.py`](../tools/resonant_identity_operators.py).
Its independent replay is
[`check_resonant_identity_operators.py`](../tools/check_resonant_identity_operators.py),
with saved evidence in
[`resonant-identity-operators.json`](../evidence/resonant-identity-operators.json).

The saved run contains 41 complete two-cut comparisons with the generic
sixteen-state Weyl engine, 120 full-pressure comparisons, and an original
physical control of both \(256\)-by-\(256\) cuts and all 8,027
endpoint coefficients. Finite search scopes belong to their separate
search records. No finite search result is needed for the obstruction
in Section 4.
