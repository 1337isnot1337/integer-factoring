# Every one-mode binary Laurent filter is dominated at zero charge

## Theorem and scope

Consider one selected identity direction in the
[prime-power Gaussian construction](prime-power-joint-identity-filters.md),
at aggregate charge \(J=0\). Allow arbitrary finite Laurent
polynomials over any fixed characteristic-two extension field, with
arbitrary nonzero endpoint values and the free diagonal endpoint repair.
Then its exact two-cut cost is at least the Gaussian cost at every
scale. Equality holds only when both occupation weights are the same
scalar Laurent monomial.

This covers arbitrary degree and arbitrary multilevel cancellations
in a single binary direction. It is stronger than a finite coefficient
screen or a monomial-filter exclusion. It does not assert a joint-mode
tensorization, a nonzero-charge comparison, or a theorem for other
digit counts. Together with the natural Gaussian pressure floor, it
rules out a strict critical witness in this entire one-mode family.

## 1. Exact scalar profile data

Let \(n\) be a power of two, let \(P\) be a positive odd integer,
and put \(m=nP\). Select any \(0\le j<n\), and write
\(\lambda=2n-1-2j>0\). Reverse the occupation indices as
\(f_0=\alpha_1\), \(f_1=\alpha_0\). The complete zero-place
transform is

\[
\beta_0=f_0,\qquad
\beta_1=f_1+(1+x^m)f_0.
\tag{1}
\]

Assume \(f_0(1),f_1(1)\ne0\). Both beta polynomials are then
nonzero as well, since \(\beta_i(1)=f_i(1)\). Define

\[
b_i=\mathrm{ord}_x\beta_i,\qquad D_i=\deg f_i.
\]

Orders and degrees may be negative. They are the exact finite Laurent
orders and complete degrees; no coefficient or denominator is omitted.
The [joint expectation formula](joint-binary-filter-synthesis.md),
specialized to one mode and \(J=0\), gives

\[
R=(\pi_0U_0+\pi_1U_1)(\pi_0V_0+\pi_1V_1),
\tag{2}
\]

where

\[
\begin{gathered}
q=e^{2t/P}>1,\qquad
U_0=q^{-b_0},\quad U_1=q^{m-b_1},\quad V_i=q^{D_i},\\
\pi_0=\frac1{1+e^{-\lambda t}},\qquad
\pi_1=\frac{e^{-\lambda t}}{1+e^{-\lambda t}},\qquad
\pi_0>\pi_1>0.
\end{gathered}
\tag{3}
\]

Here \(R\) denotes the ratio to the Gaussian cost at the same scale.
A common Laurent clearing shifts every \(D_i,b_i\) equally and
cancels between the two factors. The endpoint repair is a constant
diagonal map regular at both cuts, so it also leaves (2) unchanged.

## 2. The complete degree/order dichotomy

Always \(D_0\ge b_0\), so \(U_0V_0\ge1\).

If \(b_1\le D_1\), then \(U_1V_1\ge q^m>1\). Cauchy--Schwarz
immediately gives

\[
R\ge(\pi_0\sqrt{U_0V_0}+\pi_1\sqrt{U_1V_1})^2>1.
\tag{4}
\]

Otherwise \(b_1>D_1\). In (1), cancellation of the lowest term of
\((1+x^m)f_0\), whose order is \(b_0\), forces
\(\mathrm{ord}_x f_1=b_0\). Also \(b_1>D_1\) implies
\(b_1\le\deg((1+x^m)f_0)=D_0+m\). Consequently

\[
U_0V_1\ge1,\qquad U_1V_0\ge1,
\tag{5}
\]

in addition to \(U_0V_0\ge1\). To justify the first cancellation
claim explicitly: if the two input orders differed, their sum would
have the smaller order, which is at most \(D_1\); this contradicts
\(b_1>D_1\).

Put \(z=U_1/U_0\). Equations (5) give

\[
R\ge(\pi_0+\pi_1z)
       [\pi_0\max(1,1/z)+\pi_1].
\tag{6}
\]

For \(z\ge1\), the right side is \(\pi_0+\pi_1z\ge1\).
For \(0<z\le1\), it is
\(\pi_0^2/z+2\pi_0\pi_1+\pi_1^2z\). Its derivative is
\(-\pi_0^2/z^2+\pi_1^2<0\), so it is at least its value one at
\(z=1\). This proves \(R\ge1\) in every case.

## 3. Equality and verification

Equality in the second case requires \(z=1\),
\(D_0=D_1=b_0\), and \(b_1=b_0+m\). Thus
\(f_0=A x^{b_0}\) and \(f_1=B x^{b_0}\). Equation (1) has order
\(b_0+m\) only if \(A=B\). Conversely this common scalar monomial
has \(U_0=U_1=q^{-b_0}\), \(V_0=V_1=q^{b_0}\), and \(R=1\).

The mathematical argument is uniform in the field, degree, and
interval length. The independent finite replay checks the complete
polynomial transform, every order/degree implication, and the equality
classification over binary and quaternary coefficient catalogues.
It also evaluates representative full cost ratios with exact rational
arithmetic for independently varied biased probabilities.

```sh
python tools/check_one_mode_binary_laurent_filters.py
```

The complete counts and exact example catalogue are in
[`one-mode-binary-laurent-filters.json`](../evidence/one-mode-binary-laurent-filters.json).
The proof applies to the physical probabilities (3), while the rational
controls deliberately verify the stronger algebraic statement for any
\(\pi_0>\pi_1>0\) with \(\pi_0+\pi_1=1\).
