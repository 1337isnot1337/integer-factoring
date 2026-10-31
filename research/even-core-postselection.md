# A polynomial crossing family from a fixed even core

## Result and scope

There are explicitly generated tensors

\[
A(w),B(w),C(w)\in
\mathbb F_2[w]^{4\times4\times4\times4}
\]

of respective degrees \(17,32,14\) such that

\[
K_\lambda(w)=A(w)+\lambda B(w)+\lambda^2C(w)
\tag{1}
\]

has the strict crossing endpoint

\[
K_\lambda(1)_{xyuv}=[x=v][y=u]
\tag{2}
\]

over every field of characteristic two, for every scalar parameter
\(\lambda\). The parameter-zero member is Gaussian. Every nonzero
rational-function specialization \(\lambda=\lambda(w)\) fails at least
one of two explicit Wick relations.

This is a constructive family for further profile and type searches.
It does not establish an improved factoring exponent.

## 1. Construction by fixing the extreme modes

Start with the four-mode one-species Gaussian \(G(w)\) on the interval
\(\{0,1,2,3\}\), with ambient width \(M=4\). Use the kernel and conventions
from [the interval construction](schur-one-species-interval.md). Set

\[
c=1+w,\qquad Q=x_0y_3u_0v_3,
\qquad \mathcal T=G\left(1+\frac{Q}{c^2}\right).
\tag{3}
\]

Keep modes \(1,2\) as the two occupation bits of each physical letter.
On the extreme modes \(0,3\), use the input vector

\[
|\varnothing\rangle+|\{0,3\}\rangle
\]

on both input legs. Use the output covector

\[
\langle\varnothing|+\lambda\langle\{0,3\}|
\]

on both output legs. Both vectors have even fermionic parity.
There are two output legs, so this contraction is quadratic in
\(\lambda\).

The numerator is obtained by contracting \(c^2G+QG\).
Grouping its terms by the number of both-occupied output extreme states
gives three binary polynomial tensors. Each is divisible by \(c^2\).
Dividing them gives \(A,B,C\) in (1).

The exact generator uses only minors of an \(8\times8\) polynomial
matrix. It performs 768 exact polynomial divisions. In the four-index
coefficient order, each letter is a two-bit occupation mask for modes
\(1,2\), and the dense index is
\(x+4y+16u+64v\).

### Why the denominator cancels

The energy argument in
[the residual-code construction](conjugate-residual-phase-codes.md)
applies before the extreme-mode contraction. Put \(w=1+t^3\) and multiply
each coefficient by its output-minus-input index energy.
The resulting Gaussian and quartic tensor is regular in \(t\).

Since the extreme core is always an even state, charge conservation fixes
the energy increment of every term contributing to one logical entry:

\[
E=\frac12\Delta(n_2-n_1),\qquad |E|\le2,
\tag{4}
\]

where \(\Delta\) means the sum on output legs minus the sum on input legs.
The value of \(E\) is independent of which vacuum/both components were
chosen. On physical support it is an integer.

If the contracted rational entry had a pole at \(c=0\), its order after
the substitution would be at most \(-3+E<0\), contrary to regularity.
There are no other possible finite poles. Thus the contraction is
polynomial in \(w\). The same argument works with \(\lambda\) treated as
an indeterminate. Equivalently, the three parameter coefficients divide
by \(c^2\) separately.

The exact endpoint identities are

\[
A(1)=E_{\rm cross},\qquad B(1)=C(1)=0.
\tag{5}
\]

The checker verifies all 768 coefficient identities. They prove the
strict local endpoint (2), before taking tensor powers or restricting to
word types.

## 2. Two exact obstructions to Gaussianity

For an even tensor \(K\), and a four-element set \(S=\{a,b,c,d\}\) of
occupation variables, define its quartic Wick residual by

\[
\mathcal W_S(K)=
 K_\varnothing K_{abcd}
\quad +K_{ab}K_{cd}+K_{ac}K_{bd}+K_{ad}K_{bc}.
\tag{6}
\]

A fermionic Gaussian tensor satisfies these Pfaffian relations.
All signs are plus in characteristic two.

In the occupation ordering

\[
(x_1,x_2,y_1,y_2,u_1,u_2,v_1,v_2),
\]

the family (1) has the following exact residuals:

\[
\begin{aligned}
\mathcal W_{\{x_1,x_2,y_1,u_1\}}(K_\lambda)
 &=\lambda^3w^2c^4(w^8+\lambda),\\
\mathcal W_{\{x_1,x_2,u_1,v_1\}}(K_\lambda)
 &=\lambda^3w^2c^7(1+\lambda).
\end{aligned}
\tag{7}
\]

These are masks 23 and 83 in the stated bit order.
The checker verifies (7) as polynomial identities in both \(w\) and
\(\lambda\), by comparing all five parameter coefficients.

For any nonzero \(\lambda(w)\) in a characteristic-two rational-function
field, simultaneous vanishing of the two residuals would force
\(\lambda=w^8\) and \(\lambda=1\). This is impossible because \(w\) is
indeterminate. Thus every nonzero rational-function specialization
violates a Wick relation.

In particular, the binary tensor \(K_1\) has degree 32 and

\[
\mathcal W_{\{x_1,x_2,y_1,u_1\}}(K_1)
=w^2(1+w)^{12}\ne0.
\tag{8}
\]

Its vacuum coefficient is \(w^{16}\). A complete check of the 99 even
Pfaffian identities beyond the two-variable coefficients finds 68
failures.

For \(\lambda=0\), the vacuum coefficient is one, every odd coefficient
vanishes, and all 99 Pfaffian identities hold. This gives an exact
Gaussian control for the same construction.

## 3. Parameter paths

The endpoint identities (5) allow any polynomial substitution
\(\lambda=\lambda(w)\). Laurent-polynomial substitutions are also
available after multiplication by a sufficiently large common power of
\(w\); that power equals one at the crossing endpoint.

For example, \(\lambda=w^k\) gives a finite family to screen with integer
Smith orders. A common monomial clearing factor shifts the two places
in opposite directions and does not change their product pressure.
The same cancellation does not apply to an arbitrary non-monomial
denominator. Such a denominator must be retained in any proposed
certificate.

The value \(\lambda=w^8\) cancels the first residual in (7), but the
second remains nonzero. It is not a Gaussian exception.

The construction therefore provides endpoint-preserving deformations
whose non-Gaussian behavior is certified independently of a numerical
optimizer. Improving the pressure remains a separate task.

## 4. Reproduction

Run

~~~sh
python tools/check_even_core_postselection.py \
  --output evidence/even-core-postselection.json
~~~

Add a seed directory to export \(A,B,C\), the two binary constant
specializations, and all four constant specializations over
\(\mathbb F_4\):

~~~sh
python tools/check_even_core_postselection.py \
  --seed-directory /tmp/even-core-seeds
~~~

The [checker](../tools/check_even_core_postselection.py) uses exact
polynomial arithmetic throughout. The
[evidence](../evidence/even-core-postselection.json) records coefficient
hashes, division and endpoint counts, the two symbolic Wick identities,
and the Gaussian control. No native search library or Lean formalization
is required.
