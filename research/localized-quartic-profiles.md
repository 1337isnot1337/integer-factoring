# Localized quartic interactions: determinant gains and exact zero profiles

## Status

This note gives exact formulas for one family of non-Gaussian perturbations of
the [one-species interval crossing](schur-one-species-interval.md). It separates
a determinant improvement from the complete Smith pressure. No improved
factoring exponent is claimed.

The formulas cover Laurent coefficients and coefficients with a pole at the
crossing point. A pole at the crossing point requires a separate endpoint
argument; the [positive-energy endpoint lemma](positive-energy-quartic-endpoints.md)
provides one such argument. The final section records why a limiting
rank-deficient endpoint is not already a crossing on its natural alphabet.

## 1. The interval Gaussian and its inverse

Work over a field of characteristic two. Let \(M\ge2\) be a power of two,
\(z=w^2\), \(q=1+z\), and \(c=1+w\). For indices in \([0,M)\), put

\[
 P_{ab}=[b\subseteq a],\qquad
 Q_{aa}=q^a,\qquad H_{ab}=z^{M-1-(a\mathbin\vert b)}.
\]

Binary inclusion is denoted by \(\subseteq\). The one-particle matrix of the
one-species tensor, with rows \((v,y)\) and columns \((u,x)\), is

\[
 D=\begin{pmatrix}A&z^MP\\I&B\end{pmatrix},
 \qquad A=cQP^T,\qquad B=cH.
 \tag{1}
\]

Its scalar tensor \(G\) consists of all minors of \(D\). In particular,
the coefficient of the physical occupancy tuple \((x,y,u,v)\) is the minor
with rows \((v,y)\) and columns \((u,x)\).

The full interval identities \(P^2=I\) and
\(AB=(1+z^M)P\) give the exact rational inverse

\[
 D^{-1}=\begin{pmatrix}BP&z^MI\\P&PA\end{pmatrix}.
 \tag{2}
\]

Two entries needed below have particularly simple forms:

\[
 (BP)_{cd}=[c\subseteq d]c\,q^{M-1-d},\qquad
 (PA)_{ab}=c\,z^{a\mathbin\&b}.
 \tag{3}
\]

These identities follow either by multiplying the binary digit matrices or
by direct polynomial multiplication.

## 2. A two-creation, two-annihilation determinant identity

Let \(a_1,a_2\) be independent creation operators on an exterior algebra
with \(r\) one-particle modes. Let \(b_1,b_2\) be annihilation operators, and
write \(C_{ij}=\{b_i,a_j\}\). For

\[
 N=a_1a_2b_1b_2,\qquad \mu=\det C,
\]

the anticommutation relations give

\[
 N^2=\mu N.
 \tag{4}
\]

If \(\mu\ne0\), the image is exactly the subspace divisible by
\(a_1a_2\), which has dimension \(2^{r-2}\). Indeed, the displayed operator
acts there as multiplication by \(\mu\). Therefore

\[
 \det(I+fN)=(1+f\mu)^{2^{r-2}}.
 \tag{5}
\]

If \(\mu=0\), (4) makes \(N\) nilpotent and the determinant is one.
Thus (5) also holds in this case. This is an identity over the rational
function field; it does not assert anything about the individual Smith
orders.

Now multiply \(G\), as an exterior polynomial, by

\[
 1+f(w)\mathcal Q,\qquad
 \mathcal Q=x_a y_b u_c v_d,\qquad d=M-1,\quad a<d.
 \tag{6}
\]

For the column flattening, the contraction matrix is the submatrix of (2)
with rows \((u_c,x_a)\) and columns \((v_d,y_b)\). Equations (2)–(3) imply

\[
 \mu_C=(BP)_{cd}(PA)_{ab}+z^M[c=b]P_{ad}
       =q\,w^{2(a\mathbin\&b)}.
 \tag{7}
\]

For the row flattening, the contraction is
\(\mu_R=z^{-M}[b=c]P_{ad}=0\). Consequently its determinant is unchanged,
while the column determinant is multiplied by

\[
 \left(1+f(w)q\,w^{2(a\mathbin\&b)}\right)^{2^{2M-2}}.
 \tag{8}
\]

For disjoint \(a,b\), several illustrative choices are:

| Coefficient \(f\) | Ratio \(1+fq\) | Column determinant-order gain at infinity, divided by \(2^{2M}\) |
|---|---|---:|
| \((1+w)w^{-3}\) | \(w^{-1}+w^{-2}+w^{-3}\) | \(1/4\) |
| \(qw^{-4}\) | \(w^{-4}\) | \(1\) |
| \(w^{-2}\) | \(w^{-2}\) | \(1/2\) |
| \(1/[w(1+w)]\) | \(w^{-1}\) | \(1/4\) |
| \(q^{-1}\) | \(0\) | The column flattening loses one quarter of its generic rank. |

The last three rows do not have a locally unchanged endpoint. Endpoint
validity must be established after the physical restriction.

## 3. Exact zero-place profiles

Use the physical mode heights

\[
 h_x(a)=2a,\quad h_y(a)=2a+2-2M,\quad
 h_u(a)=2M-1,\quad h_v(a)=1-2M.
 \tag{9}
\]

The row and column internal Gaussian factors have regular coefficients:

\[
 x_a y_b:\ c w^{2(a\mathbin\&b)},\qquad
 u_c v_d:\ c q^d[d\subseteq c].
\]

They are units in their exterior algebras. Removing them is a regular
unimodular transformation on each flattened side. These factors commute
with the quartic multiplier. The remaining cross matrix is

\[
 \mathrm{diag}(w^{2a+1})P^T\ \oplus\
 \mathrm{diag}(w^{2a+1}).
\]

The constant Pascal change of the \(v\) variables makes it diagonal. Its
inverse fixes the top-index variable \(v_{M-1}\), because the last column
of \(P\) has only its diagonal entry. This detail matters: a quartic using
an arbitrary \(v_d\) generally spreads into several directions under the
same change.

Label the resulting two copies of the mode set by \(X,Y\), with energies
\(e_{X,j}=e_{Y,j}=2j+1\). Put

\[
 I=\{(X,a),(Y,b)\},\quad J=\{(X,M-1),(Y,c)\},
\]

\[
 A_0=2(a+b+1),\quad B_0=2(M+c),\quad
 r_0=\mathrm{ord}_0(f)+2(a+b+1-M).
 \tag{10}
\]

Since \(a<M-1\), the two subsets \(I,J\) differ, even when \(b=c\).
For every spectator subset \(T\) outside \(I\cup J\), the perturbation
changes exactly the entry from \(J\cup T\) to \(I\cup T\). The associated
two-dimensional block has the form

\[
 w^{e(T)}
 \begin{pmatrix}w^{A_0}&w^{r_0}u(w)\\0&w^{B_0}\end{pmatrix},
 \qquad u(0)\ne0.
\]

Its Smith orders are

\[
 e(T)+m_0,\quad e(T)+A_0+B_0-m_0,\qquad
 m_0=\min(A_0,B_0,r_0).
 \tag{11}
\]

All other subset-sum orders are unchanged. This proves the full
\(2^{2M}\)-entry profile without constructing that matrix.

For \(s>1\), let \(U_0(s)\) be the sum of \(s^{-\nu}\) over the
zero-place orders. The unperturbed sum is
\(U_0^{\rm base}=\prod_{j=0}^{M-1}(1+s^{-(2j+1)})^2\). Equation (11) gives

\[
 \frac{U_0(s)}{U_0^{\rm base}(s)}=
 1+\frac{s^{-m_0}+s^{-(A_0+B_0-m_0)}-s^{-A_0}-s^{-B_0}}
          {\prod_{i\in I\cup J}(1+s^{-e_i})}.
 \tag{12}
\]

The numerator is nonnegative, and strictly positive when
\(m_0<\min(A_0,B_0)\). Thus a determinant gain in (8) comes with a
precisely measurable zero-place spread in this family.

For the complementary choice \(a+b=M-1\), a coefficient of order
\(-D\) replaces the active orders \(2M,2M+2c\) by
\(-D,4M+2c+D\). Conversely, the low-index choice \(a=b=0\) does not make
this trade inexpensive: with \(f=w^{-2}\), its new pivot is \(-2M\).
More generally, the determinant-improving monomial
\(f=w^{-2(a\mathbin\&b)-2}\) has

\[
 r_0=-2\bigl(M-(a\mathbin\vert b)\bigr).
\]

Keeping this pole bounded requires the union of the two indices to lie
near the top of the interval. The active row energy is then also large.
This is a constraint on this localized family, not a general obstruction
to non-Gaussian constructions or other physical gauges.

## 4. The rank-deficient rational endpoint

Take \(a=c=0\), \(b=d=T=M-1\), and \(f=q^{-1}\). The quartic increases
the index sum by \(2T\) and has a double pole at \(w=1\). To compensate
one quartic, a term in a common fixed-mode type needs at least two
Gaussian internal edges. Equality is possible only for the two edges

\[
 u_Tv_0,\qquad x_Ty_0,
\]

each of which decreases the index sum by \(T\). All remaining cross edges
must have equal indices. Thus the endpoint after the global type
restriction is the corresponding type restriction of the following
four-letter tensor on the two extreme modes:

\[
 E=(1+x_0v_0)(1+x_Tv_T)(1+u_0y_0)(1+u_Ty_T)
   (1+u_Tv_0)(1+x_Ty_0)(1+x_0y_Tu_0v_T).
 \tag{13}
\]

The unaffected modes retain their crossing factors. Formula (13) is a
description of the residual endpoint, not a crossing claim. Its column
flattening has rank 12 on a 16-dimensional flattened space.

An equivalent local construction makes this endpoint literal. Introduce
\(t\) by \(1+w=t^T\), and multiply each tensor coefficient by
\(t^{\text{output index sum}-\text{input index sum}}\). Every Gaussian
edge and the quartic factor is then regular in \(t\). Its value at
\(t=0\) is (13), with ordinary crossing factors on the other modes.
On a common fixed-mode type the additional monomial gauge cancels.

Restricting each physical letter to one of the singleton states
\(\{0\},\{T\}\), encode these as bits 0 and 1. Then

\[
 E(x,y,u,v)=[x+y=u+v],
 \tag{14}
\]

where the sum is an integer sum. Every two distinct binary words of the
same weight give the unwanted endpoint entry
\(E^{\otimes m}(x,y,x,y)=1\). Already the words \((0,1)\) and
\((1,0)\) supply a two-copy negative control. Therefore merely removing
the both-occupied physical letter does not repair the endpoint. Further
linear restrictions or a separate endpoint construction are required.

## 5. Exact checks

Run:

```sh
python tools/check_localized_quartic_profiles.py --output evidence/localized-quartic-profiles.json
```

The checker verifies the polynomial inverse and all relevant minors at
widths 2, 4, and 8; tests the creation-annihilation identity on 256 exact
four-mode operator pairs; independently computes 20 raw two-mode
zero-place Smith profiles, with determinant-order checks at infinity
when the tensor remains full rank; and checks selected profiles by
block-Toeplitz ranks. It also independently replays the complete
256 by 256 raw zero-place matrices for \(M=4,a=1,b=2,c\in\{0,1\}\)
with coefficient \((1+w)w^{-3}\), constructs (13) literally, and verifies
its rank and its off-crossing common-type entry.

The evidence includes the predicted complete four-mode zero profiles for
five coefficient valuations and four choices of \(c\). An independent
native calculation during discovery matched every order for the original
\((1+w)w^{-3}\) perturbation at all four values of \(c\). That native
discovery calculation is not a dependency of the exact checker.
