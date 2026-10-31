# A fixed-size zero profile for the qutrit identity summand

## Result and scope

The first identity-block summand of a truncated symmetric Pascal
interval has an exact two-mode active space at zero, for every interval
length. The basis changes are regular at zero, even when the physical
charge gauge is unbalanced. The other modes remain independent Gaussian
spectators.

For the quadratic perturbation alone, every active Smith order has a
closed form. In particular, the apparently favorable infinity change
has a substantial zero cost: the highest zero Gaussian mode disappears
from the active spectrum. More general qutrit digit filters still reduce
to five matrices of orders at most three. This note proves the reduction;
the [identity-filter implementation](../tools/truncated_symmetric_identity_filters.py)
uses the resulting exact finite matrices. No improved factoring exponent
is claimed here.

## 1. The polynomial summand and its zero coordinates

Work over a field of characteristic three. Take the interval
\(S=\{0,\ldots,n-1\}\), with ambient parameter \(M=2n\), in the
[truncated symmetric Pascal construction](truncated-symmetric-pascal-crossings.md).
Write its one-particle matrix as

~~~math
\mathcal D(w)=\begin{pmatrix}A&B\\ I&C\end{pmatrix},
\qquad
A=\operatorname{diag}((1+w)^{i+1})P^{\mathsf T},
\quad B=w^{2n}R,
\quad R=P^{-1}.
\tag{1}
~~~

The physical row variables are \((V,Y)\), and the column variables
are \((U,X)\). Every variable has cube zero. The polynomial
identity-block summand with index zero is

~~~math
Q_0=\left(Y_0+\sum_i A_{i0}V_i\right)
    \left(U_0+\sum_i C_{0i}X_i\right).
\tag{2}
~~~

Here \(A_{i0}=(1+w)[i=0]\) and
\(C_{0i}=(-1)^i w^{2n-1-i}\). Thus (2) has no endpoint pole.

Set \(w=v^P\), \(v=-z^2\), where \(P\) is a positive odd integer.
Use the following zero heights, with rational \(J\):

~~~math
h(V_i)=PJ,\quad h(Y_i)=P(2i+1),\quad
h(U_i)=-PJ,\quad h(X_i)=P(2i-4n+1).
\tag{3}
~~~

A further fixed ramification makes all exponents integral if necessary.
Equivalently, work over a fixed Puiseux extension and normalize its
valuation by \(\operatorname{ord}(z)=1\). Ramification multiplies
all integer-valued orders by the same positive constant, so every
regularity comparison and displayed normalized order is unchanged.
The internal Gaussian factors are regular under (3): the changes on
\(V\) and \(U\) cancel, while the \(Y,X\) heights are unchanged.
Their constant terms in the nilpotent variables are one, so they can
be removed by invertible regular multiplication. Multiplication by
these factors commutes with a multiplier polynomial in the same
variables, including (2).

Make the constant substitution \(V'=R^{\mathsf T}V\), and then omit
the prime. The Gaussian cross terms become diagonal, with orders

~~~math
k_i^-=P(2i+1-J),\qquad k_i^+=P(2i+1+J).
\tag{4}
~~~

Since the first row of \((R^{\mathsf T})^{-1}=P^{\mathsf T}\)
consists of ones, (2) becomes

~~~math
Q_0=\left(z^aY_0+(1-z^{2P})\sum_{i=0}^{n-1}V_i\right)
    \left(U_0-z^{-a}\sum_{i=0}^{n-1}X_i\right),
\qquad a=P(1-J).
\tag{5}
~~~

The following regular concentration is the step that prevents the
active space from growing with \(n\).

## 2. Concentrating both linear forms by regular transformations

Let \(t=n-1\), \(k_i=k_i^+\), and

~~~math
D=\operatorname{diag}(z^{k_i}),\qquad
u(z)=\sum_{h=0}^{n-1}z^{2Ph}.
\tag{6}
~~~

The constant term of \(u\) is one. The Gaussian in these variables
is \(X^{\mathsf T}DV\). A pair of substitutions
\(X_{\rm old}=L^{\mathsf T}X_{\rm new}\) and
\(V_{\rm old}=RV_{\rm new}\), with
\(R=D^{-1}L^{-1}D\), preserves this bilinear form.

First use

~~~math
L_1=I-\sum_{i<t}E_{i,t},\qquad
R_1=I+\sum_{i<t}z^{k_t-k_i}E_{i,t}.
\tag{7}
~~~

These matrices and their inverses are regular, because
\(k_t-k_i=2P(t-i)>0\). The vector of coefficients of
\(\sum X_i\) becomes \(e_t\). To track the other form, put
\(h=D^{-1}\mathbf 1\). The vector \(L_1^{-\mathsf T}h\)
has coordinates \(z^{-k_i}\) for \(i<t\) and final coordinate

~~~math
\gamma=\sum_{i=0}^{n-1}z^{-k_i}=z^{-k_t}u(z).
\tag{8}
~~~

Next use

~~~math
L_2=I+\sum_{i<t}\frac{z^{k_t-k_i}}{u(z)}E_{t,i},\qquad
R_2=I-\sum_{i<t}\frac1{u(z)}E_{t,i}.
\tag{9}
~~~

Again all four matrices, including their inverses, are regular.
The first matrix fixes \(e_t\). Its inverse transpose sends the
vector in (8) to \(\gamma e_t\): in coordinate \(i<t\), the
subtracted term is exactly
\(z^{k_t-k_i}u^{-1}\gamma=z^{-k_i}\).

Taking \(L=L_2L_1\) and \(R=R_1R_2\) therefore gives the exact
linear-form identities

~~~math
\sum_iX_i\longmapsto X_t,\qquad
\sum_iV_i\longmapsto u(z)V_t.
\tag{10}
~~~

These are identities over the local field with regular invertible
basis changes. They do not require each diagonal entry of \(D\)
to be regular; only the positive differences \(k_t-k_i\) are used.
The truncated symmetric functor carries these linear substitutions
to regular invertible maps on the coefficient matrices.

Since \((1-z^{2P})u(z)=1-z^{2Pn}\), the exact active tensor is

~~~math
\begin{aligned}
G_{\rm act}&=\mathrm{Exp}_{<3}(z^aYU)\,
             \mathrm{Exp}_{<3}(z^bXV),\\
Q_{\rm act}&=(z^aY+(1-z^{a+b})V)(U-z^{-a}X),\\
a&=P(1-J),\qquad b=P(2n-1+J),\qquad a+b=2Pn>0.
\end{aligned}
\tag{11}
~~~

The remaining zero modes are independent spectators with orders

~~~math
P(2i+1-J),\quad 1\le i<n;\qquad
P(2i+1+J),\quad 0\le i<n-1.
\tag{12}
~~~

Each spectator contributes the three digit orders \(0,k,2k\).
No regularity at infinity of (7) or (9) is asserted or needed: this
is a zero-cut proof.

## 3. Complete zero spectrum for a pure quadratic multiplier

Let \(F(z)\) have valuation \(-h<0\), and consider
\((1+FQ_{\rm act}^2)G_{\rm act}\). Write
\(A=z^a\), \(B=z^b\), \(T=AB\), and \(c=1-T\).
The zero cut uses rows \((Y,X)\) and columns \((U,V)\).
It preserves the charge
\(q=\deg Y-\deg X=\deg U-\deg V\), with block sizes
\(1,2,3,2,1\) for \(q=-2,-1,0,1,2\).

After extracting \(A^q\), the raw Gaussian coefficient blocks are

~~~math
\begin{array}{c|c}
q& A^{-q}M_q\\ \hline
2&2+F\\[2pt]
1&\begin{pmatrix}1&2Fc\\F&2T-F(1+T)\end{pmatrix}\\[8pt]
0&\begin{pmatrix}
1&0&Fc^2\\
0&T+2Fc&Fc(1+T)\\
F&-F(1+T)&T^2+2F(1+T)^2
\end{pmatrix}\\[12pt]
-1&\begin{pmatrix}T&Fc^2\\2Fc&2T^2+Fc(1+T)\end{pmatrix}\\[8pt]
-2&2T^2+Fc^2.
\end{array}
\tag{13}
~~~

Within a block, the smaller total degree comes first. Input-factorial
normalizations change these matrices only by constant diagonal units.

Divide each block by \(A^qF\). Since \(a+b>0\), we have
\(T\equiv0\), \(c\equiv1\), and \(F^{-1}\equiv0\) modulo
the maximal ideal. The two outer blocks reduce to one. The other
three reduce respectively to

~~~math
\begin{pmatrix}0&2\\1&2\end{pmatrix},\qquad
\begin{pmatrix}0&0&1\\0&2&1\\1&2&2\end{pmatrix},\qquad
\begin{pmatrix}0&1\\2&1\end{pmatrix}.
\tag{14}
~~~

All three determinants are one in characteristic three. Every
normalized block is therefore a regular invertible matrix. This proves
the complete active spectrum

~~~math
\boxed{\quad
\{\,a(r-s)-h:0\le r,s<3\,\}.
\quad}
\tag{15}
~~~

For the originally proposed slow multiplier
\(F=(1+v)/v=1-z^{-2}\), take \(h=2\). At the balanced gauge
\(J=1\), all nine active orders are \(-2\). In particular, the
positive Gaussian order \(b=2nP\) no longer contributes to the
active profile. The full spectrum is the sum of (15) and all digit
sums of the spectators (12).

## 4. General digit filters remain finite

For digit weights \(\alpha_0=1\), \(\alpha_1=r\),
\(\alpha_2=s\), the identity-mode multiplier is

~~~math
1+(r-1)Q_{\rm act}+(2s-r+2)Q_{\rm act}^2.
\tag{16}
~~~

The same proof of (11) applies. In fact, all five charge blocks admit
explicit polynomial unimodular factorizations. Put

~~~math
S=s,\qquad D=cs-r,\qquad F=1+cr+c^2s,\qquad c=1-T.
\tag{17}
~~~

After extracting the charge scalar \(A^q\), the two outer blocks
are \(2S\) and \(2F\). For charge \(+1\), the raw coefficient
matrix is

~~~math
\begin{pmatrix}-1&c\\-1&c-1\end{pmatrix}
\begin{pmatrix}D&0\\0&S\end{pmatrix}
\begin{pmatrix}1&-c\\1&1-c\end{pmatrix}.
\tag{18}
~~~

Both outer determinants are one. For charge \(-1\), the corresponding
factorization is

~~~math
\begin{pmatrix}1&-c\\1&1-c\end{pmatrix}
\begin{pmatrix}F&0\\0&D\end{pmatrix}
\begin{pmatrix}1&-c\\1&1-c\end{pmatrix}.
\tag{19}
~~~

For charge zero, write the matrix as
\(L\operatorname{diag}(F,D,S)R\), where

~~~math
L=\begin{pmatrix}
1&c&c^2\\
1&c+1&c^2-c\\
-1&1-c&-c^2-c-1
\end{pmatrix},\qquad
R=\begin{pmatrix}
1&-c&-c^2\\
1&-c-1&c-c^2\\
1&1-c&-c^2-c-1
\end{pmatrix}.
\tag{20}
~~~

These identities follow by multiplying the displayed matrices in
characteristic three. Subtracting \(c\) times the first column from
the second column of \(L\), and \(c^2\) times the first from the
third, gives \(\det L=-1\). Subtracting the first row of \(R\)
from each remaining row gives \(\det R=1\). Thus the factors are
unimodular over the polynomial ring in \(c\), without division by
\(c\). Since \(c=1-z^{a+b}\) is regular, they preserve local
Smith orders.

Let \(\sigma=\operatorname{ord}(S)\),
\(\delta=\operatorname{ord}(D)\), and
\(\phi=\operatorname{ord}(F)\). The complete active profile is

~~~math
\boxed{\quad
\begin{array}{c|c}
q&\text{Smith orders in the charge block}\\ \hline
2&2a+\sigma\\
1&a+\delta,\ a+\sigma\\
0&\phi,\ \delta,\ \sigma\\
-1&-a+\phi,\ -a+\delta\\
-2&-2a+\phi.
\end{array}
\quad}
\tag{21}
~~~

An identically zero scalar contributes an infinite order and a zero
direction. For the admissible monomial digit filters the three
polynomials can be evaluated directly. Equal-exponent cancellation
must be performed over \(\mathbb F_3\), including on walls where
multiple terms tie. There is no series cutoff or generic-coefficient
assumption in (21).

The identity-filter implementation also evaluates the fixed polynomial
minors independently. Its direct polynomial cut controls check both
the reduction and these closed formulas in small cases. The all-size
conclusion follows from the regular transformations and exact matrix
factorizations above.
