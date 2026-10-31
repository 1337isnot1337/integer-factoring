# A bounded zero core for coupled identity-mode filters

## Result and scope

In every prime characteristic, an arbitrary joint occupation filter on
the first \(r\) identity-block directions of a Pascal interval has
an exact zero reduction to \(2r\) active modes. The interval can
have any length \(n\ge r\). All other zero modes are explicit
Gaussian spectators. At infinity only the selected \(r\) modes
are active, and their joint digit weights appear directly.

For two qutrit directions, the zero matrix therefore has order
\(3^4=81\), independently of \(n\); the infinity active matrix
has order nine. The filter can correlate the two occupations, so it
need not factor into single-direction filters. This theorem supplies
a finite constructive search family with both cuts and the endpoint
specified. It does not assert a pressure improvement.

## 1. Polynomial identity summands and joint filters

Use the initial Pascal interval \(0,\ldots,n-1\), ambient
parameter \(M=2n\), and characteristic \(p\). Write the
one-particle kernel in physical boundary order \((V,Y)\mid(U,X)\)
as

~~~math
\mathcal D(w)=\begin{pmatrix}A&B\\I&C\end{pmatrix}.
\tag{1}
~~~

For \(0\le j<r\), its polynomial identity-block summand is

~~~math
Q_j(w)=
\left(Y_j+\sum_i A_{ij}(w)V_i\right)
\left(U_j+\sum_i C_{ji}(w)X_i\right).
\tag{2}
~~~

All physical variables have \(p\)-th power zero, so \(Q_j^p=0\).
Let \(\alpha_{\mathbf d}(v)\), indexed by
\(\mathbf d\in\{0,\ldots,p-1\}^r\), be joint digit weights.
For positive odd \(P\), multiply the ramified Gaussian by

~~~math
\mathcal F_\alpha
=\left(\prod_{j=0}^{r-1}E(-Q_j(v^P))\right)
 \sum_{\mathbf d}\alpha_{\mathbf d}(v)
 \prod_{j=0}^{r-1}\frac{Q_j(v^P)^{d_j}}{d_j!},
\qquad E(H)=\sum_{d=0}^{p-1}\frac{H^d}{d!}.
\tag{3}
~~~

Every denominator is a field unit. If every digit weight is regular
at \(v=-1\) and equals one there, (3) equals one at the endpoint,
coefficient by coefficient. Thus the original typed crossing endpoint
is preserved. Polynomial weights give a polynomial tensor. Laurent
polynomial weights can be cleared by one common monomial with endpoint
value one; the two profile shifts cancel in their partition product.
Other rational denominators require their full clearing cost and a
separate endpoint check when they vanish at the endpoint.

For polynomial weights, a direct degree bound is
\(P\deg G+2nPr(p-1)+\max_{\mathbf d}\deg\alpha_{\mathbf d}\).
Indeed, each \(Q_j\) has coefficient degree at most \(2n\),
and \(Q_j^p=0\) truncates every individual exponent. The finite
alphabet and full polynomial degree remain part of the construction
cost; the small profile matrices do not make these parameters free.

## 2. The zero forms have a common Pascal basis

Set \(v=-z^2\), \(t=z^{2P}\), \(c=1-t\), and use the
aggregate zero gauge \(J\) from the
[single-direction reduction](qutrit-identity-zero-reduction.md).
Put

~~~math
a_j=P(2j+1-J),\qquad b_i=P(2i+1+J).
\tag{4}
~~~

Rational gauges are interpreted in a fixed Puiseux extension with
\(\operatorname{ord}(z)=1\), or after an equivalent fixed
ramification. These choices preserve all valuation comparisons.

Removing the regular internal Gaussian factors and applying the
constant Pascal change gives diagonal zero modes \(Y_iU_i\) of
orders \(a_i\), and \(X_iV_i\) of orders \(b_i\). Define

~~~math
f_j(i;c)=\sum_{h=0}^{j}\binom jh\binom ih c^h.
\tag{5}
~~~

Then (2) becomes exactly

~~~math
\boxed{\quad
Q_j=\left(z^{a_j}Y_j+c\sum_{i=0}^{n-1}f_j(i;c)V_i\right)
    \left(U_j-z^{-a_j}\sum_{i=0}^{n-1}f_j(i;c)X_i\right).
\quad}
\tag{6}
~~~

For the first factor, the original column of \(A\) has coefficients
\(\binom ja(1+w)^{a+1}\). The inverse signed-Pascal change on
\(V\) therefore gives \(c f_j(i;c)\).

For the second factor, the signed Delannoy polynomial satisfies

~~~math
t^{i+j}O_{ji}(-t^{-1})=f_j(i;1-t).
\tag{7}
~~~

This is an integer polynomial identity before reduction modulo \(p\).
The coefficient of \(t^\ell\) on its left is
\((-1)^\ell\binom j\ell\binom{i+j-\ell}{j}\).
Expanding the right and applying Vandermonde gives the same coefficient:

~~~math
\sum_{h\ge\ell}\binom jh\binom ih\binom h\ell
=\binom j\ell\binom{i+j-\ell}{j}.
\tag{8}
~~~

The physical monomial gauges then supply the factor \(-z^{-a_j}\)
in (6). For example, \(f_0(i)=1\) and \(f_1(i)=1+ci\).

## 3. The top consecutive minor is a unit in every characteristic

Let \(E\) be the \(n\)-by-\(r\) matrix with entries
\(E_{ij}=f_j(i;c)\). Split its rows into the first \(n-r\)
and the last \(r\) rows:

~~~math
E=\begin{pmatrix}E_L\\F\end{pmatrix}.
\tag{9}
~~~

Write \(V_{ih}=\binom ih\) for \(0\le h<r\), and let
\(P_r=(\binom jh)_{0\le j,h<r}\). Equation (5) gives

~~~math
E=V\operatorname{diag}(1,c,\ldots,c^{r-1})P_r^{\mathsf T}.
\tag{10}
~~~

The bottom consecutive Pascal matrix \(V_T\) has determinant one
over the integers. One direct proof uses
\(\binom{n-r+i}{h}=\sum_\ell\binom i\ell\binom{n-r}{h-\ell}\):
it factors \(V_T\) into two triangular integer matrices with
diagonal entries one. This proof applies even when \(r\ge p\);
it does not divide by any factorial.

Consequently

~~~math
\det F=c^{r(r-1)/2},\qquad
B:=E_LF^{-1}=V_LV_T^{-1}.
\tag{11}
~~~

The first quantity is a unit at zero because \(c(0)=1\).
The matrix \(B\) is constant over the prime field: all powers of
\(c\) and the common Pascal factor cancel. For \(r=2\), its
row indexed by \(i<n-2\) is
\((n-1-i,\ i-n+2)\), reduced in the field.

## 4. Two regular changes concentrate every form simultaneously

Let \(D=\operatorname{diag}(z^{b_i})\), split as
\(\operatorname{diag}(D_L,D_T)\) at the same indices as (9).
The \(X,V\) Gaussian is \(X^{\mathsf T}DV\). A substitution

~~~math
X_{\rm old}=L^{\mathsf T}X_{\rm new},\qquad
V_{\rm old}=RV_{\rm new},\qquad
R=D^{-1}L^{-1}D
\tag{12}
~~~

preserves that diagonal Gaussian. Its left and right matrices must
both be regular; arbitrary rational changes would not preserve Smith
orders.

First take

~~~math
L_1=\begin{pmatrix}I&-B\\0&I\end{pmatrix},\qquad
R_1=\begin{pmatrix}I&D_L^{-1}BD_T\\0&I\end{pmatrix}.
\tag{13}
~~~

These matrices and their inverses are regular. Every entry of the
off-diagonal block of \(R_1\) has a factor
\(z^{b_t-b_i}\), where \(t\) is a retained top index and
\(i\) is a discarded low index. This difference is positive.
Also \(L_1E=(0,F)^{\mathsf T}\).

Define

~~~math
W=I+D_TB^{\mathsf T}D_L^{-1}B.
\tag{14}
~~~

Every entry of the added matrix has positive valuation, so
\(W\equiv I\pmod z\). In particular, \(W\) and its inverse
are regular. Under the inverse transpose in (13), the derivative
form transforms as

~~~math
L_1^{-\mathsf T}D^{-1}E
=\begin{pmatrix}D_L^{-1}BF\\D_T^{-1}WF\end{pmatrix}.
\tag{15}
~~~

Now put

~~~math
C_*=D_TW^{-\mathsf T}B^{\mathsf T}D_L^{-1},\qquad
L_2=\begin{pmatrix}I&0\\C_*&I\end{pmatrix},\qquad
R_2=\begin{pmatrix}I&0\\-W^{-\mathsf T}B^{\mathsf T}&I\end{pmatrix}.
\tag{16}
~~~

The matrices \(L_2,R_2\) are companions as in (12). The matrix
\(C_*\) is regular because \(W^{-\mathsf T}\) is regular and
every retained diagonal order exceeds every discarded one. The
companion \(R_2\) is visibly regular as well. The inverses only
change the signs of their off-diagonal blocks, so they are regular.

The matrix \(L_2\) fixes \((0,F)^{\mathsf T}\). Its inverse
transpose kills the upper part of (15), because

~~~math
C_*^{\mathsf T}D_T^{-1}WF=D_L^{-1}BF.
\tag{17}
~~~

With \(L=L_2L_1\), \(R=R_1R_2\), the two complete form
matrices therefore become

~~~math
LE=\begin{pmatrix}0\\F\end{pmatrix},\qquad
R^{\mathsf T}E=\begin{pmatrix}0\\WF\end{pmatrix}.
\tag{18}
~~~

This proves simultaneous concentration of all \(r\) forms on
each side into the last \(r\) diagonal modes. If \(r=n\), the
discarded block is empty, \(W=I\), and the statement is immediate.
No sign or characteristic restriction was used in this argument.

## 5. The exact active tensor and spectators

Write \(t_\ell=n-r+\ell\) for \(0\le\ell<r\). The active
zero Gaussian and identity summands are exactly

~~~math
G_{0,\rm act}
=\prod_{j=0}^{r-1}E(z^{a_j}Y_jU_j)
 \prod_{\ell=0}^{r-1}E(z^{b_{t_\ell}}X_\ell V_\ell),
\tag{19}
~~~

~~~math
Q_{j,\rm act}
=\left(z^{a_j}Y_j+c\sum_{\ell=0}^{r-1}(WF)_{\ell j}V_\ell\right)
 \left(U_j-z^{-a_j}\sum_{\ell=0}^{r-1}F_{\ell j}X_\ell\right).
\tag{20}
~~~

Apply the same joint filter (3) to (19), using (20). Its zero
coefficient matrix has order \(p^{2r}\). The discarded zero
spectators have orders

~~~math
a_i\quad(r\le i<n),\qquad
b_i\quad(0\le i<n-r).
\tag{21}
~~~

Notice that the final coefficients \(F\) and \(WF\) in (20)
are polynomial in the local parameter after clearing rational gauge
exponents. The inverse of \(W\) is needed only to prove the regular
basis change; it introduces no denominator into the active tensor.

At infinity, the usual regular identity-block Schur reduction isolates
the selected summands as \(r\) independent Gaussian modes of orders

~~~math
\kappa_j=P(2j+1-2n),\qquad 0\le j<r.
\tag{22}
~~~

The joint filter changes the corresponding active digit entries to
\(\alpha_{\mathbf d}\) times their diagonal Gaussian entries.
Use the same ramified infinity coordinate as the Gaussian profiles,
and write
\(\nu_\infty(\alpha)=\operatorname{ord}_{z=0}\alpha(-z^{-2})\).
For a polynomial of degree \(d\), this is \(-2d\).
The complete infinity active orders are therefore

~~~math
\nu_\infty(\alpha_{\mathbf d})
 +\sum_{j=0}^{r-1}\kappa_jd_j,
\qquad \mathbf d\in\{0,\ldots,p-1\}^r.
\tag{23}
~~~

The infinity spectators are the remaining interval orders
\(P(2i+1-2n)\), \(r\le i<2n\). These statements concern
actual tensor factors; no spectators or physical letters are omitted
from the pressure calculation.

For \(r=1\), (14) is the geometric sum from the single-direction
proof, and (20) reduces to
\((z^aY+(1-z^{2nP})V)(U-z^{-a}X)\). Thus the construction
recovers the earlier universal active tensor exactly. For \(r\ge2\),
joint weights can introduce correlations that the single-direction
finite-difference theorem does not decide. Their finite zero matrices
must be analyzed directly or by an additional proved reduction.

The [coupled-filter implementation](../tools/coupled_identity_mode_filters.py)
constructs (19)--(20), applies arbitrary joint digit weights, and
computes both active profiles with their complete Gaussian spectators.
Its [independent checker](../tools/check_coupled_identity_mode_filters.py)
reconstructs the original physical tensors and compares both full cuts.
The [recorded evidence](../evidence/coupled-identity-mode-filters.json)
includes 72 exact Pascal/interpolation/Gram controls, 32 regressions
against the proved single-direction spectrum, and 26 original
full-tensor cut comparisons with correlated binary and qutrit weights.

## 6. A weighted criterion for terminal and mixed directions

The two-shear argument also permits retained indices in other parts of
the Gaussian interval. Let \(T\) be \(r\) retained indices and
\(I\) their complement, each listed increasingly. Let \(E\) be
the form matrix of any \(r\) selected polynomial identity summands.
Suppose \(F=E_T\) is invertible and put

~~~math
B=E_I F^{-1},\qquad K=D_I^{-1}BD_T.
\tag{24}
~~~

**Weighted concentration criterion.** If both \(B\) and \(K\)
are regular, and all the diagonal Gaussian orders are distinct, the
same two shears concentrate both complete form matrices to \(T\)
by regular changes and regular inverse companions.

To prove this, arrange discarded coordinates before retained ones,
and use (13) with \(R_1\)'s off-diagonal block equal to \(K\).
Set \(W=I+K^{\mathsf T}B\). If \(k_i\) denotes a Gaussian
order, regularity in (24) gives
\(\operatorname{ord} B_{ia}\ge\max(0,k_i-k_a)\).
A term of \((W-I)_{ab}\) therefore has valuation at least

~~~math
\max(0,k_a-k_i)+\max(0,k_i-k_b).
\tag{25}
~~~

For \(a=b\), this is \(|k_a-k_i|>0\). For \(k_a>k_b\),
it is at least \(k_a-k_b>0\). Thus \(W\) reduces to an upper
triangular matrix with diagonal entries one, and is a regular unit.
Moreover,

~~~math
W^{\mathsf T}=D_T^{-1}WD_T,\qquad
D_TW^{-\mathsf T}D_T^{-1}=W^{-1}.
\tag{26}
~~~

The second shear coefficient from (16) can consequently be written
\(C_*=W^{-1}K^{\mathsf T}\), which is regular. Its companion
has coefficient \(-W^{-\mathsf T}B^{\mathsf T}\), also regular.
Equations (17)--(18) apply unchanged. No positivity assertion about
an arbitrary quadratic form is used; the ordered valuation bound
(25) is what proves the unit property.

### Terminal directions when the interval length is a prime power

Suppose \(n\) is a power of \(p\), and write \(\tau=1-c=z^{2P}\).
For \(0\le i,j<n\), the signed reversal identity is

~~~math
f_{n-1-j}(i;c)=\tau^i f_j(i;-c/\tau).
\tag{27}
~~~

Indeed, \(\binom{n-1-j}{h}=(-1)^h\binom{j+h}{h}\) in the
field for \(h<n\). Expanding the right side of (27) and applying
Vandermonde gives the same coefficients.

Select the last \(r\) identity summands, and retain the first
\(r\) \(X,V\) indices. The interpolation matrix has the form

~~~math
B_{ih}=\tau^{i-h}(B_0)_{ih},\qquad i\ge r,\quad h<r,
\tag{28}
~~~

where \(B_0\) is a constant consecutive-Pascal interpolation
matrix. Its retained minor has determinant a constant sign times
\(c^{r(r-1)/2}\), hence is a unit. Since
\(D_i=z^{P(1+J)}\tau^i\), the matrix \(K\) in (24) is
exactly \(B_0\). Both matrices are regular, so the criterion gives
\(2r\) zero active modes. This geometry retains the infinity modes
of orders \(-P,-3P,\ldots,-(2r-1)P\).

### The two extreme identity directions

Still assume \(n\ge2\) is a power of \(p\). Select \(Q_0,Q_{n-1}\)
and retain \(X,V\) indices \(0,t\), where \(t=n-1\). The form
row is \(E_i=(1,\tau^i)\), and \(g=1-\tau^t\) is a unit.
For \(0<i<t\), interpolation gives

~~~math
B_{i0}=\frac{\tau^i-\tau^t}{g},\qquad
B_{it}=\frac{1-\tau^i}{g}.
\tag{29}
~~~

Their orders in \(\tau\) are \(i\) and zero. The corresponding
entries of \(K\) have orders zero and \(t-i\), respectively.
Thus (24) applies and again leaves only four zero active modes.

Although \(B\) has denominator \(g\), the final form coefficient
\(WF=F+D_TB^{\mathsf T}D_I^{-1}E_I\) has at most that one
unit denominator. Scaling each active \(V\) variable by \(g\)
clears it by a regular diagonal unit. The auxiliary Gaussian entries
then become \(g z^{b_i}\), for \(i\in\{0,t\}\), and must be retained in the active
matrix; this operation does not assert an unchanged monomial Gaussian
formula. The original global tensor still uses the polynomial
summands (2), so the local auxiliary denominator does not introduce
an endpoint pole into the construction.
