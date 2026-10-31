# Exact rational module orbits

This note gives an exact alternative to finite-jet module states when the
normalized row-lattice orbit is finite. It retains rank-deficient lattices,
including directions that vanish identically. A completed orbit computes the
grouped-channel word pressure exactly.

The worked example is a two-dimensional companion crossing. Its two cuts have
exact orbits with two and seven states. The resulting fixed-factorization
pressure satisfies $B(s)>10$ for every $s>0$. It therefore cannot improve
the elementary certificate $256/27$. This is a conclusion about this
factorization's word bound, not a lower bound on the actual coefficient ranks
or on every factorization of the crossing. There is no improvement to the
factoring exponent in this note.

The main result is the finite-orbit certificate and its exact checker. It
complements [finite-jet modules](cancellation-module-smith-automata.md) and
[minimal boundary representations](minimal-boundary-module-completeness.md).

## 1. Exact row lattices

Let $\mathbb F$ be a finite field, $\mathcal O=\mathbb F[[\tau]]$, and
$K=\mathbb F(\tau)$. Each channel $i$ has a finite collection of regular
matrices

\[
U_i(a)\in (K\cap\mathcal O)^{D\times D}.
\]

The physical choice $a$ is internal to a channel. A channel word is counted
once; its order is minimized over its physical choices. Starting from the free
row lattice $M_\varnothing=\mathcal O^D$, set

\[
M_{wi}=\sum_a M_wU_i(a),\qquad
v(w)=\operatorname{ord} M_w.
\tag{1}
\]

For a nonzero row lattice, its order is the minimum valuation of all entries
in any generating set. This is also the minimum entry valuation among all
physical matrix products for $w$. The zero lattice has infinite order.

Each $M_w$ has a finite basis with entries in $K\cap\mathcal O$, even though
the coefficients in its $\mathcal O$-span may be arbitrary formal series.
Its rank can be smaller than $D$. We do not insert a tail
$\tau^P\mathcal O^D$ into a missing direction.

## 2. A canonical basis with exact span witnesses

We use row Hermite reduction over the discrete valuation ring $\mathcal O$.
Process columns from left to right. Among the remaining rows with a nonzero
entry in the current column, choose an entry of minimum order $v$. Divide
the pivot row by its unit part, so the pivot becomes $\tau^v$. Every other
remaining entry in that column is a regular multiple of the pivot, so it can
be eliminated by a regular row operation.

For an earlier pivot row, reduce its entry in the current column modulo
$\tau^v$. Subtracting its polynomial remainder leaves a multiple of
$\tau^v$, so this reduction also uses a regular row operation. The resulting
basis has increasing pivot columns, pure-power pivots, and polynomial entries
of degree smaller than the pivot order above each pivot. Missing pivots stay
missing.

These operations preserve the exact row lattice. The checker verifies both
directions independently of any finite-jet comparison:

- Each returned basis row has a recorded regular linear combination of the
  original generators.
- Each original generator reduces to zero against the returned basis using
  regular coefficients.

The rational arithmetic is exact. Taking a remainder modulo $\tau^v$ uses
only a finite jet, but this is an exact remainder operation inside the Hermite
algorithm. It does not approximate the row lattice.

The normal form is canonical. One can prove uniqueness by comparing the first
pivot column of two bases for the same lattice, then passing to the submodule
with that coordinate zero and continuing. The bounds on entries above each
pivot give the same unique polynomial remainders. The finite-orbit theorem
below only needs exact equality of the represented lattices; canonical form
provides a useful way to recognize that equality.

## 3. Finite exact orbit theorem

For a normalized nonzero lattice $L\subseteq\mathcal O^D$ with
$\operatorname{ord}L=0$, define

\[
R_i(L)=\sum_a LU_i(a),\qquad
\gamma(L,i)=\operatorname{ord}R_i(L),\qquad
\Phi_i(L)=\tau^{-\gamma(L,i)}R_i(L).
\tag{2}
\]

A zero image is a dead transition. Otherwise $\Phi_i(L)$ is another
normalized exact lattice, possibly of smaller rank.

**Theorem.** Suppose the normalized lattices reachable from $\mathcal O^D$
under (2) form a finite set. For every channel word, its exact order $v(w)$
is the sum of the transition orders along the corresponding path. A word
meeting a dead transition has $v(w)=\infty$.

**Proof.** If $M_w=\tau^gL$, then

\[
M_{wi}=\tau^gR_i(L)
       =\tau^{g+\gamma(L,i)}\Phi_i(L).
\]

The last lattice has order zero. Induction, starting with $g=0$, proves the
claim. Every operation is an equality of full lattices. $\square$

Let $h_i\in\mathbb Z$ be the channel weights from a local factorization and
fix $s>0$. The exact free-boundary partition function is

\[
Q_t(s)=\sum_{|w|=t}s^{\sum_{j=1}^t h_{w_j}-v(w)},
\tag{3}
\]

where a zero word contributes zero. The finite matrix

\[
T(s)_{L,L'}=
\sum_{\substack{i:\ \Phi_i(L)=L'}}
s^{h_i-\gamma(L,i)}
\tag{4}
\]

satisfies the exact identity

\[
Q_t(s)=e_{\mathcal O^D}^{\mathsf T}T(s)^t\mathbf1.
\tag{5}
\]

All states are reachable, and the final vector is positive. Hence the
exponential limsup rate in (5) is the spectral radius of $T(s)$. This is the
usual finite nonnegative-matrix argument: Jordan normal form gives the upper
rate, and a reachable component of maximal spectral radius gives the lower
rate. If all sufficiently long paths are dead, both rates are zero.

After reachable and observable minimization, the boundary completeness theorem
identifies this free rate with the limsup rate of the exact fixed-boundary
grouped-channel partition function. The method still does not identify
cancellations between different channel words in a coefficient matrix.

There is no claim that every input has a finite orbit. The checker stops
inconclusively when a state or coefficient-degree budget is exceeded. It
never turns a budget failure into an exact certificate.

## 4. A companion crossing with a rank-deficient exact state

Use physical labels $x,y,u,v\in\{0,1\}$, and work over $\mathbb F_2$. Put

\[
C_{xyuv}=[x=u][y=v],\qquad E_{xyuv}=[x=v][y=u],\qquad N=E+C.
\]

Let $e=(1,1)^{\mathsf T}$, $f=(1,0)$, and

\[
K_{xyuv}(z)=
\begin{pmatrix}
zN_{xyuv}&C_{xyuv}\\
zE_{xyuv}&0
\end{pmatrix}.
\tag{6}
\]

For every physical tuple, $K(1)e=Ee$, and $fe=1$. Consequently

\[
fK_{x_1y_1u_1v_1}(1)\cdots K_{x_ty_tu_tv_t}(1)e
=\prod_{j=1}^t E_{x_jy_ju_jv_j}.
\tag{7}
\]

Thus (6) gives a valid ordered crossing at every length. The constant
auxiliary matrices are strictly upper triangular, so their products vanish
after two factors. The highest scalar coefficient is $N^{\otimes t}$,
whose row-cut rank is one. These favorable endpoint structures motivated the
candidate.

For the checker, conjugate by

\[
S=\begin{pmatrix}1&0\\1&1\end{pmatrix},\qquad S^2=I.
\]

This changes the boundary vectors to $e_0^{\mathsf T}$ and $e_0$. It
preserves the scalar crossing at every length.

We take the checker-generated Schur factorization with both auxiliary indices
on side zero. The infinity row cut has weights $(1,1,1,1)$; the zero column
cut has weights $(0,-1,-1,-1)$. Reachable and observable minimization leaves
dimension two on both sides.

The exact infinity orbit is

\[
\mathcal O^2,\qquad \mathcal O(1,0)+\mathcal O(0,\tau).
\]

Its weighted transition matrix is

\[
T_\infty(s)=\begin{pmatrix}0&4s\\3&s\end{pmatrix}.
\tag{8}
\]

The exact zero orbit has seven states. In the checker's order, they are the
row spans of

\[
\begin{aligned}
&\begin{pmatrix}1&0\\0&1\end{pmatrix},\quad
\begin{pmatrix}1&1\\0&\tau\end{pmatrix},\quad
\begin{pmatrix}1&0\end{pmatrix},\quad
\begin{pmatrix}1&0\\0&\tau\end{pmatrix},\\
&\begin{pmatrix}\tau&0\\0&1\end{pmatrix},\quad
\begin{pmatrix}1&1+\tau\\0&\tau^2\end{pmatrix},\quad
\begin{pmatrix}1&\tau\\0&\tau^2\end{pmatrix}.
\end{aligned}
\tag{9}
\]

The third state has rank one. Finite full-lift states always contain
$\tau^P\mathcal O^2$, so they cannot represent this exact state at any
finite precision.

States with indices $0$, $1,5$, $2,3,6$, and $4$ respectively have
identical weighted outgoing sums to those four classes. The quotient matrix is

\[
T_0(s)=
\begin{pmatrix}
s^{-1}&1&2s^{-1}&0\\
s^{-1}&0&3s^{-1}&0\\
0&1&2s^{-1}&s^{-1}\\
0&1+s^{-1}&2s^{-2}&0
\end{pmatrix}.
\tag{10}
\]

The checker verifies this quotient as an identity of Laurent polynomials in
$s$, using the exact transition orders. It is not inferred from numerical
evaluation. The quotient preserves (5), since the final vector is constant
on each class.

## 5. A uniform lower bound for this factorization

Let $h=\rho(T_\infty(s))$. Equation (8) gives

\[
h^2-sh-12s=0,\qquad s=\frac{h^2}{h+12}.
\tag{11}
\]

The middle two states of (10) give the principal submatrix

\[
\begin{pmatrix}0&3s^{-1}\\1&2s^{-1}\end{pmatrix}.
\]

Its positive eigenvalue $\ell$ satisfies

\[
s\ell^2-2\ell-3=0,\qquad
\ell=\frac{1+\sqrt{1+3s}}{s}.
\tag{12}
\]

Nonnegative walk counts give $\rho(T_0(s))\ge\ell$. Evaluate the left side
of (12) at $10/h>0$. Using (11), it becomes

\[
s\left(\frac{10}{h}\right)^2-2\frac{10}{h}-3
=\frac{-3h^2+44h-240}{h(h+12)}<0.
\tag{13}
\]

The numerator is negative for every real $h$, because its leading
coefficient is negative and its discriminant is
$44^2-4\cdot3\cdot240=-944$. Thus $\ell>10/h$, proving

\[
\boxed{\rho(T_\infty(s))\rho(T_0(s))>10>\frac{256}{27}}
\qquad(s>0).
\tag{14}
\]

This rules out an improvement from increasing jet precision for this fixed
factorization. A different channel factorization, a different crossing, or
cancellation between coefficient terms is outside the conclusion.

As an exact upper check at $s=5$, the evidence supplies positive rational
vectors solving

\[
\left(\frac{266}{25}I-T_\infty(5)\right)h_\infty=\mathbf1,
\qquad
\left(\frac{607}{500}I-T_{0,\mathrm{seven}}(5)\right)h_0=\mathbf1.
\]

They certify the strict upper bound

\[
B(5)<\frac{80731}{6250}=12.91696.
\tag{15}
\]

No floating-point eigenvalue calculation is used in the certificate.

## 6. Replay and scope

Run from the repository root:

```sh
python3 tools/check_exact_module_orbits.py --output evidence/exact-rational-module-orbits.json
python3 -m unittest discover -s tests -p test_exact_module_orbits.py
```

The checker verifies exact two-sided lattice span identities at every
transition, completed orbit closure, local endpoint identities, minimality,
all 170 channel words through length three by independent literal products,
and 300 independent finite-jet span comparisons for random rational lattices.
It also checks the symbolic four-state quotient and both rational positive
potentials.

The [evidence](../evidence/exact-rational-module-orbits.json) records every
state, every transition, all weights, and the rational potential vectors.
The [checker](../tools/check_exact_module_orbits.py) is reusable for other
regular rational matrix families whose exact normalized lattice orbit closes.
