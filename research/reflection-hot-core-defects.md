# Three active modes per hot reflection-logical variable

## Result

Let \(M\ge4\) be a power of two and \(m=M/2\). Use the
[reflection-compressed Gaussian](reflection-dual-rail-compression.md),
with \(m\) logical modes on each physical leg. For a power of two
\(k\le m\), select the hottest logical indices
\(m-1,m-2,\ldots,m-k\). Any even exterior polynomial defect
on these selected variables has a regular reduction at both natural
cuts to

\[
k+\min(m,2k)
\tag{1}
\]

active Gaussian modes. Every other mode remains an explicit Gaussian
spectator. Thus one selected mode needs three active modes and an
\(8\)-by-\(8\) matrix; two selected modes need six active modes
and a \(64\)-by-\(64\) matrix, independently of the interval length.

The defect coefficients can be arbitrary Laurent polynomials over
\(\mathbb F_4\), subject to the endpoint condition below. The
reduction does not assume that the defect is Gaussian or that its
monomials preserve occupation. It initially uses the natural gauges;
arbitrary physical state gauges require the corresponding annihilation
and creation support as well. No improved factoring exponent is claimed.

## 1. Polynomial family and endpoint

Write \(G_M(w)\) for the complete polynomial reflection tensor,
including its common scalar \(w^{M^2/2}\). Let \(P\) be a
positive integer. For an even subset \(H\) of the selected physical
variables, write \(X_H\) for its exterior monomial and choose
Laurent polynomials \(f_H(w)\). Require

\[
f_\varnothing(1)=1,\qquad f_H(1)=0\quad(H\ne\varnothing).
\tag{2}
\]

Then

\[
T(w)=G_M(w^P)\sum_H f_H(w)X_H
\tag{3}
\]

has the same literal crossing endpoint. If \(w^h\) clears all
negative powers of the coefficients, \(w^hT(w)\) is a polynomial
tensor with the same endpoint and conservative degree bound

\[
\deg(w^hT)\le 2PM^2+h+\max_H\deg f_H.
\tag{4}
\]

Multiplication by an exterior variable does not raise the degree in
\(w\). The monomial clearing shifts the complete zero profile by
\(+h\) and the infinity profile by \(-h\). These shifts cancel
in the profile product. The Gaussian scalar contributes
\(\pm PM^2/2\) in the same way. No global endpoint-vanishing
denominator is divided out.

## 2. Both cuts have the same regular cross reduction

At either place use its local parameter \(z\), set

\[
s=z^{4P},\quad c=1+z^{2P},\quad F=z^{2MP},
\quad D_i=z^{P(2M-2-4i)}\quad(0\le i<m).
\tag{5}
\]

With the natural logical heights from the compression theorem, order
the zero cut as \((X,Y)\mid(U,V)\) and the infinity cut as
\((U,X)\mid(V,Y)\). The cross matrix is the same at both places:

\[
K=\begin{pmatrix}A&D\\D&A\end{pmatrix},\qquad
A_{ij}=c z^{P(2M-2-4(i\mathbin{\&}j))}.
\tag{6}
\]

All internal Gaussian factors are regular and have regular inverses.
Put \(R_{ij}=s^{i\mathbin{\&}\neg j}\) and \(L=R^{\mathsf T}\).
Then \(D^{-1}A=cR\), \(AD^{-1}=cL\), and the regular matrices

\[
U=\begin{pmatrix}0&I\\I&cL\end{pmatrix},\qquad
V=\begin{pmatrix}I&cR\\0&I\end{pmatrix}
\tag{7}
\]

satisfy \(UKV=\operatorname{diag}(D,FD)\). Their inverses are
regular as well. Reverse the coordinate indices in each diagonal
block, writing \(r=m-1-i\). The complete one-particle orders become

\[
\kappa_{A,r}=P(2+4r),\qquad
\kappa_{B,r}=P(2M+2+4r).
\tag{8}
\]

Their union is exactly \(P(2,6,\ldots,4M-2)\).

## 3. The transported physical multiplication operators

Let \(C_{A,r},C_{B,r}\) denote creation operators in the diagonal
Fock basis, and \(A_{A,r},A_{B,r}\) their annihilators. Put
\(u=z^{2P}\). The selected physical logical index \(m-1-j\),
where \(0\le j<k\), has four multiplication operators:

\[
\begin{aligned}
\ell_{1,j}&=C_{B,j},\\
\ell_{2,j}&=C_{A,j}
 +c\sum_{r=0}^{m-1}s^{r-(j\mathbin{\&}r)}C_{B,r},\\
r_{1,j}&=(us^j)^{-1}A_{A,j}
 +\frac{c}{Fu}\sum_{r=0}^{m-1}s^{-(j\mathbin{\&}r)}A_{B,r},\\
r_{2,j}&=(Fus^j)^{-1}A_{B,j}.
\end{aligned}
\tag{9}
\]

At zero, \((\ell_1,\ell_2,r_1,r_2)=(X,Y,U,V)\); at
infinity they are \((U,X,V,Y)\). A physical input variable acts
by right multiplication by a Fock annihilator. Transporting it
through the diagonal Gaussian gives the two \(r\)-operators in
(9). Confusing this operation with an input Fock creator gives a
larger and incorrect support requirement for polynomial defects.

Formula (9) follows directly from the columns of \(U\) for
physical row multiplication and the rows of \(V\), divided by
the diagonal Gaussian coefficients, for physical input multiplication.
Internal creation factors commute with row multiplication; internal
annihilation factors commute with input multiplication. Thus they
do not add terms to (9).

## 4. Regular concentration in each residue class

Write \(r=t+k\ell\), with \(0\le t<k\), and
\(N=m/k\). Since \(j<k\) and \(k\) is a power of two,
\(j\mathbin{\&}r=j\mathbin{\&}t\). In every residue class,
the spreading creation sequence is proportional to
\((s^k)^\ell\), while the spreading annihilation sequence is
constant. Exceptional coordinate vectors occur only at \(\ell=0\).

If \(N\ge2\), retain \(\ell=0,1\). For every \(\ell\ge2\),
add \((s^k)^{\ell-1}\) times creation row 1 to row \(\ell\),
and apply the inverse transpose to annihilators. This kills the
geometric creation sequence outside the retained rows. It changes
the constant annihilation coefficient at 1 to

\[
g_N=\sum_{h=0}^{N-2}(s^k)^h,
\tag{10}
\]

a regular unit. Then add \(g_N^{-1}\) times row \(\ell\)
to row 1 for every \(\ell\ge2\). This kills the remaining
outside annihilation coordinates and does not reopen a creation
coordinate. Index 0 is untouched.

For a shear with coefficient \(f\), row \(a\), and column
\(b\), its Gaussian inverse companion has coefficient order
\(\operatorname{ord}(f)+\kappa_b-\kappa_a\). In the first
sweep this order is zero. In the second it is nonnegative, while
the forward coefficient has order zero. Thus both changes and all
Gaussian companions are regular. The first block in (8) already
uses only indices \(0,\ldots,k-1\). This proves (1).

The retained indices in the complete increasing list of orders are

\[
\{0,\ldots,k-1\}\ \cup\
\{m,\ldots,m+\min(m,2k)-1\}.
\tag{11}
\]

Every other order in \(P(2,6,\ldots,4M-2)\) is a Gaussian
spectator at both places. The regular basis changes prove an exact
tensor factorization, not just an estimate for the smallest orders.
The literal endpoint of (3) supplies full generic rank, and therefore
full rank of its active matrix.

## 5. Complete active matrices and precision

For one hot mode, the active orders are
\(2P,(2M+2)P,(2M+6)P\). The final operator coefficients are
Laurent polynomials: the unit inverse in the proof cancels from
the retained images. The same holds for every selected block in
(11). The implementation constructs every active entry exactly
before Smith elimination.

The one-mode interface is
[`reflection_hot_core_profiles.py`](../tools/reflection_hot_core_profiles.py).
Its `compile_interval(M, place, ramification=1)` returns all eight
even monomial component matrices, active and spectator orders,
and the omitted common scalar order. The two-mode extension is
[`reflection_hot_block_profiles.py`](../tools/reflection_hot_block_profiles.py),
with `core=2` and 128 even monomial components.

In either interface, coefficient arrays are complete Laurent
polynomials over \(\mathbb F_4\). The native Smith window grows
when required; a source exceeding its supported exact range is
rejected. The floating pressure helper includes all spectators and
the full physical energies
\(4P(M-1-2i)\), for \(0\le i<m\), on both alphabets.
Its optimized roots are exploratory values, not strict certificates.

The [independent checker](../tools/check_reflection_hot_core_profiles.py)
and [saved evidence](../evidence/reflection-hot-core-profiles.json)
distinguish original polynomial endpoint controls, original complete
cut profiles, and the all-size regular-support argument. The extended
saved run passes 28 regular-support cases, 322 original complete
profiles, and 1,150,976 endpoint coefficient equalities.

```sh
python tools/check_reflection_hot_core_profiles.py --extended --output evidence/reflection-hot-core-profiles.json
```
