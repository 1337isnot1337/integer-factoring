# Arbitrary state gauges on a bounded reflection hot core

## Result and scope

The [hot-block reduction](reflection-hot-core-defects.md) for polynomial
multipliers extends to the complete creation and annihilation algebra
of the selected physical states. If \(k\) is a power of two and
the \(k\) hottest logical modes are selected, both cuts have a
regular reduction with at most

\[
\min(M,11k)
\tag{1}
\]

active Gaussian modes. All remaining modes are explicit spectators.
This permits arbitrary diagonal weights on each physical core state,
including weights that make the original internal Gaussian factors
nonregular. It also permits constant coefficient masks on the full
four-leg core, provided that they preserve the literal endpoint.

For one selected logical mode, (1) gives eleven active modes for
every \(M\ge16\). The coefficient matrix has order 2048, but
conserves Fock parity and splits into two square blocks of order 1024.
The implementation keeps the earlier native interfaces unchanged and
uses a separate parity exporter for this larger active space.

This is an exact constructive family. Neither the support bound nor
the floating gauge optimizer asserts an improved factoring exponent.

## 1. The complete transported operators

Use the notation of the hot-block theorem:
\(m=M/2\), \(s=z^{4P}\), \(c=1+z^{2P}\),
\(F=z^{2MP}\), \(u=z^{2P}\), and \(B=1+F\).
Reverse the diagonal indices, so the orders in the two blocks are
\(P(2+4r)\) and \(P(2M+2+4r)\).

For selected logical index \(m-1-j\), write
\(h_{jr}=j\mathbin{\&}r\). In the diagonal Fock basis, the
complete operators are

\[
\begin{aligned}
\ell^+_{1,j}&=C_{B,j},\\
\ell^-_{1,j}&=c\sum_r s^{h_{jr}}C_{A,r}
 +c\sum_r s^{j-h_{jr}}A_{A,r}+A_{B,j},\\
\ell^+_{2,j}&=C_{A,j}+c\sum_r s^{r-h_{jr}}C_{B,r},\\
\ell^-_{2,j}&=B C_{A,m-1-j}
 +cF\sum_r s^{h_{jr}}C_{B,r}+A_{A,j},\\
r^+_{1,j}&=us^j C_{A,j}
 +\frac{B}{us^{m-1-j}}A_{A,m-1-j}
 +\frac c u\sum_r s^{h_{jr}-r}A_{B,r},\\
r^-_{1,j}&=(us^j)^{-1}A_{A,j}
 +\frac c{Fu}\sum_r s^{-h_{jr}}A_{B,r},\\
r^+_{2,j}&=cu\sum_r s^{r+j-h_{jr}}C_{A,r}
 +Fus^j C_{B,j}+\frac c u\sum_r s^{h_{jr}-r}A_{A,r},\\
r^-_{2,j}&=(Fus^j)^{-1}A_{B,j}.
\end{aligned}
\tag{2}
\]

Here \(C\) and \(A\) denote creation and annihilation operators.
The physical letter assignments are the same as in the hot-block
theorem. The right operators are the input Fock operators transported
through the diagonal Gaussian. Physical input-variable multiplication
uses \(r^-\); its occupation projector uses the full canonical pair.

For an independent derivation, let \(J\) be reflection on the
unreversed indices. The two regular internal Gaussian matrices are

\[
S=\begin{pmatrix}BJ&cH_m(s)\\cH_m(s)&BJ\end{pmatrix}.
\tag{3}
\]

For the regular cross reduction \(U,V\), direct multiplication gives

\[
US=\begin{pmatrix}cH&BJ\\0&cFH\end{pmatrix},\qquad
SV=\begin{pmatrix}BJ&cFH\\cH&0\end{pmatrix}.
\tag{4}
\]

Row creators come from the columns of \(U\); row annihilators
have creation part \(US\) and annihilation part \(U^{-1}\).
Input creators have creation part \(\Sigma V^{-1}\) and
annihilation part \(SV\Sigma^{-1}\); input annihilators have
part \(V\Sigma^{-1}\). Reversing indices gives (2).

The left operators are regular. Pulling each right operator back
through \(\Sigma\), by dividing its creation coefficients and
multiplying its annihilation coefficients by the corresponding
diagonal entry, also gives regular operators. Each side satisfies
the complete canonical anticommutation relations.

## 2. Six and five retained coordinates per residue

Write \(r=t+k\ell\), with \(0\le t<k\) and \(N=m/k\).
Because \(j<k\), the bitwise quantity \(h_{jr}\) depends
only on \(t\). Put \(v=s^k\). In each residue class, the
first diagonal block has creation support in

\[
\{e_0,e_{N-1},1,v^\ell\}
\]

and annihilation support in

\[
\{e_0,e_{N-1},1,v^{-\ell}\}.
\tag{5}
\]

The second block has the same two spreading sequences, but only
the exceptional vector \(e_0\). The reflected exceptional
index \(m-1-j\) in (2) lies in the first block and at the last
position of a residue class.

When \(N\ge8\), use four shears to concentrate (5). First
kill the geometric creation sequence using position 1. Then kill
the remaining constant sequence using a high position \(a\).
These two sweeps have coefficients

\[
v^{\ell-1},\qquad
\frac{1+v^{\ell-1}}{1+v^{a-1}}.
\tag{6}
\]

Next kill the inverse-geometric annihilation sequence using a
second high position \(b\), and then the constant annihilation
sequence using position 2. Their coefficients are

\[
v^{b-\ell},\qquad
\frac{1+v^{b-\ell}}{1+v^{b-2}}.
\tag{7}
\]

Use \((a,b)=(N-2,N-3)\) in the first block and
\((a,b)=(N-1,N-2)\) in the second. The retained positions are

\[
\begin{array}{c|c}
A&\{0,1,2,N-3,N-2,N-1\}\\
B&\{0,1,2,N-2,N-1\}.
\end{array}
\tag{8}
\]

The first two sweeps leave every creation coordinate outside their
pivots zero. The last two therefore do not reopen one. Every
denominator displayed in (6)--(7) has constant coefficient one.
The monomial coefficients pay exactly for their Gaussian inverse
companions. The unit coefficients point from a lower to a higher
Gaussian order in the companion calculation. Thus all four sweeps,
their inverses, and their Gaussian companions are regular.

There are six retained coordinates in the first block and five in
the second for each residue. If \(N<8\), it is a power of two
at most four, and keeping all \(M=2kN\) modes is smaller than
\(11k\). This proves (1). Every original operator has zero
support outside the retained coordinates after the changes.

## 3. State projectors, coefficient masks, and diagonal gauges

The complete physical core algebra is generated by the canonical
pairs in (2). Hence the same regular factorization applies to every
core occupation projector and every matrix unit. The full cut
matrix is the resulting active matrix tensored with all unused
Gaussian modes.

Regular physical core coordinates are obtained from the joint
vacuum projector. Its residue has rank \(2^{a-2k}\), where
\(a\) is the active mode count. Indeed, the \(2k\) canonical
pairs generate the full matrix algebra of order \(2^{2k}\);
its unital residue representation has multiplicity \(2^{a-2k}\).
A unit vacuum minor gives a
regular coordinate basis and regular inverse. Creation operators
extend the vacuum basis to all core states; the canonical relations
show that this complete basis is invertible. The construction is
performed separately on the two sides, pulling the input operators
back through \(\Sigma\) first.

In these coordinates, a physical core state gauge is an ordinary
row or column weight. No internal-regularity cone is imposed on
these added weights. Every cold physical bit remains in the type
partition, and every Gaussian spectator remains in the profile
partition. Rational weights can be interpreted after a fixed
ramification; floating weights are used only for exploration.

A constant coefficient mask on the four physical core labels is
also an exact entry mask in these coordinates. If
\(m(x,y,y,x)=1\), the literal endpoint is preserved. More
generally, polynomial exterior defects can be applied by their
standard core creation and input-annihilation matrix units. The
known monomial gauge factors and any common Laurent clearing are
retained explicitly. This introduces no global rational repair.

## 4. Parity export and complete precision accounting

The transported generators are odd and the vacuum projectors are
even. Their coordinate bases preserve parity, with the parity of a
coordinate given by its vacuum parity plus its core occupation
parity. Every allowed exterior defect is even, and coefficient masks
are diagonal in occupation. Therefore both active matrices split
into two equal square parity blocks. Their complete Smith spectrum
is the union of the block spectra.

The separate native exporter
[`native_fock_f4_parity.cpp`](../tools/native_fock_f4_parity.cpp)
checks the row and column parity maps and rejects a nonzero entry
outside the advertised sectors. It supports eleven active modes
without changing the earlier ten-mode interfaces. The scalar
coefficient arrays retain their individual Laurent shifts and the
common absolute precision bound.

The reported builder degree bound, actual per-entry coefficient width,
retained source window, and Smith elimination precision are separate.
All retained source coefficients are packed before elimination.
If a promised pivot is unresolved, the elimination workspace can widen
without discarding or repacking a source coefficient. An incomplete
source remains subject to the independent absolute threshold below.

The Python constructor reports whether its arrays are complete
polynomials. Otherwise, every weighted profile is accepted only
below the certified retained coefficient threshold, including the
least row and column weights. A larger requested window rebuilds
the constructor; unresolved profiles raise an error rather than
being counted as successful searches.

The implementation files are
[`reflection_hot_state_support.py`](../tools/reflection_hot_state_support.py),
[`reflection_hot_state_coordinates.py`](../tools/reflection_hot_state_coordinates.py),
and [`reflection_hot_state_gauges.py`](../tools/reflection_hot_state_gauges.py).
The coordinate matrices at the two natural cuts are identical after
the theorem's boundary permutations. The implementation checks the
two regular core bases before reusing the same numerical Gaussian
matrix with the appropriate physical labels.

## 5. A smaller core for monotone state heights

For polynomial exterior defects, a useful subfamily has only
\(\min(M,4k)\) active modes. Suppose each added state height
is monotone under changing a selected occupation from zero to one:

\[
h_{\nu,\ell}(a+2^j)\ge h_{\nu,\ell}(a)
\quad\text{when }a_j=0.
\tag{9}
\]

Every original internal Gaussian coefficient is regular in the
natural gauge. Its exterior creation terms only add occupations,
so (9) keeps it and its inverse regular after the added gauge.
Exterior polynomial multiplication commutes with these factors.
They may therefore be removed before computing the weighted profile.
This argument does not permit an arbitrary coefficient mask, which
need not commute with the internal Gaussian.

Removing the internal factors deletes the creation part of each
left annihilator and the annihilation part of each right creator
in (2). In each residue class of each diagonal block, the remaining
creators span \(\{e_0,v^\ell\}\), and the remaining
annihilators span \(\{e_0,1\}\). The two regular sweeps
from the polynomial hot-block theorem concentrate these to positions
\(0,1\). There are two positions in each of the two blocks,
for a total of \(4k\) modes. Their exact orders are

\[
\{P(2+4r):0\le r<\min(m,2k)\}
\;\cup\;
\{P(2M+2+4r):0\le r<\min(m,2k)\}.
\tag{10}
\]

All other orders are unchanged Gaussian spectators.

For one hot logical bit, after subtracting each state's constant
height, (9) is simply the positive orthant in the eight added heights.
The active matrix has order 16 for every interval \(M\), and
its parity blocks have order 8. The constructor option
`internal_regular=True` selects this branch. The evaluator checks
(9) on every call and rejects a coefficient mask in this branch.
The eleven-mode constructor remains available for arbitrary signs.

## 6. Independent controls

The [portable checker](../tools/check_reflection_hot_state_gauges.py)
and its [saved record](../evidence/reflection-hot-state-gauges.json)
reconstruct the transported operators directly from the original
alternating kernels. They also compare the complete original cut
matrices with the reduced matrices under nontrivial dyadic state
heights, including negative heights in the general branch and
nonlinear two-bit heights.

The extended replay checks 2,944 original-kernel operator coefficients,
56 complete weighted profiles, and 328,960 literal endpoint
coefficients. All 27 finite-difference derivative controls agree to
less than \(1.7\times10^{-10}\). A separate high-height control
uses complete sources of at most 14 coefficients, widens the Smith
workspace from 512 to 4096 coefficients, and agrees with two original
matrix profiles computed at precision 8192.

The search driver
[`search_reflection_hot_state_gauges.py`](../tools/search_reflection_hot_state_gauges.py)
supports the fast monotone branch and the general arbitrary-sign
branch. It records every coefficient row, the actual state heights,
all Gaussian spectators, the complete active orders, and the numerical
optimizer's status. These optimization records are finite exploratory
screens and are not proofs of global optimality.
