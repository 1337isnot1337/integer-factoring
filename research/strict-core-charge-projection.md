# Fixed-core charge projections on strict Pascal intervals

## Result and scope

A nonlinear charge depending on a fixed physical core can be projected
onto a long strict Gaussian interval without increasing the polynomial
degree or deleting physical states. If the interval length \(M\) and
core length \(k\) are powers of two, with \(k\le M\), all profiles
reduce to at most \(4k\) active fermion modes and explicit Gaussian
spectators. The bound is independent of \(M\).

For two physical core modes, each cut therefore uses a matrix of order
at most \(256\). All \(3^4\) charge functions reduce to fourteen
classes under affine changes of charge. The same theorem places the
four-mode indicator-12 extension in sixteen active modes, but the
current dense evaluator does not construct its \(65{,}536\)-order
matrix. This is an admissible family and an exact profile reduction;
no improved factoring exponent is asserted.

## 1. Strict repair simplifies the interval kernel

Write \(u=1+w\), \(t=w^2\), and index all one-particle blocks by
\(0\le a,b<M\). The strict Pascal repair gives

~~~math
\mathcal D(w)=
\begin{pmatrix}uG&w^{2M}I\\I&uH\end{pmatrix},\qquad
G_{ab}=t^{a\mathbin{\&}b},\quad
H_{ab}=t^{M-1-(a\mathbin{|}b)}.
\tag{1}
~~~

Indeed, for the binary Pascal matrix \(P\),
\(P^2=I\) and
\(P\operatorname{diag}((1+t)^a)P^{\mathsf T}=G\).
Applying \(P\) to the original \(V\) rows therefore changes the
upper-right Pascal block into the identity and gives (1). It is the
same [strict repair](strict-gaussian-relabeling-crossings.md), with
no division and no increase in the maximum entry degree.

The one-bit factors satisfy \(GH=(1+t)^{M-1}I\). Consequently

~~~math
\mathcal D(w)^{-1}
=\begin{pmatrix}uH&w^{2M}I\\I&uG\end{pmatrix},
\qquad \det\mathcal D(w)=1.
\tag{2}
~~~

At \(w=1\), the diagonal blocks of (1) vanish and both other
blocks are identities. Its full exterior coefficient tensor \(C\)
therefore has the strict local endpoint
\(C(1;x,y,u,v)=[x=v][y=u]\).

## 2. The phase projection and its exact endpoint

The core consists of physical mode indices \(0,\ldots,k-1\).
For any function \(f:\{0,1\}^k\to\mathbb Z/3\), define

~~~math
C_f(x,y,u,v)=
\mathbf1_{f(x)+f(u)=f(y)+f(v)}\,C(x,y,u,v),
\tag{3}
~~~

where each argument of \(f\) is restricted to its core coordinates.
The endpoint is unchanged. The coefficients remain binary, every
physical state is retained, and (3) does not increase the maximum
entry degree.

Let \(\omega\in\mathbb F_4\) have order three, and let \(S\)
be the physical diagonal map with entry \(\omega^{f(\mathrm{core})}\).
Equation (3) is the sum of the three phase conjugates with factors
\(S_x^jS_u^jS_y^{-j}S_v^{-j}\), for \(j=0,1,2\).
The sum is one when the charge difference vanishes and zero otherwise.
There is no averaging denominator: three equals one in characteristic
two. This also proves exact descent from \(\mathbb F_4\) to
\(\mathbb F_2\).

Every diagonal function of the core state is a polynomial in the
\(k\) physical number operators \(c_i a_i\). It suffices to
bound the support of their creation and annihilation generators at
both cuts. Bounding only one kind of generator would not suffice.

## 3. An explicit regular infinity diagonalization

Use the physical gauges below, in boundary order \((x,y,u,v)\):

~~~math
h_0(i)=(2i+1-2M,\ 2i+1,\ 0,\ 0),\qquad
h_\infty(i)=(M,\ -M,\ 2i+1-M,\ 2i+1-M).
\tag{4}
~~~

They agree with the earlier strict-crossing convention up to a gauge
that vanishes on every particle-balanced coefficient.

At infinity put \(z=w^{-1}\), \(t=z^2\), and now \(u=1+z\).
The gauged matrix (1) becomes

~~~math
\begin{pmatrix}A&D_-\\D_-&A\end{pmatrix},\qquad
A_{ab}=u z^{1-2M+2(a\mathbin{|}b)},\quad
(D_-)_{aa}=z^{1-2M+2a}.
\tag{5}
~~~

Let \(R=AD_-^{-1}\) and \(L=D_-^{-1}A\). Their entries are
\(u t^{a\mathbin{\&}\neg b}\) and
\(u t^{b\mathbin{\&}\neg a}\), respectively, so both are
regular. The one-bit product identity gives

~~~math
AD_-^{-1}A=(1+z^{2M})D_-.
\tag{6}
~~~

To see this directly, sum over each bit of the middle index. If the
two outer bits differ, the two summands cancel. For equal outer bits,
the product is \(t^a\prod_j(1+t^{2^j})\). Multiplying by
\(u^2\) changes that product into \(t^a(1+t^M)\), proving (6).

It follows that

~~~math
\begin{pmatrix}0&I\\I&R\end{pmatrix}
\begin{pmatrix}A&D_-\\D_-&A\end{pmatrix}
\begin{pmatrix}I&L\\0&I\end{pmatrix}
=\operatorname{diag}(D_-,D_+),\qquad
(D_+)_{aa}=z^{1+2a}.
\tag{7}
~~~

Both outside matrices and their inverses are regular. The complete
one-particle orders are \(1-2M,3-2M,\ldots,2M-1\).

At zero the internal \(XY\) and \(VU\) Gaussian coefficients
are both \(u t^{a\mathbin{\&}b}\), which are regular. Removing
their invertible multiplication maps leaves two copies of the diagonal
\(\sigma_a=z^{2a+1}\). These operations establish the Gaussian
background at both places before the core projection is imposed.

## 4. The conjugate operators have two sequences per residue class

The exact transported generator formulas are implemented in
[`strict_images`](../tools/strict_core_charge_profiles.py). Here is the
support information they imply. It includes both the original left
operators and the right operators transported through the Gaussian
diagonal.

Write \(a=r+kj\), with \(0\le r<k\), \(0\le j<M/k\).
Since \(k\) is a power of two and every core index \(i<k\),

~~~math
a\mathbin{\&}i=r\mathbin{\&}i,\quad
a\mathbin{\&}\neg i=kj+(r\mathbin{\&}\neg i),\quad
a\mathbin{|}i=kj+(r\mathbin{|}i).
\tag{8}
~~~

Put \(s=z^{2k}\). In each residue class of each diagonal block,
the possible creation vectors and annihilation covectors lie in:

| Place | Creation span | Annihilation span |
| --- | --- | --- |
| Infinity | \(e_0,(s^j)_j\) | \(e_0,(1)_j\) |
| Zero | \(e_0,(1)_j\) | \(e_0,(s^{-j})_j\) |

Scalar multiples may have nonzero valuations. This does not affect
the common concentration procedure. The Gaussian order increases by
\(2k\) between consecutive indices in the residue class.

For example, at infinity the transported \(x_i\) creation has
first-block coefficient
\(u z^{1-2M+2(a\mathbin{|}i)}\), while the transported \(u_i\)
annihilation has second-block coefficient
\(u z^{-1}t^{-(a\mathbin{\&}i)}\). They give the geometric
creation sequence and the constant annihilation sequence in the
first row of the table. At zero, the noncoordinate creation
coefficients are \(u t^{a\mathbin{\&}i}\), and the
noncoordinate annihilation coefficients are
\(u t^{a\mathbin{\&}i}/\sigma_a\), giving the second row.

## 5. Regular concentration at both places

For a diagonal Gaussian with orders \(\kappa_i\), a shear
\(I+fE_{ij}\) has regular inverse companion
\(I+f z^{\kappa_j-\kappa_i}E_{ij}\) provided

~~~math
\operatorname{ord}(f)\ge0,\qquad
\operatorname{ord}(f)+\kappa_j-\kappa_i\ge0.
\tag{9}
~~~

The shear adds \(f\) times creation coordinate \(j\) to
coordinate \(i\), and \(f\) times annihilation coordinate
\(i\) to coordinate \(j\). Thus it preserves every canonical
and cross pairing.

Let \(N=M/k\ge2\), and put
\(F(s)=\sum_{h=0}^{N-2}s^h\), a unit at zero.

At infinity, retain indices \(j=0,1\) in every residue class.
For each \(j\ge2\), add \(s^{j-1}\) times row 1 to row
\(j\). This kills the geometric creation sequence outside the
retained set and multiplies the constant annihilation coordinate at
1 by \(F\). Next add \(F^{-1}\) times row \(j\) to row 1
for every \(j\ge2\). This kills the remaining annihilation
coordinates and does not reopen any creation coordinate. The first
sweep pays exactly for its reverse companion; the second has a
regular forward companion.

At zero, retain \(j=0,N-1\). For each \(1\le j<N-1\), add
the last row to row \(j\). This kills the constant creation
sequence outside the retained set and multiplies the last inverse
geometric annihilation coordinate by \(F\). Then add
\(s^{N-1-j}F^{-1}\) times row \(j\) to the last row. This
kills every remaining outside annihilation coordinate. Again the
creation support stays fixed, and (9) holds for every operation.

No operation touches index zero. Each of the \(k\) residue
classes in each of two blocks therefore retains only two modes.
When \(N=1\), keep all \(2M\) modes. This proves the active
bound \(\min(2M,4k)\), with explicit spectators:

~~~math
\begin{array}{c|c}
\text{place}&\text{retained one-particle indices in each block}\\ \hline
\infty&\{0,\ldots,2k-1\}\\
0&\{0,\ldots,k-1\}\cup\{M-k,\ldots,M-1\}.
\end{array}
\tag{10}
~~~

Every matrix and its inverse companion is regular. After these
changes, every projected operator acts on the active factors only.
The full tensor is thus the active matrix tensored with all unused
Gaussian modes. Its complete profile partition is the active
partition times \(\prod_{i\text{ spectator}}(1+s^{-\kappa_i})\).
The strict endpoint supplies full generic cut rank, hence full active
rank \(2^{\min(2M,4k)}\). Rank is not inferred from missing pivots
in a truncated series computation.

## 6. Core gauges and finite searches

The core-support theorem permits arbitrary diagonal weights on the
core state of every physical leg. The regular vacuum-projector
construction from the
[core-coordinate proof](rectangular-defect-factorization.md#7-regular-coordinate-bases-for-the-core-operators)
places these number operators in ordinary coordinate form. The right
operators are first pulled back through the Gaussian diagonal; their
integrality follows from the regular inverse companions proved above.
Unit vacuum minors certify both complete coordinate bases.

Two further real bulk parameters are available. At zero, add
\((J,-J,0,0)\) per occupied mode in boundary order \((x,y,u,v)\).
At infinity, add \((0,K,0,K)\). The corresponding active row
charges are the difference of the two zero block occupations and the
total infinity occupation. Spectators shift by the same signed
charges. The paired physical energies become
\(4i-2M+2+K+J\) and \(4i-2M+2+K-J\), with the additional
core-state weights retained exactly.

The [profile evaluator](../tools/strict_core_charge_profiles.py) uses
all Gaussian spectators and every physical letter. The
[core-gauge evaluator](../tools/strict_core_charge_gauges.py) adds the
regular coordinate bases. For two core modes, subtracting \(f(0)\)
normalizes one charge to zero; multiplication by two pairs the 26
nonconstant normalized functions. Thus fourteen representatives,
including the constant function, cover every charge projection.

The [independent checker](../tools/check_strict_core_charge_profiles.py)
verifies both conjugate supports, all canonical and cross pairings,
regularity before transporting the right generators, and complete
original cut profiles. The
[recorded evidence](../evidence/strict-core-charge-profiles.json)
includes 34 original cut comparisons and two different active-space
reductions with identical complete spectator spectra. These controls
support the implementation; the uniform bound follows from (8)--(10).

The independent
[weighted-core checker](../tools/check_strict_core_charge_gauges.py)
also compares arbitrary quarter-integer physical core weights against
the original complete cut matrices. Its
[evidence](../evidence/strict-core-charge-gauges.json) contains 18
weighted cut comparisons and all 27 independent parameter derivatives.
The coordinate bases intertwine each physical core number operator, so
all charge classes can be obtained by exact entry masks of one cached
Gaussian matrix in those coordinates. The implementation uses this
identity to avoid rebuilding a fermionic program for every charge.
