# Exact Gaussian complements at the five-symbol reference

## Result and scope

The new complementary directions are built from the same strict Gaussian
whose five-symbol restriction has numerical ratio approximately
\(0.7263479477\). Their alignment is checked on the complete polynomial
coefficients, not inferred from matching dimensions or physical weights.

There are 65 complementary four-particle states. Seventeen give direct
degree-at-most-18 perturbation slabs after exactly the same monomial
normalization as the reference. Individual common-factor normalization gives
thirty admissible slabs, whose 505 logical placements and monomial shifts
have exact coefficient rank 295 over both \(\mathbf F_2\) and \(\mathbf F_4\).

These are verified construction spaces around the existing reference.
They do not establish a factoring exponent improvement. The larger
[full complement-module saturation](aligned-gaussian-complement-module.md)
uses polynomial combinations between these source rows; its scope is stated
separately.

The constructor is
[`aligned_gaussian_complements.py`](../tools/aligned_gaussian_complements.py),
with the complete
[source record](../evidence/aligned-gaussian-complement-directions-q5.json).
It uses the existing exact selected-minor engine and retains the attribution
of the [reflection compression](reflection-dual-rail-compression.md).

## 1. The literal strict Gaussian

Put \(u=1+w\) and let \(0\le i,j<8\). Define

\[
G_{ij}=w^{2(i\mathbin{\&}j)},\qquad
H_{ij}=w^{2(7-(i\mathbin{|}j))}.
\]

The strict one-particle matrix is

\[
D(w)=\begin{pmatrix}uG&w^{16}I\\I&uH\end{pmatrix}.
\tag{1}
\]

For physical occupation masks \(x,y,u_0,v\), each with four occupied
modes, take the minor with row set \(v\cup(y+8)\) and column set
\(u_0\cup(x+8)\). Here the symbol \(u_0\) denotes a physical mask,
whereas \(u=1+w\) in (1) is a polynomial.

At \(w=1\), matrix (1) is the block swap. Its selected minor is therefore

\[
[x=v][y=u_0].
\tag{2}
\]

The retained physical code, in the reference's label order, is

\[
C=(15,77,43,23,51),
\]

whose occupied mode sets are

\[
\{0,1,2,3\},\quad \{0,2,3,6\},\quad
\{0,1,3,5\},\quad\{0,1,2,4\},\quad\{0,1,4,5\}.
\]

These are the physical states corresponding to reflection labels
\(\{0,2,4,8,12\}\). Let \(G_C(w)\) be their restriction. For the
degree-18 reference \(B\), the exact identity is

\[
G_C(w)=w^{16}B(w^2).
\tag{3}
\]

The constructor checks (3) at every one of its repeated baseline entries.
The independent checker also recomputes all 625 baseline minors directly
from the closed matrix (1), using complete Python polynomial arithmetic.

The occupation grades are \((16,26,22,18,24)\), giving the natural paired
energies \((0,10,6,2,8)\) after subtracting the first grade. Matching these
grades alone would not identify an odd-variable monomial basis with the
Schur fermion basis. Equation (3) is the coefficient-level alignment gate.

## 2. Complementary columns and the common variable

For a four-element occupation set \(a\notin C\), define the slab

\[
K_a(w;y,u_0,v)
=\det D(w)_{\,v\cup(y+8),\,u_0\cup(a+8)},
\qquad y,u_0,v\in C.
\tag{4}
\]

There are \(\binom84-5=65\) choices of \(a\). Equation (2) gives
\(K_a(1)=0\) entry by entry.

Every entry of (4) has only even powers of \(w\). To see this, consider
one determinant term. If it uses \(k\) entries from the upper-left block
of (1), the block row and column counts force it to use \(k\) entries
from the lower-right block. Its factor from these blocks is

\[
(1+w)^{2k}=(1+w^2)^k.
\]

All remaining monomials in (1) already have even exponents. Therefore

\[
K_a(w)=\widetilde K_a(z),\qquad z=w^2,
\tag{5}
\]

with a genuine polynomial slab \(\widetilde K_a\). The implementation
checks the exponent condition in every literal entry before making this
substitution.

For seventeen complementary states, every entry is divisible by \(z^8\)
and the largest entry degree is at most 26. Hence

\[
\Delta_a(z)=z^{-8}\widetilde K_a(z)
\tag{6}
\]

is a polynomial endpoint-zero slab of degree at most 18. This uses exactly
the reference normalization in (3). Replacing a selected \(x\)-boundary
state by itself plus a constant multiple of \(a\), while retaining the
other boundaries, gives the tensor direction (6). The direct family thus
retains the original complementary Gaussian column, including all its
polynomial factors.

The [raw-column export](../evidence/aligned-gaussian-raw-complement-directions-q5.json)
records these seventeen slabs and all five logical \(x\)-placements.
It is reproduced by
[`aligned_gaussian_subfamilies.py`](../tools/aligned_gaussian_subfamilies.py).

## 3. Further common-factor normalization

For each nonzero \(\widetilde K_a\), let \(g_a\) be its common entry
gcd. As in the
[individual kernel-normalization construction](saturated-dual-kernel-perturbations.md),

\[
H_a(z)=(1+z)\frac{\widetilde K_a(z)}{g_a(z)}
\tag{7}
\]

is a polynomial slab with zero endpoint. Its validity is checked directly
as a polynomial tensor. It does not require the corresponding rational
physical split to be defined at the endpoint.

Thirty normalized slabs fit degree 18. For each, every logical placement
\([x=b]z^sH_a\) with nonnegative \(s\) and degree at most 18 is retained.
This gives 505 literal generator directions. Exact complete-coefficient
elimination gives rank 295; scalar extension to \(\mathbf F_4\) preserves
that rank.

A single-coordinate scan in an independent basis does not exhaust the
original 505 generator shapes. To keep that distinction explicit, the
exporter can partition the original nonzero generators into three
independent layers of sizes 295, 200, and 10. Together the layers cover
every original generator exactly once. This supports complete one-generator
tests without discarding a useful correlated combination merely because
it is linearly dependent on earlier generators.

## 4. Independent finite verification

The [independent record](../evidence/aligned-gaussian-complement-controls-q5.json)
checks:

- All 625 baseline coefficients as complete Python determinants of (1).
- Two entries from each complementary slab, for 130 independent minors.
- All 8,125 complementary endpoint-zero entries and all 8,125 exact
  normalization identities \(g_aH_a=(1+z)\widetilde K_a\).
- All 2,125 entries in the seventeen direct degree-18 slabs.
- The exact rank of all 505 candidate directions.

The native source constructor itself checks 40,625 baseline coefficient
equalities while generating the 65 rectangular tensors. Its determinant
engine has an explicit complete-degree bound before using formal series;
the independent Python controls do not call that native determinant routine.

```sh
python3 tools/aligned_gaussian_complements.py \
  --output /tmp/aligned-gaussian-complement-directions-q5.json
python3 tools/check_aligned_gaussian_complements.py \
  evidence/aligned-gaussian-complement-directions-q5.json \
  --output /tmp/aligned-gaussian-complement-controls-q5.json
python3 tools/aligned_gaussian_subfamilies.py \
  evidence/aligned-gaussian-complement-directions-q5.json \
  --raw-output /tmp/aligned-gaussian-raw-complement-directions-q5.json
```

Numerical coefficient scans and row-target searches use these literal
directions as input. Any reported favorable pressure still requires full
profile replay and a rigorous numerical inequality before it can change
the factoring bound.
