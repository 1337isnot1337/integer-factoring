# Complete crossing synthesis from nonbinary reference frames

## Result and scope

The complete-entry affine constructor now accepts a reference tensor over
the four-element field, rather than requiring a binary reference. A retained
nonbinary tensor can therefore supply the regular frames for the next exact
coefficient search. Every entry coefficient remains free, subject to the
strict crossing endpoint, the stated polynomial degree, and the selected
local row conditions.

The new implementation is
[`f4_recentered_lattices.py`](../tools/f4_recentered_lattices.py). Its independent
finite controls include binary restriction of scalars, complete Toeplitz
ranks, direct witness substitution, and agreement with the previous binary
constructor. This is a construction tool and an ordinary algebraic argument.
It does not establish a factoring exponent above one half.

The predecessor is the
[complete-entry affine construction](full-alphabet-two-cut-affine-synthesis.md),
with the [finite-precision frame method](finite-precision-smith-frames.md).
The previous coefficient/height searches could reuse selected binary planes
when those planes happened to retain the same profiles. That condition is
not required here. The construction and checks used substantial OpenAI
assistant support and retain the predecessor's source attribution.

## 1. Full coefficient space over the four-element field

Write

\[
F=\mathbf F_2[\omega]/(\omega^2+\omega+1).
\]

Let \(B(w;x,y,u,v)\) be any \(F\)-polynomial tensor of degree at most
\(D\), on an alphabet of size \(q\), with

\[
B(1;x,y,u,v)=[x=v][y=u].
\tag{1}
\]

Every other tensor of degree at most \(D\) with the same endpoint has
the unique form

\[
T=B+(1+w)Q,\qquad \deg Q_{xyuv}<D.
\tag{2}
\]

Indeed, each entry of \(T-B\) vanishes at one. Division by the monic
polynomial \(w+1\) is exact and lowers degree by one. Thus (2) has exactly
\(Dq^4\) independent unknowns over \(F\), irrespective of the
reference tensor's coefficient colors. No occupation, parity, or Gaussian
condition is imposed.

At infinity use the cut \((x,u)\mid(y,v)\); at zero use
\((x,y)\mid(u,v)\). Eight independent integer physical height lists are
retained. If \(h_e\) denotes the sum of the four heights on one entry,
the matrices used for the finite equations are

\[
A_0(T)_{xyuv}=z^{c_0+h_0(x,y,u,v)}T(z;x,y,u,v),
\]

\[
A_\infty(T)_{xyuv}
 =z^{c_\infty+h_\infty(x,y,u,v)}T(z^{-1};x,y,u,v),
\]

where

\[
c_0=-\sum_l\min_a h_{0,l}(a),\qquad
c_\infty=D-\sum_l\min_a h_{\infty,l}(a).
\tag{3}
\]

These are polynomial matrices for every choice of the unknown coefficients.
Both scalar shifts are removed from the reported profiles. The original
degree \(D\) remains explicit throughout the search.

## 2. Exact finite polynomial frames certify complete orders

Let \(A\) be an \(n\)-square polynomial matrix, where \(n=q^2\).
The frame routine computes actual polynomial matrices \(L,R\), with
invertible constant terms, such that

\[
LAR\equiv\mathrm{diag}(d_1,\ldots,d_n)\pmod {z^P},
\qquad \sigma_i=\mathrm{ord}(d_i)<P,
\tag{4}
\]

and the \(\sigma_i\) are sorted.

Every elimination uses a least-order pivot. Writing its value as
\(z^a u\), with \(u(0)\ne0\), an off-diagonal entry is eliminated by
multiplying its row or column by \(u\) and subtracting the pivot row or
column times the corresponding polynomial quotient. These operations are
regular units. Truncating their polynomial representatives modulo \(z^P\)
retains their constant terms and all congruences below \(P\). The routine
checks (4) and the ranks of \(L(0),R(0)\) explicitly.

Equation (4) certifies the complete profile of the untruncated matrix. In
fact, divide row \(i\) of the actual product \(LAR\) by
\(z^{\sigma_i}\). Every error term has order at least
\(P-\sigma_i>0\), and the remaining diagonal entries have nonzero
constant terms. The resulting matrix is a regular unit. Hence

\[
LAR=\mathrm{diag}(z^{\sigma_1},\ldots,z^{\sigma_n})U,
\qquad U(0)\text{ invertible}.
\tag{5}
\]

The \(\sigma_i\) are therefore the full Smith orders. A proposed precision
that misses any of the \(n\) pivots is rejected. The native profile engine
supplies a useful precision hint; the exact polynomial identity and (5)
are the independent complete-profile check.

For a requested row vector \(\tau\), the constructor chooses
\(P\ge\max_i\tau_i\). It then requires

\[
[z^k](L A_e(T)R)_{ij}=0
\qquad(k<\tau_i).
\tag{6}
\]

All equations in (6) are linear over \(F\) in the same coefficient
unknowns from (2). The reversed polynomial and the extra \(1+z\) factor
at infinity are included. The two places are solved simultaneously.
An out-of-window request is rejected rather than interpreting unknown
coefficients as zero.

For a feasible tensor, (6) factors the transformed matrix as
\(\mathrm{diag}(z^{\tau_i})C\), with \(C\) regular. Every
\(k\)-square minor has order at least the sum of the \(k\) smallest
row targets. This supplies the stated Smith-prefix bounds. The implementation
also recomputes the realized complete profiles; target bounds alone are not
reported as a pressure result.

## 3. Affine elimination and repeated reference changes

An equation over \(F\) is stored as two binary coefficient planes and one
field-valued right-hand side. Multiplication by \(\omega\) sends the two
planes \((a,b)\) to \((b,a+b)\). Pivot coefficients are normalized to
one by exact field inversion. The affine system is not separated into an
independent binary affine problem and a binary homogeneous problem: that
separation is generally invalid once the reference frames are nonbinary.

Common row floors are eliminated once. Every later equation is reduced
against those pivots and expressed in the remaining free coordinates.
This is an exact parameterization of the full floor fiber. The checker
compares feasibility and dimension against direct elimination in the
original coordinates.

After a feasible tensor has been reconstructed and checked, it may replace
\(B\). Its endpoint still satisfies (1), so (2) again describes the entire
degree-bounded endpoint fiber. Its new regular frames can define different
row conditions even when its two unordered profiles equal those of the old
reference. Recentered searches remain finite searches of selected frames;
they are not exhaustive searches of all tensors with a target profile.

## 4. A smaller literal source for the five-symbol search

The earlier [physical-code search](physical-code-gauge-searches.md) retained
the interval-eight reflection code with physical letters

\[
\{0,2,4,8,12\}.
\]

Its 625 literal entries all have a common factor \(w^{16}\). After removing
that factor, every exponent is even, and the largest exponent is 36.
Consequently the substitution

\[
T_5(z)=w^{-16}G_{\{0,2,4,8,12\}}(w),\qquad z=w^2,
\tag{7}
\]

gives a degree-18 strict crossing on five symbols. Equation (7) is checked
entry by entry by the search constructor. The division is an exact monomial
division, and the new variable has only integer exponents. Evaluating at one
preserves the crossing endpoint.

The natural paired energies become \((0,10,6,2,8)\), with half those
heights on the two active legs at each place. Removing a common monomial
shifts the two profiles oppositely; the variable change rescales profiles,
physical heights, and the positive evaluation scale together. The complete
pressure is unchanged. The source's numerical ratio is approximately
\(0.7263479477\), still above \(1/\sqrt2\).

The degree reduction leaves 11,250 free field coefficients in (2).
The driver
[`search_f4_recentered_crossings.py`](../tools/search_f4_recentered_crossings.py)
records the complete reference tensors, both profiles, all physical heights,
row windows, completed iteration counts, and reference changes. A finite
budget or an interrupted schedule is distinguished from an exhaustive
search or a completed requested iteration count. Its pressure optimizations
are floating point. A strict numerical candidate still requires an exact
rational pressure certificate before any bound is claimed.

## 5. Independent controls

Run:

```sh
python3 tools/check_f4_recentered_lattices.py \
  --output /tmp/f4-recentered-lattice-controls.json
```

The [recorded controls](../evidence/f4-recentered-lattice-controls.json) cover
1,665 literal endpoint entries, 39 full tensor substitutions, 32 independent
Toeplitz-rank calculations, and four genuinely nonbinary reference changes.
Twenty additional complete-fiber queries agree with the predecessor binary
implementation. Four out-of-window requests are rejected.

For the Toeplitz controls, each \(F\)-coefficient is expanded into its
literal two-by-two binary multiplication matrix. The resulting binary rank
at precision \(P\) is compared with

\[
2\sum_i\max(0,P-\sigma_i).
\]

Thus these controls do not use the native \(F\)-profile implementation to
calculate both sides of the comparison. The evidence has the source hashes
and the precise small dimensions covered. It supports the construction and
finite implementation checks, not a new asymptotic factoring exponent.
