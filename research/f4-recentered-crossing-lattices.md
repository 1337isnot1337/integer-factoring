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
\operatorname{coeff}_{z^k}\bigl(L A_e(T)R\bigr)_{ij}=0
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

For dense nonbinary frames, repeated elimination of the same pivot coordinates
was the main computational cost. The accelerated implementation expresses
each original unknown in free affine coordinates by triangular substitution.
It caches the 256 binary sums within each block of eight coordinate
expressions. The two coefficient planes use the same tables; multiplication
by \(\omega\) is applied to the second sum. The affine constant is an extra
coordinate. These are exact substitutions, with the same output equations
as direct elimination. A memory estimate selects direct elimination when
the cached tables would exceed the stated implementation budget.

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

gives a degree-18 strict crossing on five symbols. The standalone search
constructor evaluates the corresponding small Laurent-Pfaffian recurrence
directly in \(z\). Its 625 entries were compared exactly with (7), retaining
the source digest
`91e0ee1014b79e5e042053201717ccd95ad33a27ec7f1ee9a6c2144b10aa75e9`.
The division in (7) is an exact monomial
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

## 5. Completed search from a different source family

The driver also supports the five-symbol tensor

\[
R_{xyuv}=[x=y][u=v],\qquad E_{xyuv}=[x=v][y=u],
\qquad B(z)=R+z(E+R).
\tag{8}
\]

It has degree one and satisfies the strict endpoint \(B(1)=E\).
This reference is independent of the reflection Gaussian source. The
completed run allowed the whole degree-six endpoint fiber, used paired
energies \((0,12,6,0,8)\) at both places, set row floors three below the
reference orders, and set caps four above them.

| Reference | Completed steps | Distinct target vectors | Profile classes in that reference | Complete profile queries |
| --- | ---: | ---: | ---: | ---: |
| Binary Brauer source | 200 | 155 | 278 | 465 |
| Retained nonbinary F4 tensor | 200 | 8 | 16 | 24 |

The common-floor fibers had 1,770 and 1,501 free F4 coordinates. The best
recorded ratio moved from the original reference's approximately
\(0.9957096676\) to \(0.9825257604\). This improves this deliberately
different finite source only; it is far above the existing factoring
threshold and does not improve the repository's best bound.

The [complete record](../evidence/f4-recentered-brauer-five-search.json) has
both reference tensors, all physical heights, and 24 retained witnesses.
Its [independent replay](../evidence/f4-recentered-brauer-five-replay.json)
certifies 52 complete profiles by exact polynomial regular frames, verifies
48 untruncated target-frame products, and checks 16,250 literal endpoint
entries. This replay does not call the native Smith profile oracle or solve
the search's affine systems. The numerical roots are recomputed separately.

```sh
python3 tools/check_f4_recentered_walks.py \
  evidence/f4-recentered-brauer-five-search.json \
  --output /tmp/f4-recentered-brauer-five-replay.json
```

The search driver can restart a new finite schedule from a selected retained
literal witness with `--initial-record` and `--witness-index`. This is a new
schedule from a saved tensor, not a claim to resume the exact earlier
iteration state.

### Restarts from the complete allowed row floor

The driver accepts `--restart-every K`. At every positive multiple of `K`,
it begins the greedy row search from all of the declared floor targets at
once. The default `K=0` retains the earlier local walk. This releases a
possible search bottleneck: lowering only a subset of the current rows can
repeatedly recover the same reference even when other parts of the allowed
fiber are feasible. A restart still stays within the exact common-floor
fiber, the same endpoint equations, and the same degree bound. Its counter
is retained for each reference. This is a search change, not a claim of
exhaustiveness or an improved pressure bound.

A focused two-reference control uses the degree-two five-symbol Brauer
fiber, eight steps per reference, and a restart every two steps. It completes
six global restarts and thirteen complete profile queries, including a
nonbinary reference. The [saved control](../evidence/f4-global-restart-controls-search.json)
and [independent replay](../evidence/f4-global-restart-controls-replay.json)
certify thirty complete profiles, twenty-six untruncated target products,
and 9,375 literal endpoint entries. Run `make verify-f4-global-restarts`.
The numerical source-family ratio is still above 0.99 and supplies no
factoring improvement.

## 6. Independent controls

Run:

```sh
python3 tools/check_f4_recentered_lattices.py \
  --output /tmp/f4-recentered-lattice-controls.json
```

The [recorded controls](../evidence/f4-recentered-lattice-controls.json) cover
1,665 literal endpoint entries, 39 full tensor substitutions, 32 independent
Toeplitz-rank calculations, and four genuinely nonbinary reference changes.
Twenty additional complete-fiber queries agree with an independent literal
binary Smith-jet construction. Four out-of-window requests are rejected. The accelerated
coordinate substitution has another 100 literal substitution controls,
including dense triangular systems with 257 variables. Four complete
crossing fibers agree equation by equation with direct elimination.

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
