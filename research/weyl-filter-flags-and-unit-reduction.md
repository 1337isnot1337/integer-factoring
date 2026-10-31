# Polynomial filter flags and exact zero-unit reduction

## Result and scope

Two exact changes make the existing
[Weyl lattice synthesis](weyl-lattice-filter-synthesis.md) easier to
search. Constant occupation-block flags over \(\mathbf F_4\) allow
noncoordinate row-divisibility constraints. Common polynomial factors
that are units at zero can be removed before coefficient elimination.
Both changes preserve the physical endpoint construction and charge
the original alpha degrees at infinity.

These are constructive synthesis tools. The finite searches below
found no improved critical pressure. They do not exclude arbitrary
polynomial filters, variable flags, or other degree ranges.

## 1. Constant flags over a coefficient extension

Let \(M_{\boldsymbol\alpha}(x)\) be the reduced polynomial zero matrix
in the reference regular Gaussian Smith bases. It depends linearly
on the diagonal filter polynomials. Each block has fixed occupation.
For an invertible constant block matrix \(A\) over \(\mathbf F_4\),
replace the reference by \(A M_{\boldsymbol\alpha}\).
The matrix \(A\) and its inverse are regular at zero. Since all rows
in the block have the same aggregate charge, the comparison holds for
every value of that charge.

Prescribing

\[
\mathrm{ord}_x (A M_{\boldsymbol\alpha})_{ij}\ge t_i
\quad\text{for every }j
\tag{1}
\]

gives homogeneous linear coefficient equations over \(\mathbf F_4\).
For every minor of size \(k\), factor the corresponding row powers
\(x^{t_i}\). Its order is at least the sum of the \(k\) smallest
targets. Therefore the Smith order vector weakly supermajorizes the
sorted target vector. The partition-function comparison from the
Weyl lattice theorem applies to all Smith prefixes, rather than only
to the determinant.

In dimensions two and three, a complete flag can be represented by
an ordered constant basis. There are five full flags in
\(\mathbf F_4^2\) and 105 in \(\mathbf F_4^3\). The latter count is
21 projective lines times five planes through each line. The tool
enumerates those incidence pairs directly and checks invertibility.
When row targets are ordered, bases that differ by operations
preserving the associated flag impose the same filtered-row condition.
The search records still state their precise tested flag list; a
subset of flags is not reported as all flags.

The homogeneous solution space has at most eight endpoint coordinates
for the two-/three-mode binary engines. Over \(\mathbf F_{16}\),
nonvanishing of every coordinate projection is sufficient for an
everywhere-nonzero endpoint vector: eight proper hyperplanes cannot
cover a vector space over a field of size 16. The
[finite-field linear helper](../tools/f4_lattice_linear.py) constructs
the lift. The fixed physical-Y diagonal endpoint repair is regular at
both cuts, as proved in the
[repaired-endpoint note](repaired-endpoint-state-filters.md).

For the four-state case, \(\mathbf F_4\) already suffices whenever
all coordinate projections are nonzero. The implementation uses a
complete fallback at the equality boundary; it does not rely on the
greedy support-growth argument when the coordinate count equals the
field order. For eight states, the same \(\mathbf F_4\) coordinate
criterion alone would be invalid.

The implementation is
[`weyl_flag_filter_synthesis.py`](../tools/weyl_flag_filter_synthesis.py).
Its [independent checks](../evidence/weyl-flag-filter-synthesis.json)
include 36 known-Smith matrices over \(\mathbf F_{16}\), a degree-600
late-coefficient control, 30 native \(\mathbf F_4\) profile comparisons,
160 comparisons with the original identity-lattice solver, and 16
literal polynomial flag witnesses.

At prime-power interval lengths, the
[simultaneous identity-mode theorem](prime-power-joint-identity-filters.md)
provides a diagonal scalar-transform description of every joint
diagonal filter. Thus noncoordinate flags are not needed to compute
its exact profiles there. They remain a valid synthesis interface,
and can be useful when the reference family is not simultaneously
diagonal, including lengths just below a prime power.

## 2. Removing complete polynomial units

Write the reference source as

\[
M_{\boldsymbol\alpha}(x)=\sum_d\alpha_d(x)M_d(x).
\tag{2}
\]

Suppose a polynomial \(u(x)\), with \(u(0)\ne0\), divides every
source entry in one fixed row. Dividing that row in every \(M_d\) by
\(u\) is multiplication by a diagonal matrix invertible over the
local ring at zero. The same holds for a common factor in one fixed
column. These operations preserve every Smith order.

They also preserve the complete coefficient conditions (1), before
any new flag is applied. Multiplication or division by a zero unit
does not change the order of an individual entry. Equivalently,
convolution by the unit is invertible on every finite coefficient
jet. Thus each row-divisibility condition has exactly the same
solution space, not merely the same feasible dimension.

[`unit_reduced_lattice_basis.py`](../tools/unit_reduced_lattice_basis.py)
takes the greatest common divisor of all nonzero source entries in a
row or column, removes its monomial factor, and divides by the
remaining unit using complete polynomial long division. It alternates
rows and columns until no nonconstant common unit remains. Every
nontrivial division decreases the total source degree, so the process
terminates. No coefficient window or series cutoff is used.

This is a basis simplification at the zero cut. The literal filter
polynomials, their endpoint values, and their full infinity degrees
are unchanged. A new constant flag applied after this simplification
is a new valid regular reference frame; it need not coincide with
the same constant flag applied before the simplification.

The [unit-removal checker](../tools/check_unit_reduced_lattice_basis.py)
reconstructs every original source entry by restoring the recorded
row and column units. It also compares canonical rowspaces of the
full coefficient constraints and their deterministic literal lifts.
Its [record](../evidence/unit-reduced-lattice-basis.json) contains
225 complete polynomial identities, 40 full constraint-rowspace
comparisons, and 16 identical nonzero-endpoint witnesses.

| Interval and modes | Maximum source degree, before | After | Nonzero coefficient terms, before | After |
|---|---:|---:|---:|---:|
| \(n=3,r=2\), initial | 26 | 7 | 120 | 44 |
| \(n=5,r=3\), initial | 320 | 22 | 5,268 | 477 |
| \(n=31,r=3\), initial | 2,221 | 124 | 27,706 | 2,580 |
| \(n=32,r=3\), initial | 285 | 96 | 1,120 | 48 |
| \(n=4,P=3,r=2\), mixed | 87 | 24 | 88 | 14 |

## 3. Completed finite records

The [two-mode mixed flag screen](../evidence/weyl-f4-flag-mixed-n32-screen.json)
uses the two nonbinary flags in its one-particle block, 216 degree
vectors, and all 24 row-priority orders. It completed 325,704 exact
feasibility queries and retained 810 degree/flag/target combinations,
with 117 distinct actual profile pairs. All returned witnesses have
their complete profiles computed. The stated 132-point scale/charge
grid found no same-gauge gain below \(10^{-9}\) and no negative
critical pressure.

The [split three-mode lattice screen](../evidence/weyl-lattice-split-n32-screen.json)
completed 314,832 queries over 136 degree vectors, with 346 returned
targets. The [initial three-mode lattice screen](../evidence/weyl-lattice-three-mode-initial-n32-screen.json)
completed 1,082,674 queries over 512 degree vectors, with 4,270 returned
targets and 2,044 actual profile pairs. These finite grids also found
no gain. Their different scope from a complete coefficient catalogue
is explained in the [joint search note](joint-filter-search-records.md).

The latter record uses a lossless shared-profile representation.
[`compact_joint_lattice_screen.py`](../tools/compact_joint_lattice_screen.py)
restores every original row and verifies the expanded-row SHA-256
digest. Its 4,270 rows occupy about 1.62 MB; no searched row or
coefficient witness hash is discarded.

The [split \(\mathbf F_4\) monomial record](../evidence/f4-split-n32-screen.json)
contains all 369,603 stated exponent/color cases, grouped into 3,850
exact profile classes. Its complete 132-point grid found no gain.
The later all-index theorem and screens extend the tested physical
selections beyond this split geometry.

The recovered
[\(n=31\) initial three-mode screen](../evidence/weyl-lattice-three-mode-initial-n31-screen.json)
was rerun to completion using the unit-reduced basis. It tested 64
degree vectors with 16 row priorities, completed 47,524 exact
feasibility queries, and retained 192 targets. Its complete returned
witness profiles found no gain on the recorded grid. The minimum
sampled critical pressure was approximately 0.0367599061. This search
uses the initial-mode reduction, which is valid at \(n=31\); it does
not apply the prime-power all-indices theorem at a non-prime-power
length.

## 4. Reproduction and interruption scope

```sh
python tools/check_weyl_flag_filter_synthesis.py
python tools/check_unit_reduced_lattice_basis.py
python tools/search_weyl_flag_filters.py --n 32 --modes 2 --geometry mixed --flag-family nonbinary --field-size 4 --output /tmp/mixed-flag-screen.json
python tools/unit_reduced_lattice_basis.py --n 31 --modes 3 --geometry initial --degree-grid 0,1,31,32 --priority-limit 16 --output /tmp/initial-n31-screen.json
```

The last command applies the equivalent smaller source in its own
Python process and records the removed-unit count. The separate
[`search_unit_reduced_flag_filters.py`](../tools/search_unit_reduced_flag_filters.py)
entry point applies constant flags after that reduction.

An additional 24-flag \(n=31\) run was interrupted by the execution
environment outage, and its temporary output did not survive.
It is not part of the completed counts above. The accompanying
fixed-coordinate \(n=31\) screen was rerun after recovery and is
included in Section 3. No result from an
unretained temporary search is used as evidence.
