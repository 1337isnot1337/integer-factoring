# Polynomial synthesis through reduced Smith lattices

## Result and scope

Polynomial joint filters can now be synthesized with several independent
Smith-prefix targets. The method fixes a regular basis for the small
Weyl-reduced matrix, then imposes a separate valuation on each row.
These are linear coefficient conditions. Their solutions certify every
Smith prefix in each occupation block, and therefore an upper bound on
the complete zero-cut partition function.

The current implementation treats two or three binary identity
directions and permits `F4` coefficients with arbitrary nonzero
occupation endpoints. It constructs the diagonal physical endpoint
repair explicitly through the existing tensor interface. For eight
occupation states, it uses an exact endpoint-image test; it does not
assume that nonzero coordinate projections alone suffice over `F4`.

This supplies a constructive search method, not a new factoring rate.
The recorded finite searches have not found an improved critical
pressure. An infeasible row-lattice target excludes only that target in
that fixed basis. It is not an exclusion of every tensor with the same
Smith-prefix bounds.

## 1. A small polynomial matrix linear in all filter coefficients

Let `r` be two or three, `n>=r`, and let `P` be positive and odd. Use one
of the selections in
[the coupled identity theorem](coupled-identity-filter-reduction.md):
initial directions for every `n`, terminal directions when `n` is a
power of two, or the two-direction mixed selection when `n` is a power
of two. The [split-cluster theorem](split-identity-mode-reduction.md)
also permits `split:a:b` with `a+b=r` at power-of-two lengths. Write
the selected indices as `j_1,...,j_r` and put `x=-v`.

The [Weyl theorem](identity-operator-weyl-reduction.md) reduces the zero
cut to a `2^r` by `2^r` matrix, with occupation blocks of sizes
`binom(r,d)`. After multiplying by its known zero-unit denominator and
including the `J=1` row gauge, this matrix has the form

\[
M_\alpha(x)=\sum_{e\in\{0,1\}^r}\alpha_e(x)M_e(x),
\qquad M_e(x)\in\mathbb F_2[x]^{2^r\times2^r}.
\tag{1}
\]

The matrices `M_e` preserve total occupation. The unknown filters are
polynomials in `F4[x]` with

\[
\alpha_e(1)\ne0.
\tag{2}
\]

The geometric coefficients in (1) lie in `F2`, even though the unknown
filter coefficients may lie in `F4`. The common denominator clearing
has nonzero constant term, so it changes no zero order. It is an
auxiliary basis operation; it does not alter the degree of the original
physical filter.

By the [diagonal endpoint repair](repaired-endpoint-state-filters.md),
(2) gives the original crossing endpoint after a constant invertible
map on the physical `Y` alphabet. This repair preserves both weighted
Smith profiles and ordinary polynomial degree.

## 2. Compute an exact regular reference basis

Set every `alpha_e` equal to one, giving the reduced Gaussian matrix
`M_G`. Work separately in its occupation blocks. Exact polynomial
row and column operations construct matrices `U_d,V_d` such that

\[
U_d M_{G,d}V_d
=\mathrm{diag}(x^{g_i}u_i(x)),\qquad u_i(0)\ne0,
\tag{3}
\]

and `U_d,V_d` and their inverses are regular at zero.

The implementation uses a fraction-free elimination. Choose a
remaining entry with smallest valuation, move it to the pivot, and
write it as `x^h u`, where `u(0)=1`. Every other entry in the remaining
block is divisible by `x^h`. To clear a pivot-column entry `x^h a`,
replace its row by

\[
u\,\mathrm{row}_i-a\,\mathrm{row}_{\rm pivot}.
\tag{4}
\]

This is a polynomial operation with determinant `u`, a zero unit. Its
inverse is regular. The corresponding column operation clears the
pivot row. Iterate on the remaining block. The pivot orders cannot
decrease, and the nonsingularity of the Gaussian block guarantees
termination. Keeping the operations gives exact polynomial `U_d,V_d`.

The search reuses this basis for every polynomial filter. Define

\[
T_{\alpha,d}=U_dM_{\alpha,d}V_d.
\tag{5}
\]

Each entry of (5) is still linear in the unknown coefficients of the
`alpha_e`. The basis is fixed before imposing their degree caps or
valuation targets.

## 3. Row targets certify every Smith prefix

Choose an integer target `h_i` for every row of each occupation block,
and impose

\[
\mathrm{ord}_0(T_{\alpha,d})_{ij}\ge h_i
\quad\hbox{for every column }j.
\tag{6}
\]

Negative row targets are allowed. Let `h_(1)<=...<=h_(m)` be the
sorted targets in a block of size `m`, and let
`sigma_1<=...<=sigma_m` be its actual Smith orders in the variable `x`.
Factoring `x^{h_i}` from the selected rows of any minor gives

\[
\sum_{i=1}^k\sigma_i\ge\sum_{i=1}^k h_{(i)}
\qquad(1\le k\le m).
\tag{7}
\]

The row and column bases preserve the Smith orders, so these are
inequalities for the original reduced block.

Equation (7) implies the partition bound

\[
\sum_i s^{-2\sigma_i}\le\sum_i s^{-2h_i}
\qquad(s>1).
\tag{8}
\]

Here the factor two comes from `x=z^2`. One direct proof uses

\[
\sum_i(t-\sigma_i)_+
=\max_k\left(kt-\sum_{i=1}^k\sigma_i\right).
\]

By (7), this is no larger than the corresponding expression for `h`.
Integrating against `gamma^2 exp(-gamma t)` on the real line gives
(8), with `gamma=2 log s`. No equality of total determinant orders is
required.

For an auxiliary binary occupation `a`, an order in reduced occupation
`d` gives the full active zero order

\[
2\sigma_i-2P\sum_\ell j_\ell a_\ell
 +P(1-J)(d-|a|).
\tag{9}
\]

The last two terms are common shifts for that copy of the occupation
block. Thus (8) remains valid after reconstruction, at every aggregate
charge `J`. Gaussian spectator modes multiply the resulting bound.
The function `bound_blocks` performs this reconstruction with `h_i` in
place of `sigma_i`.

This is a sufficient certificate for many prefix inequalities. Direct
minor conditions would be polynomial, rather than linear, in the
unknown filter coefficients. The fixed row lattice keeps the search
linear, at the cost of potentially missing a useful matrix in a
different regular basis.

## 4. Charge the full infinity degree

Supply nonnegative degree caps `D_e` and require
`deg alpha_e<=D_e`. The active infinity order of state `e` is exactly

\[
P\sum_\ell(2j_\ell+1-2n)e_\ell-2\deg\alpha_e.
\tag{10}
\]

Replacing the actual degree by its cap gives an upper bound on the
infinity partition cost. Combining (8)--(10), the exact Gaussian
spectators, and the original physical partition gives a complete
two-cut cost bound for any feasible returned witness.

A common monomial clearing shifts the two cuts oppositely and cancels
from their product. A nonmonomial rational denominator still has its
degree-minus-zero-order cost. The algorithm searches the cleared
polynomials, so it does not omit that cost. The basis polynomials in
(3) are only used to certify the zero cut; they are not multiplied into
the physical crossing.

## 5. Exact endpoint feasibility over F4

Expand each `alpha_e` up to its cap. Setting the forbidden coefficients
in (6) to zero gives a homogeneous matrix `A` over `F2`. Let

\[
E=\{(\alpha_e(1))_e:A\alpha=0,\ \alpha\in\mathbb F_2[x]^S\}
 \subseteq\mathbb F_2^S,
\qquad S=2^r.
\tag{11}
\]

The endpoint image over `F4` is exactly `E+omega E`.

### Four states

For `r=2`, a full-support point exists in this extension if and only
if every coordinate projection of `E` is nonzero. At most four proper
coordinate hyperplanes cannot cover a vector space over a field with
four elements. Constructively, add scalar multiples of endpoint basis
vectors while retaining every already nonzero coordinate. Unless a
basis vector is already fully supported, at most three scalar values
are forbidden at the next step. The coefficient lifts give the actual
`F4` polynomials.

### Eight states

For `r=3`, the coordinate-projection test alone is insufficient.
Instead use the exact criterion

\[
\exists v,w\in E:\quad
\mathrm{supp}(v)\cup\mathrm{supp}(w)
=\{1,\ldots,8\}.
\tag{12}
\]

The equivalence is immediate from the unique expansion
`lambda=v+omega w`: a coordinate is zero exactly when both binary
coordinates are zero. The space `E` has at most 256 vectors. Enumerating
them and tabulating which supports have a superset in `E` decides (12)
exactly. The two stored coefficient lifts construct the polynomial
witness.

There is a small counterexample to replacing (12) by coordinate
support. Use the seven nonzero binary linear forms on three variables
as seven endpoint coordinates, and repeat the first form as an eighth
coordinate. Every projection is nonzero. But two vectors in the
three-dimensional binary input space have a common nonzero orthogonal
linear form, so their two images miss a coordinate. The corresponding
`F4` endpoint image has no full-support point. The checker tests all 64
`F4` inputs in this example.

### Coefficient elimination and witness construction

The implementation stores homogeneous coefficient rows as binary
integers. A coordinate projection is zero exactly when its endpoint
functional lies in the homogeneous row space. Reducing those
functionals also gives the small image (11). Coefficient tables are
cached for each degree-cap vector; changing a row target merely adds
or removes an initial segment of its coefficient equations.

Every successful synthesis is replayed against the literal transformed
polynomial entries in (5), the original degree caps, and its nonzero
endpoint values. The returned `alpha` values are actual finite
polynomials, not only a rank certificate.

## 6. A strict endpoint-feasibility enlargement

For two initial directions at `n=2,P=1`, the Gaussian row orders in
the implementation's occupation-block order are `(4,2,3,1)`. Take

\[
D=(0,0,0,1),\qquad h=(1,0,1,2).
\tag{13}
\]

The homogeneous constraints reduce to
`alpha_11(0)=0` and
`alpha_00(0)+alpha_01(0)+alpha_10(0)=0`. They admit

\[
(\alpha_{00},\alpha_{01},\alpha_{10},\alpha_{11})
=(1+\omega,\omega,1,x).
\tag{14}
\]

All endpoint values are nonzero. In contrast, requiring every endpoint
to equal one makes the first three constant polynomials equal one,
whose sum in characteristic two is one. Thus the unrepaired
endpoint-one target is infeasible. This proves an exact enlargement
of the feasible search region. The example is not a pressure witness.

## 7. Verification and finite search records

The implementation is
[the reduced-lattice synthesizer](../tools/weyl_lattice_filter_synthesis.py).
Its [first checker](../tools/check_weyl_lattice_filter_synthesis.py) and
[record](../evidence/weyl-lattice-filter-synthesis.json) include six
exact regular reference-basis controls, 41,472 complete small `F4`
polynomial tuples, and 4,864 exact degree/order comparisons with that
complete catalog. Eight returned polynomial witnesses are replayed in
both original physical cuts, with all 2,048 endpoint coefficients
checked. The checker also verifies (13)--(14).

The [endpoint-cover checker](../tools/check_weyl_lattice_endpoint_cover.py)
and [its record](../evidence/weyl-lattice-endpoint-cover.json) compare
240 endpoint maps with 70,236 complete `F4` assignments. They replay
five three-direction polynomial witnesses in the original matrices,
including `split:1:2`, through size `256` by `256`, and check all
143,360 endpoint coefficients.
Every reconstructed Smith-prefix bound is checked against the exact
active profile before the physical replay.

The [search driver](../tools/search_weyl_lattice_filters.py) raises row
targets in a recorded finite list of priority orders. For each distinct
feasible degree/row pair, it constructs a witness and computes its exact
zero and infinity profiles. This last step can detect a gain hidden by
a loose sufficient row bound.

The [mixed two-direction record at n=32](../evidence/weyl-lattice-mixed-n32-screen.json)
contains 512 degree-cap vectors, 384,776 exact feasibility queries,
and 998 distinct feasible degree/row pairs. The
[initial three-direction record at n=8](../evidence/weyl-lattice-three-mode-initial-n8-screen.json)
contains 64 occupation-class degree-cap vectors, 99,043 exact queries,
and 214 distinct feasible degree/row pairs. Neither record finds a
same-gauge gain or negative critical pressure in its returned
witnesses at tolerance `10^-9` on the documented grids. These searches
do not enumerate every feasible coefficient choice or every Pareto
target.

Reproduce the exact checks and finite searches from the repository root:

```bash
python tools/check_weyl_lattice_filter_synthesis.py \
  --output evidence/weyl-lattice-filter-synthesis.json
python tools/check_weyl_lattice_endpoint_cover.py \
  --output evidence/weyl-lattice-endpoint-cover.json
python tools/search_weyl_lattice_filters.py \
  --n 32 --geometry mixed \
  --output evidence/weyl-lattice-mixed-n32-screen.json
python tools/search_weyl_lattice_filters.py \
  --n 8 --modes 3 --geometry initial --degree-grid 0,1,8,9 \
  --output evidence/weyl-lattice-three-mode-initial-n8-screen.json
```
