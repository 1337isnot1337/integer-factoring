# Shared letter gauges for coherent physical maps

## Result and scope

Both complete cuts of a homogeneous coherent retraction admit a common
polyhedral family of letter gauges. The conditions use the actual induced
bit-deletion operators of both physical map spaces. They permit some
negative letter potentials. They also account for the intrinsic positive
valuation of a within-cut pairing, so they are broader than requiring
every deletion to decrease the potential.

The gauge changes the complete cut profiles and the physical type cost
together. The resulting pressure equation is given in section 3. No
factoring exponent improvement is asserted here. The implementation uses
real weights for discovery; a proposed improvement needs a rational gauge
and an exact pressure certificate.

The exact controls compare six rationally gauged complete profiles with
the literal ambient Pfaffian tensor. They include both a nonorthogonal
eight-letter map and a mixed-occupation nonlinear quotient map.

## 1. Induced deletion edges and the regularity cone

Let `P,L` satisfy the hypotheses of the
[homogeneous coherent construction](coherent-hot-mixed-boundary.md):
`LP=I`, paired labels have common heights `h_j`, and both map spaces are
closed under every bit deletion. For each physical mode `a`, form both
induced operators

\[
M^P_a=LD_aP,\qquad M^L_a=P^{\mathsf T}D_aL^{\mathsf T}.
\tag{1}
\]

An edge `(a,j,k)` means that at least one of these matrices has a nonzero
entry from label `j` to label `k`. Each edge satisfies

\[
h_j-h_k=2a+1.
\tag{2}
\]

Choose a power-of-two reflection interval `M >= 4s`, where `s` is the
number of selected hot modes. Then the selected within-leg reflection
links vanish. The only nonzero within-row or within-column Gaussian
pairings in the stable hot kernel have the form

\[
A_{ab}(t)=F_{ab}(t)=(1+t)t^{2(a\mathbin{\&}b)}.
\tag{3}
\]

The infinity kernel exchanges the complete second and third physical
roles. This gives the same pair coefficient for its two within-cut
factors. The two map spaces may exchange roles, which is why (1) includes
both induced deletion operators.

Give encoded letter `j` a gauge exponent `theta_j`. A sufficient condition
for both conjugated Gaussian factors in both cuts to remain regular is

\[
(\theta_j-\theta_k)+(\theta_l-\theta_m)
 +2(a\mathbin{\&}b)\ge0
\tag{4}
\]

for every pair of induced edges `(a,j,k)` and `(b,l,m)`. Repetition of an
edge is allowed. This is a sufficient cone; cancellations can make the
exact regularity region larger.

To prove (4), first consider one elementary Gaussian pairing. Its induced
coefficient is a product of the two deletion entries and (3). Conjugating
the row factor by the diagonal letter gauge adds the two potential
differences in (4). Conjugating the column factor in the opposite order
adds those same differences, because the column factor has the opposite
triangular orientation. Thus every elementary factor is polynomial in a
common ramified parameter. A product of these factors remains regular.
Each has identity diagonal and is triangular in height, so its inverse
is regular as well. This proves the full two-cut descent.

If every paired physical label also has one common occupation number,
then (4) is necessary for these four descended factors. Consider an entry
that lowers occupation by exactly one on each of its two legs. Only one
single-mode deletion on each leg can contribute. Equation (2) uniquely
identifies those two modes. Its coefficient is the nonzero product of
the two induced deletion entries and (3), so its valuation is exactly
`2(a & b)` and cannot cancel. Every combination of the two map spaces
occurs among the zero and infinity factors. Thus any violated inequality
creates a pole in one of them. For occupation-homogeneous maps, (4)
therefore describes the exact shared regularity region for these factors.
For mixed-occupation maps, higher odd deletions can share the same height,
so this necessity argument does not apply.

For a numerical scale `u > 0`, the implementation stores `f_j=u theta_j`.
Its cone is the homogeneous linear system

\[
(f_j-f_k)+(f_l-f_m)+2u(a\mathbin{\&}b)\ge0.
\tag{5}
\]

The vacuum potential is set to zero. A common shift of all potentials
cancels from the pressure equation, so this only fixes a normalization.

### Negative potentials and a necessary lower bound

Repeating an edge in (5) gives `f_j-f_k >= -u a`. A nonzero homogeneous
physical vector has a sequence of derivatives ending in the vacuum.
In the induced matrices, at least one path contributing to that sequence
has nonzero coefficients. Along a path deleting `n` modes, (2) and the
repeated-edge inequalities give

\[
u h_j+2f_j\ge un>0
\tag{6}
\]

for each nonvacuum label. The physical energies in section 3 are therefore
positive throughout the feasible cone. They need not have nonnegative
individual potentials.

A uniform occupation potential is more restrictive. If the two map
spaces retain the vacuum and the first single mode, the pairing between
mode zero on two different legs has unit coefficient `1+t`. Condition
(4) then forces the uniform occupation coefficient to be nonnegative.
The two distinct physical legs are essential: pairing a mode with its
copy on another leg is allowed.

## 2. Complete weighted profiles

Use the exact compound representatives

\[
C_0=(L\otimes L)E_B(P\otimes P),\qquad
C_\infty=(L\otimes P^{\mathsf T})E_B(L^{\mathsf T}\otimes P).
\]

Multiply both rows and columns by the paired-letter gauge. For a row or
column labeled `(j,k)`, its dimensionless weight is `f_j+f_k`. A matrix
coefficient with natural `t`-valuation `a` consequently has initial weight

\[
ua+f_j+f_k+f_l+f_m.
\]

The complete weighted Smith orders, after all cancellations, are denoted
`lambda_(0,i)` and `lambda_(infinity,i)`. There are exactly `q^2` orders
in each cut. Diagonal letter gauges preserve every compound occupation
block, so the implementation uses those smaller blocks when every label
has one occupation number. When labels mix occupations, height still
fixes parity and the two complete parity blocks are used instead.

The [oracle](../tools/coherent_weighted_gauges.py) constructs every exact
binary polynomial coefficient once. It then asks the weighted Smith
engine for every generic pivot. The polynomial endpoint proves full
generic rank independently of the series window. A missing pivot causes
a larger workspace or an explicit error; it never produces a truncated
profile presented as complete.

Within a fixed pivot region, each returned weighted order is affine in
`u` and the `f_j`. The selected row and column labels give its derivative.
The implementation checks that these labels and the integer base order
reconstruct every returned weight. At pivot ties the objective can be
nondifferentiable, so the reported derivative is only the derivative of
the selected region. Numerical optimization has no global guarantee.

## 3. Physical cost and the pressure equation

Each cut uses four letter gauges. Using the same gauge at both places
therefore adds a total physical type cost `8 sum_j p_j f_j`. The natural
height cost is `4u sum_j p_j h_j`. Define

\[
e_j=uh_j+2f_j.
\tag{7}
\]

The balanced pressure equality is

\[
\log\sum_i e^{-\lambda_{0,i}}
+\log\sum_i e^{-\lambda_{\infty,i}}
=4\rho\log\sum_j e^{-e_j/\rho}.
\tag{8}
\]

Both the profile weights and the right-hand side must change when a
gauge changes. Holding the natural physical heights fixed would omit the
cost of the gauge. Setting every `f_j` to zero recovers the established
natural pressure formula.

For clarity, let `L` denote the left side of (8), let
`p_j=exp(-e_j/rho)/sum exp(-e/rho)`, and let `H(p)` be its entropy. In a
fixed pivot region, implicit differentiation gives

\[
\frac{\partial\rho}{\partial f_j}
=\frac{\partial L/\partial f_j+8p_j}{4H(p)},
\qquad
\frac{\partial\rho}{\partial u}
=\frac{\partial L/\partial u+4\sum_jp_jh_j}{4H(p)}.
\tag{9}
\]

The vacuum normalization removes one potential coordinate from the
optimizer. The derivative in (9) is used only for discovery. An exact
result must choose rational gauge exponents, clear their denominators,
recompute all integral Smith orders, and certify a strict pressure
inequality with exact arithmetic. For a fixed finite alphabet, this is
an effective finite computation. The polynomial degrees, denominator
ramification, and all matrix entries remain part of its cost.

## 4. An exact negative-gauge symmetry

Let `c` be a nonnegative integer mask and put

\[
J_a=(a\mathbin{\mathrm{xor}}c)-a.
\]

For both the hot cross kernel and the within-cut pairing,

\[
t^{J_a+J_b}B_{ab}(t)
 =B_{a\mathbin{\mathrm{xor}}c,\,b\mathbin{\mathrm{xor}}c}(t),
\]

\[
t^{J_a+J_b}A_{ab}(t)
 =A_{a\mathbin{\mathrm{xor}}c,\,b\mathbin{\mathrm{xor}}c}(t).
\tag{10}
\]

These identities follow bit by bit. If both bits are zero, the exponent
increases by twice the bit value; if both are one, it decreases by that
amount; and if they differ, it stays fixed. The diagonal cross-side
entry transforms by `2a+1+2J_a=2(a xor c)+1`.

Whenever `sum_(a in A) J_a` is constant on each coherent letter support,
this gives an allowed encoded-letter gauge. If the physical mode set is
closed under the XOR permutation, it is an exact relabeling of the
kernel. In particular it preserves the full Fock pressure after relabeling
the physical states. This is a useful control for negative-gauge searches.
The continuous full-letter cone can contain gauges that have no additive
mode representation.

## 5. Independent controls

The [checker](../tools/check_coherent_shared_gauges.py) first verifies
3,264 entries of the literal stable zero and infinity kernels, including
(3) and the full role swap. It then uses three physical maps: a six-letter
inverse system, a nonorthogonal eight-letter code, and the nonlinear
eight-letter quotient with `phi(x_4)=z_0 z_1 z_2`.

For each map it chooses `u=1/2` and explicit potentials with denominator
32, including negative values. It reconstructs the complete ambient
Pfaffian coefficients, substitutes the ramified variable, applies the
integer letter shifts, and uses independent binary-polynomial Smith
elimination. All six complete cut profiles agree with the weighted
compound oracle: 328 Smith orders and 18,976 literal tensor entries.

At nearby strict-interior points, eighteen central finite differences
also agree with (9) to absolute error below `10^-6`. These derivative
checks are numerical; the preceding rational profile comparisons are
exact. The [retained control record](../evidence/coherent-shared-gauge-controls.json)
contains all rational potentials, complete integral orders, error values,
and source hashes.

The separate [search driver](../tools/search_coherent_shared_gauges.py)
retains the complete cone, maps, best complete profiles, every finished
optimizer result, and its explicit budget. It does not label a floating
point root as an exact rate certificate or a finite search as a global
optimum.
