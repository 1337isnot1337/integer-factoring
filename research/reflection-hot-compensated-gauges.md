# Negative hot-state heights with compensation on every cold bit

## Result

The [four-mode hot-state evaluator](reflection-hot-state-gauges.md)
extends beyond nonnegative relative heights. At each cut, assign a
common height \(b_\nu\ge0\) to **every physical logical bit on
all four legs**. Add a state height \(h_{\nu,\ell}(a)\) on
the selected \(k\)-bit physical core. It is sufficient that

\[
h_{\nu,\ell}(a+2^j)-h_{\nu,\ell}(a)\ge-2b_\nu
\quad(a_j=0).
\tag{1}
\]

The weighted cut still reduces regularly to at most
\(\min(M,4k)\) active modes, with all other modes explicit.
For one selected bit, its total added height may be negative:
\(h(1)-h(0)=-2b_\nu\) gives a total increment \(-b_\nu\).
Every unselected physical bit receives the positive height
\(b_\nu\). These cold-bit heights are included in both the
Gaussian spectators and the complete physical type partition.

The construction applies to endpoint-valid exterior polynomial
defects over \(\mathbb F_4\). It does not assert a new factoring
exponent, and it does not permit an arbitrary coefficient mask.

## 1. Retain the internal terms supported on the selected core

At each cut, split the regular internal Gaussian on each boundary
into factors whose variables are both in the selected core and the
remaining factors. An individual factor has the form
\(1+q(z)C_iC_j\) on the output, or the corresponding input
factor. The original natural gauge makes \(q(z)\) integral.

Keep all factors involving two selected variables inside the active
tensor. Remove only the other internal factors. If both variables
are cold, the added height is \(2b_\nu\ge0\). If precisely
one is hot, its added height is

\[
2b_\nu+h(a+2^j)-h(a)\ge0.
\tag{2}
\]

Thus every removed factor remains regular after the new gauges.
Its inverse is regular as well: in characteristic two the square
of each nilpotent quadratic factor is one. Products of these
factors and their inverses remain regular.

Exterior multiplication commutes with every internal Gaussian
factor. Consequently the retained factors can be multiplied into
the exterior defect before applying the cross-only reduction.
This multiplication is performed separately at the two cuts.
The resulting local exterior factors represent the same original
polynomial tensor; they are not separate global repairs.

For one hot bit, write its four physical variables as \(X,Y,U,V\).
With ramification \(P\), the retained internal factors are

\[
\begin{aligned}
E_0&=(1+w^{-4P}(1+w^{2P})XY)
      (1+(1+w^{2P})UV),\\
E_\infty&=(1+(1+w^{2P})UX)
           (1+(1+w^{2P})VY).
\end{aligned}
\tag{3}
\]

The local factors used by the evaluator are \(H E_0\) and
\(H E_\infty\), where \(H\) is the original endpoint-valid
defect. Each factor in (3) has endpoint one. For larger selected
cores, the implementation extracts the relevant internal entries
from the original alternating kernel and multiplies their exact
quadratic factors. It retains every such entry, including entries
within one physical leg.

## 2. The same bounded cross reduction applies

After the remaining internal factors have been removed, the
underlying Gaussian cross matrix is the original cross matrix.
The [cross-only state reduction](reflection-hot-state-gauges.md),
Section 5, uses only ordinary linear changes of the one-particle
variables. Its complete physical creation and annihilation support
has at most \(\min(M,4k)\) active modes.

These changes preserve total occupation. Their regular physical
core coordinate bases also preserve it: the joint vacuum projectors
preserve degree, and a core creation operator increases degree by
one. A coordinate indexed by core state \(a\) and vacuum row or
column \(r\) therefore has occupation

\[
|a|+|r|.
\tag{4}
\]

The common bulk height is exactly a row or column weight
\(b_\nu(|a|+|r|)\) in these coordinates. It need not be
inserted into a polynomial coefficient array or approximated by a
series. Arbitrary real heights are useful for exploration; rational
heights acquire their ordinary valuation interpretation after a
fixed ramification.

The noncore Gaussian factor with original order \(s_j\) has
one occupied variable on each side of a cut. Its new order is

\[
s_j+2b_\nu.
\tag{5}
\]

Equations (4)--(5), together with the relative core-state weights,
give the full weighted profile. They also show why dropping the
cold modes would give an incorrect cost.

## 3. Complete physical partition and degree accounting

Let \(g_i\) be the original paired-type energy of logical bit
\(i\). Every cold bit has new energy

\[
g_i+2(b_0+b_\infty).
\tag{6}
\]

Indeed, each paired physical letter appears on two legs at each
cut. For a selected core state \(a\), sum (6) over its occupied
bits, then add the original paired sums of the state heights:

\[
\begin{aligned}
g_X(a)&=\sum_{i\in a}[g_i+2(b_0+b_\infty)]
 +\sum_\nu[h_{\nu,X}(a)+h_{\nu,V}(a)],\\
g_Y(a)&=\sum_{i\in a}[g_i+2(b_0+b_\infty)]
 +\sum_\nu[h_{\nu,Y}(a)+h_{\nu,U}(a)].
\end{aligned}
\tag{7}
\]

All cold factors from (6) and all selected states from (7) enter
the physical partition. The numerical evaluator differentiates
both the active and spectator profile costs and these complete
physical energies when optimizing the two bulk heights.

The ordinary polynomial seed is still the original compressed
Gaussian multiplied by \(H\), with one common monomial clearing
if needed. The two local factors in (3) incur no new global degree:
they only record factors already present in that Gaussian. Its
common scalar and the clearing scalar shift the cuts oppositely
and cancel from the profile product.

## 4. Implementation and checks

### Independent boundary bulks

The same proof allows separate nonnegative bulks \(b_{\nu,L}\)
and \(b_{\nu,R}\) on the two boundaries of each cut. Condition
(1) uses the bulk of the boundary containing that physical leg.
The active matrix receives the row weight
\(b_{\nu,L}(|a|+|r|)\) and the column weight
\(b_{\nu,R}(|a|+|r|)\). Equation (5) becomes

\[
s_j+b_{\nu,L}+b_{\nu,R}.
\tag{8}
\]

Each paired physical letter contains one leg from each boundary at
each cut. Thus the bulk term in (6)--(7) becomes
\(\sum_\nu(b_{\nu,L}+b_{\nu,R})\) per occupied bit.
No cold physical bit is omitted. This refinement can compensate
negative hot heights on just one boundary at a cut.

### Reproducibility

The constructor is
[`reflection_hot_compensated_gauges.py`](../tools/reflection_hot_compensated_gauges.py).
It records the original and the two local exterior factors,
preserves the complete coefficient and precision metadata, and
rejects a height outside (1). The active matrix retains its two
parity blocks.

The independent checker
[`check_reflection_hot_compensated_gauges.py`](../tools/check_reflection_hot_compensated_gauges.py)
builds the complete original tensors, adds the bulk height to
every physical occupied bit, and compares both original weighted
cuts with the reduced calculation. It includes negative total hot
heights, field-valued Laurent defects, a nontrivial ramification,
and a nonlinear two-bit state gauge. It also compares the pressure
with a complete physical-alphabet calculation and checks the bulk,
state-height, and scale derivatives by finite differences.

```sh
python tools/check_reflection_hot_compensated_gauges.py --extended \
  --output evidence/reflection-hot-compensated-gauges.json
```

The exact control record is
[`reflection-hot-compensated-gauges.json`](../evidence/reflection-hot-compensated-gauges.json).
The checker includes both uniform and independent-boundary bulks,
along with regressions against the earlier positive cone. Its saved
record gives the exact profile and endpoint counts and the errors
from the finite-difference derivative controls.

The search driver is
[`search_reflection_hot_compensated_gauges.py`](../tools/search_reflection_hot_compensated_gauges.py).
For one hot bit it parameterizes (1) by nonnegative slacks
\(u_{\nu,\ell}\), setting \(h(1)-h(0)=u_{\nu,\ell}-2b_\nu\).
Its finite numerical searches retain the original coefficient row,
both bulk heights, all relative and total hot heights, and the
complete cold physical energies in each saved best point.
