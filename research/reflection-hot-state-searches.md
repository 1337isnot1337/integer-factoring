# Finite reflection hot-state search checkpoint

## Outcome and scope

These completed finite screens found no improvement over the existing
Gaussian control. They use the proved
[bounded physical-state construction](reflection-hot-state-gauges.md),
retain all Gaussian spectators and the complete physical alphabet,
and report the actual optimizer statuses. They do not exclude other
coefficient polynomials, untested heights, or better local minima.

The natural Gaussian controls are
\(\rho\approx0.710604404520\) at \(M=16\) and
\(\rho\approx0.716653659728\) at \(M=8\).
The target for an improved factoring exponent remains
\(\rho<1/\sqrt2\approx0.707106781187\).

## Completed screens

| Family | Completed scope | Best nonconstant value found |
|---|---|---:|
| One hot bit, M16, positive state heights | Seven selected rows, 64 random rows, and the Gaussian; three starts per row | 0.771282886066 |
| One hot bit, M16, positive state heights | All 234 representatives of the natural-profile classes; two starts per row | 0.758070131461 |
| One hot bit, M8, unrestricted state heights | Seven selected rows and the Gaussian; two starts per row | 0.851761145122 |
| Two hot bits, M16, nonlinear state heights | Ten selected single-defect tensors and the Gaussian; 36 integer height pairs each | 0.720070222525 |

In every row of this table, the Gaussian control was the best tensor
in the complete recorded screen. The numerical search budgets are
part of the scope. A natural-profile class need not remain one class
after changing the physical state heights; the 234-representative
screen therefore does not stand for all 65,536 coefficient rows at
arbitrary heights.

The records are:

- [Selected and random positive-height rows](../evidence/reflection-hot-positive-gauges-M16.json).
- [All natural-class representatives](../evidence/reflection-hot-positive-classes-M16.json).
- [Unrestricted M8 heights](../evidence/reflection-hot-unrestricted-gauges-M8.json).
- [Nonlinear two-bit heights](../evidence/reflection-hot-joint-heights-M16.json).

The one-bit driver is
[`search_reflection_hot_state_gauges.py`](../tools/search_reflection_hot_state_gauges.py).
The eight relative bit heights and logarithmic scale are all optimized.
The positive branch uses four active modes; the unrestricted M8
controls use all eight active modes and their two parity blocks.

## All coefficient rows at five positive height rays

A separate compiled sweep checked every row in the prescribed
\(\mathbb F_4^8\) coefficient grid at five fixed integer-height
choices. The row ordering and coefficient polynomials are exactly
those in
[`search_reflection_hot_core_filters.py`](../tools/search_reflection_hot_core_filters.py).
The following tuples list the relative occupied-state heights in
physical order \((X,Y,U,V)\):

| Ray | Zero cut | Infinity cut | Distinct profile pairs |
|---|---|---|---:|
| 1 | (2,2,2,2) | (0,0,0,0) | 228 |
| 2 | (0,0,0,0) | (2,2,2,2) | 249 |
| 3 | (2,2,2,2) | (2,2,2,2) | 243 |
| 4 | (2,0,0,2) | (2,0,0,2) | 424 |
| 5 | (0,2,2,0) | (0,2,2,0) | 424 |

Each candidate and ray received both complete active profiles, for
655,360 cut profiles in total. All selected heights are integers,
and the Smith orders are checked to be integers. The constructor
first replays 30 compiled profiles against the uncompiled state
evaluator. Its source matrices have complete Laurent coefficient
arrays before entering the Smith routine.

At each ray, the profile products were compared at three logarithmic
scales, equal to one half, one, and twice the natural Gaussian
minimizer. No coefficient row improved the Gaussian product at the
same heights and scale. These partition comparisons are numerical;
the record makes no assertion about other scales. The smallest
sampled excess over the natural Gaussian pressure equation was
approximately 0.1873888884, and the smallest sampled excess over the
critical \(\sqrt2\) rate equation was approximately 0.2165186287.

The complete record is
[`reflection-hot-positive-rays-M16.json`](../evidence/reflection-hot-positive-rays-M16.json).
The compiled driver is
[`search_reflection_hot_positive_rays.py`](../tools/search_reflection_hot_positive_rays.py).

## Nonlinear two-bit state penalties

The two-bit screen uses, on every physical leg, the state heights

\[
h(00)=h(01)=h(10)=0,\qquad h(11)=t_\nu,
\]

where the zero and infinity values are chosen independently from
\(\{0,1,2,4,8,16\}\). This is a nonlinear function of the
two occupation bits. It satisfies the monotonicity condition for
the eight-mode reduction at every interval length.

The selected tensors were the ten best recorded classes from the
M16 single-defect screen. Every one of them attained its best value
on this 36-point grid at \(t_0=t_\infty=0\). Each fixed
height pair has exact active profiles; the logarithmic scale was
then optimized without repeating Smith elimination.

The driver
[`search_reflection_hot_joint_heights.py`](../tools/search_reflection_hot_joint_heights.py)
also replays four complete original weighted profiles and 131,072
literal endpoint coefficients for nonlinear two-bit heights before
starting the finite screen. An independent complete-alphabet
pressure calculation agrees with those controls.

The newer cold-compensated negative-height construction is a
separate family. Its [completed searches](reflection-hot-compensated-searches.md)
are not included in the counts or conclusions above.
