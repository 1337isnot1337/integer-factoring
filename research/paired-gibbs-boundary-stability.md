# Positive Gibbs roots at extreme gauges

## Result and scope

The one-hour eight-symbol closed-auxiliary search recorded eight
`Nonpositive paired Gibbs root needs exact investigation` rejections
and one small-entropy rejection among 1,310,782 native profile queries.
That historical run did not save the rejected parameter vectors.
Those nine events cannot be individually reconstructed or classified
from the retained record.

This checkpoint constructs twelve separate comparable gauges inside
the same normalized bounds. Eight reproduce the nonpositive-root
guard even though exact interval certificates put their true Gibbs
ratios above one. In these examples, floating arithmetic rounds a
strictly positive partition remainder to zero. The guard prevents
those lost remainders from becoming false candidates.

A new, separate evaluator preserves the small remainder and the
physical entropy logarithmically. It checks the exact affine
consistency of the native selected minors, but does not claim that
these checks prove global profile attainment. Existing search
drivers, historical records, and source hashes remain unchanged.

## 1. The exact zero-rate remainder is positive

Consider a strict scalar polynomial crossing with complete local
orders `mu,nu` at finite real gauges and positive scale. Let
`G_X,G_Y` be the total physical energies, including the scale.
Define

\[
L=\log\sum_i e^{-\mu_i}+\log\sum_i e^{-\nu_i},
\qquad
C=L+\min G_X+\min G_Y.
\]

For every physical pair `(x,y)`, the polynomial

\[
f(z)=T_{xyyx}(z)
\]

is nonzero because `f(1)=1`. Each local minimum order is at most the
weighted valuation of this matrix entry. Adding the two places gives

\[
\min\mu+\min\nu
\leq G_X(x)+G_Y(y)
+t\bigl(\operatorname{ord}_0f-\deg f\bigr)
\leq G_X(x)+G_Y(y).
\]

Choose minimizing physical letters. Then

\[
\delta=\min G_X+\min G_Y-\min\mu-\min\nu\geq0. \tag{1}
\]

After selecting one minimum from each complete local profile, write

\[
\tau_\infty=\sum_{i\ne i_0}e^{-(\mu_i-\min\mu)},
\qquad
\tau_0=\sum_{j\ne j_0}e^{-(\nu_j-\min\nu)}.
\]

Both tails are strictly positive for the nontrivial alphabets here.
Therefore

\[
\boxed{C=\delta+\log(1+\tau_\infty)+\log(1+\tau_0)>0.} \tag{2}
\]

This is a statement about the exact complete profiles of a strict
polynomial crossing. It does not apply without modification to
arbitrary rational paths or an unverified numerical profile.

The historical evaluator computes large log-partition values and
physical minima before combining them. A tail smaller than the
relative precision of `1+tail` may disappear. Large positive and
negative minima can also cancel after the tail information has been
lost. Thus a floating value of zero for `C` does not contradict (2).
The old guard is evaluated at a small positive trial rate, so it is
not itself an exact test of the zero-rate limit.

## 2. Separate extreme-gauge witnesses

The [diagnostic](../tools/diagnose_paired_gibbs_boundary.py) uses the
closed auxiliary tensor of block length three with both auxiliary
masks zero. It checks all 4,096 coefficients by literal auxiliary-path
enumeration and checks every strict endpoint entry.

For each height `h` in `{3,6,12}`, use scale eight and assign the same
height vector at every physical leg and place. The first pattern is
`(0,h,h,h,h,h,h,h)`; the second is `(0,-h,0,0,0,0,0,0)`. Both local
cut orientations are tested. All normalized heights lie within the
historical bounds `[-12,12]`, and the scale equals the old upper bound.

For all twelve cases, the checker independently constructs the
complete weighted 64-by-64 polynomial matrices and certifies their
24 Smith profiles with exact regular frames. A 200-digit decimal
root calculation is diagnostic only. The pressure signs that certify
the retained root brackets use outward dyadic integer arithmetic with
denominator `2^768` and explicit exponential/logarithm remainders.

The two sign patterns give the same ratio brackets. The table lists
one copy of each height and cut.

| Normalized height | Cut | Floating guard | Certified ratio interval |
|---:|---:|---|---|
| 3 | 0 | Returns a value | `(1.721991, 1.721994)` |
| 3 | 1 | Returns a value | `(1.764327, 1.764330)` |
| 6 | 0 | Nonpositive root | `(1.270336, 1.270339)` |
| 6 | 1 | Nonpositive root | `(1.281940, 1.281943)` |
| 12 | 0 | Nonpositive root | `(1.119895, 1.119898)` |
| 12 | 1 | Nonpositive root | `(1.124414, 1.124417)` |

At height six the exact positive `C` is of order `1e-32`; at height
twelve it is of order `1e-73`. The old calculation returns zero in
all eight corresponding sign/cut cases. These are numerical boundary
effects, not low-ratio constructions.

The [saved diagnostic](../evidence/paired-gibbs-boundary-diagnostics.json)
retains the original search's rejection counts and explicitly sets
`historical_failures_individually_classified` to false. The equal count
of eight reproduced guard events does not identify them with the
eight historical events.

## 3. Exact affine checks on selected profiles

The [new evaluator](../tools/stable_paired_gibbs.py) obtains the native
attaining-prefix data. Each selected order is an integer affine form
in the physical heights and the common scale. The scale coefficient
is the difference of consecutive base-prefix valuations; the other
coefficients record the newly selected physical row and column labels.

Every finite binary floating input is an exact dyadic rational.
A common denominator therefore lets the evaluator reconstruct every
selected order and physical energy with integer arithmetic. It checks
agreement with the reported floating orders, exact nondecreasing order
within each cut, and the nonnegative endpoint gap (1). A negative gap
or inconsistent order sequence raises an error. Nothing is clipped
to force a positive result.

These checks are necessary consistency conditions on the selected
prefixes. They do not prove that floating weight comparisons chose
globally attaining minors at every rank. A candidate used in a theorem
still requires independent full-matrix profile certification. The
diagnostic witnesses above have that separate certification.

## 4. Stable root and gradient formulas

Put `R=2*rho` and define the minimum-shifted physical sums

\[
S_X(R)=\log\sum_x e^{-(G_X(x)-\min G_X)/R},
\qquad
S_Y(R)=\log\sum_y e^{-(G_Y(y)-\min G_Y)/R}.
\]

The paired Gibbs equation is

\[
C=R\bigl(S_X(R)+S_Y(R)\bigr).
\]

All quantities are positive. The evaluator solves the equivalent
logarithmic equation

\[
\log C-\log R-\log\bigl(S_X(R)+S_Y(R)\bigr)=0.
\]

Equation (2) supplies `log C` from the exact gap and logarithmic
partition tails. A stable `log(log(1+exp(x)))` routine retains tails
far below ordinary probability precision. The root search uses
`log R`, with explicit range failures if a positive root cannot be
represented safely.

The output retains both ordinary type probabilities and their
logarithms. It also retains the logarithm of the paired entropy.
A probability may round to one or underflow to zero in the ordinary
display without destroying the small entropy needed by the evaluator.

The gradient uses the same selected affine forms. If `d` denotes their
parameter derivative and one minimum is fixed on each side, the
fixed-type cost derivative is

\[
\begin{aligned}
A={}&d\delta
+\sum_x p_X(x)\bigl(dG_X(x)-d\min G_X\bigr)
+\sum_y p_Y(y)\bigl(dG_Y(y)-d\min G_Y\bigr)\\
&-\sum_i p_\infty(i)\bigl(d\mu_i-d\min\mu\bigr)
-\sum_j p_0(j)\bigl(d\nu_j-d\min\nu\bigr).
\end{aligned}
\]

The envelope gradient is `A/(2H)`. Ground-state differences vanish
exactly. Each remaining probability-to-entropy ratio is formed from
its logarithm before multiplying the integer affine coefficients.
This avoids subtracting order-one uniform-gradient terms and then
dividing their rounding error by a tiny entropy.

If the root is representable but a genuine gradient component exceeds
the floating range, the detailed API retains the positive root and
reports an unavailable gradient. The optimizer adapter raises a
separate error in that case; it does not return an invented zero
gradient.

## 5. Interfaces and checked controls

The detailed API is

```python
result = stable_paired_result(oracle, parameters, cut=1)
```

It returns the ratio, gradient status, ordinary and logarithmic types,
log entropy, exact leading-gap numerator and denominator, reconstructed
orders, and consistency flags. The four-item compatibility API is

```python
rho, gradient, px, py = stable_paired_value_gradient(oracle, parameters, cut=1)
```

The separate `paired_value_gradient_with_fallback` first uses the
historical evaluator. It replays its root and entropy guards through
the stable path. It also replays successful calls whose paired entropy
is below `1e-7`, because accuracy can degrade before a guard fires.
Existing callers do not change until they explicitly import this
adapter. Nonzero denominator degrees are rejected by this polynomial
endpoint-bound implementation.

The [controls](../evidence/stable-paired-gibbs-controls.json) passed all
twelve independently certified boundary brackets. Eight guarded
diagnostic calls and four successful but small-entropy calls were
replayed through the adapter. The smallest tested paired entropy has
logarithm about `-163.66`, while all corresponding envelope gradients
remain finite.

There are also 48 ordinary root/gradient comparisons on square and
rectangular Gaussian seeds and a closed auxiliary seed. The maximum
root difference is about `5.3e-14`, and the maximum gradient difference
is about `2.7e-15`. Forty-seven smooth finite-difference checks passed;
one nonsmooth direction was explicitly skipped. Eight boundary
directional derivatives were compared with independent 200-digit
complete-profile root calculations, with maximum error below `2.6e-9`.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/diagnose_paired_gibbs_boundary.py --output /tmp/paired-gibbs-boundary-diagnostics.json
OPENBLAS_NUM_THREADS=1 PYTHONPATH=tools python3 tools/check_stable_paired_gibbs.py --output /tmp/stable-paired-gibbs-controls.json
```

The [fixed-scale interval helper](../tools/paired_gibbs_interval.py)
is separate from the numerical evaluator. It certifies pressure signs
only when supplied exact integer profiles and exact energies. It
does not turn future floating profile outputs into exact certificates.
