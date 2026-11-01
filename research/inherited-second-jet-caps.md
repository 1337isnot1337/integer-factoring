# Nonlinear second-jet caps survive independent gradings and restriction

The existing universal second-jet adjugate theorem also strengthens the
independently graded physical-subalphabet constraints. Its nonlinear
profile information must be checked on the **two actual profiles**.
Substituting their average into the minimum function is unjustified.

An exact twelve-letter bare point passes the preceding linear constraints
and has critical pressure below `-1/100`, but fails these already necessary
nonlinear conditions. It is not a surviving tensor target. The
[uniform bare obstruction through eleven letters](multichannel-bare-eleven-letter-floor.md)
remains the established coverage; twelve letters and independently gauged
constructions remain open. The factoring score stays **1/2**.

## 1. Input and graded consequence

The input is equation (19) in the credited
[universal second-jet theorem](https://github.com/1337isnot1337/integer-factoring/blob/365f319fbb62c734bb409c93d39323d566b87c49/research/adjugate-seven-letter-pressure-floor.md#8-a-stronger-universal-second-jet-constraint).
For a regular square physical tensor with normalized crossed source
profile `s`, its crossed adjugate orders satisfy

```
eta_i >= min(s_i,2*s2).
```

This holds over every field, including singular leading channel factors.
No new adjugation identity is claimed here. Put
`F_k(s)=sum_(i<=k) min(s_i,2*s2)`. Substituting this prefix bound into the
existing graded cell and shared-line arguments gives, on `m` physical
letters with `n=m^2` and `C=(n-1)d+m*(sum X+sum Y)`,

```
k S_(n-1)(r) + F_k(s) <= k C - TopCell_k,       1<=k<=n, (1)
```

```
(l-1) S_n(r) + S_(n-l)(r) + F_l(s)
    <= l C - max(TopX_l+l*Ymax,TopY_l+l*Xmax), 1<=l<=m. (2)
```

Both cut orientations hold. The cofactor grading and negative energy lists
are those of the [graded adjugate input](https://github.com/1337isnot1337/integer-factoring/blob/365f319fbb62c734bb409c93d39323d566b87c49/research/graded-adjugate-caps.md).
Only the substituted lower profile differs from its rank-one version.

## 2. Restore restricted minima before lifting

Restrict to `m` smallest paired-energy letters on each side. Let its raw
profiles relative to the full normalization be `u,v`, with minima `a,b>=0`.
Its own normalized profiles are `u-a,v-b`, and its deficit is `d-a-b`.
For a sorted raw list beginning at `b`, the elementary inequality

```
F_k(v-b) >= F_k(v) - (2k-1)*b                          (3)
```

holds: the first summand changes by exactly `b`, and each later summand
changes by at most `2b`. The first-profile coefficients in both (1) and
(2) sum to `k(n-1)` with the corresponding row size `k`. Restoring the
local normalization therefore leaves an additional term

```
[k(n-3)+1]*b >= 0
```

on the raw left side. Dropping it proves the same inequality on the raw
restricted profiles. Individual Smith interlacing gives `u_i>=r_i` and
`v_i>=s_i`; all first-profile coefficients and every argument of `F_k`
are nondecreasing. Thus (1)--(2) transfer to the first `n` orders of the
full tensor. The argument also applies in the reverse orientation.

The right sides retain the old coordinatewise monotonicity in selected
energies. Choosing the smallest energy subsets suffices. Integer heights
give the polynomial proof and rational heights follow by ramification.
The [implementation](../tools/inherited_second_jet_caps.py) reuses the
source-pinned square right sides and replaces only the old rank-one prefix.

## 3. Averaging requires care

The function `F_k` is concave. Averaging the two oriented inequalities gives
the term `(F_k(r)+F_k(s))/2`, which can be **smaller** than
`F_k((r+s)/2)`. Consequently the latter substitution can reject a pair
without a valid argument. The checker accepts complete actual `r,s` and
evaluates each orientation separately. A balanced search is a restricted
search once these nonlinear conditions are imposed.

## 4. Exact controls and the twelve-letter point

The [retained point](../evidence/inherited-second-jet-q12-point.json) has
denominator `100000`, zero paired energies and 144 orders at both cuts.
It is obtained by rounding a numerical linear-relaxation point toward
smaller orders and a larger deficit. The
[checker](../tools/check_inherited_second_jet_caps.py) independently verifies
all preceding cell/shared, square adjugate, rectangular first-jet and
multichannel linear constraints, then certifies pressure below `-1/100`
with 100-bit outward arithmetic. The actual equal pair fails (1) and (2),
including their whole-alphabet determinant and longest-line cases.
No nonlinear feasible point or realizing tensor is supplied.

Nine independently graded literal strict tensors give 36 complete profiles
and exact nonlinear checks. Rational local-shift controls verify (3) and
the favorable restoration sign. The
[receipt](../evidence/inherited-second-jet-controls.json) retains counts,
full-alphabet rejections, source hashes and the outward pressure interval.
Run `make verify-inherited-second-jet`.

This extends and packages the preceding second-jet theorem; it does not
replace its proof. It postdates the full default source `7e1c1d8` and has
a separate focused replay. No new Lean work, independent review or improved
factoring algorithm is claimed. Substantial OpenAI assistant support was used.

## 5. A successor point survives the declared nonlinear checks

A balanced cutoff-four search gives a different
[exact twelve-letter point](../evidence/q12-second-jet-target.json).
Its integer deficit is `4691`, its paired energies are all zero, its two
144-order profiles are equal, and its thermal scale is `1/1000`. It passes
every original cell/shared cap, all 1,452 oriented square adjugate caps,
6,578 rectangular first-jet caps, 6,578 multichannel caps and all 1,452
inherited second-jet caps. The equal pair makes the other rectangular cut
orientation identical.

100-bit outward arithmetic certifies pressure **below `-1/100`**, numerically
about `-0.0101040`. This does not use the invalid average substitution:
the two complete proposed profiles are checked separately. The
[checker](../tools/check_q12_second_jet_target.py) and
[receipt](../evidence/q12-second-jet-target-controls.json) pin the point,
every constraint compiler and the exact pressure interval.

This is a point of the declared necessary relaxation, not a polynomial
crossing or an optimality theorem. The discovery searched one balanced
sector; it says nothing about the optimal asymmetric profile pair. In the
bare integer convention a proposed realization would have degree `4691`,
so feasibility and a useful construction remain substantial open tasks.
Run `make verify-q12-second-jet-target`. The factoring score remains `1/2`.
