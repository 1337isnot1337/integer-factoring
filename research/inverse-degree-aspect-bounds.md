# Inverse-degree aspect bounds and a complete binary-weight exclusion

## Result and scope

A scalar polynomial crossing with a low-degree polynomial inverse cannot
gain arbitrarily much by increasing its own degree. For physical alphabet
sizes q=7,8,9,10, the following exact constants work:

| q | A_q | Rational comparison scale a_q |
| --- | --- | --- |
| 7 | 14 | 3/2 |
| 8 | 20 | 269/200 |
| 9 | 35 | 239/200 |
| 10 | 91 | 130/121 |

Let R(w) be a polynomial operator on Mat_q over any field, with R(1)=I
and polynomial inverse of degree d>=1. Put M=deg R and construct the
crossing by reversing Phi(R). Then the bare ungauged Smith pressure
cannot improve the half-rate threshold whenever

```math
M\ge A_q d.
```

If d<=2, the stronger degree-only exclusion M>=A_q holds. Both conclusions
give rho>17/24>1/sqrt(2). These are restrictions on the stated certificate,
not on arbitrary coefficient-rank evaluations or factoring algorithms.
Independent physical gauges and rectangular auxiliary inverse problems
are not covered. No stronger factoring score is established.

The proof combines [inverse-degree compatibility](inverse-degree-crossing-compatibility.md),
the [first-two-coefficient tradeoff](nine-letter-quadratic-inverse-probe.md),
and the [rational entropy dual](leading-rank-entropy-floor.md). It uses
substantial OpenAI Codex assistance, with no new Lean formalization.

## 1. A quantitative version of the leading-rank dual

If the minimum order in either cut has multiplicity at least two, the
entropy theorem already excludes improvement. Otherwise, let the infinity
profile begin 0,k,.... Write H_q for that theorem's exact entropy constant,
alpha=(1/2,1/5,3/10), and

```math
\eta_q=\alpha_0/2-\alpha_1/(q-2),\qquad\kappa_q=2\eta_q.
```

The leading two channels now contribute at least 2s^(-k/2), instead of
the constant 2 used in the repeated-minimum proof. The rest of the
q-entry prefix has total order decreased by k. Repeating its weighted
arithmetic-geometric mean calculation therefore gives

```math
\log B(s)\ge H_q-\eta_q k\log s.
```

All remaining degree terms are still a nonnegative multiple of the
shared-prefix slack. Thus this is a global inequality, without an
assumption that s is close to one.

Also B(s)>=s^M: after reversal, each cut has a zero-order channel.
Choose the scale split s=a_q^(2/k). On the lower interval,

```math
B(s)\ge e^{H_q}a_q^{-\kappa_q};
```

on the upper interval, B(s)>=a_q^(2M/k). Exact rational comparisons prove

```math
e^{H_q}a_q^{-\kappa_q}>q^{17/6},\qquad
a_q^{A_q}>q^{17/6}.
```

Hence 2M>=A_q k excludes improvement at every scale. The inverse-degree
inequality k<=2d then proves the first assertion M>=A_q d.

## 2. Why a quadratic inverse gives a stronger cap

If d<=2 and the minimum infinity order has multiplicity one, let c and e
count orders one and two. The first-two-coefficient theorem gives

```math
e\ge\max(0,q^2-c^2-7c-10).
```

If c>0, then k=1. If c=0, then e>=q^2-10>0 for q>=7, so k=2. Thus
k<=2, and 2M>=A_q k holds whenever M>=A_q. The repeated-minimum case
was already excluded. No extra leading-rank hypothesis is needed for
the final quadratic-inverse degree caps.

## 3. Exact arithmetic certificate

Let D be the common probability denominator in the entropy dual. Then
exp(D H_q) and a_q are rational, and D kappa_q is an integer. The
[checker](../tools/check_quadratic_inverse_degree_caps.py) verifies

```math
\left[e^{D H_q}/a_q^{D\kappa_q}\right]^6>q^{17D},\qquad
a_q^{6A_q}>q^{17}.
```

These are integer comparisons after clearing denominators. The
[record](../evidence/quadratic-inverse-degree-caps.json) contains the exact
scales, penalties, and degree caps. The split itself may be an irrational
root of a_q, but the proved comparisons require no numerical root or log.

## 4. All 2^31 binary kernel-level weights on the fixed nine-letter flag

Return to U(X)=J_9(1)X^T, N=I-U and the deterministic kernel flag of the
[nine-letter construction](nine-letter-quadratic-inverse-probe.md).
Give each of its 31 active input levels either weight one or weight two.
If t levels receive weight two, the actual degree is M=31+t. The
inverse degree is at most two, and the endpoint and rank-one row constant
are unchanged.

For t>=4 the degree cap M>=35 excludes every path immediately.

For t<=3, let F be the first-coefficient matrix after quotienting by the
constant row and column spaces when all weights are one. Its exact rank
is 80. Let F_j be the change made by giving level j weight two. The 31
matrices F_j are already the first-jet directions computed by the
constructor. Their three largest ranks are 22,22,18. Therefore, for
any set J of at most three changed levels,

```math
\operatorname{rank}\left(F-\sum_{j\in J}F_j\right)
\ge80-22-22-18=18.
```

This rank is exactly the number c of order-one infinity channels. For
1<=s<=5/4, the inverse-degree bound nu_i<=M+2 gives

```math
B(s)\ge81s^{-2}(1+c/s)\ge99792/125.
```

For s>=5/4, B(s)>=s^M>=(5/4)^31. The checker raises both rational
bounds to the seventh power and verifies that they exceed 9^20. Thus
this case has rho>5/7, while the high-degree case has rho>17/24.

Together these arguments exclude all 2,147,483,648 binary assignments
on this fixed flag, and all characteristic-two extension fields of the
displayed binary matrices. The proof uses 31 exact first-jet ranks and
the abstract degree cap, not an enumeration of the assignments. It does
not cover arbitrary changes of flag, arbitrary matrix entries in S,T,
or independent gauges.

The [constructor checker](../tools/check_quadratic_inverse_paths.py)
with `--alphabet 9` reconstructs these directions and the rank-80
reference matrix, then verifies both cases. The
[saved certificate](../evidence/quadratic-inverse-nine-letter-probe.json)
separates the 31 directly checked directions from the 2^31 assignments
covered by the proof.

## 5. Constructive implication

For nine letters, a potentially useful quadratic-inverse path must have
M<=34, and it must use more freedom than these binary input-level weights.
Changing many level weights only moves it into the degree obstruction;
changing at most three leaves too many early channels. A general
noncommuting pencil can change several entries within and between levels
while preserving the endpoint. Its full polynomial inverse, early ranks,
and complete profiles still have to be controlled together.

Reproduce the exact proof checks with

```sh
python3 tools/check_quadratic_inverse_degree_caps.py
python3 tools/check_quadratic_inverse_paths.py --alphabet 9
```
