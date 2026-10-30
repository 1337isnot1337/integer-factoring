# Positive cut duals exclude seeds through fifty-two labels

## Result

Every exact crossing decomposition over every field satisfies

```math
\sum_i\sqrt{r_i c_i}>d^{\sqrt2}\qquad(2\le d\le52).     \tag{1}
```

Thus a seed passing the strict sufficient finite rank-sum criterion
needs at least **fifty-three labels**. The
[residual-budget theorem](weighted-cut-residual-budget.md) already
proves this through fifty labels. This note closes dimensions fifty-one
and fifty-two using six explicit positive rational dual certificates.
It makes no claim about realizability of a fifty-three-label seed or
about a factoring time exponent above one half.

A numerical linear program was used only to discover the six duals.
The coefficients were rounded downward and rescaled using exact
arithmetic until every required inequality held. The proof checker
uses only integers and fractions. It neither trusts nor invokes the
numerical optimizer.

## 1. An arbitrary positive cut dual

Fix `k>=2`. Choose nonnegative rational numbers

```math
u,v,\quad a_s,b_s\quad(1\le s<k).
```

Put

```math
Q(r,c)=ur+vc+
\sum_{s=1}^{k-1}\bigl(a_s\min(r,s^2c)+b_s\min(s^2r,c)\bigr).
                                                               \tag{2}
```

Suppose the following `2k-1` rational tests hold:

```math
Q(j^2,1)\le j\quad(1\le j\le k),\qquad
Q(1,j^2)\le j\quad(2\le j\le k).                     \tag{3}
```

Then

```math
Q(r,c)\le\sqrt{rc}
\quad\text{whenever }r,c>0\text{ and }k^{-2}\le r/c\le k^2.
                                                               \tag{4}
```

Indeed, divide by `c` and set `t=r/c`. The function `Q(t,1)`
is affine between consecutive knots in

```math
\{1/k^2,1/(k-1)^2,\ldots,1/4,1,4,\ldots,k^2\}.
```

Equations (3), with homogeneity at the reciprocal knots, bound its
value there by `sqrt(t)`. On each intervening interval, its affine
value is bounded by the chord between the two square-root values.
Concavity of the square root puts that chord below `sqrt(t)`.
This proves (4). No assertion that the selected dual is optimal is
needed.

## 2. Apply the dual after removing an actual extremal piece

Use the same classification as the two preceding weighted-cut notes.
The case `k=1` has all cut ranks equal and gives `S>=d^2`.
Otherwise an actual piece has ranks `(a,b)`, after globally exchanging
cuts if necessary, with

```math
a_0:=(k-1)^2b+1\le a\le a_1:=\min(k^2b,d^2),\qquad
1\le b\le\left\lfloor\frac{d^2-1}{(k-1)^2}\right\rfloor.
                                                               \tag{5}
```

All remaining pieces have rank ratio between `1/k^2` and `k^2`.
Their ordinary rank totals are at least `R_0,C_0`, as defined and
proved in equations (5)–(6) of the residual-budget note. In particular,
that definition retains the all-field swap/product refinement when
`b=1`; those classes are not silently discarded.

Subtract the selected piece's contribution from the universal
weighted-cut inequalities. The remaining pieces satisfy

```math
\begin{aligned}
\sum_{i\ne A}\min(r_i,s^2c_i)
&\ge L_s:=\max\{0,sd-\min(a_1,s^2b)\},\\
\sum_{i\ne A}\min(s^2r_i,c_i)
&\ge M_s:=\max\{0,sd-b\}.                             \tag{6}
\end{aligned}
```

The first inequality uses the upper bound `a<=a_1`; subtracting
`min(a_0,s^2b)` would have the wrong direction. For the second,
`a>b` gives `min(s^2a,b)=b`.

Sum (4) over the remaining pieces and apply their rank budgets to
the nonnegative coefficients in (2). Adding the selected piece gives

```math
\boxed{\quad
S\ge\sqrt{a_0b}+uR_0+vC_0+
             \sum_{s=1}^{k-1}(a_s L_s+b_s M_s).
\quad}                                                        \tag{7}
```

This is valid for every positive rational dual satisfying (3).
The only square root in its lower bound is replaced in the checker by

```math
\frac{\operatorname{isqrt}(10^{24}a_0b)}{10^{12}}.
                                                               \tag{8}
```

Equation (8) is rounded downward. All subsequent calculations are
exact rational calculations.

## 3. Exhaust all scalar classes in the two new dimensions

The checker enumerates every pair `(k,b)` in (5) for `d=51,52`.
It first computes the existing residual bound `J(d,k,b)`. If that
already proves `S>d^(577/408)`, no extra dual is needed. Exactly six
classes remain:

| `d` | `k` | `b` | Lower bound from the stored dual, rounded down |
|---:|---:|---:|---:|
| 51 | 34 | 2 | 261.2131818666 |
| 51 | 35 | 2 | 262.3777868004 |
| 52 | 34 | 2 | 267.8057173309 |
| 52 | 35 | 2 | 267.7598560452 |
| 52 | 36 | 2 | 268.9529053983 |
| 52 | 37 | 2 | 270.2256121366 |

The complete rational coefficient lists are in
[`weighted-cut-dual-witnesses.json`](../evidence/weighted-cut-dual-witnesses.json).
Their order is
`u,v,a_1,...,a_(k-1),b_1,...,b_(k-1)`.

For each class, the verifier checks every coefficient's sign and
denominator, every knot inequality (3), and every arithmetic input
to (7). It then verifies for its resulting rational lower bound `L`

```math
L^{408}>d^{577},\qquad577^2-2\cdot408^2=1.             \tag{9}
```

Consequently `L>d^(577/408)>d^sqrt(2)`. It also checks that the
stored dual keys match the full set of classes left by the earlier
bound, with no omissions, duplicate keys, or extra untested claims.
The equal-rank class `k=1` is verified separately.

This proves (1) using the previous result through fifty and the new
certificates in both remaining dimensions. The certificate stops at
fifty-two; no claim about dimension fifty-three follows from this
calculation.

## 4. Replay and proof status

```sh
python tools/check_weighted_cut_dual_certificate.py \
  --output evidence/weighted-cut-dual-certificate.json
```

The verifier reads the six finite rational witnesses, reconstructs
all coefficient inequalities, and enumerates all scalar classes for
the two new dimensions. It uses the standard library and does not
require SciPy. Discovery with floating-point linear programming is
separate from this exact proof replay.

The all-field argument is (2)–(7), together with the previously proved
quotient and rank-one complement theorems. No new Lean proof or
improvement to the factoring-time exponent is asserted.
The proofs and exact certificate replay were developed with substantial
OpenAI assistant support.
