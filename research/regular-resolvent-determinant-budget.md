# A determinant budget for every regular nilpotent resolvent

## Result and scope

Let a regular nilpotent operator on `Mat_q(F)` define the
[resolvent crossing](nilpotent-flag-crossings.md), over any field `F`.
Its existing bare two-place Smith-profile criterion cannot reach the
half-score threshold for

\[
 2\le q\le18,\qquad\text{or}\qquad q=21.                 \tag{1}
\]

This conclusion covers every leading row rank and every physical basis.
It also covers endpoint-preserving fractional-linear moves of any finite
nonzero, nonendpoint row place to infinity, including field extensions.
It does not assume the normalized two-transpose-channel form, binary
coefficients, a permutation source, or a particular Jordan chain.
In addition, every regular source of leading row rank at most two fails
the same criterion for **every** alphabet size `q>=2`.

The new ingredient is the total determinant-degree budget. The known
zero and order-one pivots constrain the remaining orders through that
budget, so their contribution cannot simply be discarded. Exact rational
weighted arithmetic-geometric mean certificates prove the finite-size
claims uniformly in the scale. This is a restriction on a sufficient
profile estimate, not a lower bound for coefficient ranks or factoring.
Physical types, gauges, other generator transforms, and different
coefficient decompositions remain outside this conclusion.

## 1. Keep the remaining Smith orders

Write `n=q^2`, `D=n-1`, `U=I-S`, and

\[
 H(w)=\Phi\bigl(U(I-wS)^{-1}\bigr).
\]

Here `S` is regular nilpotent, so `H` is a polynomial matrix of degree
at most `D`. Its value at one is the row crossing matrix, which is
invertible. Thus its `n` Smith orders `mu_i` at zero are finite,
nonnegative integers and

\[
 \sum_i\mu_i=v_w(\det H)\le\deg\det H\le nD.              \tag{2}
\]

Let `r=rank Phi(U)`. Exactly `r` orders are zero. The
[partial-transpose power-rank theorem](partial-transpose-power-rank.md)
forces at least

\[
 k=\max\left(0,\left\lceil\frac n{r^2}\right\rceil-2r\right) \tag{3}
\]

orders to equal one. In particular `k<=n-r`. The column orders of the
regular resolvent are `(0,n,...,n)`. Therefore its bare profile cost is

\[
 B(s)=s^D\left(\sum_i s^{-\mu_i}\right)
                 (1+(n-1)s^{-n}),\qquad s\ge1.             \tag{4}
\]

Select `k` of the order-one pivots and set `ell=n-r-k`. When `ell>0`,
their remaining sum is at most `nD-k`. The convex decreasing function
`x -> s^(-x)` and Jensen's inequality give

\[
 \boxed{\quad
 B(s)\ge
 \left[r s^D+k s^{D-1}
       +\ell s^{-(rD+k(D-1))/\ell}\right]
       (1+(n-1)s^{-n}).\quad}                             \tag{5}
\]

The statement is also valid at `s=1` by direct substitution. It does
not assume the remaining orders are equal or that a profile attaining
the right side is realizable. If `ell=0`, the first two terms already
give the full row sum.

The earlier scalar bound remains useful:

\[
 B(s)\ge nr+(n-1)k.                                       \tag{6}
\]

It settles most small `r`. Equation (5) settles the previously missing
intermediate ranks.

## 2. One zero-only certificate covers all larger ranks

Let `r_0` be the least positive integer with `2r_0^3>=n`.
Then (3) has `k=0` for every `r>=r_0`. The zero-only row relaxation is

\[
 Z_r(t)=r+(n-r)t^{nD/(n-r)},\qquad t=s^{-1}\in(0,1].        \tag{7}
\]

For fixed `t`, this is nondecreasing in `r`. Indeed put `x=n-r`
and `C=-nD log t>=0`. The derivative of
`n-x+x exp(-C/x)` with respect to `x` is

\[
 -1+e^{-C/x}(1+C/x)\le0,
\]

using `e^u>=1+u`. The case `r=n` follows by its limiting value,
or directly from the full constant rank. Thus a certificate for
`r=r_0,k=0` covers every `r>=r_0`. Only `1<=r<=r_0` needs to be listed.

## 3. Small integer weights remove the scale exactly

The certificates for (5) use positive integers `a,b,c`. If `k=0`, the
value of `c` is immaterial. Put

\[
 R_0=rb,\quad R_1=kc,\quad R_2=\ell a,\quad
 Q=R_0+R_1+R_2,
\]

and choose the integers so that

\[
 V=\frac{rD(b-a)+k(D-1)(c-a)}{n}\in\{1,\ldots,Q-1\}.       \tag{8}
\]

The three row weights are `R_i/Q`, omitting a zero `R_1` when `k=0`.
The two column weights are `(Q-V)/Q,V/Q`. Equation (8) says that the
weighted row exponent is `nV/Q`; the weighted column exponent is
`-nV/Q`. Applying weighted arithmetic-geometric mean separately to the
two factors of (5) consequently cancels the entire scale and gives

\[
 B(s)^Q\ge
 \frac{Q^{2Q}(n-1)^V}
 {b^{R_0}c^{R_1}a^{R_2}(Q-V)^{Q-V}V^V}.                  \tag{9}
\]

Every object in (8)--(9) is an integer. To certify `B(s)>K`, the
checker tests the single integer inequality obtained by multiplying
the denominator in (9) by `K^Q`. No logarithm, numerical minimizer,
or interval sampling enters this verification.

For each size in (1), the retained certificate covers every
`1<=r<=r_0` with either (6) or (9). Its integer lower bound `K`
satisfies

\[
 K^{408}>n^{577},\qquad577^2>2\cdot408^2.                 \tag{10}
\]

Thus `B(s)>n^sqrt(2)` for every `s>=1`, as required. There are
`76` rank cases in the retained table; the exact checker reports the
split between early-pivot and weighted certificates. The largest
weight denominator is only `2312`.

The weighted witnesses can be discovered by a numerical search, but
the committed witnesses are integers and the default replay does not
run that search. Their positive margins and (10) are verified exactly.

## 4. Rank two is excluded at every size

For `q>=5`, the existing partial-transpose theorem already proves
failure for leading row rank at most two. At `q=2`, (6) suffices.
The two remaining sizes have particularly simple certificates:

| `q` | `n` | `r,k` | `a,b` | `Q,V` | Proved bound |
|---:|---:|---:|---:|---:|---:|
| 3 | 9 | 2,0 | 6,51 | 144,80 | `B>27` |
| 4 | 16 | 2,0 | 9,121 | 368,210 | `B>64` |

These follow from (9) by exact integer comparisons. In both rows
the bound is `q^3>q^(2sqrt(2))`. Rank one is already covered by (6).
This closes the regular rank-two profile question over every field,
including sources that require scalar extension for a normal form.
It is stronger in scope than the separate finite binary normalized
four-label census, while that census gives much stronger actual pivot
counts in its smaller parameter space.

## 5. The relaxation itself does not exclude size 19

The first gap after the interval in (1) is real for these particular
constraints. At `q=19`, `n=361`, consider the abstract integer profile

\[
 (0^6,366^{325},367^{30}).                                \tag{11}
\]

It has `n` entries, total `n(n-1)=129960`, and satisfies (3) with
`r=6,k=0`. At the rational scale `s=73/72`, its cost from (4)
satisfies the exact comparison

\[
 B(s)^{169}<361^{239},\qquad239^2<2\cdot169^2.             \tag{12}
\]

The checker evaluates (12) as an integer numerator-denominator
comparison. Thus the determinant and first-pivot constraints alone
permit a value below the half-score threshold at this size. No actual
operator or crossing realizing (11) is supplied. Additional power
ranks, shared minors, physical restrictions, and construction costs
are independent obligations. The isolated exclusion at `q=21` comes
from the integer jumps in (3), so monotonicity in alphabet size should
not be inferred.

## 6. The reverse bare orientation also fails

At infinity the column leading coefficient is `U`, which is invertible,
so all column infinity orders are zero. The row orders at zero have
total at most `nD`, by the same determinant-degree argument. Jensen
therefore bounds the reverse product by

\[
 n s^D\sum_i s^{-\nu_i}\ge n^2.
\]

It cannot improve the threshold.

## 7. Fractional-linear moves remain in the same family

Fix `gamma!=0,1`, in the coefficient field or an extension, and set

\[
 \phi(x)=\frac{\gamma x}{x+\gamma-1},\qquad
 S'=(\gamma-1)S(\gamma I-S)^{-1}.                         \tag{13}
\]

The parameter move preserves zero and one and sends infinity to
`gamma`. Since
`S'=((gamma-1)/gamma)S+O(S^2)`, the operator `S'` is again regular
nilpotent. Its nilpotence index is `n` over every extension field.

Write the original operator polynomial as
`K_S(z)=z^n(zI-S)^(-1)(I-S)` before applying `Phi`. Direct multiplication
gives

\[
 \phi(x)I-S
  =\frac{(\gamma I-S)(xI-S')}{x+\gamma-1},
 \qquad I-S'=\gamma(I-S)(\gamma I-S)^{-1}.
\]

Consequently the cleared, endpoint-normalized transformed family is

\[
 \boxed{\quad
 \gamma^{-D}(x+\gamma-1)^D K_S(\phi(x))
    =x^n(xI-S')^{-1}(I-S')=K_{S'}(x).\quad}              \tag{14}
\]

Apply the linear reshuffling `Phi` to obtain the identical tensor
statement. Thus every such finite-place move is another regular
resolvent, and the all-rank exclusions (1) apply without change. The
point `gamma=1` has full row rank and cannot improve the bound. The
only finite column determinant root is zero, so choosing any other
column place gives the full-rank side of the elementary reverse
bound instead.

This closure does not show that rank two is preserved by the move.
For `q>=19` outside (1), the all-size rank-two conclusion concerns
the leading rank of the transformed source itself. Other generator
polynomials and other profile methods remain open directions.

## 8. Reproduction

```sh
python3 tools/check_regular_resolvent_determinant_budget.py
python3 tools/check_regular_resolvent_determinant_budget.py \
  --generate --output /tmp/regular-resolvent-determinant-budget.json
```

The first command replays every retained exact certificate, the
size-19 relaxation witness, and 27 cleared polynomial identities over
three prime fields. The second additionally searches for
integer weights before applying the same exact verification.
The [receipt](../evidence/regular-resolvent-determinant-budget.json)
contains the complete rank table, integer weights, integer hashes,
and the exact rational cost in (12). No factoring score change or
new Lean formalization is claimed.
