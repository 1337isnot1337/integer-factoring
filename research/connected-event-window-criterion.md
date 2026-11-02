# Connected event correlations give a finite-order window criterion

## Result and what is still missing

There is an explicit sufficient connected-correlation estimate for the
actual periodic [relation-change events](relation-change-event-density.md).
It gives a window linear in inverse event density, a correlation scale,
and the logarithm of a period bound. The elementary implication is proved
below. **The required arithmetic correlation estimate is not proved.**
Actual target-subgroup indices also remain open, so the complete factoring
score remains **1/2**.

The mechanism is the classical cumulant method: connected correlations
add under summation, and their factorial bound controls an even centered
moment. Compare Féray, [*Weighted dependency graphs*](https://arxiv.org/abs/1605.03836),
Sections 2 and 4, and Döring--Jansen--Schubert,
[*The method of cumulants for the normal approximation*](https://doi.org/10.1214/22-PS7).
Those sources do not establish a correlation theorem for these arithmetic
events. The explicit finite-order constants, target-cycle specialization,
and complete-cost interface below are the deductions recorded here.

The [registered scope](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19355534)
is Sol worker 1's analytic lane, based on `3ed3569` and refreshed to
`af12893`. No new orbit, cutoff search, or Lean formalization is involved.

## 1. Only finitely many ordinary cumulants are needed

Let `e:Z/h -> {0,1}` be a periodic flag with density `delta>0`. Choose
the start `x` uniformly over the **actual target cycle**, and put

```math
X_x=\sum_{j=0}^{W-1}e(x+j),\qquad \mu=\delta W.
```

For a random variable `X` on this finite probability space, define its
ordinary cumulants by the formal power-series identity

```math
\log\mathbb E\exp(tX)=\sum_{r\ge1}\kappa_r(X)t^r/r!.
```

Fix an integer `q>=1` and a real scale `A>=1`. Suppose **only** that

```math
|\kappa_r(X)|\le\mu\,r! A^{r-1}\quad(2\le r\le2q),
\qquad \mu\ge36Aq.                                      \tag{1}
```

Then

```math
\boxed{\Pr(X=0)<2^{-q}.}                                \tag{2}
```

Here is a finite-order proof, avoiding an assumption about convergence
of the full logarithmic series. Set `Z=X-mu`; its first cumulant is zero
and its other cumulants through order `2q` are unchanged. The
moment--cumulant identity and coefficientwise absolute values give

```math
\frac{\mathbb E Z^{2q}}{(2q)!}
 \le[t^{2q}]\exp\left(\frac{\mu A t^2}{1-At}\right).
                                                               \tag{3}
```

Only the first `2q` coefficients on the right enter the comparison.
This positive-coefficient majorant can be evaluated at any `0<t<1/A`.
Take `t=4q/mu`. Under (1), `At<=1/9` and

```math
\frac{\mu At^2}{1-At}\le q/2.
```

Consequently, even-moment Markov and `(2q)!<=(2q)^(2q)` give

```math
\Pr(X=0)\le\frac{\mathbb E Z^{2q}}{\mu^{2q}}
 \le4^{-q}\exp(q/2)<2^{-q},
```

as asserted. A high-order joint law or an independence hypothesis is not
silently inferred from a covariance bound.

If `h<=Hbar<2^q`, every nonzero empty-window probability is at least
`1/h>=1/Hbar>2^-q`. Thus (1) implies that **every** length-`W` window
hits. Given integer bounds `delta>=1/D0`, `A>=1`, choose

```math
q=\operatorname{bitlength}(\overline H),\qquad
W=36 A q D_0.                                           \tag{4}
```

It suffices to prove the cumulant bounds in (1) at this particular `W`.
If `W>=h`, a hit is automatic and no correlation hypothesis is needed.
This convention avoids imposing incompatible independence assumptions
on repeated complete cycles. The algorithm computes neither cumulants
nor the unknown period.

## 2. A concrete sufficient connected-correlation mechanism

For time indices `j_1,...,j_r`, define the joint ordinary cumulant by

```math
\kappa(e_{j_1},\ldots,e_{j_r})=
 \sum_{\pi\text{ a partition of }[r]}
 (-1)^{|\pi|-1}(|\pi|-1)!
 \prod_{B\in\pi}\mathbb E\prod_{i\in B}e(x+j_i).
                                                               \tag{5}
```

Repeated time indices are allowed. Multilinearity gives the exact identity

```math
\kappa_r(X)=\sum_{0\le j_1,\ldots,j_r<W}
                 \kappa(e_{j_1},\ldots,e_{j_r}).          \tag{6}
```

One strong sufficient arithmetic mechanism would be a nonnegative even
kernel `K(s)`, for `|s|<W`, satisfying, for all these tuples and
`2<=r<=2q`,

```math
|\kappa(e_{j_1},\ldots,e_{j_r})|
 \le\delta\sum_{T\text{ a labeled tree on }[r]}
           \prod_{\{a,b\}\in E(T)}K(j_a-j_b),
\qquad M:=\sum_{|s|<W}K(s).                              \tag{7}
```

Root a tree at its first vertex and eliminate leaves. Its sum over all
time tuples is at most `W M^(r-1)`, including tuples with equal indices.
Cayley's tree count and `r^(r-2)<=r! 3^(r-1)` therefore give

```math
|\kappa_r(X)|\le\mu r^{r-2}M^{r-1}
                    \le\mu r!(3M)^{r-1}.                \tag{8}
```

For the last elementary inequality, integrate `log x` to obtain
`r!>=(r/e)^r`, then use `e<3` and `r^2>=3`. Thus
`A=max(1,ceil(3M))` suffices. The **weaker signed sum bound (1)** is
also enough; absolute tree domination need not hold if cancellation
controls (6) directly.

Condition (7) is deliberately explicit and strong. It has no additional
uncontrolled `r!` prefactor. A family of constants depending on `r`,
a fixed-order asymptotic normality theorem, a pair-correlation estimate,
or a small-bias theorem for a separate prime walk does not imply it at
orders through `2 bitlength(Hbar)`. Periodic aliases must be retained.
This note neither posits nor proves a Markov model for the power orbit.

## 3. Conditional complete factoring invoice

For an input of `b` bits, the existing elementary bound `Hbar=2N` gives
`q=O(b)`. Suppose effective uniform theorems provide integer functions
`D0(b),A(b)` and establish (1) at (4), for every prime direction used
by the smooth-power selector, whenever that window is shorter than its
period. Apply this to the **row-change flag**, whose successful events
are accepted by the original recognizable predicate.

The selector can simply examine its original powers starting at `K_p`;
it need not compute the event flag to stop at any admissible hit. Its
existing diagonal-dominance and ambiguous-class extraction contracts
then apply. Prime/base enumeration, word lengths, all rejected queries,
and exact extraction are charged as in the
[moment transfer](cofactor-orbit-moment-criterion.md#4-what-uniform-estimates-would-imply-for-the-score).
In its more general cofactor-closure version the same assumptions are
needed on every reached residual direction, with at most
`B bitlength(B)` keys. The complete cost is polynomial in
`b,Y,B,W`, including elimination and decoding. Hence

```math
\log T_{\rm factor}
 =O(\log(b+1)+\log(Y+1)+\log(B+1)+\log D_0+\log A).
                                                               \tag{9}
```

For a prime orbit covered by the [small-index density theorem](small-index-smooth-density.md),
with index `<=D^kappa`, fixed `kappa<1/4`, take
`Y=ceil(exp((log b)^(alpha+1)))`. That theorem and the row-change bridge
give `log D0=O_alpha(b/(log b)^alpha+(log b)^(alpha+1))`.
If **all needed indices are covered** and **(1) holds uniformly with
`log A=O(b/(log b)^alpha)`**, affordable `B` would make (9) the desired
saving exponent `alpha`. For example polynomial `A` would suffice.
These are conditional implications, not a new bound. Ramified or
otherwise uncovered directions require their own handling; the density
theorem about unramified prime products is not silently applied to them.

This does not improve the logical strength of the previous moment
criterion for free: it replaces its arithmetic hypothesis with a
specific connected-correlation hypothesis, now through `O(b)` orders.
The [ambient minorant obstruction](uniform-ambient-minorant-obstruction.md)
does not prove or disprove (1).

## 4. Exact controls and their limitations

```sh
make verify-event-cumulant-windows
```

The [checker](../tools/check_event_cumulant_windows.py) uses rational
moment--cumulant recurrences and the partition formula independently.
It checks mixed cumulant summation through order five with repeated
indices; deleting the diagonal already fails at order two. Three
independent Bernoulli sums satisfy the finite-order even-moment bound.
A regular cyclic flag meets the every-window criterion. A flag with one
long contiguous cluster fails it (`A>=14` through order 18).
Pairwise independent bits with nonzero third cumulant `-1/8` guard
against treating pairwise control as a tree theorem.

The [receipt](../evidence/event-cumulant-windows.json) also reads the two
saved generator-5/7 flags at period 198850. Their equality to the
row-change flags is pinned to the existing event certificate, which is
a prerequisite. This uses the **actual** generator-seven stride, not
the earlier shifted-generator-five control. At the existing window 3260
and the fixed quarter-period window 49712, aggregate cumulants through
order 12 need integer scales between one and two. These are finite
checks of aggregate cumulant algebra, not pointwise tree bounds or estimates through all the
required orders. There are only 61/73 events in these cycles, less than
`36 bitlength(h)=648`; consequently this sufficient criterion cannot
certify any nonwrapping window there even with `A=1`. No cutoff is enlarged
to evade that honest limitation.

Credit: PR82's event definition/density bridge, PR87's original subgroup
density, the repository's cofactor moment/complete-cost transfer, and
the saved selector and Grok entry-cycle data. Cumulant and tree methods
are classical as cited above. Developed with substantial OpenAI Codex
assistance; no novelty priority or independent mathematical review claimed.

The next arithmetic question is to bound the signed connected sum (6)
for the actual reduced-norm predicate at (4), or to identify a genuine
arithmetic obstruction. A generic low-bias sample does not settle it.
