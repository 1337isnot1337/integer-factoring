# Balanced partial repairs fail at every profile scale

## Result and scope

The scale restriction in the preceding balanced-gauge obstructions can be
removed. Keeping the positive finite Gaussian tail, instead of discarding
it, treats the small scales that those arguments left open.

**Theorem.** There is a fixed threshold such that, for every power-of-two
interval length `M` above that threshold, each of the following families
has strictly positive critical pressure for every `ell=log s>0` and every
real physical-core gauge vector satisfying

```math
\theta_0+\theta_1=\theta_2+\theta_3=0.
\tag{1}
```

- The mixed four-by-three repair with `alpha=1`, `beta=omega`, its
  original quartic, and all 25 parameters of the full-core construction.
- The full four-by-four partial repair with `gamma=omega` or
  `gamma=omega^2`, its original quartic, and all 29 parameters.

Here `omega^2=omega+1` in `F4`. The gauges and the profile scale may vary
arbitrarily with `M`. There is no logarithmic lower bound on `M ell`,
no upper bound on `ell`, and no boundedness or convergence requirement
on the gauges.

This is a restriction on the stated pressure certificate for these
particular partial repairs. It does not exclude unbalanced aggregate
bulk charges, different repairs or perturbations, other coefficient
rank estimates, or an isolated witness below the unspecified finite
threshold. It does not improve the factoring exponent.

The algebraic inputs are the
[uniform rectangular entry certificates](rectangular-uniform-ray-obstruction.md), their
[moving-gauge consequences](rectangular-balanced-gauge-obstruction.md),
and the [full-alphabet certificates](full-alphabet-leading-obstruction.md),
all at predecessor commit
`18cec06deac47c8f142ea67a60e2278ac2234d67`.
Their exact finite-field and rational-dual calculations are reused.
The new step is the uniform analytic argument below.

## 1. Retain the Gaussian term in the existing comparison

Write `R=sqrt(2)`, `tau=M ell`, and let `F` be the critical pressure.
Let `F_(G,L)(ell)` be the balanced Gaussian interval pressure at even
length `L`. In the conventions of the paired-type theorem it is

```math
F_{G,L}(\ell)
=4\sum_{j=1}^{L}h((2j-1)\ell)
 -4R\sum_{j=1}^{L/2}h(R(2j-1)\ell),
\qquad h(x)=\log(1+e^{-x}).
\tag{2}
```

The existing entry and dual certificates give a nonnegative quantity
`T`. The predecessor's coarse comparison can be strengthened without
changing any estimate to

```math
F\ge F_{G,M}(\ell)+T-a\ell-b,
\qquad T\ge0,
\tag{3}
```

where `a,b` are fixed nonnegative constants. The earlier proof replaced
the first term on the right by zero. Keeping it gives (3): its spectator
estimate is `S>=F_(G,M)-20 log 2-4 ell` before that replacement.

For the rectangular family one may take

```math
a=2B+8,\qquad b=20\log2+R\log12,
\tag{4}
```

where `B` is its common eventual entry-offset constant. For the full
alphabet replace `log 12` by `log 16`, and use that family's offset
constant. The constants are independent of the gauges, `M`, and `ell`.

If a non-low state minimizes one of the offset-free physical core
energies, the same exact dual table gives

```math
T\ge2\tau.
\tag{5}
```

When both minima are low, the off-face duals and the endpoint offsets
give a lower bound on every actual non-low gap. Under `F<=0`, equation
(3), the known nonnegativity of `F_(G,M)`, and `T>=0` imply
`T<=a ell+b`. The resulting gap bounds are

```math
\begin{array}{ll}
2\tau-b-(2B+12)\ell,&\text{rectangular family},\\
2\tau-2b-(4B+20)\ell,&\text{full alphabet}.
\end{array}
\tag{6}
```

The full-alphabet `X`-top gap has off-face multiplier `-2`; this is
why its second line uses the larger constants. It is not silently
treated as a `-1` multiplier.

The all-low restriction is exactly the Gaussian interval of length
`M-2`, with its scalar factors and fixed physical energies cancelling
in pressure. Its Smith profile partitions are bounded above by the
full profile partitions. Thus, whenever a putative `F<=0` has forced
both offset-free minima to be low, (6) gives constants `C>0,D>=0`
such that

```math
F\ge F_{G,M-2}(\ell)-C e^{-R\tau+D\ell}.
\tag{7}
```

Explicit choices are `C=5R exp(b/R)`, `D=(2B+12)/R` for the
rectangle and `C=6R exp(2b/R)`, `D=(4B+20)/R` for the full alphabet.
These are exactly the predecessor's freezing estimates.

## 2. Bounded total scale has a linear Gaussian margin

Fix any constant `tau_0>0`. Since `h` is positive and decreasing,
a midpoint Riemann sum of step `d` differs from its integral divided
by `d` by at most `h(0)=log 2`. Applying this to both sums in (2)
gives, for `tau=M ell`,

```math
F_{G,M}(\ell)
\ge {2\over\ell}\int_{R\tau}^{2\tau}h(x)\,dx
      -4(1+R)\log2.
\tag{8}
```

In particular, throughout `0<tau<=tau_0`,

```math
F_{G,M}(\ell)
\ge 2(2-R)M h(2\tau_0)-4(1+R)\log2.
\tag{9}
```

The coefficient of `M` is a fixed positive number. Combining (3)
and (9), and using `ell<=tau_0/M`, proves `F>0` for all sufficiently
large `M`, uniformly over this entire scale range and all balanced
gauges. No freezing conclusion is needed in this case.

## 3. Large total scale forces the low states

Choose one fixed `tau_0` so large that

```math
\tau_0>2b+2,
\qquad
\exp\!\left({2-R\over2}\tau_0\right)>4.
\tag{10}
```

Increase the eventual threshold for `M` so that `M>2a` and `M>=8`.
Suppose `tau>=tau_0` and `F<=0`. Then

```math
T\le a\ell+b={a\over M}\tau+b<\tau<2\tau.
\tag{11}
```

Equation (5) forces both offset-free physical minima to be low.
The bound (7) therefore holds. It remains to show that its Gaussian
term dominates its error uniformly over every remaining `ell`.

## 4. The positive finite tail pays for freezing

Put `L=M-2`. Define

```math
\begin{aligned}
F_{G,\infty}(\ell)
 &=4\sum_{j\ge1}h((2j-1)\ell)
   -4R\sum_{j\ge1}h(R(2j-1)\ell),\\
P&=4R\sum_{j>L/2}h(R(2j-1)\ell),\\
N&=4\sum_{j>L}h((2j-1)\ell).
\end{aligned}
\tag{12}
```

Then `F_(G,L)=F_(G,infinity)+P-N`. The
[interval floor proof](schur-one-species-interval.md) proves
`F_(G,infinity)(ell)>0` for every `ell>0`.

For `0<u<=1`, use `u/(1+u)<=log(1+u)<=u` and sum geometric
series. These give

```math
\begin{aligned}
P&\ge
 {4R e^{-R(L+1)\ell}
  \over(1-e^{-2R\ell})(1+e^{-R(L+1)\ell})},\\
N&\le {4e^{-(2L+1)\ell}\over1-e^{-2\ell}},\\
{P\over N}&\ge
 {\exp(((2-R)L+1-R)\ell)
  \over1+e^{-R(L+1)\ell}}.
\end{aligned}
\tag{13}
```

For the last line also use
`1-exp(-R x)<=R(1-exp(-x))`, valid for `R>=1`.
For `M>=8`,

```math
((2-R)L+1-R)\ell
=(2-R)\tau+(R-3)\ell
\ge {2-R\over2}\tau.
\tag{14}
```

The denominator in (13) is at most two. Equations (10), (13), and
(14) therefore give `P>2N` whenever `tau>=tau_0`.

The first line of (13), together with `1-e^{-x}<=x` and
`1+e^{-R(L+1)ell}<=2`, also gives the useful scale-sensitive bound

```math
P\ge {e^{-R(M-1)\ell}\over\ell}.
\tag{15}
```

Consequently the error in (7), relative to `P/2`, satisfies

```math
{C e^{-R\tau+D\ell}\over P/2}
\le 2C\ell e^{(D-R)\ell}.
\tag{16}
```

Choose one fixed `delta` with `0<delta<=1` so small that
`2C delta exp(max(D-R,0) delta)<1/2`. For `0<ell<=delta`,
(16) is less than one half. Thus

```math
F\ge F_{G,\infty}+P-N-Ce^{-R\tau+D\ell}
 > F_{G,\infty}+P/4>0.
\tag{17}
```

This contradicts `F<=0`. The missing factor in the earlier scale
argument was `1/ell` in the positive tail (15). It is precisely what
allows the tail to dominate the physical-state correction even below
the logarithmic scale.

## 5. Scales bounded away from zero

The same interval analysis gives a constant `c_0>0` such that

```math
F_{G,\infty}(\ell)\ge c_0\min(\ell,e^{-\ell})
\quad(\ell>0).
\tag{18}
```

This follows from strict positivity and continuity on compact
intervals, the expansion `ell/6+O(exp(-k/ell))` near zero, and the
leading term `4 exp(-ell)(1+o(1))` at infinity.

For `ell>=delta`, the minimum in (18) is at least
`delta exp(-ell)`. Since `P>N`, the ratio of the error in (7) to
this lower bound is at most

```math
{C\over c_0\delta}
 \exp\!\left(-(RM-D-1)\ell\right)
\le {C\over c_0\delta}
 \exp\!\left(-(RM-D-1)\delta\right)
\longrightarrow0.
\tag{19}
```

The convergence is uniform for all `ell>=delta`. This gives the
last contradiction for all sufficiently large `M`.

Sections 2--5 cover every positive scale. Every threshold depends only
on the fixed algebraic families and their constants, never on the
gauge vector or the profile scale. This proves the theorem.

## 6. Verification and the remaining construction problem

The old exact entry and dual certificates are the algebraic hypotheses;
they are not replaced by numerical profile fits. The accompanying
[analytic checker](../tools/check_all_scale_balanced_obstruction.py)
replays the shared dual data, checks the constants in the reduction,
and independently tests the finite Gaussian sum and tail inequalities.
Its recorded output is
[all-scale-balanced-obstruction.json](../evidence/all-scale-balanced-obstruction.json).

```sh
python tools/check_all_scale_balanced_obstruction.py \
  --output evidence/all-scale-balanced-obstruction.json
```

The all-scale assertion is the ordinary proof above. Finite numerical
controls do not establish an effective value of its eventual interval
threshold. No new Lean proof is claimed.

For these partial repairs, merely changing the balanced core weights
or sending the profile scale below the earlier cutoff cannot provide
an asymptotic post-half certificate. A constructive escape must change
an assumption: for example, use aggregate bulk charges that are not
balanced, change the correction so the low Gaussian restriction or
its gap duals change, or bound actual coefficient ranks beyond the
two-place Smith pressure used here.

This extension and its verification used substantial OpenAI assistant
support. The Gaussian crossing, prior partial repairs, exact entry
certificates, and dual-gap arguments retain the attribution in their
linked source notes.
