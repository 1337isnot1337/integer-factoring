# Every digit count has positive charged Gaussian pressure

## Result and scope

Let `G_{n,q}(J,t)` be the complete interval Gaussian critical pressure
from the [charged Gaussian note](charged-gaussian-pressure-floor.md),
including both cut partitions and the paired physical partition. The
critical contraction is `rho=1/sqrt(2)`.

**Theorem.** For every integer `n>=1`, real aggregate charge `J`,
positive scale `t`, and real `q>1`,

\[
G_{n,q}(J,t)>0.
\tag{1}
\]

For fixed `n,J,t`, the pressure is strictly increasing as a function
of the real digit-count parameter `q>0`.

For integer `q`, this is the original geometric digit partition.
In particular, the theorem includes every odd prime-state interval
Gaussian construction. Real noninteger `q` is an analytic extension,
not a physical alphabet. The previous
[global binary theorem](global-charged-gaussian-floor.md) reached
power-of-two digit counts through binary factorization; the result here
removes that restriction.

This excludes a strict critical-pressure witness in the stated
aggregate-charge Gaussian family. It does not exclude coupled filters,
independent physical gauges, new crossing tensors, or a better analysis
of actual coefficient ranks. It does not improve the factoring exponent,
which remains `1/2`.

The proof is ordinary analysis plus complete integer-outward compact
certificates. No new Lean formalization accompanies it. The construction
and pressure normalization are retained from predecessor
`c2e951fa1f391edd97ba587d15ec13d3d6c8892d`; this extension used substantial
OpenAI assistant support.

## 1. Differentiate in the digit count

Put `R=sqrt(2)` and extend the digit free energy to real `q>0` by

\[
f_q(x)=
\begin{cases}
\log\!\left((1-e^{-qx})/(1-e^{-x})\right),&x\ne0,\\
\log q,&x=0.
\end{cases}
\tag{2}
\]

For positive integer `q`, this equals
`log(1+e^{-x}+...+e^{-(q-1)x})`. Define the smooth kernel

\[
h(x)=\begin{cases}x/(e^x-1),&x\ne0,\\1,&x=0.\end{cases}
\qquad
\frac{\partial f_q(x)}{\partial q}=\frac{h(qx)}q.
\tag{3}
\]

For any function `k`, write the full finite pressure functional as

\[
\begin{aligned}
\mathcal D_n(k;J,t)={}&
\sum_{i=0}^{n-1}\{k((2i+1-J)t)+k((2i+1+J)t)\}\\
&+\sum_{a=0}^{2n-1}k((2a+1-2n)t)\\
&-R\sum_{i=0}^{n-1}\{k((4i-2n+2-J)t/R)
                         +k((4i-2n+2+J)t/R)\}.
\end{aligned}
\tag{4}
\]

Thus `G_{n,q}=D_n[f_q]`. Put `D_n=D_n[h]`. Equation (3) gives

\[
\boxed{\frac{\partial G_{n,q}(J,t)}{\partial q}
       =\frac{D_n(J,qt)}q.}
\tag{5}
\]

The main work below proves `D_n(J,t)>0` for every `n,J,t`. Since
`f_1=0`, integrating (5) from one to `q` proves (1), and (5) also
proves strict monotonicity. The formula at zero arguments is included;
no zero-energy physical factor is discarded.

## 2. Properties of the derivative kernel

Let `f(x)=log(1+e^{-x})` and

\[
p(x)=f''(x)=\frac{e^x}{(1+e^x)^2}.
\]

The substitution `u=e^y/(1+e^y)` gives the exact convolution identity

\[
\int_{\mathbb R}f(x-y)p(y)\,dy
=\int_0^1\log\!\left(1+e^{-x}\frac u{1-u}\right)du
=h(x).
\tag{6}
\]

Consequently `h'=f'*p` lies strictly between `-1` and zero, and
`h''=p*p>0`. In particular, `h` is positive, decreasing, and convex
on the entire real line. Directly from its formula,

\[
h(-x)=h(x)+x,
\qquad e^{-x}\le h(x)\le(1+x)e^{-x}\quad(x\ge0).
\tag{7}
\]

The lower thermal bound is equivalent to `x>=1-e^{-x}`. The upper
bound is equivalent to `e^x>=1+x`. Also `h(-x)>=max(x,1)` for `x>=0`.

The Fourier transform convention here is
`hat p(xi)=integral p(x)e^{-i xi x} dx`. Euler's beta integral,
the gamma recurrence, and reflection formula give

\[
\widehat p(\xi)
=\Gamma(1-i\xi)\Gamma(1+i\xi)
=\frac{\pi\xi}{\sinh(\pi\xi)},
\qquad
\widehat{h''}(\xi)=\left(\frac{\pi\xi}{\sinh(\pi\xi)}\right)^2.
\tag{8}
\]

The value at zero is the continuous limit. The primary identities are
[DLMF 5.12.1](https://dlmf.nist.gov/5.12.E1),
[5.12.3](https://dlmf.nist.gov/5.12.E3),
[5.5.1](https://dlmf.nist.gov/5.5.E1), and
[5.5.3](https://dlmf.nist.gov/5.5.E3). Equation (6) also supplies
exponential decay of `h''` and its derivatives, sufficient for the
Poisson summation used next.

## 3. A positive periodic continuation

Define

\[
A(J,t)=\sum_{j\ge1}\{h((2j-1+J)t)+h((2j-1-J)t)\},
\qquad
H(J,t)=A(J,t)+A(0,t)-2R A(J/2,Rt).
\tag{9}
\]

These series converge for every fixed charge and positive scale.
Equation (7) gives `A(J+2,t)=A(J,t)+(J+1)t`; hence `H` is even and
four-periodic in `J`. It is enough to prove positivity on `0<=J<=2`.

### Small scales

Apply [Poisson summation](https://dlmf.nist.gov/1.8.E14) to `h''` on
the shifted odd lattice, using (8), and integrate twice in `J`.
The exact resulting Fourier series is

\[
\boxed{
A(J,t)=\frac{\pi^2}{6t}-\frac t{12}+\frac{tJ^2}4
-\frac{\pi^2}t\sum_{m\ge1}(-1)^m\cos(\pi mJ)
                   \frac1{\sinh^2(\pi^2m/t)}.}
\tag{10}
\]

Here is the constant-term calculation, so (10) has no undetermined
integration term. The function `A(J,t)-tJ^2/4` is even and two-periodic.
Integrating `A` on `[0,2]` concatenates one interval `[-t,t]` and
two copies of `[t,infinity)`. Using (7) gives

\[
\frac12\int_0^2 A(J,t)dJ
=\frac1t\int_0^\infty h(x)dx+\frac t4
=\frac{\pi^2}{6t}+\frac t4.
\]

The integral follows by integrating `x sum_{k>=1}e^{-kx}` and using
[zeta(2)=pi^2/6](https://dlmf.nist.gov/25.6.E1).
The mean of `tJ^2/4` is `t/3`, giving the constant in (10).
The cosine coefficients follow directly by dividing those of `A''`
by `-pi^2 m^2`. All Fourier series in this argument converge absolutely.

Substitution in (9) cancels the inverse-scale and charge-quadratic
terms. The nonoscillatory remainder is `t/6`. For `a>0`,

\[
\sum_{m\ge1}\frac1{\sinh^2(am)}
\le\frac{4e^{-2a}}{(1-e^{-2a})^3}.
\tag{11}
\]

Use `pi^2<10`, `pi>3`, and `R<3/2`. For `0<t<=1`, the absolute
Fourier error in `H` is strictly less than

\[
\frac{160}t\frac{e^{-12/t}}{(1-e^{-12/t})^3}<\frac t{900}.
\tag{12}
\]

For the last comparison, `e^3>20`, so `e^{12}>160000`;
`t^{-2}e^{-12/t}` increases on `(0,1]`. The remaining rational
comparison is

\[
\frac{160/160000}{(1-1/160000)^3}<\frac1{900}.
\]

Therefore `H(J,t)>149t/900>0` on this entire small-scale range.

### Large scales

Let `t>=16`, `0<=J<=2`, and set
`s=R(1-J/2)t`, `d=2Rt`, and `z=e^{-d}`. Monotonicity and (7) give

\[
2R A(J/2,Rt)
\le4R e^{-s}\left(\frac{s+1}{1-z}
                         +\frac{dz}{(1-z)^2}\right)
<6(s+1)e^{-s}.
\tag{13}
\]

To check the last constant, `d>44`, `z<eta=20^{-14}`, and
`dz<44eta`, since `d e^{-d}` decreases for `d>1`. The exact
inequality `1/(1-eta)+44eta/(1-eta)^2<21/20`, together with
`4R<40/7`, proves (13).

If `J<=1`, the positive term `A(J,t)` is at least
`e^{-(1-J)t}`. The ratio of (13) to this term is less than
`6(1+3t/2)e^{-2t/5}`. This decreases for `t>=16` and is at most
`150/e^6<150/400<1`.

If `J=1+v` with `0<=v<=1`, then `A(J,t)>=max(vt,1)` and
`s=(1-v)t/R`. For `v>=3/5`, this positive term exceeds nine,
whereas (13) is at most six. For `v<=3/5`, use `s>4t/15`:
the negative term is less than `6(1+4t/15)e^{-4t/15}<1` for
`t>=16`. At sixteen this follows from `e^4>50`; the expression
decreases afterward. Thus `H>0` at every large scale.

### Compact scales

The remaining rectangle `[1,16] x [0,2]` has a complete cover by
**803 positive outward-dyadic rectangles** at precision eighty.
The checker retains twelve terms in each series. Each missing
negative tail is bounded with the exact geometric sum in (13),
using its actual first argument and spacing; positive missing tails
are discarded only for a lower bound.

The minimum certified lower numerator is `4902945637190845203`,
with common denominator `2^80`. Every leaf and subdivision path is
recorded in the [certificate](../evidence/all-digit-charged-gaussian-floor.json).
The independent tree audit verifies complete coverage and reconstructs
each lower bound. Together with (12)--(13), this proves

\[
H(J,t)>0\qquad(J\in\mathbb R,\ t>0).
\tag{14}
\]

## 4. Two scalar kernels needed at exterior and edge charges

For every `u>0`,

\[
h(u)+h(3u)-R h(Ru)>0,
\qquad
h(u)+h(3u/2)-R h(3u/(2R))>0.
\tag{15}
\]

Both have elementary exponential proofs. For the first, the expansion
`h(x)=x sum_{k>=1}e^{-kx}` reduces positivity to

\[
e^{-v}+3e^{-3v}-2e^{-Rv}>0\qquad(v>0).
\]

Set `a=R-1`, `b=3-R`, so `a+b=2`. The minimum over the real line
of `e^{av}+3e^{-bv}` is `(2/b)(3b/a)^{a/2}`. It exceeds two:
use `41/100<a<42/100`, `158/100<b<159/100`, and the exact integer
comparison `(79/7)^41>(159/100)^200`. This is also the exponential
kernel in the [natural all-digit proof](truncated-symmetric-interval-frontier.md).

For the second inequality, put `c=3/(2R)`, `a=c-1`, and `b=3/2-c`.
Now `a+b=1/2`. It suffices that
`e^{av}+(3/2)e^{-bv}>3/2`. Its minimum is
`(1/(2b))(3b/(2a))^{2a}`. The bounds

\[
\frac3{50}<a<\frac{61}{1000},\qquad
\frac{439}{1000}<b<\frac{11}{25},\qquad
\frac{3b}{2a}>\frac{107}{10},\qquad
\left(\frac{107}{10}\right)^3>
\left(\frac{33}{25}\right)^{25}
\tag{16}
\]

show `(3b/(2a))^{2a}>33/25>3b`. The radical bounds in (16)
follow by squaring `106/75<R<1061/750`. All comparisons in this
paragraph are included as exact rational controls in the checker.

## 5. The three finite residual bases

The complete cancelled profiles in the
[small-base note](small-charged-gaussian-bases.md) remain valid with
`h` replacing `f`, because their derivation uses only
`h(-x)=h(x)+x`. They give

\[
D_n(J,t)=g(J)t+
\sum_{i=0}^{n-1}\{h(|2i+1-J|t)+h((2i+1+J)t)+2h((2i+1)t)\}
-2R\sum_{e\in\mathcal E}h(et/R),
\tag{17}
\]

with the exact data

| Length and charge | `g(J)` | Physical energies `E` |
|---|---|---|
| `n=1`, `0<=J<=2` | `(1-J)_+` | `J` |
| `n=2`, `0<=J<=2` | `(J-1)_+` | `2-J,2+J` |
| `n=2`, `2<=J<=4` | `(3-J)_+` | `J-2,J+2` |
| `n=3`, `0<=J<=2` | `(1-J)_+` | `J,4-J,4+J` |

At zero scale, `D_n(J,0)=2n(2-R)>8n/7`. Since `-1<h'<0`, the
original finite energy lists give derivative lower bound
`-(4n^2+8)` on these four charge ranges. Hence for `t<=1/16`,

\[
D_n(J,t)>\frac{8n}7-\frac{4n^2+8}{16}>0
\quad(n=1,2,3).
\tag{18}
\]

For `t>=16`, use the retained positive term `h(|1-y|t)`, where
`y=J` except in the third row, where `y=J-2`. In the half with
zero ground coefficient, put `a=|1-y|`, and let `b` be the least
physical energy divided by `R`. The displayed lists give
`b-a>=R-1>2/5` and `b<=R<3/2`. The ratio of the entire negative
contribution to the retained positive term is less than
`9(1+3t/2)e^{-2t/5}<=225/e^6<1`.

In the other half put `v=|1-y|`. The ground and retained term
are at least `max(vt,1)`, while the negative contribution is at
most `9(1+bt)e^{-bt}`, with `b=(1-v)/R`. If `v>=3/5`, the
positive term exceeds nine. If `v<=3/5`, the negative term is
less than `9(1+4t/15)e^{-4t/15}<1`; at sixteen use
`711/15<50<e^4`. These bounds decrease for larger scales.

The remaining four rectangles, with scale `[1/16,16]`, have the
following complete positive integer-outward covers:

| Base | Positive rectangles | Minimum lower numerator, denominator `2^80` |
|---|---:|---:|
| `n=1`, `[0,2]` | 213 | 4772837920546723763 |
| `n=2`, `[0,2]` | 467 | 4653847454689935807 |
| `n=2`, `[2,4]` | 258 | 4772837920546888911 |
| `n=3`, `[0,2]` | 907 | 474969336081336921 |

Equation (17) retains all physical and cut terms, including zero
energies. Thus the four bases have positive pressure at every scale.

## 6. Restore every finite length and charge

The previous folding, tail, and edge identities depend on the kernel
only through its reciprocal identity, monotonicity, and convexity.
Section 2 verifies all three properties for `h`.

First, the [finite-tail majorization](folded-charged-gaussian-floor.md)
restores the finite interval from (14), giving

\[
D_n(J,t)>0\quad\text{if}\quad
n\ge3,\quad |J|\le2n-2-2R.
\tag{19}
\]

Explicitly, put `I(x)=integral_x^infinity h(u)du`. This function is
decreasing and convex. The physical-tail endpoints
`(A-B,A-B,A+B,A+B)`, with `A=R(n+1)` and `B=|J|/R`, have
no larger ascending prefix sums than the cut endpoints
`(2n-|J|,2n,2n,2n+|J|)` precisely on the sufficient band (19).
The decreasing-convex comparison proves that restoring the finite
tails cannot decrease the positive continuation. The continuation
charge is `J+2(n mod 2)`; its four-periodicity covers either parity.

At `n>=4`, (19) contains `[0,2]`. Section 5 supplies this residual
charge range for lengths two and three. The exact folding inequality is

\[
D_n(2k+\delta,t)\ge
\frac{D_{n+k}(\delta,t)+D_{n-k}(\delta,t)}2,
\quad0\le k\le n-2,\quad0\le\delta\le2.
\tag{20}
\]

Its correction terms are nonnegative by the same decreasing-convex
argument in the folding note. Therefore `D_n(J,t)>0` for
`n>=2`, `0<=J<=2n-2`.

For the last interior band write `J=2n+a`, `-2<=a<=0`. The exact
energy cancellation gives

\[
D_n(2n+a,t)=\sum_{i=0}^{n-1}
\{2h((2i+1)t)+h((a+4i+1)t)+h((a+4i+3)t)
-2R h((a+4i+2)t/R)\}.
\tag{21}
\]

The first two summands are `D_2(4+a,t)>0`. For `i>=2`, put
`u=(2i+1)t` and `alpha=2+a/(2i+1)>=8/5`. Convexity bounds
the summand below by
`2[h(u)+h(alpha u)-R h(alpha u/R)]`. This increases with
`alpha`, since `h'` increases, and is positive already at
`alpha=3/2` by (15). Length one is covered by Section 5.

For the exterior range `a>=0`, the
[exterior reduction](exterior-charged-gaussian-floor.md) applies to
(21). Each summand is bounded below by
`2[h(u)+h(2u)-R h(Ru)]` for the appropriate positive `u`,
which is larger than the first expression in (15). Thus every
exterior charge is also positive. Evenness supplies negative charges.

This completes the proof `D_n(J,t)>0` for every integer `n>=1`,
real `J`, and `t>0`. Equation (5) proves the stated all-digit theorem.

## 7. Reproduction and evidence

Run

```sh
python tools/check_all_digit_charged_gaussian_floor.py \
  --output /tmp/all-digit-gaussian-replay.json
```

The [kernel bounds](../tools/all_digit_charged_gaussian_floor.py) reuse
the repository's outward-integer exponential arithmetic. At a nonzero
nonnegative dyadic argument they evaluate `x e^{-x}/(1-e^{-x})`,
with outward division; the zero value is exactly one. Monotonicity
then bounds a whole interval from its endpoints. Negative arguments
retain the exact reciprocal linear term. Insufficient precision near
zero raises an error instead of weakening the bound.

The [checker](../tools/check_all_digit_charged_gaussian_floor.py)
regenerates **2,648 positive rectangles** and independently audits
their **5,291 tree nodes**. Deliberately deleting a leaf or increasing
a stored lower endpoint is rejected. It also checks the twenty-one
rational tail and scalar inequalities, thirty-five exact signed-energy
profiles, and 144 outward overlaps between the original complete
pressure and its independently cancelled expression. Source hashes
are included in the [record](../evidence/all-digit-charged-gaussian-floor.json).

The compact certificates prove their displayed continuum domains;
the analytic arguments cover all unbounded scale, charge, length,
and digit-count parameters. No floating-point search or finite-rank
extrapolation is used in the positivity theorem.
