# Every digit count has the same interval pressure floor

## Result and scope

The truncated symmetric construction permits `p` states per mode in
characteristic `p`. Increasing this digit count does not improve the
displayed interval pressure, even when the digit count is allowed to
vary with the interval length.

**Theorem.** For every pair of integers `q>=2`, `n>=1`, and every real
`s>1`, the interval pressure at the original gauges is strictly positive
at `rho=1/sqrt(2)`. Thus this pressure certificate cannot cross the
half-rate threshold for any digit count or interval length.

Here `q` is an analytic digit-count parameter. When `q=p` is prime, the
[truncated symmetric Pascal construction](truncated-symmetric-pascal-crossings.md)
realizes the products below in characteristic `p`. The theorem is about
those products and their stated gauges. It does not exclude noninterval
modes, independent mode heights, non-Gaussian perturbations, or a better
bound on the actual coefficient ranks.

For each fixed `q`, the optimized gap also has a sharp asymptotic:

~~~math
\rho_{n,q}=
\frac1{\sqrt2}
+\frac{q\sqrt2}{16\pi^2n^2}
  \left(y_{n,q}^2+2y_{n,q}\right)+o(n^{-2}),
\qquad
y_{n,q}=W\!\left(\frac{12n^2}{q-1}\right).
\tag{1}
~~~

The error in (1) is asserted for fixed `q`. No uniform moving-`q`
asymptotic is needed for the first theorem, which already quantifies
every `q,n,s`.

The Gaussian profiles and endpoint are reused from predecessor
`18cec06deac47c8f142ea67a60e2278ac2234d67`. The binary arguments in
[the interval note](schur-one-species-interval.md) and
[the quantitative frontier](quantitative-interval-frontier.md) are
extended here to every digit count. The main factoring exponent remains
`1/2`. This extension used substantial OpenAI assistant support.

## 1. Exact products and their cancelled forms

Put

~~~math
L_q(x)=1+x+\cdots+x^{q-1},\qquad
f_q(x)=\log L_q(e^{-x}),\qquad R=4\rho.
~~~

At the natural interval gauges the two profile partitions and the
physical partition are

~~~math
\begin{aligned}
B_{n,q}(s)
 &=s^{(q-1)n^2}\prod_{j=1}^{n}L_q(s^{-(2j-1)})^4,\\
Z_{n,q,R}(s)
 &=\prod_{a=0}^{n-1}
   L_q\!\left(s^{-[8a-4(n-1)]/R}\right).
\end{aligned}
\tag{2}
~~~

A strict witness would have `B<Z^R`. Write `t=log s` and
`D=log B-R log Z`. The reciprocal identity
`L_q(x^(-1))=x^(-(q-1))L_q(x)` gives

~~~math
D_{n,q}(t,R)=
\begin{cases}
(q-1)t+4\displaystyle\sum_{j=1}^{n}f_q((2j-1)t)
-R\log q-2R\displaystyle\sum_{j=1}^{(n-1)/2}f_q(8jt/R),
&n\text{ odd},\\[3pt]
4\displaystyle\sum_{j=1}^{n}f_q((2j-1)t)
-2R\displaystyle\sum_{j=1}^{n/2}f_q(4(2j-1)t/R),
&n\text{ even}.
\end{cases}
\tag{3}
~~~

The function `f_q` is positive, decreasing, and convex on the positive
axis. Convexity follows directly from the log-sum-exp representation:
its second derivative is the variance of the digits `0,...,q-1` under
their Gibbs weights.

The critical value is `R0=2sqrt(2)`. Since

~~~math
\frac{\partial}{\partial R}
  \bigl(R\log Z_{n,q,R}(e^t)\bigr)
=\sum_a H(\pi_a)>0,
\tag{4}
~~~

positivity at `R0` also excludes every smaller `R`.

## 2. A periodic-series monotonicity lemma

For every integer `q>=2`, define

~~~math
a_r=1-q\,\mathbf1_{q\mid r}.
~~~

Its partial sums are `r mod q`, hence nonnegative. More explicitly,
if `b_r` is positive and strictly decreasing to zero, then

~~~math
\sum_{r\ge1}a_rb_r
=\sum_{k\ge0}\sum_{j=1}^{q-1}
  \bigl(b_{qk+j}-b_{qk+q}\bigr)>0
\tag{5}
~~~

whenever the series converges. This replaces the binary alternating
series without requiring `q` to be a power of two.

Set

~~~math
A_q(t)=\sum_{j\ge1}f_q((2j-1)t)
=\frac12\sum_{r\ge1}\frac{a_r}{r\sinh(rt)}.
\tag{6}
~~~

The second equality follows from the absolutely convergent logarithmic
series. Define

~~~math
h(x)=\frac{x\cosh x-\sinh x}{\sinh^2x}>0.
~~~

The function `h(x)/x` strictly decreases. In fact, its derivative is
`-F(x)/(x^2 sinh^3 x)`, where

~~~math
F(x)=(x^2-1)\sinh^2x-x\sinh x\cosh x+2x^2.
~~~

All lower coefficients in `2F` vanish, and the coefficient of `x^(2r)`
for `r>=2` is

~~~math
\frac{2^{2r-2}(4r^2-6r-4)}{(2r)!}.
\tag{7}
~~~

This is zero at `r=2` and positive for every `r>=3`. Therefore
`h(rt)/r` is a strictly decreasing positive sequence. Applying (5) gives

~~~math
\frac{d}{dt}\bigl(tA_q(t)\bigr)
=-\frac12\sum_{r\ge1}a_r\frac{h(rt)}r<0.
\tag{8}
~~~

The infinite even continuation of (3) at `R0` is consequently

~~~math
4\bigl(A_q(t)-\sqrt2 A_q(\sqrt2t)\bigr)>0.
\tag{9}
~~~

## 3. The odd continuation has an explicitly positive remainder

Let

~~~math
P(x)=\sum_{j\ge1}\log(1-e^{-jx}),\qquad
J_q(x)=-P(x/q)+P(x).
~~~

The function `J_q` is positive and decreasing: it is the sum of
`-log(1-exp(-jx/q))` over positive `j` not divisible by `q`.
The eta product and its modular transformation give

~~~math
\sum_{j\ge1}f_q(jt)
=\frac{I_q}{t}-\frac{\log q}{2}
 +\frac{(q-1)t}{24}-J_q(4\pi^2/t),
\qquad I_q=\frac{\pi^2}{6}\left(1-\frac1q\right).
\tag{10}
~~~

This follows by subtracting the eta transformation for `P(t)` from
the one for `P(qt)`. The primary identities are
[DLMF 27.14.12](https://dlmf.nist.gov/27.14.E12) and
[DLMF 23.18.5](https://dlmf.nist.gov/23.18.E5), already used by the
binary interval proof.

Replacing both finite sums in the odd line of (3) by infinite sums
and applying (10) yields

~~~math
\begin{aligned}
D_{\infty,q}^{\rm odd}(t,R)
={}&\frac{(8-R^2)I_q}{4t}+\frac{(q-1)t}{6}\\
&+4\bigl(J_q(2\pi^2/t)-J_q(4\pi^2/t)\bigr)
 +2R J_q(\pi^2R/(2t)).
\end{aligned}
\tag{11}
~~~

At `R=R0`, each remaining term is positive.

## 4. Finite tails, and the two small intervals

For either parity, the missing positive tail in (3) is bounded above by

~~~math
\sum_{j=n+1}^{\infty}f_q((2j-1)t)
\le\frac1{2t}\int_{2nt}^{\infty}f_q(x)\,dx.
\tag{12}
~~~

Convexity gives this midpoint bound. The missing negative tail,
multiplied by `R0/2=sqrt(2)`, is bounded below by

~~~math
\frac1{2t}\int_{\sqrt2(n+1)t}^{\infty}f_q(x)\,dx,
\tag{13}
~~~

using decreasing rectangles of the appropriate step. For `n>=3`,
`sqrt(2)(n+1)<=2n`, so restoring the finite truncation only increases
the critical pressure. Equations (9) and (11) prove the theorem for
every `n>=3`.

### Two modes

The even line of (3) gives

~~~math
\frac14D_{2,q}(t,R_0)
=f_q(t)+f_q(3t)-\sqrt2 f_q(\sqrt2t)
=\sum_{r\ge1}\frac{a_r}{r}g(rt),
\tag{14}
~~~

where `g(x)=exp(-x)+exp(-3x)-sqrt(2)exp(-sqrt(2)x)`.
This `g` is positive and strictly decreasing. To check the derivative,
put `a=sqrt(2)-1`, `b=3-sqrt(2)`, so `a+b=2`. Then

~~~math
-e^{\sqrt2x}g'(x)=e^{ax}+3e^{-bx}-2.
~~~

The minimum of the first two terms over the real line is
`(2/b)(3b/a)^(a/2)>2`. The strict comparison is elementary; one
fully rational certificate uses

~~~math
\frac{41}{100}<a<\frac{42}{100},\qquad
\frac{158}{100}<b<\frac{159}{100},\qquad
\left(\frac{79}{7}\right)^{41}
 >\left(\frac{159}{100}\right)^{200}.
\tag{15}
~~~

Taking logarithms in (15) gives `a log(3b/a)>2log b`.
Thus `g'<0`, and `g` tends to zero at infinity, proving its positivity.
The sequence `g(rt)/r` satisfies (5), so (14) is strictly positive.

### One mode

The one-mode pressure is `(q-1)t+4f_q(t)-R0 log q`.
For any probability distribution `pi` on the digits with mean
`(q-1)/4`, the Gibbs inequality gives

~~~math
(q-1)t+4f_q(t)\ge4H(\pi).
\tag{16}
~~~

Take integer weights

~~~math
w_j=\binom{q-j+1}{2},\quad 0\le j<q,
\qquad W=\binom{q+2}{3},\qquad \pi_j=w_j/W.
~~~

Their mean is `(q-1)/4`. Their entropy is strictly greater than
`(3/4)log q`. For `q=2,3,4,5`, the exact equivalent inequalities are

~~~math
W^{4W}>q^{3W}\prod_j w_j^{4w_j}.
\tag{17}
~~~

For `q>=6`, use Shannon entropy at least Renyi entropy of order two:

~~~math
H(\pi)\ge
\log\frac{W^2}{\sum_jw_j^2}
=\log\frac{5q(q+1)(q+2)}{3(3q^2+6q+1)}
>\frac34\log q.
\tag{18}
~~~

The last inequality is the positivity of

~~~math
625q(q+1)^4(q+2)^4-81(3q^2+6q+1)^4.
\tag{19}
~~~

After substituting `q=6+x`, every coefficient is a positive integer;
the checker records the complete coefficient list. Equations (16)--(18)
give `D_(1,q)>(3-R0)log q>0`. This completes the all-size theorem.

## 5. The sharp fixed-digit gap

Let `R_(n,q)(t)` be the unique root of `D_(n,q)(t,R)=0`, and set
`rho_(n,q)=(1/4)inf_(t>0)R_(n,q)(t)`. Uniqueness follows from (4),
with opposite signs at the small- and large-`R` limits.

For either parity, the eta transformations give, uniformly for `R` in
a fixed neighborhood of `R0`,

~~~math
D_{\infty,q}(t,R)
=\frac{(8-R^2)I_q}{4t}+\frac{(q-1)t}{6}
 +O(e^{-c_q/t}).
\tag{20}
~~~

For fixed `q`, `f_q(x)=exp(-x)+O(exp(-2x))` as `x` tends to infinity.
When `nt` tends to infinity, the same geometric tail calculation as
in the binary quantitative proof therefore gives

~~~math
\begin{aligned}
D_{n,q}(t,R)=D_{\infty,q}(t,R)
&+\frac{R^2}{4t}e^{-4nt/R}
  \bigl(1+O_q(t+e^{-4nt/R})\bigr)\\
&-\frac2t e^{-2nt}
  \bigl(1+O_q(t+e^{-2nt})\bigr).
\end{aligned}
\tag{21}
~~~

The endpoint arguments are unchanged: a root approaching `R0` must
have `t->0` and `nt->infinity`. If `nt->c` is finite and positive,
the Riemann-sum limit of `tD` is
`2 integral_(sqrt(2)c)^(2c) f_q(x) dx>0`; if `nt->0`, then
`D/n->(4-R)log q`. On scales bounded away from zero, the strict
infinite floor and the large-`t` leading terms exclude convergence to
`R0`. These arguments also confine minimizing sequences for large `n`.

Using a trial scale proportional to `log n/n` in (20)--(21), then
expanding `R=R0+delta`, gives throughout the possible minimizing range

~~~math
\frac\delta4
=\frac{q-1}{12R_0 I_q}t^2
 +\frac{R_0}{8I_q}e^{-\sqrt2nt}+o(n^{-2}).
\tag{22}
~~~

The two displayed terms form a strictly convex function. Its derivative
vanishes when

~~~math
(\sqrt2nt)e^{\sqrt2nt}=\frac{12n^2}{q-1}.
~~~

Substitution proves (1), and every minimizing scale satisfies
`t=y_(n,q)/(sqrt(2)n)+o(n^(-1))`. The positive leading gap is
`q(log n)^2/(2sqrt(2)pi^2 n^2)`. Thus, for each fixed digit count,
this displayed family approaches the same threshold from above.

## 6. Reproduction

~~~sh
python tools/check_truncated_symmetric_interval_frontier.py \
  --output evidence/truncated-symmetric-interval-frontier.json
~~~

The standard-library checker compares the original and cancelled
products using exact rational arithmetic, verifies the periodic block
coefficients, checks (7), proves the rational comparisons in (15) and
(17), and expands (19) independently. These are finite implementation
controls. The all-parameter pressure floor and asymptotic gap are the
ordinary proofs above; no new Lean claim is made.
