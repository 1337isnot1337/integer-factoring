# A charge interval and loss cap for every one-mode binary Laurent filter

## Results

The [zero-charge Laurent theorem](one-mode-binary-laurent-filter-floor.md)
extends to a substantial interval of nonzero aggregate charges.
Use the same construction, arbitrary finite Laurent occupation weights,
fixed characteristic-two coefficient field, and nonzero endpoint repair.
For selected index \(0\le j<n\), put \(\lambda=2n-1-2j\).

**Theorem.** The ratio \(R\) of the filtered two-cut cost to the
Gaussian cost at the same charge and scale satisfies

\[
R\ge1\qquad\text{whenever}\qquad |J|\le2n.
\tag{1}
\]

Within this charge interval, equality holds exactly for a common scalar
Laurent monomial on both occupation states. The statement holds at every scale and degree,
including arbitrary multilevel cancellations. The interval is uniform
over the selected index. Its negative-charge part uses an additional
coefficient constraint beyond the degree/order dichotomy.

There is also a charge-independent strict loss cap:

\[
\boxed{R>\frac1{1+e^{-\lambda t}}\quad(J\in\mathbb R,\ t>0).}
\tag{2}
\]

The lower bound in (2) is sharp as an infimum over charges and finite
Laurent filters. Relative improvements do occur outside (1); exact
examples are given below. A stronger binary exterior Gaussian estimate
also gives an absolute conclusion:

**Exterior corollary.** Every filter in this one-mode Laurent family
has strictly positive full critical pressure whenever \(|J|\ge2n\),
for every admissible interval length and every scale. This includes the explicit
relative-improvement examples.

## 1. Exact probabilities and three profile regimes

Write \(m=nP\), where \(P\) is a positive odd integer, and
\(f_0=\alpha_1,f_1=\alpha_0\). As before,
\(\beta_0=f_0\), \(\beta_1=f_1+(1+x^m)f_0\).
A common scalar monomial cancels from the cost, so normalize
\(\mathrm{ord}_x f_0=0\). Let

\[
b=\mathrm{ord}_x\beta_1,\quad D_i=\deg f_i,\quad
\eta=e^{2t/P},\quad H=\eta^m=e^{2nt},\quad z=\eta^{m-b}.
\tag{3}
\]

The zero-cut and infinity probabilities are, respectively,

\[
p_1=\frac1{1+e^{(\lambda+J)t}},\qquad
q_1=\frac1{1+e^{\lambda t}},\qquad
p_0=1-p_1,\quad q_0=1-q_1.
\tag{4}
\]

All four probabilities are positive, and \(q_0>q_1\). The exact
cost is

\[
R=(p_0+p_1z)(q_0\eta^{D_0}+q_1\eta^{D_1}).
\tag{5}
\]

Always \(D_0\ge0\) and \(D_1\ge\min(0,b)\): if the order of
\(f_1\) differed from zero, the smaller input order would survive
in \(\beta_1\). Also
\(b\le\max(D_1,D_0+m)\). Define

\[
A(z)=(p_0+p_1z)(q_0/z+q_1),\qquad
B(z)=(p_0+p_1z)(q_0+q_1H/z).
\tag{6}
\]

The full degree/order consequences are:

| Beta order | Range of \(z\) | Lower bound |
|---|---|---|
| \(b\le0\) | \(z\ge H\) | \(R\ge B(z)\) |
| \(0\le b\le m\) | \(1\le z\le H\) | \(R\ge p_0+p_1z\ge1\) |
| \(b>m\), \(D_1\ge b\) | \(0<z<1\) | \(R\ge B(z)\) |
| \(b>m\), \(D_1<b\) | \(0<z<1\) | \(R\ge A(z)\) |

In the last row, cancellation gives \(\mathrm{ord}_x f_1=0\)
and \(D_0\ge b-m\). The table uses complete degrees at infinity;
there is no uncharged high-degree term.

## 2. Minimizing the two elementary bounds

The unique stationary point of \(B\) is

\[
z_B=\sqrt{\frac{p_0q_1H}{p_1q_0}}
    =e^{(J+2n)t/2}.
\tag{7}
\]

The function decreases before this point and increases afterward.
If \(J\le2n\), then \(z_B\le H\), so on \(z\ge H\),
\(B(z)\ge B(H)=p_0+p_1H>1\).
If \(J\ge-2n\), then \(z_B\ge1\), so on \(0<z\le1\),
\(B(z)\ge B(1)=q_0+q_1H>1\).

On \(0<z\le1\),

\[
A'(z)=-p_0q_0/z^2+p_1q_1.
\]

Thus \(A(z)\ge A(1)=1\) whenever
\(p_0q_0\ge p_1q_1\). By (4), this is exactly
\(J\ge-2\lambda\). This already proves domination on
\([-2\min(n,\lambda),2n]\). The following coefficient argument
fills the remaining negative-charge range.

In the last row of the table, \(b>m\) and both input orders are zero.
Necessarily

\[
\max(D_0,D_1)\ge m.
\tag{7a}
\]

Indeed, if \(D_0<m\), the coefficient of \(x^m\) in
\((1+x^m)f_0\) is the nonzero constant coefficient of \(f_0\).
Since \(\beta_1\) has order greater than \(m\), \(f_1\) must
cancel that coefficient, forcing \(D_1\ge m\).

If \(D_0\ge m\), the lower bound is strengthened to

\[
R\ge(p_0+p_1z)\,[q_0\max(H,1/z)+q_1].
\tag{7b}
\]

For \(z\ge1/H\), this expression increases with \(z\). For
\(z\le1/H\), it equals \(A(z)\) and decreases with \(z\),
because its stationary point
\(\sqrt{p_0q_0/(p_1q_1)}=e^{(2\lambda+J)t/2}\)
exceeds \(1/H\) when \(J\ge-2n\). The bound is therefore at least
\(A(1/H)>1\); the strict inequality is equivalent to
\(Hp_0q_0>p_1q_1\), and follows from
\(2n+2\lambda+J>0\).

If \(D_1\ge m\), we instead obtain

\[
R\ge(p_0+p_1z)(q_0/z+q_1H)
\ge(\sqrt{Hp_0q_1}+\sqrt{p_1q_0})^2.
\tag{7c}
\]

When \(J\ge-2n\), the ratio \(Hp_0q_1/(p_1q_0)\) equals
\(e^{(2n+J)t}\ge1\). In the only previously unproved range
\(J<-2\lambda\), one has \(p_1>q_0>1/2\). Thus (7c) is at
least \(4p_1q_0>1\). Together with the earlier \(A\) bound this
proves (1) throughout the full interval \([-2n,2n]\).

Equality in the middle regime requires \(z=1\) and
\(D_0=D_1=0\). The outside regimes controlled by \(B\) are strict.
For the \(A\) regime, \(A(z)>1\) for \(z<1\), even at the
charge boundary where \(A'(1)=0\). Hence equality again forces
\(f_0=f_1\) constant after normalization, which is the common scalar
Laurent monomial stated in the theorem.

For reference, unrestricted continuous minimization gives

\[
\min_{z>0}B(z)=
(\sqrt{p_0q_0}+\sqrt{Hp_1q_1})^2,
\tag{8}
\]

and, when \(p_1+q_1>1\),

\[
\min_{0<z\le1}A(z)=
(\sqrt{p_0q_1}+\sqrt{p_1q_0})^2.
\tag{9}
\]

These are valid lower envelopes. Actual Laurent orders impose a
discrete set of \(z\), and endpoint cancellation can impose further
restrictions. No assertion that every continuous minimum is attained
by a physical filter is needed.

## 3. A universal strict loss cap

For \(z\le1\),

\[
A(z)=q_0(p_0/z+p_1)+q_1(p_0+p_1z)>q_0.
\]

For every \(z>0\),

\[
B(z)>p_0q_0+p_1Hq_1>q_0,
\]

because \(Hq_1/q_0=e^{(2n-\lambda)t}>1\). The middle regime has
\(R\ge1>q_0\). This proves (2), or equivalently

\[
\log R>-\log(1+e^{-\lambda t}).
\tag{10}
\]

For sharpness, take \(f_0=1,f_1=x^{-d}\). Its cost is
\((p_0+p_1H\eta^d)(q_0+q_1\eta^{-d})\).
First send \(J\) to positive infinity at fixed \(d\), then send
the finite positive integer \(d\) to infinity. The ratio tends to
\(q_0\). This is an infimum argument; each member still has finite
ordinary degree after monomial clearing and nonzero endpoint values.

## 4. Absolute exclusions and exact outside examples

Combining (1) with the
[charged Gaussian floor](charged-gaussian-pressure-floor.md) excludes
a critical witness throughout

\[
|J|\le C_n,\qquad n\ge3.
\]

It also excludes every even integer charge in \([-2n,2n]\),
using the all-digit even-charge Gaussian theorem, for every selected
index and every admissible binary interval length.

Relative domination cannot be extended to all real charges. With
\(n=2,P=1,t=\log2\), the following endpoint-one examples have
exactly computed ratios:

| Selected index and charge | Laurent weights \((f_0,f_1)\) | Exact ratio |
|---|---|---|
| \(j=0,J=8\) | \((1,x^{-1})\) | \(1936/2049<1\) |
| \(j=1,J=-12\) | \((1+x^2+x^4,1)\) | \(513/683<1\) |

They illustrate the scope of the relative theorem. Section 5 proves
that both have positive full critical pressure.

The exact rational replay checks the envelopes across physical charge
bands and outside-charge loss-cap controls, and replays both complete
polynomial examples:

```sh
python tools/check_charged_one_mode_filters.py
```

Its counts are recorded in
[`charged-one-mode-filters.json`](../evidence/charged-one-mode-filters.json).

## 5. Positive full pressure at every exterior charge

The [exterior Gaussian identity](exterior-charged-gaussian-floor.md)
gives, for \(J=2n+a\) with \(a\ge0\),

\[
G_n(J,t)=\sum_{i=0}^{n-1}\left[
2f((2i+1)t)+f((a+4i+1)t)+f((a+4i+3)t)
-2\sqrt2 f((a+4i+2)t/\sqrt2)\right],
\tag{11}
\]

where \(f(x)=\log(1+e^{-x})\). Each summand is strictly positive
by the two-mode Gaussian floor. We sharpen the first one using

\[
\boxed{2f(u)+f(3u)>2\sqrt2 f(\sqrt2u)\qquad(u>0).}
\tag{12}
\]

Here is a complete analytic-plus-rational proof of (12). Let its
left side minus right side be \(g(u)\). For \(0<u\le1/32\),
\(g(0)=(3-2\sqrt2)\log2>2/21\), using
\(\sqrt2<10/7\) and \(\log2>2/3\). Also \(g'(u)\ge-5/2\),
since the two negative derivative terms have combined magnitude at
most \(5/2\), while the third derivative term is positive. Thus
\(g(u)>2/21-5/64>0\).

For \(u\ge5/2\), use
\(2f(u)\ge2(\log2)e^{-u}\) and
\(2\sqrt2 f(\sqrt2u)\le2\sqrt2 e^{-\sqrt2u}\).
The ratio of the second bound to the first is below
\((9/4)e^{-2u/5}\le9/(4e)<1\), since \(e>5/2\).
The remaining interval \([1/32,5/2]\) has a complete 70-interval
outward-dyadic certificate, using the already reviewed
[`Dyadic`](../tools/charged_gaussian_floor.py) arithmetic. Every
subinterval has a positive exact lower bound; the smallest is

\[
\frac{259466609285209734441}{1208925819614629174706176}>0.
\tag{13}
\]

The checker constructs the full subdivision tree and verifies every
bound with integer arithmetic. Its compact paths and bounds are stored
with the rational probability controls. There is no truncated infinite
series in this scalar certificate.

To apply (12), set \(x=(a+2)t\). For fixed \(x\), the function
\(f(v)+f(x-v)+f(x+v)\) is strictly decreasing on
\(0<v\le x/2\), because its derivative is
\(f'(v)-f'(x-v)+f'(x+v)<0\). Therefore the first summand in (11),
after subtracting \(f(t)\), is at least
\(g(x/2)>0\). All remaining summands are positive. Evenness in
\(J\) now yields the stronger Gaussian bound

\[
\boxed{G_n(J,t)>f(t)\qquad(|J|\ge2n).}
\tag{14}
\]

For a one-mode filter, its full pressure is
\(F_n=G_n+\log R\). Since \(\lambda\ge1\), (10) and (14) give

\[
F_n(J,t)>f(t)-f(\lambda t)\ge0
\qquad(|J|\ge2n).
\tag{15}
\]

This proves the exterior corollary, including \(n=1\). The middle
charge region between the continuous Gaussian band and the exterior
range remains outside this absolute theorem except for the even
integer charges covered above.

The subsequent [global Gaussian theorem](global-charged-gaussian-floor.md)
closes the remaining middle-charge gap. Combined with the relative
interior theorem and the exterior corollary above, it excludes this
one-mode family at every real aggregate charge. The filter construction
retains its admissible power-of-two interval lengths.
