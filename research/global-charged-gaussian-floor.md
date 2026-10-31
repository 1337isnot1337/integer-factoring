# The binary interval Gaussian has positive pressure at every real charge

## Result and scope

For every integer length \(n\ge1\), every real aggregate charge
\(J\), and every scale \(t>0\), the binary interval Gaussian has
strictly positive complete critical pressure:

\[
\boxed{G_n(J,t)>0.}
\tag{1}
\]

The critical parameter is \(\rho=1/\sqrt2\). This is the same
two-cut pressure, including the paired physical partition, as in the
[charged Gaussian theorem](charged-gaussian-pressure-floor.md). It
does not assert that zero charge minimizes the pressure; that stronger
statement is false. Nor does it cover arbitrary independent physical
state gauges. It resolves the aggregate-charge family at every length,
with no prime-power restriction on this analytic pressure statement.

Combined with the
[charged one-mode Laurent theorem](charged-one-mode-binary-laurent-filters.md),
(1) also excludes every admissible binary one-mode Laurent occupation
filter in that theorem, at every real aggregate charge and every scale.
This is an obstruction within those families, not an improved factoring
bound or an obstruction to unrelated constructions.

The proof combines exact charge identities, analytic inequalities, and
complete outward-integer certificates for compact scalar domains.
No floating-point sign decision is used.

## 1. Normalization and the established ingredients

Write \(R=\sqrt2\) and \(f(x)=\log(1+e^{-x})\), for real \(x\).
The full pressure is

\[
\begin{aligned}
G_n(J,t)={}&\sum_{i=0}^{n-1}
 [f((2i+1-J)t)+f((2i+1+J)t)]\\
&+\sum_{a=0}^{2n-1}f((2a+1-2n)t)\\
&-R\sum_{i=0}^{n-1}
 [f((4i-2n+2-J)t/R)+f((4i-2n+2+J)t/R)].
\end{aligned}
\tag{2}
\]

It is even in \(J\), so it suffices to consider \(J\ge0\).
The proof uses three previously checked facts:

1. The [exterior floor](exterior-charged-gaussian-floor.md) proves
   \(G_n(J,t)>0\) for \(J\ge2n\), every \(n\ge1\).
2. The [folding and tail theorem](folded-charged-gaussian-floor.md)
   proves positivity for \(|J|\le2n-2-2R\), \(n\ge3\), and
   the exact inequality
   \[
   G_n(2k+\delta,t)\ge
   \frac{G_{n+k}(\delta,t)+G_{n-k}(\delta,t)}2,
   \quad 0\le k\le n-2,\quad0\le\delta\le2.
   \tag{3}
   \]
3. The [small-base theorem](small-charged-gaussian-bases.md) proves
   positivity for length one on \([0,2]\), length two on \([0,4]\),
   and length three on \([0,2]\), at every scale.

For every length \(m\ge4\), the interval \([0,2]\) lies inside
the tail band because \(2m-2-2R\ge6-2R>2\). Consequently

\[
G_m(\delta,t)>0\qquad(m\ge2,\ 0\le\delta\le2).
\tag{4}
\]

## 2. Folding reaches within two units of the exterior boundary

For \(n\ge2\) and \(0\le J\le2n-2\), choose
\(J=2k+\delta\) with \(0\le k\le n-2\) and
\(0\le\delta\le2\). At the right endpoint use
\(k=n-2,\delta=2\); this avoids the inadmissible lower length one.
Both lengths in (3) are at least two. Equation (4) therefore gives

\[
G_n(J,t)>0\qquad(0\le J\le2n-2).
\tag{5}
\]

Only the edge interval \(2n-2\le J\le2n\) remains.

## 3. A scalar kernel for the edge

For every \(u>0\),

\[
\boxed{h(u):=f(u)+f(3u/2)-R f(3u/(2R))>0.}
\tag{6}
\]

For \(0<u\le1/4\), use \(R<10/7\) and \(\log2>2/3\):
\(h(0)=(2-R)\log2>8/21\). On nonnegative arguments,
\(-1/2\le f'\le0\), and the differentiated negative term in
\(h\) is nonnegative. Hence \(h'(u)\ge-5/4\), giving

\[
h(u)>8/21-5/16>0.
\tag{7}
\]

For \(u\ge20\), concavity of \(\log(1+z)\) on \([0,1]\)
gives \(f(u)\ge(\log2)e^{-u}\). Also \(f(x)\le e^{-x}\).
The ratio of the negative term in (6) to this retained positive term
is therefore less than

\[
\frac{15}{7}e^{-u/20}\le\frac{15}{7e}<1,
\tag{8}
\]

using \(3/(2R)>21/20\), \(R/\log2<15/7\), and
\(e>1+1+1/2+1/6>15/7\). The second positive term may be discarded.

The compact interval \([1/4,20]\) is covered by 250 positive
outward-dyadic intervals at precision 80. The smallest certified lower
numerator is \(476138769828760\), with common denominator \(2^{80}\).
The complete subdivision tree and every integer lower bound appear in
the [certificate](../evidence/global-charged-gaussian-floor.json).
The arithmetic is the already reviewed integer Taylor implementation
in [`charged_gaussian_floor.py`](../tools/charged_gaussian_floor.py).
Thus (7)--(8) and the finite certificate prove (6) for every positive
real argument.

## 4. The edge is a two-mode base plus positive summands

Write \(J=2n+a\). The following identity is valid for every real
\(a\), including negative values:

\[
G_n(2n+a,t)=\sum_{i=0}^{n-1}g_i(a,t),
\tag{9}
\]

where

\[
g_i(a,t)=2f((2i+1)t)+f((a+4i+1)t)+f((a+4i+3)t)
 -2R f((a+4i+2)t/R).
\tag{10}
\]

For completeness, apply \(f(-x)=f(x)+x\) to the first zero-cut
list in (2), then reverse its indices. The two zero lists concatenate
into \(f((a+2i+1)t)\), \(0\le i<2n\). Their linear term is
\((n^2+na)t\). The infinity list contributes
\(2\sum_{i<n}f((2i+1)t)+n^2t\). Reversing one physical list
makes both physical lists equal to the last list in (10), and their
linear term is \((2n^2+na)t\). All linear terms cancel. The
reciprocal identity holds for every real argument, so no positivity
assumption on \(a\) was made.

Now let \(-2\le a\le0\). The first two terms of (9) are exactly
\(G_2(4+a,t)>0\), by the small-base theorem. For \(i\ge2\), put

\[
b=2i+1\ge5,\qquad u=bt>0,\qquad
\alpha=2+a/b\ge8/5>3/2.
\tag{11}
\]

Convexity gives

\[
g_i(a,t)\ge
2\left[f(u)+f(\alpha u)-R f(\alpha u/R)\right].
\tag{12}
\]

For fixed \(u>0\), the bracket increases with \(\alpha\), because
its derivative is

\[
u[f'(\alpha u)-f'(\alpha u/R)]\ge0.
\tag{13}
\]

Here \(f'\) increases and \(\alpha u\ge\alpha u/R\).
Equations (6) and (11)--(13) prove every remaining summand positive.
This proves the edge floor for \(n\ge2\). Length one is already
covered by its small-base theorem and the exterior floor.
Together with (5) and evenness, the proof of (1) is complete.

## 5. Absolute exclusion of the one-mode Laurent family

Let \(\mathcal R\) be the ratio of filtered to Gaussian two-cut cost
in the [charged one-mode theorem](charged-one-mode-binary-laurent-filters.md).
Its complete filtered critical pressure is
\(G_n(J,t)+\log\mathcal R\).
For \(|J|\le2n\), that theorem gives \(\mathcal R\ge1\), and
(1) makes the full pressure strictly positive.

For \(|J|\ge2n\), the same theorem gives the stronger exterior
bound \(G_n(J,t)>f(t)\) and the strict loss cap
\(\log\mathcal R>-f(\lambda t)\), where
\(\lambda=2n-1-2j\ge1\). Since \(f\) decreases,
\(f(\lambda t)\le f(t)\), so the full pressure is again positive.
This covers every real charge, arbitrary finite Laurent support,
every selected mode, and every fixed characteristic-two coefficient
field allowed by that construction. The filter corollary retains the
construction's admissible power-of-two interval lengths; the analytic
Gaussian theorem itself holds at every integer length. The existing
endpoint and degree conditions remain in force.

## 6. Analytic power-of-two digit counts

For a digit count \(q=2^r\), \(r\ge1\), the exact factorization

\[
1+x+\cdots+x^{2^r-1}=\prod_{h=0}^{r-1}(1+x^{2^h})
\tag{14}
\]

gives \(f_{2^r}(u)=\sum_{h=0}^{r-1}f(2^h u)\), for every real
\(u\). Every term in the full pressure is linear in this function.
Therefore

\[
G_{n,2^r}(J,t)=\sum_{h=0}^{r-1}G_n(J,2^h t)>0
\tag{15}
\]

for every \(n\ge1\), real \(J\), and \(t>0\). This is an
analytic digit-count corollary. It does not assert that each such
digit count supplies a separate prime-state tensor construction, and
it makes no claim about odd digit counts.

## 7. Replay and exact scope of the certificate

```sh
python tools/check_global_charged_gaussian_floor.py --output /tmp/global-gaussian-replay.json
```

The [checker](../tools/check_global_charged_gaussian_floor.py) regenerates
all 250 positive scalar intervals (499 tree nodes) using only integer
outward arithmetic. It also replays 65,151 exact folding-range
controls, 1,275 exact edge-reduction controls, and ten exact polynomial
digit-factorization identities. Those finite algebra
controls check the transcription; the inequalities above supply the
uniform proof. The complete earlier folding, small-base, and exterior
certificates are explicitly linked in the evidence. No unrecorded
charge grid, floating optimization, or limiting numerical experiment
is used to establish positivity.
