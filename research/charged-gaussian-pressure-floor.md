# Positive critical Gaussian pressure at nonzero aggregate charge

## Results and scope

The binary interval Gaussian has strictly positive critical pressure on
a continuous charge interval whose width grows linearly with its length.
This extends the [natural-gauge floor](truncated-symmetric-interval-frontier.md)
to charges that can occur in the occupation-filter searches.

**Theorem 1.** Put

\[
C_n=(6-4\sqrt2)n+2-2\sqrt2.
\tag{1}
\]

For every integer \(n\ge3\), every \(t>0\), and every real
\(|J|\le C_n\), the binary Gaussian pressure at
\(\rho=1/\sqrt2\) satisfies \(G_n(J,t)>0\).

The analytic statement holds for all integer lengths, including odd
lengths. A prime-power length remains necessary for the particular
all-indices filter construction. The bound is sufficient; it is not
claimed to be the largest positive charge interval.

**Theorem 2.** For every digit count \(q\ge2\), integer \(n\ge1\),
and even integer charge \(J=2k\) with \(|k|\le n\), the same
critical Gaussian pressure is strictly positive at every scale.

The proof of Theorem 1 combines analytic estimates with a finite,
outward-rational certificate. It does not assume that zero charge
minimizes pressure. In fact that minimization statement is false: the
complete finite-pressure certificate below gives
\(G_{32}(323/5,3)<G_{32}(0,3)\), with both values positive.
No conclusion is asserted here for every real charge outside (1), and
no improved factoring exponent follows from these positive floors.

## 1. Normalization and an exact even-charge identity

Write

\[
R=\sqrt2,\qquad
f_q(x)=\log\!\left(\sum_{d=0}^{q-1}e^{-dx}\right).
\]

The critical pressure in
[`gaussian_pressure`](../tools/prime_state_filter_pressure.py) is

\[
\begin{aligned}
G_{n,q}(J,t)={}&
\sum_{i=0}^{n-1}\left[f_q((2i+1-J)t)+f_q((2i+1+J)t)\right]\\
&+\sum_{a=0}^{2n-1}f_q((2a+1-2n)t)\\
&-R\sum_{i=0}^{n-1}\left[
 f_q((4i-2n+2-J)t/R)+f_q((4i-2n+2+J)t/R)
\right].
\end{aligned}
\tag{2}
\]

For the binary case the subscript \(q=2\) is omitted. This
normalization includes both cuts and the complete paired physical
partition. The scale in a ramified construction is \(t=P\log s\).

Put

\[
S_{m,q}(t)=\sum_{i=0}^{m-1}f_q((2i+1)t),\qquad S_{0,q}=G_{0,q}=0.
\]

At \(J=2k\), \(|k|\le n\), the physical energy multiset in (2)
is the union of the natural physical energy multisets for lengths
\(n+k\) and \(n-k\). The reciprocal identity
\(f_q(-x)=f_q(x)+(q-1)x\) gives

\[
\boxed{
G_{n,q}(2k,t)=\frac{G_{n+k,q}(0,t)+G_{n-k,q}(0,t)}2
 +2S_{n,q}(t)-S_{n+k,q}(t)-S_{n-k,q}(t).
}
\tag{3}
\]

For completeness, the ground contribution from the two cuts is
\((q-1)(n^2+k^2)t\). Their remaining thermal energy multiset consists
of the positive odd energies for lengths \(n+k,n-k,n,n\). These
facts prove (3) without a limiting argument. Since the increments of
\(S_{m,q}\) strictly decrease with \(m\), its last three terms are
nonnegative. At least one length on the right of (3) is positive.
The natural all-digit floor therefore proves Theorem 2. The checker
replays all three multiset identities for every \(1\le n\le64\).

## 2. A positive periodic infinite continuation

From now on let \(f(x)=\log(1+e^{-x})\). Define the convergent series

\[
A(J,t)=\sum_{j=1}^{\infty}
 \left[f((2j-1+J)t)+f((2j-1-J)t)\right]
\]

and

\[
H(J,t)=A(J,t)+A(0,t)-2R A(J/2,Rt).
\tag{4}
\]

The reciprocal identity gives
\(A(J+2,t)=A(J,t)+(J+1)t\). Thus \(H\) is even and
four-periodic in \(J\). It is enough to consider \(0\le J\le2\).
We prove

\[
\boxed{H(J,t)>0\quad(J\in\mathbb R,\ t>0).}
\tag{5}
\]

### Small scales: \(0<t\le1\)

Set \(Q=e^{-\pi^2/t}\), \(V=e^{-\pi^2/(Rt)}\), and
\(\mathcal P(y)=\sum_{j\ge1}\log(1-y^j)\). The Jacobi product,
theta transformation, and eta transformation give

\[
A(J,t)=\frac{\pi^2}{12t}-\frac t{12}+\frac{tJ^2}4
 -\mathcal P(Q^2)+\log\theta_3(\pi J/2,Q).
\tag{6}
\]

The primary identities are
[DLMF 20.5.3](https://dlmf.nist.gov/20.5.E3),
[DLMF 20.7.32](https://dlmf.nist.gov/20.7.E32), and
[DLMF 23.18.5](https://dlmf.nist.gov/23.18.E5).
For example the product first gives
\(A=\log\theta_3(iJt/2,e^{-t})-\mathcal P(e^{-2t})\);
applying the two transformations yields (6).
Substitution cancels the inverse-scale and charge-quadratic terms:

\[
\begin{aligned}
H(J,t)={}&\frac t6-2\mathcal P(Q^2)+2R\mathcal P(V^2)\\
&+\log\theta_3(\pi J/2,Q)+\log\theta_3(0,Q)
-2R\log\theta_3(\pi J/4,V).
\end{aligned}
\tag{7}
\]

Let \(u=e^{-3/t}<1/20\). Since \(\pi>3\) and \(R<3/2\),
\(Q<u^3\) and \(V<u^2\). The elementary bounds

\[
1-\frac{2y}{1-y^3}\le\theta_3(v,y)
\le1+\frac{2y}{1-y^3},\qquad
\mathcal P(y)\ge-\frac{y}{(1-y)^2}
\]

hold in the ranges used here. Also \(\theta_3(0,Q)>1\).
Use \(\log(1-x)\ge-x/(1-x)\) and
\(\log(1+x)\le x\), and discard the positive term
\(-2\mathcal P(Q^2)\). The three remaining negative errors in (7)
are bounded above by \(3u^3,4u^4,7u^2\), respectively. Their sum
is less than \(8u^2\), so

\[
H(J,t)>\frac t6-8e^{-6/t}>\frac{11t}{75}>0.
\tag{8}
\]

The second inequality uses that \(e^{-6/t}/t\) increases on
\((0,1]\), and \(e^3>20\). The latter follows already from the
Taylor sum through degree eight and is checked with exact fractions.

### Large scales: \(t\ge6\)

For \(0\le J\le2\),

\[
A(J/2,Rt)\le
\frac{2e^{-R(1-J/2)t}}{1-e^{-2Rt}}.
\tag{9}
\]

First suppose \(0\le J\le1\). Concavity of \(\log(1+y)\) on
\([0,1]\) gives
\(A(J,t)\ge(\log2)e^{-(1-J)t}\). The ratio of the upper bound
for \(2R A(J/2,Rt)\) to this lower bound is at most

\[
\frac{4R}{\log2}
\frac{e^{-[R-1+(1-R/2)J]t}}{1-e^{-2Rt}}
<\frac{9e^{-2t/5}}{1-e^{-14t/5}}<1.
\tag{10}
\]

Here \(7/5<R<3/2\) and \(\log2>2/3\). At \(t=6\), the
numerator in the last ratio is below \(9/11\), while its denominator
is above \(10/11\); the ratio decreases thereafter.
The exact Taylor sum through degree eleven proves \(e^{12/5}>11\).

For \(J=1+v\), \(0\le v\le1\), we instead have
\(A(J,t)\ge\max(vt,\log2)\). The minimum over \(v\in[0,1]\) of

\[
\max(vt,\log2)e^{(1-v)t/R}
\]

is \(\min\{t,(\log2)e^{(t-\log2)/R}\}=t\) for \(t\ge6\).
On the first branch the function decreases; on the second its
logarithm is concave, so only the endpoints are needed. The second
displayed endpoint is at least \(t\): its ratio to \(t\) increases
for \(t\ge6\), and at six this follows from
\(\log2>2/3\), \(\log2<1\), \(R<3/2\), and \(e^3>20\).
Consequently

\[
A(J,t)\ge t e^{-(1-v)t/R},\qquad
2R A(J/2,Rt)\le
\frac{4R}{1-e^{-2Rt}}e^{-(1-v)t/R}.
\]

Since \(R<10/7\) and \(e^{2Rt}>21\), the last coefficient is
strictly below six. This proves (5) throughout \(t\ge6\), also
using the positive term \(A(0,t)\).

### The compact interval: an integer certificate

The remaining rectangle \([1,6]\times[0,2]\) is covered by 290
closed rectangles with strictly positive outward-rational lower
bounds for (4). The reproducible implementation is
[`charged_gaussian_floor.py`](../tools/charged_gaussian_floor.py);
the complete compact record is
[`charged-gaussian-floor.json`](../evidence/charged-gaussian-floor.json).

Every interval endpoint is an integer divided by \(2^{80}\), and
each arithmetic operation rounds outward using integer floor and
ceiling. The interval for \(\sqrt2\) is obtained with an integer
square root. To bound \(e^{-x}\), the code first reduces the exponent
to at most \(1/2\), uses the positive exponential series through
degree 24 with a remainder below twice the first omitted term, squares
back, and takes the reciprocal. To bound \(\log(1+y)\), it uses
32 terms of the atanh series with \(z=y/(2+y)\), and the remainder

\[
\frac{2z^{65}}{65(1-z^2)}.
\]

Twenty terms of each paired series in (4) are included. Positive
omitted tails are discarded. The negative omitted tail is bounded by

\[
\frac{4R e^{-Rt(41-J/2)}}{1-e^{-2Rt}}.
\tag{11}
\]

The checker reconstructs the entire rectangle partition from binary
subdivision paths. It verifies coverage, every exact lower numerator,
and the record hash. The smallest certified lower bound is

\[
\frac{53861972282886535215}{1208925819614629174706176}>0.
\tag{12}
\]

No floating-point sign is used in this proof.

## 3. Restoring the finite interval, with both parities

Let \(\epsilon=n\bmod2\). Repeated use of
\(f(-x)=f(x)+x\) gives the exact identity

\[
G_n(J,t)=H(J+2\epsilon,t)+T_{\rm phys}-T_{\rm cuts},
\tag{13}
\]

where

\[
\begin{aligned}
T_{\rm phys}
&=2R\sum_{\ell\ge0}\sum_{\sigma=\pm1}
 f\!\left([R(n+1+2\ell)+\sigma J/R]t\right),\\
T_{\rm cuts}
&=\sum_{\ell\ge0}\left[
 f((2n+1+2\ell+J)t)+f((2n+1+2\ell-J)t)
 +2f((2n+1+2\ell)t)\right].
\end{aligned}
\tag{14}
\]

Here is the parity calculation explicitly. The infinity-cut sum is
\(n^2t+A(0,t)\) minus its indicated tail. If \(n\) is even, the
physical sum after multiplying by \(R\) is
\(n^2t+2R A(J/2,Rt)-T_{\rm phys}\).
If \(n\) is odd, the physical energies include zero and paired
multiples of four; their sum is
\((n^2-1-J)t+2R A(1+J/2,Rt)-T_{\rm phys}\).
Finally \(A(J+2,t)=A(J,t)+(J+1)t\) gives (13) in the odd case.
These identities do not require any of the displayed shifted energies
to be positive.

Convex midpoint bounds for the four cut tails, and decreasing rectangle
bounds for the two physical tails, imply

\[
T_{\rm cuts}\le\frac2t\int_{(2n-|J|)t}^{\infty} f(x)\,dx,
\qquad
T_{\rm phys}\ge\frac2t
\int_{[R(n+1)+|J|/R]t}^{\infty} f(x)\,dx.
\tag{15}
\]

The second integral starts no later than the first exactly when

\[
R(n+1)+|J|/R\le2n-|J|,
\]

which is \(|J|\le C_n\). In this range (13) is at least
\(H(J+2\epsilon,t)>0\). The range is nonempty for \(n\ge3\).
This proves Theorem 1 for both parities.

## 4. What this closes for homogeneous occupation filters

The [homogeneous Boolean-filter theorem](homogeneous-boolean-filter-bounds.md)
gives two complementary estimates. Above its explicit threshold
\(t_n^*\), every such filter has larger cost than the Gaussian at
the same real charge. Below \(t_n^*\), for fixed Boolean degree
\(d\) and bounded \(|J|\le J_0\), its full critical pressure is
positive for sufficiently large binary prime-power \(n\), uniformly
over the other stated filter parameters.

For fixed \(J_0\), (1) eventually contains \([-J_0,J_0]\).
Theorem 1 therefore also makes the large-scale branch absolutely
positive. Thus, for every fixed \(d,J_0\), the homogeneous family
has strictly positive full critical pressure at **every** scale for
all sufficiently large binary prime-power lengths. This conclusion
includes arbitrary selected subsets, odd ramification, and the full
polynomial-degree costs under that theorem's hypotheses. It does not
cover arbitrary multilevel polynomial filters.

## 5. Replay and a genuine limitation of zero-charge optimization

Run the complete certificate with

```sh
python tools/check_charged_gaussian_floor.py
```

The command verifies 40 exact arithmetic controls, 4,224 complete
even-charge multiset identities, 128 physical parity identities with
symbolic real charge, all 290 positive cover rectangles, and a complete
finite-pressure comparison. To rebuild the same record:

```sh
python tools/check_charged_gaussian_floor.py --generate
```

The finite comparison uses \(n=32,t=3,J=323/5\), and gives outward
dyadic intervals with the common denominator \(2^{80}\). Their
numerical locations are

\[
G_{32}(0,3)\simeq0.1141192364,\qquad
G_{32}(323/5,3)\simeq0.0942770143.
\]

The upper endpoint of the difference interval is strictly negative,
while both individual lower endpoints are positive. These strict signs
are certified with the same integer arithmetic, using every finite
factor in (2). This excludes a global zero-charge minimization claim;
it does not give negative critical pressure.
