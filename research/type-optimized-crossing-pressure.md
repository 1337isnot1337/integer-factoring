# Exact type optimization for crossing pressure

## Result and scope

A crossing certificate can optimize the common external letter frequencies,
as well as the physical degree gauges. The optimization over frequencies has
an exact finite form. Let \(B_{\rm raw}(s)\) be the product of the two
gauged Smith profile sums, or of two valid module-rate upper bounds. Let
\(G_a\) be the sum of the physical heights for letter \(a\), across all
four legs and both cuts. For positive integers \(A,B\), put \(R=A/B\) and choose
\(s=t^A\) with rational \(t>1\). Then

\[
\boxed{\quad
B_{\rm raw}(t^A)^B
 <\left(\sum_a t^{-BG_a}\right)^A
\quad}                                                       \tag{1}
\]

certifies a rational common letter type with width parameter
\(\rho<R/4\). All quantities in (1) are rational. The optimizing type
is explicit and rational too. In particular, a certificate with
\(2\le R<2\sqrt2\) would improve the factoring exponent through the
existing transfer theorem.

No such post-half certificate is supplied here. An exact gauge for the
[four-label Schur Gaussian seed](schur-gaussian-digit-reduction.md) gives
\(\rho<61/80\), improving that seed's ungauged sufficient bound but
remaining above \(1/\sqrt2\). The strongest complete factoring exponent
is still one half.

This extends the uniform-type physical-gauge theorem in
[PR #30 at `28105165`](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/balanced-physical-gauge-transfer.md).
The new step is arbitrary rational type selection and its exact entropy
optimization. The proof and checks used substantial OpenAI assistant support.

## 1. A common rational type preserves the crossing

Let the local alphabet have size \(q\), and choose positive rational
frequencies \(\pi_0,\ldots,\pi_{q-1}\) with sum one. For each admissible
length \(m\), restrict every external word to have exactly \(m\pi_a\)
copies of each letter \(a\). Use the same type on all four legs.
If the unrestricted polynomial coefficients sum to the crossing, their
restrictions sum to the crossing on

\[
D_m=\frac{m!}{\prod_a(m\pi_a)!}
\]

labels. Write \(H(\pi)=-\sum_a\pi_a\log\pi_a\), with natural logarithms.
The elementary type bounds are

\[
\frac{e^{mH(\pi)}}{(m+1)^q}\le D_m\le e^{mH(\pi)}.           \tag{2}
\]

For completeness, under the product distribution \(\pi\), the probability
of this type is \(D_m e^{-mH(\pi)}\). It is a most likely count vector:
relative to counts \(k_a=m\pi_a\), increasing a count multiplies its
factor by successive ratios \(k_a/(k_a+j)\), and decreasing it gives
ratios \((k_a-j)/k_a\), all at most one. There are at most \((m+1)^q\)
count vectors. Its probability is therefore between \((m+1)^{-q}\) and
one, proving (2).

## 2. Gauge convention and the pressure at a fixed type

The sign convention at infinity matters. Choose four integer height
functions independently in each cut. The two local weighted tensors are

\[
\begin{aligned}
\widetilde K_\infty(\tau;\ell_1,\ldots,\ell_4)
 &=\tau^{\sum_j h_{\infty,j}(\ell_j)}K(\tau^{-1};\ell_1,\ldots,\ell_4),\\
\widetilde K_0(\tau;\ell_1,\ldots,\ell_4)
 &=\tau^{\sum_j h_{0,j}(\ell_j)}K(\tau;\ell_1,\ldots,\ell_4).
\end{aligned}                                                 \tag{3}
\]

Auxiliary rational gauges may also be included with the usual fixed
boundary constants. Both auxiliary indices must stay on one side of each
cut, as in the existing transfer.

Let \(\kappa_{\infty,i}\) and \(\kappa_{0,i}\) be the restored local
Smith orders, with negative values allowed. For \(s>1\), put

\[
U(s)=\sum_i s^{-\kappa_{\infty,i}},\quad
V(s)=\sum_i s^{-\kappa_{0,i}},\quad
B_{\rm raw}(s)=U(s)V(s),
\qquad
G_a=\sum_{j=1}^4\bigl(h_{\infty,j}(a)+h_{0,j}(a)\bigr).       \tag{4}
\]

Thus \(G\) is an unsigned sum in the convention (3). Heights applied to
the original variable before substituting \(z=\tau^{-1}\) would have the
opposite infinity sign.

The same coefficient argument as in the predecessor gives

\[
\begin{aligned}
r(A_{m,k}^{\pi})&\le C_\infty(m+1)
 s^{m\sum_{a,j}\pi_a h_{\infty,j}(a)-k}U(s)^m,\\
c(A_{m,k}^{\pi})&\le C_0(m+1)
 s^{m\sum_{a,j}\pi_a h_{0,j}(a)+k}V(s)^m.
\end{aligned}                                                 \tag{5}
\]

The original coefficient index cancels in the product. Consequently

\[
\sum_k\sqrt{r(A_{m,k}^{\pi})c(A_{m,k}^{\pi})}
 \le C(m+1)^2 B_\pi(s)^{m/2},
\qquad B_\pi(s)=B_{\rm raw}(s)s^{\sum_a\pi_aG_a}.              \tag{6}
\]

Valid module-rate bounds can replace \(U,V\): their coefficient estimates
have the same shifts in (5). The [minimal-boundary completeness theorem](minimal-boundary-module-completeness.md)
can be used separately within either chosen factorization.

Equations (2) and (6) give

\[
\rho(\pi)=\frac{\log B_{\rm raw}(s)+(\log s)\sum_a\pi_aG_a}
                 {4H(\pi)}.                                \tag{7}
\]

For fixed rational type, admissible lengths are multiples of a fixed
denominator. Ranking and unranking the words use exact multinomial counts.
Choosing the least admissible length with \(D_m\ge2^\lambda\) gives
\(mH(\pi)=\lambda\log2+O(\log(\lambda+2))\). The complete coefficient
tables and their rank factorizations retain the existing \(2^{O(\lambda)}\)
construction cost. The constants may depend on the fixed certificate.

## 3. Eliminate the type by a finite identity

Fix a positive real \(R\), and define

\[
Z_R(s)=\sum_a s^{-G_a/R},\qquad
\pi_a^*=\frac{s^{-G_a/R}}{Z_R(s)}.                            \tag{8}
\]

For any type \(\pi\), direct substitution gives

\[
\log B_{\rm raw}+(\log s)\sum_a\pi_aG_a-RH(\pi)
=\log B_{\rm raw}-R\log Z_R
 +R\sum_a\pi_a\log\frac{\pi_a}{\pi_a^*}.                    \tag{9}
\]

The last sum is nonnegative. One elementary proof applies
\(-\log u\ge1-u\) to \(u=\pi_a^*/\pi_a\), then sums; zero
frequencies follow by continuity. Equality holds at \(\pi=\pi^*\).
Hence some type satisfies \(\rho(\pi)<R/4\) exactly when

\[
B_{\rm raw}(s)<Z_R(s)^R.                                    \tag{10}
\]

For real parameters, a strict inequality remains true for sufficiently
close rational types, by continuity. There is no numerical optimization
over a simplex left in (10). In the
rational setting of the introduction, \(s=t^A\) and \(R=A/B\), so

\[
Z_R=\sum_a t^{-BG_a},\qquad
\pi_a^*=\frac{t^{-BG_a}}{\sum_b t^{-BG_b}}.                   \tag{11}
\]

Every type coordinate is rational. Raising (10) to the positive integer
power \(B\) gives precisely (1). If (1) holds and \(R<2\sqrt2\), then
the [effective crossing transfer](pieces-to-factoring.md) gives score at
least \(\log_2(4/R)>1/2\). The present evidence does not satisfy such
a post-half inequality.

## 4. Exact gauge improvement for the four-label Gaussian seed

The seed uses labels empty, plus, minus, and both. At infinity take
the cut \((x,u)\mid(y,v)\), with plus/minus heights at the four legs

\[
(1,1),\ (1,1),\ (0,0),\ (2,2).
\]

At zero take the other cut, with heights

\[
(0,0),\ (0,0),\ (1,1),\ (-1,-1).
\]

The height of the both label is the sum of its two species heights;
the empty label has height zero. In the convention (3), both restored
order lists are

\[
(0,1,1,2,3,3,4,4,4,4,5,5,6,7,7,8).                        \tag{12}
\]

Thus

\[
B_{\rm raw}(s)=(1+s^{-1})^4(1+s^{-3})^4,
\qquad G=(0,4,4,8).                                        \tag{13}
\]

Take independent plus/minus density \(p=2/9\). The resulting rational
four-letter type is

\[
\pi=(49,14,14,4)/81,
\quad H(\pi)=2H(2/9),\quad \sum_a\pi_aG_a=16/9.
\]

At \(s=5/2\), define the rational numbers

\[
J=(5/2)^{16}(1+2/5)^{36}(1+8/125)^{36}=B_\pi(5/2)^9,
\qquad Q=\frac{9^9}{2^2 7^7}=e^{9H(2/9)}.
\]

Exact integer comparison gives \(J^{10}<Q^{61}\). Therefore

\[
\rho(\pi)=\frac{\log B_\pi}{8H(2/9)}<\frac{61}{80}.          \tag{14}
\]

This is an all-length sufficient certificate for that fixed rational
type. It is weaker than the known half-rate construction. A floating
search near \(p=0.22551,s=2.55489\) gives about \(0.76023\), but (14)
uses only the displayed rational data and does not assert a global optimum.

## 5. Reproduction and controls

Run:

```sh
python3 tools/check_type_optimized_pressure.py \
  --output evidence/type-optimized-pressure.json
```

The checker reconstructs (12) from the complete seed table and independently
from block-Toeplitz ranks. It checks (14) exactly, evaluates a rational
instance of (1), and records a failing post-half control for the same seed.
It also verifies invariance under a uniform height shift: forgetting the
mean term would falsely improve the certificate. These are finite controls
for the separate uniform proof above.

Files: [checker](../tools/check_type_optimized_pressure.py) and
[exact evidence](../evidence/type-optimized-pressure.json).
