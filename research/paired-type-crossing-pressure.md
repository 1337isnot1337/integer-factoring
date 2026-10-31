# Independent paired types and rectangular crossing pressure

## Result and scope

The two crossing pairs can use different physical alphabets and different
letter distributions. A rectangular path with alphabets
\(\mathcal X,\mathcal Y,\mathcal Y,\mathcal X\) can be combined with its
pair-swapped copy to produce an ordinary square crossing. The coefficient
pieces retain two separate indices, so their rank costs multiply exactly.

Let \(B(s)\) be the product of valid two-place profile bounds, and let
\(G_X,G_Y\) be the total physical heights for the two crossing pairs.
For \(R=a/b>0\), with positive integers \(a,b\), and rational \(t>1\), the
exact inequality

\[
\boxed{
 B(t^a)^b<
 \left[
   \left(\sum_{x\in\mathcal X}t^{-bG_X(x)}\right)
   \left(\sum_{y\in\mathcal Y}t^{-bG_Y(y)}\right)
 \right]^a
}                                                               \tag{1}
\]

certifies a square crossing family with \(\rho<R/2\). In particular,
\(R<\sqrt2\) would pass the post-half threshold. This note supplies the
transfer and exact optimization formula; it does not assert that the
rectangular seeds currently pass that threshold.

The construction also applies when both alphabets have the same size but
the two crossing pairs use different types. It extends
[type-optimized crossing pressure](type-optimized-crossing-pressure.md),
which uses one common type on all four legs. The
[coordinate-core construction](coordinate-core-crossings.md) gives
examples whose endpoint is valid under independent paired types.

## 1. Endpoint and coefficient hypotheses

Write the local tensor as \(K(w;x,y,u,v)\), where
\(x,v\in\mathcal X\), \(y,u\in\mathcal Y\),
\(|\mathcal X|=A\), and \(|\mathcal Y|=B\). The two cost cuts are

\[
 (x,y)\mid(u,v),\qquad (x,u)\mid(y,v).
\]

Both matrices have size \(AB\) by \(AB\). The desired endpoint is
\(E(x,y,u,v)=[x=v][y=u]\).

It is enough that the following statement holds for every pair of rational
types \(\pi_X,\pi_Y\) and all admissible lengths \(m\): restrict the
\(x,v\) words to type \(\pi_X\), and the \(y,u\) words to type
\(\pi_Y\). The resulting tensor \(K_m^{\pi_X,\pi_Y}(w)\) is a Laurent
polynomial of degree range \(O(m)\), and its value at \(w=1\) is the
crossing on those type classes. A scalar monomial can clear the finite
negative degree range.

Thus the local tensor itself may be rational and may fail to have a
literal crossing endpoint before the type projection. The stated
polynomiality and endpoint hypotheses must be proved for the particular
family. They do not follow merely from the rectangular alphabet sizes.

The dimensions of the two type classes are

\[
 D_X(m)=\frac{m!}{\prod_x(m\pi_X(x))!},\qquad
 D_Y(m)=\frac{m!}{\prod_y(m\pi_Y(y))!}.
\]

With \(H_X=H(\pi_X)\) and \(H_Y=H(\pi_Y)\), the usual multinomial
bounds give

\[
 \log D_X=mH_X+O(\log(m+1)),\qquad
 \log D_Y=mH_Y+O(\log(m+1)).                                \tag{2}
\]

## 2. Gauges and the coefficient cost

Use the infinity convention
\(\widetilde K_\infty(\tau)=\tau^{\sum h_{\infty,j}}K(\tau^{-1})\),
and the zero convention
\(\widetilde K_0(\tau)=\tau^{\sum h_{0,j}}K(\tau)\).
Let the restored Smith orders in complementary cost cuts be
\(\kappa_{\infty,i}\) and \(\kappa_{0,i}\). Negative orders are allowed.
Set

\[
 U(s)=\sum_i s^{-\kappa_{\infty,i}},\qquad
 V(s)=\sum_i s^{-\kappa_{0,i}},\qquad B(s)=U(s)V(s).            \tag{3}
\]

Any valid module bounds with the same coefficient-shift convention can
replace these profile sums. Define

\[
\begin{aligned}
 G_X(x)&=\sum_{e\in\{0,\infty\}}
                 \bigl(h_{e,x}(x)+h_{e,v}(x)\bigr),\\
 G_Y(y)&=\sum_{e\in\{0,\infty\}}
                 \bigl(h_{e,y}(y)+h_{e,u}(y)\bigr).
\end{aligned}                                                \tag{4}
\]

These are unsigned sums in the displayed infinity convention. The same
coefficient argument as in the common-type theorem yields, for the
coefficient pieces \(A_{m,k}\) of the restricted path,

\[
 C_m:=\sum_k\sqrt{r(A_{m,k})c(A_{m,k})}
 \le \operatorname{poly}(m)\,
 \left[B(s)s^{\langle\pi_X,G_X\rangle+
                     \langle\pi_Y,G_Y\rangle}\right]^{m/2}.  \tag{5}
\]

The coefficient degree cancels between the two rank bounds. The gauge
shifts are constant on each relevant type class, which is the only place
where a common type is needed in that argument.

For a tensor represented as \(N(w)/(1+w)^d\), use its actual rational
Smith orders. Relative to the numerator, the zero orders are unchanged,
and every infinity order increases by \(d\). Thus \(B(s)\) is divided
by \(s^d\). Apply this correction once. Multiplication by a common power
of \(w\) shifts the two order lists in opposite directions and leaves
\(B(s)\) unchanged.

## 3. Convert the rectangle to an equal alphabet

Define the pair-swapped tensor

\[
 K^{\rm sw}(x',y',u',v')=K(u',v',x',y'),
\]

whose alphabets are \(\mathcal Y,\mathcal X,\mathcal X,\mathcal Y\).
This permutation transposes the first cost flattening and permutes rows
and columns of the second. Hence it preserves both ranks of every piece.
Its type dimensions are \(D_Y,D_X,D_X,D_Y\).

Tensor the original restricted crossing with the swapped one. Every
physical alphabet then has size

\[
 D_m=D_X(m)D_Y(m).
\]

Keep the coefficient pieces indexed separately:

\[
 A_{m,k}\otimes A^{\rm sw}_{m,l}\qquad\text{for every pair }(k,l).
\]

Their sum is the square crossing. The ranks in either cut multiply under
tensor products, giving the exact cost identity

\[
\begin{aligned}
 \sum_{k,l}\sqrt{
  r(A_{m,k}\otimes A^{\rm sw}_{m,l})
  c(A_{m,k}\otimes A^{\rm sw}_{m,l})}
 &=\left(\sum_k\sqrt{r(A_{m,k})c(A_{m,k})}\right)^2\\
 &=C_m^2.
\end{aligned}                                                \tag{6}
\]

Equivalently, use two independent polynomial variables. There is no need
to merge coefficients having the same sum \(k+l\), and no rank-cost
inequality for such a merge is assumed.

Combining (2), (5), and (6) gives the sufficient width parameter

\[
 \boxed{
 \rho(\pi_X,\pi_Y)=
 \frac{\log B(s)+(\log s)
       (\langle\pi_X,G_X\rangle+\langle\pi_Y,G_Y\rangle)}
      {2(H_X+H_Y)}.
 }                                                          \tag{7}
\]

At least one entropy must be positive. Taking admissible lengths through
multiples of the two type denominators, and then restricting the square
crossing to the required number of labels, retains the usual polynomial
factors and effective construction bounds.

## 4. Exact independent Gibbs optimization

Fix \(R>0\). Define

\[
 Z_X=\sum_x s^{-G_X(x)/R},\quad
 Z_Y=\sum_y s^{-G_Y(y)/R},\qquad
 \pi_X^*(x)=\frac{s^{-G_X(x)/R}}{Z_X},\quad
 \pi_Y^*(y)=\frac{s^{-G_Y(y)/R}}{Z_Y}.                        \tag{8}
\]

The two entropy identities add to

\[
\begin{aligned}
 &\log B+(\log s)
   (\langle\pi_X,G_X\rangle+\langle\pi_Y,G_Y\rangle)
   -R(H_X+H_Y)\\
 &\quad=\log B-R\log(Z_XZ_Y)
  +R\left[D(\pi_X\|\pi_X^*)+D(\pi_Y\|\pi_Y^*)\right].
\end{aligned}                                                \tag{9}
\]

Both relative entropies are nonnegative. Therefore some pair of types
satisfies \(\rho<R/2\) exactly when

\[
 B(s)<(Z_XZ_Y)^R.                                           \tag{10}
\]

For integer heights, \(R=a/b\) and \(s=t^a\) with rational \(t>1\),
both Gibbs types are rational, and (10) is exactly the rational comparison
(1). A strict inequality at real parameters can instead be transferred to
nearby rational types by continuity.

For a square alphabet with \(\pi_X=\pi_Y\), (7) becomes the original
common-type formula. If also \(G_X=G_Y=G/2\), (10) is precisely its Gibbs
criterion with the common-type parameter \(2R=4\rho\). Independent paired
types can improve that constrained optimization; they do not require any
new endpoint hypothesis when the local endpoint is already the crossing.

## 5. Reproduction and limits

The exact checker verifies the pair-swap cut conventions, the tensor-product
rank identity on a rectangular decomposition, the multinomial dimensions,
and the rational Gibbs specialization. It also checks the constant
rectangular crossing control: uniform paired types give exactly \(\rho=1\).

```sh
python3 tools/check_paired_type_crossing_pressure.py \
  --output evidence/paired-type-crossing-pressure.json
```

Floating profile or type optimization is exploratory. A proposed
post-half result still requires exact Smith or rank bounds, exact endpoint
and polynomiality proofs, and a strict certified instance of (1).
