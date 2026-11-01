# Asymmetric supporting lines for crossing pieces

## Result and scope

A complete crossing family can be tested against an asymmetric linear
bound on its two factor widths. This extends
[independent paired-type pressure](paired-type-crossing-pressure.md).
It uses the same coefficient construction and retains both indices when
the rectangular alphabet is paired with its swapped copy.

Let the normalized widths of a square crossing piece be

\[
 r=\frac{\log \operatorname{rank}_R}{2\log D},\qquad
 c=\frac{\log \operatorname{rank}_C}{2\log D},
\]

where every physical alphabet has size \(D\). For fixed positive
\(a,b\) with \(a^2+b^2=1\), the criterion below gives the uniform bound

\[
 ar+bc\le R+o(1).                                             \tag{1}
\]

If \(R<1\), a fixed rectangular aspect ratio gives a strict improvement
over the leading exponent \(1/2\) in the
[rectangular grid recurrence](rectangular-grid-frontier.md). This is a
geometric transfer statement. No currently certified seed in this note
has \(R<1\). The companion
[complete rectangular transfer](asymmetric-modular-transfer.md) proves
the uniform evaluator and factoring costs under an effective all-width
crossing hypothesis.

## 1. Two complementary coefficient bounds

Take a polynomial tensor \(K(w;x,y,u,v)\) over a fixed finite field, with
alphabets \(\mathcal X,\mathcal Y,\mathcal Y,\mathcal X\), whose value at
\(w=1\) is \([x=v][y=u]\). The two cuts are
\((x,y)\mid(u,v)\) and \((x,u)\mid(y,v)\). Choose one cut at infinity
and the other at zero. The names \(R,C\) can be interchanged throughout.

Use physical gauges \(h_{e,j}\) at \(e\in\{\infty,0\}\), with the raw
infinity convention

\[
 \widetilde K_\infty(\tau)=\tau^{\sum_jh_{\infty,j}}K(\tau^{-1}),
 \qquad
 \widetilde K_0(\tau)=\tau^{\sum_jh_{0,j}}K(\tau).
\]

Let their complete Smith orders be \(\kappa_{\infty,i}\) and
\(\kappa_{0,i}\). These are raw Laurent orders; an entry of degree \(d\)
has infinity valuation \(-d\), without an additional homogenizing
degree. Rational heights and orders can be made integral by a fixed
ramification. Negative orders are allowed.

Restrict the \(x,v\) words to a rational type \(\pi_X\), and the
\(y,u\) words to a rational type \(\pi_Y\). Write

\[
 g_e=\sum_x\pi_X(x)(h_{e,x}(x)+h_{e,v}(x))
     +\sum_y\pi_Y(y)(h_{e,y}(y)+h_{e,u}(y)).
\]

Let \(A_{m,k}\) be coefficient \(w^k\) in this restricted \(m\)-fold
tensor power. Choose one intrinsic scale \(t>0\). The Smith-module
coefficient bound, used at scales \(t/a\) and \(t/b\), gives

\[
\begin{aligned}
 a\log\operatorname{rank}_R A_{m,k}
 &\le ma\log\sum_i e^{-t\kappa_{\infty,i}/a}
       +t(mg_\infty-k)+O(\log(m+1)),\\
 b\log\operatorname{rank}_C A_{m,k}
 &\le mb\log\sum_i e^{-t\kappa_{0,i}/b}
       +t(mg_0+k)+O(\log(m+1)).
\end{aligned}                                                  \tag{2}
\]

Zero pieces can be omitted. The constants depend on the fixed path,
gauges and positive parameters, and are uniform over all coefficient
indices. One proof expands the two invertible local Smith matrices as
power series. A channel with total order at most the requested
coefficient index contributes only \(O(m+1)\) rank-one terms. The
exponential bound on these channels is the displayed profile sum raised
to the \(m\)-th power. Type projections cannot increase matrix rank.

Adding (2) cancels the coefficient index exactly. Put

\[
\begin{aligned}
 \Phi(a,b)&=a\log\sum_i e^{-t\kappa_{\infty,i}/a}
           +b\log\sum_i e^{-t\kappa_{0,i}/b},\\
 G_X(x)&=t\sum_e(h_{e,x}(x)+h_{e,v}(x)),\\
 G_Y(y)&=t\sum_e(h_{e,y}(y)+h_{e,u}(y)),\\
 F(\pi_X,\pi_Y)&=\Phi(a,b)+\langle\pi_X,G_X\rangle
                              +\langle\pi_Y,G_Y\rangle.
\end{aligned}                                                  \tag{3}
\]

Thus \(a\log\operatorname{rank}_R A_{m,k}
+b\log\operatorname{rank}_C A_{m,k}\le mF+O(\log(m+1))\).
The same intrinsic \(t\) is essential; independently choosing the two
Chernoff scales would leave a coefficient-dependent term.

## 2. Equal alphabets and the exact normalization

The pair swap \(K^{\rm sw}(x',y',u',v')=K(u',v',x',y')\) preserves
both cut ranks. It transposes the first flattening and permutes rows and
columns of the second. Pair the two restricted paths and retain the
separate pieces \(A_{m,k}\otimes A^{\rm sw}_{m,l}\).

Their sum is a square crossing with physical dimension
\(D_m=D_X(m)D_Y(m)\). Writing
\(H=H(\pi_X)+H(\pi_Y)>0\), multinomial bounds give
\(\log D_m=mH+O(\log(m+1))\). Both ranks multiply, so each paired piece
satisfies

\[
 ar+bc\le \frac{F(\pi_X,\pi_Y)}{H}+O\!\left(\frac{\log(m+1)}m\right).
                                                               \tag{4}
\]

The numerator of the rank bound doubles, while the normalization is
\(2\log D_m\). These factors cancel. There are polynomially many
paired coefficient pieces. Keeping the two indices avoids any
unsupported bound for a merged coefficient.

For \(a=b=1/\sqrt2\), divide all parameters of the usual paired-type
pressure, including the intrinsic scale, by \(\sqrt2\). Formula (4)
then gives \(R=\sqrt2\rho\), with \(\rho\) in the convention of the
paired-type note. A constant crossing has \(r=c=1\), and therefore
has optimal supporting-line value \(a+b\).

## 3. Gibbs optimization and rational certificates

For a proposed \(R>0\), define

\[
 Z_X(R)=\sum_xe^{-G_X(x)/R},\qquad
 Z_Y(R)=\sum_ye^{-G_Y(y)/R},
 \qquad \pi_X^*(x)=\frac{e^{-G_X(x)/R}}{Z_X(R)},
\]

and similarly for \(\pi_Y^*\). The exact entropy identity is

\[
 F(\pi_X,\pi_Y)-RH
 =\Phi-R\log(Z_XZ_Y)
  +R\bigl(D(\pi_X\Vert\pi_X^*)+D(\pi_Y\Vert\pi_Y^*)\bigr).
                                                               \tag{5}
\]

Consequently some types satisfy \(F/H<R\) if and only if

\[
 \Phi(a,b)<R\log(Z_X(R)Z_Y(R)).                              \tag{6}
\]

Strict real parameters can be replaced by nearby rational parameters
and types. Positive rational points on the unit circle are dense in its
positive arc. The numerical
[support-line implementation](../tools/asymmetric_crossing_pressure.py)
finds the Gibbs root; it is not an exact certificate.

Here is a specialization using rational arithmetic only. Let
\(a=u/c\), \(b=v/c\), where positive integers satisfy
\(u^2+v^2=c^2\). Write the target \(R=e/f\), choose rational \(z>1\),
and set \(t=uve\log z\). In the following display let
\(g_X=G_X/t\), \(g_Y=G_Y/t\) be the **unscaled** total physical heights,
and put
\(U(s)=\sum_i s^{-\kappa_{\infty,i}}\),
\(V(s)=\sum_i s^{-\kappa_{0,i}}\). Then (6) is exactly

\[
 \boxed{
 U(z^{vce})^{uf}\,V(z^{uce})^{vf}
 <\left[
    \left(\sum_xz^{-uvf g_X(x)}\right)
    \left(\sum_yz^{-uvf g_Y(y)}\right)
   \right]^{ce}.}
                                                               \tag{7}
\]

For integral heights and orders all terms are rational. If their common
denominator is \(L\), apply (7) to the integral lists
\(L\kappa_\infty,L\kappa_0,Lg_X,Lg_Y\). This corresponds to choosing
\(t=Luve\log z\) for the original lists, or to the fixed ramified path
at scale \(uve\log z\). Every power can then be checked by integers.
A global monomial
\(w^d\) subtracts \(d\) from every infinity order and adds \(d\) to
every zero order. Its two factors in (7) cancel exactly.

## 4. The geometric recurrence

Suppose every crossing piece satisfies \(ar+bc\le R\). Fix
\(0<\theta<1\) and propose exponent \(0<\alpha<1\). A rectangular
split with \(\log p/\log K=\theta\) has normalized orientation costs
\(r/\theta^\alpha\) and \(c/(1-\theta)^\alpha\). Then

\[
 \min\left(\frac r{\theta^\alpha},
           \frac c{(1-\theta)^\alpha}\right)
 \le \frac{R}{a\theta^\alpha+b(1-\theta)^\alpha}.              \tag{8}
\]

Indeed, if both quotients exceeded the right side, their positive
weighted sum would contradict \(ar+bc\le R\).

At \(\alpha=1/2\), choose \(\theta=a^2\). The denominator is exactly
\(a^2+b^2=1\). Thus **every strict gap \(R<1\) permits some
\(\alpha>1/2\)** by continuity. No second support inequality is needed
for this leading, fixed-aspect recurrence. A linear bound also holds for
averages of profiles at distinct rectangles.

For a specified \(\alpha\), the largest possible denominator is

\[
 D_\alpha(a,b)=
 \left(a^{1/(1-\alpha)}+b^{1/(1-\alpha)}\right)^{1-\alpha},
 \quad
 \theta=\frac{a^{1/(1-\alpha)}}
 {a^{1/(1-\alpha)}+b^{1/(1-\alpha)}}.                        \tag{9}
\]

This follows by differentiating the strictly concave two-term function,
or by the two-term Hölder inequality. For balanced weights,
\(D_\alpha=2^{1/2-\alpha}\), so \(R=\sqrt2\rho<D_\alpha\) is exactly
the familiar condition \(\rho<2^{-\alpha}\).

For example, consider the **conditional geometric** line with
\(a=3/5\), \(b=4/5\), \(R=9/10\). At
\(\theta=9/25\), \(\alpha=3/5\), exact fifth-power comparisons give

\[
 (9/25)^{3/5}>27/50,\qquad (16/25)^{3/5}>19/25.
\]

The denominator in (8) therefore exceeds \(233/250\), and the
orientation cost is less than \(225/233<1\). The checker also gives
rational intervals for the smaller gap \(R=99/100\),
\(\alpha=51/100\), at the same aspect ratio. These examples verify the
geometric implication; they do not assert that such a crossing seed has
been constructed.

## 5. Connect the geometry to a full evaluator

The strict geometric margin has to dominate the finite grid terms.
For a \(p\) by \(d\) grid, the exact normalized volumes are

\[
 \frac{U_R}{N}=\frac{1-1/p}{d}
       +(1-1/p)(1-1/d)\,\bar r,
 \quad
 \frac{U_C}{N}=\frac{1-1/d}{p}
       +(1-1/p)(1-1/d)\,\bar c.
\]

Both straight-edge terms tend to zero when the chosen fixed aspect has
both sides growing. A complete uniform recursion must additionally
prove integer grouping and destination bounds, all-width type padding,
width ceilings, piece and enumeration counts, effective rank factors,
their bit costs, and termination with a controlled small-parameter base
case. The [complete transfer proof](asymmetric-modular-transfer.md)
discharges the grouping, padding, enumeration, circuit, and factoring
costs under an effective all-width crossing hypothesis. Equation (7)
alone does not supply that hypothesis. A numerical supporting-line
value below one still needs an exact path/profile certificate, common
types, all-width padding, and effective factors before it can change
the factoring exponent.

## 6. Exact controls

```sh
python3 tools/check_asymmetric_support_line.py \
  --output evidence/asymmetric-support-line.json
```

The checker verifies the rational Gibbs identity, the constant rectangular
crossing normalization, cancellation of a common monomial, the exact
half-exponent equalizer, and rational-power intervals in Section 4.
It records no new crossing or factoring bound.

The [numerical oracle controls](../tools/check_asymmetric_pressure_oracle.py)
separately compare balanced normalization against the existing paired-type
oracle, check every gradient coordinate at four generic points, and test
raw-infinity monomial cancellation in both cuts. Their
[record](../evidence/asymmetric-pressure-oracle.json) validates the search
interface and is not an exact pressure certificate.
