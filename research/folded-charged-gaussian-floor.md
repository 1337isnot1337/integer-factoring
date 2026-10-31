# A wider continuous Gaussian pressure floor

## Result

For every integer \(n\ge3\), every \(t>0\), and every real
\(|J|\le2n-2-2\sqrt2\), the binary Gaussian critical pressure is strictly
positive:

\[
G_n(J,t)>0.
\tag{1}
\]

This follows from a sharper finite-tail comparison with the previously
proved [positive periodic continuation](charged-gaussian-pressure-floor.md).
The note also proves an exact continuous charge-folding inequality.
Combining (1) with the
[exterior floor](exterior-charged-gaussian-floor.md) leaves only the
two bands \(2n-2-2\sqrt2<|J|<2n\) outside these continuous conclusions.
The separate even-charge theorem still applies inside those bands.
No all-real-charge floor is claimed here.

## 1. An exact identity for every digit count

Let \(q\ge2\), \(R=\sqrt2\), and
\(f=f_q=\log(1+e^{-x}+\cdots+e^{-(q-1)x})\). Write
\(S_m(t)=\sum_{i=0}^{m-1}f((2i+1)t)\).
For integers \(0\le k\le n-2\), put

\[
M=n-k,\qquad L=n+k,\qquad J=2k+\delta,
\qquad 0\le\delta\le2.
\tag{2}
\]

For \(0\le\ell<k\), define \(e_\ell=2(M+1+2\ell)\) and

\[
\begin{aligned}
C_\ell={}&\frac12\sum_{s=\pm1}
 \bigl[f((e_\ell+s+\delta)t)-f((e_\ell+s-\delta)t)\bigr]\\
&-R\bigl[f((e_\ell+\delta)t/R)-f((e_\ell-\delta)t/R)\bigr].
\end{aligned}
\tag{3}
\]

The complete pressure satisfies

\[
\boxed{
G_{n,q}(2k+\delta,t)
=\frac{G_{L,q}(\delta,t)+G_{M,q}(\delta,t)}2
 +2S_n(t)-S_L(t)-S_M(t)+\sum_{\ell=0}^{k-1}C_\ell.
}
\tag{4}
\]

To derive it, shift the positive and negative zero-cut odd-energy
lists by \(2k\). After the reciprocal identity, their thermal sums
become
\(\sum_{i<L}f((2i+1+\delta)t)
+\sum_{i<M}f((2i+1-\delta)t)\).
The linear term is \((q-1)(k^2+k\delta)t\).
For the physical energy list, the interval shifted by \(2k\)
consists of the centered length-\(M\) list and the \(k\) extra
positive energies \(e_\ell\). The centered length-\(L\) list
adds both signs of these extra energies. Applying the reciprocal
identity to this comparison cancels the term
\((q-1)k\delta t\) and gives precisely (3)--(4).
The infinity-cut ground contribution completes the remaining
\((q-1)(n^2+k^2)t\) term on both sides.

## 2. Every correction is nonnegative

Since \(f\) strictly decreases, the increments of \(S_m\) decrease
with \(m\). Therefore

\[
2S_n-S_L-S_M
=\sum_{i=0}^{k-1}
 \bigl[f((2(M+i)+1)t)-f((2(n+i)+1)t)\bigr]\ge0.
\tag{5}
\]

For (3), put \(h=-f'\). The function \(h\) is positive and
decreasing because \(f\) is decreasing and convex. Differentiating
under the finite integral gives

\[
C_\ell=t\int_{-\delta}^{\delta}\left[
h((e_\ell+x)t/R)
-\frac{h((e_\ell-1+x)t)+h((e_\ell+1+x)t)}2
\right]dx.
\tag{6}
\]

For every integrand point, \(e_\ell+x\ge2M\ge4\). Hence

\[
\frac{e_\ell+x}{R}\le e_\ell+x-1,
\tag{7}
\]

because \(R/(R-1)=2+R<4\). The first argument of \(h\) in
(6) is no larger than either of the other two. Monotonicity proves
\(C_\ell\ge0\). Thus, for every digit count,

\[
\boxed{
G_{n,q}(2k+\delta,t)
\ge\frac{G_{n+k,q}(\delta,t)+G_{n-k,q}(\delta,t)}2.
}
\tag{8}
\]

This inequality is exact and uniform in the scale. Positivity at the
two smaller-charge instances is its only additional requirement.

## 3. Apply the known binary continuous-charge base

The existing binary theorem proves positivity at length \(m\ge3\)
for \(|\delta|\le C_m\), where

\[
C_m=(6-4R)m+2-2R.
\]

The bound increases with \(m\), and \(C_9>2\). The last strict
inequality is equivalent to \(27>19R\), or the integer comparison
\(729>722\).

For \(n\ge9\) and \(0\le J\le2n-16\), choose (2) with
\(0\le k\le n-9\) and \(0\le\delta\le2\). This includes
the upper endpoint by taking \(k=n-9,\delta=2\).
Both lengths in (8) are at least nine, so their pressures are
strictly positive. This proves the separate folding corollary
\(G_n(J,t)>0\) for \(n\ge9\), \(|J|\le2n-16\).
Evenness covers negative charges. No new numerical sign certificate
is needed beyond the already reviewed continuous-charge theorem.

## 4. Preserve all tail endpoints for a sharper band

The exact continuation identity from the charged-floor theorem is
\(G_n(J,t)=H(J+2(n\bmod2),t)+T_{\rm phys}-T_{\rm cuts}\),
where \(H>0\) at every real charge and every scale. Put

\[
A=R(n+1),\qquad C=2n,\qquad B=|J|/R,\qquad D=C-A,
\qquad I(x)=\int_x^\infty f(u)\,du.
\tag{9}
\]

Keeping the two physical and four cut endpoints separately in the
rectangle and convex-midpoint tail bounds gives

\[
\begin{aligned}
2t(T_{\rm phys}-T_{\rm cuts})\ge{}&
2I((A-B)t)+2I((A+B)t)\\
&-I((C-|J|)t)-2I(Ct)-I((C+|J|)t).
\end{aligned}
\tag{10}
\]

Compare the ascending vectors

\[
v=(A-B,A-B,A+B,A+B),\qquad
w=(C-|J|,C,C,C+|J|).
\tag{11}
\]

The differences between their successive prefix sums, with the
prefix of \(v\) subtracted from that of \(w\), are

\[
D-|J|+B,\quad 2D-|J|+2B,\quad
3D-|J|+B,\quad4D.
\tag{12}
\]

All four are nonnegative provided
\(D\ge|J|(1-1/R)\). The second is automatic from \(D\ge0\),
because \(2B=R|J|\), and the third follows from the first.

The function \(I\) is decreasing and convex. The prefix comparison
implies \(\sum_i I(v_i t)\ge\sum_i I(w_i t)\). One direct proof
uses
\(\sum_i(s-v_i)_+=\max_k(ks-\sum_{i\le k}v_i)\): each maximum
for \(v\) is at least the corresponding maximum for \(w\).
Integrate this inequality against the nonnegative second derivative
of the decreasing convex function. Its value and derivative tend
to zero at positive infinity, so no affine remainder is present.

Equation (10) is therefore nonnegative whenever

\[
|J|\le\frac{(2-R)n-R}{1-1/R}
      =2n-2-2R.
\tag{13}
\]

For \(n\ge3\) this interval is nonempty. Restoring the strictly
positive continuation proves the leading result (1), for either
parity of \(n\). The comparison does not assert positivity in the
remaining edge bands.

## 5. Independent controls

The [checker](../tools/check_folded_charged_gaussian_floor.py) represents
the pressure as two exact signed-energy coefficient tables and a
rational ground coefficient. It applies
\(f_q(-x)=f_q(x)+(q-1)x\), then compares both sides of (4) as
formal expressions. This verifies the identity independently of the
digit count and without evaluating logarithms.

The [record](../evidence/folded-charged-gaussian-floor.json) includes
those exact identities, rational checks of the monotonicity
hypotheses, the integer \(C_9>2\) comparison, and exact
\(\mathbb Q(\sqrt2)\) prefix comparisons for the sharper tail
bound. It also independently evaluates the full pressure formula.
The floating comparisons check transcription; all sign conclusions
follow from the analytic proof and the existing continuation theorem.

```bash
python3 tools/check_folded_charged_gaussian_floor.py --output /tmp/folded-gaussian-replay.json
```
