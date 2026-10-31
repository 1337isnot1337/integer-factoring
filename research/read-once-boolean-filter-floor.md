# Read-once Boolean filters cannot improve the interior charged Gaussian

## Result and scope

Consider the homogeneous Boolean occupation filters from
[the exact filter construction](homogeneous-boolean-filter-bounds.md).
Suppose the nonzero homogeneous polynomial is represented by a formula
using XOR and AND, with each selected Boolean variable appearing once.
For every admissible interval length, every scale, and every aggregate
charge in the full interval `|J| <= 2n`, its complete two-cut cost is
strictly larger than the Gaussian cost at the same parameters.

This statement has no bound on the number of selected modes or on the
Boolean degree. Combining it with the
[global charged Gaussian floor](global-charged-gaussian-floor.md)
excludes a strict pressure witness from this family throughout that
charge interval.

For degree one, the polynomial is the XOR of any nonempty subset of
selected occupations. The exclusion extends to **every real charge**,
including exterior charges, at every scale and interval length allowed
by the construction. Thus correlating an arbitrary number of occupations
through a single parity does not escape the earlier one-mode obstruction.

These are ordinary mathematical obstructions to the stated filter family.
They establish no new factoring bound. Polynomials whose representations
reuse a variable, arbitrary multilevel filters, and the exterior-charge
case of degree at least two remain outside the general theorem below.

## 1. The exact cost ratio

Let `n` be a power of two, `P` a positive odd integer, and `m=nP`.
Select `r <= n` distinct indices `j_i` and put

\[
\lambda_i=2n-1-2j_i\in\{1,3,\ldots,2n-1\},\qquad t>0.
\]

Let `g` be a nonzero homogeneous multilinear polynomial of degree `d`
over the binary field. For an integer `L>md`, the filter is

\[
f_k(x)=1+(x^{md}+x^L)g(k),\qquad k\in\{0,1\}^r.
\tag{1}
\]

Every endpoint value is one. The complete polynomial identity in the
preceding construction gives the ratio

\[
\mathcal R=(1-q+qa)(1-p+pH/a),
\tag{2}
\]

where

\[
\begin{aligned}
p&=\Pr(g(X)=1),& \Pr(X_i=1)&=(1+e^{t\lambda_i})^{-1},\\
q&=\Pr(g(Y)=1),& \Pr(Y_i=1)&=(1+e^{t(\lambda_i+J)})^{-1},\\
a&=e^{-2t(L-md)/P}\in(0,1),& H&=e^{2ntd}>1.
\end{aligned}
\tag{3}
\]

All coordinates in each product distribution are independent. The
physical paired partition cancels in this ratio because the comparison
uses the same charge, scale, and physical alphabet. The full filtered
critical pressure is `G_n(J,t)+log R`.

## 2. A log-odds bound for read-once formulas

For a probability `u` in `(0,1)`, write
`ell(u)=log(u/(1-u))`. A read-once formula here has Boolean variables
at its leaves and XOR or AND at its internal nodes. Its two child
formulas use disjoint variable sets. This disjointness makes their
outputs independent under every product measure.

**Lemma.** Let `g` be a nonconstant read-once formula of algebraic degree
`d` over the binary field. Change each input log-odds by the same real
number `theta`. If `p(theta)` is the output-one probability, then

\[
\left|\frac{d}{d\theta}\ell(p(\theta))\right|\le d.
\tag{4}
\]

Consequently, for every pair of real shifts,

\[
|\ell(p(\theta_1))-\ell(p(\theta_0))|
\le d|\theta_1-\theta_0|.
\tag{5}
\]

The lemma does not require homogeneity. Homogeneity is needed in (1)
for the exact two-cut identity (2).

**Proof.** A leaf has derivative one. We prove the induction step for
two child probabilities `u,v` and their output probability `w`.

For AND, `w=uv`, and the derivatives with respect to the child log-odds
are

\[
\frac{\partial\ell(w)}{\partial\ell(u)}
=\frac{1-u}{1-uv}\le1,\qquad
\frac{\partial\ell(w)}{\partial\ell(v)}
=\frac{1-v}{1-uv}\le1.
\tag{6}
\]

Both are nonnegative. The chain rule bounds the absolute derivative
by the sum of the two child degrees. The algebraic degree of their
AND is that sum because their variables are disjoint.

For XOR, `w=u+v-2uv`. Put `A=|1-2u|` and `B=|1-2v|`. Both lie in
`[0,1)`. Direct differentiation gives

\[
\begin{aligned}
\left|\frac{\partial\ell(w)}{\partial\ell(u)}\right|
+\left|\frac{\partial\ell(w)}{\partial\ell(v)}\right|
&=\frac{(1-A^2)B+(1-B^2)A}{1-A^2B^2}\\
&=\frac{A+B}{1+AB}\le1.
\end{aligned}
\tag{7}
\]

The final inequality is equivalent to `(1-A)(1-B)>=0`. The chain rule
therefore bounds the derivative by the larger child degree. That is
the degree of the XOR: the disjoint nonconstant polynomials have no
common nonconstant monomial to cancel. This proves (4); integration
proves (5). All formula probabilities lie strictly between zero and
one, so the differentiations are valid. \(\square\)

## 3. Every interior charge is excluded

The change from `X` to `Y` in (3) shifts every input log-odds by `-Jt`.
The lemma implies

\[
\ell(q)-\ell(p)\le d|J|t\le2ntd=\log H
\qquad(|J|\le2n).
\tag{8}
\]

Equivalently,

\[
q(1-p)\le Hp(1-q).
\tag{9}
\]

Regard (2) as a function of the real parameter `a` on `(0,1]`:

\[
\mathcal R(a)=(1-q)(1-p)+pqH+q(1-p)a+pH(1-q)/a.
\]

Equation (9) gives

\[
\mathcal R'(a)=q(1-p)-pH(1-q)/a^2\le0.
\]

Therefore

\[
\boxed{\mathcal R(a)\ge\mathcal R(1)=1-p+pH>1.}
\tag{10}
\]

The actual filter has `0<a<1`; permitting `a=1` only enlarges the
minimization domain. The strict inequality follows from `p>0` and
`H>1`, including the charge endpoints. The global Gaussian theorem
now proves strictly positive full critical pressure for every
homogeneous read-once filter in this interior charge interval.

Examples covered by (10) include a parity of disjoint equal-degree
monomials, a product of parities on disjoint variable blocks, and any
homogeneous nesting of those operations. Overlapping monomials do
not automatically satisfy the read-once hypothesis.

## 4. Parity filters are excluded at all exterior charges too

For degree one, write `g(X)` as the XOR of its nonempty selected
subset. With `u_i=Pr(X_i=1)` from (3), independence gives

\[
1-2p=\prod_i(1-2u_i).
\tag{11}
\]

Here every `u_i` lies in `(0,1/2)`. The product in (11) is no larger
than any one factor, so `p>=u_i` for every selected index.

At every charge, the first factor of (2) is at least `a`, and hence

\[
\mathcal R>a(1-p)+pH>pH.
\tag{12}
\]

For any selected `lambda_i` and `d=1`,

\[
pH\ge\frac{e^{2nt}}{1+e^{\lambda_i t}}
=\frac{e^{(2n-\lambda_i)t}}{1+e^{-\lambda_i t}}
>\frac1{1+e^{-t}},
\tag{13}
\]

because `1<=lambda_i<=2n-1`. Thus the universal strict loss cap is

\[
\log\mathcal R>-\log(1+e^{-t}).
\tag{14}
\]

The earlier
[exterior Gaussian estimate](charged-one-mode-binary-laurent-filters.md#5-positive-full-pressure-at-every-exterior-charge)
proves `G_n(J,t)>log(1+e^{-t})` for `|J|>=2n`. Combining it with
(14) makes the full filtered pressure positive in the entire exterior.
Section 3 covers the interior, proving the all-charge parity result.

Relative gains at exterior charges remain possible. The conclusion
concerns the full pressure with the physical partition retained, not
the false claim that every filter ratio is at least one everywhere.

## 5. Verification and remaining question

Run the independent exact checker with

```sh
python3 tools/check_read_once_boolean_filter_floor.py
```

It checks the gate derivative identities using rational probabilities,
enumerates homogeneous read-once formulas on small variable sets,
compares their recursive probabilities with direct truth-table sums,
and checks the common-shift derivative bound and the interior odds
inequality. It also reconstructs the literal binary polynomial filter,
checks its full finite-difference transform and endpoint, and verifies
the exact ratio and exterior loss cap on rational exponential scales.
The counts and a strict exterior relative-gain example are retained in
[the replay record](../evidence/read-once-boolean-filter-floor.json).
The infinite-parameter conclusions follow from the proofs above.

The open generalization is an appropriate probability or profile bound
for homogeneous polynomials that reuse variables. An unsupported
extension of (4) to all binary polynomials is not used here. The
exterior region for degree at least two also needs a stronger absolute
comparison than the degree-one loss cap.

This result builds on the cited exact filter construction and Gaussian
floor at repository revision `c2e951fa1f391edd97ba587d15ec13d3d6c8892d`.
The new log-odds argument and checker were developed with substantial
OpenAI assistant support. No Lean source was changed.
