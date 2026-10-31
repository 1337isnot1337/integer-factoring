# Homogeneous Boolean occupation filters: exact profiles and pressure bounds

## Result and scope

There is an explicit correlated polynomial filter for every nonzero
homogeneous Boolean polynomial on the selected occupation digits.
Its exact two-cut cost reduces to two support probabilities. An
elementary biased support bound then gives three conclusions:

1. At aggregate charge \(J=0\), every nonconstant filter in this
   family costs strictly more than the Gaussian at every scale.
2. At every real \(J\), the same strict comparison holds above an
   explicit scale threshold, uniformly over the selected indices,
   polynomial degree, and Boolean pattern.
3. Below that threshold, the full critical pressure is positive for
   bounded \(|J|\), fixed Boolean degree, and sufficiently large
   interval length. This uses the exact positive finite Gaussian tail.

The remaining charged large-scale absolute floor depends on the
charged Gaussian itself. No theorem that \(J=0\) minimizes Gaussian
pressure over all real charges is assumed. The conclusions concern
this homogeneous one-level family; arbitrary multilevel polynomial
filters and nonhomogeneous cancellations remain outside its scope.

## 1. A correlated family with a strict endpoint

Let \(n\) be a power of two, \(P\) a positive odd integer,
\(m=nP\), and select \(r\le n\) distinct physical identity indices.
Write the occupation states as \(k\in\{0,1\}^r\). Let

\[
g(k)=\sum_{|A|=d} c_A\prod_{i\in A}k_i,
\qquad c_A\in\mathbf F_2,
\qquad 1\le d\le r,
\tag{1}
\]

be a nonzero homogeneous multilinear Boolean polynomial. Choose an
integer \(L>md\), and reverse the alpha indices by setting

\[
f_k(x)=\alpha_{\mathbf1-k}(x)
       =1+(x^{md}+x^L)g(k).
\tag{2}
\]

Every endpoint value is one. Thus the
[prime-power identity construction](prime-power-joint-identity-filters.md)
has its original endpoint without an additional scalar correction.
The full degree is \(L\) on the support of \(g\) and zero off that
support. This degree is charged at infinity.

Put \(C=1+x^m\). The exact binary transform is

\[
\beta_k=\sum_{\ell\le k}C^{|k|-|\ell|}f_\ell.
\tag{3}
\]

For a monomial \(k_A\) of degree \(d\), its transform is zero unless
\(A\subseteq\mathrm{supp}(k)\), and otherwise equals
\((1+C)^{|k|-d}\). Homogeneity makes that power the same for all
monomials in (1). Consequently

\[
\boxed{
\beta_k(x)=
\begin{cases}
x^{m|k|},&g(k)=0,\\
x^{m(|k|-d)+L},&g(k)=1.
\end{cases}}
\tag{4}
\]

The second case only occurs when \(|k|\ge d\). These are complete
polynomial identities, not leading-term approximations.

## 2. The exact two-cut ratio

For the selected indices \(j_1,\ldots,j_r\), put

\[
\lambda_i=2n-1-2j_i,
\qquad t=P\log s,
\qquad
\mu_\lambda(k)\propto e^{-t\lambda\cdot k}.
\tag{5}
\]

Define

\[
p_g=\Pr_{\mu_\lambda}(g=1),\qquad
q_g=\Pr_{\mu_{\lambda+J\mathbf1}}(g=1),\qquad
a=e^{-2t(L-md)/P},\qquad H=e^{2ntd}.
\tag{6}
\]

Here \(0<p_g,q_g<1\), \(0<a<1\), and \(H>1\).
Substituting (4) and the full alpha degrees into the
[joint expectation identity](joint-binary-filter-synthesis.md) gives

\[
\boxed{R=(1-q_g+q_g a)(1-p_g+p_gH/a).}
\tag{7}
\]

This is the ratio to the Gaussian cost at the **same** charge and
scale. The physical partition therefore cancels from the ratio.

At \(J=0\), the two probabilities coincide. Cauchy--Schwarz gives

\[
R\ge(1-p_g+p_g\sqrt H)^2>1.
\tag{8}
\]

This proves strict domination at every scale for the natural charge.
Together with the
[natural Gaussian pressure floor](truncated-symmetric-interval-frontier.md),
it excludes a threshold improvement from this family at \(J=0\),
for every admissible interval length.

## 3. A biased Boolean support lemma

**Lemma.** Let \(g\) be a nonzero multilinear polynomial over
\(\mathbf F_2\) of degree at most \(d\). Under independent Boolean
digits with probabilities \(0<u_i\le1/2\),

\[
\Pr(g=1)\ge\prod_{i=1}^{d}u_{(i)},
\tag{9}
\]

where \(u_{(1)}\le\cdots\le u_{(r)}\) are the sorted probabilities.
For degree zero the empty product is one.

**Proof.** It is enough to use the actual degree; a smaller degree
gives a larger product lower bound. For positive degree, choose a
variable in a monomial of maximum degree and write
\(g=g_0+k_i g_1\), with \(g_1\ne0\) and degree \(d-1\).
Whenever \(g_1=1\), exactly one of the two values obtained by varying
\(k_i\) equals one. Its conditional probability is at least \(u_i\).
Induction gives \(u_i\) times the product of the \(d-1\) smallest
remaining probabilities. This is a product of \(d\) distinct
probabilities, hence at least the right side of (9). \(\square\)

The lemma applies to nonhomogeneous Boolean polynomials as well.
Homogeneity is needed for the exact cancellation identity (4).
It is not asserted for arbitrary extension-field coefficient patterns.

## 4. Uniform comparison over every aggregate charge

The first factor of (7) is at least \(a\). Therefore

\[
R\ge a(1-p_g)+p_g H>p_gH.
\tag{10}
\]

The physical infinity probabilities are
\(u_i=(1+e^{t\lambda_i})^{-1}<1/2\). Distinct physical indices imply
that their \(d\) smallest probabilities are bounded below by the
probabilities at energies
\(2n-1,2n-3,\ldots,2n-(2d-1)\). Combining (9) and (10),

\[
\boxed{
\log R>
td^2-\sum_{i=0}^{d-1}
\log(1+e^{-t(2n-1-2i)}).
}
\tag{11}
\]

Using only \(\lambda_i\le2n-1\) gives the simpler bound

\[
\boxed{
\log R>d\,[t-\log(1+e^{-(2n-1)t})].
}
\tag{12}
\]

Let \(t_n^*>0\) be the unique solution of

\[
t_n^*=\log(1+e^{-(2n-1)t_n^*}),
\quad\text{equivalently}\quad
e^{2nt_n^*}(1-e^{-t_n^*})=1.
\tag{13}
\]

The bracket in (12) is strictly increasing. Hence \(R>1\) for
every real \(J\) whenever \(t\ge t_n^*\). This conclusion is uniform
over \(r,d,g,L,P\) and the selected indices, within the construction's
stated hypotheses.

There can be relative gains below this scale at extreme charges.
For example, \(n=32,P=1,r=d=1,j_1=0,L=266,J=-1000,t=0.01\)
has \(R<1\). Its full critical pressure remains strongly positive.
The exact profiles and floating comparison are replayed in the
[checker record](../evidence/homogeneous-boolean-filter-bounds.json).
Thus (13) is a scope restriction with mathematical content, not an
omitted all-scale relative theorem.

## 5. The small-scale branch at bounded charge

Let \(G_n(J,t)\) denote the complete binary Gaussian critical
pressure at \(\rho=1/\sqrt2\), with the same normalization as
[`gaussian_pressure`](../tools/prime_state_filter_pressure.py).
The filtered pressure is
\(F_n(J,t)=G_n(J,t)+\log R\).

For binary prime-power \(n\ge4\), the positive infinite continuation
and finite-tail bounds in the natural interval theorem imply

\[
G_n(0,t)>
\frac2t\int_{\sqrt2(n+1)t}^{2nt}\log(1+e^{-x})\,dx
\ge
\frac{2\log2}{t}
\left(e^{-\sqrt2(n+1)t}-e^{-2nt}\right).
\tag{14}
\]

The final inequality uses concavity of \(\log(1+y)\) on
\(0\le y\le1\), which gives \(\log(1+y)\ge y\log2\).

The Gaussian zero-cut log partition is convex in \(J\). Its physical
partition has \(2n\) binary factors, enters with coefficient
\(-\sqrt2\), and uses scale \(t/\sqrt2\). Since the second derivative
of \(\log(1+e^{-x})\) is at most \(1/4\),

\[
\frac{\partial^2 G_n}{\partial J^2}
\ge-\frac{nt^2}{2\sqrt2}.
\]

The pressure is even in \(J\). Integrating this lower bound twice
therefore gives, for every real \(J\),

\[
G_n(J,t)\ge G_n(0,t)-\frac{nt^2J^2}{4\sqrt2}.
\tag{15}
\]

Put

\[
T_n=\frac{\log(2n)}{2n-1}.
\tag{16}
\]

Then \(t_n^*<T_n\), because
\(e^{(2n-1)T_n}(e^{T_n}-1)>2nT_n>1\).
Also

\[
\frac{e^{-\sqrt2(n+1)t}-e^{-2nt}}t
=\int_{\sqrt2(n+1)}^{2n}e^{-ut}\,du
\]

decreases with \(t\). Equation (12) gives
\(\log R>-d\log2\). Hence for
\(0<t\le t_n^*\) and \(|J|\le J_0\),

\[
\boxed{
F_n(J,t)>
\frac{2\log2}{T_n}
\left(e^{-\sqrt2(n+1)T_n}-e^{-2nT_n}\right)
-d\log2-\frac{nT_n^2J_0^2}{4\sqrt2}.
}
\tag{17}
\]

This explicit bound is uniform over the selected subset, the number
of selected modes, the homogeneous pattern, \(L\), and \(P\).
Its first term is asymptotic to

\[
2\log2\,
\frac{(2n)^{1-1/\sqrt2}}{\log(2n)},
\tag{18}
\]

which diverges. For fixed \(d,J_0\), the last term tends to zero.
Thus (17) proves strictly positive full critical pressure throughout
the small-scale branch for all sufficiently large binary prime powers.

Combined with Section 4, any hypothetical new charged threshold
crossing from this family, at those interval lengths and bounded
charges, would have to occur where the Gaussian at the same charge
already has negative critical pressure. This is a relative reduction
to the charged Gaussian question. It is not a proof of an absolute
charged-Gaussian floor.

The subsequent [charged Gaussian floor](charged-gaussian-pressure-floor.md)
closes that gap for bounded charge. For fixed Boolean degree and a fixed
charge bound, this entire family has positive full critical pressure at
every scale for all sufficiently large binary prime-power lengths.

## 6. Exact replay

The [implementation](../tools/homogeneous_boolean_filter_bounds.py)
constructs the literal binary alpha polynomials and both profiles.
Run:

```sh
python tools/check_homogeneous_boolean_filter_bounds.py
```

The [record](../evidence/homogeneous-boolean-filter-bounds.json) checks
68,031 exact biased-support inequalities, including every nonzero
Boolean polynomial in up to four variables; 8,892 complete polynomial
transform identities for every homogeneous pattern in up to five
variables; and 276,520 endpoint equalities. Eighty independent
complete-profile cost comparisons agree with (7), with maximum
floating discrepancy below \(8\cdot10^{-15}\).

The analytic inequalities are proved above. Floating evaluations of
(17) in the record illustrate the scale of the bound; they are not
used as interval-certified cutoffs in the theorem.
