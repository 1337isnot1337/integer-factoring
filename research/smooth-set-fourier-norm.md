# The recognizable smooth set can have exponentially large Fourier norm

The new prime-form character estimates do not turn the recognizable
smooth set into a low-complexity Fourier test. Along an explicit family
of quadratic orders, its Fourier algebra norm is exponential in the
discriminant bit length, even though its density is subexponential.

More precisely, let `A` run through odd primes, put

$$\Delta=1-4A^4,\qquad D=|\Delta|,\qquad b=\mathrm{bits}(D),$$

and use

$$R=A^2-1,\qquad
 Y=\exp\bigl(\sqrt{\log D\log\log D}\bigr).$$

Let `S_R` consist of reduced primitive forms with leading coefficient
at most `R`, coprime to `Delta`, and `Y`-smooth. With normalized Fourier
coefficients on the proper class group `G`,

$$\|1_{S_R}\|_{\mathcal A}
 :=\sum_{\chi\in\widehat G}|\widehat{1_{S_R}}(\chi)|
 \ge 2^{b/128-o(b)}.\tag{1}$$

A varying norm-gap family in Section 5 sharpens the lower bound to
`2^(b/64-o(b))`.

This uses the same pinned finite-order Hecke zero-free theorem as the
[prime-form construction](polylog-prime-form-generators.md). It is an
obstruction to a specific deterministic hitting argument. It is not a
lower bound for factoring or for other algorithms that use the set.

## 1. A missed translate forces a Fourier norm lower bound

Let `G` be any finite abelian group, let `S` have density `delta>0`, and
let `mu` be a probability measure on `G`. Suppose

$$|\widehat\mu(\chi)|\le\epsilon\qquad(\chi\ne1),$$

where `0<epsilon`, and suppose a translate of its support misses `S`:

$$t\,\mathrm{supp}(\mu)\cap S=\varnothing.$$

Use the conventions

$$\widehat f(\chi)=\frac1{|G|}\sum_{g\in G}
 f(g)\overline{\chi(g)},\qquad
 \widehat\mu(\chi)=\sum_{g\in G}\mu(g)\chi(g).$$

Fourier inversion gives

$$0=\sum_u\mu(u)1_S(tu)
 =\delta+\sum_{\chi\ne1}\widehat{1_S}(\chi)
       \chi(t)\widehat\mu(\chi).$$

Taking absolute values proves

$$\delta\le\epsilon\sum_{\chi\ne1}
 |\widehat{1_S}(\chi)|,$$

and hence

$$\boxed{\ \|1_S\|_{\mathcal A}
 \ge\delta\bigl(1+\epsilon^{-1}\bigr).\ }\tag{2}$$

Equivalently, a proposed upper bound `M` for this norm gives a uniform
single-set hitting theorem from bias `epsilon` only when
`epsilon (M-delta)<delta`. A small density alone is insufficient.

## 2. An exponentially small bias supported on a missed translate

Set `X=A-1`. Give both orientations of every split odd prime `p<X`,
coprime to `Delta`, weight

$$\log(p)\frac pX\left(1-\frac pX\right),$$

and normalize. The
[smoothed prime estimates](explicit-prime-form-cutoffs.md) show that
the unnormalized character sums satisfy

$$A_\chi(X)=\delta_\chi X/6+
 O\bigl(X^{31/32}\log(D+3)
       +\sqrt X\log X+\log(D+3)\log X\bigr),\tag{3}$$

uniformly in `Delta` and the class-group character. Here `X` tends to
infinity with `A`, and `log D=4 log A+O(1)`. Thus the principal sum is
positive for every sufficiently large `A`, and the normalized measure
has all nonprincipal biases at most

$$\epsilon_A=C\frac{\log A}{A^{1/32}}\tag{4}$$

for an absolute effective constant `C`. Increasing `C` if necessary
absorbs the smaller error terms. For the asymptotic lower bound it is
enough to use the exact real weights in this mathematical measure; no
algorithm is being asked to enumerate its exponentially large support.

The target form

$$Q=(A,1,A^3)$$

misses `S_R` after multiplication by every form in this support.
Indeed, the [exact prime-translate calculation](prime-translates-miss-smooth-forms.md)
shows that the reduced product with a prime form of norm `p<A` has
leading coefficient exactly `Ap`. For sufficiently large `A`, we have
`Y<A`, so that coefficient is not `Y`-smooth. The radius and coprimality
conditions both hold; the failure is precisely smoothness.

Consequently (2) applies to this actual arithmetic set with `t=[Q]`
and the bias bound (4).

## 3. The density loses only a subexponential factor

For completeness, the density input needed here is much weaker than
an asymptotic formula for smooth forms. A direct lower bound suffices.

The principal case of (3), now at cutoff `Y`, gives at least

$$M\ge cY/\log Y$$

split odd primes at most `Y` and coprime to `Delta`, once `A` is large
enough. Put

$$k=\left\lfloor\frac{\log R}{\log Y}\right\rfloor.$$

Choose `k` distinct such primes and one of the two orientations at each
prime. Their product is a primitive ideal of squarefree norm
`a<=Y^k<=R`. Centering its middle coefficient gives a reduced form,
because `4R^2<D`. Different prime subsets give different norms;
different choices of orientations give different middle coefficients
modulo `2a`, by the Chinese remainder theorem. Thus

$$|S_R|\ge 2^k\binom Mk.\tag{5}$$

The elementary quadratic-order class-number upper bound is

$$|G|\le C_0\sqrt D\,\log^2(D+3).\tag{6}$$

One may even use a single logarithm. Write `Delta=f^2 D_K` and apply
the order class-number formula, including its unit index
([Daniels--Lozano-Robledo, Theorem 2.2](https://alozano.clas.uconn.edu/wp-content/uploads/sites/490/2014/01/Daniels_and_Lozano-Robledo_2.pdf),
quoting Cox, Theorem 7.24). Combining it with the imaginary quadratic
[class-number formula, Theorem 16.1](https://people.reed.edu/~jerry/361/lectures/iqclassno.pdf),
we obtain

$$h(\Delta)=\frac{w_\Delta\sqrt D}{2\pi}L(1,\chi_\Delta),$$

where `w_Delta<=6` and `chi_Delta` is the possibly imprimitive
Kronecker character modulo `D`. Its partial sums have absolute value
at most `D`, and its sum over a full period is zero. Splitting its
Dirichlet series after the first period and applying partial summation
to the tail gives `L(1,chi_Delta)<=log D+2`. This proves the stronger
bound `h(Delta)<=(3/pi)sqrt(D)(log D+2)` and, in particular, (6).

It follows, using `binom(M,k)>=(M/k)^k`, that

$$\begin{aligned}
 \log\frac{|S_R|}{|G|}
 &\ge k\log Y-\tfrac12\log D
       -O\bigl(k\log\log D+\log\log D\bigr)\\
 &\ge-O\left(\log Y+
             \frac{\log D}{\log Y}\log\log D\right).
\end{aligned}\tag{7}$$

Here `log R=(1/2) log D-log 2+o(1)`, and the floor defining `k` costs
at most `log Y`. Our specified `Y` therefore yields

$$\delta_A:=\frac{|S_R|}{|G|}
 \ge\exp\bigl(-O(\sqrt{\log D\log\log D})\bigr)
 =2^{-o(b)}.\tag{8}$$

This argument works for the nonmaximal orders as well. It neither
assumes that `D` is fundamental nor requires a factorization of `D`.

## 4. Combining the bounds

By (2), (4), and (8),

$$\begin{aligned}
 \log\|1_{S_R}\|_{\mathcal A}
 &\ge \log\delta_A+\log\epsilon_A^{-1}\\
 &\ge \frac1{32}\log A
       -O\bigl(\sqrt{\log D\log\log D}\bigr).
\end{aligned}$$

Since `log D=4 log A+O(1)` and `b log 2=log D+O(1)`, this is (1).
In particular, there is no uniform bound

$$\|1_{S_R}\|_{\mathcal A}\le2^{o(b)}$$

for these recognizable smooth sets. Combining a subexponential prime
sample with such a proposed bound cannot supply the missing universal
single-set extractor.

The positive correlation kernel remains different in a precise way.
For

$$C_S(g)=\frac1{|G|}\sum_x1_S(x)1_S(g^{-1}x),$$

its Fourier coefficients are `abs(fhat(chi))^2`. They are nonnegative,
and their sum is exactly `delta`, by Parseval. This is why the existing
[positive-kernel pair theorem](positive-relation-kernel.md) avoids the
large Fourier norm in (1). Finding its promised pair efficiently
remains the separate arithmetic task.

The proof is asymptotic and does not infer (1) from finite numerical
Fourier calculations. The underlying exact prime-translate identities
are replayed by `tools/check_prime_translate_obstruction.py`.

## 5. A larger norm gap doubles the exponent

The same argument gives a stronger family without new analytic input.
Let `A` again run through odd primes, and choose

$$r=\left\lfloor\sqrt{\log_2 A}\right\rfloor\ge2.$$

Now define

$$\Delta=1-4A^{2r+2},\qquad
 Q=(A,1,A^{2r+1}),\qquad
 R=A^{r+1}-1,\qquad X=A^r-1.$$

Use the same smoothness bound
`Y=exp(sqrt(log D loglog D))`. Since
`log D=(2r+2)log A+O(1)` and `r=Theta(sqrt(log A))`, we have `Y<A`
for sufficiently large `A`.

Form the weighted prime distribution at cutoff `X` as before, but
delete both orientations at `p=A`. For each remaining prime `p<X`,
the product with `Q` has primitive norm `Ap`. Its centered middle
coefficient gives a reduced form: `Ap<R` and `4R^2<D`. Thus its
leading coefficient is exactly `Ap`, and it lies outside `S_R`.
Deleting `p=A` is necessary here because the inverse orientation at
that prime would cancel the target.

Each deleted orientation has unnormalized weight at most
`(log A) A/X`. Dividing by the principal mass, which is asymptotic to
`X/6`, and renormalizing changes every character expectation by at
most

$$O\left(\frac{A\log A}{X^2}\right).$$

For `r>=2`, this is smaller than `X^(-1/32)log D`. The resulting
measure therefore still has nonprincipal character biases at most

$$\epsilon=O\bigl(X^{-1/32}\log D\bigr).$$

The proof of (8) is unchanged: `4R^2<D` and
`log R=(1/2)log D-log 2+o(1)`. Finally,

$$\log X=r\log A+o(1)
 =\left(\tfrac12-o(1)\right)\log D.$$

Using (2) gives

$$\boxed{\ \|1_{S_R}\|_{\mathcal A}
          \ge 2^{b/64-o(b)}.\ }\tag{9}$$

The sample in this proof is a mathematical certificate of a large
Fourier norm. Constructing its support would take exponential time;
no efficient sampler or extractor is claimed at that cutoff.

Substantial OpenAI assistant support was used. No Lean formalization
or stronger factoring score is claimed.
