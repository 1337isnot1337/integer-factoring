# Uniform charge bands for every homogeneous Boolean filter

## Result and scope

The [homogeneous one-level occupation filters](homogeneous-boolean-filter-bounds.md)
have strictly positive full critical pressure in each of these ranges:

1. Every exterior charge `|J|>=2n`, at every admissible interval length
   and every positive scale.
2. Every charge `|J|<=2n-2-2sqrt(2)` when `n>=17`, at every positive scale.

These statements permit every nonzero homogeneous multilinear Boolean
polynomial, including overlapping monomials. The degree and number of
selected modes can grow all the way to `n`. There is no fixed-degree
or fixed-charge qualification.

The physical construction requires power-of-two interval length, so
the second result applies to every such `n>=32`. At those lengths the
remaining charge region consists of two intervals of fixed width
`2+2sqrt(2)` immediately inside `J=+/-2n`. The
[read-once theorem](read-once-boolean-filter-floor.md) already excludes
that region for homogeneous read-once formulas, and excludes degree
one at every charge. General overlapping degree at least two remains
outside those additional conclusions.

The argument gives stronger ordinary obstructions within this filter
family. It provides no strict crossing witness and no improved
factoring exponent.

## 1. A degree-independent absolute loss cap

Retain the construction and notation of the earlier homogeneous note:
`n` is a power of two, `P` is positive and odd, `g` is a nonzero
homogeneous Boolean polynomial of degree `1<=d<=r<=n`, `L>nPd`, and
the selected physical energies are distinct members of
`{1,3,...,2n-1}`. Write

\[
f(x)=\log(1+e^{-x}),\qquad t>0.
\]

The exact relative cost is

\[
\mathcal R=(1-q+qa)(1-p+pH/a),\qquad
H=e^{2ntd},\qquad 0<a<1,
\tag{1}
\]

where `p` is the support probability under independent occupation
probabilities `(1+e^{t lambda_i})^{-1}` and `q` is the shifted support
probability. All charge dependence is contained in `q`.

The previous biased-support lemma and distinct physical energies give

\[
p\ge\prod_{i=0}^{d-1}
  \frac1{1+e^{(2n-1-2i)t}}.
\tag{2}
\]

The first factor of (1) is at least `a`, so `R>pH`. Multiplying (2)
by `H` and rearranging each factor gives

\[
\begin{aligned}
\mathcal R>pH
&\ge\prod_{i=0}^{d-1}
 \frac1{e^{-(2i+1)t}+e^{-2nt}}\\
&>(1+e^{-2nt})^{-d}
\ge(1+e^{-2nt})^{-n}.
\end{aligned}
\tag{3}
\]

Therefore the uniform strict loss cap is

\[
\boxed{\log\mathcal R>-n f(2nt).}
\tag{4}
\]

It holds for every real aggregate charge and every homogeneous pattern
in the stated construction. It includes the full polynomial degree at
infinity; no numerator or endpoint cost was dropped.

It remains to find charge ranges where the complete Gaussian pressure
is larger than the same quantity.

## 2. A stronger Gaussian floor at every exterior charge

Let `R=sqrt(2)` in this section; the filter ratio continues to be
written `mathcal R`. The
[exterior identity](exterior-charged-gaussian-floor.md) writes, for
`J=2n+c` with `c>=0`,

\[
G_n(2n+c,t)=\sum_{i=0}^{n-1}g_i(c,t),
\tag{5}
\]

where, with `b=2i+1`,

\[
g_i(c,t)=2f(bt)+f((c+2b-1)t)+f((c+2b+1)t)
         -2R f((c+2b)t/R).
\tag{6}
\]

Put `u=bt` and `alpha=2+c/b>=2`. Convexity gives

\[
g_i(c,t)\ge2[f(u)+f(\alpha u)-R f(\alpha u/R)].
\tag{7}
\]

The bracket increases with `alpha`, because its derivative is
`u[f'(alpha u)-f'(alpha u/R)]>=0`. Thus (7) is at least

\[
k(u)=2[f(u)+f(2u)-R f(Ru)].
\tag{8}
\]

The existing
[strict scalar inequality](charged-one-mode-binary-laurent-filters.md#5-positive-full-pressure-at-every-exterior-charge)
proves `2f(u)+f(3u)>2R f(Ru)` for all `u>0`. Since `f(2u)>f(3u)`,

\[
k(u)-f(2u)=2f(u)+f(2u)-2R f(Ru)>0.
\tag{9}
\]

Summing (7)--(9) and using convexity once more gives

\[
\boxed{
G_n(J,t)>\sum_{i=0}^{n-1}f((4i+2)t)
          \ge n f(2nt)\qquad(|J|\ge2n).
}
\tag{10}
\]

The last step is Jensen's inequality: the mean of `4i+2` is `2n`.
Evenness supplies the negative exterior. This is stronger than the
earlier exterior floor `G_n>f(t)` in the particular direction needed
to absorb an arbitrary homogeneous filter's loss.

Combining (4) and (10) proves strictly positive filtered pressure at
every exterior charge, with no restriction on the Boolean degree.
The scalar inequality used here already has its full analytic and
outward-integer proof; this argument adds no unproved sign estimate.

## 3. The finite-tail lower bound increases across its charge band

The exact
[finite-tail comparison](folded-charged-gaussian-floor.md#4-preserve-all-tail-endpoints-for-a-sharper-band)
uses

\[
A=R(n+1),\qquad C=2n,\qquad
I(x)=\int_x^\infty f(u)\,du.
\]

Define the retained lower numerator

\[
\begin{aligned}
E(J)={}&2I((A-J/R)t)+2I((A+J/R)t)\\
&-I((C-J)t)-2I(Ct)-I((C+J)t).
\end{aligned}
\tag{11}
\]

The positive periodic continuation gives

\[
G_n(J,t)>E(|J|)/(2t).
\tag{12}
\]

This statement retains the continuation charge shift for odd `n`;
its all-charge positivity makes that shift irrelevant to (12).

Put `h=-f'`. This function is positive and decreasing. Direct
differentiation of (11), using `A=(C+2)/R`, gives

\[
E'(J)=t^2\int_{C-J}^{C+J}
 \left[h((v+2)t/R)-h(vt)\right]dv.
\tag{13}
\]

If `0<=J<=2n-2-2R`, every integration point satisfies
`v>=2+2R=2/(R-1)`. Hence `(v+2)/R<=v`, and the integrand in
(13) is nonnegative. Therefore

\[
E(J)\ge E(0)=4\int_{At}^{Ct}f(u)\,du
\quad(0\le J\le2n-2-2R).
\tag{14}
\]

This improves a sign-only majorization argument to an explicit
lower bound throughout the whole band. Since `f` decreases,
(12)--(14) give

\[
\boxed{
G_n(J,t)>2[(2-R)n-R]f(2nt)
\quad(|J|\le2n-2-2R).
}
\tag{15}
\]

The coefficient in (15) exceeds `n` whenever
`(3-2R)n>2R`. This holds for every integer `n>=17`:
at seventeen the difference is `51-36R>0`, as follows from
`51^2=2601>2592=2*36^2`; the coefficient then increases with `n`.
Consequently

\[
G_n(J,t)>n f(2nt)
\quad(n\ge17,\ |J|\le2n-2-2R).
\tag{16}
\]

Equations (4) and (16) prove the second full-filter exclusion.

## 4. Exact controls and remaining scope

Run

```sh
python3 tools/check_homogeneous_boolean_charge_bands.py
```

The checker enumerates every nonzero homogeneous Boolean pattern
through four variables and verifies the support-to-loss calculation
using rational exponential scales. It includes patterns with
overlapping monomials. It also checks the exact physical-energy mean,
the scalar comparison relations, the charge-band endpoint ordering,
and the sharp integer coefficient threshold. The
[record](../evidence/homogeneous-boolean-charge-bands.json)
retains the completed counts and the exact proof dependencies.

The uniform theorem comes from (1)--(16), with the already established
Gaussian continuation and scalar inequality. The finite controls do
not replace those proofs. The small admissible lengths below thirty-two
retain their earlier results; no unchecked extension of the middle-band
claim to those lengths is asserted.

For lengths at least thirty-two, a possible witness in this homogeneous
one-level family must now use a polynomial outside the read-once theorem
and lie in

\[
2n-2-2\sqrt2<|J|<2n.
\]

The previous large-scale homogeneous bound also continues to apply in
that remaining region. Arbitrary multilevel occupation filters and
off-diagonal operators are separate families.

The proof uses the existing construction and pressure normalization at
revision `c2e951fa1f391edd97ba587d15ec13d3d6c8892d`. The new loss cap,
stronger Gaussian floors, and checker used substantial OpenAI assistant
support. No Lean source was changed.
