# Slow endpoint corrections after ramification

## Result and scope

A polynomial crossing seed can be ramified by a large fixed integer while
its endpoint correction has degree one. This gives a constructive family
whose zero and infinity profiles have exact eventual formulas

\[
\kappa_i(P)=P\alpha_i+\beta_i.
\tag{1}
\]

The formula has an explicit threshold, obtained from the degree in the
correction parameter. Its finite `1/P` term can improve the rank cost at
a specified scale. An exact characteristic-three example below does so.
It does not improve the optimized interval bound or the factoring
exponent.

The generic limit cannot improve a retained background at fixed gauges.
That fact makes the finite correction term, and its sign, the relevant
search target. The note gives an exact way to compute that term without
treating the ramification degree as free input-dependent overhead.

## 1. A polynomial endpoint construction

Let `C(w)` be a finite four-boundary polynomial tensor over a field `K`.
Suppose its ordinary tensor powers, after the required matching physical
types, have the crossing endpoint at `w=-1`. This includes the truncated
symmetric Pascal tensors, and polynomial seeds with an untyped strict
crossing endpoint.

Let `B(w,b)` be a polynomial tensor in two scalar variables, with the
same physical alphabets, satisfying

\[
B(-1,0)=C(-1).
\tag{2}
\]

For a positive integer `P` such that `(-1)^P=-1`, set

\[
T_P(v)=B(v^P,1+v).
\tag{3}
\]

Then `T_P(-1)=C(-1)` entry by entry. Every tensor power has exactly the
same typed crossing endpoint. In characteristic two every positive
integer `P` is allowed; in odd characteristic it suffices that `P` is
odd. Taking `P=p^r` in characteristic `p` additionally gives
`1+v^P=(1+v)^P`, but this identity is not needed for (3).

For example, in the finite truncated algebra one can take

\[
B(w,b)=G(w)\left(1+\sum_j b^{e_j}w^{k_j}H_j\right),
\qquad e_j>0,\quad k_j\ge0.
\tag{4}
\]

Here the `H_j` are arbitrary nonconstant monomials. No sign restriction
on their occupied-index excess is needed: every correction vanishes
before the type restriction at `v=-1`. A negative-index monomial might
fail a separate *local index-gauged regularity* argument. That failure
does not invalidate the unweighted polynomial endpoint proof above.

Laurent powers of `v` can also be used if a single monomial clears their
poles. Its value at `v=-1` is a field unit, so a constant normalization
restores the endpoint. A scalar monomial shifts the two profiles by
opposite constants and cancels from their partition product. This
observation does not permit arbitrary rational denominators: additional
finite poles require their own cancellation proof.

### Gauges and finite construction

Fix any rational physical-letter gauges at zero and infinity. Under the
substitution `w=v^P`, multiply these gauges by `P`. This applies to
nonlinear letter weights as well as sums of one-particle weights. It
preserves the background's normalized pressure when the evaluation
scale is changed from `exp(t)` to `exp(t/P)`.

If the entries of `B` have degrees at most `d` in `w` and `s` in `b`,
the entries of (3) have degree at most

\[
Pd+s.
\tag{5}
\]

The `m`-copy degree is at most `m(Pd+s)`. For each fixed `P`, the local
tensor and all coefficient tables are finite and explicit, and the
usual fixed-seed transfer applies after a strict pressure certificate.
If `P` is allowed to depend on the factoring input, its construction
cost and coefficient count must be charged separately. This note does
not supply such a uniform variable-seed transfer.

## 2. Separated degree bands give exact minor valuations

Consider either polynomial cut matrix `M(w,b)` of `B`, of generic rank
`r` over `K(w,b)`. Assume every entry has `b`-degree at most `s`. For a
nonzero `k`-minor write

\[
\Delta(w,b)=\sum_{a=a_-}^{a_+}w^a f_a(b),
\qquad f_{a_-}f_{a_+}\ne0.
\tag{6}
\]

Every `f_a` has degree at most `ks`. If `P>ks`, the intervals
`[Pa,Pa+ks]` for distinct integer `a` are disjoint. Consequently the
substituted minor is nonzero, and

\[
\begin{aligned}
\mathrm{ord}_0\Delta(v^P,1+v)
 &=Pa_-+\mathrm{ord}_{b=1}f_{a_-}(b),\\
\mathrm{ord}_\infty\Delta(v^P,1+v)
 &=-Pa_+-\deg_b f_{a_+}(b).
\end{aligned}
\tag{7}
\]

The first secondary term lies in `[0,ks]`, and the second lies in
`[-ks,0]`. Equation (7) also proves that the specialized matrix has
generic rank exactly `r` when `P>rs`: an `r`-minor remains nonzero, and
every larger minor was already the zero polynomial.

These are exact polynomial statements. They do not depend on a sampled
value of `b` or on a finite series window.

## 3. A uniform threshold for the complete weighted profiles

Let `h_i` and `c_l` be the fixed row and column weights for this cut.
Let `L` be a common positive denominator, so all `Lh_i,Lc_l` are
integers. Give the substituted matrix row weights `Ph_i` and column
weights `Pc_l`.

For each nonzero minor, let `a_\sigma(\Delta)` and
`b_\sigma(\Delta)` denote the leading and secondary terms in (7), for
`\sigma=0,\infty`. Define

\[
\begin{aligned}
A_k^\sigma
 &=\min_{|I|=|J|=k}
 \left(a_\sigma(\Delta_{I,J})+
       \sum_{i\in I}h_i+\sum_{j\in J}c_j\right),\\
B_k^\sigma
 &=\min\{b_\sigma(\Delta_{I,J}):
       (I,J)\text{ attains }A_k^\sigma\}.
\end{aligned}
\tag{8}
\]

Zero minors are omitted, and `A_0=B_0=0`. The numbers `A_k` are the
weighted determinantal minima of `M(w,b)` over the coefficient field
`K(b)`, where `b` is a unit for the `w`-valuation.

**Affine profile theorem.** For every integer

\[
P>Lrs,
\tag{9}
\]

the actual weighted determinantal minima after substitution are

\[
\mu_k^\sigma(P)=PA_k^\sigma+B_k^\sigma.
\tag{10}
\]

Thus the complete ordered profile is exactly

\[
\kappa_i^\sigma(P)=P\alpha_i^\sigma+\beta_i^\sigma,
\qquad
\alpha_i^\sigma=A_i^\sigma-A_{i-1}^\sigma,\quad
\beta_i^\sigma=B_i^\sigma-B_{i-1}^\sigma.
\tag{11}
\]

**Proof.** All leading weighted minor values lie in `(1/L) Z`. A
nonminimal leading value therefore exceeds `A_k` by at least `1/L`.
The entire secondary range in (7) has width at most `ks`. Under (9),
no secondary advantage can compensate for this leading loss. First
minimize the leading value and then its secondary term, proving (10).
Successive differences of determinantal minima give (11). The rank
statement from Section 2 supplies every nonzero order. `\square`

The theorem is a finite search interface. One can enumerate minors for
a small seed, or use an exact weighted-minor oracle that retains its
attaining minors. Alternatively, once the threshold is proved, two
exact complete profiles beyond it determine the affine functions in
(11). Numerical fits before an established threshold do not do this.

The simple bound `|\beta_i^\sigma|\le rs` follows from the secondary
ranges above. Additional constant offsets in physical gauges can also
be retained: add their row and column sums to `b_\sigma` and increase
the threshold by the resulting finite secondary range.

## 4. The generic background floor and the finite opportunity

Suppose, more strongly than (2), that

\[
B(w,0)=C(w)
\tag{12}
\]

as an entire polynomial tensor, and that the two background cuts are
generically invertible. Every minor of `M(w,b)` has constant coefficient
in `b` equal to the corresponding background minor. Therefore its
generic `w`-valuation is at most the background valuation at zero, and
also at infinity. After adding the same weights and minimizing,

\[
A_k^\sigma\le \mu_{C,k}^\sigma.
\tag{13}
\]

This implies the partition inequality

\[
\sum_i e^{-t\alpha_i^\sigma}
\ge\sum_i e^{-t\kappa_{C,i}^\sigma},\qquad t>0.
\tag{14}
\]

For completeness, prefix-sum comparison implies
`\sum_i(u-\alpha_i)_+\ge\sum_i(u-\kappa_{C,i})_+`, since either sum
is the maximum over `k` of `ku` minus its first `k` orders. Integrating
against `t^2 e^{-tu}` proves (14). Thus the generic limit at fixed
gauges and fixed scale cannot beat its retained background.

The finite correction can nevertheless have the favorable sign. Write
the critical pressure with the common-type normalization as

\[
D_P(R,t)=
\sum_{\sigma=0,\infty}\log
 \sum_i e^{-t(\alpha_i^\sigma+\beta_i^\sigma/P)}
 -R\log\sum_a e^{-t g_a/R},
\tag{15}
\]

where `g_a` is the unscaled total physical energy. The paired-type
version replaces the last term by its two physical partition terms.
Let `\pi_i^\sigma` be the Gibbs weights proportional to
`e^{-t\alpha_i^\sigma}`. Then

\[
D_P=D_{\rm gen}
 -\frac{t}{P}\sum_{\sigma,i}\pi_i^\sigma\beta_i^\sigma
 +O(P^{-2}).
\tag{16}
\]

For fixed `t`, the absolute remainder is at most
`t^2 r^2 s^2/P^2` when both cuts have rank at most `r`. This follows
from the second derivative of log-sum-exp: it is the variance of
`\beta`, at most `r^2s^2` per cut.

A positive weighted sum of the secondary orders therefore decreases
the pressure at first order. It must overcome any generic background
gap. Near an interval whose background gap tends to zero, this finite
term remains a concrete search opportunity; (13) does not exclude it.

## 5. An exact qutrit control with a finite cost decrease

Work over `F3` in the algebra where each of `V,Y,U,X` has cube zero.
Let `G(w)` be the factorial-normalized truncated symmetric kernel of

\[
D(w)=\begin{pmatrix}1+w&w^2\\1&w\end{pmatrix},
\qquad H=V^2U^2.
\tag{17}
\]

For odd `P`, take the polynomial tensor

\[
N_P(v)=-\left(vG(v^P)+(1+v)H G(v^P)\right).
\tag{18}
\]

It has the original typed crossing at `v=-1`. This is (3) with
`B(w,b)=-((b-1)G(w)+bHG(w))`, so `B(w,0)=G(w)` and `s=1`.

Use the one-mode natural gauges from the truncated symmetric
construction, multiplied by `P`, and then make `v=-z^2`. The exact
zero profile is

\[
\bigl(0,\ P+2,\ P+2,\ 2P+2,\ 2P+2,\ 2P+2,
       3P+2,\ 3P+2,\ 4P+4\bigr),
\tag{19}
\]

and the infinity profile is the sorted list

\[
\bigl(-2P-2,\ -P-2,\ -P-2,\ -2,\ -2,\ 0,
       P-2,\ P,\ 2P\bigr).
\tag{20}
\]

These formulas for every positive odd `P` follow from the exact finite
certificate and the affine theorem. Before `v=-z^2`, the gauges have
denominator two, the rank is nine, and `s=1`; hence `P>18` suffices.
The checker proves (19)--(20) for every odd `P` through 21. The values
19 and 21 determine the affine profiles beyond the threshold; the
smaller values cover the remaining cases. Sorting is retained in (20)
for the small values of `P`.

Put `x=s^{-P}`, `r=s^2`, and `L=1+x+x^2`. After canceling the harmless
opposite scalar shifts from clearing the Laurent expression, the two
partition sums are

\[
\begin{aligned}
Z_0&=r+2x+3x^2+2x^3+x^4/r,\\
Z_\infty&=x^{-2}+2x^{-1}+2+r^{-1}+x+x/r+x^2/r.
\end{aligned}
\tag{21}
\]

The unperturbed product is `x^{-2}L^4`. At the exact rational choice
`P=3`, `s=11/10`, the ratio is

\[
\frac{Z_0Z_\infty}{x^{-2}L^4}
=\frac{99612598272277477054122841}
       {101096444552703336598632100}
<1.
\tag{22}
\]

This is about a `1.47%` decrease in the rank cost at those same gauges
and that same scale. It is not an improvement over the Gaussian's
optimized scale, and neither tensor here has pressure below the
half-rate threshold.

The secondary weighted sum in (16), for (18), is

\[
\frac{2(-1+x^2+x^3+2x^4)}{(1+x+x^2)^2},
\qquad x=e^{-t}.
\tag{23}
\]

It is positive for sufficiently small positive `t`, explaining the
finite improvement in (22). This example verifies that the finite
term is substantive even though its generic leading profile is the
Gaussian profile.

More precisely, (23) is

\[
\frac23-\frac{14}{9}t+O(t^2)
\tag{24}
\]

as `t` tends to zero. At the one-mode Gaussian optimum, however,
`x_*=(\sqrt{13}-1)/6` satisfies `3x_*^2+x_*-1=0`. Substituting this
relation into (23) gives exactly `-1/2`. Thus the pressure correction
at that optimum is `+t_*/(2P)+O(P^{-2})`, where `t_*=-\log x_*`.
The favorable small-scale correction is not a favorable correction
at the one-mode optimum. The Gaussian critical pressure there is
already strictly positive by the every-digit interval theorem.

## 6. Reproduction and search scope

Run

```sh
python tools/check_slow_ramification_crossings.py \
  --output evidence/slow-ramification-crossings.json
```

The checker independently enumerates bivariate minors of small random
matrices, verifies their separated valuations, compares complete
weighted profiles with power-series elimination, and verifies the
qutrit endpoint and profiles used in (22). All claimed comparisons are
over finite fields or rational numbers. Floating searches of additional
quartics and opposite-charge cubic pairs are separate experiments.

The reusable F4 constructor is
[`ramified_polynomial_crossings.py`](../tools/ramified_polynomial_crossings.py).
Strict polynomial Gaussian seeds and relabeling mixtures are supplied
by [`relabelled_gaussian_sums.py`](../tools/relabelled_gaussian_sums.py).
The background constructions are
[truncated symmetric Pascal crossings](truncated-symmetric-pascal-crossings.md)
and [every-digit interval frontier](truncated-symmetric-interval-frontier.md).
