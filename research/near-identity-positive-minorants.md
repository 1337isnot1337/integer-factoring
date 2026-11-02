# Path C: a positive Fourier minorant can bound the raw power search

**Status.** A valid sufficient near-identity hitting criterion and its
conditional factoring transfer. The proposed uniformly subexponential
**ambient mean** is now [excluded on an arithmetic semiprime family](ambient-minorant-norm-gap.md)
under the same pinned analytic input. The one-sided negative-spectrum
repair is excluded there as well. Subgroup-specific certificates remain
outside that obstruction. The complete factoring score stays **1/2**.

The [moment criterion](cofactor-orbit-moment-criterion.md) controls every
translated window of a target cycle. The actual
[smooth power selector](smooth-leading-power-selection.md) only needs
one successful exponent beyond the small least-power threshold. This
weaker requirement permits a different positivity argument. In contrast
to the [positive pair kernel](positive-relation-kernel.md), its positive
values must imply membership in the **original recognizable smooth set**,
not merely a quotient of two smooth classes.

The supported-function version is the classical finite-group Turan
extremal framework; see Kolountzakis and Revesz,
[*Turan's extremal problem for positive definite functions on groups*](https://arxiv.org/abs/math/0312218).
The argument here also permits nonpositive values outside the success set.
It is a proposed arithmetic interface, not a novelty claim for Fejer
positivity or positive-definite functions.

## 1. An exact finite-group hitting lemma

Let `G` be a finite abelian group, `S` an inversion-symmetric subset, and
`f:G->R` satisfy

\[
 f(1)=1,\qquad \widehat f(\chi)\ge0\quad\hbox{for every character},
 \qquad f(x)\le0\quad(x\notin S).
 \tag{1}
\]

Normalize Fourier coefficients by `|G|`, and let
`tau=widehat f(1)=|G|^{-1} sum_x f(x)>0`. Fourier inversion gives
`sum_chi widehat f(chi)=1`, hence `|f(x)|<=1`. Reality and nonnegative
Fourier coefficients imply `f(x^{-1})=f(x)`.

For `g in G` and integer `L>=0`, the exact Fejer identity is

\[
 \begin{aligned}
 \sum_{j=-L}^{L}\left(1-\frac{|j|}{L+1}\right)f(g^j)
 &=\sum_\chi\widehat f(\chi)\,
       \frac{\left|\sum_{u=0}^{L}\chi(g)^u\right|^2}{L+1}\\
 &\ge (L+1)\tau.
 \end{aligned}
 \tag{2}
\]

The identity follows by expanding the squared geometric sum and counting
the `L+1-|j|` pairs with difference `j`. Every Fourier contribution is
nonnegative. The principal character alone gives the last line; characters
trivial on the target subgroup only strengthen it. No subgroup-density
estimate or mixing assumption is used.

Fix `K>=1` and `L>=K`. If `g^j` is outside `S` for every `K<=j<=L`,
inversion symmetry gives the same assertion for negative exponents.
By (1), those summands in (2) are nonpositive. The remaining at most
`2K-1` summands each have value at most one. Therefore

\[
 \boxed{(L+1)\tau>2K-1\quad\Longrightarrow\quad
          g^k\in S\ \hbox{for some }K\le k\le L.}
 \tag{3}
\]

In particular, if an integer `D>=1` satisfies `tau>=1/D`, then
`L=(2K-1)D` works. The strict inequality is retained:
`(L+1)tau >= 2K-1+1/D`. The statement does not control arbitrary
translated windows; its center is the identity.

## 2. Conditional factoring consequence without evaluating f

Take `G` to be the proper class group, and `S` the classes with an exactly
decodable factor-base-smooth reduced leading norm at most the existing
reduced-form height `H`. This set is inversion symmetric. Useful proper
gcds encountered by the selector terminate the factoring step separately.

Suppose a uniform theorem supplies a bound `D(b)>=1` and guarantees that
**for every preprocessed b-bit input**
there exists an `f` satisfying (1) on its class group with
`tau>=1/D(b)`. The function itself need not be computed or evaluated.
Its existence proves that the unchanged raw selector, for every target
prime `p`, succeeds among

\[
 K_p\le k\le(2K_p-1)D(b),\qquad
 K_p=\min\{k:p^k>H\}.
 \tag{4}
\]

The algorithm tests consecutive powers until it recognizes a successful
leading norm, and decodes its row. It needs to compute neither the kernel
nor `D(b)` or the endpoint in (4); the theorem bounds its stopping time.
The existing weighted
diagonal-dominance argument and ambiguous-class extractor apply. Under
their inherited generation, preprocessing and extraction contracts, the
complete cost is

\[
 m^{O(1)}D(b)\operatorname{poly}(b,\log D(b),\log Y),
 \tag{5}
\]

including all rejected tests and the polynomial least-power thresholds.
For fixed `c`, a base cutoff `Y=2^{O((log b)^c)}` and a bound
`D(b)=exp(O(b/(log b)^alpha))` would therefore give saving exponent
`alpha`. The bound must be uniform over the inputs: measured inverse
means on finitely many cycles do not supply it. A fixed unknown absolute
constant in such a theorem does not obstruct this search-until-success
algorithm, which never uses it as a stopping parameter.

This removes the need to assume a lower-tail moment at every phase.
It replaces that obligation with a stronger structural property of the
original arithmetic success set. Neither property is proved here.

## 3. Density alone does not give this minorant

There is an elementary sharp counterexample even to a large global mean.
Let `G=Z/(2m)`, `m>=2`, and let `S` consist of zero and all odd residues.
Its density tends to `1/2`. For any function satisfying (1), write `E`
and `O` for its sums on even and odd residues. The off-set sign condition
gives `E<=f(0)=1`. Nonnegativity of the parity Fourier coefficient gives
`E-O>=0`. Consequently

\[
 \sum_x f(x)=E+O\le2E\le2,\qquad\tau\le1/m.
\]

Equality is attained by `f=1_{0}+(1_{1}+1_{-1})/2`, whose unnormalized
Fourier transform is `1+cos(2 pi k/(2m))>=0`. Thus constant density can
coexist with an inverse minorant mean as large as the group order.
This is a finite-group counterexample to a density-only implication,
not an arithmetic counterexample under the project's stronger contracts.

The known half-radius autocorrelation is a valid supported minorant:
if `A A^{-1} subset S`, take
`f(x)=|A intersect xA|/|A|`. Its mean is exactly `|A|/|G|`.
Substituting the larger density `|S|/|G|` for that mean is unjustified.
The existing construction supplies no lower bound of the strength
required in (5). The pair kernel built from the whole smooth set has
larger support `S S^{-1}`; using it in (3) would prove only a pair hit,
whose fast arithmetic extraction remains open.

## 4. A one-sided spectral bound would suffice

There is a concrete special case of the open minorant condition. Write
the real, unnormalized character sums as

\[
 F_S(\chi)=\sum_{x\in S}\overline{\chi(x)}.
\]

Suppose `B>=0` bounds only their negative part:
`F_S(chi)>=-B` for every character. Since `1 in S`, the function

\[
 f(x)=\frac{\mathbf1_S(x)+B\mathbf1_{\{1\}}(x)}{B+1}
\]

satisfies (1), with exact mean

\[
 \tau=\frac{|S|+B}{|G|(B+1)}\ge\frac{\delta}{B+1},
 \qquad\delta=\frac{|S|}{|G|}.
\]

Indeed its normalized Fourier coefficients are
`(F_S(chi)+B)/(|G|(B+1))>=0`; all other requirements are immediate.
Thus uniform subexponential bounds of the desired saving scale for both
`1/delta` and `B+1` would suffice through (3)--(5).

This asks for the **negative** spectrum, not a bound on absolute character
sums. The [absolute-value Parseval obstruction](leading-power-window-barrier.md)
does not itself lower-bound `B`: large positive coefficients can carry
that energy. Reusing an absolute bound would, however, retain its cost.
The [arithmetic norm-gap obstruction](ambient-minorant-norm-gap.md#4-large-global-density-and-the-negative-spectrum-consequence)
now forces an exponential negative bound on a semiprime family with
subexponential inverse density. The general ambient minorant condition
is excluded uniformly on that family too; this conclusion does not
assume that the minorant uses this identity correction.

The uniform ambient-minorant program is stopped by the
[arithmetic obstruction](ambient-minorant-norm-gap.md). A possible
continuation would need a large mean on the target subgroup, another
representation, or the original density/moment or pooled-rank obligations.
The exact criterion above remains valid and can still be useful for a
restricted input family when its hypotheses are separately established.

The [clean combined replay](../evidence/pathc-positive-minorant-full-verification.json)
passes the default gate and the focused character-lift gate at exact
`dd1b393235daa2c25f8529c0b3b405479d3be0aa`. Later changes add the
one-sided corollary and verification records in prose; source code is
unchanged. Documentation links are rechecked separately. No Lean or
independent mathematical review is claimed.
