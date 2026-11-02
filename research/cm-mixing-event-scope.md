# CM mixing and the arithmetic event-correlation gap

## Audit result

The available CM mixing theorems do not directly supply the arithmetic
correlations required by the [static p-free flag](event-word-boundary.md)
or [finite-order event criterion](connected-event-window-criterion.md).
Their averaging sets, hypotheses, error scales and test functions must
all be transferred. The diagnostic was
[registered before this audit](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19360397).
This extends the earlier [literature routing](window-lemma-character-duality.md#4-literature-searched-2026-11-01),
which already distinguishes sparse suborbits from random-walk expansion.
No complete factoring improvement follows; the score remains **1/2**.

There is also an elementary, explicit test-function obstruction: on
arbitrarily large order discriminants, any geometric function agreeing
with the ambient truncated smooth flag at all CM points must have
Lipschitz constant at least `sqrt(D/20)`. Thus a subexponential
test-function budget cannot be assumed. This does not exclude arithmetic
weighted trace formulae, functions constructed only on a particular
target orbit, or stronger estimates that charge the actual norm loss.

## 1. Exact primary-source scopes

**Khayutin, Joint equidistribution of CM points, arXiv1710.04557v3,
Theorem 1.3** ([versioned primary text](https://arxiv.org/html/1710.04557v3#S1.SS2.SSS1)).
The averaging set is a full toral packet, paired with a class-group
translate. Discriminants are fundamental; two fixed primes split in the
CM fields, there is no exceptional Landau--Siegel zero, and the minimum
ideal norm of the shifts tends to infinity. The conclusion is weak-*
equidistribution, with component restrictions at general level. Its
fixed-product extension does not state uniform growing-order arithmetic
cumulant bounds. This is not a theorem about short consecutive powers
inside an arbitrary prime-generated cyclic subgroup.

**Blomer--Brumley--Khayutin, The mixing conjecture under GRH,
arXiv2212.06280v2, Theorems 2.1 and 4.2**
([versioned primary text](https://arxiv.org/html/2212.06280v2#S4.SS4)).
For fundamental discriminants and almost maximal level, the effective
same-packet statement assumes GRH and, in the indefinite case, Ramanujan.
The main compact quaternionic setting uses a non-split algebra. If `q`
is the minimum ideal norm of the shift, the displayed error envelope is

```math
h_\varepsilon(q,D)=
\begin{cases}
q^{-1/2+\varepsilon}+q^{1/2}D^{-1/4+\varepsilon},
 &q\le D^{1/2-\varepsilon},\\
(\log D)^{-1/18+\varepsilon},&q\ge D^{1/2-\varepsilon}.
\end{cases}                                                   \tag{1}
```

Theorem 4.2 tests the discrete spectrum, and its implied constant depends
on the test function through a Sobolev norm. Section 2.3 describes an unconditional small-shift range using
subconvexity; that does not remove the hypotheses of the global theorem.
The repository's [fixed Hecke zero-free strip](polylog-prime-form-generators.md)
does not imply GRH or Ramanujan. We have not audited the imported
automorphic proofs independently.

**Blomer--Brumley--Radziwill, Joint Linnik problems,
arXiv2603.05609v1, Theorem 2.5 and Section 3.1**
([versioned primary text](https://arxiv.org/html/2603.05609v1#S2.SS3)).
Here fundamental-discriminant packets average the whole class group,
with coordinates acted on by `[a]` and `[a]^nu`. For `nu=1` the two
quaternion algebras must be non-isomorphic; for `nu=2` the first is
non-split. The quantitative quadratic-character zero-free hypothesis is
compatible with a fixed strip by choosing, for example,
`psi(D)=sqrt(log D)`. The errors are `psi(D)^(-rho)` and
`(log psi(D))^(-1/2)`, respectively, with test-function dependence;
Section 3.1 charges polynomial dependence on spectral parameters.
This removes important hypotheses for distinct Linnik problems, but is
not same-speed mixing on one modular curve or an arbitrary prime cycle.
It is not a uniform theorem for all order discriminants.

## 2. Rates and averaging must match the event problem

Suppose the density scale under study is
`delta_b=exp(-c b/(log b)^alpha)`, with fixed positive `c,alpha`,
and `b=bits(D)`. For every fixed `A>0`,

```math
\frac{b^{-A}}{\delta_b^2}
 =\exp\left(\frac{2cb}{(\log b)^\alpha}-A\log b\right)
 \longrightarrow\infty.                                      \tag{2}
```

Therefore the logarithmic error in (1), even with bounded test-function
norm, does not certify an error below the product-density scale in this
regime. The displayed rate is an upper envelope, not a lower bound on
the true error: (2) stops this direct use of the theorem, not the
possibility of stronger arithmetic cancellation. It also does not claim
that every actual target density is as small as the currently guaranteed
scale. The newer joint theorem's displayed errors have the same
insufficiency at that scale, even before its different ambient setting.

Restriction to an actual target subgroup `U<=G` needs its own work.
For any class functions `f,g` and shift `c in U`, orthogonality gives

```math
\mathbb E_{x\in U} f(x)g(cx)
 =\sum_{\psi\in U^\perp}\mathbb E_{x\in G}\psi(x)f(x)g(cx).
                                                               \tag{3}
```

There is no missing reciprocal-index factor in this normalized identity.
Even a uniform twisted ambient error `epsilon` could lose a factor
`[G:U]` on summing (3). Untwisted packet equidistribution alone does not
control it. Finally, two-point mixing does not imply the connected
estimates through order `O(log |G|)` required by the event criterion.

## 3. Adjacent CM points force a large interpolation norm

Fix any `alpha>0`. For large integers `k`, set

```math
b=2k+5,\quad X=2^k,\quad
Y=\lceil\exp((\log b)^{1+\alpha})\rceil.
```

Here `Y=o(X)`. Choose a prime `q` in `(Y,2Y]`, using Bertrand's
postulate. The next multiple of `q` after `X` is not `Y`-smooth and
lies at most `2Y` steps away. Since `X` is smooth, there is an integer
`a in [X,X+2Y)` for which `a` is smooth and `a+1` is not. One may
find the first such boundary by trial division through `Y`; this has
cost polynomial in `Y,b`, with no general factoring oracle.

Set

```math
D=20a(a+1)-1,\quad
Q_a=(a,1,5(a+1)),\quad Q_{a+1}=(a+1,1,5a).
```

For all sufficiently large `k`, `bits(D)=b`, so the chosen cutoff is
self-consistent. Both forms are primitive strict reduced forms of
discriminant `-D`, with leading norms coprime to `D`. Also
`a+1<=floor(sqrt(D)/4)` for `a>=5`. Their truncated smooth flags are
therefore one and zero at the same safe radius used by the event note.

The associated points in the standard fundamental region are

```math
z_a=\frac{-1+i\sqrt D}{2a},\qquad
z_{a+1}=\frac{-1+i\sqrt D}{2(a+1)}.
```

They stay in a fixed compact region around height `sqrt(5)`. The exact
hyperbolic distance formula gives, for their distance in the upper
half-plane,

```math
\cosh d(z_a,z_{a+1})-1
 =\frac{D+1}{2Da(a+1)}=\frac{10}{D},\qquad
d(z_a,z_{a+1})\le\sqrt{20/D}.                          \tag{4}
```

The quotient distance is no larger. A geometric function approximating
both flags to error at most `eta<1/2` must consequently satisfy

```math
\operatorname{Lip}(F_D)\ge(1-2\eta)\sqrt{D/20}.          \tag{5}
```

The same lower bound, up to a fixed constant, holds for any chosen
Sobolev norm controlling the `C^1` norm on this fixed compact region.
We do not assume without inspection which exact norm a source theorem
uses. Equation (5) rules out a free subexponential Lipschitz interpolation
of the ambient arithmetic flag.

This is a full-interpolation statement. It does not rule out a one-sided
geometric minorant that discards these particular smooth points; the
retained mass and norm of such a minorant would need separate proofs.
The [ambient Fourier-positive obstruction](uniform-ambient-minorant-obstruction.md)
has its own positivity hypotheses and is not automatically a theorem about
every geometric minorant.

Under the existing principal split-prime band estimate, this example
also applies to an ambient **p-free** flag for some target prime `p<=Y`.
The band has `>>Y/log Y` unramified split primes, while at most
`O(b/log Y)` of them divide `a`. Since `Y/b` tends to infinity, choose
one not dividing `a`. The first point then has p-free flag one, while
the rough second point still has flag zero. The target prime may vary
with `D`. This does not put both points on its cyclic subgroup, and does
not exclude a different interpolant defined only on that orbit.

The constructed discriminants are order discriminants; no infinite
fundamental-discriminant or hard-composite subfamily is claimed here.
Thus (5) is an obstruction to a uniform all-order interpolation step,
not a contradiction to the fundamental-discriminant mixing theorems.

## 4. Evidence and next mechanism

```sh
make verify-cm-smooth-flag-geometry
```

The [checker](../tools/check_cm_smooth_flag_geometry.py) verifies five
fixed adjacent-form constructions, exact reducedness, coprimality,
safe-radius membership, opposite flags, and the rational distance
identity. Four controls have an eligible p-free target; a small-cutoff
control has none, guarding the large-cutoff analytic qualification.
The [receipt](../evidence/cm-smooth-flag-geometry.json) pins the source.
There is no class-cycle, discriminant-range or cutoff sweep. These
controls do not reverify any external automorphic theorem.

The remaining potentially useful route is an **arithmetic weighted**
period or correlation estimate that handles the smooth-word coefficients
directly, charges subgroup restriction, and reaches the required orders
and rarity scale. A bounded-complexity geometric mixing import is stopped
by the audited gaps. No new conditional cumulant theorem is asserted.

The cited results retain their authors' credit. Hyperbolic distance,
character orthogonality and Sobolev embedding are classical; Sol worker 1
developed the explicit flag-budget diagnostic with substantial OpenAI
assistant support. This is author verification, without independent
mathematical review, new Lean formalization or a stronger complete bound.
