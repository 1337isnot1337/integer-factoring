# Original smooth-set density in small-index subgroup cosets

## Theorem and remaining limits

Under the pinned uniform analytic and effective-estimate contract of the
[prime-form theorem](polylog-prime-form-generators.md), the original
recognizable smooth set has uniformly large density in every coset of a
subgroup of sufficiently small index.

Fix `0<kappa<1/64`. For every sufficiently large negative order discriminant
`Delta`, set `D=|Delta|`, `b=bits(D)`, `G=Pic(O_Delta)` and

```math
Y=\lceil\exp((\log b)^2)\rceil.
```

Let `S_Y` be the reduced primitive classes whose leading coefficient is
`Y`-smooth and coprime to `Delta`. For every subgroup `U<=G` of index
`d<=D^kappa` and every coset `tU`,

```math
\boxed{\frac{|S_Y\cap tU|}{|U|}
 >=\exp\left(-O_\kappa\left(\frac b{\log b}+\log^2 b\right)\right).}
                                                               \tag{1}
```

All constructed leading coefficients are below `sqrt(D)/4`, so they
satisfy the full reduced-height recognition predicate. Adding useful
proper-gcd or cofactor successes only increases density. For a covered
prime-power orbit, the [row-change lemma](relation-change-event-density.md)
also gives change-event density for its fully smooth rows at this same scale.

This proves a part of the original Path C density obligation under an
explicit subgroup-index condition. It does not bound the indices of the
actual target subgroups, give a short hit window, or establish the needed
event moments. The complete factoring exponent stays `1/2`.

The ingredients are the existing split-prime character estimate, an
elementary probabilistic partition, quotient Fourier mixing and the
small-norm reduction geometry. No novelty is claimed for the concentration
or finite-group mixing principles. This is not a positive Fourier
minorant of the smooth set; the
[ambient-minorant obstruction](uniform-ambient-minorant-obstruction.md)
remains valid.

## 1. A positive prime-band measure

Use the nonnegative band weight `psi` in
[smooth-band-density.md](smooth-band-density.md#1-uniform-split-primes-in-a-short-multiplicative-interval).
For each unramified split odd prime `p<Y`, let

```math
w_p=2\log(p)\psi(p/Y),\qquad W=\sum_p w_p,\qquad
w_{max}=\max_p w_p.
```

The two orientations of `g_p` receive half this weight each. The existing
uniform estimates, applied to the same band combination for every
character, give

```math
W\asymp Y,\qquad w_{max}<=\tfrac12\log Y,
```

```math
\left|\frac1W\sum_p w_p\operatorname{Re}\chi(g_p)\right|
 <=eta:=C bY^{-1/32}\quad(chi\ne1),                       \tag{2}
```

for one fixed sufficiently large constant `C`. The smaller prime-power
and conductor-removal errors are absorbed at this choice of `Y`.
In particular `eta->0`. The normalized measure is symmetric, so these
character averages are real.

Put

```math
R=\lfloor\sqrt D/4\rfloor,\qquad
L=\lfloor\log R/\log Y\rfloor.
```

We need `L` independent steps **with distinct underlying primes**.
Deleting repeated-prime tuples after mixing would lose too much compared
with the probability `1/d` when `d` is exponentially large. Partitioning
the primes before taking products avoids that loss.

## 2. A partition preserving all character bounds

Assign each underlying prime independently and uniformly to one of `L`
colors. For character `chi`, including the principal one, define

```math
T_{j,chi}=\sum_{p:\,color(p)=j}w_p\operatorname{Re}\chi(g_p).
```

Its expectation is the corresponding full weighted sum divided by `L`.
For fixed `(j,chi)`, the independent summands have range lengths at most
`w_p`. The elementary Hoeffding bound and `sum w_p^2<=w_max W` give

```math
\Pr\left(\left|T_{j,chi}-\mathbb ET_{j,chi}\right|
                        >eta W/L\right)
 <=2\exp\left(-\frac{2eta^2 W}{L^2w_{max}}\right).         \tag{3}
```

For completeness, the bound follows by writing the logarithmic moment
generating function of a bounded variable. Its second derivative is a
variance under an exponential tilt, at most one quarter of the squared
range. Integrating twice bounds the centered log moment by the squared
range times the parameter squared divided by eight. Multiplying these
bounds for independent summands and optimizing the exponential Markov
inequality gives (3), including both tails.

There are `L|G|` pairs. Since `|G|<=D`, `L<=b`, `W>=cY` and
`eta=C bY^(-1/32)`, the exponent in (3) is at least a constant times
`Y^(15/16)/log Y`, which eventually dominates `log(2L|G|)=O(b)`.
A union bound therefore has total probability less than one. Thus some
coloring satisfies all these inequalities simultaneously.

Fix such a coloring. Its color masses `W_j=T_(j,1)` satisfy
`(1-eta)W/L<=W_j<=(1+eta)W/L`. On each color normalize the original
positive weights and orientations to a probability measure `mu_j`.
Equations (2)--(3) imply

```math
|\widehat{mu_j}(chi)|<=\frac{2eta}{1-eta}<=4eta=:beta
\quad(chi\ne1),                                        \tag{4}
```

for sufficiently large `b`. Every color is nonempty. A single oriented
raw prime choice in color `j` has probability at most
`w_max/(2W_j)<=Lw_max/W`.

This is an existence proof for a density count. No algorithm is asked to
find the coloring or evaluate unknown class characters.

## 3. Mix on the quotient and count distinct actual classes

Choose one oriented prime independently from each color. Orthogonality
of the `d` characters trivial on `U` and (4) give, for every `tU`,

```math
\Pr(g_{p_1}^{s_1}\cdots g_{p_L}^{s_L}\in tU)
 >=1/d-beta^L>=1/(2d),                                  \tag{5}
```

provided `d beta^L<=1/2`. Translation changes the character terms by
unit-modulus phases and leaves this bound unchanged.

The primes in a raw tuple are distinct because their colors differ.
Their oriented ideal product is primitive, coprime to `Delta`, and has
squarefree norm

```math
a=\prod_jp_j<Y^L<=R.
```

Centering its middle coefficient `c` gives the form
`(a,c,(c^2+D)/(4a))`. Its trailing coefficient exceeds `a` since
`a<sqrt(D)/2`, so it is already reduced. Distinct raw tuples give distinct
classes: the underlying prime set determines the norm, the orientations
determine distinct roots modulo `2a` by CRT, and a primitive strict
reduced form has a unique representative with centered middle coefficient.
Coprimality excludes the boundary ambiguities. Fixing prime colors also
removes any permutation multiplicity.

Every raw tuple has probability at most `(Lw_max/W)^L`. Thus (5) and
injectivity prove the explicit counting inequality

```math
|S_Y\cap tU|>=\frac1{2d}\left(\frac W{Lw_{max}}\right)^L,
\qquad
\frac{|S_Y\cap tU|}{|U|}
 >=\frac1{2|G|}\left(\frac W{Lw_{max}}\right)^L.           \tag{6}
```

These are original smooth reduced classes, not quotients of smooth
classes and not unextracted represented-norm witnesses.

## 4. The asymptotic index and density estimates

From the definition of `L` and `beta=4C bY^(-1/32)`,

```math
L\log(1/beta)
 =\frac{\log D}{64}
   -O\left(\frac b{\log b}+\log^2 b\right).
```

Consequently `d<=D^kappa` with fixed `kappa<1/64` ensures
`d beta^L<=1/2` for all sufficiently large `D`, uniformly over `U`.
The index restriction is a condition of this method, not a density
obstruction for larger indices.

Use the existing elementary class-number upper bound
`|G|=O(sqrt(D)log D)`, as in the
[global smooth-density proof](smooth-set-fourier-norm.md#3-the-density-loses-only-a-subexponential-factor).
Since `W/w_max>=cY/log Y` and `L log Y>=log R-log Y`, taking logarithms
in (6) gives

```math
\log\frac{|S_Y\cap tU|}{|U|}
 >=-O\left(\log Y+\frac b{\log Y}\log b\right),
```

which is (1). The familiar global count is the special case `U=G`;
the quotient argument supplies the subgroup-coset extension.

## 5. Consequence for original prime orbits

For a target prime class `g_p` with `[G:<g_p>]<=D^kappa`, its original
smooth-hit density is at least `exp(-O(b/log b))`. The fixed row-multiplicity
threshold `K_p=min{k:p^k>H_red}` is `O(b)`. Therefore the existing uniform
row-change bridge gives event density at least

```math
kappa_p>=delta_p/(2K_p-1)>=exp(-O(b/log b)).
```

No raw-cluster-length conjecture is needed for this corollary. Actual
target-subgroup index bounds and quantitative short-window or moment
estimates remain open. Density by itself does not bound the first hit
in a short window, so this partial theorem is not a complete factoring
improvement or a supplied selection algorithm.

## Exact controls and review scope

```sh
make verify-small-index-smooth-density
```

The [checker](../tools/check_small_index_smooth_density.py) verifies a
quotient-mixing count on 3375 raw tuples, for every coset of four subgroups,
and compares 64 distinct prime-product classes at one fixed discriminant
against an independent literal CRT construction. These controls exercise
the quotient count and the critical injection below the reduction
threshold. The [receipt](../evidence/small-index-smooth-density.json)
records exact forms and rational probabilities. They do not infer the
uniform analytic estimate, probabilistic partition or asymptotic density
from finite data. The ordinary proofs and external contract are explicit;
no new Lean or independent mathematical review is claimed.
