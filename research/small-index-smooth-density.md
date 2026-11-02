# Original smooth density on subgroups from even prime products

## Theorem and remaining limits

The original recognizable smooth set has uniformly large density on
subgroups of index below the quarter-power scale. A separate character
estimate extends the conclusion to every coset in a smaller index range.
Both are ordinary deductions from the identified analytic input, with
explicit limits.

Let `Delta` be a sufficiently large negative order discriminant, put
`D=|Delta|`, `b=bits(D)`, `G=Pic(O_Delta)` and

```math
Y=\lceil\exp((\log b)^2)\rceil.
```

Let `S_Y` consist of reduced primitive classes whose leading coefficient
is `Y`-smooth and coprime to `Delta`. Under the principal split-prime
band estimate in [smooth-band-density.md](smooth-band-density.md), for
every fixed `0<kappa<1/4` and every subgroup `U<=G` of index `d<=D^kappa`,

```math
\boxed{\frac{|S_Y\cap U|}{|U|}
 >=\exp\left(-O_\kappa\left(\frac b{\log b}+\log^2 b\right)\right).}
                                                               \tag{1}
```

Under the full pinned [prime-character estimate](polylog-prime-form-generators.md),
the same density bound holds in **every coset `vU`** for every fixed
`0<kappa<1/64` and `d<=D^kappa`. The subgroup conclusion uses only the
principal abundance estimate; the all-coset conclusion additionally
uses the nonprincipal estimate. The repository's pinned analytic
contract supplies both.

All constructed leading coefficients are below `sqrt(D)/4`, so they
satisfy the full reduced-height predicate. Adding proper-gcd or cofactor
successes only increases density. The [row-change lemma](relation-change-event-density.md)
then gives change-event density at the same scale for covered prime orbits.

Actual target-subgroup indices, hit gaps and event moments remain open.
This partial theorem does not improve the complete factoring exponent
beyond `1/2`. It constructs no large ambient Fourier minorant and is
consistent with the [ambient obstruction](uniform-ambient-minorant-obstruction.md).
The concentration and finite-group Fourier principles are classical;
the arithmetic deduction and its stated scope are what is recorded here.

## 1. A symmetric prime-band measure and an even norm budget

Use the nonnegative band kernel `psi` from
[smooth-band-density.md](smooth-band-density.md#1-uniform-split-primes-in-a-short-multiplicative-interval).
For each unramified split odd prime `p<Y`, define

```math
w_p=2\log(p)\psi(p/Y),\quad W=\sum_p w_p,\quad
w_{max}=\max_p w_p.
```

The principal estimate gives `W asymp Y`, and
`w_max<=log(Y)/2`. Give the two orientations half the weight each.
For every class character, its normalized average is real:

```math
a_\chi=\frac1W\sum_p w_p\operatorname{Re}\chi(g_p)\in[-1,1].
```

No smallness assumption on `a_chi` is needed for (1). The full
nonprincipal estimate, used only for the all-coset result, gives

```math
|a_\chi|<=eta:=C bY^{-1/32}\quad(chi\ne1)                 \tag{2}
```

for one fixed sufficiently large constant `C`.

Set

```math
R=\lfloor\sqrt D/4\rfloor,\qquad
L=2\left\lfloor\frac{\log R}{2\log Y}\right\rfloor.
```

Thus `L` is even, `Y^L<=R` and `L log Y=log R-O(log Y)`.
We will choose one prime from each of `L` disjoint colors, avoiding
repeated-prime cancellations while keeping the same norm budget.

## 2. A simultaneous high-precision color partition

Assign each underlying prime independently and uniformly to one of `L`
colors. For every character, including the principal one, set

```math
T_{j,chi}=\sum_{p:\,color(p)=j}w_p\operatorname{Re}\chi(g_p).
```

Its expectation is `Wa_chi/L`. Put `t=b^2Y^(-1/2)`, which tends to zero.
For fixed `(j,chi)`, the independent summands have range lengths at most
`w_p`, so the elementary Hoeffding bound gives

```math
\Pr\left(|T_{j,chi}-Wa_\chi/L|>tW/L\right)
 <=2\exp\left(-\frac{2t^2W}{L^2w_{max}}\right).           \tag{3}
```

One proof bounds the centered logarithmic moment-generating function.
Its second derivative is a variance under an exponential tilt, at most
one quarter of the squared range. Integrating twice, multiplying over
independent summands and optimizing exponential Markov gives (3), with
both tails. We used `sum_p w_p^2<=w_max W`.

There are `L|G|` pairs. Since `|G|<=D`, `L<=b` and `W>=cY`, the exponent
in (3) is at least a constant times `b^2/log Y`, which eventually exceeds
`log(2L|G|)=O(b)`. The union bound is less than one. Fix a coloring for
which all deviations are at most `tW/L`.

Its color masses satisfy `(1-t)W/L<=W_j<=(1+t)W/L`. Normalize the
original positive weights and orientations in each color to a measure
`mu_j`. For every character,

```math
|\widehat{mu_j}(chi)-a_\chi|<=\frac{2t}{1-t}<=4t.         \tag{4}
```

Indeed, both its weighted numerator and its mass differ from their
expectations by at most `tW/L`, and `|a_chi|<=1`. Every color is nonempty.
A single oriented raw prime choice has probability at most
`w_max/(2W_j)<=Lw_max/W`. This is an existence argument for a count;
an algorithm need not compute the coloring or the class characters.

## 3. Even products give a one-sided Fourier bound

Let `mu=mu_1*...*mu_L`. All factors in each Fourier coefficient are real.
If `|a_chi|>4t`, equation (4) makes them all have the same sign. Their
product is nonnegative because `L` is even. Otherwise every factor has
absolute value at most `8t`. Hence

```math
\widehat\mu(chi)>=-\gamma,\qquad\gamma=(8t)^L,
                                                               \tag{5}
```

for **every** character, without using (2).

For a subgroup of index `d`, orthogonality of its `d` trivial-on-`U`
characters has phase one at the subgroup itself. Therefore

```math
\Pr_\mu(g\in U)>=1/d-\gamma>=1/(2d)
\quad\hbox{if }d\gamma<=1/2.                            \tag{6}
```

For an arbitrary coset, its phases can change the positive terms' signs.
Here use (2): since `t/eta->0`, equation (4) gives
`|widehat{mu_j}(chi)|<=2eta=:beta` for nonprincipal characters and large
`b`. Thus the separate bound is

```math
\Pr_\mu(g\in vU)>=1/d-\beta^L>=1/(2d)
\quad\hbox{if }d\beta^L<=1/2.                           \tag{7}
```

The distinction between (6) and (7) is essential. The quarter-power
statement concerns the subgroup itself.

## 4. Count distinct original smooth reduced classes

Choose the oriented primes independently by their colors. Their
underlying primes are distinct. Their ideal product is primitive,
coprime to `Delta`, and has squarefree norm

```math
a=\prod_jp_j<Y^L<=R.
```

Center its middle coefficient `c`. The form
`(a,c,(c^2+D)/(4a))` is already reduced, since `a<sqrt(D)/2` forces its
trailing coefficient above `a`. Distinct raw tuples give distinct
classes: the prime set determines the norm, orientations determine
distinct roots modulo `2a` by CRT, and a primitive strict reduced form
has a unique centered representative. Coprimality excludes boundary
ambiguities. Fixed prime colors remove permutation multiplicities.

Every raw tuple has probability at most `(Lw_max/W)^L`. Applying (6),
or (7) for a permitted coset, gives

```math
|S_Y\cap vU|>=\frac1{2d}\left(\frac W{Lw_{max}}\right)^L,
\qquad
\frac{|S_Y\cap vU|}{|U|}
 >=\frac1{2|G|}\left(\frac W{Lw_{max}}\right)^L.           \tag{8}
```

In the subgroup case take `v=1`. These are recognizable smooth reduced
leading coefficients, rather than quotients or unextracted norm witnesses.

## 5. The two index ranges and the common density scale

For `t=b^2Y^(-1/2)` and even `L` as above,

```math
\log\gamma=-\frac{\log D}{4}
              +O\left(\frac b{\log b}+\log^2b\right).
```

Thus `d<=D^kappa`, fixed `kappa<1/4`, ensures the subgroup condition in
(6), uniformly for sufficiently large `D`. For the separate all-coset
condition, (2) gives

```math
\log(\beta^L)=-\frac{\log D}{64}
              +O\left(\frac b{\log b}+\log^2b\right),
```

so fixed `kappa<1/64` ensures (7). Neither threshold asserts that density
fails beyond the range of this proof.

Use the existing elementary upper bound `|G|=O(sqrt(D)log D)`, as in the
[global density count](smooth-set-fourier-norm.md#3-the-density-loses-only-a-subexponential-factor).
From `W/w_max>=cY/log Y` and `L log Y>=log R-2log Y`, equation (8) gives

```math
\log\frac{|S_Y\cap vU|}{|U|}
 >=-O\left(\log Y+\frac b{\log Y}\log b\right).
```

This proves (1) and its smaller-range all-coset extension. The global
count is the special case `U=G`; the even-product argument is what
extends the subgroup range beyond absolute character mixing.

## 6. Consequence for original prime orbits

For a target prime class `g_p` with `[G:<g_p>]<=D^kappa`, fixed
`kappa<1/4`, its original smooth-hit density is at least
`exp(-O(b/log b))`. The canonical-row threshold
`K_p=min{k:p^k>H_red}` is `O(b)`. The existing row-change bridge gives

```math
kappa_p>=delta_p/(2K_p-1)>=exp(-O(b/log b)).
```

No raw-cluster-length conjecture is needed. Actual target-subgroup
indices and quantitative hit-window or moment bounds remain open.
Density alone does not supply a short first-hit window, so no complete
factoring improvement or new selection algorithm follows here.

## Exact controls and review scope

```sh
make verify-small-index-smooth-density
```

The [checker](../tools/check_small_index_smooth_density.py) checks 50625
raw tuples for every coset of four finite subgroups. Its four color
distributions have eight genuinely negative product Fourier coefficients;
the exact one-sided lower bound and the subgroup count still hold.
Separately, 64 distinct arithmetic classes at one fixed discriminant
are compared with a literal CRT construction. These three-color geometric
controls check injectivity, whose proof does not require parity.
The [receipt](../evidence/small-index-smooth-density.json) records exact
forms and rational quantities. Finite data do not establish the analytic
input, simultaneous color partition or asymptotic theorem. Those
ordinary proofs and external assumptions are stated explicitly;
no new Lean or independent mathematical review is claimed.
