# An arithmetic obstruction to uniformly large ambient positive minorants

## Result and scope

The [near-identity criterion](near-identity-positive-minorants.md) is a valid
sufficient hitting lemma, but its proposed uniformly subexponential
**ambient mean** cannot hold for all inputs. Under the same pinned analytic
input as the [prime-form sample theorem](polylog-prime-form-generators.md),
there is an infinite family of semiprimes for which every admissible
ambient positive minorant has mean at most

```math
tau=O\left(b\,2^{-b/5376}\right).
```

These semiprimes survive the stated trial division through
`2^ceil(sqrt(b))`, primality and perfect-power preprocessing. The success
set is the **full reduced-height smooth-or-proper-gcd set**, rather than a
smaller radius restriction. For a quasipolynomial smoothness cutoff, its
global density is still `2^-o(b)` under that analytic input. The associated
one-sided negative-spectrum repair must also pay an exponential cost.

This stops the proposed uniform ambient-minorant and uniform one-sided
negative-spectrum routes. It is not a factoring lower bound or an
obstruction to original target-cycle density, power-orbit moments,
subgroup-specific minorants, a different multiplier choice, or smooth-pair
extraction. The complete factoring exponent stays `1/2`.

The mechanism combines the existing
[arithmetic missed-translate calculation](prime-translates-miss-smooth-forms.md)
with Fourier positivity. The earlier
[Fourier algebra norm obstruction](smooth-set-fourier-norm.md) used the
same analytic prime estimate with an absolute Fourier bound. The argument
here controls **every** admissible positive minorant and also applies when
the function may be negative outside the success set.

## 1. A missed low-bias translate bounds every positive minorant

Let `G` be finite abelian, `S` a subset, and `mu` a probability measure
with nonprincipal absolute character bias at most `epsilon`. Suppose
`t supp(mu)` is disjoint from `S`. Let a real function `f` satisfy

```math
f(1)=1,\qquad \widehat f(chi)>=0,\qquad f(x)<=0\quad(x\notin S),
```

using Fourier coefficients normalized by `|G|`. Write
`tau=widehat f(1)`. Fourier inversion at one gives
`sum_chi widehat f(chi)=1`, so

```math
0>=\mathbb E_{u\sim mu}f(tu)
 >=tau-epsilon(1-tau).
```

The first inequality uses the missed translate. For the second, each
nonprincipal character expectation has absolute value at most `epsilon`,
translation has unit-modulus phase, and all Fourier coefficients are
nonnegative. Consequently

```math
\boxed{tau<=epsilon/(1+epsilon).}                         \tag{1}
```

This bound is sharp: on a group of order `h`, take `S={1}`, let `mu` be
uniform on the other elements, and use `f=1_{\{1\}}`. Then
`epsilon=1/(h-1)` and `tau=1/h`.

## 2. The full arithmetic success set still has a missed translate

Take an odd prime `A`, and let `P,R` be distinct primes with

```math
P=1\pmod{4A^2},\qquad R=-1\pmod{4A^2P^2}.
```

Set `N=PR`, `Delta=-N`, `C=(N+1)/(4A)` and `Q=(A,1,C)`.
Then `N=3 mod 4`, `Delta` is a fundamental negative discriminant, and
`Q` is primitive and reduced. In fact,
`P>=4A^2+1` and `R>=4A^2P^2-1`, so `C>A^3` and `N>4A^4`.

For every unramified split odd prime `p<A`, compose `Q` with either
orientation of its prime form. Since `p` and `A` are distinct and
coprime to `N`, the product has primitive norm `Ap`. Center its middle
coefficient `B` so `|B|<=Ap`. Its trailing coefficient is

```math
(B^2+N)/(4Ap)>Ap,
```

because `N>4A^4>4(Ap)^2`. Hence it is already reduced and its leading
coefficient is exactly `Ap`.

Let `H=floor(sqrt(N/3))`, and let `S` contain every class whose reduced
leading coefficient is `Y`-smooth, together with every class whose leading
coefficient has a proper gcd with `N`. For any `Y<A`, the translated
coefficient `Ap` is at most `H`, coprime to `N` and divisible by the prime
`A>Y`. It fails precisely smoothness. Thus **every** such prime translate
misses the full success set `S`.

Use the real weighted prime distribution of the
[existing prime estimate](polylog-prime-form-generators.md), at cutoff
`X=A-1`. Under its pinned analytic and effective-estimate contract, all
nonprincipal class-character biases satisfy

```math
epsilon_A=O\left(X^{-1/32}\log(N+3)
 +X^{-1/2}\log X+X^{-1}\log(N+3)\log X\right).
                                                               \tag{2}
```

This is a mathematical probability measure, not an algorithm required to
enumerate exponentially many primes. Its entire support has the missed
translate above. No finite character computation substitutes for (2).

## 3. Polynomial-size semiprimes and the preprocessing check

Choose `P` and `R` to be the least primes in their displayed progressions.
Only a coarse form of Linnik's theorem is needed: the least prime in any
coprime progression of modulus `q` is `O(q^6)`, with an absolute constant.
For a primary source, Xylouris,
[*On Linnik's constant*, arXiv:0906.2749](https://arxiv.org/abs/0906.2749),
states the stronger admissible exponent 5.2. No claim about the latest
record or a new least-prime proof is made here.

It follows successively that

```math
P=O(A^{12}),\qquad R=O(A^{156}),\qquad N=O(A^{168}).       \tag{3}
```

Thus `b=bits(N)<=168 log_2 A+O(1)`, while `log N=Theta(log A)`.
The smaller factor `P>=4A^2+1` eventually exceeds
`2^ceil(sqrt(b))`. The two distinct prime factors make `N` neither prime
nor a perfect power. The stated preprocessing therefore leaves this
family. Also `R/P` tends to infinity, so the example does not depend on
the trivial `4A^4-1` difference-of-squares family used earlier.

Equations (1)--(3) give

```math
tau<=O\left(b A^{-1/32}\right)
    <=O\left(b\,2^{-b/5376}\right).                     \tag{4}
```

The examples are a mathematical obstruction to a proposed certificate.
Known factors used to specify them do not imply that they are hard inputs,
and no hardness assumption is needed for (4). The statement concerns the
discriminant `-N`; an algorithm allowed another multiplier requires its
own analysis.

## 4. Large global density and the negative-spectrum consequence

For a concrete cutoff, use `Y=ceil(exp((log b)^2))`, which is smaller than
`A` eventually along (3). The existing
[smooth-density counting proof](smooth-set-fourier-norm.md#3-the-density-loses-only-a-subexponential-factor)
applies here with `R_0=floor(sqrt(N)/2)` and
`k=floor(log(R_0)/log(Y))`. The pinned prime estimate supplies at least
`M>=cY/log Y` eligible split primes below `Y`. Products of `k` distinct
such primes give distinct reduced leading coefficients at most `R_0`;
`4R_0^2<N` guarantees reduction. With
`|G|=O(sqrt(N) log N)`, the binomial count gives

```math
delta=|S|/|G|>=exp\left(-O\left(\log Y+
              \frac{b\log b}{\log Y}\right)\right)=2^{-o(b)}.
                                                               \tag{5}
```

Adding the proper-gcd success classes only increases density, and their
addition did not repair the missed translate. Thus this is an actual
arithmetic large-density example, not just the abstract density-only
counterexample in the earlier minorant note.

Let `B>=0` bound the negative part of the unnormalized real sums
`F_S(chi)`, so `F_S(chi)>=-B`. Inversion symmetry of the reduced leading
coefficient makes these sums real. Since the identity is smooth, the
standard identity-corrected function

```math
f=(1_S+B\,1_{\{1\}})/(B+1)
```

is an admissible minorant with mean
`tau=(|S|+B)/(|G|(B+1))>=delta/(B+1)`. Applying (1) therefore forces

```math
B+1>=delta(1+epsilon_A)/epsilon_A
    >=2^{b/5376-o(b)}.                                  \tag{6}
```

This obstructs a uniformly subexponential **one-sided negative** bound,
without using the absolute-value Parseval obstruction. Large positive
Fourier coefficients cannot remove (6).

## 5. What can still matter for Path C

The obstruction uses an ambient mean and an arbitrary translated prime
sample. It does not supply a long rough interval on a target prime-power
orbit. The actual selector starts near the identity on that orbit.
A positive function on the target subgroup, or a global positive function
with a large **subgroup** mean, is a different possible certificate. The
mean over that subgroup includes every ambient character trivial on it,
not just the ambient principal character. Equation (1) does not bound
that subgroup mean unless its own hypotheses hold inside the subgroup.

Uniform original target-cycle density, event-moment estimates and the
pooled-rank approach remain open. A continuation should address those
arithmetic obligations or a subgroup-specific certificate, rather than
optimize the now-excluded uniform ambient minorant.

## Exact controls and review scope

```sh
make verify-ambient-minorant-gap
```

The [checker](../tools/check_ambient_minorant_gap.py) verifies the sharp
finite Fourier lemma on `(Z/2)^4`, including a function negative outside
the set, and checks 26 oriented prime translates for six actual
semiprimes against the full-height recognition predicate. It verifies
the factor primality by trial division and the finite preprocessing
inequalities. The [receipt](../evidence/ambient-minorant-gap.json) records
the exact forms. These controls do not prove Linnik's theorem, the uniform
analytic estimate, the infinite family, or independent mathematical review;
those ingredients and the ordinary deduction are identified above.
