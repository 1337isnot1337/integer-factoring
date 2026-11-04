# Every large class group obstructs the ambient-minorant program

## Result and scope

Under the existing pinned analytic and effective-estimate contract,
the [ambient positive-minorant obstruction](ambient-minorant-norm-gap.md)
can be made **uniform over every sufficiently large quadratic-order
discriminant**. No special semiprime family or least-prime progression
construction is needed.

Let `D=|Delta|`, `b=bits(D)`, and let the smoothness cutoff `Y` and residual
cofactor cap `B` satisfy `Y,B=2^o(b)`. Use the full reduced-height
smooth/cofactor/proper-gcd success predicate. Every real function on the
proper class group with `f(1)=1`, nonnegative Fourier coefficients and
`f<=0` outside that success set has ambient mean

```math
\boxed{tau<=2^{-b/64+o(b)}.}                              \tag{1}
```

The result also holds if ambiguous-class early successes are added.
Consequently changing a permitted discriminant multiplier does not repair
this ambient-mean certificate when the cutoffs remain subexponential in
the resulting bit length. For quasipolynomial smoothness cutoff, the
original smooth set has density `2^-o(b)` and the one-sided negative
Fourier bound must be at least `2^{b/64-o(b)}`.

The uniformly large ambient-mean program is stopped. The exact
near-identity hitting lemma remains valid. Original target-cycle density,
power-orbit moments, subgroup-specific means and efficient smooth-pair
extraction remain separate open obligations. No factoring lower bound,
new Lean proof or improved complete exponent is claimed.

## 1. Analytic inputs already in use

Use the [uniform prime character estimate](polylog-prime-form-generators.md)
and the nonnegative band kernel from
[smooth-band-density.md](smooth-band-density.md#1-uniform-split-primes-in-a-short-multiplicative-interval):

```math
\psi(t)=w(t)-\frac{21}{16}w(4t/3)+\frac{27}{64}w(16t/9),
\qquad w(t)=t(1-t)\mathbf1_{[0,1]}(t).
```

It is supported in `[9/16,1]`, nonnegative, and has integral `259/6144`.
For every class character, apply this same linear combination to the
prime estimates at `x,3x/4,9x/16`. The principal main term is
`(259/6144)x`; every nonprincipal sum is

```math
O\left(bx^{31/32}+\sqrt{x}\log x+b\log x\right).           \tag{2}
```

All constants are uniform in the order and character under the identified
contract. The sum of the absolute combination coefficients is
`175/64<3`, so the original uniform error envelope changes only by an
absolute constant. Normalizing the positive band weights therefore gives
nonprincipal character bias

```math
epsilon(x)=O\left(bx^{-1/32}+x^{-1/2}\log x+
                              bx^{-1}\log x\right).      \tag{3}
```

The existing principal band estimate also supplies an unramified split
odd prime in `(9x/16,x)` whenever
`x>=T=(10^6(b+1))^32`. This explicit threshold is inherited from that
note; the finite controls below do not establish it. The argument does
not need to compute the probability measure or enumerate its support.

## 2. A missed translate for an arbitrary discriminant

Set

```math
Z=\max(Y,B,T),\qquad X=\left\lfloor\frac{\sqrt D}{8Z}\right\rfloor.
                                                               \tag{4}
```

Apply the split-prime band estimate at `2Z` to choose an unramified split
prime `A` with `9Z/8<A<2Z`. Let `Q` be either prime-form class of norm `A`.
For sufficiently large `D`, subexponential `Z` ensures `X>=T` and
`9X/16>2Z`. Take the positive prime sample at cutoff `X` from Section 1.
Every supported prime `p` is distinct from `A`, coprime to `Delta`, and

```math
Ap<2ZX<=\sqrt D/4.                                      \tag{5}
```

The product of `Q` with either orientation of a supported prime form has
primitive norm `Ap`. Centering its middle coefficient `c` gives

```math
\left(Ap,c,\frac{c^2+D}{4Ap}\right),\qquad |c|<=Ap.
```

Its trailing coefficient exceeds `Ap` by (5). It is already reduced, so
its actual reduced leading coefficient is **exactly `Ap`**. This is the
only norm-gap fact needed; the middle coefficient of `Q` can be arbitrary.

The coefficient `Ap` lies below the full height `floor(sqrt(D/3))`, is
coprime to `Delta` and hence to any input `N` dividing `Delta`, and contains
the prime `A>Z>=Y,B`. Stripping every base prime at most `Y` leaves a
cofactor at least `A>B`. Thus every supported translate fails precisely
the smooth/cofactor requirement and yields no proper gcd.

It is not an ambiguous class either. An ambiguous reduced form has middle
coefficient zero or of absolute value its leading coefficient, or equal
leading and trailing coefficients. The last case is excluded by (5).
Either middle-coefficient case would force a prime divisor of `Ap` to
ramify, contradicting coprimality with `Delta`. Adding ambiguous-class
early successes therefore does not change this missed translate.

## 3. Fourier consequence

The [missed-translate lemma](ambient-minorant-norm-gap.md#1-a-missed-low-bias-translate-bounds-every-positive-minorant)
gives `tau<=epsilon(X)/(1+epsilon(X))` for every admissible minorant.
Since `Z=2^o(b)`, (4) gives

```math
\log_2X=b/2-o(b),\qquad epsilon(X)=2^{-b/64+o(b)},
```

proving (1). More explicitly, the leading term is
`O(b D^(-1/64) Z^(1/32))`; the other two terms in (3) are asymptotically
smaller on this scale. No numerical Fourier optimization or new uniform
smooth-power theorem enters the argument.

For `Y=ceil(exp((log b)^2))`, the existing distinct-split-prime product
count gives global smooth density
`delta>=exp(-O(log^2 b+b/log b))=2^-o(b)`, uniformly in `Delta`.
The full success set is inversion symmetric and contains the identity.
If `L` bounds the negative unnormalized character sums, the
identity-corrected minorant has mean at least `delta/(L+1)`. Applying (1)
forces

```math
L+1>=2^{b/64-o(b)}.                                     \tag{6}
```

This one-sided conclusion uses positivity and the missed translate;
it does not rely on an absolute-value Parseval lower bound.

## 4. Scope of possible continuations

An ambient mean is the principal Fourier coefficient over the full class
group. A mean on one target subgroup includes every ambient character
trivial on that subgroup, so it is a different quantity. In particular,
(1) does not establish a rough interval on the target power orbit.

There is a precise limitation even for a subgroup certificate. Extend a
positive-definite function on a subgroup `H` by zero outside `H`.
Its ambient Fourier coefficients are the nonnegative subgroup coefficients
divided by `[G:H]`. Hence its ambient mean is `tau_H/[G:H]`, and (1) implies

```math
tau_H<=[G:H]2^{-b/64+o(b)}.
```

Thus a subgroup of subexponential index cannot escape this particular
minorant obstruction. No bound on the index of the actual prime-generated
cyclic subgroup is supplied here. Subgroup means and the original moment
or pooled-rank questions require their own arithmetic arguments.

## Exact controls and review scope

```sh
make verify-uniform-ambient-minorant-gap
```

The [checker](../tools/check_uniform_ambient_minorant_gap.py) verifies the
band polynomial pieces and integral exactly, and checks 20 oriented prime
compositions at five fixed fundamental/nonmaximal order discriminants.
The middle coefficient is unrestricted. Every product stays reduced,
fails the full-height cofactor predicate and is nonambiguous. Its
[receipt](../evidence/uniform-ambient-minorant-gap.json) records exact forms.
These small geometric cutoffs are not claimed to meet the analytic
threshold `T`; uniform existence and bias use the stated external input.
The earlier semiprime and sharp finite Fourier controls remain intact.

The [clean full default replay](../evidence/uniform-ambient-minorant-full-verification.json)
passes at exact source `bfb5f74b487b725b4e3ccf9ab2819b3f94151e98`
in 773.91 seconds. The archive follow-up changes review prose and receipts.
