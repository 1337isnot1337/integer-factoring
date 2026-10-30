# The verified degree-twelve factoring bound

## Statement

Let `n` be the binary length of the integer being factored. The closed theorem
`Factoring.K6.factoring_real_rate` asserts that there are one deterministic
finite-alphabet machine `M` and a constant `C > 0` such that

```math
T_M(n)\le C\exp\!\left(\frac{Cn}{(\log n)^{1/30}}\right)
\qquad\text{for every }n\ge2.
\tag{1}
```

Correctness on every input is part of `M : FactoringMachine`. The type is
the reference's `TM2ComputableInTime encodeNat encodeNat Nat.minFac`.
Its `FiniteAlphabets` condition applies to every stack alphabet. The actual
machine steps, including compilation and arithmetic costs, are bounded.

In terms of the integer `N`, the asymptotic bound is equivalently

```math
\exp\!\left(O\!\left(
\frac{\log N}{(\log\log N)^{1/30}}
\right)\right).
\tag{2}
```

For the exact integer statement, set
`r(n) = Nat.nthRoot 30 (Nat.log 2 n)`. These are natural-number operations;
in particular `Nat.log 2 0 = 0`. The theorem supplies one finite-alphabet
`M : FactoringMachine` and one natural constant `C`, valid at every natural
input length:

```math
\exists M\;\exists C\in\mathbb N\;\forall n\in\mathbb N,\quad
T_M(n)\le C\,2^{C\left(
\left\lfloor\frac{n}{r(n)+1}\right\rfloor+1
\right)}.
\tag{3}
```

Lean checks (1) and (3), as well as their connection, in
[`Certified.lean`](../variants/k6/FactoringK6/Certified.lean).
The [variant instructions](../variants/k6/README.md) state the current build
status and reproduction command. The [source audit](k6-proof-audit.md) checks
the theorem's meaning and the changed program.

## Comparison target

The pinned reference is Isaac Rudich and Louis-Martin Rousseau's
[`SubexponentialDeterministicFactoring@c91cc3f`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
Its published endpoint is one machine satisfying `N^{o(1)}`.

The official BigO
[`integer-factoring/deterministic-bit-model/v1` data](https://github.com/bigospeedrun/site/blob/5f428921611f2f0b740e4553a1c2cf92cd6165b2/data/integer-factoring.json)
assigns that reference score 0. A positive score `alpha` requires, for every
fixed `epsilon > 0`, eventual time at most

```math
2^{n/(\log n)^{\alpha-\epsilon}}.
```

Equation (1) meets this criterion for `alpha = 1/30`. To include the fixed
prefactor explicitly, put `a = 1/30` and take any fixed `epsilon > 0`. Then

```math
\frac{\log_2\!\left(C\exp(Cn/(\ln n)^a)\right)}
     {n/(\ln n)^{a-\epsilon}}
=\frac{\ln C}{\ln2}\frac{(\ln n)^{a-\epsilon}}{n}
 +\frac{C}{\ln2\,(\ln n)^\epsilon}
\longrightarrow0.
```

The ratio is therefore at most 1 for all sufficiently large `n`, which gives
the required bound. Multiplying by a fixed polynomial in `n` only adds a
constant times `log n` to the numerator, and also leaves the limit zero.
Changing the logarithm base has the same harmless constant-factor effect.

The machine and the constant in (1) are chosen once, before `epsilon`; only
the eventual input threshold depends on `epsilon`. This meets the official
requirement of one algorithm for every input. This is a mathematical score
comparison; no submission or acceptance by the site is asserted.

The [baseline analysis](baseline-audit.md) extracts the more explicit exponent

```math
\alpha_{\rm ref}=\frac{\log(12/11)}{\log16}
=0.0313827205209647\ldots.
```

The new exponent is also strictly larger than this number. Lean proves

```math
\frac{\log(12/11)}{\log16}<\frac1{30},
```

using the exact integer inequality `12^30 < 16 * 11^30`. Taking logarithms
gives `30 log(12/11) < log 16`, and `log 16 > 0` permits the division.
For any fixed positive constants, the new main exponent is asymptotically
smaller, because

```math
\frac{n/(\log n)^{1/30}}{n/(\log n)^{\alpha_{\rm ref}}}
=(\log n)^{\alpha_{\rm ref}-1/30}\longrightarrow0.
```

The closed comparison theorem combines `alpha_ref < 1/30` with the new
machine bound. The unchanged reference's rate at `alpha_ref` is the
separate handwritten baseline derivation; it is not an additional machine
theorem contained in `certified_improvement`. The two programs use the same
factoring specification and machine model, but the new program has the
changed schedule described below. This is a comparison of proved upper-bound
expressions. It does not prove a lower bound on the reference program's
runtime, or that the new program is faster at each input.

## Change to the uniform algorithm

The original recursion uses

```math
G_d=3K_d^8,\qquad K_{d+1}=G_d^2=9K_d^{16}.
```

The existing second-moment argument admits a stronger grouping threshold:

```math
g>K^2+4D^2K^4.
\tag{4}
```

The original accuracy parameter is `D = K`. For every `K >= 2`,

```math
K^2+4K^6<5K^6.
```

Therefore the same grouping guarantees and value-preservation arguments work
with

```math
K_0=48,\qquad G_d=5K_d^6,\qquad K_{d+1}=25K_d^{12}.
\tag{5}
```

The patch changes the actual program expressions and their cost proofs.
The program still searches deterministically for the first valid grouping.
The use of a counting argument to prove existence adds no random choices
to the executable algorithm. The original crossing factors, modular zero
test, prefix search, output encoding, and TM simulator are retained.

## Why the logarithmic exponent improves

The input length determines a maximal admissible depth `ell`, with the
positive-depth cutoff

```math
K_{\ell}^{40}\le n<K_{\ell+1}^{40}.
\tag{6}
```

The degree-twelve schedule satisfies

```math
2K_d\le(2K_0)^{12^d}.
```

Thus `log_2 n` is bounded by a constant times `12^ell`. The smaller growth
degree permits more useful recursion levels at a fixed input length.

The reference's exact volume inequality is

```math
12K N_{\rm next}\le(11K+24)N+24K G^4.
\tag{7}
```

For `K >= 600`, the coefficient of `N` is at most `23/25`, since

```math
\frac{11}{12}+\frac2K\le\frac{23}{25}.
```

The machine still starts at `K_0 = 48`. The proof explicitly handles its
single exceptional bottom transition using nonexpansion. Every positive
schedule height has at least 600 groups. Losing that one contraction changes
only a fixed constant.

The exact certificate

```math
12\cdot23^{30}<25^{30}
\tag{8}
```

turns (6) and the contraction into the checked inequality

```math
\bigl(r(n)+1\bigr)23^{\ell}
\le3\cdot25^{\ell}.
\tag{9}
```

The proof bounds the top construction, maximum exponent over every recursion
level, path counts, modular budgets, additive padding, and bookkeeping with
the same weight. For the actual uniform engine exponent `Zu(n)`, it obtains

```math
Zu(n)\le A\left(
\left\lfloor\frac{n}{r(n)+1}\right\rfloor+1
\right)
\tag{10}
```

for a fixed natural `A`. The finite depth-zero range is bounded explicitly by
a finite sum; it is not omitted from the theorem.

## Full machine cost

`MachineBridge.lean` reconstructs the actual counting specification at value
bound

```math
V(n)=2^{12Zu(n)+186}.
```

The number 186 accounts for `25^40 <= 2^186`. The source command cost is a
fixed polynomial in `V(n)+1`; the checked compiler and simulator add fixed
polynomial overhead. Their composition therefore stays within

```math
C\,2^{C(Zu(n)+1)}.
```

Combining this with (10) proves (3). The real conversion proves, for
`n >= 2`, that its integer denominator `D` satisfies

```math
(\log n)^{1/30}\le D\le n,
\qquad
\lfloor n/D\rfloor+1\le\frac{2n}{(\log n)^{1/30}}.
```

This yields (1) with all rounding and prefactors included.

## Verification and scope

Run `make verify-k6`. The gate checks the exact public reference revision,
applies the source patch in its own Git root, checks protected definitions,
and rebuilds both the modified `Solution` and all new proof modules.
It then runs the endpoint source again and checks its axiom reports.

The final statements require only `propext`, `Classical.choice`, and
`Quot.sound`. The reference's intentionally unfinished challenge statement
is preserved as a specification template; it is not a proof dependency of
the certified result.

The verified bound has enormous constants. It establishes no practical
factoring benchmark or usable crossover point. Exponent `1/30` is smaller
than the reference's conjectured `1/2`; the
[Appendix B integration route](appendix-b-factoring-route.md) records a
separate possible approach to that stronger bound, with its formal adapters
and machine composition still open.
