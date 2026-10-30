# A verified deterministic factoring bound with exponent 1/30

Reference revision:
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).

## Result

The modified uniform factoring program satisfies

```math
\exists M\;\exists C>0\;\forall n\ge2,\qquad
T_M(n)\le C\exp\!\left(\frac{Cn}{(\log n)^{1/30}}\right).
```

The machine `M` has finite alphabets and computes `Nat.minFac` on every
natural-number input. It retains `Nat.minFac 0 = 2` and `Nat.minFac 1 = 1`.
The input and output use the reference's binary encoding. The variable `n`
denotes input length. Time denotes the reference's `TM2` machine-step bound.
The existential quantifier for `M` precedes every input-length quantifier.

The final module also proves the exact integer bound

```math
\exists M\;\exists c\in\mathbb N\;\forall n\in\mathbb N,\qquad
T_M(n)\le c\,2^{c(\lfloor n/(R(n)+1)\rfloor+1)}.
```

Here `R(n) = Nat.nthRoot 30 (Nat.log 2 n)`. Thus `R(0) = 0`, and for `n >= 1`,

```math
R(n)=\left\lfloor\sqrt[30]{\lfloor\log_2 n\rfloor}\right\rfloor.
```

The integer theorem covers the finite exceptions, including lengths below two.
The real theorem uses `n >= 2`, so its logarithm is positive.

The declarations appear in
[`Certified.lean`](../variants/k6/FactoringK6/Certified.lean):

- `Factoring.K6.factoring_root_rate`: the exact integer bound and finite alphabets.
- `Factoring.K6.factoring_real_rate`: the real-logarithm bound for the same machine.
- `Factoring.K6.certified_improvement`: the real bound and the strict exponent comparison.

These are closed theorems. They have no unproved grouping, runtime, or asymptotic
hypotheses. Their printed axioms are exactly `propext`, `Classical.choice`, and
`Quot.sound`. The [verification record](../evidence/k6-verification.json) identifies
the source hashes and reproduction command. The
[final axiom output](../evidence/k6-certified.log) records the three declarations.

## What improves

The reference's published theorem gives one machine with time `2^{o(n)}`.
The [baseline audit](baseline-audit.md) extracts the more explicit rate

```math
\beta=\frac{\log(12/11)}{\log16}=0.0313827205\ldots.
```

That extraction remains a handwritten derivation. The new machine's full
`1/30` bound and the inequality `beta < 1/30` are Lean theorems.
The comparison follows from the exact integer inequality

```math
12^{30}<16\cdot11^{30}.
```

A larger logarithmic saving exponent gives a smaller asymptotic upper-bound
form. Thus `1/30 = 0.033333...` improves the explicit reference bound.
This statement concerns upper bounds. It gives no lower bound on the reference
machine's actual runtime or the best possible analysis of that machine.

The source's square-root-logarithm conjecture asks for exponent `1/2`.
This result does not establish that conjecture. The stronger candidate near
`0.07551` in [integration.md](integration.md) also remains open.

## The actual program change

The reference chooses

```math
K_0=48,\qquad G_d=3K_d^8,\qquad K_{d+1}=9K_d^{16}.
```

The variant chooses

```math
K_0=48,\qquad G_d=5K_d^6,\qquad K_{d+1}=25K_d^{12}.
```

The stronger grouping lemma retains the reference's accuracy and label-width
conclusions. Its sufficient threshold is `g > K^2 + 4 D^2 K^4`.
The original accuracy uses `D = K`. For the relevant `K`, the choice
`g = 5 K^6` exceeds that threshold.

The [patch](../variants/k6/reference.patch) changes both mathematical parameters
and the runtime expression constructors. The top grouping and recursive grouping
use the new schedule. The corresponding specifications and cost proofs compile
against those changed constructors.

The crossing identity, counting reduction, target function, and machine model
retain their reference definitions. The compiler and simulator also retain
their reference source. The [source audit](k6-source-audit.md) records their
byte comparisons and identifies every modified source file.

Some old eighth-power expressions remain as conservative upper bounds.
They do not select the runtime schedule. The patched proof includes
`Gs_le_eighth` for those uses.

## Proof of the quantitative rate

### 1. A useful contraction at every positive height

The exact width formulas imply, for `K >= 600`,

```math
25N'\le23N+50G^4.
```

Here `N` and `N'` bound the label volume before and after one transition.
The statement gives contraction `23/25`, plus an additive width allowance.

The actual initial parameter remains `K_0 = 48`. Every positive height satisfies
`K_d >= 48^2 >= 600`. The proof uses the reference's nonexpansion bound for
the single bottom transition. For positive depth `ell >= 1`, the principal
bottom term is

```math
\left\lfloor\frac{23^{\ell-1}N_{\rm top}}{25^{\ell-1}}\right\rfloor.
```

The missing contraction changes a fixed constant. It does not change the
logarithmic saving exponent. `BottomException.lean` proves this assertion for
the actual base parameter, without a large replacement base.

### 2. Depth supplies a thirtieth-root saving

The runtime chooses its depth from the test `K_ell^40 <= n`.
The failed next test and the degree-twelve recurrence give

```math
\lfloor\log_2 n\rfloor\le5280\cdot12^\ell.
```

The exact certificate

```math
12\cdot23^{30}<25^{30}
```

then gives

```math
R(n)23^\ell\le2\cdot25^\ell,
\qquad (R(n)+1)23^\ell\le3\cdot25^\ell.
```

The proof handles the floors with natural-number arithmetic.
It does not substitute floating-point approximations for roots or logarithms.
The bottom exception adds a safe fixed factor of `25` to the weighted bound.

### 3. Every level and every auxiliary cost fits the bound

A small bottom table alone does not bound the algorithm.
`ChainEnvelope.lean` bounds all intermediate table exponents and the bottom
table in one expression. `BottomException.lean` specializes that expression
to the actual base parameter `48`.

`BudgetBounds.lean` bounds the top table, path enumeration, modulus sizes,
logarithmic terms, and depth allowances. These estimates use the same
parameters as the compiled program. `Quantitative.lean` combines them into

```math
\exists A\in\mathbb N\;\forall n,\qquad
Z_u(n)\le A\left(\left\lfloor\frac{n}{R(n)+1}\right\rfloor+1\right).
```

The positive-depth proof uses the runtime's size conditions.
A finite constant covers the zero-depth inputs. No asymptotic exception remains
outside the theorem's quantified bound.

### 4. The compiled machine satisfies the claimed time bound

`MachineBridge.lean` supplies a common value allowance `2^(12 * Zu(n) + 186)`.
The constant `186` covers the new recurrence coefficient through `25^40 <= 2^186`.
The source's uniform counting program retains its verified fixed polynomial
cost, including degree `2408`.

The proof carries that allowance through the prefix search, compiler, and
finite-alphabet simulation. Each polynomial cost contributes only a fixed
factor to the final exponent. The resulting theorem gives one machine with

```math
T_M(n)\le c\,2^{c(Z_u(n)+1)}.
```

Substitution of the bound on `Zu` proves `factoring_root_rate`.
The theorem does not choose a new machine for each input length or desired rate.

### 5. Lean verifies the conversion to real logarithms

`RateInterpretation.lean` proves, for `n >= 2`,

```math
(\log n)^{1/30}\le R(n)+1\le n.
```

It then absorbs the outer coefficient, the added one, and the base-two
exponential into a fixed positive real constant. This proves the displayed
real bound with the same machine. The logarithmic-power conversion is part
of the formal proof.

## Files and dependencies

| Module | Role |
| --- | --- |
| `RateArithmetic.lean` | Natural roots, depth bounds, and exact exponent certificates |
| `ChainEnvelope.lean` | Contraction and all-level table bounds |
| `BottomException.lean` | The actual base `48` and its single exceptional transition |
| `BudgetBounds.lean` | Top, path, modulus, logarithm, and depth costs |
| `Quantitative.lean` | A bound on the complete source exponent `Zu` |
| `MachineBridge.lean` | Value requirements and the compiled machine time bound |
| `RateInterpretation.lean` | Real logarithms, exponentials, and strict comparison |
| `Certified.lean` | Closed final theorems with no remaining hypotheses |

The proof modules live in [variants/k6/FactoringK6](../variants/k6/FactoringK6).
The patch applies to the pinned upstream source. The build script registers
the additional `FactoringK6` library in a fresh temporary source tree.

## Reproduction

The build requires Lean 4.35.0, Lake, Python 3.11 or newer, Git, Bash, and `flock`.
It also requires network access for the public reference and the Mathlib cache.
The pinned Mathlib revision is `56a005855ea17e0b0958f73ec801ae140b0a1679`.

From a checkout of this repository, run:

```sh
make verify
make verify-k6
```

For the older local lemmas and reference envelope, run:

```sh
make verify-formal
```

`make verify-k6` creates a fresh archive of the pinned reference.
It initializes an isolated Git repository and applies the patch there.
It verifies all 26 patched files and 43 protected files against their audited hashes.
Then it copies the repository's tracked proof sources.
It also verifies the target definitions against Git and the dependency revisions.
Then it builds `Solution` and all `FactoringK6` modules.
Finally, it runs `Certified.lean` again and verifies every final axiom report.

The script rejects `sorryAx`, failed proofs, and unapproved final axioms.
It records source hashes, build output, and final axiom reports under `build/k6/`.
Project proof outputs are isolated from the original reference build.
Only pinned dependency packages share a cache.

To retain the temporary source tree after a successful build, run:

```sh
K6_KEEP_BUILD=1 make verify-k6
```

The GitHub workflow includes a separate `k6-rate` job for this reproduction.
It pins the Lean release archive by SHA-256.

The [local verification record](../evidence/k6-verification.json) distinguishes
the completed local proof checks from this fresh CI build. It also records
the first archive attempt, which rejected an unapplied patch at the final theorem.
The isolated Git root and source-hash gate correct that packaging error.

## Scope and credit

The result verifies a stronger deterministic asymptotic upper bound for the
reference's factoring problem and machine model. The algorithm has very large
constants. This work contains no practical factoring benchmark.

Isaac Rudich and Louis-Martin Rousseau supplied the reference development.
Their project credits Claude for its generated ideas and Lean code.
The k6 extension used substantial OpenAI assistant support.
The files retain the Apache-2.0 license and identify adapted source.

Lean verification establishes the formal statements relative to the recorded
toolchain and standard axioms. Independent review of the construction and its
significance remains welcome. No leaderboard acceptance or publication review
is claimed.
