# A route from Appendix B to the square-root logarithm factoring bound

## Status

The fully assembled clique theorem in `anthropics/formal-math/eth` already
provides the strong recursive evaluator that the factoring reference does
not use. Combining that evaluator with the reference's top multiplication
network gives a mathematical route to

```math
 T(n)=2^{O(n/\sqrt{\log n})},
```

where `n` is the input bit length. This is the rate of the reference's
`SqrtLogLogConjecture`, up to harmless integer rounding and constant factors.

The source theorems cited below are proved Lean declarations in the pinned
upstream project. This repository has not integrated and rebuilt them in
its factoring development. The new connection between the two projects is
**not yet a checked factoring-machine theorem**. In particular, a
value-preserving table adapter and the machine composition described in
Section 6 still need proofs. The `alpha = 1/2` factoring rate in this note
remains conditional on those obligations.

The separate [certified k6 result](certified-k6.md) has exponent `alpha = 1/30`
and uses `G = 5K^6`, `K_next = 25K^12`. Its
[closed machine theorem](../variants/k6/FactoringK6/Certified.lean) does not
use Appendix B. Its verification does not certify the stronger route below.

Sources inspected:

- Factoring reference: `IsaacRudich/SubexponentialDeterministicFactoring`,
  commit `c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.
- Appendix A/B formalization: `anthropics/formal-math`,
  commit `16a5ce41304fcdaea4e0be7c132e1c361dc26da5`.

The pinned projects use different Lean/mathlib versions. The factoring
reference uses Lean 4.35.0; the inspected ETH project uses Lean 4.33.1.
Import compatibility must be checked before treating their theorems as one
compiled development.

## 1. The existing evaluator is a machine theorem for written tables

The useful interface is `AppRunningTime.PartB`, rather than the weaker
interface stated only in counted table operations. It asserts the existence
of a function

```lean
run : U × U × PInst → ℤ
```

and a natural-valued bound `need l lam B`, with these properties:

1. `PolyT run`.
2. If `lam ≥ 1`, `B ≥ 1`, `WFe (K x₀ l) lam B I`, and the pad has at least
   `2 ^ need l lam B` units, then
   `run (pad, List.replicate l (), I) = zval (K x₀ l) I`.
3. The logarithm of the required pad satisfies

```math
 \operatorname{need}(\ell,\lambda,B)
 \le e_\ell K_\ell^2\lambda
   +C_0K_\ell^2\log_2\bigl(K_\ell(\lambda+B)\bigr).
 \tag{1}
```

Here `U = List Unit`,
`PInst = U × List (List ℤ)`, and a well-formed instance has `K` rows, each
containing exactly `2 ^ ((K - 1) * lam)` integers. Its first component is the
label length in unary. Each entry has absolute value at most `2 ^ B`.

This interface charges for written tables. It does not give the algorithm
free query access to a table supplied by an arbitrary function.

`PolyT` is defined using actual finite-symbol, four-stack machine fragments:
for some fragment and constant `c`, the input word is replaced by the output
word in at most `c * (len input + 1) ^ c` steps. The fragment preserves the
rest of the stacks, so such computations can be composed. For the concrete
input and output types above, the encoding is the ordinary list/product
encoding, with an integer represented by its sign and binary absolute value.
The general freedom to choose an encoding for an abstract state type does
not erase any input tables or integer output in this interface.

The assembled theorem `AppRunningTime.thmRunningTime` proves `PartB` after
the stated pieces hypotheses. The assembled
`AppRunningTime.thmHalfPieces` supplies those hypotheses for

```math
 \rho=1/\sqrt2.
```

Choose the absolute constants `c`, `C₀`, and `x₀` once, with `x₀ ≥ 18` and
`2 ^ x₀ ≥ c + 5`. None depends on the factoring input.

Primary source paths:

- [The exact PartA/PartB/half-pieces statements](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppRunningTimeStatement.lean)
- [Assembly with all hypotheses discharged](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Main.lean)
- [PInst, table indexing, zval, and WF](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/PaddedDescent/F/Data.lean)
- [WFe and the schedule solver interface](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Alg/Interfaces.lean)
- [The machine definition of PolyT](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/machines/ToTuring/PolyT.lean)
- [The concrete data encodings](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/machines/ToTuring/Enc.lean)

## 2. The schedule and a computable pad

The schedule is

```math
 K_0=2^{x_0},\qquad
 g_\ell=2^{3\ell+x_0}K_\ell,\qquad
 K_{\ell+1}=g_\ell^2.
```

The proved PartA conclusions include

```math
\begin{aligned}
 \log_2K_\ell
   &=(3x_0+6)2^\ell-6\ell-2x_0-6,\\
 K_\ell&<K_{\ell+1}\le K_\ell^{10},\\
 \tfrac12\rho^\ell&\le e_\ell\le0.85\rho^\ell.
\end{aligned}
\tag{2}
```

In particular, `log K_l = Θ(2^l)` and
`e_l = O(1 / sqrt(log K_l))`.

There is a small computability issue in the existential PartB statement:
it does not assert that its chosen `need` function is computable. Do not
make the factoring algorithm call that function. Instead set

```math
 r_\ell=2^{\lfloor\ell/2\rfloor},\quad P=K_\ell,
```

and use the explicit integer budget

```math
 q=\left\lfloor\frac{P^2\lambda}{r_\ell}\right\rfloor+1
    +C_0P^2\,\operatorname{size}\bigl(P(\lambda+B)\bigr).
 \tag{3}
```

`Nat.size` is binary bit length. Since
`e_l ≤ rho^l ≤ 1/r_l` and `log₂ t ≤ Nat.size t` for positive `t`, (1)
implies `need ≤ q`. The algorithm writes a pad of exactly `2 ^ q` units.
Every operation in (3) is integer arithmetic. No real power, real logarithm,
noncomputable choice of an input-dependent number, or advice string is
needed at run time.

For sufficiently large `n`, select the largest `ell` such that
`K_ell ^ 40 ≤ n`. Then

```math
 P^{40}\le n<K_{\ell+1}^{40}\le P^{400}.
 \tag{4}
```

The closed form in (2) gives

```math
 \log_2 n<80(3x_0+6)2^\ell.
```

As `2^ell ≤ 2 r_l²`, this implies

```math
 \sqrt{\log_2 n}+1\le C r_\ell
 \tag{5}
```

for a fixed constant `C`. Also `r_l ≤ P`. These are elementary natural
power and square-root inequalities after replacing the real logarithm by
`Nat.log 2 n`. They provide a way to prove the final integer-valued target
without a general real-asymptotic library argument.

The schedule is defined in
[LevelsStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/LevelsStatement.lean).
The existing `ClauseBound.Gen.StPick` interface is proved by
[`ClauseBound.Gen.Pick.stPick`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Gen/Pick.lean),
which supplies a polynomial-time selection routine from a unary upper
bound. An explicit loop can instead test the cutoff (4).

## 3. Reuse the reference's top multiplication network

Let `N ≥ 2` have bit length `n`. Put `s = sOf n` and `S = 2 ^ s`.
The reference supplies

```math
 n\le4S^4,\qquad S^4\le4n+1.
 \tag{6}
```

Write `A = 2S` and `fw = 8S² + 1`, as in the source. Retain its original
top grouping size

```math
 G=3P^8.
```

This top grouping is separate from the faster recursive schedule in
Section 2. There is no need to change the source's universal grouping
lemma or its multiplication construction.

The source constructs `clusterNet s ℤ N k x G hG`. Its value is

```math
 \#\{d:d\mid N,\;\lfloor d/2^k\rfloor=x\}.
 \tag{7}
```

The valid range is `1 ≤ N < 2 ^ (pL s * pL s)` and `0 ≤ k ≤ n`.
The choice of `s` ensures the upper bound on `N`.

`good_exists` proves that a valid grouping of the `G` clusters into `P`
classes exists. The source's `firstGood` search and `firstGood_spec` then
select a valid grouping deterministically. They guarantee at most
`floor(2G/P)` clusters per class and a common pair-label capacity

```math
 L=\left\lfloor
   \frac{P(10A^2fw)+G\,\operatorname{clP}(s,G)}{P^3}
   \right\rfloor.
```

Take `lambda = max 1 L`. Build the source's `Net.toClique lambda` from
that same grouping. Its padding enforces zero in unused channel bits, so
the extra capacity does not multiply the value. The value remains (7).

The existing source volume estimate is

```math
 P^2L\le720S^4+306GS^3.
 \tag{8}
```

The rounding to positive `lambda` adds at most `P²`.

There is a useful exact consequence of (4) and (6). Since `P ≥ 2`,

```math
 P^{40}\le4S^4\quad\Longrightarrow\quad P^9\le S.
```

Indeed, `P^4 ≥ 4` first gives `P^36 ≤ S^4`, and fourth powers are
strictly increasing on nonnegative integers. It follows that

```math
 PGS^3=3P^9S^3\le3S^4=O(n).
 \tag{9}
```

Equations (8)--(9) give `P² lambda = O(n)`. This argument uses integer
powers rather than the informal estimate `G n^(3/4) = O(n/P)`.

Primary source paths:

- [Multiplication gate network and local weights](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean)
- [Cluster value, width, and total bounds](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Cluster.lean)
- [First valid grouping and preservation of value](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineGroup.lean)
- [The existing top-clique construction](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineTop.lean)
- [Top volume estimates](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineCostRate.lean)

## 4. Table entries and materialization cost

All original gate weights and half-edge factors are 0 or 1.
`Fac.eval_numTest` shows that contraction of the two halves of one original
edge is again 0 or 1. The source's `Net.coarsen_F_eval` identifies a
coarsened edge value with the product of its original edge values.
Therefore the first and second coarsenings retain 0-or-1 local weights and
0-or-1 half-edge factors. A final clique-table entry counts compatible
settings in one group.

The checked `cluster_s_le` theorem bounds a cluster's number of settings by
`2 ^ Wc`, where

```math
 W_c=8\bigl(\lfloor A/G\rfloor+2\bigr)Afw
        +2\bigl(\lfloor A^2/G\rfloor+2\bigr)fw.
```

Thus one can take

```math
 B=1+W_c\lfloor2G/P\rfloor
 \tag{10}
```

as a bound on the bit length of every table entry. The extra 1 makes
`B ≥ 1` and permits a strict modulus bound if needed.

Using `fw ≤ 9S²`, the nonconstant part of (10) satisfies

```math
 P(B-1)\le720S^4+648GS^3=O(n).
 \tag{11}
```

One way to prove (11) with no division of real quantities is to multiply
the formulas by `G`, use `G floor(A/G) ≤ A`, then multiply the group-size
inequality `P floor(2G/P) ≤ 2G`.

Materializing all `P` rows has `P 2^((P-1)lambda)` entries. Each entry can
be computed by enumeration of the local group settings and the internal
channel sums already present in the reference. The latter sums add at
most the existing top cost exponent. Searching all groupings costs
`P ^ G` times a polynomial. The source's conservative top requirements are
controlled by the same terms `n/P`, `G S³`, and fixed powers through `P^38`,
plus logarithmic or polynomial factors. By (4) and (9), all these
exponents are `O(n/P)`.

There are two implementation choices:

1. Compute the nonnegative integer table entries directly with binary
   arithmetic, using the integer bound (10).
2. Reuse the reference's modular `GroupEx.tab`, with modulus
   `m = 2 ^ B`. Every entry is strictly smaller than `m`, so its residue
   is its exact nonnegative integer value. The reference's unary-style
   value-bound transfer then still costs only `2 ^ O(n/P)`, because
   `log₂ m = B = O(n/P)`.

The local bound in (10) is essential for the second choice. Taking a crude
global bound `B = O(n)` would make the reference's bound polynomial in
the *numeric* modulus as large as `2 ^ O(n)`. A claim that such a choice
preserves the desired rate would be incorrect.

The exact top-cost extraction and integer-entry lemma remain new formal
obligations. Citing the full existing `topE` cost is insufficient because
that expression also evaluates the old recursive chain; it must be
stopped after its top-table constructor.

Relevant source paths:

- [Equality factor contraction](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Gates.lean)
- [Coarsening formulas and coarsen_F_eval](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Net.lean)
- [Wc, cluster_s_le, and TopReq](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineCostTop.lean)
- [Grouped settings bound and table-expression bound](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineCostGroup.lean)
- [Why the old transfer pays polynomially in numeric V](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/UniGlue.lean)

## 5. Conditional complexity calculation

The first term of (3) is `O(n/r_l)` by (8). The second is

```math
 O(P^2\log n)=O(n/r_\ell),
```

using `n < P^400`, `P^40 ≤ n`, and `r_l ≤ P`. The table input itself has
length `2 ^ O(n/P)`, which is at most `2 ^ O(n/r_l)`. The explicit pad
has length `2 ^ O(n/r_l)`. Consequently the `PolyT run` machine has
running time `2 ^ O(n/r_l)` as well: its fixed polynomial overhead merely
changes the absolute constant in the exponent.

After the constructor and machine-composition obligations in Section 6 are
proved, the top-table constructor has the same or smaller bound. One exact
divisor-prefix count would then take `2 ^ O(n/r_l)` steps.

Subtract the contribution from divisor 1. The reference's prefix search
then uses at most `n` such counts to recover `Nat.minFac N`. This only
adds a polynomial factor. Inputs 0 and 1 and all bit lengths below the
fixed first cutoff use a finite fallback. Formula (5) gives the desired
bound

```math
 T(n)\le C\,2^{C n/(\operatorname{Nat.sqrt}(\operatorname{Nat.log}_2n)+1)}
```

after increasing the absolute constant. Integer divisions and the finite
small-input range need to be handled in the final formal statement.

## 6. Minimal formal bridge and verification gates

The following are genuinely new obligations; none is discharged merely by
the existence of the reference's old factoring machine.

1. **Shared compatible dependencies.** Import the assembled Appendix B
   theorem and the pinned factoring network into one Lean project, and
   check the theorem's axioms. Do not import a statement-only file in place
   of its assembly proof.
2. **Clique codec.** The source clique accepts a function
   `Fin P → Bits lambda` at a vertex and forces the self-coordinate to
   zero. The Appendix B table codec concatenates the `P-1` off-diagonal
   strings, using `[]` for the self-coordinate. Prove the insertion/removal
   of the diagonal and the block codec preserve the summed clique value.
   The result must be `WFe P lambda B I` and equality of `zval P I` with
   the count in (7).
3. **Bounded top constructor.** Extract the top search and table expression
   before `deepEx`. Prove its uniform computation bound, the entry bound
   (10), and its compatibility with the table codec. If reusing the modular
   interpreter, prove residues equal the integer entries with the local
   modulus from Section 4.
4. **Actual machine composition.** Build the pad in (3), call the finite
   fragment supplied by `PolyT`, read its encoded integer, and compose it
   with the source prefix search. A theorem about counted operations alone
   is insufficient. The final result must be a single finite-alphabet
   machine computing `Nat.minFac` under the reference's encoding.
5. **Uniform rate.** Prove the natural inequalities (4)--(5), (8)--(11),
   absorb the fixed polynomial machine overhead, and discharge the finite
   fallback. Then state the reference's square-root-logarithm target and
   check the final theorem and its axioms in CI.

The stronger decomposition, bounded-edge grouping, recursion, and
finite-symbol table evaluator already exist in the pinned upstream project.
Connecting them to factoring still requires the formal integration above.
Until those obligations are checked, this note supports a conditional
`alpha = 1/2` route. The repository's separate certified machine exponent
remains `1/30`.
