# The Appendix B evaluator: machine interface and square-root saving

## Scope

This note audits the existing Appendix B clique evaluator in
[`anthropics/formal-math`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth),
commit `16a5ce41304fcdaea4e0be7c132e1c361dc26da5`. It proves an ordinary
mathematical composition lemma for using that evaluator after a top-table
constructor. It also proves a limitation of the declared Appendix B channel
widths. It does not claim to improve the Appendix B evaluator itself.

The source theorem supplies one finite-alphabet machine fragment, not an
oracle that evaluates arbitrary functions. Its input contains all clique
tables in a fixed binary encoding. The fragment and all constants can be
chosen once, independently of the factoring input and recursion depth.

The composition and width-limitation arguments below are handwritten. The
pinned source was subsequently rebuilt, and
[ExplicitPad.lean](../formal-appendix/ExplicitPad.lean) checks the assembled
evaluator with an explicit pad under Lean 4.33.1. Compatibility with the
factoring reference's Lean 4.35.0 development and the complete machine
composition remain formalization tasks.

## 1. The theorem that is actually available

The statements are in
[`AppRunningTimeStatement.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppRunningTimeStatement.lean).
The proofs are assembled in
[`GeneralPieces/Main.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Main.lean).
Using the statement file alone would not discharge its hypotheses.

`AppRunningTime.thmHalfPieces` supplies an absolute integer `c` and the
programs satisfying `IsPieces (1 / sqrt 2) c ...`. In particular, its
factorizations, coefficient bounds, and running-time certificates are
included. The assembled proof supplies the ring identities called
`AppHalfPieces.pieceFacts`; they are not left as assumptions.

`AppRunningTime.thmRunningTime.b` states

```text
for every c, there exists C0 such that for every rho and x0,
  Exist rho c, x0 >= 18, and c + 5 <= 2^x0
  imply PartB rho x0 C0.
```

Choose the half-pieces family once, choose its `C0` once, and fix any
`x0 >= max(18, ceil(log2(c+5)))`. Then `PartB` supplies functions

```text
run  : U × U × PInst -> integer
need : natural × natural × natural -> natural
```

with `PolyT run`. For every level `l`, positive integers `lambda,B`,
well-formed `K_l`-node instance `I`, and pad of at least
`2^need(l,lambda,B)` units, `run(pad, unary(l), I)` equals the clique value.
The bound is

```math
\mathrm{need}(\ell,\lambda,B)
\le e_\ell K_\ell^2\lambda
   +C_0K_\ell^2\log_2\bigl(K_\ell(\lambda+B)\bigr).                 \tag{1}
```

The schedule and its estimates are

```math
\begin{aligned}
K_0&=2^{x_0}, & g_\ell&=2^{3\ell+x_0}K_\ell,
&K_{\ell+1}&=g_\ell^2,\\
\log_2 K_\ell&=(3x_0+6)2^\ell-6\ell-2x_0-6,\\
K_\ell&<K_{\ell+1}\le K_\ell^{10},
&\tfrac12\rho^\ell&\le e_\ell\le0.85\rho^\ell,
&\rho&=1/\sqrt2.
\end{aligned}                                                        \tag{2}
```

These are `thmRunningTime.a`, with the definitions in
[`LevelsStatement.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/LevelsStatement.lean).

### The existential statement comes from a named uniform recursion

[`GeneralPieces/Rec/PartB.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Rec/PartB.lean)
proves the stronger named interface `GeneralPieces.Tied.PartB` and then
forgets the names. The source chain is:

| Required fact | Assembled proof |
| --- | --- |
| Correct descent through all levels | `solveG2`, using `Descend.stSolveG2` and `Level.stExpandG2` |
| A sufficient finite pad | `Need.stNeedG` |
| Rate of that pad | `stRate2` |
| Polynomial-time finite-symbol implementation | `certSolveG`, then `runB_polyT` |
| Correct schedule chosen from table data | `pick_rows` |
| All clauses together | `partB_tied`, then `partB_of_tied` |

`Tied.runB` reads the number of rows from the written table and chooses
the matching schedule. Its depth argument is not used to conjure a long
schedule from a short input. Every recursive level meets case II of the
source grouping lemma. The bound includes signed arithmetic, intermediate
instances, factors, and the complete descent, through the source's pad
and `PolyT` certificates.

The conversions from the half-pieces theorem to the recursion are in
[`GeneralPieces/ThmRunningTime.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/ThmRunningTime.lean).
`par_of_isPieces` and `exist_par_one` perform the permitted trimming of
factorizations to at most `2 lambda` bits. The final `thmRunningTime_b`
passes their conclusions to `partB_tied`; it does not assume a solver for
the new factoring problem.

## 2. Exact input and output encodings

[`PaddedDescent/F/Data.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/PaddedDescent/F/Data.lean)
defines

```text
PInst = U × List (List integer).
```

The first component is `lambda` in unary. A `K`-node well-formed instance
has `K` rows of length `J=2^((K-1)lambda)`, with every entry `z` satisfying
`abs(z) <= 2^B`. The exact-width predicate `WFe` additionally equates the
stored unary length to `lambda`.

The ordinary instances of `Enc` in
[`ToTuring/Enc.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/machines/ToTuring/Enc.lean)
encode a natural number as its binary digits and an integer as its sign
paired with its binary absolute value. Products and lists use bracketed
binary trees. In
[`TM/Tree.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/machines/ToTuring/TM/Tree.lean)
the alphabet has exactly four symbols: two bits and two brackets.

Write `s(z)` for the bit length of `abs(z)`, with `s(0)=0`. Directly from
these encodings:

```math
\begin{aligned}
|\mathrm{enc}(U_m)|&=4m+2,\\
|\mathrm{enc}(z)|&=3s(z)+5,\\
|\mathrm{enc}(\mathrm{pad}_F, U_\ell,I)|
 &=14+4(F+\ell+\lambda+K)+7KJ+3\sum_z s(z).
\end{aligned}
\tag{3}
```

The sum ranges over every table entry, including zeros. Since the bound
on entries is inclusive, `s(z) <= B+1`, not always `B`. Thus

```math
|\mathrm{enc}(\mathrm{pad}_F,U_\ell,I)|
\le14+4(F+\ell+\lambda+K)+(3B+10)K2^{(K-1)\lambda}.               \tag{4}
```

`PolyT run`, defined in
[`ToTuring/PolyT.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/machines/ToTuring/PolyT.lean),
means that a single fragment takes at most `a*(L+1)^a` steps on an
encoded input of length `L`, for an absolute integer `a`. It returns the
actual encoded output integer and preserves the surrounding stack contents.
The definition also bounds the length of the written output.

The fragment has finite control, four work stacks, and finite local state;
see `Frag`, `Fr`, and `Sym` in the pinned machine files. Each instruction
is a fixed finite statement. Expanding these fixed statements to ordinary
finite-alphabet Turing-machine transitions changes time by a fixed factor.
A usual simulation between a fixed number of stacks and tapes has fixed
polynomial overhead. Both effects preserve the logarithmic saving exponent.

Some *internal state types* elsewhere in the source permit custom encodings.
That flexibility does not shorten the explicit `PInst` input here. Its
integer, list, and product encodings are fixed as above and determine the
tables and returned integer. There is no free table-query oracle.

## 3. An explicit pad; no input-dependent choice oracle

The public existential `PartB` does not promise computability of its
chosen function `need`. The algorithm must not call that existential
function without another argument. No such call is necessary.

Put `P=K_l`, `r=2^floor(l/2)`, and let `size(t)` be binary bit length.
Use the integer

```math
q=\left\lfloor\frac{P^2\lambda}{r}\right\rfloor+1
  +C_0P^2\mathrm{size}\bigl(P(\lambda+B)\bigr).             \tag{5}
```

From (2), `e_l <= rho^l <= 1/r`. Also `size(t) >= log2(t)` for
positive integers. The floor plus one is at least the unfloored quotient.
Consequently (1) gives `need(l,lambda,B) <= q`. A pad of exactly `2^q`
units therefore suffices.

All operations in (5) are integer operations. The same fixed constants
`c,C0,x0` and the same fragment are used on every input. Selecting absolute
constants once from their existence proofs is permitted in a theorem
asserting existence of one uniform machine. It is not a supply of advice
depending on input length. No explicit small numerical values or practical
implementation of those constants is claimed here.

## 4. Composition lemma for a top constructor

Suppose that for an `n`-bit factoring input and each requested divisor
prefix, a single deterministic constructor produces a written instance
`I` with these properties:

1. Choose the largest level `l` with `P^40 <= n`, where `P=K_l`.
2. `I` satisfies `WFe P lambda B`, with `lambda,B >= 1`.
3. Its value is exactly the requested divisor-prefix count.
4. `P^2 lambda <= A n` and `B <= A n/P` for a fixed constant `A >= 1`.
5. Construction takes `2^O(n/P)` bit operations, uniformly in the input
   and prefix.

Then one divisor-prefix count can be computed in

```math
2^{O(n/\sqrt{\log n})}                                          \tag{6}
```

bit operations by one fixed deterministic finite-alphabet machine.

### Proof: input length and pad

Maximality of `l` and (2) give

```math
P^{40}\le n<K_{\ell+1}^{40}\le P^{400},\qquad r\le P.            \tag{7}
```

The bound `r<=P` follows directly from the schedule: `K_l` grows at least
as fast as `2^l`. Moreover,

```math
P(\lambda+B)\le A n/P+A n\le2An.
```

Therefore the second term of (5) is `O(P^2 log n)`. It is
`O(n/r)`: from (7),

```math
rP^2\log_2 n\le400P^3\log_2P\le400P^4\le400n.                 \tag{8}
```

Thus `q=O(n/r)`. In (4), the table exponent satisfies
`(P-1)lambda <= A n/P`; the table-entry bit bound and all prefactors
are polynomial in `n`. Since `r<=P`, the whole encoded input, including
the written pad, has length `2^O(n/r)`.

Writing the pad by repeated doubling or a binary counter takes time
polynomial in its written length. Constructing the bracket encoding of
the tables and reading the signed binary output also take polynomial
time in their written lengths. The fixed polynomial bound of `PolyT`
therefore yields `2^O(n/r)` steps. The top constructor takes no more.

### Proof: level versus input length

The closed schedule formula in (2), used at `l+1`, gives

```math
\log_2n<80(3x_0+6)2^\ell\le160(3x_0+6)r^2.                    \tag{9}
```

It follows that `1/r=O(1/sqrt(log n))`, proving (6). The floor in
`r=2^floor(l/2)` changes only the factor two in (9).

The level can be found using successive integer powers and comparisons.
The next candidate `K_(l+1)` has `O(log n)` bits by (7), so this does not
require superpolynomial preprocessing. All lengths below the fixed first
cutoff can use a fixed deterministic fallback. There are finitely many
such inputs, so a larger absolute constant covers them.

Finally, at most `n` prefix queries find the least divisor greater than
one. Subtract divisor one's contribution from each count and descend to
the lexicographically least surviving fixed-width binary prefix. That
divisor is prime, since a composite least divisor would have a smaller
nontrivial divisor. The factor `n` and binary bookkeeping are absorbed
by (6). This last argument needs the top constructor's exact value
identity; the evaluator alone does not supply it.

## 5. Why the declared half-pieces widths stop at exponent one half

This section concerns only the channel counts used in the published
Appendix B construction. It is not a lower bound on tensor rank, all
clique algorithms, or factoring.

Let `p(h)` be the partition number, and let `h=h(lambda)` be the least
integer at least two with `p(h)>=2^lambda`. The central piece `s=0`
has the declared counts

```math
R_{\rm row}(h,0)=R_{\rm col}(h,0)=\sum_{j=0}^{2h}p(j).           \tag{10}
```

They appear explicitly in
[`AppHalfPiecesStatement.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
The source proves, for the same positive constant
`kappa=pi*sqrt(2/3)/log(2)`,

```math
\log_2p(E)=\kappa\sqrt E+O(\log(E+1)).                          \tag{11}
```

This follows from the two bounds in
[`ClauseBound/Part/Bounds.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Part/Bounds.lean).
Minimality of `h` gives
`log2 p(h-1) < lambda <= log2 p(h)` for positive `lambda`.
Together with (11) and `h=O(lambda^2)`, proved in
[`GeneralPieces/Ring/Bits.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Ring/Bits.lean),
this implies

```math
\kappa\sqrt h=\lambda+O(\log(\lambda+1)).                       \tag{12}
```

Monotonicity of the partition numbers gives

```math
p(2h)\le\sum_{j=0}^{2h}p(j)\le(2h+1)p(2h).
```

Using (11)--(12) and then rounding up to bit widths,

```math
\left\lceil\log_2R_{\rm row}(h,0)\right\rceil
+\left\lceil\log_2R_{\rm col}(h,0)\right\rceil
=2\sqrt2\,\lambda+O(\log(\lambda+1)).                           \tag{13}
```

If these same widths satisfied the pieces condition with a fixed
`rho<1/sqrt(2)`, it would require the left side of (13) to be at most
`4 rho lambda + c log2(2 lambda)`. The positive linear gap
`(2 sqrt(2)-4 rho)lambda` cannot be absorbed by `O(log lambda)`.
Consequently no such `rho` works for these declared widths. Trimming a
width to `2 lambda` does not help this central piece asymptotically,
because each of its widths is `sqrt(2)lambda+O(log lambda)`.

The general recursion gives exponent `log2(1/rho)`, so its present
half-pieces interface reaches exactly `1/2`. Better constants, a different
fixed initial level, or tighter rounding in the existing counts cannot
yield a larger fixed exponent. A larger exponent would require smaller
effective factorizations, a different pieces construction, or a different
recursion argument. Equation (13) does not rule out those changes.

## 6. Independent review of the written top constructor

The companion [top-constructor proof](half-top-constructor.md) supplies
the five properties used in Section 4, with `A=6769`. This review checked
its construction against the clean factoring-reference commit
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
It found no additional mathematical premise in the interface.

In particular, `MulNet.mS` is a family of explicit finite tuples. Input
gates store bounded binary digits. The middle and inverse gates store
tuples of residues modulo `2^(8S^2)+1`. Carry gates store one residue and
one bounded carry. `EngineMul.gcode` and `gsz` express the same sets by
mixed-radix codes. Their coordinates can be enumerated with ordinary
binary counters. The cardinality bound therefore corresponds to an
effective enumeration, rather than merely a finite set whose elements
lack an algorithmic representation.

The source `mW` is explicitly zero or one. Its `mV` messages are finite
sums, products, divisions, and powers of bounded integers or residues.
The largest message has `8S^2+1=O(sqrt(n))` bits. The functions `faOut`,
`midOut`, and `ibOut` have polynomially many summands and polynomially
bounded exponent indices. Repeated squaring computes each modular power
in polynomial bit time. The inverse-transform multiplier `tInv` is
explicitly `2^((2*pm-1)*(2s+2))` modulo the same ring modulus; no search for
an inverse is required. The prefix and carry tests also have polynomial
index ranges. A fixed setting can thus be checked in `n^c` bit operations
for one absolute `c`.

The central cost point survives expansion of the actual definitions.
`Fac.numTest` accepts a channel only when it equals the endpoint's
deterministic message. The messages fit their channel widths. Each full
internal factor is consequently exactly a zero-or-one equality test.
`Net.coarsen` groups original settings and incorporates those equality
tests into the node weight. It does not require enumerating independent
internal channel assignments to compute that weight. Its `orient`
definition also selects the correct original endpoint when group order
reverses an edge. Using one direct grouping of the original gates avoids
any dependence on the block order of an intermediate grouped table.

The grouping theorem's hypotheses match `G=3P^8`, `P>=2`, the proved
cluster degree bound, and the ordered total-width bound. The claimed
entry exponent `H`, width bound `P^2 lambda<=4717n`, and entry-bit bound
`B<=6769n/P` follow with the stated floor rounding. The exhaustive search
cost `P^G` and the full entry-by-entry enumeration are both included in
`2^O(n/P)`. `MulCount.mulNet_value` supplies the exact requested count
over the integers, not merely its residue modulo a bound-dependent
modulus.

Finally, the source's zero self block and the target's empty self block
give the same off-diagonal flattening order. Zero padding gives exactly
one extension of a legal crossing label. I independently replayed
`python3 -m unittest tests.test_half_clique_codec -v`; all seven tests
passed. They include an intentional one-endpoint block-order error, which
changes the value as expected. These finite tests check the codec and
grouping identities; the source inspection and estimates above justify
the uniform asymptotic construction.

## 7. Audit scope and source hashes

This audit traced the assembled theorem, the named recursion, the machine
and data definitions, and the top-constructor properties above. The new
composition is an ordinary written mathematical proof. The later Lean check
certifies the evaluator and explicit pad, not the composition or factoring
endpoint.

The following are SHA-256 hashes of the inspected source files at the
pinned commit. Paths are relative to `eth/proof/`.

| Source | SHA-256 |
| --- | --- |
| `statements/AppRunningTimeStatement.lean` | `07725ec86e202a1b23b7621f46cf6830d42bb9df6e000e1ae9ead7848c6c64ba` |
| `rates/GeneralPieces/Main.lean` | `6e3bc1a617e2a74445d280ddc4374d13d2c4fe5d9f396b812111ed08c586376c` |
| `rates/GeneralPieces/Rec/PartB.lean` | `69aa82cd09d7ebaca0ee485c4f79b9c1f50f7a7ca196ac83d60fd254d2d4318d` |
| `machines/ToTuring/PolyT.lean` | `28925257b2c11e91b40c499deea5630b6613b4d8466d3954a09be15923cc67a7` |
| `machines/ToTuring/Enc.lean` | `66aa4b980b5c1d3191489a04f6fe5e2243c9a31e6c04c569e95b23ba2d45ba38` |
| `rates/PaddedDescent/F/Data.lean` | `08f9ef710e1229529285c9df75e5189efe935555b4f6c12ac1b1cea2497a9f9d` |
| `rates/ClauseBound/Part/Bounds.lean` | `292273a74d9897da78e530b95920a8f486b33d7c7f2422bf3d6c064b09c7278d` |

The following additional SHA-256 hashes identify key files inspected in
the pinned factoring reference. Paths are relative to its repository root.

| Source | SHA-256 |
| --- | --- |
| `FactoringAux/Gates.lean` | `c93e99043116c086748626882de1694635a3186f8545c713b0a5be80739a450c` |
| `FactoringAux/MulNet.lean` | `5f01187574f1bc08a442a6c9f8fde654a15d25224f95f410f1c3cd85c9df0b67` |
| `FactoringAux/EngineMul.lean` | `23e9b8a165596171f5e55fc38b922c7a6ae8125c707573f038621393d525a3f7` |
| `FactoringAux/Net.lean` | `4fc0003f2560e30f5998ac52f65058210a3df6d574382bf68078d7e7ab4b76e6` |
| `FactoringAux/Transform.lean` | `c62e78903c6ab7c9ca26e201dbcd5c5928a4d2128d0eb47af22c2be0629e7e68` |
| `FactoringAux/EngineGroup.lean` | `f190262f74b1e92dd2abab68c9d7c3fc48972e78c72cfac6baf2ceea13f83ebe` |
| `FactoringAux/MulCount.lean` | `e334f81b9c7558c36e3151287f1d9c3a883a86e35b411e6898a7ac4e18d1074f` |

The evaluator, grouping algorithm, half-pieces construction, partition
bounds, and finite-symbol implementation are prior work from the pinned
Anthropic development. This note's contribution is the explicit pad and
encoding accounting, the stated top-constructor composition lemma, and
the limitation of these particular declared channel counts.
