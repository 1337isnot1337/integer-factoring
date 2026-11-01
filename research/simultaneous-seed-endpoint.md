# Simultaneous time and workspace at the finite-seed endpoint

## Result and dependencies

The best limiting score supplied by an enumerable stream of finite crossing
certificates is attainable by one machine with the same saving in workspace.
The seed generator may be arbitrarily slow, the good widths may have arbitrary
gaps, and the coefficient fields may vary. The conclusion concerns the
challenge's limiting exponent. It does not assert a bounded multiplicative
constant in the exponent at the endpoint.

For a seed on `d=2^a` labels over an explicitly presented finite field, write

```math
E_d(x,y,u,v)=[x=v][y=u]=\sum_{i=1}^{M} A_i(x,y,u,v).
```

Discard zero tensors. Let `R_i,C_i` be their exact ranks in the cuts
`(x,y)|(u,v)` and `(x,u)|(y,v)`. Put

```math
w=\max_i\bigl(\lceil\log_2R_i\rceil+\lceil\log_2C_i\rceil\bigr),
\qquad c=\lceil\log_2M\rceil,\qquad Q(S)=\frac{w+2c}{4a}.
\tag{1}
```

Fix a computable stream of complete finite seed records, allowing invalid
records. Include the one-piece crossing over `F_2` as a permanent fallback,
and define

```math
\rho_* =\inf_{S\text{ valid in the stream or the fallback}}Q(S).
\tag{2}
```

**Theorem.** If `rho_*<1`, one deterministic finite-alphabet machine computes
the least prime factor of every binary integer. With
`alpha_*=log_2(1/rho_*)`, for every fixed `epsilon>0`, its worst-case
running time and working space on `n`-bit inputs satisfy

```math
T(n)\le 2^{\,n/(\log_2 n)^{\alpha_*-\epsilon}},\qquad
S(n)\le \frac{n}{(\log_2 n)^{\alpha_*-\epsilon}}
\tag{3}
```

for all sufficiently large `n`. The same machine works for every epsilon.
Input is read-only and output is streamed. The conclusion also holds for
complete prime factorization, sorted divisors, and exact divisor-prefix or
interval counts with read-only endpoints.

Constants and eventual thresholds may depend on the fixed enumeration program
and its infimum. The program, the evaluator schedule, and the local routines
are fixed independently of the input length and of epsilon.

This is an ordinary transfer theorem, not a new numerical factoring bound.
The known unconditional exponent remains `1/2`. It adds workspace to the
existing endpoint selection; it does not supply a crossing seed below the
half-rate threshold.

The proof uses these immutable research dependencies:

- [Packed time and space transfer, Sections 1–2 and 4](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/simultaneous-time-space-transfer.md): streamed top construction, grid recurrence, type-word channels, and packed assignments.
- [Gap-free finite-seed transfer](https://github.com/1337isnot1337/integer-factoring/blob/023e4fc4482d72f25738f6b1d0f9cd31583baed1/research/gap-free-seed-score-transfer.md) and [field selection](https://github.com/1337isnot1337/integer-factoring/blob/023e4fc4482d72f25738f6b1d0f9cd31583baed1/research/universal-field-seed-selection.md): the time-only endpoint and its rational quality selector.
- [Recovery in varying characteristic](https://github.com/1337isnot1337/integer-factoring/blob/69ee5a3aa8f575846a94be5c92a51bd8af8cb3a5/research/characteristic-independent-recovery.md): short verified divisor descriptors over one explicitly represented finite field per input.

The underlying evaluator is Anthropic's Appendix A/B at
[`16a5ce4`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth),
and the multiplication-network reduction is Rudich and Rousseau's at
[`c91cc3f`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
No new Lean proof is claimed.

## 1. A bounded global selector fits sublinear space

Choose a fixed self-delimiting encoding for field and tensor records. A field
record consists of a prime `p`, a positive degree `e`, and a monic irreducible
polynomial over `F_p`; elements have their polynomial-basis encoding. Tensor
entries are explicit. There are no compressed tensor or rank oracles.

For the original input length `n`, put

```math
B(n)=\max\{1,\lfloor\log_2\log_2(n+4)\rfloor\},\qquad
u(n)=\lfloor n^{1/64}\rfloor.
\tag{4}
```

These choices are implemented with bit lengths and integer roots; changing
either displayed floor by an additive constant has no effect below. Simulate
the fixed stream program for `u(n)` steps and retain its output prefix. A
direct simulation uses `O(u(n))` working cells, since the program and its
alphabet are fixed. A scan-based implementation has polynomial time in
`u(n)`, which is harmless here. Parse complete records and consider only those
whose entire encoding and each of `a,M,ceil(log_2(p^e))` are at most `B(n)`.
Reject malformed headers against the actual recorded length before allocating
any asserted array. The fallback is considered regardless of these filters.

Validate each candidate exactly. Primality by trial division and testing all
monic possible polynomial factors suffice: `p^e<=2^B=O(log(n+4))` for an
accepted candidate. Check its crossing sum, eliminate zero pieces, and compute
both ranks and rank factorizations over its field. The arrays have polynomial
size in the explicit record length. Thus all this validation uses
`poly(log(n+4))` space and time per candidate. Retain a minimum-quality valid
record, using exact rational comparisons in (1), with a fixed tie rule. Erase
the simulation and unselected output before evaluating the factoring network.

All preprocessing has polynomial time in `n` and workspace

```math
O(n^{1/64}+\operatorname{polylog}n).
\tag{5}
```

Its polynomial exponent is irrelevant to (3); the displayed workspace bound
is the reason to limit simulation by `n^{1/64}` rather than by `n`.
Store only the selected tables, their rank factorizations, and the finite-field
arithmetic data, of size `poly(B)`.

Write the selected field and quality as `F_n,Q_n`. Use this same field and seed
for every level and every weighted query on the input. Every fixed valid
record eventually appears within `u(n)` steps and passes all the size filters.
Therefore

```math
\rho_*\le Q_n\le1,\qquad Q_n\longrightarrow\rho_*.
\tag{6}
```

The algorithm never computes `rho_*` or a convergence modulus. Invalid
records cannot change correctness or the infimum because they are rejected.
The fallback and exact grid identities give termination on all inputs,
including the finite range before any useful record is selected.

For every valid seed, rank subadditivity gives `sum R_i>=d^2`, hence some
`R_i>=d^2/M`. Consequently `w>=2a-log_2M`, and (1) implies `Q>=1/2`.
Thus `rho_*>=1/2` and `alpha_*<=1`.

## 2. Amplification with explicit channels, without dense expanded matrices

Fix the selected seed for one input. At a recursive label width `lambda`,
write `lambda=ka+b`, where `0<=b<a`. For each count vector
`nu=(nu_1,...,nu_M)` with `sum nu_i=k`, form

```math
B_\nu=\left(\sum_{\operatorname{type}(v)=\nu}
 A_{v_1}\otimes\cdots\otimes A_{v_k}\right)\otimes E_{2^b}.
\tag{7}
```

The empty tail is omitted. Summing over types gives the full crossing exactly,
by distributivity; no factorial is inverted in the coefficient field.
Keep zero pieces as zero pieces, so deciding whether a large tensor vanishes
is unnecessary. The number of count vectors is at most `(k+1)^M`.

Let `D_nu=k!/prod nu_i!`. An explicit row channel consists of a word of type
`nu`, one seed row channel at each position, and one full tail channel. Its
cardinality, before power-of-two padding, is

```math
\mathcal R_\nu=D_\nu\prod_i R_i^{\nu_i}\,2^{2b};
\qquad
\mathcal C_\nu=D_\nu\prod_i C_i^{\nu_i}\,2^{2b}
\tag{8}
```

for the column cut. Products of the corresponding left and right seed factors
give the two factorizations. Padding the two channel sets once, after forming
their complete cardinalities, gives total width at most

```math
2\log_2D_\nu+\sum_i\nu_i(\log_2R_i+\log_2C_i)+4b+2
\le4Q_n\lambda+4a+2.
\tag{9}
```

Here `D_nu<=M^k` and (1) were used. If either padded width exceeds `2lambda`,
replace that cut by the identity matrix times the complete piece matrix. The
complete-entry routine below makes this cap effective, and capping only lowers
the total width.

Because `Q_n<=1`, equation (1) gives `c<=2a`, hence

```math
\log_2M\le2a,\qquad k\lceil\log_2M\rceil\le2\lambda.
\tag{10}
```

This uniform bound is essential. It pays for a word in `[M]^k` even when the
selected seed changes with `n`.

To compute a factor entry, unrank its type word using first-letter counts
`D_nu*nu_i/k`, then decode the remaining mixed-radix local channel index.
The product of local channel cardinalities depends only on the type, so this
is a bijective numbering. Multiply the seed factor entries sequentially.
For a capped factor or complete-piece entry, enumerate all words in `[M]^k`,
keep those of the required type, and add their products of seed entries.
Equation (10) bounds this by `2^{2lambda}` words and `O(lambda)` word bits.

Avoid storing `k!`, which could take `Theta(k log k)` bits. Instead compute
`D_nu` as a product of binomial coefficients. Replace each lower argument
`v` by `min(v,m-v)`, and use the integer recurrence
`z_(j+1)=z_j*(m-j)/(j+1)` only through that smaller argument. The partial
binomial coefficients are then increasing. Each partial result is
at most `D_nu<=M^k<=2^{2lambda}`; a numerator multiplication adds only
`O(log(k+1))` bits. The same bounded arithmetic computes the first-letter
counts. These are integer indexing operations, not divisions in `F_n`.

At each active recursive call the packed-grid cutoff gives

```math
\lambda\ge n^{1/4}/8,\qquad \lambda\le 2^{16}(n+2)^3.
\tag{11}
```

The upper bound is the crude separate-width-cap envelope and does not assume
the sharper contraction. Thus `M log(lambda+1)=O(B log(n+2))=o(lambda)`.
The selected raw tables occupy `poly(B)=o(lambda)` space. Exact field
arithmetic has polylogarithmic bit cost, and only a constant number of field
elements is live in a local query. Consequently all expanded-family queries
and streamed enumeration take `2^{O(lambda)}` bit time and `O(lambda)`
working bits, with constants independent of the chosen seed. No matrix of
order `2^{2lambda}` is materialized or eliminated.

The type count, descriptor size, and piece-ordinal width obey the common bound

```math
f_n(\lambda)=O(M\log(\lambda+1))=O(B\log(n+2)).
\tag{12}
```

## 3. Packed recursion with an input-dependent limiting contraction

Use the fixed source schedule

```math
K_0=2^{x_0},\quad g_j=2^{3j+x_0}K_j,\quad
K_{j+1}=g_j^2,\quad \zeta_j=2^{-j-5},
\tag{13}
```

the streamed top construction, and the `n^{3/4}` early-volume cutoff from
the packed transfer. The fixed constant `x_0` covers the uniform local
algorithms above. That construction supplies top volume `V_ell=O(n)`,
`2^ell=Theta(log n)`, and (11) at every call performing a grid step.
If the input is too short for a scheduled top, use ordinary trial division;
this is one fixed finite range. A valid top of level zero is an ordinary
finite leaf and is evaluated by complete assignment enumeration. Neither
case requires evaluating the displayed asymptotic bounds.

The grouping estimate applied to (9), with child-width rounding included,
has the form

```math
V_{j-1}\le(\beta_n+3\zeta_{j-1})V_j,
\qquad
\beta_n=Q_n+C(B+1)n^{-1/4}
\tag{14}
```

for one fixed sufficiently large `C`. This `beta_n` is an analytical upper
bound, not a parameter queried by the machine. Since `beta_n -> rho_*<1`
and `beta_n>=1/2`, the source product estimate is uniform for large `n`:

```math
V_j\le D\beta_n^{\ell-j}n,\qquad \sum_jV_j=O(n).
\tag{15}
```

Indeed, divide the factors in (14) by `beta_n` and bound their product by
`exp(6 sum_j zeta_j)`. Eventually `beta_n` is also bounded above by one
fixed constant below one, which gives the second inequality.

Along an evaluation path, packed type ordinals and width metadata occupy

```math
O\!\left(\sum_jK_j^2 B\log(n+2)\right)
=O(n^{3/4}B\log(n+2))
\tag{16}
```

bits, using `K_j^2=V_j/lambda_j`, (11), and (15). The same expression bounds
the logarithm of the piece-assignment branching count. Raw descriptors need
not be retained at every crossing site: the ordinal regenerates the type.

Grouped assignments, upward query arguments, and one active local scratch
area are bounded by

```math
\sum_{j=1}^{\ell}O(V_j/K_{j-1})
\le O(n\beta_n^\ell)
       \sum_{j\ge1}\frac{2^j}{K_{j-1}}
=O(n\beta_n^\ell).
\tag{17}
```

The series converges by the fixed schedule. Terminal assignments cost
`O(n beta_n^ell+n^{3/4})` bits. The top construction and grouping metadata
are `O(n/P+n^{3/4} polylog n)`, with `P=Omega(n^{1/40})` in this streamed
parameter choice. Seed tables are held once. The scalar registers in the
`O(loglog n)` active frames add only `polylog n` bits. Together with (5),
these give weighted-evaluator workspace

```math
O\!\left(n\beta_n^\ell+n^{39/40}
             +n^{3/4}\operatorname{polylog}n+n^{1/64}\right).
\tag{18}
```

The packed transfer charges the same terms to logarithmic running time:
(16) for piece choices, (17) for repeated local queries, and the explicit
leaf and top costs. Equivalently, after inlining the halting finite-tape
routines, the configuration-count argument in that transfer gives time at
most `poly(n,S) 2^{O(S)}` from this actual workspace bound. The one-way
write-only output does not affect future transitions and is omitted from
the internal configuration.
The termination needed here is established by the algorithms: recursion
strictly decreases the integer level, every local enumeration is finite,
and preprocessing has its explicit step budget. During recovery below,
successor scans strictly increase a divisor and therefore have at most
`tau(N)+1` passes. No output-dependent state is hidden on the output tape.

## 4. Recovery and the limiting exponent

Apply the pinned varying-characteristic recovery theorem in the selected
field `F_n`. In characteristic two, it uses only the degree-two coefficient
algebra `F_n[W]/(W^2+W+1)` for ternary characters. In odd characteristic,
binary characters already lie in the field. The same selected seed extends
to this algebra by its exact sum and product identities. All field and
extension representation costs are included in the polylogarithmic bounds.

That recovery enumerates short descriptors, regenerates candidate divisor
bits, and checks the factor pair by streamed multiplication. A modular count
of one is never accepted as proof of isolation. The elementary bound
`tau(N)<=2^{O(n/log n)}` supplies descriptors of `O(n/log n)` bits. It adds
`O(n/log n)` workspace and `2^{O(n/log n)} poly(n)` query overhead.

Put `alpha_n=log_2(1/beta_n)`. Then `alpha_n -> alpha_*` and the fixed
schedule gives

```math
\beta_n^\ell=\frac1{(\log n)^{\alpha_*-o(1)}}.
\tag{19}
```

All power-saving terms in (18) are smaller than `n/(log n)^t` for every
fixed `t`. Recovery also fits the claimed limiting bound because
`alpha_*<=1`. Thus both logarithmic time and workspace are
`n/(log n)^{alpha_*-o(1)}`. Any fixed multiplicative constants are absorbed
by an arbitrarily small positive logarithmic power, proving (3).

Arbitrarily slow `Q_n -> rho_*` is enough. A stronger statement with
`exp(O(n/(log n)^alpha_*))` time and `O(n/(log n)^alpha_*)` space does
not follow from this argument without quantitative convergence.

## 5. Rank-sum and all-field versions

Close a computable seed stream under every finite type amplification, using
dovetailing so every finite construction is eventually emitted. For one seed,
write `s(S)=sum_i sqrt(R_i C_i)`. The established finite-type estimate gives

```math
Q(S^{(k)})\le\frac{\log_2 s(S)}{2a}
       +\frac{2+2\lceil M\log_2(k+1)\rceil}{4ak}.
\tag{20}
```

Conversely `s(S)<=2^{w/2+c}`, so
`log_2(s(S))/(2a)<=Q(S)`. On an amplification-closed stream the two infima
are therefore equal. The theorem attains their common endpoint in both
time exponent and workspace, while selecting by the rational number (1).
No exact real logarithm or square-root comparison is an instruction.

For a universal enumeration over all explicit finite fields, the
[finite-field specialization theorem](https://github.com/1337isnot1337/integer-factoring/blob/023e4fc4482d72f25738f6b1d0f9cd31583baed1/research/finite-field-rank-specialization.md)
preserves every finite exact rank profile from an arbitrary coefficient
field. Hence the same endpoint includes all coefficient fields. This is an
existence and resource statement, not a computed numerical value of the
infimum. It strengthens the workspace conclusion of seed selection; it does
not improve the separate universal verified-search principle for time alone.

## 6. Exact finite checks

Run:

```sh
python3 tools/check_simultaneous_seed_endpoint.py \
  --output evidence/simultaneous-seed-endpoint.json
```

The checker constructs small exact crossing seeds, computes their ranks and
factorizations independently, checks both direct type-word factorizations,
tests mixed-radix word/channel unranking and tails, and verifies bounded
integer multinomial arithmetic. Delayed-selection and scalar-convergence
checks cover only those interfaces, not the existence of better seeds.
The infinite-input theorem follows from the proof and pinned dependencies,
not from finite rank extrapolation. This note and checker used substantial
OpenAI assistant support.
