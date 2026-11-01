# Streamed modular factors in every fixed characteristic

## Result

For every fixed prime `p`, the original normalized quotient crossing
family admits scalar-factor queries in `2^O(lambda)` bit time and
`poly(lambda)` space. The constants may depend on `p`.
This extends [the characteristic-two routine](modular-factor-space.md)
to repeated parts and odd characteristic. In particular it supplies the
local interface at `p=3`.

The underlying family is from
[modular-crossing.md](modular-crossing.md): retained state partitions
have no part divisible by `p`, and every multiplicity is less than `p`.
The original dual polynomial is retained as an operator parameter.
The scalar normalization is
`b^mu=z_mu^(-1)t^mu`, where `z_mu` is a unit modulo `p`.
The proof below does not divide by a factorial divisible by `p`.

This is a local time-space theorem. Its application to factoring also
needs a streamed grid evaluator and a nonzeroness algorithm over suitable
coefficient rings. It supplies no new contraction exponent by itself.

## 1. Repeated columns and differentiation have different factors

Let an operator monomial have `m` parts. Write the input multiplicities
as `sigma_k` and the output multiplicities as `pi_k`.
Choose unaffected multiplicities

```math
0\le\eta_k\le\min(\sigma_k,\pi_k),\qquad
2\sum_k\eta_k=|\sigma|+|\pi|-m,
\tag{1}
```

where the bars here denote particle counts. Put

```math
a_k=\sigma_k-\eta_k,\qquad c_k=\pi_k-\eta_k.
```

Make `a_k` identical annihilation columns and `c_k` identical creation
columns of mode `k`. There are exactly `m` signed columns. For operator
part `j`, their scalar entries are

```math
C(j,k)=\binom{k-1}{j-1},\qquad
D(j,k)=(-1)^{j-1}\binom{k+j-1}{j-1},
\tag{2}
```

with `C(j,k)=0` when `j>k`. The mode multiplier `k` is deliberately
excluded from `D` and included in the weight below.

The source normal-order coefficient is exactly

```math
\sum_\eta w_\eta\operatorname{per}M_\eta,
\qquad
w_\eta=\prod_k
 \frac{\sigma_k!\,k^{a_k}}
 {\eta_k!\,a_k!\,c_k!}.
\tag{3}
```

To prove (3), first assign each operator position to its signed mode.
Normal ordering then gives the derivative multiplier
`prod_k k^(a_k) sigma_k!/eta_k!`.
The permanent enumerates bijections to individually named columns.
It repeats each assignment exactly `prod_k a_k! c_k!` times.
Dividing by this overcount gives (3). Identical operator rows are
retained: their assignments already appear with the multiplicities in
the source's ordered recursion. No additional operator factorial is
introduced.

All of `sigma_k,pi_k,eta_k,a_k,c_k` are less than `p`.
Every denominator in (3) is therefore a unit in `F_p`.
This is the reason the division remains valid in odd characteristic.
There is no cancellation of a formal `p!` or division of zero by zero.

For the bottom factor, multiply (3) by `z_mu^(-1)`. This is the only
family normalization: the top dual remains the original `aFam(mu)`.

## 2. Fold the dual polynomial and every angle coefficient at once

For a label `mu` of weight `h`, its top dual is

```math
a^\mu=\prod_{n\in\mu}
 (-1)^{n-1}\sum_{j=1}^{n}\binom nj t_j.
\tag{4}
```

Do not enumerate its monomials. For each original operator row `n`
and signed mode `k`, use these polynomial entries instead:

```math
\begin{aligned}
F^+_{n,k}(S)
 &=(-1)^{n-1}\sum_{j=1}^{n}
   \binom nj C(j,k)(1+S)^j,\\
F^-_{n,k}(S)
 &=(-1)^{n-1}\sum_{j=1}^{n}
   \binom nj D(j,k)(1+S)^j.
\end{aligned}
\tag{5}
```

The full top-factor profile is

```math
(1+S)^h\sum_\eta w_\eta\operatorname{per}F_\eta(S).
\tag{6}
```

Its coefficient of `S^i` is the source angle-`i` factor. Indeed, for
each term of (4), the row factors contribute
`(1+S)^(sum j)`, and the additional factor in (6) gives
`binom(h+sum j,i)`. The permanent's multilinearity accounts for all
lowered operator choices, including repeated lowered parts.

The calculation keeps every `j` in (5). Terms whose operator monomial
is excluded from the state quotient are evaluated by the same complete
normal-order formula; there is no unsupported reduction of the operator
parameter. Thus only the already-established state-space quotient is
needed for this proof.

Each row polynomial has degree at most `n`; the permanent has degree
at most `h`, and (6) has degree at most `2h`. All these bounds are
independent of the number of monomials in the expanded dual.
The separate sign `(-1)^i` in a complete crossing piece is applied
after extracting the angle coefficient.

## 3. Stream all subset and multiplicity loops

Over `F_p[S]`, Ryser's formula is

```math
\operatorname{per}M
=\sum_{J\subseteq\{1,\ldots,m\}}
 (-1)^{m-|J|}\prod_{r=1}^{m}\sum_{j\in J}M_{rj}.
\tag{7}
```

Enumerate `J` sequentially with an `m`-bit counter. Keep one product,
one row sum, and one accumulator, truncated at the known degree bound.
Enumerate the unaffected multiplicities in (1) with a mixed-radix
counter, and discard the current permanent before proceeding to the
next choice. No table indexed by all subsets is needed.

With multiplicity at most `p-1`, a partition with `r` parts has weight
at least `r^2/(2(p-1))`. Thus the operator label of weight `h` has
`O_p(sqrt(h))` parts, and the input and output states of weight at most
`4h` have the same bound. There are at most
`p^O_p(sqrt(h))` unaffected-multiplicity choices and `2^O_p(sqrt(h))`
Ryser subsets per choice. The polynomial matrix (5), the counters,
the factorial units, and the coefficient arrays occupy `poly_p(h)`
space. Direct finite-field polynomial arithmetic gives

```math
\text{time}=2^{O_p(\sqrt h)}\operatorname{poly}_p(h),
\qquad
\text{space}=\operatorname{poly}_p(h).
\tag{8}
```

For the effective partition family `h=O_p(lambda^2)`, this is the
claimed `2^O(lambda)` time and polynomial space.
The exact half endpoint in
[rough-input-transfer.md, Section 6](rough-input-transfer.md) also
works for every fixed prime. Replace `A_2` there by
`A_p=(pi^2/6)(1-1/p)^2` in the same moment and cumulative-count bounds.
This gives width sum `2sqrt(2)lambda+O_p(lambda^(3/4))`, so the enlarged
transfer applies at the exact value `rho=1/sqrt(2)`.
The retained partition counts and unranking use a polynomial-size
bounded-multiplicity dynamic program through weight `4h`, exactly as
in the characteristic-two case. Its entries have `O_p(sqrt h)` bits.

An entire complete piece can also be evaluated in polynomial space by
summing its channel factorization sequentially. There are only
`2^O(lambda)` channel indices, and each factor query has the bound (8).
This preserves the same exponential-time order and supplies the direct
identity factorization used to cap an individual width at `2lambda`.

All factor values lie in `F_p`. For an ambient commutative `F_p`-algebra,
the local routine returns the same prime-field scalar for the circuit to
embed. The routine does not materialize a table over the coefficient ring.

## 4. Reproduction and controls

Run:

```sh
python3 tools/check_modular_factor_space_odd.py \
  --output evidence/modular-factor-space-odd.json
```

The checker compares the streamed formula with the independent literal
source recursion at primes three and five, including repeated input,
output, and operator modes. It checks every angle coefficient in its
finite cases and verifies complete normalized crossing sums in small
weights. Negative controls separately detect missing repeated-column
division, missing derivative multiplicity, and an omitted dual sign.
These checks verify the formula's implementation; the uniform space
bound is the streaming proof above.

The source operators are prior work by Anthropic at
[`16a5ce4`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
This extension and its checks used substantial OpenAI assistant support.
No new Lean formalization is claimed.
