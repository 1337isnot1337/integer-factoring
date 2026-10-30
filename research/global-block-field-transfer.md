# One global block length permits growing coefficient fields

## Theorem and exact interface

The coefficient field of a crossing construction may vary with its label
length, and its degree may be exponential in a fixed power of that length. The evaluator
can choose one raw block length from the original input length and use
that same length throughout its recursion.

Fix `1/2 <= rho < 1`, an integer `h>=1`, and fixed positive constants.
For every positive integer `t`, assume uniform data with these properties:

- A monic irreducible polynomial `f_t in Z[X]` describes
  `K_t=Q[X]/(f_t)`, of degree `e_t>=1`, in the coordinate basis
  `1,X,...,X^(e_t-1)`.
- Polynomially many pieces sum to a nonzero scalar times the crossing
  on `t`-bit labels. Their displayed factorizations have separate binary
  widths at most `2t` and total width at most
  `4rho t+C t/log(t+2)`.
- The polynomial, scale, complete piece enumeration, widths, and any
  requested scalar factor entry are generated in `2^{O(t^h)}` **bit time**.
  Field elements are explicit lists of rational coordinates with binary
  numerators and nonzero binary denominators. All constants are uniform
  in `t`.

Then, with `alpha=log_2(1/rho)`, one deterministic finite-alphabet
machine factors every `n`-bit integer `N>=2` with

```math
T(n)=\exp\!\left(O\!\left(\frac n{(\log(n+2))^\alpha}\right)\right),
\qquad
S(n)=O\!\left(\frac n{(\log(n+2))^\alpha}\right).        \tag{1}
```

Input is read-only and output is streamed. The conclusion includes the
least prime factor, complete factorization, sorted divisors, and exact
prefix or interval counts with read-only endpoints, including `alpha=1`.

The representation permits `e_t=2^{O(t^h)}` and coefficient heights with
`2^{O(t^h)}` bits. A cost bound with any fixed real power at least one
is covered by rounding that power upward to an integer `h`. The width
condition includes every fixed power-saving error `O(t^theta)` with
`theta<1`. No fixed number field or polynomial degree bound is
assumed. An unspecified algebraic oracle does not meet the interface.
This is a conditional transfer: a new family with a smaller `rho`
remains necessary for an unconditional score above one half.

The proof extends the fixed-field case of the
[rational simultaneous transfer](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/rational-simultaneous-transfer.md)
at `1476525`, using the existing
[packed grid and type grouping](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/simultaneous-time-space-transfer.md).
The new argument chooses a global raw length and restricts its tail.
No new Lean formalization is claimed.

## 1. Compile one raw family into one finite ring

At length `t`, normalize the crossing scale to one inside `K_t`.
The extended Euclidean algorithm modulo `f_t`, or rational linear
algebra on the multiplication matrix, computes its inverse in `2^{O(t^h)}`
bit time. To see the height bound explicitly, an `e_t` by `e_t`
integer matrix with `H`-bit entries has determinants and cofactors with
`O(e_t(H+log(e_t+1)))` bits. Clearing input denominators and using
fraction-free elimination therefore preserves the `2^{O(t^h)}` bit bound.

Enumerate every coordinate of every normalized entry in all four
factor tables. A table has at most `2^(4t)` entries, since it uses
two external `t`-bit labels and at most `2t` channel bits. Polynomial
piece count and the entry-generation bound give `2^{O(t^h)}` total output
bits. Let `D_t` be the product of all positive coordinate denominators.
Its bit length and construction time are `2^{O(t^h)}`. Repetitions are
allowed; complete piece entries are sums of factor products and need
no further denominator collection.

Use

```math
J_t=2D_t+1,\qquad
R_t=(\mathbb Z/J_t\mathbb Z)[X]/(f_t(X)).               \tag{2}
```

Every collected denominator is a unit modulo `J_t`. Monicity gives
a unique remainder of degree below `e_t`, so `R_t` is free of rank
`e_t` over `Z/J_t Z`. Constants embed faithfully, zero differs from
one, and two is a unit. Reducibility modulo `J_t`, composite moduli,
and zero divisors cause no difficulty.

Reduction from `Z[1/D_t][X]/(f_t)` is a ring homomorphism. The
normalized crossing and both factorizations hold in that localized
ring because their monic remainders are zero coordinatewise. They
therefore survive in (2). All normalization occurred before reduction;
no arbitrary algebraic element is inverted in the finite ring.

An element is an array of `e_t` residues. Direct convolution and monic
reduction implement arithmetic. Even storing the full multiplication
tensor, with `e_t^3` residues, and the entire normalized raw family
costs `2^{O(t^h)}` bits. The `2^(4t)` factor-table size and the
`e_t^3` multiplication-tensor size both fit this bound because `h>=1`.
Include all coefficient bit lengths and sequential
tape indexing. There is consequently one fixed integer `C_*>=1`
such that preprocessing, stored data, one embedded raw query, and one
ring operation have time and space at most

```math
G(t)=2^{C_*(t^h+1)}.                                   \tag{3}
```

This bound pays for the full field degree and multiplication data.
The supplied family fixes `C_*`; no algorithm to infer complexity
bounds from arbitrary code is required.

## 2. A single global length and its restricted tail

Fix an integer `B>=128C_*`. Let `iroot_h(v)` be the largest nonnegative
integer `u` satisfying `u^h<=v`, and choose

```math
L=\lceil\log_2(n+2)\rceil,\qquad
m=m(n)=\max\left(1,
 \operatorname{iroot}_h(\lfloor L/B\rfloor)\right).
```

Compile only this raw family, over only `R_m`. Integer bit length,
division, and a fixed-power integer-root search compute the choice.
In particular `m=Theta((log(n+2))^(1/h))`. Since
`m^h<=max(1,L/B)`, equation (3) gives

```math
G(m(n))=O((n+2)^{1/128}).                              \tag{4}
```

Use the explicit top clique and `n^(3/4)` early-leaf cutoff of the
packed grid. Every node that performs a crossing step satisfies

```math
\lambda_j\ge n^{1/4}/8.                               \tag{5}
```

For sufficiently large inputs this is greater than `m`. The finite
remaining input range can use the original exhaustive procedure.
At an active node write `lambda=bm+q`, with `0<=q<m`. Every full
block uses the length-`m` family.

For the tail, inject `q`-bit labels into `m`-bit labels by appending
zeros and restrict every raw tensor and both factorizations to these
labels. The sum remains exactly the `q`-bit crossing. Cap each
restricted cut at `2q` bits using the identity matrix and the
complete restricted piece matrix. A complete entry is a sum over at
most `2^(2m)` raw channels, so its cost is `2^{O(m^h)}`. The tail has
at most the same number of pieces and total width at most `4q`.
It calls no length-`q` field or family.

For a preliminary label envelope, the separate width caps give
`V_(j-1)<=(1+3zeta_(j-1))V_j+K_(j-1)^2<=4V_j`, because
`zeta<=1/16` and `K_(j-1)^2<=V_j`. Top volume is at most
`1296n`. The schedule has `log_2 K_ell>=18*2^ell`, while the
selected top satisfies `K_ell<=A` and `log_2 A<=n`. Thus
`4^ell<=n^2/324`, and every visited volume is at most `4n^3`.
In particular,

```math
\lambda_j\le\Lambda(n)=2^{16}(n+2)^3.                 \tag{6}
```

These bounds use only the width caps and fixed schedule. They do not
assume the sharper contraction proved below.

## 3. Type grouping and explicit indexing

Let the compiled raw pieces be `A_1,...,A_k`, with `k=O(m^a)`.
For a type `nu` with `sum nu_i=b` and one restricted tail piece
when `q>0`, use

```math
B_{\nu,a}=\left(\sum_{\operatorname{type}(w)=\nu}
A_{w_1}\otimes\cdots\otimes A_{w_b}\right)\otimes A_a^{(q)}.
                                                               \tag{7}
```

Omit the tail for `q=0`. Summing these tensors gives the crossing.
Each factorization channel chooses a type word and an original
channel at every block. Both cuts therefore pay the full multinomial
factor `M_nu=b!/prod nu_i!`; no word orders are discarded.

Since `M_nu<=k^b`, the total binary width satisfies

```math
r_{\nu,a}+c_{\nu,a}\le4\rho\lambda+
O\!\left(\frac{\lambda}{\log(m+2)}
 +\lambda\frac{\log(m+2)}m+m\right).                  \tag{8}
```

The last term pays for the tail and rounding. Cap each full cut at
`2lambda` with a complete-piece query if necessary. The number of
pieces is at most `binom(b+k-1,k-1)*max(1,k)`. Their entropy and
explicit type-descriptor length are bounded by

```math
f_n(\lambda)=O(k\log(b+1)+\log(k+1))
=O((\log(n+2))^{a/h+1}),                              \tag{9}
```

uniformly over (5)--(6). This local family depends on the original
input length through `m`, which is permitted by the direct grid
algorithm. We do not invoke a fixed local-family theorem outside its
hypotheses.

Compute multinomials as products of binomial coefficients, using the
exact multiplicative recurrence and the smaller lower argument for
each binomial. Every partial binomial and product is at most `M_nu`;
an intermediate multiplication adds only `O(log b)` bits. The
first-letter counts for word unranking are `M_nu*nu_i/b`. Continue
with the chosen count and decremented type. These are exact integer
indexing divisions, never divisions of coefficient-ring values.

Word, multinomial, and channel counters use `O(lambda)` bits, and
(9) bounds the smaller type vector. Factor queries visit blocks
sequentially. Complete-piece and capped-factor queries enumerate type
words and add their products. Their cost is at most `k^b` times a
polynomial in `n`, hence `2^O(lambda)` by (5)--(6). Complete
piece enumeration also fits this bound by (9). The globally compiled
data and a constant number of ring accumulators add only (4), which
fits the same local space bound at every executed crossing step.

## 4. Uniform contraction and complete grid accounting

The relative slack in (8), together with child-width rounding, is
uniformly bounded by

```math
\varepsilon_n=O\!\left(\frac1{\log\log n}
 +\frac{\log\log n}{(\log n)^{1/h}}
 +\frac{(\log n)^{1/h}}{n^{1/4}}+n^{-1/2}\right).    \tag{10}
```

There are `ell=O(loglog n)` levels, so `ell*epsilon_n=O(1)`.
The resulting product of relative errors is bounded by
`exp(O(ell*epsilon_n))=O(1)`. A limit of one is unnecessary.
The summable schedule errors and the direct grouping recurrence give

```math
V_j\le C\rho^{\ell-j}V_\ell,\qquad
\sum_j V_j=O(n),\qquad
\rho^\ell=\Theta((\log n)^{-\alpha}).                 \tag{11}
```

The piece-assignment entropy over one worst path is at most

```math
\sum_j O(K_j^2 f_n(\lambda_j))
 =O(n^{3/4}(\log n)^{a/h+1}),                         \tag{12}
```

by (5), (9), and (11). The same quantity bounds the stored piece
ordinals and widths. It is smaller than `n/(log n)^alpha` for
every fixed `alpha`. The established parent-query series is

```math
\sum_{j=1}^{\ell} V_j/K_{j-1}=O(n\rho^\ell).           \tag{13}
```

Local generation and arithmetic have the explicit polynomial time
cost (4), already absorbed into the fixed
`2^O(V_j/K_(j-1))` query cost. Early leaves have volume at most
`n^(3/4)`; final leaves have volume `O(n rho^ell)`. The source
grouping searches retain their smaller cost. Thus the weighted grid
has the time in (1).

Its packed assignments, upward query arguments, and structural data
have the space in (1). There are `O(loglog n)` structural/query
frames, with a constant number of ring registers each. Their total
scalar storage is `O(G(m(n))*loglog(n+2))`, with the compiled global
data held only once. The streamed top adds `O(n/P+n^(3/4))`
binary bits and the same small registers. All are lower order than
(1); no unspecified polynomial-space term is hidden in the claim.

The arithmetic part of this proof only needs effective scale-one
pieces over one nonzero finite commutative ring per raw length, with
the bit-cost bound (3). Its characteristic may vary with that length.
The number-field construction (2) additionally makes two invertible,
which is the property used by the following recovery argument.
For explicit finite fields with varying characteristic, the
[characteristic-independent recovery](characteristic-independent-recovery.md)
supplies the corresponding conclusion using a constant-degree coefficient
extension when necessary.

## 5. Recovery and direct exact counts

Use the paired-descriptor recovery from the rational simultaneous
transfer. Binary character phases `+1,-1` and normalization by
`2^k` work over (2). Hash the first factor's bits, tag one bit of
either factor at a time, and verify virtual candidate factors by
streamed binary multiplication. Zero and one are distinct constant
vectors, and a singleton bucket gives the true factor bits even if
the coefficient ring has zero divisors.

The uniform divisor bound gives `O(n/log(n+2))`-bit descriptors
and `exp(O(n/log(n+2)))` recovery overhead. This is absorbed in
(1) for `alpha<=1`. Successor scans and the same short descriptors
give sorted divisors and complete prime factorization.

For direct exact counts use the small-modulus generalized CRT method
from the
[direct-count theorem](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/direct-sublinear-divisor-counts.md).
For each eligible `j` with `gcd(j,2D_m)=1`, replace `J_m` in
(2) by `j`, retaining the same monic polynomial and normalized raw
family. The integer count is the constant coordinate of the returned
array, with every other coordinate zero.

Let `q_N=O(n/log(n+2))` be the proved count-bit bound and
`b_D=bitlength(D_m)<=G(m(n))`. The supply limit
`H=4(q_N+(b_D+1)^2+4)` is `O(n/log(n+2))`. Each call uses
`O(log n)`-bit residues per basis coefficient. Recompiling its
arithmetic data costs at most `G(m(n))*polylog(n)` bits and time,
which is `O(n^(1/64))` for large `n` by (4). Across active frames,
scalar storage remains lower order than (1). The outside CRT state
has `O(q_N)` bits and is held once. At most `H` queries and
polynomial extra work therefore give (1), including `alpha=1`,
without an extra `loglog n` workspace factor.

All queries use the same original input length and raw field. The
strongest unconditional exponent remains `1/2`; this result removes
a coefficient-field restriction from a route to a future improvement.

## Verification and source scope

The companion checker tests exact field normalization, both crossing
factorizations, composite moduli, zero-padded tail restrictions and
capped factors, type-word indexing, and explicit global bit budgets.
Its additional tensor-quotient and repeated-field zero-divisor controls
are not needed by the single-field proof.

```sh
python3 tools/check_growing_number_field_transfer.py \
  --output evidence/growing-number-field-transfer.json
```

Finite checks validate these interfaces; the asymptotic bounds follow
from Sections 1--5. They do not establish a smaller contraction or
execute the enormous asymptotic factorer.

The source evaluator is Claude/Anthropic's Appendix B development at
`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`; the multiplication
reduction is Isaac Rudich and Louis-Martin Rousseau's work at
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`. This extension was
prepared with substantial OpenAI assistant support. No new Lean
formalization is asserted.
