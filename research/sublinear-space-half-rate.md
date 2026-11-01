# Half-rate factoring with sublinear working space

## Theorem and machine model

For every `n`-bit integer `N>=2`, a deterministic finite-alphabet machine
computes its least prime factor with

```math
\boxed{\quad
T(n)\le C\exp\!\left(\frac{Cn}{\sqrt{\log(n+2)}}\right),
\qquad
S(n)\le\frac{Cn}{\sqrt{\log(n+2)}}.
\quad}                                                        \tag{1}
```

The input is on a read-only tape and the output is streamed to a
write-only tape. `S(n)` counts working bits, so this is a sublinear
**workspace** bound, not a claim that an `n`-bit input takes fewer than
`n` bits to represent. The same bounds support complete prime
factorization, exact divisor-prefix counts, and sorted divisor output.
Prefix inputs have `0<=k<=n` and `0<=x<2^(n-k)` and are also read-only.

The time exponent remains one half. This note sharpens
[the linear-space theorem](linear-space-half-rate.md) by avoiding
materialized candidate integers and by charging the actual lengths of
the grid's local counters. It uses the original polynomial-piece regular
crossing family, its [linear-space factor queries](near-linear-factor-space.md),
and [streamed top queries](streamed-top-space.md).

The recovery records both factor inputs of the multiplication network.
A column-by-column multiplication check verifies their virtual bit
streams with logarithmic arithmetic workspace. No new Lean
formalization is asserted.

Write `F(n)=n/sqrt(log(n+2))` for the rest of the note.

## 1. A uniform short description of every divisor bucket

We need a divisor bound on every input, without first storing a
factored cofactor. For `n>=4`, put

```math
L=\lfloor\log_2 n\rfloor,\qquad s=\lfloor L/2\rfloor.
```

The elementary small/large-prime split from the linear-space theorem
gives

```math
\log_2\tau(N)\le2^s\log_2(n+1)+n/s.
```

Here `s>=L/3`, `2^s<=sqrt(n)`, and `log_2(n+1)<=L+1`. The exact
inequality

```math
L(L+1)\le6\,2^{L/2}\quad(L\ge2)                       \tag{2}
```

shows that the first term is at most `6n/L`; the second is at most
`3n/L`. To check (2), square it at `L=2,...,6`. Increasing `L` by
two multiplies its left side by `(L+2)(L+3)/(L(L+1))<=2` for
`L>=5`, matching the right-side factor. This proves all remaining
integers by parity induction.

Consequently the computable integer

```math
q=\min\left(n,\left\lceil\frac{9n}{L}\right\rceil\right),
\qquad K=2^q                                             \tag{3}
```

is a valid divisor-count bound. The alternative `q=n` is valid because
`tau(N)<=N<2^n`. Handle `n<4` directly. We have
`q=O(n/log(n+2))`.

Use the ternary evaluation hashes from the linear-space theorem, with
the least `k>=1` such that `Q=3^k>n(K-1)`. Thus

```math
k=O(n/\log(n+2)),\qquad Q=2^{O(n/\log(n+2))}.            \tag{4}
```

The fixed sixteen-element coefficient ring is unchanged. Tag bits of
both input factors, one at a time as specified in Section 3. For each
`t,b in F_(3^k)`, the character sums define a scalar `s(t,b)` and two
possible `n`-bit words. Validate that the scalar is one and all `2n`
coordinates are zero or one, by sequential queries. If validation
succeeds, denote the words by `d(t,b),e(t,b)`. They are **virtual
integers**: store only the pair `(t,b)`, and recompute any requested bit.

The control-field modulus and a constant number of such descriptors
use `O(k)` bits. Each bit query costs `Q` weighted evaluator calls and
polynomial extra work. Every actual divisor of `N` is represented by
at least one valid singleton-bucket descriptor, together with its
cofactor. A valid-looking colliding bucket may represent a false pair,
so arithmetic verification is still required.

## 2. The weighted grid evaluator uses only `O(F(n))` working bits

Choose the original top schedule of
[polynomial-top-queries.md](polynomial-top-queries.md), with largest
`P=K_ell<=A` and the `n^(3/4)` early-leaf cutoff. At every level that
performs a grid step,

```math
\lambda_j\ge n^{1/4}/8,
\qquad V_j=K_j^2\lambda_j
\le D\rho^{\ell-j}V_\ell,
\quad V_\ell=O(n),\quad \rho=1/\sqrt2.                 \tag{5}
```

The established parent-query series gives

```math
\sum_{j=1}^{\ell}\frac{V_j}{K_{j-1}}
=O(n\rho^\ell)=O(F(n)).                                \tag{6}
```

We now use the actual descriptor sizes, instead of bounding every item
at a level by the full `V_j`.

### Piece and grouping metadata

There are `O(K_j^2)` crossing sites. The original regular family has
polynomially many pieces and `O(lambda_j)`-bit channels, so a piece
index and its widths use `O(log(lambda_j+2))` bits. Equations (5) and
the geometric sum of volumes give

```math
\sum_j O\!\left(K_j^2\log(\lambda_j+2)\right)
=O(n^{3/4}\log(n+2))=o(F(n)).                           \tag{7}
```

Grouping codes, crossing positions, and group-width records fit within
`O(P^2 log(n+2) log log(n+2))`, which is also `o(F(n))` since
`P=O(n^(1/4))`. They can instead be generated from indices. Neither
implementation stores an array of ring elements.

### Local assignments and query arguments

A grouped-table query at a grid step enumerates assignments to
`O(V_j/K_(j-1))` binary local-setting coordinates in the source grouping
construction. Its one current setting, and the arguments of one
parent-node query, use `O(V_j/K_(j-1))` bits. Products visit their factors sequentially;
there is no cache of all vertex arguments.

There is one structural frame per level and at most one additional
upward table-query frame at that level. Equation (6) bounds the total
length of all these local counters and arguments by `O(F(n))`.
Use packed strings and formula-derived offsets. All ring accumulators
have constant size because the coefficient ring has sixteen elements.

The local crossing routines require `O(lambda_j)` scratch bits. Since
`lambda_j=V_j/K_j^2<=V_j/K_(j-1)`, this fits the same bound. Queries are
sequential, so only one local scratch workspace is active. A bottom
leaf stores one full edge assignment of volume `O(n rho^ell)`; an early
leaf stores at most `n^(3/4)` bits. Both are `O(F(n))`.

### The top query

The [streamed top-query construction](streamed-top-space.md) uses

```math
O(n/P+Af+\log(n+2))
```

working bits, apart from the supplied hash-control registers. It stores
the current boundary and one short interpolation polynomial. It
recomputes input, middle, and inverse outputs through the constant
number of transform layers, and handles carries in a forward loop.

Here `Af=O(n^(3/4))`. The largest scheduled `P<=A` satisfies
`A<K_(ell+1)<=P^10`, so `P=Omega(n^(1/40))` and
`n/P=O(n^(39/40))`. Both quantities are `o(F(n))`. The on-demand
hash phase adds `O(k+log n)=o(F(n))` bits by (4).

The entire weighted evaluation therefore uses `O(F(n))` workspace.
Its established time `exp(O(F(n)))` is unchanged: local factor time
is still `2^O(lambda)`, top recomputation is polynomial, and the
parent-query costs were already charged by (6).

## 3. Verify both virtual factors by streamed multiplication

The settings-to-divisors bijection fixes both factor inputs. Give their
bits separate variable families. The weighted network then computes

```math
W_N(z,w)=\sum_{ab=N}\prod_{j=0}^{n-1}z_j^{a_j}w_j^{b_j}. \tag{8}
```

This changes only input-gate weights. It preserves the graph, all local
crossings, and the unique-extension top queries. The hash phases still
depend only on the first input: set
`z_j=omega^(v dot coords(t^j))` and `w_j=1` before adding a tag.

To query bit `j` of the first factor, multiply `z_j` by `1+epsilon`.
To query bit `j` of the second factor, instead multiply `w_j` by
`1+epsilon`. Apply the same character projection to bucket `b`. If that
bucket isolates a divisor `a`, the two returned bit streams are exactly
`a` and `N/a`. The coefficient ring still has sixteen elements. Each
bit requires a separate query, so no `2n`-coordinate algebra is stored.

For possibly colliding buckets, verify the proposed pair by multiplication.
Here is an elementary small-space multiplication comparator for any
three virtual `n`-bit integers `a,b,c`. Set `carry=0`. For
`k=0,...,2n-1`, compute

```math
u_k=\operatorname{carry}
       +\sum_{\substack{0\le i<n\\0\le k-i<n}}a_i b_{k-i}.
```

Compare `u_k mod 2` with `c_k`, where `c_k=0` for `k>=n`, and then set
`carry=floor(u_k/2)`. Require a final zero carry. Each column has at
most `n` summands. Inductively `carry<=n`, so the carry, column sum,
and indices use `O(log(n+2))` bits. The check uses `O(n^2)` bit queries
and is exactly ordinary binary multiplication with its columns streamed.

For bucket validation take `c=N`, read directly from the input. Accept
the descriptor only if `d(t,b)e(t,b)=N`. Every accepted first coordinate
is then a positive true divisor; every true divisor is accepted at a
singleton bucket. Compare candidates or test a read-only prefix by
scanning their requested bits. The same comparator also tests products
between three stored bucket descriptors for Section 5.

All bit requests ultimately call the weighted evaluator on the original
read-only `N`. They are made sequentially while retaining only the
logarithmic multiplication or comparison state and the constant number
of descriptors. Verification therefore has time `exp(O(F(n)))` and
workspace `O(F(n))`. It does not call a factoring routine on a virtual
new input and does not create a recursive verification dependency.

## 4. Least-factor recovery and exact enumeration

Scan all `Q^2` bucket descriptors. Validate their coordinates, verify
their represented integers, and retain a descriptor of the smallest
verified divisor greater than one. Only the best and current descriptors
are needed. The divisor `N` itself is always isolated somewhere, so a
best descriptor always exists. Since every true divisor is isolated,
its value is the least prime factor, including when `N` is prime.

Output its bits by querying the saved descriptor, from most to least
significant, suppressing leading zeros. No `n`-bit candidate register
is needed. Equations (4)--(6) and the multiplication comparator prove
(1) for least-factor output. The `poly(n,K)` number of calls adds
`exp(O(n/log(n+2)))` time, which is absorbed by `exp(O(F(n)))`.

For sorted divisor enumeration, retain the last output descriptor and
repeat the scan for the smallest verified divisor strictly greater
than it. Start with a zero sentinel and include one. There are at most
`tau(N)+1<=K+1` scans. This is still only a polynomial factor in `K`.
It handles duplicate hash representations without storing a set.

For an exact prefix or interval count, test each newly selected divisor
against the read-only endpoint data and increment a counter. The counter
has at most `q+1=O(n/log(n+2))` bits. No divisor list or factorization
list is retained. The same procedure can stream all divisors instead.

## 5. Complete prime factorization without a cofactor register

One may also use the same immutable input `N` throughout complete
factorization. Visit its divisors in increasing order using Section 4.
A divisor `d>1` is prime exactly when it is not the product of two
nonunit divisors `a,c` of `N`. If `d` is composite, choose any proper
factorization `d=ac`; both factors divide `N`, so both have accepted
descriptors.

For each `d`, scan pairs of bucket descriptors, validate their factor
pairs as in Section 3, require `a,c>1`, and compare their product with
`d` by the streamed comparator. If a pair passes, skip `d`; otherwise
`d` is prime. There are only a constant number of nested scans and live
descriptors. Their storage is `O(k+log n)` and their number is
`poly(n,K)`.

For a prime `d`, retain a descriptor for the current power, initially
`d` itself, and set the exponent to one. Scan accepted divisor handles
for an integer equal to `current_power*d`, using the same comparator.
If one is found, retain its descriptor and increment the exponent.
If none is found, the next power does not divide `N`, because every
divisor is represented. The loop has at most `n` successful powers.
It therefore returns the exact prime multiplicity without constructing
the power as an integer. Output the prime and exponent, or stream that
prime repeatedly.

This procedure stores neither the cofactor nor a prime-factor list.
All arithmetic input bits ultimately come from the original read-only
`N` and a constant number of bucket descriptors. Its time is
`poly(n,K)` weighted-evaluator calls and its workspace is `O(F(n))`.
This proves the complete-factorization claim.

## Verification and limits

The companion checks cover the exact bound (2), finite divisor-count
bounds, bitwise descriptor comparisons, verification, duplicate removal,
sorted output, prefix counts, and prime/exponent extraction. The top
query and local factor checkers test the two new streaming interfaces.

The full enormous evaluator is not materialized in the finite tests. The
claims in (1) are ordinary asymptotic proofs supported by exact finite
checks and independent internal review; they are not new Lean results
or a practical factoring benchmark.

```sh
python tools/check_sublinear_space_recovery.py \
  --output evidence/sublinear-space-recovery.json
```
