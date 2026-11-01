# Half-rate deterministic factoring in linear space

## Theorem and scope

There is a deterministic finite-alphabet machine that, on every `n`-bit
integer `N>=2`, returns its least prime factor with

```math
\boxed{\quad
T(n)\le C\exp\!\left(\frac{Cn}{\sqrt{\log(n+2)}}\right),
\qquad S(n)\le Cn.
\quad}                                                        \tag{1}
```

Here `S(n)` counts work-tape bits. The machine can also stream the
complete prime factorization or compute **exact divisor-prefix counts**
within the same bounds. Constants are independent of the input.
Finitely many small inputs use trial division.

This improves the earlier polynomial-space bound. It **does not**
improve the time exponent beyond one half. The new exact-count result
uses factorization followed by divisor enumeration. It does not assert
linear space for the earlier integral-Schur counting implementation.

The proof combines four explicit ingredients:

- The existing half-rate crossing and streamed grid evaluator.
- [Linear-space local factors](near-linear-factor-space.md), including
  label and channel unranking, widths, and capped piece queries.
- The [unique-extension top constructor](polynomial-top-queries.md).
- The deterministic isolation proof in
  [bit-tagged witness recovery](bit-tagged-witness-recovery.md), evaluated
  one tag coordinate at a time.

Only the last three ingredients require sharper storage accounting.
The source operators and recursion remain the prior constructions
identified in those notes. This is an ordinary mathematical proof,
without new Lean formalization.

## 1. Recover one bit at a time over a fixed sixteen-element ring

Let `E=F_4`, with an element `omega` satisfying `omega^2+omega+1=0`.
Use the fixed coefficient ring

```math
R=E[\epsilon]/(\epsilon^2).
```

An element is a pair `(a,b)` of elements of `E`, with

```math
(a,b)+(c,d)=(a+c,b+d),\qquad
(a,b)(c,d)=(ac,ad+bc).                                  \tag{2}
```

There are exactly sixteen elements. Ring operations can therefore be
implemented by fixed tables and require constant space. The local
crossing factors are bits, so their extension to `R` is immediate.

For a witness set `S subset {0,1}^n`, the weighted network computes

```math
C_S(z_0,\ldots,z_{n-1})
=\sum_{d\in S}\prod_{l=0}^{n-1}z_l^{d_l}.                \tag{3}
```

Fix a target bit `j` and scalar phases `beta_l in E`. Substitute

```math
z_l=\beta_l(1+[l=j]\epsilon).
```

The returned pair is exactly

```math
\left(\sum_{d\in S}\prod_l\beta_l^{d_l},\quad
       \sum_{d\in S}d_j\prod_l\beta_l^{d_l}\right).       \tag{4}
```

This follows by expanding each monomial: its sole possible tag is the
one at bit `j`. Running (4) separately for `j=0,...,n-1` recovers the
same coordinates as the former `(n+1)`-dimensional tag algebra. It
costs a factor of at most `n` in time. The evaluator now has only
constant-size coefficient-ring elements, including at every recursive
level. Keep the decoded candidate in an `n`-bit register.

## 2. Stream the isolation hash and its phases

After removing small prime factors, suppose the remaining divisor set
has size at most `K=2^q`. Choose the least positive integer `k` with

```math
Q=3^k>n(K-1),\qquad H=\mathbb F_{3^k}.                  \tag{5}
```

Then `k=O(q+log(n+2))`. The control field `H` is used only to generate
the phases. The evaluator's coefficient ring remains the fixed ring
`R` of Section 1.

For each `t in H`, define

```math
H_t(d)=\sum_{l=0}^{n-1}d_lt^l\in H.
```

For each fixed witness `d`, fewer than `n(K-1)` choices of `t` can
collide with another witness. Indeed, each difference polynomial is
nonzero and has degree at most `n-1`. Equation (5) therefore guarantees
that some hash isolates each witness.

Use the polynomial-basis coordinate vectors in `F_3^k`. For a bucket
`b in H`, character `v in F_3^k`, and bit query `j`, set

```math
\beta_l=\omega^{v\cdot\operatorname{coords}(t^l)}.
```

Multiply the result of (4) by
`omega^(-v dot coords(b))` and sum over `v`. The usual coordinatewise
geometric sums give

```math
\sum_{v\in\mathbb F_3^k}
\omega^{v\cdot(\operatorname{coords}(a)-\operatorname{coords}(b))}
=[a=b].                                                    \tag{6}
```

The unnormalized scalar is `3^k`, which equals one in characteristic
two, so no extra factor is missing. The resulting pair is

```math
\left(\sum_{d\in S:H_t(d)=b}1,\quad
       \sum_{d\in S:H_t(d)=b}d_j\right).                  \tag{7}
```

A singleton bucket returns scalar one and the witness's `j`th bit.
Combine the bit queries for this same `t,b` into one candidate. Require
scalar one and zero-one bit coordinates, then verify the candidate by
integer division and endpoint checks. A colliding bucket may produce a
candidate, but verification prevents a false factor. Every witness
appears correctly in some singleton bucket.

### Storage and time of the control computations

Construct `H` by enumerating monic degree-`k` ternary polynomials and
testing division by every monic polynomial of degree at most `k/2`.
Keep only the current candidate, divisor, and polynomial remainder.
All have `O(k)` coefficients in a fixed alphabet. The first irreducible
polynomial gives the field. Construction costs `3^O(k)` time and `O(k)`
space. Polynomial-basis field addition, multiplication, division, and
powering also have implementations with `O(k+log(n+2))` storage.

Generate a requested phase by computing `t^l` with binary powering.
Do not store the list of all `n` powers. Keep `t,b,v`, the field modulus,
the current bit index, and a constant number of field registers. These
use `O(k+log(n+2))` bits. Sequentially scan all `t,b,v` and all target
bits. There are at most `n Q^3` evaluator calls. All extra work costs

```math
2^{O(q+\log(n+2))}
```

times the evaluator's cost and uses `O(n+q)` storage, including the
candidate and least-candidate registers.

For factoring, apply this to the proper-divisor polynomial: subtract
the known monomials for one and `N` from the network's divisor
polynomial. Each endpoint monomial is evaluated by scanning its bits
and multiplying the on-demand weights. Keep the smallest verified
proper divisor over the whole hash scan. If one exists, it is the least
prime factor because every proper divisor is eventually isolated. If
none exists, `N` is prime and the answer is `N`.

## 3. The original top query needs only linear storage

Use the original arithmetic parameters and the explicit grouping in
[polynomial-top-queries.md](polynomial-top-queries.md):

```math
A=2S,\quad h=S,\quad d=2S^2,\quad f=4d+1,
\qquad A^2f\le144n.
```

The arithmetic modulus has `f=O(sqrt(n))` bits. Across all original
gates, there are only `O(A^2)` stored arithmetic words:

- `O(A)` input gates, each with `O(A)` digits and outputs.
- `O(A)` middle and inverse gates, each with `O(A)` incoming and
  outgoing words.
- `O(A^2)` carry gates with a constant number of words each.

Thus even storing all recovered local settings takes `O(A^2 f)=O(n)`
bits. An `h`-by-`h` interpolation workspace uses another `O(h^2 f)=O(n)`
bits. Ordinary modular arithmetic, Euclidean inversion, and polynomial
interpolation use a constant number of such workspaces. The unit
denominators are established in the top-query proof, so this is valid
over its composite arithmetic modulus.

The gate metadata have `O(A^2)` records with `O(log(A+2))`-bit indices.
They also fit in `O(n)` bits, or can be generated from the gate formulas.
One queried boundary has at most the top information volume, which is
`O(n)`. Padding and local checks are streamed.

After determining the unique valid setting, scan its first-factor bits
and multiply the weights from Section 2. Each weight is generated on
demand in `O(k+log n)` additional space, and its value is in the fixed
ring `R`. The resulting top-table query has polynomial bit time and
`O(n)` workspace. No table is materialized.

## 4. Linear storage for the full grid recursion

Use the larger top schedule and the `n^(3/4)` volume cutoff proved in
the top-query note. At level `j`, write

```math
V_j=K_j^2\lambda_j.
```

The existing half-rate analysis gives a fixed `D` with

```math
V_j\le D\rho^{\ell-j}V_\ell,
\qquad V_\ell=O(n),\quad \rho=1/\sqrt2,
\qquad \ell=O(\log\log(n+2)).                            \tag{8}
```

In particular, the **sum** of the volumes on any live ancestor path is
`O(n)`, not merely each separate volume. Small width errors and rounding
are already included in the constant `D`.

At one active level store only the following data:

| Data | Bit bound |
| --- | ---: |
| Piece choices and their two widths | `O(K_j^2 log(lambda_j+2))=O(V_j)` |
| Current grouping, orientation, and grouped-width counters | `O(V_j)` |
| Current edge labels and local-setting indices | `O(V_j)` |
| Table-query boundary and loop positions | `O(V_j)` |
| Ring sum and product accumulators | Constant per loop frame |

For the first row, the characteristic-two family has polynomially many
pieces and each width is `O(lambda_j)`. Use packed label strings and
compute offsets from dimensions, instead of allocating a pointer for
each bit. Grouping codes and grouped-width arrays fit within the same
bound. For example, a grouped edge-width array has at most
`K_(j-1)^2` entries of `O(log(V_j+2))` bits, while the square-growth
schedule has `K_(j-1)^2<=K_j`. Thus its size is `O(V_j)`.

A structural recursive call has one frame per level. While answering
a table query, there is at most one additional parent-table query frame
at each ancestor level: each product queries its factors sequentially.
Consequently (8) bounds the stored descriptors by `O(n)`. Return
positions and the `O(log log n)` constant-size ring accumulators are
smaller than that bound.

A local crossing-factor query at level `j` uses `O(lambda_j)` bits by
the [local factor theorem](near-linear-factor-space.md). Queries are
sequential, so only one such scratch workspace is live. It is `O(n)`.
At a base or early leaf, enumerate one complete edge assignment at a
time. Its length is at most its information volume `V_j=O(n)`, and its
ring accumulators have constant size. The top query uses `O(n)` space
by Section 3. All these workspaces can be reused.

This proves that each weighted evaluator call over `R` uses `O(n)`
bits. The time is still the established half-rate time. The parent-query
cost was already charged by the convergent series

```math
\sum_{j=1}^{\ell}\frac{V_j}{K_{j-1}}
=O(n\rho^\ell)=O\!\left(\frac n{\sqrt{\log(n+2)}}\right).
```

The sharper local factors retain their original `2^O(lambda_j)` time
bound. Recomputing hash powers and using packed sequential-tape data
adds only polynomial overhead. No memoized circuit values or dense
factor tables are assumed.

## 5. Complete factoring time and storage

Take `q=ceil(sqrt(n))` and `s=ceil(n/q)`. Ascending trial division
through `2^s` takes `exp(O(sqrt(n)))` time and `O(n)` storage. Every
factor discovered in ascending order is prime after smaller factors
have been removed. The remaining cofactor has at most `2^q` divisors,
by the previously proved rough-input estimate.

Apply Sections 1--4 to that cofactor, retaining the original `n` as a
uniform bound. Field construction, character scans, separate bit queries,
and verification multiply the evaluator cost by `exp(O(sqrt(n)))`.
This is absorbed by the half-rate exponent in (1). Their storage is
`O(n+q)=O(n)`.

For a complete factorization, return the least prime, divide it out,
and repeat on the cofactor. There are at most `n` calls. Stream the
prime factors to the output tape; there is no list of pending recursive
cofactors or all recovered witnesses. The time bound and linear working
space are unchanged. This proves (1).

## 6. Exact divisor-prefix counts in the same linear space

Write the prime factorization as `N=prod_i p_i^{e_i}`. Store the distinct
primes and positive exponents in a self-delimiting representation.
The total description length is `O(n)` because

```math
\sum_i\log_2p_i\le\log_2N<n,
\qquad
\sum_i\log_2(e_i+1)\le\sum_i e_i\le\log_2N<n.
```

Length delimiters increase this bound by only a constant factor. One
may generate this list by the repeated least-factor procedure in
Section 5, dividing out each prime fully before recording its exponent.

Enumerate every vector `(a_i)` with `0<=a_i<=e_i` in mixed radix.
It represents the unique divisor `d=prod_i p_i^{a_i}`. The exponent
vector has `O(n)` bits, and the current product has at most `n` bits.
For a prefix query with `0<=k<=n` low bits free and high prefix
`0<=x<2^(n-k)`, increment
an integer counter exactly when `floor(d/2^k)=x`. Include or exclude
the divisor one according to the requested convention. Each divisor
appears exactly once by unique prime factorization.

The number `tau(N)` of divisors is small enough for this enumeration.
Here is an elementary bound sufficient for the time comparison.
For any integer `s>=1`, split the prime factors at `2^s`. There are at
most `2^s` small primes, each with exponent at most `n`. For large
primes, `sum e_i<=n/s`, and `e_i+1<=2^{e_i}`. Hence

```math
\log_2\tau(N)\le2^s\log_2(n+1)+n/s.                    \tag{9}
```

For large `n`, choose `s=floor((log_2 n)/2)`. Then the first term is
`O(sqrt(n) log n)=o(n/log n)`, and the second is `O(n/log n)`.
Absorb the finitely many smaller cases into a constant. Thus

```math
\tau(N)=\exp(O(n/\log(n+2))).                           \tag{10}
```

Generating and checking one divisor costs polynomial time in `n` and
linear space. The count register also has `O(n)` bits. Equation (10)
is lower order than the time in (1). Exact prefix counting, interval
counting, and streaming all divisors therefore fit the same half-rate
time and `O(n)` workspace after factorization. No table of all divisors
is retained.

## Verification boundary

The companion checker compares coordinatewise sixteen-element-ring
queries with the full tag algebra and with direct bucket sums. It
checks on-demand phases, singleton isolation, candidate verification,
prime endpoints, and least-factor recovery on finite instances. It
also checks mixed-radix divisor enumeration and exact prefix counts
against independent divisibility tests. The
local factor checker independently verifies the compressed partitions
and the exact factor values used by Section 4.

These tests do not instantiate the asymptotically huge clique evaluator.
The uniform time and linear-space statements rely on the proofs above
and the cited evaluator and local-factor arguments.

```sh
python tools/check_linear_space_recovery.py \
  --output evidence/linear-space-recovery.json
```
