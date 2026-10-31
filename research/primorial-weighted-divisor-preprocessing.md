# Faster exact divisor-count preprocessing by primorial weights

## Result and scope

For an input-length bound `n >= 1`, let

\[
D(n)=\max_{1\le N<2^n}\tau(N),\qquad
q(n)=\lceil\log_2D(n)\rceil.
\]

These exact parameters and a least integer attaining the maximum can be
computed deterministically in

\[
\boxed{2^{O(\sqrt{n/\log(n+2)})}\operatorname{poly}(n)}
\tag{1}
\]

bit time. A version with read-only input and streamed output uses

\[
\boxed{O(\sqrt{n\log(n+2)})}
\tag{2}
\]

working bits. It stores the current and best exponent descriptors, and
regenerates product bits when comparing them. It does not store `D(n)` or
the maximizing integer in its working memory.

This improves the earlier `2^{O(sqrt(n))}` time and
`O(sqrt(n) log(n+2))` workspace estimates in the
[exact divisor-count parameter](https://github.com/1337isnot1337/integer-factoring/blob/dd28d375ea65863356c67b051f5d4f217562475d/research/exact-divisor-count-parameter.md).
The existing fast ordinary-integer implementation already has the sharper
time bound (1). The smaller-space variant below changes its enumeration.

The number of canonical prime-exponent patterns below `2^n` is itself

\[
2^{\Theta(\sqrt{n/\log(n+2)})}.
\tag{3}
\]

Thus the order of the exponent in (1) is tight for methods that visit
every such pattern. This is an auxiliary preprocessing improvement. It
does not change the proved factoring score `1/2`, and does not lower-bound
factoring or every algorithm for `D(n)`.

No prime-number theorem, partition asymptotic, or new Lean proof is used.

## 1. Exponent differences expose the correct weights

Write the primes as `p_1 < p_2 < ...` and their primorials as

\[
P_j=\prod_{i=1}^j p_i.
\]

Every divisor count attained below `2^n` has a representative

\[
M=\prod_{i=1}^k p_i^{e_i}<2^n,
\qquad e_1\ge e_2\ge\cdots\ge e_k\ge1.
\tag{4}
\]

Replacing the prime bases by the first `k` primes and moving larger
exponents to smaller primes cannot increase the integer or change its
divisor count. This also preserves the search for the least maximizing
integer, as proved in the pinned predecessor.

Set `e_(k+1)=0` and

\[
m_j=e_j-e_{j+1}\ge0.
\]

Then

\[
M=\prod_{j\ge1}P_j^{m_j},\qquad
e_i=\sum_{j\ge i}m_j.
\tag{5}
\]

The map between finite nonnegative vectors `m` and nonincreasing
exponent vectors is a bijection. The actual budget is
`sum_j m_j log_2(P_j) < n`, rather than just `sum_i e_i < n`.

For a cheap integer lower weight, define

\[
\ell_j=\lfloor\log_2(j+1)\rfloor,
\qquad
w_j=\max\{j,\lfloor j\ell_j/4\rfloor\}.
\tag{6}
\]

These weights can be computed from integer bit lengths using
`O(log(j+2))` bits. They satisfy

\[
j\le w_j\le\log_2P_j,
\qquad
w_j\ge\frac{j\log_2(j+1)}{16}.
\tag{7}
\]

For the upper inequality, use `p_i >= i+1`. At least `j/2` factors
of `P_j` are at least `(j+2)/2`; therefore

\[
\log_2P_j\ge\frac j2\log_2\frac{j+2}{2}
\ge\frac j4\log_2(j+1),
\]

where the last inequality follows from `(j+2)^2 >= 4(j+1)`.
Also `P_j >= 2^j`. These prove both parts of the maximum in (6)
are valid lower weights. For the second bound in (7), if `j>=4`
then `j ell_j/4 >= 2`, so flooring loses at most a factor two;
`ell_j >= log_2(j+1)/2` loses at most another factor two. The
three cases `j=1,2,3` satisfy it directly.

Every feasible vector in (5) consequently obeys the integer budget

\[
\sum_jm_j w_j\le n-1.
\tag{8}
\]

Vectors satisfying (8) can include infeasible integers. The algorithm
tests the exact product; it does not accept the lower-weight test as
proof of `M < 2^n`.

## 2. An elementary weighted-partition bound

Here and in this section's estimates, `log` is the natural logarithm.
Let `w_j` be any positive sequence with `w_j >= j` and
`w_j >= c j log(j+1)` for a fixed `c>0`. Define

\[
A(B)=\#\{(m_j):m_j\in\mathbb Z_{\ge0},
                  \quad\sum_jm_jw_j\le B\}.
\]

Then

\[
\log A(B)=O\!\left(\sqrt{B/\log(B+2)}\right).
\tag{9}
\]

All vectors have finite support, since `w_j >= j`. For `beta>0`,
the product of geometric series gives

\[
A(B)\le e^{\beta B}Z(\beta),\qquad
\log Z(\beta)=\sum_{a\ge1}\frac1a
                         \sum_{j\ge1}e^{-a\beta w_j}.
\tag{10}
\]

The series converge by `w_j >= j`, so these expansions and the
interchange of nonnegative sums are valid.

First, for `0<t<=e^(-4)`, put `L=log(1/t)` and
`J=ceil(1/(tL))`. If `j>=J`, then
`log(j+1)>=L-log L>=L/2`. Consequently

\[
\sum_{j\ge1}e^{-tw_j}
\le J-1+\sum_{j\ge J}e^{-ctLj/2}
\le\frac{C_c}{t\log(1/t)}
\tag{11}
\]

for a fixed constant `C_c`. The geometric-series estimate
`1/(1-e^(-u)) <= 1+1/u` supplies the last step. This handles the
small indices explicitly; their contribution is included in `J-1`.

Now take `0<beta<=e^(-8)`, put `L=log(1/beta)` and
`H=floor(beta^(-1/2))`. For `a<=H`, equation (11) applies to
`t=a beta`, with `log(1/t)>=L/2`. Since
`sum_(a>=1) a^(-2)<=2`, this part of (10) is at most
`4C_c/(beta L)`. For `a>H`, the simpler bound `w_j>=j` gives

\[
\sum_{a>H}\frac1a\sum_j e^{-a\beta w_j}
\le\frac1\beta\sum_{a>H}\frac1{a^2}
\le\frac{2}{\sqrt\beta}
=O\!\left(\frac1{\beta\log(1/\beta)}\right).
\]

Thus

\[
\log Z(\beta)=O\!\left(\frac1{\beta\log(1/\beta)}\right).
\tag{12}
\]

For sufficiently large `B`, choose
`beta=(B log B)^(-1/2)`. Both `beta B` and the right side of
(12) are `O(sqrt(B/log B))`, proving (9). Enlarging the fixed
constant covers the remaining finite range and `B=0`.

Applying (9) to (6), using (7), proves the upper bound in (3) and
bounds the entire surrogate search (8), including rejected candidates.

## 3. A matching count of feasible patterns

The elementary prime bound `p_j <= 4j^2` is proved in the pinned
predecessor using the central binomial coefficient. It follows that,
for `j<=k`,

\[
\log_2P_j\le j\log_2(4k^2).
\]

Choose

\[
k=\left\lfloor\sqrt{
       \frac{n-1}{4\log_2(n+1)}}\right\rfloor.
\]

For `k>=1`, take every vector with `m_1,...,m_k` independently
equal to zero or one and all later entries zero. Its logarithmic
product cost is at most

\[
\sum_{j=1}^k\log_2P_j
\le k(k+1)\log_2(2k)
\le 2k^2\log_2(n+1)
\le(n-1)/2<n.
\tag{13}
\]

Here `2k<=n+1` follows from the displayed choice of `k`.
The `2^k` vectors give distinct integers by unique prime
factorization and (5). Therefore there are at least
`2^{Omega(sqrt(n/log(n+2)))}` feasible canonical patterns. Together
with Section 2, this proves (3).

This count restriction concerns complete pattern enumeration only.
An algorithm that computes the maximum without visiting every
pattern is not covered.

## 4. The fast implementation already has the sharper time bound

The pinned
[`exact_divisor_maximum.py`](https://github.com/1337isnot1337/integer-factoring/blob/dd28d375ea65863356c67b051f5d4f217562475d/tools/exact_divisor_maximum.py)
extends nonincreasing prime-exponent prefixes. It checks a prefix
only while its exact product is below `2^n`; the first rejected
extension stops that exponent loop. Each recursive node is a
different feasible canonical integer, including the empty prefix.

Thus its number of recursive nodes is exactly the quantity counted
in (3). It performs at most polynomially many `O(n)`-bit arithmetic
operations per node, generates polynomially many small primes, and
has only polynomial initialization. Equations (3) and (9) prove (1)
for that existing implementation. No timing extrapolation is used.

The implementation stores full products and therefore does not
establish the smaller workspace claim (2). The next section does.

## 5. Enumeration and comparison in smaller workspace

Enumerate the vectors in (8) by descending distinct support indices.
At a node with remaining budget `B` and a previous support index,
loop over smaller indices `j<=B`, and over positive
`m<=floor(B/w_j)`. Append `(j,m)` and recurse with the smaller
budget. Visit the current vector before its extensions, including
the empty root. This enumerates every vector in (8) exactly once.
Every partial vector is itself a vector in (8), so the number of
nodes is `A(n-1)`. The index and multiplicity loops add only a
polynomial factor in `n` per node.

If the support has `r` distinct indices, (7) implies

\[
n>\sum_{j:m_j>0}w_j
\ge\frac1{16}\sum_{i=1}^r i\log_2(i+1)
=\Omega(r^2\log(r+2)).
\]

Consequently `r=O(sqrt(n/log(n+2)))`. Each pair `(j,m)` and each
active enumeration frame uses `O(log(n+2))` bits. The current
descriptor, its stack, and a retained best descriptor therefore
fit (2).

Exact feasibility and objective comparisons must not store large
products. The budget provides a uniform bound even for rejected
candidates:

\[
\sum_i e_i=\sum_j j m_j\le\sum_j w_jm_j\le n-1.
\tag{14}
\]

Thus `M` is a product of at most `n-1` primes, each at most
`4n^2`, with an `O(n log(n+2))`-bit bound. Its divisor count is a
product of at most `n-1` values `e_i+1`, and has at most `n` bits
because `e_i+1<=2^(e_i)`. The callback for any factor regenerates
`e_i` from (5) by scanning the descriptor. A requested prime is
regenerated by trial scanning integers through `4n^2`, with trial
divisors through `2n`, using `O(log(n+2))` bits and polynomial time.

For completeness, a product bit can be regenerated with balanced
recursion over its factors. At a multiplication of child products
`A,B`, stream columns zero through the requested position `k`:

\[
c_j=\operatorname{carry}_j+\sum_{i=0}^j A_iB_{j-i},
\quad \operatorname{bit}_j=c_j\bmod2,
\quad \operatorname{carry}_{j+1}=\lfloor c_j/2\rfloor.
\tag{15}
\]

The common product-bit bound is `O(n log(n+2))`, so the carry and
all counters use `O(log(n+2))` bits per frame. The balanced
recursion has depth `O(log(n+2))`. It uses `O(log^2(n+2))`
additional bits and `2^{O(log^2(n+2))}` bit time per requested
product bit, including leaf regeneration. It stores no partial
product and uses no memoization. This is the primitive checked in
the predecessor; it is not a claim that ordinary Python integers
use sublinear storage.

Check `M < 2^n` by regenerating every bit from position `n` through
the common bound. Compare candidate and best divisor counts by
regenerating their bits from high to low. On a tie, compare `M`
and the best integer in the same way. Store only the winning
descriptor. The output integer and `D(n)` can likewise be streamed
one bit at a time. To output `q(n)`, find the highest set bit of
`D(n)` and check whether any lower bit is set; this distinguishes
an exact power of two without storing it.

The total bit time of this variant is

\[
2^{O(\sqrt{n/\log(n+2)}+\log^2(n+2))}
=2^{O(\sqrt{n/\log(n+2)})},
\]

because `log^2 n = o(sqrt(n/log n))`. Its extra arithmetic
workspace fits inside (2). Finite exceptional lengths are handled
directly. Both algorithms are uniform and deterministic.

## 6. Verification

Run:

```sh
python3 tools/check_primorial_weighted_preprocessing.py \
  --output evidence/primorial-weighted-preprocessing.json
```

The checker independently verifies the surrogate weights against
actual primorials, compares weighted-pattern counts with a generating
function dynamic program, checks the exponent/difference bijection,
and compares exact maxima with both a divisor sieve and an independent
unrestricted exponent search. It also checks the product-bit primitive
with descriptor-regenerated factors and verifies lower-bound pattern
families from Section 3 using integer comparisons.

These finite checks audit the algebra and enumeration. Equations
(9), (12), and (14) give the asymptotic time and workspace bounds.
The practical checker and oracle use full integers; their measured
memory does not implement the smaller-space specification above.

This note and its checker used OpenAI assistant support.
