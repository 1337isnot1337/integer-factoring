# Arbitrary-bound divisor maxima and logspace record membership

## Result and verified implementation

The [larger checkpoint follow-up](large-divisor-checkpoints.md) records independent exact agreement through 32,768 bits and fixes the CLI decimal-conversion limit. The mathematical search and certificate algorithm below are unchanged.

For an integer bound `X>=2`, define

\[
D_X=\max_{1\le m<X}\tau(m),\qquad
M_X=\min\{m<X:\tau(m)=D_X\}.
\]

**Theorem.** One explicit deterministic algorithm computes `D_X`, `M_X`,
the prime exponents of `M_X`, and `ceil(log2 D_X)` in time polynomial
in the binary input length of `X`. A version with read-only input and
streamed output uses `O(log n)` working bits, where
`n=ceil(log2 X)`. The algorithm has no numerical constants imported
from an asymptotic highly-composite estimate.

The practical Python implementation uses ordinary integers and polynomial
workspace. It extends the
[adaptive benefit oracle](adaptive-benefit-divisor-preprocessing.md)
with arbitrary integer bounds, certified pruning during the height
search, and a computed final radius. The polynomial-time and small-space
proofs retain the dependencies stated there: Nicolas's benefit and
uniform gap estimates, and uniform iterated multiplication for the
logarithmic-space version.

**Corollary.** The language of highly composite positive integers,
written in binary, belongs to deterministic logarithmic space. The
convention is a strict divisor record: `N` is highly composite when
`tau(m)<tau(N)` for every positive `m<N`. Under this convention, `1`
is highly composite.

The auxiliary bound has improved; the complete factoring exponent remains
`1/2`. Maximizing divisor counts below a bound does not factor an arbitrary
integer supplied as input.

## 1. Passing from powers of two to an arbitrary bound

Set

\[
n=\lceil\log_2 X\rceil=\operatorname{bitlen}(X-1).
\]

Then `2^(n-1)<X<=2^n`. For `X>=3`, the rational superior-parameter
search in the predecessor applies unchanged. Its lower endpoint has
superior product at least `2^n>=X`; its upper endpoint has product
`2<X`. It therefore finds adjacent mesh values satisfying

\[
N_\epsilon<X\le N_{\epsilon-1/Q},\qquad Q=\Theta(n^2).
\tag{1}
\]

All candidate prime and exponent bounds remain valid because feasible
integers are below `2^n`. In the stopping certificate, use the actual
logarithm of `X`:

\[
B_\epsilon(D_0)=
\epsilon\log X-\log D_0+
\log\tau(N_\epsilon)-\epsilon\log N_\epsilon.
\tag{2}
\]

The same exact interval arithmetic gives a certified upper bound for
(2). Every integer `m<X` with `tau(m)>=D0` has benefit at most that
bound. The proof of correctness is otherwise unchanged.

The runtime proof also transfers. The last strict record `A<X` and the
uniform next-record gap give `log A=Theta(log X)=Theta(n)` for large
inputs. The mesh error remains `O(n/Q)=O(1/n)`, and the terminating
benefit radius remains a fixed negative power of the largest prime of
`A`. Hence there are polynomially many descriptors, with `O(log n)`
bits per compact descriptor. The separate case `X=2` has `M_X=D_X=1`.

The precise source statements are identified in the predecessor:
Nicolas, *Répartition des nombres hautement composés de Ramanujan*,
Can. J. Math. **23** (1971), 116–130,
[primary PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/4E3B60B2BC31651C55D01FFC9ECB4497/S0008414X00053037a.pdf/repartition_des_nombres_hautement_composes_de_ramanujan.pdf).
The per-record gap proved on pages 128–129 is used; a bound only on
the total number of records would not suffice for this runtime argument.

## 2. Nonnegative prime costs allow safe pruning

Let `a_i=v_(p_i)(N_epsilon)` be the maximizing exponent at prime `p_i`.
Its local loss at another exponent `e>=0` is

\[
c_i(e)=\log(a_i+1)-\log(e+1)
       +\epsilon(e-a_i)\log p_i\ge0.
\tag{3}
\]

The total benefit is `sum_i c_i(e_i)`. The implementation computes
an integer lower bound for each cost, scaled by `Q*2^F`, using the
already proved outward logarithmic intervals. If a rounded lower bound
is negative, replacing it by zero is valid because of (3).

For clarity, suppose `ell(v)<=log v<=ell(v)+err(v)`. A valid unscaled
lower bound for (3) is

\[
\ell(a_i+1)-\ell(e+1)-\operatorname{err}(e+1)
+\epsilon(e-a_i)
\begin{cases}
\ell(p_i),&e\ge a_i,\\
\ell(p_i)+\operatorname{err}(p_i),&e<a_i.
\end{cases}
\tag{4}
\]

The positive reference logarithm is rounded downward. The negative
candidate logarithm is rounded upward. The sign of `e-a_i` determines
which endpoint is used for `log p_i`. Thus every omitted bit of precision
weakens pruning rather than excluding a permitted candidate.

Consider a partial height profile `h_1,...,h_j`.

- Primes beyond the selected support `h_1` already have final exponent
  zero.
- A prime which left the active set between heights `h_(r-1)` and `h_r`
  already has final exponent `r-1`.
- The first `h_j` primes remain active and will have exponents at least
  `j`.

The first two groups contribute their exact local-cost lower bounds.
For an active prime, the minimum possible future cost is zero if
`a_i>=j`; if `a_i<j`, it is `c_i(j)`. The latter fact follows because
all marginal gains after the maximizing exponent are nonpositive, so
its loss is nondecreasing for `e>=j`.

Let `F_j` be the sum of lower costs of frozen exponents, and let
`A_j(h_j)` be the sum of these active-prime minimum lower costs. Then

\[
F_j+A_j(h_j)\le\operatorname{ben}_\epsilon(m)
\tag{5}
\]

for every completion of that prefix. A branch with
`F_j+A_j(h_j)>b` can therefore be skipped during a radius-`b` search.
This does not exclude any candidate promised by the terminal certificate
(2), including a smaller witness tied for the best divisor count.

The implementation stores prefix sums of the local costs. When the next
height is `h<=h_j`, exactly the primes with indices in `(h,h_j]` become
frozen at exponent `j`; their prefix-sum difference updates `F`.
The active lower-cost sum is another prefix lookup. Each node performs
only exact integer operations. At the last height, all remaining exponents
are fixed by the global exponent cap.

These additional bounds are optional for the theorem's correctness: the
unpruned adaptive predecessor remains an independent algorithm. They
substantially reduce the finite search in the recorded larger cases.

## 3. A computed last radius avoids a factor-two overshoot

The predecessor doubles the benefit radius after every failed certificate.
This successor retains that schedule until the computed upper bound makes
a smaller final step available.

Keep all radii on the dyadic grid

\[
\Delta=\frac1{2^{16}R},\qquad
R=\text{least power of two strictly greater than }n+1.
\]

The initial radius is `b0=1/R`. If the current certified upper bound
for (2) is `U>b`, set

\[
b_{\rm next}
=\min\{2b,\;\Delta(\lceil U/\Delta\rceil+1)\}.
\tag{6}
\]

The second choice is strictly larger than `U` and has an additional grid
unit for independent verification. If that choice is smaller than doubling,
the next stage must finish: it searches a cover of every possible improver,
and retaining a larger divisor count can only lower the required upper
bound. The fixed-point lower logarithm used for the divisor count is
nondecreasing. This follows from monotonicity of every positive atanh-series
operation and agreement of range reduction at powers of two.

There are still only `O(log n)` stages. The chosen radius never exceeds
the corresponding doubling step, and the extra grid size is `O(1/n)`.
The polynomial descriptor bound and `O(log n)` counter size from the
predecessor are therefore preserved. No convergence assumption or floating
comparison decides whether to stop.

## 4. Strict-record membership in logarithmic space

For every positive integer `N`,

\[
\boxed{\quad
N\text{ is highly composite}\quad\Longleftrightarrow\quad
M_{N+1}=N.\quad}
\tag{7}
\]

If `N` is a strict record, it uniquely attains the maximum divisor count
on `1,...,N`, so it is the least maximizer. Conversely, if some earlier
integer has at least as many divisors as `N`, then either `N` is not a
maximizer or an earlier maximizer wins the least-integer tie rule. In
both cases `M_(N+1)!=N`.

For `N=1`, the comparison is over an empty set of earlier positive
integers, and `M_2=1`; hence (7) accepts it. Zero and negative integers
are outside the positive-integer convention and the provided predicate
returns false for them. A tied divisor count is not a strict record:
for example, `18` is rejected because `12` is an earlier integer with
the same divisor count.

For the space bound, let the binary input length be `n`. The bits of
`N+1` can be generated from the read-only input with an `O(log n)`-bit
index and carry computation. Run the logarithmic-space oracle on that
virtual input. Retain its compact winning descriptor, regenerate each
bit of `M_(N+1)`, and compare it with the corresponding bit of `N`.
There is no need to retain an `n`-bit witness. Uniform iterated
multiplication and the constant number of virtual-input layers preserve
polynomial time and `O(log n)` workspace. This proves membership in the
usual deterministic class `L`.

The uniform arithmetic dependency is Hesse, Allender, and Barrington,
*Uniform Constant-Depth Threshold Circuits for Division and Iterated
Multiplication*, JCSS **65** (2002), 695–716,
[primary paper](https://people.cs.rutgers.edu/~allender/papers/division.pdf),
Corollary 6.3. This is a proved dependency of the small-space version;
the Python implementation uses normal integers and caches.

## 5. Resource accounting for the new features

The input itself is the binary bound `X`. No factorization of `X` is part
of the interface. The superior integer is constructed from its maximizing
prime exponents, and the winning integer is represented by its height
profile.

For the small-space implementation, local-cost prefix sums in Section 2
are recomputed when needed. Their counters and fixed-point values use
`O(log n)` bits. The prefix recursion can be replaced by the same compact
Cartesian descriptor enumeration as before; pruning is an optional test
on that descriptor. The prime-exponent output is streamed by recomputing
`e_i=#{j:h_j>=i}`. Its total output length is polynomial in `n`.

The ordinary Python search stores the tables and partial integer products.
Its space is polynomial in `n`, and it visits a subset of the polynomial
candidate box proved in the predecessor. Its time remains polynomial.
These are asymptotic statements with fixed absolute exponents; the source
estimates do not supply small numerical runtime exponents here.

## 6. Reproduction and independent checks

Use the command-line interface as follows:

```sh
python3 tools/certified_divisor_oracle.py --bits 8192
python3 tools/certified_divisor_oracle.py --below 100000
python3 tools/certified_divisor_oracle.py --is-record 83160
python3 tools/check_certified_divisor_oracle.py --large
```

The arbitrary bound `100000` returns least witness `83160` with `128`
divisors. The evidence in
[`certified-divisor-oracle.json`](../evidence/certified-divisor-oracle.json)
includes exact large outputs and the following independent checks:

- A divisor sieve inspects all `65,535` positive integers below `2^16`.
  It verifies 384 arbitrary-bound queries, including neighbors of every
  strict record and a spread across the full range.
- There are 260 strict-record predicate controls, including `N=1`, zero,
  and tied divisor counts.
- Exact rational logarithmic intervals verify 3,365 superior marginal
  inequalities, the adjacent parameter brackets, and the terminal benefit
  certificates.
- A separate search regenerates terminal windows using rational intervals
  and enumerates in prime-exponent order. It verifies the least witnesses
  and divisor counts through `X=2^8192`. The main search uses height order
  and different fixed-point lower bounds. Their profile counts need not
  agree because their safe pruning strengths differ.
- Five additional arbitrary bounds, including `10^12+39` and `2^128-17`,
  receive full independent benefit-certificate replays.

The independent prime-order checker collapses a suffix of exponents forced
to one into one exact product and divisor-count update. This avoids a deep
recursion through a long list of forced primes; it does not omit any choice.

| Bit bound `n` | Divisor parameter `ceil(log2 D(n))` | Unpruned adaptive profiles | New pruned profiles |
| ---: | ---: | ---: | ---: |
| 128 | 29 | 4,466 | 184 |
| 256 | 47 | 3,321 | 155 |
| 512 | 80 | 3,077 | 132 |
| 1,024 | 137 | 298,964 | 1,355 |
| 2,048 | 239 | 176,133 | 596 |
| 4,096 | 426 | — | 935 |
| 8,192 | 768 | — | 6,762 |

Profiles are summed across the adaptive stages. The JSON separately
records recursion nodes and the terminal independent enumeration. These
finite counts support reproducibility; the polynomial theorem rests on
the proof, not on a fitted timing or growth curve.
