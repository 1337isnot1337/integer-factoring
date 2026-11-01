# An integer divisor-count bound with leading coefficient one

## Result

The uniform divisor bound used by the current recovery and direct-count
arguments can be made substantially smaller by optimizing its existing
small/large-prime split. No prime search or factorization is needed.

For an integer `n >= 1`, define

```math
L=\lfloor\log_2 n\rfloor,\qquad
h=\lceil\log_2 n\rceil,
```

and, for `1 <= s <= L`, define the integer

```math
q_s(n)=2^{s-1}h+\left\lfloor\frac{n-1}{s}\right\rfloor.
\tag{1}
```

Let

```math
q_*(n)=\min\left\{
n-1,\ \left\lceil\frac n2\right\rceil+1,
\ q_s(n):1\le s\le L
\right\}.
\tag{2}
```

The last set is empty when `n=1`. Then, for every positive integer
`N < 2^n`,

```math
\tau(N)\le2^{q_*(n)},
\qquad
q_*(n)=(1+o(1))\frac n{\log_2 n}.
\tag{3}
```

The integer `q_*(n)` is computed using `O(log n)` work bits and
polynomial time in `log n`, given `n` in binary. The algorithm scans
only `floor(log2 n)` candidate splits. It does not enumerate integers
up to `2^s` or compute the primes below that bound.

This improves an auxiliary bound, not the factoring score. Replacing
the existing `q=min(n,ceil(9n/floor(log2 n)))` by `q_*(n)` in the
proved hash, cyclic-polynomial, and direct-count constructions changes
their constants while leaving the unconditional time exponent `1/2`.
The new leading coefficient in the universal bound is one instead of
nine. No optimality claim about all possible finite bounds for `tau`
is needed for (3).

The predecessor is Section 1 of
[sublinear-space-half-rate.md at 1476525](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/sublinear-space-half-rate.md).
The applications are the
[direct-count construction](https://github.com/1337isnot1337/integer-factoring/blob/147652574e7ce1ff68b2c378712846cbf438da34/research/direct-sublinear-divisor-counts.md)
and [sparse-divisor nonzeroness](sparse-divisor-nonzeroness.md).
The proof below is self-contained elementary arithmetic.

## 1. Every candidate is a valid upper bound

Write `N=product_p p^a_p`. If `N=1`, its divisor count is one and
all the claimed bounds are immediate. Otherwise

```math
\tau(N)=\prod_p(a_p+1).
```

Fix `s >= 1`. There are at most `2^(s-1)` primes at most `2^s`:
prime two is one of them, and every other such prime is an odd integer
between three and `2^s-1`. Every exponent satisfies `a_p <= n-1`,
since `2^a_p <= p^a_p <= N < 2^n`. Hence each small-prime factor
`a_p+1` is at most `n <= 2^h`. Their total contribution is at
most

```math
2^{2^{s-1}h}.
\tag{4}
```

Let `Omega_large=sum_(p>2^s) a_p`. Because each such prime is
strictly greater than `2^s`,

```math
s\Omega_{\rm large}<\log_2N<n.
```

The left side is an integer, so
`Omega_large <= floor((n-1)/s)`. The elementary inequality
`a+1 <= 2^a` for positive integers gives the large-prime
contribution at most

```math
2^{\lfloor(n-1)/s\rfloor}.
\tag{5}
```

Multiplying (4) and (5) proves `tau(N) <= 2^q_s(n)`.
The prime factorization appears only in this proof; the algorithm
computing (1) does not know it.

The two additional candidates in (2) are also valid. Every proper
divisor of `N>=2` is at most `floor(N/2)`, so

```math
\tau(N)\le\lfloor N/2\rfloor+1\le2^{n-1}.
```

Pairing each divisor with its complementary divisor also gives
`tau(N) <= 2sqrt(N) < 2^(n/2+1)`, which implies the second
candidate in (2). Taking the minimum preserves validity.

## 2. Exact computation

For a positive Python integer, `bit_length(n)-1` equals
`floor(log2 n)` and `bit_length(n-1)` equals `ceil(log2 n)`.
Thus the complete parameter algorithm is:

```python
def divisor_count_bits(n):
    if n < 1:
        raise ValueError("n must be positive")
    floor_log = n.bit_length() - 1
    ceil_log = (n - 1).bit_length()
    bound = min(n - 1, (n + 1) // 2 + 1)
    for split in range(1, floor_log + 1):
        candidate = (1 << (split - 1)) * ceil_log + (n - 1) // split
        bound = min(bound, candidate)
    return bound
```

Every candidate has `O(log n)` bits: `2^s <= n` and
`h <= log2 n+1`. A fixed number of such registers suffices.
There are `O(log n)` iterations, each using integer arithmetic on
`O(log n)`-bit operands. Ordinary schoolbook arithmetic gives, for
example, `O((log(n+2))^3)` bit time and `O(log(n+2))` space.

The function returns zero at `n=1`; this correctly bounds the sole
positive integer `N<2`. For the factoring task, `N>=2` implies
`n>=2`, and the returned bound is positive.

## 3. Leading coefficient one

Suppose `L=floor(log2 n) >= 16`, and put

```math
r=\lceil\log_2(L+1)\rceil,\qquad s=L-3r.
\tag{6}
```

Then `s>=1`. To see this without a real-number estimate, write
`2^k <= L < 2^(k+1)`, where `k>=4`. Here `r=k+1`, and
`3(k+1)<2^k<=L`. The strict inequality holds at `k=4` and
is preserved by induction.

Since `h<=L+1`, `2^L<=n`, and `2^r>=L+1`, the chosen
candidate satisfies

```math
\begin{aligned}
q_*(n)
&\le 2^{L-3r-1}h+\frac n{L-3r}\\
&\le\frac{n}{2(L+1)^2}+\frac n{L-3r}.
\end{aligned}
\tag{7}
```

As `r=O(log L)`, this yields

```math
q_*(n)\le\frac nL
 \left(1+O\!\left(\frac{\log L}{L}\right)\right).
\tag{8}
```

For the matching statement about the numerical parameter `q_*`, each
candidate in (1) is at least `floor((n-1)/L)`. When `L>=2`,
the two other candidates in (2) are at least that integer as well.
Therefore

```math
q_*(n)\ge\left\lfloor\frac{n-1}{L}\right\rfloor.
\tag{9}
```

Equations (8)--(9), and `L<=log2 n<L+1`, prove the asymptotic
equality for `q_*` in (3). This is a statement about the explicit
bound, not an assertion that every integer has that many divisors.

### Comparison with the existing coarse bound

For `n>=4`, choose `s=floor(L/2)`. Then `s>=L/3`, and
the elementary estimate used by the predecessor,

```math
L(L+1)\le6\,2^{L/2}\qquad(L\ge2),
```

gives

```math
2^{s-1}h\le3n/L,\qquad (n-1)/s\le3n/L.
```

Consequently `q_*(n)<=6n/L`, as well as `q_*(n)<=n-1`.
In particular it never exceeds the old integer bound
`min(n,ceil(9n/L))`. The exact parameter scan is substantially
better than this comparison bound for large `n`.

## 4. Consequences for existing transfer arguments

Each application uses only a proved upper bound `tau(N)<=2^q`,
so substitution is immediate after retaining any required `+1` margin:

| Existing operation | Valid parameter after substitution |
| --- | --- |
| Direct exact count over a power-of-two ring | `Z / 2^(q_*+1) Z` |
| Rational direct count with denominator product `D` | `Z / (2^(q_*+1) D+1) Z` |
| Sparse cyclic nonzeroness | support cap `K=2^q_*` |
| Ternary singleton-bucket recovery | least `k` with `3^k>n(2^q_*-1)` |
| Small-modulus CRT endpoint | retain `q_*+O(log n)` bits for its accumulated residue and modulus |

For example, the ternary control parameter has

```math
k\le\left\lceil\frac{q_*(n)+\log_2 n}{\log_2 3}\right\rceil+1
 =\left(\frac1{\log_2 3}+o(1)\right)\frac n{\log_2 n}.
```

The strict comparison in its definition is preserved even when the
upper estimate happens to be integral. The actual algorithm still
chooses `k` using successive integer powers of three, not logarithms.

These changes reduce auxiliary parameters. They do not change the
crossing contraction, number of grid levels, or leading factoring
score. The new bound can be used without changing any existing
formal theorem or asserting that this note has been formalized.

## 5. Reproduction

Run:

```sh
python3 tools/check_optimized_divisor_bound.py
```

The checker compares the bound against independently sieved divisor
counts, verifies every split candidate on that range, tests bit-length
boundaries, checks the exact leading-coefficient inequalities on large
integer inputs, and exercises the downstream modulus and bucket
parameter inequalities. Its evidence file identifies the checked ranges
and checker hash. The proof above supplies the all-input assertion.
