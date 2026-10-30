# A counterexample to whole-root-interval parity recovery

## Exact finite result

The following integer is composite and is not a perfect power:

```math
N=61\cdot71^2\cdot380117^2\cdot175531457527^2
 =1368960876177490538878515064474553061781.             \tag{1}
```

All four displayed factors are prime. Nevertheless, for every integer
`k>=2`,

```math
\#\{d\mid N:1<d\le\lfloor N^{1/k}\rfloor\}
\equiv0\pmod2.                                          \tag{2}
```

Thus a procedure that searches only for an odd **whole** root-interval
count cannot always start parity witness recovery, even after excluding
perfect powers. This does not rule out a more general algorithm that
queries other intervals, individual binary-prefix blocks, or adaptively
chosen prefixes. It is not a lower bound for factoring.

The [power-divisor construction](parity-power-divisor-recovery.md) avoids
this obstruction by counting `d` such that `d^(2^j)|N`, which is different
from counting divisors satisfying only `d^(2^j)<=N`.

## Integer certificate

There are exactly `2*3^3=54` divisors. Enumerate each once as

```math
d=61^{e_0}71^{e_1}380117^{e_2}175531457527^{e_3},
\quad e_0\in\{0,1\},\quad e_1,e_2,e_3\in\{0,1,2\}.
```

For every divisor `d>1`, let `m(d)` be the largest integer such that
`d^m<=N`. The divisors with `m(d)>=2` have this exact histogram:

| `m(d)` | Number of divisors |
| ---: | ---: |
| 2 | 10 |
| 3 | 6 |
| 4 | 2 |
| 5 | 2 |
| 7 | 2 |
| 10 | 2 |
| 21 | 2 |

Every histogram entry is even. The count in (2) is the sum of the
entries with `m(d)>=k`, proving (2). In particular, the counts for
`k=2,...,10` are

```text
26, 16, 10, 8, 6, 6, 4, 4, 4.
```

They equal two for `k=11,...,21` and zero for `k>=22`.
The last assertion also follows directly from `61^22>N`, since 61
is the least nonunit divisor.

The integer in (1) is not a perfect power because the prime 61 has
exponent one. Its divisor count is 54, so the ordinary square-root
parity is even, as the complement-pairing formula predicts.

## Reproduction and provenance

Run

```sh
python3 tools/check_root_threshold_counterexample.py
```

The checker verifies primality of the four stated factors by exhaustive
trial division, multiplies them to obtain (1), enumerates all 54 divisors,
and proves every histogram bin using the two integer inequalities
`d^m<=N<d^(m+1)`. It uses no floating-point logarithms, probabilistic
primality test, or numerical-rank tolerance.

The example was found by a mixed-integer search over normalized logarithmic
weights. That search is not part of the certificate or a dependency of
the checker. The exact four primes and integer inequalities above are
sufficient to reproduce the result with Python's standard library.
The discovery and checker used substantial OpenAI assistant support.
