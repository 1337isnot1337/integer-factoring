# Factoring from one prime field with a chosen trial budget

## Result and scope

The [precision-one Bezout coefficient](prime-bezout-counting.md) can be
used in one chosen odd prime characteristic. Chinese remaindering and
an upper bound on the total number of divisors are unnecessary. A short
trial-division stage instead bounds the number of prime factors with
multiplicity. Ordinary divisor counts become exact on intervals already
known to contain only prime divisors.

Let `n=bit_length(N)` and let `ell` be an odd prime. Suppose the
globally normalized scalar in the linked note can be supplied for every
requested ordinary residue rectangle in bit time `C(n,ell)`, including
its construction and arithmetic. There is a deterministic reduction with

\[
 T(N)\le
  2^{\lceil n/\ell\rceil}\operatorname{poly}(n)
  +\operatorname{poly}(n,\ell)\bigl(1+C(n,\ell)\bigr).
 \tag{1}
\]

Every query uses **the same `ell`, precision one, and one leading
coefficient**. It is enough to choose `ell>n` for polynomial ordinary
overhead and no exponential preprocessing. More generally, a target
saving exponent `alpha` needs only a suitable polylogarithmic
characteristic at least `(log n)^alpha`, provided the provider itself
has the target bound at that characteristic.

The [implementation](../tools/one_prime_bezout_recovery.py) receives that
provider explicitly. It supplies no faster global coefficient evaluator.
The unconditional factoring score is still `1/2`.

## 1. A roughness bound controls prime counts

Set

\[
 b=\lceil n/\ell\rceil,\qquad B=2^b.                 \tag{2}
\]

Trial-divide `N` by all integers below `B`, and also through `ell`.
Stop at its integer square root. The first divisor found is already the
smallest prime factor. If the square-root limit is reached without a
divisor, return `N`, which is prime.

In the remaining case, `ell` does not divide `N`, and every prime
factor is at least `B`. Write `Omega(N)` for the number of prime factors
with multiplicity; this notation is only used in the proof. Since
`N<2^n`,

\[
 B^{\Omega(N)}\le N<2^n,
 \qquad
 \Omega(N)\le\left\lfloor\frac{n-1}{b}\right\rfloor<\ell.
 \tag{3}
\]

Only `O(B+ell)` trial divisions are charged. The algorithm neither
computes `Omega(N)` nor assumes the entire divisor set has fewer than
`ell` elements. That larger set can have zero cardinality modulo `ell`
despite containing many factors.

## 2. Certified prime-only intervals have exact counts

Suppose `L>=2` is a certified lower bound on every prime factor of
`N`. Every divisor `d` of `N` in

\[
 [L,L^2)                                                   \tag{4}
\]

is prime. Otherwise two of its prime factors, with multiplicity, are
both at least `L`, forcing `d>=L^2`.

The number of divisors in (4) is therefore at most `Omega(N)<ell`.
Its residue modulo `ell`, represented between zero and `ell-1`, is
its **exact nonnegative integer value**. The same is true of every
subinterval of (4).

Let

\[
 D_q=\#\{a:a\mid N,\ q\le a,\ N/a\ge2\}\pmod\ell.
 \tag{5}
\]

For `L^2<=N`, the linked prime-adic theorem gives the count in (4)
as `(D_L-D_(L^2)) mod ell`. Each suffix in (5) is a sum of `ell-1`
normalized coefficients, one for each nonzero row residue modulo
`ell`. All their ordinary row and column intervals extend past `N`.
They partition the relevant actual pairs without overlap.

If this interval count is zero, there is no prime factor below `L^2`.
The algorithm may therefore replace `L` by `L^2` while preserving its
certificate. If it is positive, binary search within (4) locates the
least divisor there: every prefix count is exact by the same bound.
There are no prime factors below `L`, so the recovered value is the
smallest prime factor of `N`.

If instead `L^2>N`, the certificate already proves `N` prime. Every
composite integer has a prime divisor at most its square root, which
would be less than `L`. The boundary `L^2=N` is queried; it can contain
the prime factor `L` of a prime square and must not be declared empty.

Starting at `L=B`, the lower bounds are `2^b,2^(2b),2^(4b),...`.
There are `O(log n)` empty intervals before success or primality, and
at most `n` binary-search steps in the successful interval. Caching
shared suffix endpoints gives `O(n)` suffix queries, or `O(ell n)`
individual coefficient queries.

For example, `N=11*13*17*19=46189` and `ell=7` have fourteen proper
divisors. The total proper-divisor count vanishes modulo seven. Here
`b=3`, and `[8,64)` has exactly the four prime divisors; its count is
nonzero and exact. The wrapper therefore recovers eleven even though
a single global proper-divisor zero test would fail.

## 3. A polynomial default field needs no prime-distribution theorem

One may supply `ell` directly. The default implementation finds the
first odd prime above `n` by literal primality tests. This search is
polynomial in `n`; a deliberately loose elementary bound suffices.

Put `X=16(n+1)^4`, let `R_X` be the product of primes at most `X`,
and set `M=X/2`. The central-binomial argument in the linked note gives

\[
 \log_2 R_X\ge \frac{M}{\log_2 X}.
 \tag{6}
\]

For `n>=2`, `log_2 X=4+4log_2(n+1)<=4(n+1)`, so the right side
is at least `2(n+1)^3`. If every prime at most `X` were at most `n`,
their product would divide `n!` and have logarithm at most
`n log_2 n`, a contradiction. Thus a prime `n<ell<=X` exists.
Testing all candidates and their possible divisors has polynomial cost.

At this default field, `b=1`; preliminary trial division goes only
through `ell`. Every output coefficient has `O(log n)` bits. The
query degrees and normalization exponents remain binary integers with
`O(n+log ell)` bits, exactly as in the prime-adic contract.

## 4. The uniform cost target is smaller

The wrapper performs only polynomially many bounded integer operations
apart from the explicit `O(2^b+ell)` trial divisions and its coefficient
calls. The promised uniform provider cost therefore gives (1).

For a fixed `alpha>0`, choose a prime `ell` at least
`(log n)^alpha` and bounded by a fixed power of `log n`. Such a prime
can be found by the elementary selection argument in Section 3 applied
to `ceil((log n)^alpha)`; rational target powers can use ordinary
integer comparisons. Then

\[
 2^{\lceil n/\ell\rceil}
   =\exp\!\left(O\!\left(\frac{n}{(\log n)^\alpha}\right)\right).
 \tag{7}
\]

Consequently a provider with that time bound at the selected
polylogarithmic field transfers the same saving exponent. Its dependence
on `ell` must be included in the hypothesis. An algorithm supplied
separately for each fixed prime, with uncontrolled constants, does not
meet it. A large-degree product or a reciprocal-matching tree also
cannot be enumerated inside an uncharged coefficient call.

This strengthens the precision-one recovery interface: polynomially
many distinct characteristics are not required, and a scored target
does not force the chosen characteristic to have polynomial magnitude
in `n`. It does not reduce the still-open cost of obtaining even one
of the specified globally normalized coefficients.

## 5. Exact controls

The [checker](../tools/check_one_prime_bezout_recovery.py) uses explicitly
supplied finite divisor sets as its oracle. It checks each claimed
prime-only interval against literal divisors, including the zero-total
example above, several empty growth stages, and a prime-square input.
All actual queries in one run use one identical field. Prime selection
and the binary trial-budget contract are also checked independently.

The receipt records 1,541 wrapper runs, including 219 that use the
supplied oracle, 4,226 finite-list coefficient queries, 291 certified
prime-only intervals, and 128 explicit prime selections. The 1,001-bit
control constructs only its trial-budget contract.

The large-input controls construct contracts only. Neither those controls
nor the finite oracle implement a fast coefficient provider. The existing
prime-adic checker separately verifies the normalized scalar against
complete rational Bezout sums.

```sh
python3 tools/check_one_prime_bezout_recovery.py \
  --output /tmp/one-prime-bezout-recovery.json
cmp /tmp/one-prime-bezout-recovery.json evidence/one-prime-bezout-recovery.json
```

This companion uses the prime-adic normalization and residue partition
proved in the linked note. Its additional result is the one-field
prime-interval recovery and the fully charged trial-budget transfer.
