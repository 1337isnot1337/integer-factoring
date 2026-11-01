# Rough inputs make every local operand directly observable

## Result

A small-factor pretest within the existing construction budget makes the
prime-modulus factor-search constructor entirely elementary. Every local
operand is then recovered by a single Vandermonde interpolation. No
free-digit enumeration, cyclotomic norm argument, lattice reduction, or
closest-vector algorithm is needed for its Boolean table entries.

This is a factor-search-specific refinement. It searches only the proper
divisor lengths that can contain the least prime factor, after explicitly
removing all smaller factors. For each fixed `epsilon>0`, it retains the
[width-sensitive schedule](width-sensitive-prime-constructor.md):

\[
 P=\Theta_\epsilon(\sqrt{n/\log n}),\qquad
 P^2\lambda<(12+\epsilon)n,\qquad B_{\rm tab}=1.           \tag{1}
\]

Each unweighted table entry has polynomial bit cost, and all written
tables plus the pretest cost `exp(O_epsilon(sqrt(n log n)))`. Weighted
entries are uniformly generated polynomial-size products over the
chosen coefficient algebra. Its operation cost is still charged.

The same observation applies to the earlier explicit narrow-prime
schedule. It changes no crossing exponent and makes no new score claim.

## 1. The pretest and the retained lengths

Let `N` have exactly `n>=2` bits. Use the prime constructor's dimensions
`T=A^2`, radix `b=2^w`, and `P` common-index groups. Put

\[
 K=\lceil A/P\rceil,\qquad 2\le K\le A/2,
 \qquad u=A(K-1)w.                                      \tag{2}
\]

The chosen digit width satisfies

\[
 n+1\le w(T-A-1).                                      \tag{3}
\]

Test integers in increasing order from two through
`min(2^u,floor(sqrt(N)))`, by exact division. A discovered first divisor
is prime. If the search reaches `floor(sqrt(N))` without finding one,
return `N` as prime. In particular, this completes the algorithm whenever
`u>=ceil(n/2)`.

Otherwise no prime factor of `N` is at most `2^u`. Retain only exact
first-operand bit lengths

\[
 u+1\le\ell\le\lceil n/2\rceil.                         \tag{4}
\]

Every composite survivor has a least prime factor `p` with `p>2^u`
and `p<=sqrt(N)`, so its bit length is in (4). Conversely, a divisor
in (4) is proper: for the relevant `n>=3`, its bit length is less than
`n`. Thus an empty retained divisor family means that `N` is prime.
There is no need to keep or subtract the divisor `N` in this interface.

## 2. Both operand heights fit the exterior observations

In the instance with exact length `ell`, the input digit capacities are

\[
 d_p=\lceil\ell/w\rceil,\qquad
 d_q=\lceil(n-\ell+1)/w\rceil.                           \tag{5}
\]

The stored digits at any first-transform residue index form an initial
segment, of height at most `ceil(d_p/A)` or `ceil(d_q/A)`. Both heights
are at most `A-K`.

For the second operand, (3)--(4) give

\[
 n-\ell+1\le n-u
 \le w(T-A-1)-1-A(K-1)w
   =w(T-AK-1)-1<w(T-AK).                                \tag{6}
\]

Hence `d_q<=T-AK=A(A-K)`. For the first operand,

\[
 \ell\le\lceil n/2\rceil\le(n+1)/2
 \le w(T-A-1)/2<wT/2\le w(T-AK),                         \tag{7}
\]

where the last inequality uses `K<=A/2`. This proves the same bound
for `d_p`.

A group of size `k<=K` has `A-k>=A-K` external middle recipients.
Therefore every internal first-transform polynomial has at least as
many exterior observations as coefficients. Strip off its nonzero
twist and interpolate from the first required number of observations.
Over the prime modulus `Q`, the selected evaluation points are distinct,
so the Vandermonde determinant is nonzero and invertible. The empty
polynomial case is checked directly. All remaining observations and
actual bounded digit domains are then checked.

This gives at most one input setting for each internal first-transform
gate. The middle gates, inverse gates, and every carry run are forced
from those settings and the boundary labels exactly as in the
[observable-boundary decoder](observable-boundary-constructor.md).
Check the target digits, carry bounds, outgoing equalities, prefixes,
leading digits, and padding. An entry is consequently zero or one, and
is computed in polynomial bit time using only elementary modular
arithmetic. This uniqueness proof uses only (6)--(7) and interpolation.

Weighted entries are zero or the product of the local weights on that
unique setting. All decisions and the circuit topology are independent
of the coefficient values, so the construction remains valid over the
commutative algebras required by the recursive evaluator.

### The observation does not require matching index blocks

More generally, take any partition of the original arithmetic gates in
which each group contains at most `H` middle-recipient indices, with
`2<=H<=A/2`. If both operand heights are at most `A-H`, the same proof
gives polynomial-time Boolean local entries. Every internal first gate
has at least `A-H` observed outputs outside its group, and the directed
transform and carry dependencies then force all other internal settings.
The first, middle, and inverse families need not use matching index sets
for this conclusion; contiguity is also unnecessary. Replacing `K` by
`H` in the pretest (2) supplies the required height bound.

The common contiguous blocks are still used for the stated channel-width
estimate (1). This broader decoding lemma concerns partitions of the
original arithmetic gates. It does not assert that recursively produced
crossing pieces retain the same observation interface.

## 3. The pretest fits the same budget

The exact choice of `K` gives `K-1<A/P`, so

\[
 u<\frac{A^2w}{P}.                                      \tag{8}
\]

The width-sensitive schedule has `A^2w<(1+epsilon/100)n`, and the earlier
explicit schedule has `A^2w<1.001n`. Thus the pretest costs
`2^O(n/P) poly(n)`. It fits the existing written-table scale.

There is also a direct comparison with the channel volume. Since
`lambda>=6K^2f`, `f>2w`, and `PK>=A`,

\[
 P\lambda>12PK^2w\ge12AKw>12u.                          \tag{9}
\]

Every one of the `P 2^((P-1)lambda)` table labels now needs only
polynomial work. All operand lengths in (4) and all prefix queries add
only polynomial factors. This proves (1) and the full written cost.

## 4. Complete proper-factor recovery

Attach the usual monomial weight to the first input and sum the disjoint
length instances (4). For a queried binary prefix, the resulting
polynomial is

\[
 H(X)=\sum_{\substack{d\mid N,\ u+1\le\operatorname{bits}(d)
                    \le\lceil n/2\rceil\\
                    d\text{ has the queried prefix}}}X^d.\tag{10}
\]

Each divisor has one exact length and one global network setting, so its
coefficient is one in every fixed prime characteristic. The supports
are distinct even when different local assignments elsewhere in a
network description would have had the same unweighted count.

If `N` survived the pretest, every prime factor exceeds `2^u`. Writing
`Omega(N)` for its number of prime factors with multiplicity gives

\[
 \Omega(N)<n/u,\qquad
 \tau(N)\le2^{\lfloor(n-1)/u\rfloor}.                    \tag{11}
\]

The [anchored cyclic-polynomial test](sparse-divisor-nonzeroness.md),
Sections 2--3, applies with that smaller support bound. Its family size,
ring dimension, and preparation cost are `exp(O(n/u+log n))`. This is
the same pretest principle recorded independently in
[PR 31](https://github.com/1337isnot1337/integer-factoring/pull/31);
the new use here is to force the local height inequalities (6)--(7).

Moreover, (3) implies `w>n/A^2`, and `A<=PK`, so

\[
 n/u<A/(K-1)\le PK/(K-1)\le2P.                          \tag{12}
\]

Hence the cyclic recovery overhead is only `exp(O(P+log n))` on these
survivors. For the schedules in (1), it is below the written-table cost.

Encode every divisor as an `n`-bit binary word with leading zeroes.
Test the empty prefix first. If it is empty, return `N` as prime by
Section 1. Otherwise descend those binary prefixes, keeping the left child
when its polynomial is nonzero. This returns the least divisor in the
retained family, which is the least prime factor of `N`. Divide it out
and repeat, as in the existing complete-factorization transfer. The
number of prefixes and factors is polynomial in the input length.

This supplies a complete arithmetic interface, while the evaluator's
crossing exponent remains a separate hypothesis.

## 5. Exact replay

The companion checker
[`check_rough_observable_factor_constructor.py`](../tools/check_rough_observable_factor_constructor.py)
checks the retained-length completeness on small integers, verifies
(6)--(7) at exact schedule boundaries, and checks direct exterior
interpolation with no free digit on explicit prime-modulus networks.
Run:

```sh
python tools/check_rough_observable_factor_constructor.py --output evidence/rough-observable-factor-constructor.json
```

The committed [evidence](../evidence/rough-observable-factor-constructor.json)
records 5,493 rough-input support/completeness cases, 2,950 least-prime
recoveries using 32,995 zero-padded prefix decisions, 18,236 small-factor
returns, and 841 complete trial-prime returns. It also records 2,626
exact height-boundary checks, 119 pretest-budget checks, 160 direct
interpolations with no free digit, and 30 local entries checked against
10,828 complete local assignments. All comparisons passed.
