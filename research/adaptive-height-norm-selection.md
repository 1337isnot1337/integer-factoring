# Adaptive target powers in bounded-height norm selection

## Main result

The logarithm-weighted relation method permits the norm selector to
choose its target exponent. For a prime-form generator `g_p`, an
acceptable answer is any positive exponent `k` and primitive smooth
represented norm `n` satisfying

$$n\le H,\qquad n<p^k,\qquad
 1\le k\le K_p:=\min\{j\ge1:p^j>H\}. \tag{1}$$

The previous height-only contract fixed `k=K_p`. Allowing all exponents
in (1) preserves its independence proof, bit-size bounds, and existence
guarantee. It expands the search region without requiring a stronger
existence theorem.

The accompanying finite search also tests leading coefficients before
searching coordinates, and interleaves the remaining target powers.
At the same input, factor base, and height as the earlier two large
height-only fixtures, it finds every required row and recovers both
factors with fewer ellipse coordinate visits. The extra power and
leading-coefficient work is recorded separately. No uniform speed
bound or stronger complete factoring score follows from these finite
comparisons.

## 1. Adaptive powers preserve weighted independence

Let `Delta=-4N`, with `N>1` odd. Fix prime forms `g_1,...,g_m` at
distinct unramified split primes `p_1,...,p_m`, and let
`phi:Z^m -> C_Delta` be the associated class-group map.

Suppose distinct target indices `i` receive answers satisfying (1).
The [primitive norm-to-word converter](smooth-norm-short-words.md)
gives a signed exponent vector `v_i` with

$$\phi(v_i)=g_i^{k_i},\qquad
 \sum_j |v_{i,j}|\log p_j=\log n_i. \tag{2}$$

Set `r_i=k_i e_i-v_i` and `w_j=log p_j`. Then

$$\begin{aligned}
 |r_{i,i}|w_i
 &\ge k_iw_i-|v_{i,i}|w_i\\
 &>\sum_{j\ne i}|v_{i,j}|w_j
  =\sum_{j\ne i}|r_{i,j}|w_j. \tag{3}
\end{aligned}$$

Thus the same logarithm-weighted strict dominance holds. Restricting
to the designated target columns preserves it, so every accepted
partial matrix has independent rows. The
[determinant-free partial extractor](determinant-free-ambiguous-extraction.md)
therefore applies immediately. It can succeed before every target has
received an answer. A failed partial extraction still has only the
torsion-image meaning stated in that theorem.

The upper bound on `k_i` keeps the same exact product guard valid:

$$\begin{aligned}
 \sum_j |r_{i,j}|\log_2 p_j
 &\le k_i\log_2p_i+\log_2 n_i\\
 &\le K_i\log_2p_i+\log_2H\\
 &\le 2\log_2H+\log_2p_i. \tag{4}
\end{aligned}$$

Consequently, the existing
[`weighted_norm_relations.py`](../tools/weighted_norm_relations.py)
wrapper accepts these rows with the original height parameter `H`.
Its pre-exponentiation bit guard and integer-product comparisons need
no modification.

The [uniform norm-existence theorem](powered-prime-smooth-norms.md)
covers every class, including `g_i^K_i`. That old target remains in
the adaptive schedule. Hence the new contract retains the old
existence guarantee under the same pinned analytic input. Smaller
powers need not have acceptable norms; the selector may skip them.

## 2. A complete finite schedule

For a reduced positive form `Q=(A,B,C)`, with `|B|<=A<=C`,

$$Q(x,y)\ge A(x^2-|xy|+y^2)\ge A$$

for every nonzero integer vector. Therefore a power `g_p^k` cannot
contribute an acceptable answer when its reduced leading coefficient
exceeds `min(H,p^k-1)`.

The implemented schedule does the following for each base prime.

1. Compute every reduced power through `K_p`, and discard powers
   failing this elementary minimum test.
2. Test the remaining leading coefficients for a proper gcd with `N`.
   Such a gcd already completes the requested splitting step.
3. In descending exponent order, test every coprime leading coefficient
   for factor-base smoothness. A successful coefficient uses the
   primitive vector `(1,0)`.
4. For every remaining power, begin with the cap
   `min(H,p^k-1,4A)`. Search its complete primitive, coprime ellipse.
   Interleave the powers, multiplying each exhausted cap by four until
   it reaches `min(H,p^k-1)`.
5. Convert the first accepted norm to a relation and independently
   verify its class identity. If all powers exhaust their caps, record
   that target as failed and continue with the other generators.

The [exact norm sieve](complete-smooth-norm-sieve.md) is complete at
each supplied cap. Its multiplicity budget is `bits(H)`, which imposes
no additional restriction on a positive integer at most `H`. Thus an
exhausted target certifies that no primitive, coprime, base-smooth norm
satisfies (1) for that prime. This conclusion concerns the supplied
finite base and height only.

The modular root schedule is cached separately for each target power.
Repeated caps still revisit coordinates; those visits are counted
again. The output also records the power multiplications, leading
division tests, ellipse rows, primitive candidates, modular strikes,
prime divisions, and number of target root schedules constructed.
Base construction, word conversion, relation checks, and extraction
are actual executed steps. No class-group enumeration, determinant,
or supplied annihilator is used.

There are at most `bits(H)+1` candidate powers per prime and
polynomially many caps. This bounds the number of scheduled searches,
but does not bound their coordinate regions by the desired factoring
time. The finite scheduler is a concrete baseline for the remaining
selection problem.

## 3. Exact finite comparison

The comparison uses the same input integers, cutoffs, and heights as
the [fixed-power unsigned-norm fixtures](height-only-selected-integers.md).
The old and new search schedules differ in three ways: adaptive powers,
the leading-coefficient pass, and the initial `4A` cap. The table
compares their combined effect; it does not attribute the improvement
to one change alone.

| Input bits | Prime cutoff | Height | Accepted rows | Fixed-power ellipse visits | Adaptive ellipse visits | Extra leading tests | Power multiplications |
|---|---:|---:|---:|---:|---:|---:|---:|
| 60 | 1,021 | `2^42` | 92 | 53,238 | 31,237 | 139 | 553 |
| 92 | 4,093 | `2^64` | 284 | 8,926,750 | 6,053,115 | 484 | 1,939 |

The 60-bit input is `1000000016000000063`, and the recovered factor
is `1000000007`. Of its 92 rows, 33 use an exponent smaller than
`K_p`; 10 norms are found by the leading-coefficient pass. The
accepted exponents range from 4 to 27.

The 92-bit input is `4951760154835678088235319297`, and the recovered
factor is `2147483647`. Of its 284 rows, 89 use a smaller exponent;
one norm is found by the leading-coefficient pass. The accepted
exponents range from 5 to 16.

The adaptive searches perform 12,706 and 137,877 leading division
tests, respectively. They construct 185 and 779 target root schedules,
containing 17,020 and 221,236 prime-root records. These costs must be
included when assessing elapsed time; an ellipse-visit comparison
alone is not a runtime theorem. Of the accepted rows, 66 and 209 fail
ordinary unweighted diagonal dominance, while all pass the exact
weighted guard and relation checks.

The checker independently reconstructs target powers, coordinates,
norms, words, and all relation identities. It also checks 13 smaller
cases, including early extraction, direct gcd success, an empty base,
and unresolved prime inputs. For failed small targets it uses an
independent square enumeration, justified by
`Q(x,y)>=A(x*x+y*y)/2`, to exhaust all 27 relevant powers. Five
malformed requests are rejected. In total, 393 accepted rows are
independently replayed.

## 4. Consequence for supplied unsigned candidates

The adaptive contract also allows a supplied primitive norm `n` to be
tested against all powers `g_p^k`, `1<=k<=K_p`. The complete orientation
decoder can be called for each candidate target; retain a match only
when `n<p^k`. If the prime `p` is not supplied, all base primes may be
tried as well.

For a polynomial factor base and `log H=O(b)`, this adds only a
polynomial number of target tests. The existing unsigned decoding
overhead `2^O(b/log b)` is unchanged. The selector may therefore
return a stream of integers with the required coverage property,
instead of a preassigned norm for each fixed-power target. Finding
such a stream at the required total cost remains open.

For an imprimitive represented integer, the admissibility test is
instead on its primitive divisor `m=n/c^2`: require `m<p^k` and
`m<=H`. Different representations can have different contents. A
complete implementation must continue through possible matches if
the first representation fails this inequality; rejecting only that
first match would be insufficient. The
[arbitrary-norm decoder](arbitrary-factored-norm-representation.md)
now implements this gate through `maximum_primitive_norm=J`. Set
`J=min(H,p^k-1)`. For each left-half class it retains the least
primitive partial norm, then tests every right-half allocation. If
any represented primitive divisor meets the bound, this search finds
one with the same mixed-radix meet-in-the-middle complexity. The
current adaptive selector produces primitive coordinate witnesses;
the gate makes processing arbitrary supplied integers a separate
complete option.

## Files and replay

- [`adaptive_height_norm_factor.py`](../tools/adaptive_height_norm_factor.py):
  finite selector and weighted extraction.
- [`check_adaptive_height_norm_factor.py`](../tools/check_adaptive_height_norm_factor.py):
  fresh searches, independent failed-target controls, and certificate replay.
- [`adaptive-height-norm-factoring.json`](../evidence/adaptive-height-norm-factoring.json):
  the exact records and source hashes.

```sh
python3 tools/check_adaptive_height_norm_factor.py \
  --output evidence/adaptive-height-norm-factoring.json
```

Use `--small-only` to retain the small controls and the 60-bit comparison
while omitting the longer 92-bit search. The complete supported
factoring score is unchanged.
