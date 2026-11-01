# Promote prime residuals without a primality test

## Constructive result

At residual cap `B<=Y^2`, every nontrivial residual accepted by the
[cofactor closure](cofactor-ideal-closure.md) is automatically prime.
Its arithmetic root and the complete initial split-prime base prove
this fact; the program need not factor the residual or test primality.
The residual can therefore enter the saturation product used by all
later norm queries.

The resulting algorithm uses the same FIFO generator order, threshold
power, first-admissible selection rule, and extraction schedule as the
fixed-base closure. Its accepted set grows as more primes are discovered.
The weighted-dominance and integral-elimination proofs remain valid.
Definitions now form a directed acyclic graph, and a sparse reverse
substitution avoids storing an expanded word for every generator.

The implementation is
[prime_cofactor_saturation.py](../tools/prime_cofactor_saturation.py).
Its matched finite controls reduce the number of norm queries on the
92-, 124-, and 144-bit inputs below. The larger saturation product and
larger relation coefficients are charged explicitly. Fewer queries
do not by themselves imply a faster running time. A uniform improved
factoring exponent remains unproved.

## 1. A residual prime certificate

Use the same discriminant `Delta`, height bound `H`, and initial
prime cutoff `Y` as the closure. Let `P_0` be the product of all
initial base primes. The base contains every unramified split prime
at most `Y`; it also contains the applicable prime-two class.

Suppose a residual key has norm `R>1` and root `C`, with

\[
 R\le B\le Y^2,\qquad
 \gcd(R,2\Delta P_0)=1,\qquad
 C^2\equiv\Delta\pmod{4R}.
 \tag{1}
\]

**Lemma.** These conditions imply that `R` is prime.

For any prime divisor `q` of `R`, coprimality makes `q` odd and
unramified. The root congruence makes `Delta` a nonzero square
modulo `q`, so `q` is a split prime. If `q<=Y`, it occurs in `P_0`,
contrary to (1). Therefore every prime factor of `R` is strictly
larger than `Y`. If `R` were composite, including a repeated prime
factor, it would be greater than `Y^2`. This contradicts (1).

In particular, `R>Y`. This proof refers to a hypothetical prime
divisor only logically; it is not a factorization subroutine.
The complete initial prime enumeration, gcd, residual cap, and
root congruence make (1) an independently checkable certificate.

The cap matters. At `Y=5`, for example, the residual norm `49`
from the earlier composite-ideal controls is allowed by a much
larger cap but exceeds `Y^2`. This implementation rejects such a
cap instead of treating a composite ideal as a prime ideal.

## 2. Increasing saturation and exact partial identities

Maintain a product `P` initially equal to `P_0`. At every queried
reduced leading norm `A`, compute

\[
 S=\gcd\bigl(A,P^{\mathrm{bitlength}(A)}\bigr),
 \qquad R=A/S.
 \tag{2}
\]

Modular exponentiation is used for (2). Every prime in `P` is known
either from the original sieve or from the certificate (1). The
exponent in (2) removes its full multiplicity from `A`.

If `R>B`, reject the query. If `R=1`, the leading norm is completely
smooth over the current known-prime set. If `1<R<=B`, (1) proves
that `R` is a new prime. Saturation has already removed every
previously promoted prime, so it cannot be one of those primes.

The middle coefficient supplies the root of the new prime form and
its orientation. The existing partial-word routine gives the exact
class identity

\[
 g_i^k=\prod_{p\text{ currently known}}g_p^{v_p}
              g_R^\epsilon,
 \qquad\epsilon\in\{-1,+1\},
 \qquad\prod_p p^{|v_p|}=S.
 \tag{3}
\]

Omit the residual term when `R=1`. On a new residual, append its
generator and multiply `P` by `R`. The new prime is then available
in the smooth part of every later query. Neither a primality test
nor trial division by an unknown candidate is introduced.

For a fixed comparison run, set `promote_residuals=False`. The
known records still contain all queued residual prime forms, but
the saturation product remains `P_0`. Repeated residuals then
produce the same cancellation relations as the original fixed-base
closure. This mode is checked against its complete published
augmented and eliminated relation matrices.

### Selective promotion limits the product size

The optional `promotion_cutoff=Z` promotes only discovered primes
`R<=Z`. Larger accepted residuals retain the original key-cancellation
behavior. The default is `Z=B`, which promotes every accepted prime.
The same proof (1) certifies all residuals, whether promoted or retained.

This option can keep the saturation product much smaller. A retained
large residual can recur and give a collision relation; a promoted
prime is instead extracted into the smooth word. Both cases satisfy
(3)--(6), and the definition DAG still has smaller-index dependencies.
No search-order change is required. The fixed-base predicate remains
a subset for every choice of `Z`.

## 3. Dominance on the whole augmented matrix

Every generator has the fixed norm weight `w_p=log p`. Start its
search at `K_p=min{k:p^k>H}` and scan consecutive powers. The raw
relation from (3) is

\[
 r_i=k e_i-\sum_pv_p e_p-\epsilon e_R.
 \tag{4}
\]

The right-side word has weighted absolute size at most

\[
 \sum_p|v_p|\log p+\log R
 =\log A\le\log H<k\log p_i.
 \tag{5}
\]

The new known primes are included in this sum. Collecting repeated
columns can only decrease it. Consequently the diagonal of (4)
is positive, and

\[
 (r_i)_i w_i>\sum_{j\ne i}|(r_i)_j|w_j.
 \tag{6}
\]

This includes a target prime occurring in its own smooth word,
negative smooth exponents, and either residual orientation.
Removing unprocessed columns preserves (6), so every selected
prefix is independent. A closed square relation matrix has full
rational rank without a rational-rank search.

## 4. Unit definitions form a DAG

Let the initial base size be `m`. Each new prime first appears in
a row with coefficient `-epsilon`, an integer unit. If `W_i`
denotes its parent's word on the original `m` generators, the new
definition is

\[
 W_R=\epsilon\left(kW_i-\sum_pv_pW_p\right).
 \tag{7}
\]

Every prime on the right was known before `R` was introduced.
Thus introduction order makes (7) acyclic. A definition may depend
on several earlier residual generators, unlike the fixed-base
closure's one-parent forest.

Fresh definitions supply elimination rows. Every selected row
without a fresh generator produces one row on the original base.
After `s` selected rows and `C` definitions, exactly `s-C` original
rows remain. Integer unit elimination preserves independence, so
there are at most `m` such rows. At full closure there are exactly
`m` independent rows, and the existing relation extractor still
has dimension `m`.

The implementation stores each definition in augmented-generator
coordinates. To eliminate an emitted row, process its largest
remaining noninitial index and replace that coordinate using (7).
All dependencies have smaller indices. A priority queue visits
each required definition at most once for that row; stale queue
entries are discarded. The independent checker instead scans the
indices in descending order and performs the same integer
presentation substitutions without a priority queue.

This avoids a table of all `C` expanded length-`m` generator words.
Large records store the compact definitions and raw partial
identities. Replay reconstructs every emitted integer row. A
SHA-256 digest records the reconstructed matrix compactly; the
mathematical correctness comes from the identities and integer
substitutions, rather than from a cryptographic assumption.

### Coefficient growth and arithmetic cost

Put `V=ceil(log_2 H)` and let `K` be the largest selected exponent.
Every smooth correction has `sum_p |v_p|<=V`. Therefore, along a
definition DAG of depth `D`,

\[
 \|W_R\|_\infty
 \le(K+V)\max_{p\text{ in its definition}}\|W_p\|_\infty,
 \qquad
 \|W_R\|_\infty\le\max(2,K+V)^D.
 \tag{8}
\]

Since `D<=C` and `K<=O(n)+W`, the emitted rows and substitution
coefficients have

\[
 O\bigl((C+1)\log(n+W+1)\bigr)
 \tag{9}
\]

bits. The same bound follows by summing path weights during reverse
substitution: each definition's coefficient one-norm is at most
`K+V`, and indices strictly decrease.

There are at most `O(Cn)` definition edges. Across at most `m`
emitted rows, reverse substitution uses
`O(m C n log(C+2))` integer operations, on the charged bit lengths
(9). This is a polynomial bound, not a claim of constant-size
expanded words.

The saturation product can also grow substantially. Its bit length
is at most

\[
 O(m\log(Y+1)+C\log(B+1)).
 \tag{10}
\]

Reading and reducing this integer in each norm query is charged.
The implementation records its maximum size. It is not replaced
by a unit-cost abstract product oracle.

## 5. Uniform interface and its remaining gap

There is at most one residual generator for each prime in `(Y,B]`:
the two orientations give the same canonical prime key. Hence
`C<=pi(B)<=B`. If every reached generator has a successful power
in its length-`W` window, the algorithm uses at most `(m+C)W`
norm queries and closes after at most `m+C` selections.

For every fixed target form, increasing `P` can only decrease the
residual `R` in (2). The original fixed-base predicate is therefore
a subset of the promoted predicate. A uniform fixed-base
admissibility theorem imposed on every allowed prime generator
remains sufficient for this dynamic construction. In particular,
the [subgroup lower-tail moment criterion](cofactor-orbit-moment-criterion.md)
can be used unchanged; a new independence assumption for a moving
predicate is not needed.

The stronger all-allowed-generators scope matters. Two algorithms
with different promoted predicates may reach different generators,
so a guarantee only for one particular fixed-base traversal's
reached set does not automatically transfer to the other traversal.

The total bit cost is polynomial in `n,Y,B,W`, including (9)--(10)
and the existing extractor. Thus subexponential bounds on these
parameters would transfer as before. This note proves the
prime-promotion mechanism and the exact cost accounting, not the
still-missing uniform window bound.

## 6. Matched finite comparisons and replay

The compared runs use the same inputs, `Y`, cap `B=Y^2`, generator
queue rule, threshold powers, first-admissible selection, and
extraction schedule. The only policy change is whether discovered
residual primes enter later saturation products. Disabling that
change reproduces the original closure's selected rows and query
counts exactly.

| Input bits | `Y` | Fixed queries | Promoted queries | Fixed generators | Promoted generators |
|---:|---:|---:|---:|---:|---:|
| 92 | 1021 | 81430 | 54584 | 1062 | 1373 |
| 124 | 4093 | 1391826 | 807221 | 5321 | 6386 |
| 144 | 8191 | 3395784 | 2246996 | 10245 | 13265 |

Both versions find proper factors of the three same inputs. The
promotion changes later selected powers and therefore the subsequent
generator list; it need not reduce that list's size. Nor does it
guarantee an earlier partial relation factor. These effects explain
why the query improvements are not a monotone per-prefix theorem.
For an explicit counterexample to monotonicity of the total count,
the small input `N=163*167=27221` at `Y=5` takes eight queries in the
fixed mode and nine in the promoted mode. Both return a proper factor.
This control is included in the evidence.

The [small controls](../evidence/prime-cofactor-saturation.json)
compare six arithmetic inputs against the original implementation,
verify each literal partial ideal identity and strict weighted
inequality, check exact rational ranks, and reconstruct the integer
DAG substitutions. Budget failures and an oversized cap are also
checked.

The separate [92-bit](../evidence/prime-cofactor-saturation-92.json),
[124-bit](../evidence/prime-cofactor-saturation-124.json), and
[144-bit](../evidence/prime-cofactor-saturation-144.json) records retain
all dynamic witnesses. Each matched fixed run is compared with the
published original augmented matrix and the reconstructed original
base matrix. Large replay checks every selected power, residual
prime certificate, partial identity, dominance inequality, and
integer elimination, without a dense augmented rank computation.

The large records also report generator depth, coefficient bits,
saturation-product bits, and observed wall times. Their timings are
finite measurements on the available machine, not asymptotic claims.

For example, the fully promoted 144-bit run uses a 278420-bit
saturation product. In the matched production run it took about
207 seconds versus about 51 seconds for the fixed product, despite
using fewer norm queries. The larger product and relation work must
both be charged; the norm-query count alone understates that work.
Selective promotion is provided to examine this tradeoff.

The [selective 144-bit control](../evidence/prime-cofactor-saturation-selective-144.json)
uses `Z=16Y=131056`. It takes 2305811 queries, introduces 11602 total
generators, and finds the same factor at original-base rank 512.
Its saturation product has 28932 bits. The matched run took about
55 seconds versus 49 seconds for its fixed-mode control. Thus this
choice retains most of the query reduction while substantially
reducing the overhead of promoting every discovered prime. It still
does not establish a wall-time speedup over the matched fixed mode.

```sh
python3 tools/check_prime_cofactor_saturation.py \
  --output evidence/prime-cofactor-saturation.json
python3 tools/check_prime_cofactor_saturation.py --benchmark 92 \
  --output evidence/prime-cofactor-saturation-92.json
python3 tools/check_prime_cofactor_saturation.py --replay \
  evidence/prime-cofactor-saturation-92.json
```

The other benchmark choices are `124` and `144`.
Add `--promotion-cutoff 131056` to the 144-bit benchmark to replay
the selective control. The benchmark record can also be passed to
`--replay` for its exact witness verification alone.
