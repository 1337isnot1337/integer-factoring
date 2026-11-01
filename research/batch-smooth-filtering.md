# Batched smooth-part filtering with charged early stopping

## Result and scope

The [fixed-base power selector](smooth-leading-power-selection.md) and
[prime-cofactor saturation](prime-cofactor-saturation.md) now have an opt-in
`smooth_batch_size` parameter. The default remains one. A larger value uses a
product/remainder tree to reduce the common prime product modulo several
candidate norms at once.

This implements the classical algorithm of D. J. Bernstein,
[*How to find smooth parts of integers*](https://cr.yp.to/papers.html#smoothparts),
2004, Algorithm 2.1. The smooth-parts method and its product-tree acceleration
are not new mathematical results. The contribution here is its integration
with the two exact class-group selectors, preservation of their first-success
witnesses, explicit ahead-work and budget accounting, and matched finite
measurements. The source baseline is
`941cc1d8f56e6e80ec45938f1d660532de2f5046`.

The change does not bound the number of rejected powers. It therefore does
not improve the complete factoring score, which remains `1/2`.

## 1. Exact smooth parts without factoring a candidate

Let `P>=1` be the product of the allowed primes. It may contain repeated prime
factors; only its support matters. Let the positive candidate norms be
`A_1,...,A_m`, each with at most `b` bits. Put

$$E=2^{\lceil\log_2 b\rceil}.$$

Then `b<=E<2b`, except that `E=b=1` at the smallest boundary. For each candidate,
compute

$$r_i=P\bmod A_i,\qquad S_i=\gcd(A_i,r_i^E\bmod A_i).$$

Every prime valuation of `A_i` is less than or equal to `b`. If a prime `q`
divides `P`, its valuation in `P^E` is at least `E` and therefore at least its
valuation in `A_i`. If `q` does not divide `P`, its valuation in `P^E` is zero.
Thus `S_i` is exactly the largest divisor of `A_i` supported on the allowed
primes. In particular,

$$A_i\text{ is supported on }P\quad\Longleftrightarrow\quad
r_i^E\equiv0\pmod {A_i}.$$

The convention for `A_i=1` gives remainder zero, smooth part one, and a true
smoothness flag. The empty prime base has `P=1` and accepts only candidate one.
There is no coprimality requirement between different candidates. Repeated
norms, shared prime factors, and arbitrarily high prime powers are allowed.

The power-of-two exponent avoids the extra multiplications in an arbitrary
binary exponent. This changes neither the smooth part nor the accepted set.

### A bounded gcd shortcut

The implementation first sets `S=gcd(A,r)` and `T=A/S`. Every allowed prime
that divides `A` already divides `S`. If `gcd(S,T)=1`, the supported part is
therefore complete and no modular power is needed. Otherwise replace
`S` by `S gcd(S,T)` and divide `T` by the same gcd. The product `ST=A` is
preserved, and all primes in `S` remain allowed.

For a supported prime with input valuation `e`, one such replacement sends
its current valuation `v` in `S` to `min(2v,e)`. The code makes at most two
replacements. If they do not finish the smooth part, every supported prime
has valuation at least `min(4,e)` in `S`. The exact fallback is consequently

$$\gcd\bigl(A,S^{\lceil E/4\rceil}\bmod A\bigr).$$

This still saturates every input valuation because `4 ceil(E/4)>=E>=b`.
The bounded shortcut adds only a constant number of gcds and divisions;
it retains the leaf bound `O(M(b) log(b+1))`. It is the bounded version of
the gcd/powering combinations discussed in Section 3 of Bernstein's paper,
not a smoothness heuristic. The large prime-power controls exercise the
fallback as well as the early exits.

## 2. One remainder tree for a whole batch

Build a balanced binary product tree whose leaves are the candidate norms.
At each internal node store the product of its children. Reduce `P` modulo
the root product. Descend the tree, reducing each parent's remainder modulo
each child product.

At every node with product `Q`, the stored value is `P mod Q`. A child's
product divides its parent's product, so reducing the parent remainder gives
the same result as reducing `P` directly. This proves the invariant by
induction. The leaf remainders are precisely the `r_i` above.

This is ordinary integer division. It does not use CRT inversion and does
not require pairwise coprime norms. The implementation never writes `P^E`.
It performs only the small modular powers after the common reductions.

### Bit costs

Write `R=bits(P)`, `b_i=bits(A_i)`, and `L=sum_i b_i`. Let `M(t)` be a regular
fast multiplication bound, with division of `t`-bit integers charged within
`O(M(t))`. The product and remainder trees cost

$$O\bigl(M(R+L)+M(L)\log(m+1)\bigr).$$

The leaf modular powers and gcds cost

$$O\left(\sum_i M(b_i)\log(b+1)\right).$$

Consequently the total is bounded by

$$O\bigl(M(R+L)+M(L)(\log(m+1)+\log(b+1))\bigr).$$

The direct stored-tree implementation uses `O(L log(m+1))` additional bits
for its products, plus the original `R`-bit prime product. Holding `m`
quadratic-form states separately costs `O(m n)` bits for an `n`-bit input.
All prime-product construction, promotions, form compositions, successful
norm decoding, and relation extraction remain charged by the existing
selectors.

These are fast-arithmetic bounds for the stated algorithms. Python's integer
division implementation is not asserted to attain them. The finite timings
below concern the actual Python implementation.

The benefit is largest when the common product is much longer than each
candidate. The old scalar path repeatedly reads and reduces that same long
integer. A small base can make the tree overhead larger than the work saved;
there is no unconditional wall-time improvement for every input or batch size.

## 3. Preserve first-success selection and charge every generated norm

For the fixed-base selector, consecutive powers are generated in blocks of
at most the requested size. The flags in a block are inspected in the same
order as the scalar routine. The first accepted norm and every later relation
therefore agree with the scalar run when the search is not stopped by a
different global budget.

The promoted selector starts each generator's blocks at sizes
`1,2,4,...`, capped by `smooth_batch_size`. This avoids constructing a full
large block for targets that succeed immediately. The saturation product is
fixed throughout one generator's search. It changes only after an accepted
norm; the remainder of that final batch is discarded before the next
generator. No stored smooth part is reused after a promotion.

Suppose a search first succeeds after `q` inspected candidates, and the cap
is `B`. If the final batch ends at candidate `g`, then

$$q\le g\le q+B-1.$$

The same bound holds if an involution or gcd ends the search. An exhausted
complete window has no ahead work. Across `J` started generator searches,
the extra candidate count is at most `J(B-1)`.

The geometric promoted schedule has a stronger bound independent of `B`.
Before a block of length `t`, at least `t-1` candidates have already been
inspected. This follows from the preceding geometric sum and remains true
after reaching the cap or clipping the final block. If `u` candidates
precede the stopping block, then `q>=u+1` and

$$g=u+t\le2u+1\le2q-1.$$

Thus each nonempty generator search constructs fewer than twice as many
norms as its scalar counterpart, regardless of the requested batch cap.
This statement concerns the geometric saturation routine; the fixed-size
blocks in the fixed-base routine have only the preceding additive bound.

The output fields distinguish this work:

* `norm_queries` counts every constructed candidate, including discarded
  ahead work.
* `scanned_norm_queries` counts candidates inspected in scalar first-success
  order.
* `generated_ahead_norms` is their difference.
* `modular_smoothness_tests` counts the filters applied to all constructed
  norms, including the uninspected tail of the final batch.
* The saturation routine's composition counter includes every composition
  used to construct that tail.

For a fixed prime-base run, each attempt's `queries` field counts generated
norms, while its `scanned_queries` field records the first-success position.
The original fields and values are retained when the option is left at one.

The promoted selector also clips a batch to the remaining `max_queries`
budget. It never exceeds that cap. A capped batched run may stop earlier in
the scalar witness sequence because ahead work consumes part of the budget.
Its accepted witness list remains a prefix of the uncapped scalar list.
No equality of completion behavior is claimed at the same insufficient cap.

### A size-balanced amortized filtering bound

For an arbitrary geometric cap `B`, a search with prime-product length
`R`, norm bit bound `b`, and first stop `q>=1` has

$$\widetilde O\left(R+qb+\frac{qR}{B}\right)$$

filtering cost with quasi-linear arithmetic. There are
`O(log(B+1)+q/B)` blocks, their total leaf bit length is below `2qb`,
and each block costs soft-O of `R` plus its leaf bit length. The extra
tree storage is `O(Bb log(B+1))` bits. This is a time-space tradeoff for
the filtering component; state generation and the other selector work
remain charged separately.

The geometric schedule also gives a useful asymptotic consequence without
any promise on where the first acceptable norm occurs. Consider one nonempty
search with a frozen `R`-bit prime product, candidate norms of at most `b`
bits, and scalar stopping position `q`. Choose the cap

$$B=\max\{1,\lceil R/b\rceil\}.$$

This is an analytic choice for the sequence helper; the selector's practical
option remains an explicit integer cap. No machine-independent timing
threshold is inferred from this choice.

There are `O(log(B+1)+q/B)` blocks before stopping. Indeed, at most
`O(log(B+1))` blocks precede the cap, and the remaining blocks have size `B`
except for a possible clipped final block. The already proved inequality
`g<2q` bounds their total size. Every block has total leaf bit length at most

$$Bb\le R+b.$$

The preceding remainder-tree analysis therefore gives the explicit bound

$$O\left(\left(\log(B+1)+\frac qB\right)
M(R+b)\left(\log(B+1)+\log(b+1)\right)\right).$$

For quasi-linear multiplication and division this is

$$\widetilde O(R+qb)$$

bit operations. To see the linear terms, note that
`(R+b)/B<=2b` both when `R>=b` and when `R<b`. The stored tree needs
`O((R+b) log(B+1))` bits, independently of `q`. If a quadratic-form state
occupies `O(n)` bits, the additional state storage is `O(Bn)`.
For the selectors here one may use the leading-norm height's bit length,
which is `Theta(n)`. The tree and batch states then occupy
$\widetilde O(R+n)$ additional bits.

For `J` searches whose products remain fixed within each search, write
`R_j` for their product lengths and `Q=sum_j q_j` for the total scalar
inspected count. Applying the cap separately to each search gives

$$\widetilde O\left(\sum_{j=1}^{J}R_j+Qb\right)$$

filtering work while generating fewer than `2Q` norms. In contrast, separately
reducing the entire current product for every candidate has the direct
bound $\widetilde O(\sum_j q_jR_j+Qb)$. The improvement concerns
smooth-part filtering only. Initial powers, successor form compositions, construction of
the prime products, accepted-word decoding, and relation extraction must
still be added. In particular this analysis does not improve the existing
bounds for `Q` or `J` and does not improve the complete factoring exponent.

## 4. Exact controls and measurements

Run the default controls with:

```sh
python3 tools/batch_smooth_check.py \
  --output evidence/batch-smooth-controls.json
```

The retained controls check 2,152 simultaneous remainders, 6,020 smooth parts
and flags, six invalid-input cases, and 250 finite sequence/batch cases.
They also check the generated-work bound at all 6,125 possible stopping
positions in those sequence controls.
Direct trial division supplies the independent smooth-part reference. The
moduli include repeats, common factors, one, and large prime powers.

There are 80 matched small selector runs over batch caps `2,8,32,128`.
They compare the complete selected witnesses, generator definitions, relation
matrices, factor, and exact extraction result with the scalar path. For
saturation, the full augmented matrix and generator definitions determine
the fixed rows by exact reverse substitution; the fixed-row digests also
agree. A further
100 capped runs check the hard query budget and witness-prefix property.
The existing independent ideal-identity and integer-elimination verifier is
used on the reference records; exact payload equality transfers that replay
to the corresponding batched runs.

The matched large benchmark command is:

```sh
python3 tools/batch_smooth_check.py --benchmark 144 --kind saturation \
  --batches 32 --output evidence/batch-smooth-saturation-144.json
```

`--kind fixed` selects the fixed-prime-base algorithm. The checker records
all generated and inspected candidates and a SHA-256 digest of the complete
mathematical payload. It compares the payload itself in memory; the digest
is an evidence identifier, not a substitute for the comparison.

The retained complete runs use batch size 32. They compare the combined
remainder-tree and bounded-gcd implementation with the original scalar
modular-power filter. Times include the complete selector run and its
relation extraction.

| Fixture | Scalar wall seconds | Batched wall seconds | Generated / inspected norms | Fixed relation rows |
| --- | ---: | ---: | ---: | ---: |
| [92-bit fixed base](../evidence/batch-smooth-fixed-92.json) | 1.088 | 1.069 | 220,064 / 218,738 | 82 |
| [124-bit fixed base](../evidence/batch-smooth-fixed-124.json) | 23.417 | 20.560 | 3,654,400 / 3,650,183 | 274 |
| [144-bit promoted base](../evidence/batch-smooth-saturation-144.json) | 198.180 | 129.698 | 2,430,743 / 2,246,996 | 513 |

The 144-bit run is about `1.53` times as fast while charging all 183,747
generated-ahead norms. Its CPU times are `197.924` and `129.672` seconds,
so the measured improvement is also present in CPU time. Both traversals
return `4108225714090758818561` and have the same complete mathematical
payload. The 124-bit run returns `2305843009213693951` in both modes and
charges all 4,217 extra norms; its CPU times are `23.400` and `20.560` seconds.

The 92-bit times are close enough that a single run does not establish a
stable improvement at the smaller size. These finite observations support
an opt-in implementation, not a universal timing threshold or automatic
batch-size policy.

### Matched existing cofactor policies

The fully promoted 144-bit comparison above concerns one traversal policy.
The repository also has faster unpromoted and selectively promoted policies
on that input. The [separate matched policy control](../evidence/batch-smooth-cofactor-policies-144.json)
checks these with the same scalar and batch-32 filters:

| 144-bit cofactor policy | Scalar wall seconds | Batched wall seconds | Generated-ahead norms |
| --- | ---: | ---: | ---: |
| Unpromoted product | 53.327 | 53.064 | 149,926 |
| Promotion up to `16Y` | 60.848 | 51.755 | 163,680 |

The unpromoted pair has essentially equal times. The selective pair is
about `1.18` times as fast; its CPU times are `60.735` and `51.745` seconds.
Within each pair, every selected witness, augmented relation, generator
definition, and extraction result agrees. The scalar inspected counts are
3,395,784 and 2,305,811 respectively, and all ahead work is charged.

Thus the fully promoted speedup is not an overall best-policy runtime
record. The selective improvement brings that policy close to the
unpromoted run on this fixture. These measurements further support an
explicit option rather than automatic policy selection.

```sh
python3 tools/batch_smooth_policy_check.py \
  --output evidence/batch-smooth-cofactor-policies-144.json
```

The large retained records report their own run times and source hashes.
Timings depend on the machine and current load. The arithmetic equalities,
first-success comparisons, and counts are the reproducible claims.

## Remaining task

This change reduces the cost of testing supplied consecutive powers. It
does not prove that a suitable power occurs in a short window or improve
the existing bounds on promoted generators. The necessary uniform search
bound remains an obligation in the arithmetic route. The complete factoring
exponent is unchanged.

This implementation and analysis used substantial OpenAI assistant support.
