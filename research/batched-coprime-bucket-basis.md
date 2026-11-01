# Prepared bucket remainders for an exact oriented coprime basis

## Result

The supplied oriented ideal keys can be processed with local product
trees and prepared remainders without changing their order, indices,
signed atom words, or ideal identities. On the complete retained 96-bit
prefix, the matched median basis-construction CPU time falls from
4.8112 seconds for the global quotient tree to 1.7320 seconds. The
complete mathematical output remains identical. This is an exact
implementation improvement; it provides no new relation-supply bound
and no improvement to the factoring score exponent `1/2`.

The implementation consists of a [local bucket forest](../tools/coprime_bucket_forest.py)
and a [prepared-batch adapter](../tools/batched_coprime_bucket_basis.py).
Both retain the local oriented refinement algorithm from the
[incremental basis](incremental-oriented-coprime-basis.md).

## 1. A forest of coprime products

The active atom norms are pairwise coprime. Store them in fixed-size
buckets, each with its own binary product tree. The bucket capacity is
a power of two. Empty slots contain 1; a global free-slot heap preserves
the same smallest-available-slot convention as the original tree.
Capacity grows by adding buckets. No root product over all buckets is
maintained during each insertion or deletion.

Within a bucket, changing a leaf from `a` to `b` updates each ancestor
by exact division by `a` and multiplication by `b`. The proof of this
operation is the same as for the
[small-leaf quotient tree](quotient-update-oriented-basis.md).
The complete product is reconstructed by a balanced product tree only
when an exported result requests it; that cost is counted separately.

Suppose a query norm `x` has known common divisor
`g = gcd(x, L R)` with the product of two child subtrees. Coprimality
gives the exact identity

$$g_L=\gcd(g,L),\qquad g_R=g/g_L.$$

Thus one child gcd and one exact division suffice. The implementation
chooses the child with the shorter product. Empty subtrees require no
gcd. The procedure reports exactly the leaves sharing a prime with the
query; it does not require prime factorization of any atom.

After processing one bucket, its common divisor is removed from the
remaining query. Prime powers left in the query cannot occur in a later
bucket, since different bucket products are coprime. This removal
therefore preserves all later overlap decisions.

## 2. Freeze remainders while allowing exact updates

For a batch of already supplied query norms `x_1,...,x_l`, form a
product tree of those norms. For each current bucket product `P_0`,
propagate its remainder down that tree. This produces `P_0 mod x_i`
for every query. The query norms need not be pairwise coprime: this is
ordinary repeated reduction, not a Chinese remainder reconstruction.

The bucket is allowed to change during the batch. Track two integers:

* `D`, the product of the frozen leaves that have since been removed;
* `E`, the product of the current leaves inserted after the snapshot.

The invariant is

$$P_{\rm current}=(P_0/D)E,\qquad
  \gcd(P_0/D,D)=\gcd(P_0/D,E)=1.$$

The first coprimality follows because `D` is a product of whole frozen
leaves. The second follows from the current pairwise-coprime leaf
invariant. Removing a newly inserted leaf divides `E` exactly. Removing
a frozen leaf multiplies `D` and removes its slot from the frozen-slot
map. Consequently, reusing a slot cannot make a later inserted leaf
look like a frozen leaf.

Let `a` be the remaining part of the current query, so `a | x_i`. Its
old common divisor is available from the prepared small remainder:

$$g_0=\gcd(a,P_0\bmod x_i)=\gcd(a,P_0).$$

The current common divisor is then

$$\boxed{\quad
\gcd(a,P_{\rm current})=
\frac{g_0}{\gcd(g_0,D)}\,\gcd(a,E).
\quad}$$

To prove this formula, split the frozen product into the two coprime
parts `D` and `P_0/D`. The displayed quotient removes exactly the
`D` part of `g_0`. The remaining part and `gcd(a,E)` are coprime,
so their product is the gcd with the current product. The formula
also covers refinements that lower prime-power exponents. A newly
allocated bucket has frozen product 1, and hence old common divisor 1.

The resulting common divisor is routed through the current local tree.
There is no additional gcd against the current bucket root.

## 3. Buffering preserves the search prefix

`BufferedBucketOrientedBasis` uses batch and bucket sizes of 512 by
default. Its public key list includes pending keys immediately, and
`add` returns their original sequential indices. A full buffer is
processed in exactly that order. Every `result` call flushes pending
keys before exporting the mathematical result.

The prepared batch contains only keys already supplied by the caller.
It does not generate later powers or inspect later search candidates.
It therefore cannot improve a stopping-query count by using uncharged
future witnesses. A caller may flush at any prefix. Small or frequent
flushes are correct, although they can reduce the measured benefit.

The complete output equality required by the checker covers the
discriminant, input keys, atoms, sparse signed words, norm-splitting
count, orientation-splitting count, and input bit bound. Only tree
strategy labels and operation counters differ. All output ideal
identities are checked independently after each retained fixture.

## 4. Charged costs

The forest avoids multiplying every changed leaf into a large global
root. It instead updates the `log_2 B` ancestors in a bucket of capacity
`B`. If every atom has at most `b` bits, a bucket product, frozen removed
product, or current inserted product has at most `B b` bits.

Preparing a batch builds its product tree once and shares that tree
across buckets. Remainder preparation is explicitly counted. Queries
still scan bucket roots, and deltas still incur gcds, multiplications,
and exact divisions. The implementation does not treat these operations
as free, and the receipt records both counts and operand-bit sums.
The retained counters are measurements of this Python implementation,
not a proof that Python's arithmetic realizes a particular fast-integer
asymptotic bound.

The snapshot also stores current slot identities and the prepared
remainders. For `K` nonempty buckets and a batch of `l` nonunit norms,
there are at most `K l` stored leaf remainders. The largest retained
96-bit batch stores 19,968 such entries. Rebuilding the complete product
on export is charged inside the measured basis-construction interval.

This component refines a supplied family of norms. Classical
[coprime factorization](https://cr.yp.to/lineartime/dcba-20040404.pdf)
already has essentially linear algorithms in the total integer input
length. The present implementation study does not claim a new general
coprime-factorization complexity theorem. Its purpose is an exact,
replayable acceleration of this repository's oriented streaming basis.

## 5. Matched retained measurements

The [receipt](../evidence/batched-coprime-bucket-basis.json) uses three
repetitions of each method and rotates their execution order. The
measured interval includes all insertions and result export. Every
method's complete payload is compared outside that interval, and the
common payload receives independent literal-ideal verification. Other research processes were
active, so CPU measurements are more useful than cross-run wall times.

| Complete supplied prefix | Global quotient tree CPU | Local buckets CPU | Prepared buckets CPU |
| --- | ---: | ---: | ---: |
| 80-bit, 9,776 keys and 10,079 final atoms | 1.273018 s | 0.737970 s | 0.571367 s |
| 96-bit, 19,705 keys and 20,402 final atoms | 4.811191 s | 2.757884 s | 1.732033 s |

The corresponding wall-time medians are 2.851536, 1.676700, and
1.283544 seconds for 80 bits, and 10.888209, 6.210762, and 3.877935
seconds for 96 bits. These are component measurements, not complete
factoring runtimes.

For the 96-bit prefix, total tree-gcd operand bits fall from
11,819,337,591 to 308,434,845. Preparation additionally processes
388,203,658 remainder operand bits, and delta arithmetic remains
separately visible in the receipt. The final global product has
633,029 bits. All three methods make exactly 28,993 sparse usage writes
and return the same 27,584 final usage entries.

## 6. Verification

The focused checker covers 1,476 prepared overlap queries with 2,239
interleaved mutable-leaf transitions. It checks every bucket node, the
current-product invariant, and coprimality of retained and added leaves
after every query. The cases include prime powers,
slot reuse, removal of frozen and new leaves, growth, and batch sizes
that differ from bucket sizes.

Eight oriented ideal banks cover both input orders, repeated keys, unit
norms, and 44 complete prefix comparisons. They verify 1,134 literal
ideal identities. The two full fixtures add 29,481 literal ideal checks,
pairwise-coprime atom checks, and complete mathematical payload equality
against the original global quotient tree.

```sh
python3 tools/check_batched_coprime_bucket_basis.py \
  --output /tmp/batched-coprime-bucket-basis-replay.json
```

The retained receipt records hashes of the checker, both new modules,
and all three shared basis modules. The implementation and analysis
used substantial OpenAI assistant support.
