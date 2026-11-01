# Sparse incremental refinement of oriented residual ideals

## Result and scope

The [implementation](../tools/incremental_oriented_coprime_basis.py) maintains
the oriented coprime basis from [the earlier theorem](oriented-coprime-basis.md)
as residual keys arrive. A product tree finds the old atoms that intersect
each new norm. Only those atoms enter the local refinement. Sparse usage
dictionaries replace the dense matrix of all keys against all atoms.

Every supplied ideal retains an exact integral ideal factorization. Thus
the same generator substitution and torsion-image preservation theorem
apply. The improvement concerns the cost of processing supplied witnesses.
It supplies no uniform bound on how many witnesses a factoring search needs.

On the retained 40-bit selector prefix, the new implementation reproduces
all 249 atoms and every signed word for the 250 residual keys. It uses
2719 decision gcds instead of 406191 for one complete old rebuild, or
711903 across the five old checkpoints. It stores 361 nonzero exponents
instead of 62250 dense entries. The independently verified factor 1000003
is therefore preserved on the same 320 queries.

Factoring integers into pairwise coprime factors is classical. Bernstein's
[primary bibliography](https://cr.yp.to/coprimes.html) and
[2005 paper](https://cr.yp.to/lineartime/dcba-20040404.pdf) describe much
stronger batch complexity results. This small incremental implementation
does not claim that essentially-linear cost. It combines a product-tree
index with the existing exact oriented splits and sparse provenance.

## 1. State and local replacement

Fix a negative quadratic-order discriminant Delta. Input keys `(R_i,B_i)`
have odd positive norms, `gcd(R_i,Delta)=1`, and
`B_i^2 = Delta (mod 4 R_i)`. Roots are normalized modulo `2 R_i`.

The current atoms `(d_j,C_j)` have pairwise coprime norms greater than one.
Each atom has a sparse dictionary `input_index -> signed_exponent`.
For each old key, the dictionaries give the literal identity

$$I(R_i,B_i)=\prod_{w_{ij}>0}I(d_j,C_j)^{w_{ij}}
             \prod_{w_{ij}<0}\overline{I(d_j,C_j)}^{-w_{ij}}, \tag{1}$$

where `I(R,B)=(R,(-B+sqrt(Delta))/2)`. In particular,

$$R_i=\prod_j d_j^{|w_{ij}|}. \tag{2}$$

For a new norm R, let H be the old atoms with `gcd(R,d_j)>1`.
Refine only H together with `(R,B)`, using the proved norm split

`(a,u),(b,v) -> (g,u+v),(a/g,u),(b/g,v)`, where `g=gcd(a,b)`,

followed by the root-orientation split from the earlier theorem. Local
exponent vectors are sparse dictionaries too. No primality test or
factorization oracle is called.

Every unchanged old atom is coprime to every local input norm. It remains
coprime to every new local atom. Hence the replacement preserves global
pairwise coprimality.

If the local substitution writes an old atom as signed word V, substitute
V into every old usage of that atom. This preserves (1): a negative old
usage conjugates all local factors. The new input gets its own local word.
Because distinct old atoms have coprime norms, a new atom cannot occur in
the words of two distinct old atoms. Consequently the sparse substitution
does not introduce cancellation between their norm valuations.

An old atom whose norm and canonical root remain unchanged keeps its
usage dictionary. Adding the new key then requires just one new entry.
When an atom splits, only the usages of that atom are rewritten. This
avoids copying the users of unrelated atoms or repeatedly copying the
entire dictionary of a frequently reused unchanged atom.

Induction proves (1) and (2) after every insertion. The independent
verifier also checks them by generic two-dimensional ideal HNF products.
The class-word map therefore commutes with evaluation, so transporting
genuine relation rows retains every old evaluated torsion class, as in
the earlier theorem. This statement concerns evaluated torsion images;
it does not assert that refining the generators preserves the original
integer presentation.

## 2. Product-tree discovery

The tree has a power-of-two number C of leaf slots. A vacancy has value
one. Every internal node stores the product of its descendants.
Slots remain stable when capacity doubles.

At a node with product P, calculate `g=gcd(R,P)`. If `g=1`, no descendant
can overlap the incoming norm. Otherwise recurse into the children using
g as the new gcd argument. Since every child product divides P, this
replacement finds exactly the same overlapping leaves. No leaf is
missed. Unit-valued subtrees require no gcd.

If h leaves overlap and the tree height is `H=log2(C)`, the query uses at
most

$$1+2hH \tag{3}$$

decision gcds. At each depth there are at most h matching nodes, and only
the children of matching nodes are examined at the next depth. When
`h=0`, a nonempty tree stops at its root.

Insertion and removal update one root path. Each change uses at most H
integer multiplications. Capacity doubling rebuilds the product tree;
these multiplications are also counted. Unchanged atom norms need no
tree update.

Equation (3) is a count of gcd operations, not a claim that each gcd has
constant bit cost. The root product can have as many bits as all current
atoms together. The receipt separately records gcd operand bits,
multiplication operand bits, leaf updates, and local refinement work.

## 3. Sparse space and the complete component cost

Let `S=sum_i bit_length(R_i)`, let a be the maximum number of live atoms,
and let n be the bit length of the discriminant. Equation (2) gives

$$\#\{(i,j):w_{ij}\ne0\}\le S. \tag{4}$$

Each atom norm is greater than one, and every usage consumes at least one
bit of an input norm's logarithm. Each exponent has magnitude at most S.
The stored key indices and signed exponents therefore use
`O(S log(S+1))` bits. Atom norms and roots have total bit length `O(S)`.
The reduced atom forms additionally use `O(a n)` bits.

At each tree level the nonunit product bit lengths sum to `O(S+a)`.
The complete tree uses `O(S log(a+1)+1)` bits. Capacity is at most twice
the larger of one and the maximum live atom count. The local sparse refinement has the same
usage bound for its local input norms, whose total bit length is `O(S)`.
Temporary old and replacement dictionaries change these bounds only by
a constant factor. Thus the maintained basis, local refinement workspace,
and sparse output use

$$O(S\log(S+1)+(a+1)n) \quad\hbox{bits}. \tag{5}$$

Requesting the optional dense legacy output adds its explicitly written
`number_of_keys * number_of_atoms` entries. The default output remains
sparse. The bound does not include arbitrary downstream relation matrices
or their presentation certificates.

Time consists of all tree gcds, every product update and rebuild, all
local refinement operations, form reduction, and every rewritten usage
entry. The old elementary local refinement has only a polynomial bound;
it is not replaced by Bernstein's fast batch algorithm here. In
particular, repeated root gcds and updates can process large products.
An operation-count improvement therefore does not imply a global
essentially-linear bit bound. These costs remain charged even when a
later search obtains a factor from only a small part of the stored data.

## 4. Matched finite checks and timing

The [checker](../tools/check_incremental_oriented_coprime_basis.py) verifies
12336 independent product-tree overlap queries, including slot reuse and
five capacity doublings. Each query is compared with a literal scan of
all live leaves, and its visited-node count is checked against (3).

It checks 623 literal ideal identities across 43 matched batches,
including three input orders, 29 mixed-orientation batches, repeated
powers, a retained composite power atom of norm 169, unit norms, empty
input, invalid inputs, and a corrupted orientation certificate. The
legacy dense verifier and the new sparse verifier both replay the finite
ideal controls.

Every checkpoint from the four retained selector prefixes is also
compared with a fresh old rebuild. The final atoms and signed words are
identical, so the previously checked transported rows and factors remain
the same. These comparisons do not assert uniqueness of the atom basis
for arbitrary inputs.

| Input N | Keys | Atoms | Sparse exponent entries | Dense entries | Preserved factor |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1022117 | 7 | 8 | 9 | 56 | 1009 |
| 1040399 | 11 | 9 | 15 | 99 | 1019 |
| 100160063 | 38 | 34 | 43 | 1292 | 10007 |
| 1000036000099 | 250 | 249 | 361 | 62250 | 1000003 |

Seven matched construction-and-export runs on the 250-key fixture
alternate execution order. The retained median wall times are
0.05479 seconds for one old rebuild and 0.005392 seconds for incremental
insertion from empty, approximately a factor of 10.2. Median CPU times
are 0.05479 and 0.005395 seconds. Both outputs agree exactly after
conversion to the same sparse representation.

A separate `tracemalloc` pass records peak Python allocations of 1167924
and 247672 bytes, respectively. These are Python allocations during the
component call, not whole-process resident memory. The timed region
excludes input power generation and independent certificate verification.
It therefore does not establish an end-to-end factoring speedup.

The [receipt](../evidence/incremental-oriented-coprime-basis.json) retains
all seven timings, structural counts, source hashes, and checkpoint
comparisons. Replaying it will produce new machine-dependent timings:

```sh
python3 tools/check_incremental_oriented_coprime_basis.py \
  --output evidence/incremental-oriented-coprime-basis.json
```

This implementation and analysis used substantial OpenAI assistant support.
