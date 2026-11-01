# Essentially linear construction and verification of compact oriented ideals

## The improved supplied-data bound

Compact norm/root output permits a sharper bound than
`S (n+log(S+2))^O(1)` for the earlier
[oriented construction](oriented-basis-linear-charges.md).
The discriminant need not be processed separately at full length for
every key, and reduced quadratic-form triples need not be emitted.

Let `Delta` be a negative quadratic-order discriminant. Supply a finite
list of odd unramified ideal keys `(R_i,B_i)`, allowing repeated norms
and repeated keys. Roots may be arbitrary signed integers on input.
Let `B` be the total binary input length, counting `Delta` once and
counting every norm and raw root occurrence.

**Theorem.** There is a deterministic algorithm which validates the keys,
constructs pairwise coprime compact atoms `(A_j,C_j)`, and writes every
exact signed ideal word in

$$ (B+2)^{1+o(1)} $$

bit operations. The atoms satisfy `A_j>1`, and their roots are canonical
modulo `2A_j`. For every supplied key, the output means the literal
integral-ideal identity

$$
I(R_i,B_i)=\prod_{w_{ij}>0}I(A_j,C_j)^{w_{ij}}
\prod_{w_{ij}<0}\overline{I(A_j,C_j)}^{-w_{ij}},
\qquad I(R,B)=\left(R,\frac{-B+\sqrt\Delta}{2}\right).
\tag{1}
$$

A separate arithmetic verifier checks these identities in essentially
linear time in the combined input and written certificate size. The
constructor's valid output has `O(B log(B+2))` bits, so its verification
also satisfies the displayed bound. Reading an arbitrarily oversized
candidate certificate is charged to that certificate's own length.

The algorithm uses the classical coprime construction, the
[two-pass orientation theorem](circuit-box-orientation-batch.md), and
batched discriminant arithmetic. It does not find the supplied keys or
prove a stronger relation-supply bound. The complete factoring score
exponent remains `1/2`.

## 1. Compact representation and input normalization

All atoms refer to the one shared discriminant. Writing a norm and a
root specifies the ideal in (1), without writing the trailing coefficient
`(C_j^2-Delta)/(4A_j)` of its quadratic form.

Set

$$n=\mathrm{bits}(|\Delta|),\qquad S=\sum_i\mathrm{bits}(R_i).$$

Normalizing every raw root modulo `2R_i` has essentially linear cost in
the raw root and norm input lengths. This normalization preserves the
congruence modulo `4R_i`, since changing a root by `2R_i` changes its
square by a multiple of `4R_i`. After normalization, all root data have
`O(S)` bits.

Repeated norms can be deduplicated inside an integer construction, but
every original key and every output word still counts. Unit ideals
remain explicit empty words, so an arbitrarily long list of units is
also charged to its input length.

## 2. Validate the large discriminant once per batch

Build a balanced product tree for the moduli `4R_i`. Its total leaf
length is `O(S)`. Reduce `Delta` once modulo the root product, then
propagate remainders down the tree. This computes

$$d_i=\Delta\bmod4R_i$$

in essentially linear time in `n+S`. At every level after the initial
reduction, the total modulus and remainder lengths are `O(S)`.
The moduli need not be coprime, so repeated norms cause no exception.

Each local check now uses only `O(bits(R_i))` bits:

$$B_i^2\equiv d_i\pmod{4R_i},\qquad\gcd(R_i,d_i)=1.\tag{2}$$

The second equality is equivalent to unramifiedness because
`d_i=Delta mod R_i`. Thus the large discriminant contributes once to
the batch reduction, instead of once to every local gcd and congruence.
The product tree, all remainders, all local squares, and all gcds are
included in the cost.

## 3. Two coprime passes have only linear integer data

The first pass constructs pairwise coprime parent norms and sparse words

$$R_i=\prod_j A_j^{e_{ij}}.$$

Bernstein's classical
[coprime-basis and sparse-factorization algorithms](https://cr.yp.to/lineartime/dcba-20040404.pdf),
Theorems 18.2 and 21.3, have essentially linear integer-data cost with
fast multiplication, exact division, and gcd. The
[executable integer engine](classical-oriented-basis-engine.md) retains
those algorithms explicitly.

The incidence inequality

$$\sum_{i,j:e_{ij}>0}\log A_j\le\sum_i\log R_i\tag{3}$$

bounds all required parent/root data. Reduce each key root at its
incident parents by a product/remainder tree. Choose one reference
root at each parent, form each normalized half-difference, and replace
a nonzero half-difference by its gcd with that parent.

As proved in the
[unitary preprocessing note](circuit-box-oriented-batch-controls.md),
these gcds are unitary parent divisors and have total operand length
`O(S)`. The second pass receives the parents and these divisors, also
of total length `O(S)`. Its natural coprime basis contains no unrelated
components, and every parent-to-child exponent is one.

Sparse membership in those second-pass words determines every root
sign. Each output incidence uses at least one bit of its input norm;
there are `O(S)` incidences. Parent-to-child traversal therefore has
essentially linear output cost. No key-by-atom rectangle is built.
The raw-difference option also has linear total input data, but the
preferred unitary option avoids its irrelevant auxiliary components.

At no stage in this compact construction is a reduced quadratic form
computed. Every occurrence of the full discriminant is confined to
the batched validation and constant many input checks.

## 4. A complete arithmetic verifier

The verifier first binds the record to the original input by normalizing
and comparing the supplied keys. It validates the input and atom roots
with the same batched remainder method. Atom roots must be canonical.
For atom norms, set `P=product_j A_j`. Product/remainder trees compute
`P mod A_j^2`, and the test

$$\gcd\left(A_j,\frac{P\bmod A_j^2}{A_j}\right)=1\tag{4}$$

for every `j` proves pairwise coprimality. The quotient is integral and
is congruent to `P/A_j` modulo `A_j`. All atom data are charged.

For each sparse row, the verifier rejects repeated columns, invalid
indices, and zero exponents. Before taking any powers, it checks

$$
\sum_j |w_{ij}|\bigl(\mathrm{bits}(A_j)-1\bigr)
\le\mathrm{bits}(R_i)-1.
\tag{5}
$$

Every valid norm factorization satisfies (5). It prevents a huge
binary-encoded exponent from requesting an enormous intermediate
power. Since `A_j>1`, the total bit lengths of the requested powers
are `O(bits(R_i))` after this guard. Binary exponentiation followed by
a **balanced** product tree then verifies

$$R_i=\prod_j A_j^{|w_{ij}|}.\tag{6}$$

Using a growing sequential product in place of this balanced tree
would not justify the sharp bound. The compact implementation disables
the older norm engine's sequential verification path and performs the
complete checks here instead.

One remainder tree per sparse row verifies all root signs:

$$B_i\equiv\mathrm{sign}(w_{ij})C_j\pmod{2A_j}.
\tag{7}$$

Equation (3), now applied to the output words, bounds the total moduli
and root input lengths. These checks do not repeatedly scan the full
discriminant.

To prove that (2), (4), (6), and (7) imply (1), consider a prime power
`p^v` in an atom norm. The prime is odd and does not divide `Delta`.
Its two quadratic root branches lift uniquely. Equation (7) selects
the same branch as the input ideal, or its conjugate according to the
sign. Raising that atom ideal to `|w_ij|` supplies the whole power
`p^(v |w_ij|)` on this branch. Distinct atoms have disjoint prime
support by (4), and (6) accounts for the entire input norm. Local
ideal factors therefore agree everywhere, proving the literal integral
identity, including higher prime powers.

Prime factors are used only in this proof, not in the algorithm.
Negative exponents in the sparse encoding denote conjugate integral
factors in (1); the verifier does not identify conjugation with an
integral inverse.

## 5. Cost model and excluded output work

The integer input to both coprime passes, the product/remainder trees,
and the sparse reconstruction has total length `O(B)`. A logarithmic
number of tree levels and the classical algorithm's logarithmic
overheads remain essentially linear. Column indices and exponent
encodings contribute `O(B log(B+2))` output bits.

For a deterministic bit bound, sorting and dictionaries can be
implemented with comparison trees or tries; insertion order can be
retained by a linked sequence. No expected-time hashing assumption is
needed. The readable Python reference uses native integers and
dictionaries. Its observed timings are not asserted to realize the
fast-arithmetic asymptotic theorem under adversarial hash behavior.

The output restriction matters. For many small atom norms, their
reduced form triples can each contain a discriminant-sized trailing
coefficient. For example, choose distinct odd primes with product `P`
and `Delta=1-4P*2^h`; the keys `(p,1)` are valid. For sufficiently
large `h`, their forms are already reduced, and their trailing
coefficients have order `h+log P` bits each.
Writing all those form triples can exceed essentially linear total
compact input size. Producing them, reducing arbitrary classes, and
performing later class-group operations are separate tasks.

The theorem does include verification of the literal identities (1).
It uses root and norm arithmetic instead of explicitly multiplying
the ideal HNF matrices. Independent HNF controls below test the
mathematical implication on finite examples.

## 6. Executable controls

The implementation is
[compact_oriented_basis.py](../tools/compact_oriented_basis.py).
Its two-pass body is adapted from the frozen earlier constructor;
the changes are batched validation, compact atoms, and the separate
balanced-product verifier. The original sources remain unchanged.

The [checker](../tools/check_compact_oriented_basis.py) and
[receipt](../evidence/compact-oriented-basis.json) retain 17 small
batches, including translated raw roots, unit ideals, both orientation
variants, and 256 mixed prime-power orientations up to 186 bits.
They compare the complete compact payload with the earlier constructor
and independently check 3,146 literal HNF identities in total.
There are also 572 independent small prime-expansion comparisons.

The complete 64-bit fixture checks all 1,806 keys and 1,854 atoms.
The 80- and 96-bit controls check their first 128 keys, producing 137
and 129 atoms. They do not claim full replay of those larger key sets.
Nine source hashes bind the retained controls.

Five invalid inputs and nine corrupted certificates are rejected.
One corruption supplies an exponent whose binary representation has
100,001 bits; it is rejected by (5) before exponentiation.

Three constructed validation controls each contain 2,048 keys and
discriminants of 4,218, 65,658, and 1,048,698 bits. The batch validator
matches the independent old per-key validator exactly, and the complete
compact output has the expected 24 atoms and signed words. At the
largest size, the recorded single CPU observations are 0.598641 seconds
for old validation and 0.022091 seconds for batch validation. These
are validation observations on constructed inputs, not factoring runs
or a stable runtime benchmark.

The large discriminant contributes once to the initial remainder in
each validation batch. At the largest size, the entire validation tree
records 1,426,122 remainder operand bits, plus separately charged
products and local arithmetic. The 2,048 local unramified gcds use
13,299 operand bits in total.

```sh
python3 tools/check_compact_oriented_basis.py \
  --output /tmp/compact-oriented-basis-replay.json
```

The sharper bound concerns supplied ideals and their exact compact
certificate. Relation generation remains the unresolved cost relevant
to a stronger factoring score.

This implementation and analysis used substantial OpenAI assistant support.
