# Decode an arbitrary represented smooth norm

The unsigned height-only selector can be relaxed again: its answer need
not have a primitive representation. It may also contain even powers of
odd inert primes. A complete decoder recovers a primitive represented
divisor and the corresponding signed word. Its worst-case cost remains
`2^O(b/log b)` for a norm with `O(b)` bits.

This is a fixed-norm decoding result. It enlarges the acceptable outputs
of a future selector; it does not find the needed norms uniformly fast.

## 1. Remove the primitivity promise

Let `Q` be a primitive positive integral binary quadratic form of negative
discriminant `Delta`. The supplied positive integer `n` is factored and
coprime to `Delta`. Begin with split prime factors

$$n=\prod_{j=1}^r p_j^{e_j}.$$

For each `p_j`, the factor-base data supply a certified prime and a root
`0<=b_j<2p_j`, with `b_j^2=Delta mod 4p_j`. Let `g_j` be the associated
proper prime-form class.

A primitive representation of a divisor `m` of `n` corresponds to a word
in which each prime has only one orientation. An arbitrary representation
`Q(x,y)=n` has content `c=gcd(x,y)` and primitive norm

$$m=Q(x/c,y/c)=n/c^2.$$

Consequently, the possible signed exponents at a split prime are exactly

$$-e_j,-e_j+2,\ldots,e_j-2,e_j. \tag{1}$$

Equivalently, allocate `a_j` copies to one prime ideal and `e_j-a_j`
to its conjugate, for `0<=a_j<=e_j`. The class contribution is
`g_j^(2a_j-e_j)`, and the removed integer content is
`p_j^min(a_j,e_j-a_j)`.

Thus `Q` represents `n`, allowing imprimitive coordinates, if and only if

$$[Q]=\prod_j g_j^{2a_j-e_j}
\quad\hbox{for some }0\le a_j\le e_j. \tag{2}$$

For the reverse direction, the existing
[primitive word-to-norm converter](smooth-norm-short-words.md) constructs
a primitive representation of

$$m=\prod_jp_j^{|2a_j-e_j|}.$$

Multiplying its coordinates by
`c=product_j p_j^min(a_j,e_j-a_j)` gives a representation of `n`, because
`m c^2=n`. This proves completeness using the already established
primitive correspondence; no class order or class-group enumeration is
needed.

The fixed primitive decoder is still useful when its stronger promise is
available. Its two choices per prime can be substantially cheaper than
the `e_j+1` choices in (1).

## 2. Balance the mixed radices

Put

$$T=\prod_{j=1}^r(e_j+1),\qquad M=\max_j(e_j+1).$$

Assign each radix `e_j+1` to the currently smaller product. Starting with
two products equal to one, their ratio stays at most `M`: after a radix
`u<=M` is assigned to the smaller side, the new ratio is at most the
maximum of the old ratio and `u`.

The two final products `T_1,T_2` therefore satisfy

$$T_1T_2=T,\qquad
\max(T_1,T_2)\le\sqrt{MT}. \tag{3}$$

Enumerate all class products on one side and sort their canonical reduced
forms. Enumerate the other side and binary-search for its complement in
the target class. This is a deterministic meet-in-the-middle algorithm.
It uses sorting, rather than an expected-time hash-table assumption.
All allocations at a given prime are produced by a progression with
class-group step `g_j^2`.

Let `B=bits(n)`. Every exponent is at most `B`, so `M<=B+1`.
The total number of choices is the divisor count `tau(n)`. Thus the
decoder costs

$$\sqrt{(B+1)\tau(n)}\,\mathrm{poly}(b+B) \tag{4}$$

in time and at most this much space, where `b` includes the written
form and discriminant sizes. The reconstruction of coordinates uses
two-dimensional ideal HNF and the retained integral reduction map.
For a reduced target, its written size is `O(log |Delta|)`.

### Complete filtering by the primitive norm

Some applications require the recovered primitive divisor `m` to meet
a smaller bound `J`. Testing only the first representation is not
complete: the same supplied `n` can have representations with different
integer contents, and therefore different primitive divisors.

The optional `maximum_primitive_norm=J` parameter handles this exactly.
For each left-table class, retain its allocation of smallest primitive
partial norm. Sorting by `(class, partial_norm, allocation)` makes this
the first record for that class. For every right allocation, multiply
its primitive partial norm by the minimum left norm of the complementary
class. Accept when this product is at most `J`; otherwise continue with
the next right allocation.

If any matched pair satisfies the bound, the minimum left choice for
that same right allocation also satisfies it. This proves completeness
across all contents without increasing either table size or changing
(4). All computed partial-norm products have at most the supplied
norm's bit length; the written size of `J` is also charged.

For an adaptive target power `g_p^k`, setting `J=p^k-1` therefore tests
the strict weighted-dominance gate `m<p^k` correctly. It is insufficient
to reject a candidate just because the first unbounded representation
has `m>=p^k`.

## 3. An elementary bound for the divisor choices

No prime-distribution theorem is required for the useful worst-case
constant. Write `z=log_2 B` and, for sufficiently large `B`, split the
prime factors at

$$U=\frac{B}{z^3}.$$

There are at most `U` small primes, and each contributes a factor at
most `B+1` to `tau(n)`. For the large primes, use `e+1<=2^e` and
`sum e log_2 p<=B`. These bounds give

$$\begin{aligned}
\log_2\tau(n)
&\le U\log_2(B+1)+\frac{B}{\log_2U}\\
&=\frac{B\log_2(B+1)}{(\log_2B)^3}
  +\frac{B}{\log_2B-3\log_2\log_2B}\\
&=(1+o(1))\frac{B}{\log_2B}. \tag{5}
\end{aligned}$$

Combining (3)--(5), including the polynomial imbalance factor, gives

$$2^{(1/2+o(1))B/\log_2B}\,\mathrm{poly}(b). \tag{6}$$

This has the same worst-case leading coefficient as the height-only
primitive decoder. It does not assert equal costs on individual norms:
a large prime power has two primitive orientations but many possible
contents under the relaxed contract.

## 4. Odd inert primes have forced content

Suppose an odd prime `p` is coprime to `Delta` and `Delta` is a quadratic
nonresidue modulo `p`. The quadratic form modulo `p` is anisotropic:
`Q(x,y)=0 mod p` forces `x=y=0 mod p`. Indeed, a nonzero solution would
give a root of the quadratic polynomial and make its discriminant a
square; its leading coefficient cannot vanish when that discriminant
is a nonresidue.

It follows that an inert prime in a represented norm has even valuation.
Its entire contribution is forced content `p^(e/2)`. Strip all such
factors, run (2) on the split part, and restore this content when returning
the full coordinates. An odd inert exponent certifies nonrepresentation.
This only decreases the choice count in (4).

The implementation accepts `root=None` for an odd inert prime and checks
the nonresidue condition. Prime certification remains part of the
supplied factor-base data. It deliberately excludes an inert prime two
from this shortcut. In the factoring application `Delta=-4N`, the relaxed
selector can simply be required to return an odd norm coprime to `N`.

## 5. Consequence for the height-only selector

Choose `K_p` as the least exponent with `p^K_p>H`, as in the
[weighted relation theorem](height-only-smooth-norm-relations.md).
Suppose the selector returns any odd `Y`-smooth integer `n<=H`, coprime
to `N`, represented by `g_p^K_p`. Neither primitivity nor splitting data
need be supplied.

Trial division over the primes up to `Y` factors `n`. The decoder above
recovers a primitive norm `m=n/c^2<=H`, supported on split primes, and
a word `v` satisfying

$$\phi(v)=g_p^{K_p},\qquad
\sum_q|v_q|\log q=\log m\le\log H<K_p\log p.$$

The resulting relation is logarithm-weighted dominant exactly as before.
With a polynomial base and `log_2H=(C+o(1))b`, the processing overhead is

$$2^{(C/2+o(1))b/\log_2b}.$$

Uniform existence still follows from the stronger primitive powered-prime
norm theorem. The new claim concerns the weaker promise needed from a
successful selector. Its search cost remains unresolved.

## 6. Exact implementation and replay

The implementation is
[`arbitrary_norm_representation.py`](../tools/arbitrary_norm_representation.py),
with API

```python
represent_arbitrary_factored_norm(form, prime_factors)

represent_arbitrary_factored_norm(
    form, prime_factors, maximum_primitive_norm=bound,
)
```

It returns the full representation, the primitive representation, their
integer content ratio, and the signed split-prime word. It accepts
unreduced input forms and returns coordinates in the supplied form.

The [checker](../tools/check_arbitrary_norm_representation.py) compares
13,745 decisions against complete integer ellipse enumeration across
21 discriminants, including nonfundamental orders. It finds 1,168
represented norms; 746 of them have no primitive representation at the
supplied norm. It also checks 360 unreduced-form cases, balanced-radix
partitions, and a larger supplied-word norm. The split-square control
`Q=(1,1,6), n=4` requires content two; the inert-square control `n=25`
requires content five.

The bounded mode has 845 additional complete-ellipse controls and competing
content examples. At `Q=(7,-6,16)`, the norm `9*7^16` first gives primitive
divisor `7^14` in the unbounded search, but the bounded search finds
`7^4`; a bound one below `7^4` rejects every allocation. A two-prime
control also requires choosing the least-norm record among repeated
left-table classes. These cases distinguish complete bounded decoding
from filtering a single returned witness.

Run:

```sh
python3 tools/check_arbitrary_norm_representation.py
```

Results and source hashes are in
[`arbitrary-norm-representation.json`](../evidence/arbitrary-norm-representation.json).

## 7. Integer-only factoring consumer

The separate wrapper
[`arbitrary_selected_norm_factoring.py`](../tools/arbitrary_selected_norm_factoring.py)
accepts only `(target_prime, n)` pairs through

```python
factor_from_arbitrary_selected_norms(
    integer, cutoff, height_bound, selected_norms,
    extract_each=True, keep_certificates=False,
)
```

It trial-divides over all primes through the cutoff, classifies the
coprime factors using the constructed prime-form base, and runs the
complete decoder. A successful candidate supplies a primitive norm
divisor to the existing height-only weighted consumer. Re-decoding that
divisor adds a second instance of a cost already bounded by (6), and
preserves the deterministic overall estimate. This wrapper leaves the
earlier primitive-only API and its evidence unchanged.

Even candidates are rejected by this interface. A proper candidate gcd
still returns a factor directly. Failed candidates, retries, duplicate
successful targets, and partial relation extraction are supported.
Written candidate bits and both decoding passes are charged explicitly.

The [consumer checker](../tools/check_arbitrary_selected_norm_factoring.py)
starts with the finite variable-power norms from the height-only replay.
Whenever the height allows it, the checker multiplies a selected norm
by the square of a known odd inert prime. The resulting norm has no
primitive representation, while the target class and height stay fixed.
Only the modified unsigned integers are supplied to the wrapper.

| Input bits | Accepted target rows | Inputs forced to be imprimitive | Recovered factor |
| --- | --- | --- | --- |
| 60 | 92 | 63 | `1000000007` |
| 92 | 284 | 252 | `2147483647` |

Every full and primitive coordinate pair, the content ratio, every signed
word, and every resulting relation are checked independently. A further
control with `N=103`, `H=7^20`, and selected norm `9*7^16` recovers content
21 and primitive norm `7^14`: both forced inert content and optional
split-prime content are removed in the same successful decoding.

Run:

```sh
python3 tools/check_arbitrary_selected_norm_factoring.py
```

The source-pinned finite results are in
[`arbitrary-selected-norm-factoring.json`](../evidence/arbitrary-selected-norm-factoring.json).
The padding construction is a validation device. It is not a new method
for choosing the original target-dependent norms.
