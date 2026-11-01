# Height-only smooth norms suffice for a score-one interface

## Main improvement

The arithmetic selector can be asked for **any primitive smooth norm
below a stated height**, without a bound on its number of prime factors
counted with multiplicity. Different prime-form generators use different
target exponents. Logarithmic column weights then establish the required
independence of the relation rows.

For an unsigned norm with `O(b)` bits, the existing orientation decoder
already costs `2^O(b/log b)`, regardless of multiplicity: its exponential
factor depends on the number of *distinct* prime factors. Therefore this
relaxed unsigned-norm contract suffices for any target score at most one.
If the selector also supplies coordinates, processing has
polynomial overhead for a polynomial factor base and preserves any fixed
target score. Both versions enlarge the set of acceptable selector
outputs. The fast selector itself remains open.

## 1. Choose target powers from the height

Let `N>1` be odd, let `Delta=-4N`, and let `g_1,...,g_m` be the chosen
prime forms at distinct unramified split primes `p_1,...,p_m`.
Fix an integer norm-height bound `H>=1`. For each generator choose

$$K_i=\min\{k\ge1:p_i^k>H\}. \tag{1}$$

No logarithms are needed to compute (1): multiply successive powers
until they exceed `H`.

An acceptable answer for target `i` is a positive integer `n_i<=H`
supported on the factor base, coprime to `Delta`, and primitively
represented by the form for `g_i^K_i`. There is no restriction on
`Omega(n_i)`.

The primitive norm-to-word theorem, or the complete orientation decoder,
produces a signed prime word `v_i` such that

$$\phi(v_i)=g_i^{K_i},\qquad
 \sum_j |v_{i,j}|\log p_j=\log n_i\le\log H. \tag{2}$$

The equality in (2) uses the primitive representation and coprimality:
for each prime power in `n_i`, exactly one of the two prime-ideal
orientations occurs. The existing
[word converter](smooth-norm-short-words.md) and
[factored-norm decoder](factored-norm-orientation-search.md) provide
this exact signed word.

Set

$$r_i=K_i e_i-v_i.$$

These are genuine class-group relations. More importantly, with
`w_j=log p_j>0`, equations (1)--(2) give

$$\begin{aligned}
 |r_{i,i}|w_i
 &\ge K_iw_i-|v_{i,i}|w_i\\
 &>\sum_{j\ne i}|v_{i,j}|w_j
  =\sum_{j\ne i}|r_{i,j}|w_j. \tag{3}
\end{aligned}$$

Thus the relation matrix is strictly diagonally dominant after
multiplying column `j` by `w_j`.

The same argument applies to any subset of distinct target indices.
Restrict the corresponding rows to their designated target columns;
discarding the other columns only decreases the right side of (3).
The resulting square minor is nonsingular. To see this directly,
a nonzero kernel vector would have a coordinate of maximal absolute
value, and the corresponding row would contradict strict diagonal
dominance. Multiplication by the positive column weights is invertible,
so the original rational minor is also nonsingular.

Consequently, the
[partial relation-halving theorem](determinant-free-ambiguous-extraction.md)
applies without change. A failed attempt has only its stated
torsion-image meaning. If the generators span the full class group
and all targets receive acceptable norms, the square relation matrix
gives a proper factor whenever `N` has at least two distinct prime
factors. The rows need not generate the full relation kernel.

## 2. Exact weighted validation has polynomial cost

The new wrapper validates (3) on the designated minor by the equivalent
integer comparison

$$p_i^{|r_{i,i}|}>
 \prod_{j\ne i}p_j^{|r_{i,j}|}. \tag{4}$$

There is no numerical approximation to `log p_j`.

For arbitrary binary-encoded exponents, writing the integers in (4)
could require exponential time. The wrapper therefore checks a size
bound before exponentiation. Let `B=bits(H)` and
`W=max_j bits(p_j)`. For a row constructed above,

$$\sum_j |r_{i,j}|\log_2 p_j
 \le K_i\log_2p_i+\log_2 n_i
 \le2\log_2H+\log_2p_i.$$

Because `bits(p)<=2 log_2 p` for every integer `p>=2`, this implies

$$\sum_j |r_{i,j}|\,\mathrm{bits}(p_j)\le4B+2W. \tag{5}$$

The implementation rejects rows violating (5) before forming any
power product. All products in (4) then have `O(B+W)` written bits.
The target exponents also satisfy `K_i<=B+1`, and can be computed
with integers of at most `B+W+1` bits. Thus all extra validation is
polynomial in the written input dimensions and height.

The API in
[`weighted_norm_relations.py`](../tools/weighted_norm_relations.py) is:

```python
least_power_above(base, height_bound)

factor_from_weighted_relations(
    integer, generators, relation_matrix, weight_bases,
    pivot_columns=None, *, height_bound,
    verify_relations=True, keep_updates=False,
)
```

For the prime-form application, `weight_bases[j]=p_j`. Positive integer
bases suffice for the independence check itself; primality is supplied
by the separate factor-base constructor. As before, relation identities
can be omitted from this wrapper's verification only when the caller
has already established them.

## 3. An unsigned norm still has a score-one decoding cost

Suppose the selector returns only `n`, without a factorization,
orientations, or representing coordinates. Trial division over the
constructed factor base recovers its factorization or rejects it.
The complete orientation decoder then uses at most

$$2^{\lceil\omega(n)/2\rceil}\,\mathrm{poly}(b+B+\log Y) \tag{6}$$

operations, where `Y` is the prime cutoff and `omega` counts distinct
prime factors. Repeated powers are kept together in each orientation
choice. They do not increase the exponential term in (6).

If `n` has `k` distinct prime factors, their product is at least `k!`,
so `k!<=n<=H`. The elementary bound

$$\log(k!)\ge k\log k-k+1$$

therefore gives

$$\omega(n)=O\!\left(\frac{B}{\log(B+2)}\right). \tag{7}$$

In particular, if `B=O(b)`, equation (6) is `2^O(b/log b)` even when
`Omega(n)` is much larger. For a sharper constant, if
`log_2 H<=(C+o(1))b` with fixed `C>0`, then

$$\omega(n)\le(C+o(1))\frac b{\log_2 b},$$

and the orientation factor is at most

$$2^{(C/2+o(1))b/\log_2 b}. \tag{8}$$

This relaxation trades a larger possible decoding constant for a less
restricted selector. The multiplicity-bounded interface remains useful
when one wants `2^O(sqrt(b))` processing or its smaller explicit
constants. Equation (8) suffices for the score-one goal and every
smaller score. It does not by itself supply the same guarantee for
scores larger than one.

## 4. Uniform existence survives the changed powers

Under the pinned fixed-strip analytic input, the
[powered-prime norm theorem](powered-prime-smooth-norms.md) supplies
smooth primitive norm witnesses for **every** class. Hence it applies
to the altered targets `g_i^K_i` in (1).

For example, its polynomial-base specialization gives

$$Y=b^{64+o(1)},\qquad
 H=2^{(3/2+o(1))b}.$$

The construction provides acceptable witnesses below `H`, and the
prime forms generate the full class group for sufficiently large `b`.
The selector need not meet the original multiplicity guarantee; any
primitive factor-base-smooth represented norm below `H` is acceptable.

There are polynomially many targets, all target powers have polynomial
cost, and weighted extraction has polynomial cost. By (8), each supplied
norm can be decoded in

$$2^{(3/4+o(1))b/\log_2 b}\,\mathrm{poly}(b).$$

Thus a deterministic selector with total cost `2^O(b/log b)`, allowing
that many candidate integers and charging their written lengths, would
complete a score-one factoring procedure through this height-only
interface. Polynomial primality and perfect-power preprocessing and
recursive splitting are handled as in the
[arithmetic framework](arithmetic-framework.md).

Uniform existence and verification of supplied answers are established.
The required selector is still missing, so the complete supported
factoring score remains unchanged.

### Polynomial-base height tradeoff

More generally, fix an integer `d>=2` and set

$$q=\lceil\log_2 b\rceil,\qquad v=32dq,\qquad Y=2^v.$$

In the exact powered-prime parameters, `w=dq` and
`a=ceil(log_2(10000(b+w+2)))=log_2 b+O(1)`. Consequently,

$$\begin{aligned}
k&\sim\frac b{64d\log_2 b},&
h&\sim(d-1)\log_2b,\\
\log_2\mathrm{Den}&\sim\frac b{64d},&
t&\sim\frac{k}{d-1}.
\end{aligned}$$

Taking `H=Y^(2k+t)` gives

$$Y=b^{32d+o(1)},\qquad
 \log_2H=\left(1+\frac1{2(d-1)}+o(1)\right)b. \tag{9}$$

Every fixed `d>=2` satisfies the eligibility conditions for sufficiently
large `b`. The unsigned orientation exponent in (8) is then

$$\left(\frac12+\frac1{4(d-1)}+o(1)\right)
  \frac b{\log_2 b}.$$

Increasing the fixed polynomial degree of the base therefore approaches
the height coefficient one and the unsigned decoding coefficient one-half.
The case `d=2` is the specialization displayed above.

### When the selector supplies coordinates

If a candidate is supplied as integers `(x,y)`, compute its target norm,
check `0<Q(x,y)<=H` and coprimality, and trial-divide over the factor
base. The coordinates need not initially be primitive: divide both by
`d=gcd(x,y)`, replacing the norm by `Q(x,y)/d^2<=H`. The normalized
pair is primitive, and its prime support remains in the same base.
The direct norm-to-word converter determines the prime orientations in
polynomial time; no meet-in-the-middle search is needed.
For reduced positive forms and `H=2^O(b)`, every valid coordinate has
`O(b)` bits, by the elementary ellipse bounds.

With any fixed `d` in (9), base construction, candidate verification,
word conversion, variable target powers, and weighted relation extraction
are all polynomial in `b`. Thus a coordinate selector with total cost

$$\exp\!\left(O\!\left(\frac b{(\log b)^\alpha}\right)\right)$$

would preserve that bound for every fixed `alpha>0`. The score-one
qualification above belongs to the **unsigned-norm decoding overhead**;
it is not a limitation of logarithm-weighted relations. This remains a
conditional interface statement, with the selector cost unresolved.

## 5. Exact finite controls

The checker
[`check_weighted_norm_relations.py`](../tools/check_weighted_norm_relations.py)
includes 33,280 exact least-power comparisons, 140 independently checked
nonsingular weighted matrices, and explicit rejection of invalid or
oversized product comparisons.

It also searches bounded norms for six actual discriminants, with
`H=2N`, and verifies every resulting relation identity in the class
group before extraction. Among the 146 rows, **31 fail ordinary
unweighted diagonal dominance**. All pass the logarithmic condition,
and all six inputs yield proper factors. This confirms that the
weighted contract accepts examples excluded by the previous wrapper.
The small cutoff and height choices in these finite runs are not
asserted to be uniform existence parameters.

The largest control is the previous 60-bit input
`N=1000000016000000063`, with cutoff `1021`. It uses 92 rows,
target exponents from 7 to 39, and 48,591 recorded coordinate points.
Twenty-three rows fail ordinary dominance; weighted extraction returns
the proper factor `1000000007`. Nonconsecutive partial row subsets
are checked separately with their designated target columns.

Run:

```sh
python3 tools/check_weighted_norm_relations.py \
  --output evidence/weighted-norm-relations.json
```

The [evidence](../evidence/weighted-norm-relations.json) includes the
exact norms, coordinates, signed words, target exponents, relation rows,
and extraction results.
