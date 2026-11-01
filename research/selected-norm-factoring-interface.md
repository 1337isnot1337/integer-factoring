# Factoring from selected integers without supplied representations

## Executable selector contract

The remaining arithmetic search problem can be exposed through a small
interface. A selector receives a target form and returns a positive
integer. It need not return its prime factorization, prime orientations,
or a representing vector. The exact consumer below reconstructs and
checks all of that data before using the resulting relation.

Let `N>1` be odd, `Delta=-4N`, and let `Y` be a prime cutoff. Construct
one chosen prime-form orientation `g_p` for every split prime `p<=Y`
coprime to the discriminant, using the
[deterministic root constructor](prime-form-root-construction.md).
Fix a nonnegative multiplicity budget `L` and set `K=L+1`.

For any target index `p`, an acceptable selector answer is an integer
`n` with all three properties:

- `n` is coprime to `Delta` and supported on the factor base;
- its total prime multiplicity is at most `L`;
- the form representing `g_p^K` primitively represents `n`.

The consumer can check these conditions using only `n`. It trial-divides
over the constructed base, then applies the complete
[factored-norm orientation decoder](factored-norm-orientation-search.md).
The decoder either proves that this fixed norm has no primitive
representation, or returns a signed prime word `v` and primitive
coordinates. In the successful case,

$$\phi(v)=g_p^K,\qquad \|v\|_1=\Omega(n)\le L<K.$$

Thus `K e_p-v` is a genuine strictly dominant relation. Different
accepted target indices give independent rows, including when some
other targets have no successful answer yet.

The implementation accepts an iterator of pairs `(p,n)`. It rejects
unsuitable candidates, permits a later retry for the same target, and
retains at most one relation per target. After any successful answer
it can invoke [partial relation halving](determinant-free-ambiguous-extraction.md).
Every returned divisor is checked exactly. A failed partial attempt has
only the torsion-image meaning proved in that extraction theorem.

If a candidate or a small discriminant prime already gives a proper
gcd with `N`, the consumer returns that factor directly. This is an
additional exact success path, not a hypothesis about selector quality.

## Full guarantee and its cost

Suppose the constructed prime forms generate the full class group,
and an acceptable norm is eventually supplied for every target. If
`N` has at least two distinct prime factors, the consumer returns a
proper divisor. Strict dominance gives a full-rank relation matrix;
the determinant-free extraction theorem supplies the conclusion.
It requires no promise that the supplied relations generate the full
kernel of the class map.

Under the pinned fixed-strip analytic input, the
[powered-prime theorem](powered-prime-smooth-norms.md) supplies explicit
uniform values of `Y` and `L` for which every target has an acceptable
answer. It also supplies full generation. The existence theorem does
not provide the selector or charge its running time by answer length.

Let `m` be the number of factor-base records, let `s` be the number of
candidate integers examined, and let `B` bound their written bit lengths.
After the base has been constructed, trial division takes

$$s\,(m+L+1)\,\mathrm{poly}(B+\log Y)$$

bit operations: a candidate is stopped as soon as its multiplicity
exceeds `L`, and every other unsuccessful prime test costs one division.
For candidates that pass these checks, orientation matching costs

$$s\,2^{\lceil L/2\rceil}\,
  \mathrm{poly}(\log |\Delta|+B+\log Y).$$

The row construction, target powers, and all extraction attempts are
polynomial in `m,L,log|Delta|`. Trying each accepted prefix adds at
most a further factor `m` to the extraction cost. The implementation
can instead extract once after the supplied sequence is consumed.

At the balanced powered-prime parameters, `m=2^O(sqrt(b))`,
`L=O(sqrt(b))`, and valid answers have `B=O(b)`, where
`b=bits(|Delta|)`. Thus any `2^O(sqrt(b))`-time selector producing
at most `2^O(sqrt(b))` candidates would complete a factoring procedure
at that scale. Without such a selector, this remains a checked
consumer of supplied answers. No stronger complete score follows.

## Implementation and independent replay

The consumer is
[`selected_norm_factoring.py`](../tools/selected_norm_factoring.py).
The checker is
[`check_selected_norm_factoring.py`](../tools/check_selected_norm_factoring.py).

It reuses the previously searched 60-bit and 92-bit norm fixtures but
passes only each target's prime index and represented integer into the
consumer. It withholds every old factorization, orientation, coordinate,
relation, determinant, and annihilator. All prime roots are constructed
again, so their chosen orientations can differ from the old fixture.
An integer represented by a form is also represented by its inverse,
which makes the unsigned norm fixture valid under that orientation change.

Every recovered vector is checked against its target form, including
primitivity and coprimality. Every recovered relation is independently
evaluated in the class group. The two proper factors are then checked
by integer division. Separate controls exercise unsuccessful candidates,
retries, a completed-target duplicate, early partial extraction, and
the direct gcd success paths.

Run:

```sh
python3 tools/check_selected_norm_factoring.py
```

The exact counts and sample recovered certificates are recorded in
[`selected-norm-factoring.json`](../evidence/selected-norm-factoring.json).
This replay measures decoding and extraction from selected integers.
The original cost of finding those integers remains a separate input.

A [separate larger replay](../tools/check_selected_norm_factoring_124_bit.py)
applies the same interface to the 124-bit partial-search fixture. It
decodes all 942 supplied integers and independently checks every recovered
relation before confirming the factor `2305843009213693951`. Its
[evidence](../evidence/selected-norm-factoring-124-bit.json) records the
source fixture hash and contains a sample reconstructed certificate.
This is a test of the same consumer at a larger size, with the norm
selection cost still charged to the original search.

## A polynomial factor base already suffices for a score-one interface

The balanced choice of parameters targets `2^O(sqrt(b))` total time.
There is also a simpler specialization with a factor base of fixed
polynomial size. In the powered-prime construction, choose

$$q=\lceil\log_2 b\rceil,\qquad v=64q,\qquad Y=2^v.$$

Then `Y=b^(64+o(1))`. Its explicit eligibility inequalities hold for
all sufficiently large `b`. With the notation of that construction,

$$k=(1+o(1))\frac{b}{128\log_2 b},\qquad
 h=(1+o(1))\log_2 b.$$

Indeed, `w=v/32=2q` and
`a=ceil(log_2(10000(b+w+2)))=log_2 b+O(1)`, so `h=w-a` has
the stated size. The density denominator satisfies

$$\begin{aligned}
\log_2\mathrm{Den}
 &=\log_2(b+2)+\lceil b/2\rceil-kv
     +k\log_2(32vk)\\
 &=(1+o(1))\frac b{128}.
\end{aligned}$$

Here `ceil(b/2)-kv=O(log b)` and `32vk=(16+o(1))b`.
Consequently the walk length and the final multiplicity obey

$$t=(1+o(1))k,\qquad
 L=2k+t=\left(\frac3{128}+o(1)\right)
       \frac b{\log_2 b}.\tag{1}$$

Every target has an acceptable norm with

$$\log_2 n\le Lv=\left(\frac32+o(1)\right)b.\tag{2}$$

Both the factor-base setup and all integer linear algebra now have
polynomial cost in `b`. A supplied valid norm can be decoded within

$$2^{(3/256+o(1))b/\log_2 b}\,\mathrm{poly}(b).\tag{3}$$

The polynomial number of target rows does not change that exponent.
Thus a deterministic selector with total cost `2^O(b/log b)`, allowing
at most that many candidates and charging their written sizes, would
give a complete score-one factoring algorithm through this interface.
Polynomial primality and perfect-power preprocessing, followed by
recursive splitting, are accounted for in the
[arithmetic framework](arithmetic-framework.md).
Equations (1)--(3) are guarantees about existence and processing of
supplied answers. The required selector is still open.
