# An integer-only selector interface with no multiplicity promise

The [weighted relation theorem](height-only-smooth-norm-relations.md)
permits a simpler selector contract. Choose a norm-height bound `H`.
For each factor-base prime `p`, set `K_p` to the least positive integer
with `p^(K_p)>H`. A selector need only supply a positive smooth integer
`n<=H` that the target `g_p^(K_p)` primitively represents. It supplies
no factorization, orientation, coordinates, or multiplicity bound.

The new consumer
[`height_selected_norm_factoring.py`](../tools/height_selected_norm_factoring.py)
implements that contract. It constructs the prime roots, factors each
candidate over the base, and runs the complete orientation decoder.
A successful answer supplies a signed word `v`, with

$$\sum_q |v_q|\log q=\log n< K_p\log p.$$

The resulting row `K_p e_p-v` is strictly dominant after the column
for prime `q` is scaled by `log q`. The weighted extraction wrapper
checks this with bounded integer products. It then applies the same
determinant-free relation halving as before. All returned divisors are
checked exactly.

The program permits failed candidates, later retries, and partial
extraction after each accepted target. At most one relation is retained
per target prime. A failed partial extraction retains the torsion-image
scope of the relation theorem.

## Processing cost

Write `B=bits(H)` and let `r=omega(n)` be the number of distinct
prime factors of an accepted norm. Since its distinct primes are
at least `2,3,...,r+1`, one has `r!<=n<=H`. Consequently,

$$r=O\!\left(\frac B{\log(B+2)}\right).$$

The decoder uses `2^(ceil(r/2))` times a polynomial in the written
data sizes. Repeated copies of a prime are handled by powering and
do not add orientation choices. Trial division over a base of size
`m` costs `(m+B)` polynomial-size divisions per candidate.

Thus, with `B=O(b)` and a polynomial factor base, all processing
after norm selection fits `2^O(b/log b)`. The
[polynomial-base powered parameters](selected-norm-factoring-interface.md)
prove that every altered target has an acceptable norm at such a
height. A selector at that running-time scale would therefore give
score one. This implementation does not supply that missing uniform
search bound.

## Fresh finite searches

The [checker](../tools/check_height_selected_norm_factoring.py) computes
new target exponents from `H`, searches the corresponding exact ellipses,
and passes only each selected integer to the consumer. It uses a sieve
multiplicity parameter `bits(H)`, which is vacuous because `2^Omega(n)<=n`.
The old fixed-exponent representations are not inputs to these searches.

It independently checks each recovered vector, the least-power condition,
the weighted dominance inequality, and every relation identity. The
evidence records all coordinate visits, modular strikes, and prime
divisions, as well as target powers and selected integers. These search
costs are included explicitly and are not replaced by witness length.

| Input bits | Height `H` | Target rows | Exponent range | Rows requiring logarithmic weights | Coordinate visits |
| --- | --- | --- | --- | --- | --- |
| 60 | `2^42` | 92 | 5--27 | 65 | 53,238 |
| 92 | `2^64` | 284 | 6--18 | 196 | 8,926,750 |

Every target in both finite bases received an acceptable norm. The
integer-only consumer recovered factors `1000000007` and `2147483647`,
respectively. Its largest distinct-prime supports were seven and eight.
The independently verified rows requiring logarithmic weights fail the
ordinary unweighted diagonal-dominance inequality. They provide concrete
controls for the expanded contract, beyond merely choosing a smaller
target exponent for an already accepted row.

A separate control uses `N=103`, `H=7^20`, and the target `g_7^21`.
The consumer accepts `n=7^16`, recovers primitive coordinates, and
uses only one orientation bit despite multiplicity sixteen. This
checks the distinction between the removed multiplicity promise and
the distinct-prime decoding cost.

Run:

```sh
python3 tools/check_height_selected_norm_factoring.py
```

The exact results are in
[`height-selected-norm-factoring.json`](../evidence/height-selected-norm-factoring.json).
They test the relaxed contract and its full processing chain; they
establish no stronger complete factoring exponent.
