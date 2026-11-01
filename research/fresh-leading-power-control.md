# A fresh 124-bit consecutive-power factoring control

The fixed-base algorithm from
[smooth leading-power selection](smooth-leading-power-selection.md)
also factors a 124-bit semiprime generated from a fixed public seed:

\[
\begin{aligned}
 N&=16690589162031198237853308486068400277,\\
 N&=4509576403369997671\cdot3701143448763469987.
\end{aligned}
\]

Both factors have recursively verified prime certificates. They are
unrelated to the Mersenne-based controls used elsewhere in the repository.
The search function receives only \(N\), the cutoff, and the window;
the factors are used to construct and verify the test input.

## Reproducible input

The seed is the UTF-8 string
`integer-factoring-leading-power-fresh-2026-11-01`.
For each suffix `:left` and `:right`, interpret its SHA-256 digest as a
big-endian integer \(h\), and start at

\[
 a=2^{61}\mathbin{\mathrm{or}}
       (h\mathbin{\mathrm{and}}(2^{61}-1))
       \mathbin{\mathrm{or}}1.
\]

This gives starting candidates 4509576403369997599 and
3701143448763469981. The first primes at or above them are the two factors
displayed above, after 36 and three increments by two, respectively.
The checker verifies all 39 rejected odd candidates have a small-prime
or strong probable-prime-test compositeness witness. It verifies the
chosen primes with the complete recursive Lucas certificates, not with
a probable-prime acceptance test.

Input construction is separate from the factoring experiment. The exact
search is called as:

```python
factor_with_smooth_leading_powers(
    16690589162031198237853308486068400277, 4093, 100000)
```

## Exact search result

| Quantity | Value |
|---|---:|
| Input bit length | 124 |
| Prime cutoff | 4093 |
| Window per target | 100,000 |
| Prime-form generators | 262 |
| Accepted independent relations | 248 |
| Tested leading norms | 2,301,916 |
| Modular smoothness tests | 2,301,916 |
| Known-prime divisions after the filter | 48,870 |
| Largest tested exponent | 51,101 |
| Largest tested norm bit length | 63 |
| General integer-factoring calls in the search | 0 |
| New primes introduced | 0 |
| Final presentation deficiency | 14 |
| Involutions evaluated in the final extraction | 2 |
| Identity saturations in the final extraction | 0 |
| Returned factor | 3701143448763469987 |

Every target reached before the factor has a smooth leading norm inside
its window. The algorithm factors before reaching the remaining 14 base
generators. It performs no determinant computation or coordinate search.

The partial relation matrix, all selected powers and factorizations, and
the complete operation counters are stored in
[fresh-leading-power-control.json](../evidence/fresh-leading-power-control.json).
The checker reruns the full deterministic search, including rejected
candidates, and independently verifies every accepted class relation and
the weighted independence hypotheses.

```sh
python3 tools/check_fresh_leading_power_control.py
```

This is an additional finite exact control. The uniform short-window
bound required for an improved factoring exponent remains open. The
per-candidate and full-search bit costs are given in the main power-scan
note and include every rejected candidate.
