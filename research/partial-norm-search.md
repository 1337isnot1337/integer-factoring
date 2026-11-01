# Continuing norm search after a target fails

The finite norm-search prototype can use the relations it has already
found, even when another target has no witness within its current
search bounds. The new implementation skips that target and continues.
It periodically applies
[determinant-free relation halving](determinant-free-ambiguous-extraction.md)
to the independent rows collected so far.

This change can reduce the number of successful norm selections needed
for a particular input. It does not give a new bound for uniform norm
selection, and it does not improve the complete factoring score.

## Why a partial set of norm witnesses is enough to try extraction

Let `g_1,...,g_m` be the chosen prime forms of discriminant `-4N`.
Fix a multiplicity budget `L` and set `K=L+1`. Suppose norm search
succeeds for distinct target indices in a set `I`, producing words

$$\phi(v_i)=g_i^K,\qquad \|v_i\|_1\le L\quad (i\in I).$$

The resulting relation row is `r_i=K e_i-v_i`. Restrict these rows to
the columns indexed by `I`. In each row the diagonal entry has
absolute value strictly larger than the sum of the other absolute
entries, since

$$K-|v_{i,i}|>\sum_{j\ne i}|v_{i,j}|.$$

Thus the designated square minor is nonsingular, and the full rows
are independent over the rationals. No determinant computation is
needed to establish this fact.

Write `R` for the rectangular relation matrix, `L_R` for its integer
row lattice, and

$$L_R^{\mathrm{sat}}=(\mathrm{span}_{\mathbb Q}L_R)\cap\mathbb Z^m.$$

The halving extractor searches for useful involutions in the image
of `L_R^sat/L_R`. It either returns a proper factor or certifies that
this particular torsion image has no useful involution. The latter
outcome does not obstruct adding a new independent relation row.
The relevant torsion images can only increase as rows are added.

If the prime forms generate the full class group, the deficiency
corollary gives an additional guarantee: with `r` accepted rows,
failure implies that `N` has at most `m-r+1` distinct prime factors.
The finite prototype does not certify full generation at its small
chosen cutoffs. Its returned factors are nevertheless exact, because
each relation and each final divisor is checked directly.

## Finite search and extraction policy

[`arithmetic_partial_norm_factor.py`](../tools/arithmetic_partial_norm_factor.py)
uses the same complete ellipse sieve as the
[original search prototype](complete-smooth-norm-sieve.md).
For a reduced target `(A,B,C)`, its first norm cap is `4 max(A,C)`;
each unsuccessful expansion multiplies that cap by four.

When the configured number of expansions is exhausted, the program
records that target and proceeds to the next generator. A successful
norm is converted to a prime word and evaluated again in the class
group before its row is accepted. Extraction runs after a configured
number of additional successful rows and once at the end.

The program records every visited target, cap, witness, failed search,
and extraction attempt. It does not compute a relation determinant,
an annihilator, generator orders, or a complete class group. The
original cost of searching the ellipses is retained in the output.

For example:

```sh
python3 tools/arithmetic_partial_norm_factor.py 10807 \
  --smooth-bound 31 --maximum-expansions 0 \
  --extraction-interval 1 --keep-relations
```

This input has four chosen prime forms. Search fails at two target
caps, but the two successful rows still return the proper factor `101`.
The partial matrix has pivot columns `1,3`, using zero-based indices:

$$R=\begin{pmatrix}2&15&1&0\\0&0&0&16\end{pmatrix}.$$

The two failed ellipses contain no primitive coprime smooth witness
within their respective caps, as independently checked by enumeration.
The second successful norm is `23`, and its row already supplies the
cyclic relation that yields a useful involution after halving.

## Exact replay

[`check_partial_norm_factoring.py`](../tools/check_partial_norm_factoring.py)
checks both fresh finite searches and stored larger certificates.
It recomputes target powers, checks the integer norm equations and
coprimality, refactors each represented norm over the independently
rebuilt base, reconstructs the oriented words, evaluates their classes,
and reconstructs the dominant relation rows. It then repeats partial
extraction and checks the returned divisor by integer division.

The small failed-target controls also enumerate their empty ellipses
independently of the sieve. The supplied larger certificates replay
the successful representations and extraction; they do not repeat
the original coordinate search or re-establish its reported counters.

The fresh controls include these results:

| Input | Prime cutoff | Accepted rows | Skipped targets | Coordinate points | Proper factor |
| --- | ---: | ---: | ---: | ---: | ---: |
| `1022117` | 127 | 4 | 0 | 36 | `1009` |
| `100160063` | 127 | 8 | 0 | 259 | `10007` |
| `10002200057` | 127 | 19 | 0 | 744 | `100003` |
| `10807` | 31 | 2 | 2 | 24 | `101` |
| `13493` | 13 | 1 | 1 | 7 | `103` |

The first three runs extract every four accepted rows. The last two
extract after each accepted row and allow only the initial norm cap.
Separate controls cover an unresolved prime input, an empty finite
factor base, a zero multiplicity budget, and a direct factor-base gcd.

### A 124-bit searched certificate

The larger stored record uses

$$N=21267647932558653899591465697288388633,$$

with prime cutoff `16381`, multiplicity budget `124`, and 942 prime
forms. Its first pass allowed eight norm-cap expansions per target.
That pass produced 930 rows and twelve exhausted targets, after
163,494,260 recorded coordinate points. The corresponding partial
torsion image had no useful involution.

Only the twelve failed target searches were then extended. Each
found a witness at its next cap. The successful target searches are
independent, so the original generator traversal and extraction policy
were replayed using the cached search results. The final record is
the deterministic policy with the enlarged expansion budget, including
its original per-target caps and counters.

All 942 completed rows yield

$$N=2305843009213693951\cdot9223372036854775783.$$

The final stored search visits 168,843,556 coordinate points. Every
represented norm has at most 83 bits and at most thirteen prime
factors counted with multiplicity. Relation halving finds a useful
involution on its first kernel pass, with 4,811 group multiplications
after the row identities have already been established. It uses no
determinant, annihilator, or full class-group enumeration.

The replay checks all 942 primitive representations and reconstructs
all 942 class identities before repeating extraction. This extends
the earlier 60-bit and 92-bit finite examples. It is a concrete
searched certificate, not a bound for selecting norms on arbitrary
inputs of that length.

Run the checker with:

```sh
python3 tools/check_partial_norm_factoring.py
```

The output is stored in
[`partial-norm-factoring.json`](../evidence/partial-norm-factoring.json).
Any larger supplied search records are stored separately in
[`partial-norm-factoring-certificates.json`](../evidence/partial-norm-factoring-certificates.json).
The complete algorithmic bottleneck remains a uniform deterministic
method for selecting enough suitable norms within the desired time.
