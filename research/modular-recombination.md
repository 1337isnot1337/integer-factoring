# Effective recombinations of the complete modular pieces

## Result and scope

On the three-label characteristic-two domain at weight 12, there is an
explicit 17-piece crossing decomposition whose maximum row-rank times
column-rank is **18**. This improves the maximum 27 in PR #22's earlier
17-piece decomposition. The new pieces are sums of at most two symmetry
images of the original pieces. The exact coefficient lists and both ranks
are recorded in [modular-recombination.json](../evidence/modular-recombination.json).
Every listed piece is reconstructed, both ranks are checked by independent
elimination, and the complete crossing sum is replayed.

This is a finite algebraic improvement, **not a factoring exponent above
1/2**. The labels number three, rather than `2^lambda`, and there is no
proved favorable profile as the label dimension grows. The five-label
case improves from 210 to 154 in the same search. The eight-label case
in the larger suite improves from 2703 to 1935; even this is too large
to suggest a better leading coefficient on its own.

The search also gives a uniform effective closure operation: polynomially
many bounded-size recombinations can be searched, ranked, and assembled
into a crossing family within the existing `2^O(lambda)` construction cost.
Thus a future uniform width theorem for the selected pieces would suffice;
there is no extra factor-generation obligation hidden behind that theorem.

The underlying crossing operators are Claude (Anthropic)'s Appendix B,
pinned to [formal-math@16a5ce4](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The polynomial profile evaluator is reused unchanged from
[PR #22, commit 735b3dd](https://github.com/1337isnot1337/integer-factoring/blob/735b3ddb37d2e63dc469cf833687355e42097dd4/research/modular-pfaffian.md).
That note proves the six-block Pfaffian formula and its polynomial bit cost.
This work adds the recombination search, its effective closure argument,
exact witnesses, and a different family of structured-label experiments.
No new Lean formalization is involved. The research and checker used
substantial OpenAI assistant support; the existing Apache-2.0 attribution
is retained.

## 1. Symmetries and allowed candidates

Write the crossing as `E(x,y,u,v)=[x=v][y=u]`. Its eight coordinate
symmetries are exactly the permutations preserving the unordered matching
`{{0,3},{1,2}}`. Each symmetry fixes `E`, although it generally changes
an individual piece. Every symmetry either preserves or interchanges the
two flattenings, possibly transposing their matrices. Thus the product
of their ranks is unchanged by a symmetry.

For a domain with `D` labels, construct the complete matrices `A_s` from
the Pfaffian profiles. Flatten them in row order `(x,y)` and column order
`(u,v)`. The column flattening exchanges coordinates `y,u`.
Take all eight symmetry images of every nonzero piece, discard zeroes
and duplicate tensors, then optionally add every sum of two distinct
images. All sums and ranks are over `F_2`.

The original decomposition is included in the candidate list, so the
crossing is in its span. A selected subset of candidates need not be a
partition of the original summands: an exact linear equation determines
which candidates sum to `E`. Characteristic two makes each coefficient
zero or one.

## 2. Exact optimization within the candidate list

Give a candidate `F` cost `rank_row(F)*rank_col(F)`. Sort candidates by
cost and insert their flattened tensor vectors into a binary linear basis.
At the first cutoff whose span contains `E`, solve for a decomposition.
Record each chosen candidate as one or two original indices and symmetry
permutations, together with its two ranks. Reconstruct those matrices
from the original profiles and check their sum entry by entry.

This gives the **minimum feasible maximum rank product within the declared
candidate list**. The checker separately inserts every candidate of strictly
smaller cost and verifies that their span does not contain `E`. This does
not establish minimality over all combinations of three or more images,
all elements of the full span, or arbitrary tensor decompositions.

A nonzero binary matrix of rank `r` has an effective rank-`r` factorization.
Select independent columns `J` and rows `I` with invertible
`B=A[I,J]`; then `A=A[:,J] B^-1 A[I,:]`. Padding with zero channels gives
width `ceil(log2 r)`. A zero matrix uses one zero channel, of width zero.
Consequently a rank-product cutoff `Q` gives width sum at most
`log2 Q+2`. Rank product is the search objective; it does not equal the
sum of the two separately rounded binary widths.

## 3. Uniform construction cost

Suppose that `h(lambda)=poly(lambda)` and a fixed algorithm produces
`D=2^lambda` distinct-odd-part labels of that common weight within
`2^O(lambda)` time. There are `O(h)` source pieces. The profile formula
computes all their matrices using `D^4 poly(h)=2^O(lambda)` bit operations.
There are `O(h)` single-image and `O(h^2)` pair candidates. Each tensor
has `D^4` coefficient bits. Forming them, computing both matrix ranks,
sorting, solving the span problem, and obtaining rank factorizations all
have fixed polynomial cost in `D` and `h`, hence `2^O(lambda)`.

The selected decomposition uses at most the number of candidates, which
is polynomial in `lambda`. Its coefficients and factors belong to `F_2`,
and its exact crossing multiplier is one. All identities extend to every
commutative `F_2`-algebra. Selection and elimination happen over the base
field before evaluating an input table; no coefficient-algebra value is
inverted. These observations verify every effective modular-pieces
condition except the desired uniform width bound.

More precisely, if the selected cutoff satisfies
`log2 Q(lambda)<=4 rho lambda+O(log lambda)`, then the extra two rounding
bits fit the logarithmic remainder and the
[fixed-field transfer](scale-one-modular-transfer.md) gives factoring
exponent `log2(1/rho)`. The source family is included, so recombination
cannot worsen its rank-product cutoff. This does not prove a stricter
uniform coefficient. Taking tensor powers of a finite witness introduces
exponentially many choices; polynomial-piece compression and its ranks
must be analyzed separately.

## 4. Structured binary labels

The other experiment uses `lambda` independent equal-sum blocks. In block
`j`, choose `{a_j+1,a_j+7}` or `{a_j+3,a_j+5}`, where
`a_j=spacing*j+offset`, with `spacing>=8` divisible by eight and
`offset>=0` divisible by eight. The blocks are disjoint, all parts are
odd, and every label has common weight
`spacing*lambda*(lambda-1)+(2*offset+8)*lambda`.
For fixed spacing and offset this is quadratic in `lambda` and yields
exactly `2^lambda` labels. The profile formula therefore makes this a
fully effective family at every width.

The favorable one-bit width sum of two does not persist in the checked
cases. For the five tested spacings and offsets, width two has maximum
width sum seven or eight. Width three has maximum width sum twelve,
with maximum rank product 2500 or 3136. These are exact finite failures
of a tempting extrapolation, not an asymptotic obstruction to every label
encoding or to later recombinations.

## 5. Reproduce and evidence limits

```sh
python3 tools/check_modular_recombination.py --output evidence/modular-recombination.json
python3 tools/check_modular_recombination.py --full --output evidence/modular-recombination-full.json
make verify
```

The default suite independently compares 528 complete-piece entries with
the original normal-order recurrence at weight eight, replays the crossing
identities, checks the single-image and pair candidate cutoffs at weights
8, 12, and 16, and checks ten structured-label cases. The full suite adds
weights 20 and 24 and five three-bit structured cases. Removing a nonzero
selected piece is a negative control. Source hashes bind the evidence to
the profile evaluator and checker used.

The ordinary bound established by the repository remains `alpha=1/2`.
The next constructive questions are whether larger recombinations improve
the growing-domain profiles and whether an explicit family admits a
uniform rank-width proof. Neither is answered by the finite witnesses.
