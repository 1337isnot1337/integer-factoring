# Binary code restrictions: an affine-code obstruction

## Result and scope

Full-word barriers do not cover arbitrary subsets of word labels. One
natural family uses all affine Boolean functions `f_(a,b)(t)=a dot t+b`
on `F_2^r` as words of length `m=2^r`. There are `D=2^(r+1)` labels, so
the binary label width is `lambda=r+1`. Restrict the all-one-weight
[23-term Pfaffian profile](hamming-tensor-power.md) to these words.
It retains scale-one crossing, has a constant bound on the number of
coefficient positions, and is constructible in `2^O(lambda)` time despite
its exponentially large degree.

It nevertheless fails strongly. Every decomposition of the crossing
inside its complete coefficient symmetry span, over any characteristic-two
extension field, has a piece with **both** ranks at least

```math
\frac{(2^r-1)(2^r-2)}6-6D.                              (1)
```

For `r>=9`, this is at least `D^2/48`, forcing leading contraction
`rho>=1`. This particular code gives no leading contraction in the
existing evaluator. Arbitrary collective recombinations and any piece
count are allowed. This is a construction obstruction, not a factoring
lower bound. Other codes and local tensors are outside it. Substantial
OpenAI assistant support was used; independent review is pending and no
new Lean formalization is added. The factoring exponent remains `1/2`.

## 1. Uniform affine-flat histograms

For four affine words, the joint evaluation map has image an affine flat
`A` in `F_2^4`. Its dimension `q` is the rank of the four linear parts;
each image point has `2^(r-q)` preimages. A local term contributes

```math
L_M(z) A_M(0000;z)
 \left(\prod_{p\in A} A_M(p;z)\right)^{2^{r-q}}.           (2)
```

This is the existing polynomial restricted to the code, not a projection
of the original source operators. The last power is a Frobenius
substitution. Each local entry has degree at most four, `|A|<=16`, and
the fixed low/leading factor has degree at most twelve. For `r>=4`, all
exponents belong to

```math
\{b+a2^{r-q}:0<=b<=12,\ 0<=a<=64,\ 0<=q<=4\}.
```

Thus at most 4,225 coefficient positions occur, independently of `r`.
Generating all words, all `D^4` profiles, and exact dense factorizations
still costs `2^O(lambda)` bit operations. This effectiveness statement
alone supplies no favorable rank bound.

## 2. Clear the rank-three support

The 23 local support sets have sizes four, five, or eight. None contains
the full four-dimensional flat. Of the 30 affine hyperplanes, only two
are contained in any local support:

```text
y=v: {0,2,5,7,8,10,13,15}   (terms 0 and 5)
u=v: {0,3,4,7,8,11,12,15}   (terms 13 and 14)
```

The index bits are `x,y,u,v`; term order comes from `local_terms()`.
The checker enumerates the complete constant table. Consequently the
profile vanishes on joint linear rank four. On rank three it is

```math
G_r(z)[y=v]+H_r(z)[u=v].                                 (3)
```

Equality means equality of complete affine functions. A rank-three image
has exactly one defining affine equation. If this is a coordinate
equality, its histogram is the entire corresponding hyperplane with
uniform multiplicity; hence the surviving profile depends only on which
equality occurs. Complement equalities and other equations give zero.

After a crossing-preserving coordinate permutation, these supports lie
among the six coordinate-equality hyperplanes. For each pair `i<j`, fix
a rank-three tuple whose only equality is `i=j`. Evaluation there defines
a linear functional `g_ij` on the entire coefficient symmetry span.
Every member `F` has a corrected tensor

```math
F'=F-\sum_{i<j}g_{ij}(F)[label_i=label_j]                 (4)
```

vanishing on every tuple of joint linear rank at least three. This is a
well-defined map on actual tensors. An equality tensor has rank at most
`D` in either cut: its rank is one if both tested coordinates are on the
same side, and `D` otherwise. The correction has both ranks at most `6D`.
All `g_ij` values on the crossing vanish, so the crossing is unchanged.
The correction is essential; uncorrected profiles need not have the
block support used next.

## 3. Many separate two-dimensional blocks

Changing variables by `GL_r(F_2)` permutes free bit positions. Thus every
profile coefficient, its coordinate symmetries, and (4) are invariant.
In the row cut retain only pairs `(x,y)` with independent linear parts,
and similarly independent column pairs `(u,v)`. Index each pair by its
two-dimensional span `W`. Different row/column spans give joint rank at
least three, so the retained corrected matrix is block diagonal in `W`.
Each block has 24 rows and columns: six ordered bases and four choices of
constant parts. The same holds in the other cut, on pairs `(x,u),(y,v)`.

Fix independent linear parts `a,b`, with constant parts zero. The common
functional

```math
ell(F)=F'(a,b,b,a)                                       (5)
```

equals one on the crossing. Some piece of every decomposition is
detected by (5). Transitivity on ordered independent pairs makes its
corrected block nonzero for every `W`, in both cuts. The same entry has
independent pairs in both arrangements; no unrelated separators are used.

There are `(2^r-1)(2^r-2)/6` such subspaces, each having six ordered bases.
Nonzero blocks have rank at least one. Subtracting the allowance `6D`
proves (1). Its right side equals `D^2/24-25D/4+1/3`, which is at least
`D^2/48` for `r>=9`. Both widths are therefore at least
`2lambda-log2(48)`, giving `rho>=1`.

## 4. Replay and the next domain question

```sh
python3 tools/check_affine_code_barrier.py
python3 tools/check_affine_code_barrier.py --output evidence/affine-code-barrier.json
make verify
```

The replay checks the full hyperplane table, correction identities,
independent source-profile comparisons, and block support in a complete
small code. The arbitrary-size conclusion uses the proof above.

Quadratic Boolean evaluation codes remain outside this obstruction.
They have `1+r+binom(r,2)` independent coefficient bits and length `2^r`.
For arbitrary width `lambda`, take the first `lambda` multilinear
monomials at the least `r` with enough monomials. Evaluation is injective
by Boolean polynomial inversion. Then `r=O(sqrt(lambda))`; raw degree
and piece count are `2^O(sqrt(lambda))`. Their logarithmic piece cost
fits the [subexponential-piece transfer](sparse-width-interpolation.md)
pinned there to PR #23. Generation and dense factorization remain
`2^O(lambda)`; a uniform contraction theorem is still missing. Exhausting
their entire coefficient span is not justified by that budget.

There is nevertheless a growing-support search within the budget. From
`M<=2^(A sqrt(lambda))` coefficient/symmetry images, enumerate all sums of
at most `K=floor(c sqrt(lambda))` images, for any fixed `c`. There are at
most `(M+1)^K=2^O(lambda)` candidates. Compute both exact ranks, sort their
integer products, and insert candidates into a linear basis until the
crossing enters the span. Dense factorization and all membership tests
cost `2^O(lambda)`. At most `M` independent candidates are retained, so
the resulting witness still has `2^O(sqrt(lambda))` pieces. This finds the
exact best cutoff within the specified candidate set, not the whole
coefficient span. If the crossing never enters, it supplies no witness.
Thus quadratic-code experiments can search combinations with an unbounded
number of source images without violating the uniform construction cost.
