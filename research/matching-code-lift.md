# A multigraded crossing on uniform set partitions

## Status and provenance

This note gives an effective family of exact crossing decompositions over
`F_2` with `exp(O(lambda/log lambda))` coefficient pieces. The new step is a
division-free lift of binary mode positions to independent variables, using
uniform set partitions as labels. It does **not** prove a contraction below
`1/sqrt(2)`. The exact finite rank tables below show why the smallest example
cannot supply that conclusion.

The source operators and their original crossing identity are from
[Anthropic Appendix B](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5),
revision `16a5ce41304fcdaea4e0be7c132e1c361dc26da5`. The six-block folded
Pfaffian for the original translated duals is prior work in
[PR22](https://github.com/1337isnot1337/integer-factoring/pull/22), revision
`c12da897fcc6081685e8c663217ce9083db87b4a`,
`research/modular-pfaffian.md`. Its one-variable formula is used here as an
input, not claimed as a new result. The accompanying checker also verifies
the folding by explicitly expanding the translated dual and applying the
independent [monomial Wick evaluator](fermionic-piece-polynomial.md).

## 1. Labels with factorial growth and short binary descriptions

Fix a block size `k >= 2`. Let `a = k ell`. A label is an unordered partition
of `{0,...,a-1}` into `ell` blocks of size `k`. Encode a block `e` by its bit
mask and use the odd mode

\[
 n_e=1+2\sum_{i\in e}2^i.
\]

The label monomials have distinct odd parts and common weight

\[
 h=\ell+2(2^a-1).
 \tag{1}
\]

Their number is

\[
 d=\frac{a!}{(k!)^{\ell}\ell!},\qquad
 \log_2 d=(1-1/k)a\log_2a-O_k(a).
 \tag{2}
\]

For `k=2`, these are perfect matchings. The numerical energy in (1) is large,
but its binary description and the number of particles are `O(a)`.

## 2. The lifted Pfaffian

Use variables `L,Z_0,...,Z_(a-1)`, and write `Z^e=prod_(i in e) Z_i`.
For labels `x,v,u,y`, make one Wick vertex for each block of each label, in
group order `U,A,B,Y`. The six upper off-diagonal blocks are

| Groups | Entry, for the blocks `e,f` in group order |
|---|---|
| `U,A` | `L(1+L) Z^(e union f)` |
| `U,B` | `[e=f] (1+L) prod_(i in e)(1+Z_i)` |
| `U,Y` | `[e=f] L Z^e` |
| `A,B` | `[e=f] L Z^e` |
| `A,Y` | `(1+L) Z^(e intersection f)` |
| `B,Y` | `[e intersection f is empty] L(1+L) Z^f prod_(i in e)(1+Z_i)` |

The diagonal blocks are zero. In characteristic two, fill the lower blocks
by symmetry. Let `P_(x,y,v,u)` be the Pfaffian, computed by its division-free
sum over perfect matchings of the `4 ell` Wick vertices.

Substituting

\[
 L=z,\qquad Z_i=z^{2^{i+1}}
 \tag{3}
\]

gives the folded source polynomial. Indeed, equal-size bit blocks satisfy
`e subset f` if and only if `e=f`; the two relevant Lucas binomial conditions
therefore become equality and disjointness. The exponents `n_e | n_f` and
`(n_e & n_f)-1` give precisely the union and intersection in the table.
The factor `(1+z)^(n_e)` becomes `(1+L)prod_(i in e)(1+Z_i)`.

More directly, the coefficient sum is an exact crossing without relying on
any specialization argument. At

\[
 L=Z_0=\cdots=Z_{a-1}=1,
\]

only the `A,B` and `U,Y` blocks survive. Their equality matrices have
Pfaffian

\[
 P_{x,y,v,u}(1,\ldots,1)=[x=v][y=u].
 \tag{4}
\]

Thus all coefficient tensors of `P` sum to the identity crossing. No
division, limiting process, or nonzero-denominator assumption occurs at
this evaluation point.

## 3. Piece count and effective construction

Each physical bit position occurs exactly once in each external label.
Inspection of the six blocks shows that an occurrence contributes at most
one power of its `Z_i`. Consequently

\[
 \deg_{Z_i}P\le4.
 \tag{5}
\]

Each Wick matching has `2 ell` edges, each of low degree at most two, so

\[
 \deg_L P\le4\ell.
\]

The number of coefficient pieces is therefore at most

\[
 (4\ell+1)5^a=\exp(O_k(a)).
 \tag{6}
\]

A subset Pfaffian recurrence has `2^(4 ell)` states. Its polynomials have
at most the number of monomials in (6). Even naive multiplication within
this degree box costs `exp(O_k(a))`. This gives a uniform entry algorithm
with that running time; it does not enumerate permutations of all bit
positions. The separate permutation-orbit optimizer in the checker is only
for the finite tables below.

To obtain exactly `2^lambda` labels, take the least multiple `a` of `k` with
`d >= 2^lambda` and restrict to a deterministic prefix. The next-step ratio
in (2) is polynomial in `a`, so

\[
 a=O_k(\lambda/\log\lambda),\qquad
 \log_2\#\text{pieces}=O_k(\lambda/\log\lambda).
 \tag{7}
\]

Full tensor tables and Gaussian rank factors can be generated in
`2^(O_k(lambda))` time. This is an admissible piece budget for the
[subexponential piece transfer](subexponential-piece-transfer.md), with
entropy exponent `delta=1`. Applying that transfer still requires a
uniform joint row/column width bound, which this note does not supply.

## 4. Exact finite ranks

Permuting the `a` bit positions permutes the labels and high variables.
Thus ranks depend only on the low degree and the sorted list of high
degrees. Both flattenings use row/column label order consistent with (4).

| `a` | `k` | Labels | Nonzero degree orbits | Maximum rank product | A worst rank pair |
|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 1 | 5 | 1 | `(1,1)` |
| 4 | 2 | 3 | 51 | 18 | `(6,3)` |
| 6 | 2 | 15 | 542 | 32,574 | `(183,178)` |
| 6 | 3 | 10 | 101 | 2,116 | `(46,46)` |

For `a=4,k=2`, a worst degree is
`L^4 Z_0^2 Z_1^3 Z_2^3 Z_3^3`, up to permutation. There are 307 nonzero
coefficient pieces before passing to degree orbits. For `a=6,k=2`, a worst
degree is `L^6` with high degrees `(2,2,2,3,3,3)`. For `a=6,k=3`, the worst
degree is `L^4 prod_i Z_i^2`.

The first nontrivial example has a small rank product, but the next two
tests do not meet the below-half target. These are exact finite facts.
They neither establish nor rule out a favorable all-size rank bound of
the form `d^c exp(O(a))`.

The univariate source specialization also illustrates why coefficient
collisions matter: for `a=6,k=2`, some source coefficients have both ranks
`225`, whereas the lifted maximum is smaller. This statement is a separate
finite observation; the companion evidence records the multigraded ranks.

## 5. The homogeneous pairing on matching labels

For `k=2`, the homogeneous vacuum pairing is the determinant over `F_2` of
the edge-disjointness matrix between two matchings. If their union has
alternating cycles with `r_1,...,r_t` edges from each matching, let `A_r`
denote the relation of a single nontrivial `r`-cycle and otherwise common
edges. Then the pairing matrix is exactly

\[
 G_\ell=[\ell\text{ even}]I+
          \sum_{\substack{3\le r\le\ell\\r\text{ odd}}}A_r.
 \tag{8}
\]

To prove (8), write the disjointness matrix as `J+C`. A common-edge block
of `C` is `[1]`; a nontrivial cycle block is `I+S`, with `S` its cyclic
permutation. With two nontrivial cycles, `C` has nullity at least two and
the rank-one update `J` cannot make it invertible. With exactly one such
cycle, its adjugate is the all-ones matrix. The determinant update is its
cycle length modulo two. With no nontrivial cycle, the matrix is `I+J`,
whose determinant is `1+ell` modulo two. This proves every case.

The full finite Gram ranks for `a=2,4,6,8` are respectively `0,3,4,105`.
These finite checks establish the small cases only. The companion
[Jack-integrality proof](matching-gram-unipotence.md) now establishes
invertibility for every even `ell`, and gives an explicit finite polynomial
inverse. That all-size conclusion uses the proof rather than an extrapolation
from this table. The standard association-scheme framework is described in [Srinivasan (2020)](https://doi.org/10.5802/alco.104).

## 6. Reproduction

Run

```sh
python tools/check_matching_code_lift.py --output evidence/matching-code-lift.json
```

The checker verifies all 175 quartet-orbit representatives in the table
against explicit, unfurled translated-dual monomial Wick sums. It compares
all canonical coefficient matrices for `a<=4` with direct enumeration,
checks 58 scalar source-recursion entries in both operator orders, and
checks 11,260 pairing entries against independent determinants. All ranks
are exact Gaussian ranks over `F_2`.

The implementation uses mixed-radix bit packing only after the coordinate
degree bounds (5) have been proved. It copies symmetry-equivalent entries
with logical OR; treating those copies as new summands and XORing them
would be incorrect.
