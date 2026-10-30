# What would improve the square-root logarithm exponent

## Scope

This note studies changes beyond the Appendix B rate `alpha = 1/2`.
It proves restricted obstructions and records finite rank experiments.
It does not prove a factoring bound with `alpha > 1/2`, or a lower bound
on factoring algorithms.

The source is Claude (Anthropic), *A 2^o(n) Time Algorithm for 3-SAT*,
in `anthropics/formal-math` at
`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`.
The [paper](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/paper.pdf),
especially Appendices A, B and H, supplies the operator construction and
the recursion. The [evaluator audit](half-evaluator-audit.md) identifies
the assembled source theorem used by the factoring composition.

## 1. The two parameters that determine the exponent

For the existing recursive evaluator, the number of clique vertices satisfies

```math
K_{\ell+1}=g_\ell^2,\qquad
g_\ell=2^{3\ell+x_0}K_\ell.
```

Consequently `log K_l = Theta(2^l)`. Its leading cost coefficient contracts
by `rho` at each level. With a polynomial cutoff in the input length, the
logarithmic saving exponent is

```math
\alpha=\frac{\log(1/\rho)}{\log 2}.
```

Thus obtaining a larger exponent requires a smaller leading contraction,
a recursion with a genuinely smaller growth base, or a different method
for representing and evaluating the factoring instance. Reducing fixed
constants in the existing schedule does not change this equation.

The source's schedule and rate are stated in
[LevelsStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/LevelsStatement.lean)
and [AppRunningTimeStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppRunningTimeStatement.lean).

### A counting restriction on slower grid growth

Consider a rectangular `g` by `h` layout, with either orientation available
for the next recursive call. Suppose the row network is partitioned into
`K` groups, each of size at most `2g/K`, and the column network into `K`
groups, each of size at most `2h/K`. A nonempty network has a nonempty
group. Its size is at least one. Therefore

```math
g\ge K/2,\qquad h\ge K/2,\qquad gh\ge K^2/4.
```

This argument permits empty groups. Replacing the square grid by a
rectangle cannot reduce the growth base below two while retaining this
same two-orientation, same-child-size interface. More precisely, if
`K_(l+1) >= K_l^2/4` and `K_0 > 4`, then induction gives

```math
\log K_\ell\ge 2^\ell\log(K_0/4)+\log4.
```

The conclusion is restricted to that interface. It does not exclude a
different recursion whose two orientations produce different child
parameters, or one that exploits a property of the factoring tables.

## 2. Changing energy parameters does not improve the declared channels

The Appendix B space is

```math
V=\mathbb Q[t_1,t_2,\ldots],\qquad \deg t_j=j,\qquad
\dim V_E=p(E),
```

where `p(E)` is the partition number. Its logarithm grows as
`c sqrt(E) + O(log E)` for a positive absolute constant `c`.
The source encodes each `lambda`-bit string by a partition of the least
integer `h` with `p(h) >= 2^lambda`.

For the central piece, the declared row and column channel counts are both

```math
R_0=\sum_{j=0}^{2h}p(j).
```

Since `p(2h) <= R_0 <= (2h+1)p(2h)`, their total width is
`2 sqrt(2) lambda + O(log lambda)`. Keeping these channels cannot give
any fixed `rho < 1/sqrt(2)`. This is a statement about the declared
factorization. It is not a lower bound on the actual matrix ranks of the
pieces, which could have redundant channels.

An asymmetric choice of operator weight `h` and state weight `E` does
not fix this count. This is a direct extension of the source's `LemmaQ`,
which quantifies arbitrary `E`. For term `i`, the two operators have
shifts `s-h+i` and `h-s-i`. Applying the two possible orders to `V_E`
therefore gives channel counts

```math
R^{\rm row}(E,h,s)=\sum_{i=0}^{2h}p(E+s-h+i),\qquad
R^{\rm col}(E,h,s)=\sum_{i=0}^{2h}p(E+h-s-i),
```

where `p(j)=0` for negative `j`. Use partitions of `h` for the operator
labels `x,v` and partitions of `E` for the state labels `y,u`.
`LemmaQ.c`, summed over `|s| <= E+h`, and the unchanged dual pairing
give the crossing identity. To encode equally many operator labels and
state labels, both `p(h)` and `p(E)` must be at least `2^lambda`.
Let `h_lambda` be the least allowed weight. At `s=0` the displayed two
sums coincide, contain `p(E+h)`, and are at most `(2h+1)p(E+h)`.
Thus the leading normalized total width is at least

```math
\frac{2c\sqrt{E+h}}{4c\sqrt{h_\lambda}}
\ge\frac1{\sqrt2},
```

provided the weights remain polynomial in `lambda`, as required for the
same polynomial piece count. Equality in the leading term occurs at
`E = h = h_lambda`. The optional trivial width cap `2 lambda` does
not improve this minimum: a cap that becomes active has a larger leading
width than the balanced optimum.

## 3. Sparse partitions are not a drop-in replacement

Appendix H observes that a suitable operator construction with
`log dim V_E = Theta(E^a)` would suggest

```math
\rho=2^{a-1},\qquad \alpha=1-a.
```

This suggests using partitions into squares or other sparse sets of parts,
whose generating functions grow more slowly. The existing operators do
not preserve those smaller spaces.

For a single generator `t_n`, the source's definition gives exactly

```math
t_n[r]=\binom{r-1}{n-1}\,M_{t_r}\qquad(r\ge n),
```

where `M_(t_r)` denotes multiplication by `t_r`. There is only one mode
index in the defining sum, so it must be `-r`; the binomial coefficient
is nonzero over the rationals.

Therefore a subspace containing `1` and closed under every positive mode
of any one fixed `t_n` contains every `t_r` for `r >= n`. Repeated
applications also give every monomial in those variables. In particular,
its weight-`nq` component contains all monomials obtained by replacing
every part `j` of a partition of `q` by `nj`. Its dimension is at least
`p(q)`. This rules out a uniform growth estimate `exp(O(E^a))` with
`a < 1/2` for such a closed space.

For square parts, failure is immediate:

```math
t_1[2]1=t_2.
```

The right side is outside `Q[t_1,t_4,t_9,...]`. Keeping the old operator
identities while substituting square-partition dimensions would count
channels that do not span the intermediate results.

### Nonlinear regrading also breaks a required identity

Giving `t_j` the new weight `j^2` does not retain the source's property
that every `a[r]` has a definite weight shift. Direct expansion gives

```math
(t_1^2)[2]t_1=t_1^3+2t_3.
```

Both terms have original weight three. Their new weights would be three
and nine. The two intermediate-energy bounds used in the crossing
factorizations would consequently be false under that regrading.

These arguments exclude the indicated substitutions using the unchanged
operators. They do not exclude a new algebra, different operators,
input-dependent restrictions with a proved invariant space, or explicit
rank compression after the operator products are formed.

### Projecting away high modes breaks the full crossing identity

A different suggestion is to retain `Q[t_1,...,t_L]` and project each
intermediate result back into that ring. The smaller dimensions cannot
be substituted without checking the crossing identity. It already fails
for `h=L=2`, where the label domain has the two partitions `(1,1)` and
`(2)`.

Take all four labels to be `(1,1)`. Then `a=b=t_1^2`,
`a^<i>=binom(4,i)t_1^2`, and `h!=2`. Let `P_L` delete every monomial
containing a mode above `L`. Define the projected version of the original
operator order by

```math
\widetilde A_s=[t_1^2]\sum_{i=0}^4(-1)^i\binom4i
(t_1^2)[s-2+i]P_2(t_1^2)[2-s-i]t_1^2.
```

The finite exact calculation is:

| `s` | -4 | -3 | -2 | -1 | 0 | 1 | 2 | 3 | 4 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Original `A_s` | 2 | -8 | 28 | -56 | 70 | -56 | 28 | -8 | 2 |
| Projected `A_s` | -36 | 80 | -64 | 8 | 46 | -56 | 28 | -8 | 2 |

The original values sum to two, as required. The projected values sum
to zero. This is a counterexample to that projected construction, with
no floating-point calculation or asymptotic assumption.

For a direct reproduction of the nine entries, put

```math
F(r)=[t_1^2](t_1^2)[r](t_1^2)[-r]t_1^2.
```

For `r=-6,...,2`, the unprojected values are
`166,120,84,56,36,16,16,0,2`; inserting `P_2` gives
`0,0,8,24,12,16,16,0,2`. Values with `r>2` vanish by weight.
Taking `sum_i (-1)^i binom(4,i) F(s-2+i)` gives the table.
The [exact projection checker](../tools/check_half_projection.py) computes
these numbers from the operator recurrence, without using the expected
table to form its answer.

There is also a general obstruction to retaining all the source identities
after a finite cutoff. Assume `h >= 1` and `b` is homogeneous of weight
`h`. Write `a = sum_(j=0)^h a_j`, where `a_j` is homogeneous of weight
`j`, and fix input and output monomials `y,u` in the retained variables.
Let `E` be the weight of `y`. Put

```math
F_j(r)=[u]a_j[r]P_Lb[-r]y.
```

Each `F_j` has finite support. Indeed, `b` has at most `h` factors,
its differentiations occur before its multiplications, and every surviving
created mode has weight at most `L`. Its positive weight shift is at most
`hL`; its negative shift cannot exceed the input weight. Thus
`P_L b[-r] y` can be nonzero only for `-hL <= r <= E`.
The projected `Q_s` coefficient is

```math
\sum_{j=0}^{h}\sum_{i=0}^{2h}
(-1)^i\binom{h+j}{i}F_j(s-h+i).
```

All sums can now be interchanged. Summing over every integer `s` gives

```math
\sum_{j=0}^{h}
\left(\sum_{i=0}^{2h}(-1)^i\binom{h+j}{i}\right)
\left(\sum_{r\in\mathbb Z}F_j(r)\right)=0,
```

because `1 <= h+j <= 2h` and each binomial sum is `(1-1)^(h+j)`.
Therefore this projection cannot retain both the source's support bound
`|s| <= E+h` and its nonzero crossing sum. This observation applies to
every finite `L`, even if `L` grows with the input parameter. A different
construction could introduce boundary corrections or another operator
identity; those would require new proofs and cost bounds.

The exact definitions are in
[AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean)
(`term`, `opMon`, `op`, `O1`), with their proofs assembled in
[GeneralPieces/Main.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Main.lean).

## 4. More grid directions do not improve the suggested exponent by themselves

There is no higher-dimensional crossing identity supplied by this note.
The following checks the proposed parameter calculation only.

If a `d`-orientation generalization groups every orientation down to
the same `K`, with the same factor-two class-size allowance, its parent
has at least `(K/2)^d` vertices. Thus its growth base is at least `d`.
If its hypothetical energy calculation gave `rho_d = d^(a-1)`, then

```math
\frac{\log(1/\rho_d)}{\log d}=1-a.
```

With ordinary partition growth `a=1/2`, this is still `1/2`.
The larger local saving is exactly offset by the larger recursion growth
base in that proposed calculation. A useful higher-dimensional method
must beat one of these premises, not simply substitute a larger `d`.

## 5. A universal restriction within polynomially many crossing pieces

Suppose `M(lambda) = lambda^O(1)` pieces sum to a nonzero scalar multiple
of the crossing on `lambda`-bit words. The crossing has row matrix rank
`2^(2 lambda)` over the rationals. Subadditivity of rank implies that
some piece has row rank at least

```math
\frac{2^{2\lambda}}{M(\lambda)}.
```

Its row width is at least `2 lambda - O(log lambda)`. Column width is
nonnegative. Consequently a uniform sum-width estimate
`c_row + c_col <= 4 rho lambda + O(log lambda)` forces `rho >= 1/2`.
For the base-two recursion, this means that this particular interface
cannot certify `alpha > 1`.

This leaves the interval `1/2 < alpha <= 1` open. It does not say that
`alpha = 1` is attainable, and it places no restriction on factoring
algorithms outside this crossing interface.

## 6. Finite tests of actual half-piece ranks

A separate Python experiment constructs the pieces from the integer
operator recurrence (O2), rather than assuming their channel counts are
their ranks. It forms `a^mu`, the binomially weighted `a^<i>`, and
`b^nu = (h!/z_nu)t^nu`, and performs Gaussian elimination over a prime
field. Each pairing identity `a^mu[-h] b^nu = h! [mu=nu]` is checked
as an exact integer equality. Each entry of `sum_s A_s = h! cross`
is checked modulo the selected prime.

The [rank checker](../tools/check_half_piece_ranks.py) uses Python's
standard library; NumPy is optional and accelerates integer matrix
products. Its default weights are `2,3,4`. Run the larger recorded suite
with `python3 tools/check_half_piece_ranks.py --full`. The matrix products
have an explicit overflow bound when NumPy is used.
The [recorded rank results](../evidence/half-piece-ranks.json) contain all
13 cases and the checker hash. They cover 221 piece matrices, 500 integer
duality equalities, and 39,896 modular entries of the summed crossing
identity. The default three cases also passed without NumPy and matched
the accelerated output.

On the full partition domains, the central piece has these ranks modulo
65521:

| Weight `h` | Number of labels `p(h)` | Row rank of `A_0` | Column rank of `A_0` |
| ---: | ---: | ---: | ---: |
| 2 | 2 | 4 | 4 |
| 3 | 3 | 9 | 9 |
| 4 | 5 | 25 | 25 |
| 5 | 7 | 49 | 49 |
| 6 | 11 | 121 | 121 |

These are full ranks. A nonzero determinant modulo a prime is also
nonzero over the rationals, so each displayed full-rank conclusion is
exact for its finite matrix. It supplies no asymptotic lower bound.
The full partition domain is also distinct from the source's first
`2^lambda` lexicographic partitions when `p(h)` is not a power of two.
For the exact source domains at `lambda=1,2,3`, corresponding to
`h=2,4,6`, the central ranks are respectively `4,16,64` in both
directions, again full. A second prime, 1000003, gives the same full
central ranks for the five full-partition cases. The central column
factorization is also constructed independently from the opposite
operator order and checked modulo the selected prime against the
reshaped row factorization.

No smaller fixed-block contraction follows from these small cases.
Searching for a useful rank reduction remains a valid direction, but
small full-rank cases cannot establish that the large cases are full
rank or asymptotically attain the declared channel counts.

## 7. Special structure of the factoring tables

The [multiplication-network structure note](half-top-structure.md) proves
a quantitative restriction for the completed constructor. If `H` is the
number of free `faQ` bits in its largest group and `V=P^2 lambda` is the
channel volume, then every grouping with separate original-wire channels
satisfies `H+V/4>=n`. Literal local enumeration costs at least `2^H`, and
the specific explicit pad used here costs at least `2^(V/r)`.
The resulting complete construction costs
`2^Omega(n/sqrt(log n))` for the fixed Appendix B schedule.

The proof permits arbitrary adaptive grouping and separate rational
factorizations of the original wires. It does not apply to joint boundary
compression, a local solver that avoids Cartesian enumeration, or a
smaller sufficient pad proved for special inputs. It is not a lower bound
on factoring or on every use of the evaluator.

## 8. Restricting the labels without changing the operators

Keeping the original operators and encoding labels by partitions with
parts at most `L` preserves the crossing identity. This differs from the
projection tested in Section 3. However, the [restricted-label note](half-restricted-labels.md)
gives an exact counterexample to the proposed intermediate-span bound
`p_{<=L}(2h)`: at `h=L=2`, four input pairs have four independent images,
while the proposed low-mode product space has dimension three. A
four-by-four integer minor has determinant four.

This refutes the literal dimension substitution. It does not exclude
asymptotic compression, and independence of a factor block does not imply
the same rank after multiplication by the other factor.

## 9. Other published evaluator and remaining mathematical target

The
[OpenAI ETH manuscript](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/build/main_full_version.tex)
at commit `4eaee87548250b2682affda7d05d947c05554bae`
states an exact-counting runtime with logarithmic saving exponent
`1/2`, including `2^O(n/sqrt(log n))` for a fixed linear clause bound.
Its main theorem therefore does not supply a stronger exponent for this
factoring effort by direct substitution. This is a comparison of stated
bounds, not an audit of that separate proof.

An improvement beyond the square-root rate still needs one concrete
new result. Examples are a uniform crossing decomposition with polynomially
many pieces and `rho < 1/sqrt(2)`, a useful rank compression of the actual
Appendix B tensors, or a different recursion that avoids the grid
growth restriction proved above. Each must include the sizes of its
integers, an effective construction, and its complete evaluation cost.
No such result is established here.
