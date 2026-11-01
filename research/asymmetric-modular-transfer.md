# A weighted two-cut bound gives a complete rectangular factoring transfer

## Result and scope

A fixed supporting line for the two crossing widths is sufficient for a
uniform rectangular recursion. The two orientations may use different child
clique sizes. This note supplies the grouping, padding, construction, and
bit-cost arguments that are separate from the geometric calculation in
[asymmetric-support-line-crossings.md](asymmetric-support-line-crossings.md).

**Conditional transfer theorem.** Fix positive rational numbers `a,b,R`.
Suppose one fixed finite field has an effective, scale-one crossing family
with polynomially many pieces at every positive bit width `lambda`. Let
`r_i,c_i` be the two explicit factorization widths. For one fixed constant
`c`, assume

```math
0\le r_i,c_i\le 2\lambda,\qquad
a r_i+b c_i\le 2R\lambda+c\log_2(2\lambda),             \tag{1}
```

and assume all piece and factor queries can be computed in `2^(c lambda)`
bit operations. The pieces sum exactly to the crossing. Their field and
algorithms are fixed independently of the input.

If rational `0<alpha<1` and `0<theta<1` satisfy

```math
R<a\theta^\alpha+b(1-\theta)^\alpha,                  \tag{2}
```

then a deterministic finite-alphabet machine computes the smallest prime
factor of every `n`-bit integer in

```math
\exp\!\left(O\!\left(\frac{n}{(\log n)^\alpha}\right)\right).
                                                               \tag{3}
```

The implied constant is independent of the input. This is an ordinary
conditional theorem. No crossing family satisfying (2) with `alpha>1/2`
is supplied here, and the established factoring score remains `1/2`.

Optimizing over `theta` makes the right side of (2)

```math
D_\alpha(a,b)=
\left(a^{1/(1-\alpha)}+b^{1/(1-\alpha)}\right)^{1-\alpha}.
                                                               \tag{4}
```

Strict inequality permits rational choices of `alpha` and `theta`.
In particular, `R<sqrt(a^2+b^2)` implies (3) for some rational
`alpha>1/2`. The normalization `a^2+b^2=1` is convenient for searches
but is unnecessary for the theorem.

## 1. Exact rectangular networks and their capacities

The rectangular grid identity is the same finite sum-product identity as
in the square construction. On a `g` by `d` grid, put `k'=gd` and
`N'=(k')^2 lambda`. There are
`binom(g,2) binom(d,2)` crossing sites. For each assignment of a piece
to every site, its term has both a row-network and a column-network
presentation, with totals

```math
U_r=g(g-1)d\lambda+2\sum_q r_{i_q},\qquad
U_c=d(d-1)g\lambda+2\sum_q c_{i_q}.                  \tag{5}
```

The row network has `g` nodes, each storing at most `g d^2 lambda`
setting bits. The capacity between two nodes is at most

```math
d\lambda+\binom d2(2\lambda)=d^2\lambda.
```

Thus the grouping lemma applies with

```math
P_r=g d^2\lambda=N'/g.
```

Every node has at most `P_r` setting bits and total incident capacity,
and every edge has capacity at most `P_r/g`. The column statements
hold with `g,d` interchanged and `P_c=N'/d`. In particular,
`g P_r=d P_c=N'` exactly. None of these bounds depends on the chosen
piece assignment.

Using (1) in (5), and discarding a factor at most one in the crossing
count, gives

```math
aU_r+bU_c
\le\left(R+\frac a d+\frac b g\right)N'
  +\frac c2(k')^2\log_2(2\lambda).                    \tag{6}
```

The precise straight-edge terms are smaller than the displayed ones.
Linearity is the reason a bound on each piece controls an arbitrary
assignment of different pieces to the grid.

The pinned source definitions are
[`ClauseBound/Alg/Data.lean`, `capG` and `totalG`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Alg/Data.lean).
The proof uses their rectangular formulas, not a square-grid asymptotic
substitution. The scale-one identity and both factorizations use only
addition and multiplication, so these network presentations hold over
every commutative algebra over the fixed field.

## 2. Choose two child sizes and pad without adding free labels

Choose positive rational `beta_r,beta_c` strictly below `theta^alpha`
and `(1-theta)^alpha`, respectively, such that

```math
D'=a\beta_r+b\beta_c>R.                              \tag{7}
```

Choose a small positive rational `zeta<1` with
`R+2 zeta(a+b)<D'`. Leave a further fixed positive margin in this
inequality. All later choices use only these fixed constants.

For a parent clique with `k` vertices, set `B=ceil(A log k)` for a
sufficiently large fixed rational `A`, and define

```math
K_r=\left\lceil\frac{k^\theta}{B}\right\rceil,
\quad g=BK_r,
\quad K_c=\left\lceil\frac{k}{B^2K_r}\right\rceil,
\quad d=BK_c.                                       \tag{8}
```

Rational powers are computed by integer root comparisons. For example,
the least integer above `k^(u/v)/B` is found by comparing
`(B t)^v` with `k^u`. The other operations in (8) are integer
arithmetic. An integer-valued computable logarithmic upper approximation
may replace `ceil(A log k)` without changing any estimate below.

As `k` grows,

```math
\begin{aligned}
&gd=k(1+o(1)),\qquad g,d\longrightarrow\infty,\\
&\log K_r=\theta\log k-\log B+o(1),\\
&\log K_c=(1-\theta)\log k-\log B+o(1).
\end{aligned}                                                    \tag{9}
```

In particular, both child sizes tend to infinity, remain smaller than
`k`, and satisfy, for all sufficiently large `k`,

```math
(\log K_r)^\alpha\ge\beta_r(\log k)^\alpha,
\qquad
(\log K_c)^\alpha\ge\beta_c(\log k)^\alpha.           \tag{10}
```

Pad the input clique from `k` to `k'=gd` vertices. Each added vertex
has local tensor one exactly when every incident edge label is the zero
word, and zero otherwise. Original vertices use their original tables
and ignore coordinates for added vertices. Each original edge assignment
has exactly one contributing extension. Thus padding preserves the
value over every coefficient ring. Assigning constant-one tables to
dummy vertices would fail this property and is not the construction.

The algorithm answers a query to a padded table with at most one query
to the original table and polynomial additional work. It does not
materialize all padded tables before beginning the grid expansion.

## 3. Group each orientation with the exact case (ii) hypotheses

The grouping lemma's case (ii) requires a network of `h` nodes, a
destination of `K` groups, `K|h`, maximum capacity at most `P/h`, and

```math
2\zeta^2h\ge K\log(2hK).                            \tag{11}
```

It supplies a grouping whose group sizes are at most `2h/K` and whose
maximum between-group width is

```math
L=\left\lceil\frac{U+2\zeta hP}{K^2}\right\rceil.
                                                               \tag{12}
```

These are the hypotheses and conclusion proved in the pinned source's
[`GeneralPieces/Group/CaseII.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Group/CaseII.lean)
and [`GeneralPieces/Group/Lemma.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Group/Lemma.lean).
The existence proof concerns capacities only; enumerating all `K^h`
groupings and checking the displayed bounds makes the choice effective.

For the row network, `h=g=BK_r`; for the column network,
`h=d=BK_c`. The divisibility condition is exact. Section 1 gives the
capacity condition, and (11) becomes

```math
2\zeta^2B\ge\log(2hK).
```

Since `h,K<=2k` for all sufficiently large `k`, a fixed large `A`
makes this true for both orientations. Use positive child widths, taking
the maximum with one if necessary. Define

```math
W_r=U_r+2\zeta N'+K_r^2,\qquad
W_c=U_c+2\zeta N'+K_c^2.                             \tag{13}
```

Then the grouped child volumes satisfy `K_i^2 L_i<=W_i`.
The extra `K_i^2` covers the ceiling and the positive-width convention.
It is `o(N')` uniformly for `lambda>=1`.

Choose the orientation minimizing `W_i/beta_i`. This comparison uses
only rational arithmetic and integer widths; no numerical comparison
of real logarithms controls the algorithm. For nonnegative `W_r,W_c`,

```math
\min\left(\frac{W_r}{\beta_r},\frac{W_c}{\beta_c}\right)
\le\frac{aW_r+bW_c}{a\beta_r+b\beta_c}.              \tag{14}
```

Equations (6), (7), (9), and (13), with the reserved strict margin,
give fixed `eta<1`, `C_1`, and an effective threshold `k_0` such that
the selected child satisfies

```math
\frac{K_i^2L_i}{(\log K_i)^\alpha}
\le\eta\frac{k^2\lambda}{(\log k)^\alpha}
  +C_1 k^2\log_2(2\lambda),\qquad k\ge k_0.          \tag{15}
```

The threshold is uniform in `lambda>=1`. Indeed, every error outside
the logarithmic-width term is bounded by an explicit function of `k`
tending to zero: the padding ratio, `1/g`, `1/d`, and
`(K_r^2+K_c^2)/(k')^2`. Elementary bounds for logarithms by positive
powers make an upper threshold computable from the fixed constants.
The last term in (15) will enter the polynomial factor of the circuit
bound; it is not silently discarded.

## 4. Charge the complete recursive circuit

Let a table query count as one input operation. We prove, by strong
induction on `k`, that one uniformly generated arithmetic circuit has
size at most

```math
S(k,\lambda)=
2^{Ck^2\lambda/(\log(k+2))^\alpha}
\bigl(k(\lambda+1)+2\bigr)^{Dk^2},                 \tag{16}
```

for fixed sufficiently large `C,D`. Changing between `log k` and
`log(k+2)` only changes the margin and the finite threshold. Below
`k_0`, enumerate all edge-label assignments. There are only finitely
many such vertex counts, so choosing `C` large proves every base case
uniformly in `lambda`.

For larger `k`, enumerate all piece assignments. Their number is at
most `(2lambda)^(c(k')^2/4)`, which is a fixed power of
`k(lambda+1)+2` with exponent `O(k^2)`. For each assignment, compute
the two totals, select its orientation by (14), and search for a good
grouping. The grouping search examines at most `K_i^h` possibilities.
Its logarithm is `O(k log k)`, since `h<=2k`; this is
`o(k^2 lambda/(log k)^alpha)` uniformly in `lambda>=1`.

A grouped table query enumerates at most `2^(2hP/K_i)` settings.
Here `hP=N'`, so this is `2^(2N'/K_i)`. Its summands contain a
polynomial number of original table queries and piece-factor queries.
For two network nodes inside one group, their kernel can be evaluated
directly from the original piece values and straight-edge equalities.
Alternatively, summing its channel labels costs `2^O(N'/K_i)` and
fits the same bound.
Generating any factor costs at most `2^(c lambda)`. Thus replacing
one child input query by its defining circuit multiplies circuit size
by at most

```math
2^{O(N'/K_i+c\lambda)}
\bigl(k(\lambda+1)+2\bigr)^{O(1)}.                  \tag{17}
```

Both `K_r,K_c` grow as a fixed positive power of `k` divided by a
logarithm. The exponent in (17) is consequently
`o(k^2 lambda/(log k)^alpha)`, uniformly in `lambda`.
Together with grouping and index construction, it fits in the fixed
margin `1-eta` of (15), after increasing the finite base threshold.

By (15), the leading exponential term of the induction hypothesis
is at most the parent leading exponential term, with that margin,
times `2^(C C_1 k^2 log_2(2lambda))`. Choose `D` large enough to
absorb this latter factor and the number of piece assignments.
The child polynomial factor also fits: `L_i<=O(k^2 lambda)`, so

```math
\log\bigl(K_i(L_i+1)+2\bigr)
\le C_2\log\bigl(k(\lambda+1)+2\bigr),
\qquad K_i^2=o(k^2).                                \tag{18}
```

For sufficiently large `k`, its exponent is at most `Dk^2/2` after
the fixed factor `C_2` is included. Reserve the other half of `Dk^2`
for the nonrecursive polynomial factors just listed. This closes the
induction. The choices are noncircular: first fix the geometric margin
and a threshold for its exponential errors and for (18). The inequality
`C_2K_i^2<=k^2/2` is independent of `C,D`. Then choose `C>=1` for
the finite base cases, and finally choose `D` to cover the logarithmic
width error and the remaining polynomial factors.

The output is a circuit of additions and multiplications over the fixed
field algebra. The algorithm performs no division by input table
values. Piece coefficients and grouping choices are independent of those
values. Circuit generation, addresses, and loop counters have fixed
polynomial overhead in (16). This proves a uniform circuit bound over
every commutative algebra over the field, including rings with zero
divisors.

## 5. Compose with the multiplication constructor and sparse recovery

Choose `P=floor(n^(1/80))`, adjusting finitely many small inputs by
trial division. Then `P^40<=n`, `log P=Theta(log n)`, and the
existing multiplication constructor supplies exact weighted divisor
prefix instances with `P` vertices,

```math
P^2\lambda=O(n),\qquad
\text{construction cost}=2^{O(n/P)}.                \tag{19}
```

The older constructor already suffices. The
[operand-length construction](operand-length-transform-constructor.md)
also applies; summing its disjoint length classes adds only a polynomial
number of calls. Its smaller constants do not enter the exponent claim.

Equations (16) and (19) bound the weighted arithmetic circuit by

```math
\exp\!\left(O\!\left(
\frac{n}{(\log n)^\alpha}+P^2\log(n+2)+\frac nP
\right)\right)
=\exp\!\left(O\!\left(\frac{n}{(\log n)^\alpha}\right)\right).
                                                               \tag{20}
```

To obtain an ordinary factoring machine from a fixed-field circuit, use
the [sparse divisor nonzeroness proof](sparse-divisor-nonzeroness.md).
It evaluates the divisor polynomial in
`exp(O(n/log n))` cyclic polynomial rings, each of dimension
`exp(O(n/log n))`, and detects a nonempty divisor prefix in at least
one. The construction and the explicit bijection between valid network
settings and divisors are given in Sections 4--6 of
[scale-one-modular-transfer.md](scale-one-modular-transfer.md).

The same argument works for any fixed finite coefficient field:
its characteristic supplies the prime-field embedding, and each field
operation has constant bit cost. Alternatively, a fixed extension-field
crossing can be descended to its prime field with a fixed rank factor,
absorbed by the `O(log(2lambda))` allowance in (1).

Represent cyclic-ring elements by coefficient lists. Addition and
schoolbook multiplication have polynomial cost in the ring dimension.
Generating and simulating the entire circuit on a finite alphabet also
adds only a fixed polynomial factor. The combined exponent is

```math
O\!\left(\frac{n}{(\log n)^\alpha}+\frac n{\log n}\right).
```

Since `alpha<1`, the second term is absorbed by the first. A polynomial
number of prefix decisions returns the least divisor greater than one,
which is prime. This proves (3), with all ring and construction costs
included.

## 6. What remains a discovery problem

This theorem turns a strict supporting-line certificate below the unit
quarter circle into a sufficient route past `1/2`. It does not produce
that certificate. A numerical value must still be replaced by exact
rank bounds, a scale-one reconstruction, a common physical type, and
effective pieces at every label width. The new search code checks full
local profiles, but its floating pressure values are exploratory.

The rectangular grid and grouping lemmas are prior work of Anthropic;
the multiplication network is due to Isaac Rudich and Louis-Martin
Rousseau. This note extends the repository's ordinary fixed-field
transfer to unequal child sizes. No Lean formalization is added.
