# Sum-of-squares identities and the two crossing ranks

## Result and scope

A low sum-of-squares or partial-transpose rank does not by itself give
the two simultaneous flattening bounds needed by the factoring evaluator.
This note makes that distinction quantitative and records an explicit
many-mask example with a complete rank profile.

- A direct Gram/partial-transpose pair always has rank product at least
  `m(m²-m)` on an `m`-label crossing. It requires `rho>=3/4` asymptotically.
  The same restriction holds for a finite-order degeneration of this pair.
- A determinant multiplication construction has individual orbit terms
  approaching `rho=1/2`, but exponentially many such terms. Grouping them
  into polynomially many transition coefficients has the exact rate
  `rho=log2(12)/4=0.896240625...`.

The latter is a uniform effective family over `F_5`; its associated
factoring exponent is `log2(4/log2(12))=0.158041971...`, below the already
established `1/2`. No improved factoring rate is claimed.
These restrictions concern the specified conversions. They do not rule
out all uses of sums of squares, several unrelated partial transposes,
or other decompositions of the crossing.

The literature motivation is Hrubeš,
[*A Subquadratic Upper Bound on Sum-Of-Squares Composition Formulas*](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CCC.2024.12),
CCC 2024, and Rossman--Zhu,
[*Multiquadratic Sum-of-Squares Lower Bounds Imply VNC1 != VNP*](https://arxiv.org/abs/2512.01227),
Sections 3 and 6. Their partial-transpose complexity permits a separate
transpose mask for each summand. Our interface additionally needs both
specified crossing flattenings of every retained piece to be small.
The proofs below use their algebraic setting but derive the required
crossing profiles directly.

## 1. Exact conventions and the direct two-piece conversion

The crossing tensor is `E(x,y,u,v)=[x=v][y=u]`.
Its row matrix, with rows `(x,y)` and columns `(v,u)`, is `I_(m²)`.
For a matrix `A` in this indexing, define

```math
(\Gamma A)_{(x,y),(v,u)}=A_{(x,u),(v,y)}.
\tag{1}
```

This is exactly the other crossing flattening, up to the names of its
indices. Thus a tensor represented by `A` has rank pair
`(rank A, rank Gamma A)`.

Suppose a bilinear sum-of-squares identity over a field of characteristic
different from two is

```math
\left(\sum_x X_x^2\right)\left(\sum_yY_y^2\right)
=\sum_{a=1}^{s}\left(\sum_{x,y}L_{a,x,y}X_xY_y\right)^2.
```

Let `A=L^T L`. Comparing the polarized coefficients gives

```math
A+\Gamma A=2I_{m^2},\qquad \operatorname{rank}A\le s.
\tag{2}
```

The two crossing pieces are `A/2` and `(Gamma A)/2`. Their rank pairs
are `(r,r')` and `(r',r)`. If `s<=m^(c+o(1))`, the immediate guarantee
is `rho<=(2+c)/4`, since the other rank is at most `m²`.
Hrubeš's construction gives a fixed exponent below `1.616`, not a
coefficient below `1/sqrt(2)` under this conversion.
For a fixed odd finite field its coefficient algorithm is polynomial
in `m`: it enumerates subsets of `O(log m)` coordinates, constructs the
explicit Clifford matrices, and solves the `O(log m)` binomial system
from the paper's Appendix A. All field coefficients have fixed size.
Consequently arithmetic effectivity is not the missing hypothesis here.

There is a stronger restriction independent of the particular SoS formula.

**Paired partial-transpose theorem.** If a matrix over a field of
characteristic different from two satisfies (2), then, for `m>=2`,

```math
r r'\ge m(m^2-m).
\tag{3}
```

Indeed, fixing the second row and column index to the same value makes
`Gamma` act trivially. Equation (2) makes that `m`-by-`m` block of both
matrices equal to `I_m`. Hence `r,r'>=m`. Rank subadditivity also gives
`r+r'>=m²`. For two numbers with these constraints and upper bound `m²`,
their product is at least `m(m²-m)`. Taking logarithms proves the
asymptotic requirement `rho>=3/4` for these two pieces.

This argument also covers finite-order degenerations. Suppose

```math
A(t)+\Gamma A(t)=2t^\kappa I_{m^2}+O(t^{\kappa+1}).
\tag{4}
```

Work over the Laurent-series field. Each fixed-index diagonal block has
leading term `t^kappa I_m` and is invertible there. The sum in (4) is
also invertible there. Therefore the same two rank inequalities, and
(3), hold for the generic ranks. Allowing a finite coefficient jet does
not remove this particular obstruction.

## 2. A many-mask construction from determinant multiplication

Here is a different, explicit use of many partial transposes. It avoids
assuming that their individual low ranks survive aggregation.

For each of `d` positions let `X_j` be a two-by-two matrix of variables.
The elementary identity

```math
\det(X_1\cdots X_d)=\prod_{j=1}^{d}\det X_j
\tag{5}
```

provides a multiquadratic composition. Write physical indices at position
`j` as pairs `(a_j,b_j)` in `{0,1}²`, and let
`epsilon=[[0,1],[-1,0]]`. The determinant quadratic form has matrix

```math
Q_{(a,b),(c,e)}=\tfrac12\epsilon_{a,c}\epsilon_{b,e}.
```

The coefficient matrix for the left side of (5) is

```math
D_{(a_1,b_1),\ldots,(a_d,b_d);,(c_1,e_1),\ldots,(c_d,e_d)}
=\tfrac12\epsilon_{a_1,c_1}\epsilon_{b_d,e_d}
 \prod_{j=1}^{d-1}[b_j=a_{j+1}][e_j=c_{j+1}].
\tag{6}
```

For `eta in {0,1}^d`, let `D_eta` partially transpose the positions
where `eta_j=1`. The matrix `D` is symmetric, so complementary masks
give the same matrix. Choose the representatives with `eta_d=0`.
Polarizing (5) separately in its `d` variable blocks proves

```math
\sum_{\eta:\eta_d=0}D_\eta=2^{d-1}Q^{\otimes d}.
\tag{7}
```

For an ordinary crossing with identity metric, work over `F_5`, with
`i=2`, so `i²=-1`. The substitution

```math
X=\begin{pmatrix}
x_0+i x_1&x_2+i x_3\\
-x_2+i x_3&x_0-i x_1
\end{pmatrix}
```

has determinant `x_0²+x_1²+x_2²+x_3²`. It is an invertible linear
change of coordinates, applied equally at both ends of every position.
It changes the metric in (7) to the identity, preserves every relevant
flattening rank, and commutes with the partial-transpose operations.

Take even `d` and let `J` be the alternating positions. Group the
positions in `J` into one crossing label and the remaining positions
into the other. Each label has dimension `4^(d/2)=2^d`; hence its
binary length is `lambda=d`. Its column flattening is `Gamma_J`.

## 3. Exact factorization of the transition polynomial

Set `N=d-1`. The transition set of a representative `eta` consists
of the positions `j` with `eta_j != eta_(j+1)`. Every subset of
`{1,...,N}` is the transition set of exactly one representative.
Let

```math
A_t=\sum_{\eta:\eta_d=0,\,\#\mathrm{transitions}(\eta)=t}D_\eta,
\qquad A(z)=\sum_{t=0}^{N}z^t A_t.
\tag{8}
```

There are only `d` coefficient pieces. Alternating partial transpose
complements all transition sets, so

```math
\Gamma_J A_t=A_{N-t}.
\tag{9}
```

The following factorization proves their ranks exactly. Regroup the
two-bit row coordinates as the end pair `(a_1,b_d)` and the `N`
internal pairs `(b_j,a_(j+1))`, and do the same for columns.
At a transition, the local matrix in (6) is the four-by-four swap `S`.
At a nontransition it is `F=vv^T`, where `v=(1,0,0,1)^T`.
The end factor is the invertible matrix `(epsilon tensor epsilon)/2`.
The end sign is `(-1)^(eta_1)=(-1)^t` because `eta_d=0`.

Let `D_*` be the matrix for the all-transition mask; it is invertible.
Since `FS=F`, right multiplication by `D_*^(-1)` gives

```math
A(z)D_*^{-1}
=I_4\otimes\bigotimes_{j=1}^{N}(R+zI_4),
\qquad R=-vv^T,\quad R^2=-2R.
\tag{10}
```

This equality includes the sign above: every missing transition
contributes one factor `-F`. It is a direct tensor factorization,
without a simultaneous-diagonalization assumption about the original
four-index pieces.

Over every field of characteristic different from two, `R` has a
one-dimensional `-2` eigenspace and a three-dimensional zero eigenspace.
On the tensor subspace with exactly `k` zero eigenspaces, the polynomial
in (10) acts by `z^k(z-2)^(N-k)`. Therefore

```math
\boxed{\quad
r_t:=\operatorname{rank}A_t
=4\sum_{k=0}^{t}\binom Nk3^k
 \left[\binom{N-k}{t-k}\ne0\text{ in the field}\right].\quad}
\tag{11}
```

The column rank is `r_(N-t)`. In characteristic zero the indicator is
always one. In a fixed odd characteristic some summands disappear,
but the summand `k=t` never does.

## 4. The complete asymptotic profile and its cost

The surviving boundary terms in (11), at a middle value of `t`, give

```math
\max_t r_t r_{N-t}
\ge16\binom N{\lfloor N/2\rfloor}^{\!2}3^N
\ge\frac{16\,12^N}{(N+1)^2}.
\tag{12}
```

For the matching upper bound, omit the indicator and use Vandermonde:

```math
\begin{aligned}
r_t r_{N-t}
&\le16\sum_{i\le t,\,j\le N-t}\binom Ni\binom Nj3^{i+j}\\
&\le16\sum_{k=0}^{N}\binom{2N}{k}3^k
\le16(N+1)12^N.
\end{aligned}
\tag{13}
```

Thus the exact leading total width is

```math
\max_t\bigl(\lceil\log_2r_t\rceil+
                 \lceil\log_2r_{N-t}\rceil\bigr)
=\lambda\log_2 12+O(\log(\lambda+2)).
\tag{14}
```

Divide the pieces by `2^N` and apply the coordinate change of Section 2.
Their sum is the ordinary scale-one crossing. Equation (14) gives
`rho=log2(12)/4`, with each separate width at most `2lambda`.
An odd label length can use the preceding even length together with
one fixed one-bit crossing; its constant width cost is absorbed in (14).

The individual orbit terms explain why grouping mattered. A term with
`t` transitions has ranks `4^(t+1)` and `4^(N-t+1)` by (6).
Their product is `4^(N+2)`, which approaches the ideal leading parameter
`rho=1/2`. There are `2^N` such terms. In the polynomial family (8),
their exact aggregated rank product grows as `12^N` up to polynomial
factors, as proved above. No extrapolation from the individual terms
can omit that cost.

All generation requirements are effective. Formula (6) computes any
orbit entry in polynomial time in `lambda`. Summing the `2^N` masks
computes any coefficient entry in `2^O(lambda)` time; the local coordinate
change has the same allowed cost. Equation (11) computes the widths
using integer binomials and modular nonvanishing tests in polynomial
time. Exact finite-field elimination on the `2^(2lambda)`-square
flattenings computes complete row and column factors in
`2^O(lambda)` time. Complete enumeration, all factor tables, and index
metadata fit the same bound. There are `lambda+O(1)` pieces and a
nonzero, computable scale. The fixed-field transfer therefore applies,
with the weaker exponent stated at the start.

## 5. Checks and what remains open

Run:

```sh
python3 tools/check_sos_partial_transpose.py
```

The checker constructs (6) directly, compares the regrouped tensor
factorization, verifies (7), checks the two flattenings, and compares
literal finite-field ranks with (11). Separate integer checks verify
the bounds (12)--(13). Its output is
[`sos-partial-transpose.json`](../evidence/sos-partial-transpose.json).

These proofs close the direct paired conversion and this transition
coefficient family, including their stated arithmetic requirements.
They leave open regroupings outside (8), constructions with several
unrelated transpose masks, and genuinely different rank profiles.
The rank restrictions are not lower bounds on factoring.

The literature results retain their attribution. The translation,
explicit profile analysis, and checker used substantial OpenAI
assistant support. No new Lean formalization is claimed.
