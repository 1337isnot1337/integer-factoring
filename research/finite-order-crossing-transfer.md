# Exact crossings from finite-order algebraic approximations

## Result

A finite-order approximation to a crossing can supply exact crossing
pieces. The conversion costs only a factor `kappa+1` in each matrix
rank, where `kappa` is the coefficient order that contains the crossing.

Let `K` be the rationals or one fixed finite field. For each binary
label width `lambda`, suppose tensors `B_i(z)` have entries in `K[[z]]`
and satisfy

```math
\sum_i B_i(z)=z^\kappa E+O(z^{\kappa+1}),\qquad
E(x,y,u,v)=[x=v][y=u].
\tag{1}
```

Suppose the two flattenings of `B_i(z)` have ranks `r_i,c_i` over the
formal Laurent-series field `K((z))`. Define exact tensors

```math
A_i=[z^\kappa]B_i(z).
```

Then

```math
\sum_i A_i=E,\qquad
\operatorname{rank}_{\rm row}A_i\le(\kappa+1)r_i,\qquad
\operatorname{rank}_{\rm col}A_i\le(\kappa+1)c_i.
\tag{2}
```

This is useful for algebraic degenerations whose individual summands
have poles. Multiplication by a common power of `z` makes them regular
and puts the crossing at a positive coefficient order, as in (1).
The common pole order and the finite coefficient-generation costs must
be controlled. An unspecified limiting decomposition is insufficient.

## 1. Coefficient rank over a discrete valuation ring

Let `B(z)` be any rectangular matrix over `K[[z]]`, of rank `r` over
`K((z))`. We first show that

```math
B(z)=\sum_{j=1}^r L_j(z)R_j(z),
\qquad L_j\in K[[z]]^{a\times1},\quad
R_j\in K[[z]]^{1\times b}.
\tag{3}
```

If `B=0`, the claim is immediate. Otherwise choose a nonzero entry
of minimum `z`-valuation. Permute its row and column to position `(1,1)`.
Every ratio `B_i1/B_11` lies in `K[[z]]`, by minimality of the valuation.
Take `L_i=B_i1/B_11` and `R_j=B_1j`.
The residual matrix `B-LR` is regular and has zero first row and column.
Ordinary elimination over `K((z))` lowers its rank by exactly one.
Induction proves (3).

For one regular outer product,

```math
[z^\kappa](L_jR_j)
=\sum_{a=0}^{\kappa}([z^a]L_j)([z^{\kappa-a}]R_j).
```

Each summand has rank at most one. Summing over `j` proves

```math
\operatorname{rank}_K[z^\kappa]B\le(\kappa+1)r.
\tag{4}
```

Apply (4) independently to the two flattenings to obtain (2).
No generic-point evaluation, numerical limit, or division by a zero
element of a truncated polynomial ring is used.

The factor `kappa+1` is necessary in general. Set
`L_i=z^i` and `R_j=z^(kappa-j)` for `0<=i,j<=kappa`.
The matrix `LR` has generic rank one, while its coefficient of
`z^kappa` is the identity matrix of order `kappa+1`.

## 2. Effective transfer and the allowed approximation order

Assume the piece set and each required coefficient in (1) are computable
in `2^O(lambda)` bit operations, including complete index enumeration.
Assume its piece entropy obeys

```math
\log_2\max(2,|I_\lambda|)
=O\!\left(\frac{\lambda}{(\log\lambda)^\delta}\right)
```

for a fixed `delta>0`. Suppose also

```math
\begin{aligned}
\log_2\max(1,r_i)+\log_2\max(1,c_i)
&\le4\rho\lambda+
 O\!\left(\frac{\lambda}{\log\log\lambda}\right),\\
\log_2(\kappa+1)
&=O\!\left(\frac{\lambda}{\log\log\lambda}\right),
\qquad \tfrac12\le\rho<1.
\end{aligned}
\tag{5}
```

The exact coefficient tensors then satisfy the width requirement in
[subexponential-piece-transfer.md](subexponential-piece-transfer.md).
Indeed, (2) adds at most `2log2(kappa+1)` to the width sum.
Ceilings add at most two further bits. Each separate width can be capped
at `2lambda` by the full-input factorization.

The complete matrices have order `2^(2lambda)` in each flattening.
Their entries and exact Gaussian factorizations are computable in
`2^O(lambda)` bit operations over a fixed finite field.
Over the rationals, the
[rational-height theorem](rational-height-transfer.md) gives the same
conclusion, even for exponentially long coefficient numerators and
denominators. Thus no effective construction of the infinite series
factors in (3) is needed by the factoring machine.

Consequently (1) and (5) give factoring exponent
`min(log2(1/rho),delta,1)` over a fixed finite field, and
`min(log2(1/rho),delta)` over the rationals.
The fixed-field result includes its sparse-divisor recovery cost.
Polynomial approximation orders are well within the allowed range.
Orders `2^(O(lambda^beta))`, for fixed `beta<1`, also satisfy (5).

This is a conditional interface theorem. It supplies no approximation
with `rho<1/sqrt(2)` and proves no factoring exponent above `1/2`.

## 3. Regularity, poles, and finite input data

If tensors `T_i(z)` lie in `K((z))`, a known common pole bound `b`
makes `B_i=z^b T_i` regular. When
`sum_i T_i=E+O(z)`, equation (1) holds with `kappa=b`.
The order bound in (5) applies to this common pole bound.

Only the coefficient of order `kappa` is needed for the exact pieces.
A routine that outputs the complete finite jet through that order is
also sufficient, if its total bit cost meets the stated bound.
The generic ranks in the hypothesis refer to the original formal-series
matrices. They cannot be certified by a low rank at one parameter value.

For example, the one-entry matrix `(z)` has rank zero at `z=0` but
generic rank one. Its coefficient of order one is nonzero.
Also, truncating entries can change generic rank: the rank-one matrix

```math
\begin{pmatrix}1&z\\z&z^2\end{pmatrix}
```

has a first-order entry truncation whose determinant is `-z^2`.
Generic rank and the rank of a coefficient are different quantities.
Equation (4) is the precise relation used here.

## 4. Reproduction

Run:

```sh
python3 tools/check_finite_order_crossing_transfer.py
```

The checker uses exact polynomial coefficients and matrix ranks over
several fields. It checks the coefficient-rank inequality, sharp examples,
a crossing identity with cancellation at lower order, and the parameter
evaluation and truncation controls. These cases supplement the proof.

The eventual evaluator and cost transfers use the pinned Anthropic and
factoring constructions credited in the linked notes. This coefficient-rank
argument and its application used substantial OpenAI assistant support.
No new Lean formalization was run.
