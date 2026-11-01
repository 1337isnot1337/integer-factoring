# A sharp rank-sum budget with one distinguished piece

## Result

Suppose a tensor `E` has matrix rank `N` in each of two flattenings,
and

```math
E=B_1+\cdots+B_k
```

is a decomposition into nonzero tensors. Write `r_i,c_i` for their
two ranks. Assume that one specified summand has

```math
r_1,c_1\ge M,\qquad 1\le M\le N-2.
```

Then

```math
\boxed{\quad
\sum_i\sqrt{r_ic_i}\ \ge
\min\left\{
N,\quad \sqrt{M(N-1)}+\sqrt{N-M},\quad
M+2\sqrt{N-M-1}
\right\}.
\quad}                                                       \tag{1}
```

This is the **exact optimum of the numerical rank-budget relaxation**:
positive ranks in each cut, total ranks at least `N` in each cut, and
the distinguished lower bound `M`. It need not be attainable by an
actual tensor decomposition. No simultaneous rank realization is
asserted.

For the full four-point permutation seed, the cycle-block theorem gives
`N=24^2=576` and `M=72`. Thus every decomposition in that coefficient
span satisfies

```math
\sum_i\sqrt{r_ic_i}\ge72+2\sqrt{503}>116.                 \tag{2}
```

The previous finite-threshold proof only needed the weaker lower bound
96. Equation (2) does not establish the experimentally observed optimum
325 in the full span, and it does not exclude new collective
recombinations of tensor powers of a fixed seed.

## 1. Reduce to exact rank budgets

Rank subadditivity gives

```math
\sum_i r_i\ge N,\qquad\sum_i c_i\ge N.
```

Every nonzero tensor has positive rank in both cuts. If
`M+k-1>N`, positivity gives

```math
\sum_i\sqrt{r_ic_i}\ge M+k-1>N,
```

so (1) already holds. Otherwise `k<=N-M+1`. Reduce coordinates,
without crossing their respective lower bounds, until each total is
exactly `N`. The objective is increasing in every coordinate, so this
can only lower it. The same reduction is possible with integer ranks.

For fixed `k`, it is therefore enough to minimize over the two
simplices

```math
\begin{split}
r_1\ge M,\quad r_i\ge1\ (i>1),\quad\sum_i r_i=N,\\
c_1\ge M,\quad c_i\ge1\ (i>1),\quad\sum_i c_i=N.
\end{split}                                                   \tag{3}
```

The following argument also permits real coordinates, so it proves
the desired lower bound for integer ranks automatically.

## 2. A minimum occurs at two simplex vertices

For a fixed column vector `c`, the objective
`sum_i sqrt(c_i) sqrt(r_i)` is concave in `r`. Express any `r` as a
convex combination of vertices of its simplex. Concavity implies that
some vertex has objective no larger than that of `r`. Replace `r` by
such a vertex. With that `r` fixed, apply the same argument to `c`.

Thus there is a minimizer with both vectors at simplex vertices. Put

```math
R=N-M-k+1\ge0.
```

The base vector is `(M,1,...,1)`. Each vertex adds the entire residual
`R` to one coordinate. There are only the following cases.

| Locations receiving the two residuals | Objective |
| --- | --- |
| Both at the distinguished coordinate | `N` |
| Both at the same ordinary coordinate | `N` |
| One distinguished, one ordinary | `sqrt(M(N-k+1))+sqrt(N-M-k+2)+k-2` |
| Two different ordinary coordinates | `M+2sqrt(N-M-k+2)+k-3` |

The third case requires `k>=2`, and the fourth requires `k>=3`.
Every listed vector is integral, so the real and integer relaxations
have the same minimum.

## 3. Minimize over the number of summands

For the third row of the table, define

```math
g(k)=\sqrt{M(N-k+1)}+\sqrt{N-M-k+2}+k-2.
```

This is concave on `2<=k<=N-M+1`; each square root is concave and
its argument is affine. Its minimum on the interval is at an endpoint.
The two endpoint values are

```math
g(2)=\sqrt{M(N-1)}+\sqrt{N-M},\qquad
g(N-M+1)=N.                                             \tag{4}
```

Likewise

```math
h(k)=M+2\sqrt{N-M-k+2}+k-3
```

is concave for `3<=k<=N-M+1`. This interval is nonempty because
`M<=N-2`. Its endpoint values are

```math
h(3)=M+2\sqrt{N-M-1},\qquad h(N-M+1)=N.                 \tag{5}
```

Equations (4)--(5), and the aligned value `N`, prove (1).

All three candidates can be attained in the numerical relaxation:

- `k=1`, with ranks `(N,N)`.
- `k=2`, with row vector `(N-1,1)` and column vector `(M,N-M)`.
- `k=3`, with row vector `(M,N-M-1,1)` and column vector
  `(M,1,N-M-1)`.

This proves the claimed sharpness of the relaxation. When `M=N-1`
or `M=N`, a separate direct check gives the relaxed optimum `N`.
Those endpoints are excluded from (1) only to avoid an empty third
case in the displayed formula.

## 4. Application to permutation coefficient spans

The [cycle obstruction](permutation-span-cycle-obstruction.md) proves
that a decomposition of the full `S_b` identity-kernel crossing has
some summand with

```math
r,c\ge M_b=\frac{b!(b-1)!}{2},\qquad N=(b!)^2.
```

For `b=4`, the third expression in (1) is the smallest. Indeed,

```math
72+2\sqrt{503}<118<203<\sqrt{72\cdot575}
\le\sqrt{72\cdot575}+\sqrt{504},
```

and `118<576`. The strict lower bound in (2) follows from `503>22^2`.
For `b=3`, the corresponding bound is `6+2sqrt(29)`.

Tensoring the block argument alone gives the weaker distinguished
threshold `M_b^t` at dimension `N^t`. Applying (1) there does not
automatically preserve a single-seed obstruction against the
`(b!)^(sqrt(2)t)` target. In particular, for `b=4`, `72<24^sqrt(2)`.
This numerical lemma makes no new exclusion of collective tensor-power
recombinations.

## Verification

[check_distinguished_rank_budget.py](../tools/check_distinguished_rank_budget.py)
enumerates all positive integer rank vectors with exact budgets for
small `N`. It compares their minima with (1) using integer enclosures
of square roots, including the two endpoint cases `M=N-1,N`.
It also verifies the 24-label exact inequalities.

```sh
python3 tools/check_distinguished_rank_budget.py \
  --output evidence/distinguished-rank-budget.json
```

The finite checker concerns only the numerical relaxation. It does
not search for or assert tensor decompositions realizing those rank
vectors.
