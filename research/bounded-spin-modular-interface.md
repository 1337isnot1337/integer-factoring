# The remaining factoring obligation for the coherent spin families

## Conditional conclusion

The repository already proves a
[scale-one fixed-characteristic transfer](scale-one-modular-transfer.md),
using [sparse divisor nonzeroness](sparse-divisor-nonzeroness.md).
The spin constructions can use that theorem directly over `F_5`.
Lifting their ranks to characteristic zero is not required by this route.

For either the [original pencil](bounded-spin-crossing.md) or the
[enlarged six-vertex family](bounded-spin-six-vertex.md), suppose its
explicit coefficient pieces, or an effectively specified recombination,
satisfy at every binary width `lambda`

```math
\sum_i\sqrt{r(A_i)c(A_i)}
 \le \operatorname{poly}(\lambda)
       2^{(\sqrt2-\epsilon)\lambda},\qquad \epsilon>0.       \tag{1}
```

Here the pieces must sum to the crossing and the recombination algorithm,
if needed, must cost `2^O(lambda)`. Then the existing ordinary transfer
gives a deterministic smallest-prime-factor machine with saving exponent
strictly greater than `1/2`. For
`0<epsilon<=sqrt(2)-1`, its parameter is

```math
\rho=\frac{\sqrt2-\epsilon}{2},\qquad
\alpha=\log_2\frac2{\sqrt2-\epsilon}.                       \tag{2}
```

A stronger bound can instead use `rho=1/2` and the existing score-one
endpoint. These are conditional conclusions. Neither spin family is
currently known to satisfy (1).

Dense rank factorization already supplies the permitted factor algorithms.
For the original pencil, even selection of coefficient recombinations can
be made automatic by an exponential search of the permitted size. The
larger six-vertex space has a different selection-cost obligation, detailed
below. The purpose here is to separate the actual remaining rank problem
from effects already handled by the transfer theorem.

## 1. Verify the effective-family interface

Fix `lambda` and use the first `2^lambda` return paths at the least even
length `n>=4` with sufficiently many paths. The constructor proves
`n=O(lambda)`, polynomial-time path unranking, and paired encoder norms
nonzero modulo five. The endpoint is exactly the crossing, with scale one.

The original family has `O(lambda)` coefficient tensors. The enlarged
family has `O(lambda^2)`. Entries of either family are uniformly computable
in polynomial time in `lambda`, including coefficient extraction. A
specified recombination with `2^O(lambda)` generation and entry cost
has the same exponential bound. The number of retained pieces must be
polynomial in `lambda`.

For a piece `A`, each cut matrix has side `2^(2lambda)`. Write all its
entries and perform ordinary exact Gaussian elimination over the fixed
field. Naive elimination and matrix storage take at most `2^O(lambda)`
bit time. For rank `R`, select `R` independent original columns and
express every column in their basis. This gives a factorization

```math
A=UV,\qquad U\in F_5^{2^{2\lambda}\times R},\quad
V\in F_5^{R\times2^{2\lambda}}.
```

Pad its channel with zero entries to `2^ceil(log2 R)`, and repeat for the
other cut. Rank-zero pieces can be dropped. A factor-value query may
reconstruct this whole elimination each time; its cost remains
`2^O(lambda)`, which is precisely the allowed interface bound. A faster
factor algorithm may be useful but is not a missing theorem here.

The two channel bit widths obey

```math
r_i+c_i\le\log_2(r(A_i)c(A_i))+2\le4\lambda+2,
```

and individually they are at most `2lambda` because each cut matrix
has only `2^(2lambda)` columns. From (1), each individual product satisfies

```math
r_i+c_i\le(2\sqrt2-2\epsilon)\lambda+O(\log\lambda)
         =4\rho\lambda+O(\log\lambda).
```

Thus all four hypotheses of section 1 of the existing modular transfer
are met. Its arbitrary fixed constant `c` absorbs dense elimination,
piece count, indexing, padding and generation costs. The transfer, rather
than this note, supplies the weighted multiplication network, arbitrary
commutative-ring evaluator, sparse nonzeroness family and finite-alphabet
machine simulation.

The integer companions are optional for this argument. Their scale
`12^n` remains a unit modulo five and could also be normalized, but a
modular rank saving is never treated as a rational rank saving.

## 2. Original coefficient-space selection is automatically effective

Let `V_lambda` be the coefficient span of the original pencil after
restriction to binary labels. Its dimension is `O(lambda)`. Fix a
computable width budget

```math
b_\lambda=4\rho\lambda+O(\log\lambda)
```

with a fixed `rho<1/sqrt(2)`, rounded down to an integer. Suppose the
crossing belongs to the span of those tensors in `V_lambda` whose
two channel widths sum to at most `b_lambda`.

This purely algebraic containment has an automatic uniform algorithm:

1. Enumerate all `5^(n+1)=2^O(lambda)` linear combinations of the
   original coefficient tensors.
2. Write both cut matrices, compute their ranks, and retain the pieces
   meeting the width budget.
3. Perform linear elimination on their flattened tensors, retaining a
   basis of their span and solving for the crossing in that basis.
4. Scale each selected basis piece by its solution coefficient, discard
   zeros, and compute the two rank factorizations as in section 1.

At most `n+1=O(lambda)` pieces survive. Nonzero scalar multiplication
does not change rank. Each enumeration step and every matrix operation
costs `2^O(lambda)`; the total has that same form, with a fixed constant.
The algorithm need not guess a promising coefficient basis or be given
the low-rank directions as advice.

In particular, an existential decomposition satisfying (1) over `F_5`
is enough for this original space. Every piece has the individual
rank bound needed for the containment above; replacing its collection
by a spanning basis costs only an additional polynomial factor. This
statement does not establish that the containment holds. The exact
eight-site certificate proves a finite failure for the full thirteen-path
alphabet at its stated threshold.

## 3. What changes for the enlarged space

The six-vertex constructor is already effective, and the rank factors
of any explicit selected pieces are still automatic. However, its
formal coefficient count is `O(lambda^2)`. Enumerating all combinations
by that bound costs `5^O(lambda^2)`, which exceeds the allowed
`2^O(lambda)` family-generation cost.

Therefore a theorem for its raw univariate coefficients, or a uniform
`2^O(lambda)` recipe for favorable recombinations, would close the
interface. A bare existential cheap-span assertion does not yet do so
unless one separately improves the dimension bound or the selection
algorithm. Finite optimization at one fixed length does not resolve
this uniform cost issue.

## 4. Bounds on whole path alphabets also suffice

One may prove rank estimates first for all `L_n=F_(n-1)` paths, rather
than a power-of-two selection. For the least even `n` chosen above,

```math
2^\lambda\le L_n<3\cdot2^\lambda
```

apart from an immaterial finite initial case, since consecutive allowed
path counts have ratio at most three. Restricting all four legs to the
selected paths preserves the crossing identity and cannot increase
either cut rank. Replacing `L_n` by `2^lambda` in a uniform power bound
therefore changes only a constant factor. For the original pencil,
the search in section 2 then finds the required restricted pieces.

No new Lean formalization or improved exponent is asserted here. The
conditional transfer and sparse recovery are prior repository work;
this note supplies the spin-family interface and its selection-cost
distinction. The derivation used substantial Codex AI assistance.
