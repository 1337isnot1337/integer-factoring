# Near-all-ones labels and the partition-defect grading

The existing crossing identity remains valid after restricting its labels
to partitions near the all-ones partition. This is a possible change of
label family, so it deserves a separate support and rank analysis.

The natural full fixed-defect family does not improve the leading
contraction by lossless compression of the existing factor vectors.
Its creation terms already retain the square-root partition growth.
The defect cutoff also fails to be invariant under the actual operators.
More structured label subsets remain a separate question.

This note uses the Appendix B construction in
`anthropics/formal-math` at
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The operator recursion is in
[`ClauseBound/Boson/Defs.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Boson/Defs.lean).
The dual families and crossing pieces are in
[`AppHalfPiecesStatement.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
All new claims below are ordinary proofs or explicitly finite checks.

## 1. Encode labels by their defect

For a partition `mu`, define

```math
\operatorname{def}(\mu)
=\operatorname{wt}(\mu)-\operatorname{length}(\mu)
=\sum_j(\mu_j-1).
```

Fix an energy `H` and a defect `D` with `H>=2D`.
Every partition `nu` of `D` gives a partition of `H`:

```math
E_H(\nu)
=(\nu_1+1,\ldots,\nu_\ell+1,
  1^{\,H-D-\ell}).
\tag{1}
```

This is a bijection onto the weight-`H` partitions of defect `D`.
Its inverse subtracts one from every part greater than one and removes
the ones. Thus there are exactly `p(D)` such labels.
The assumption `H>=2D` guarantees nonnegative padding for every
partition of `D`.

When `H/D` grows, the fraction of parts equal to one tends to one
uniformly over this family. Yet the information content is determined
by `D`: `log2 p(D)=kappa sqrt(D)+O(log D)`.
Increasing the common all-ones padding does not create more labels.

This is a different label map from the source's canonical prefix at
weight `H`. The algebraic crossing identity remains valid on the
selected labels. Any algorithm using this change must separately bound
label generation, coefficient sizes, and evaluation costs.

## 2. The actual operator increases defect

Let `xi,eta` be partitions of `H` with defects `d_x,d_y`.
Consider a nonzero normal-order term in `(t^xi)[r]t^eta`.
Suppose it contains `A` annihilations.
Each annihilation removes one operator part and one input part.
The output has length

```math
2H-d_x-d_y-2A.
```

Its weight is `H+r`. Therefore its defect is exactly

```math
\operatorname{def}(\text{output})
=r-H+d_x+d_y+2A.
\tag{2}
```

In the highest channel `r=H`, this becomes `d_x+d_y+2A`.
Only the creation-only terms preserve the additive defect bound.
Annihilation can increase defect by an amount comparable to `H`.

There is an explicit example even when both input defects vanish:
for `H>=2`,

```math
[t_1t_{2H-1}]\,(t_1^H)[H]t_1^H=H\cdot H!\ne0.
\tag{3}
```

Indeed, the output has two parts, so exactly `H-1` factors annihilate.
There are `H` choices for the one creating factor. The annihilations
contribute the falling factorial `H!`, and the creation coefficient is one.
The output defect is `2H-2`.

Thus restricting all intermediate coordinates to defect at most
`d_x+d_y` destroys an exact coefficient. This is a support statement.
A vector with large support can still belong to a small-dimensional
family, so (3) alone is not a rank obstruction.

The inhomogeneous source vector `a^xi` does not repair this cutoff.
Every monomial in `a^xi` still has `length(xi)` parts.
In the highest row channel, the angle operator selects exactly its
weight-`H` part, namely `(-1)^(d_x)t^xi`.

## 3. The natural fixed-defect family still has a large factor rank

The embedding (1) preserves multiset addition:

```math
E_H(\nu)\mathbin{\cup}E_H(\omega)
=E_{2H}(\nu\mathbin{\cup}\omega).
\tag{4}
```

Let `D=D(lambda)` be the least integer at least two with
`p(D)>=2^lambda`. Choose the first `2^lambda` partitions of `D`
in descending-list lexicographic order. Embed them using (1), with
any `H>=2D`. Denote the resulting label set by `S`.

Consider the highest operator map

```math
U_{H,D}:\mathbb Q^{S\times S}\longrightarrow V_{2H},
\qquad
e_{(\xi,\eta)}\longmapsto(t^\xi)[H]t^\eta.
```

Every input pair has the same total number of parts, `2H-2D`.
Its creation-only term is `t^(xi union eta)` with coefficient one.
All annihilation terms have fewer parts. Different union monomials
therefore give linearly independent images.

The bounded-part construction in
[`binary-prefix-top-channel.md`](binary-prefix-top-channel.md)
counts unions of the actual prefix labels at defect weight `D`.
Its shifted threshold cutoff gives at least
`p(2D-m)/2` different unions, where

```math
m=\left\lfloor
\frac{\sqrt D\ln D}{a\sqrt2}+B\sqrt D
\right\rfloor,
\qquad a=\pi/\sqrt6,
```

for some fixed `B` and all sufficiently large `lambda`.
Equation (4) preserves their distinctness. The expansion in that proof gives

```math
\log_2\operatorname{rank}U_{H,D}
\ge\sqrt2\,\lambda
 -(3-2\sqrt2)\log_2\lambda-O(1).
\tag{5}
```

This bound is uniform over `H>=2D`: changing the padding does not
change the set of distinct union monomials.
There is no assertion that (5) is an upper bound.
Additional annihilation terms may increase the rank.

Both existing factor-vector families contain this map as a projection,
up to nonzero input scalings. The row projection is at `i=2H`.
The column projection is at `i=0` and uses
`b^nu=(H!/z_nu)t^nu`.
The same separate lossless-compression argument as in the binary-prefix
note therefore gives a width sum of at least
`2sqrt(2) lambda-O(log lambda)`.
Within that kind of compression, a bound
`w_r+w_c<=4rho lambda+o(lambda)` requires
`rho>=1/sqrt(2)`.

This is not a rank lower bound for the complete central piece.
Discarding directions invisible to the opposite factor or finding
a joint refactorization can still change the conclusion.

## 4. Structured multiplicity subsets are a different question

The full fixed-defect family is not the only possible choice.
A multiplicity box can have much smaller growth under unions.
For example, take

```math
S_H=\{(1^{H-2k},2^k):0\le k\le\lfloor H/2\rfloor\}.
```

If `M=|S_H|`, there are only `2M-1` distinct unions from `S_H×S_H`.
That count alone looks substantially better than the full-partition case.
However, exact finite calculations show that the highest channel
already has full pair rank in the following cases:

| `H` | Labels `M` | Distinct unions | Highest-channel rank |
| ---: | ---: | ---: | ---: |
| 2 | 2 | 3 | 4 |
| 3 | 2 | 3 | 4 |
| 4 | 3 | 5 | 9 |
| 5 | 3 | 5 | 9 |
| 6 | 4 | 7 | 16 |
| 7 | 4 | 7 | 16 |
| 8 | 5 | 9 | 25 |

Each rank equals `M^2` modulo both primes `65521` and `1000003`.
Since the matrix has `M^2` rows, either certificate proves the stated
rank over the rationals. These are finite statements only.
They show that union counts can substantially underestimate the
existing operator-family rank.

This two-mode family has another limitation: `M=O(H)`, so encoding
`lambda` bits requires exponential energy `H`.
Higher-dimensional multiplicity boxes can have exponentially many labels
at polynomial energy. They would need a new rank or evaluation theorem;
neither the union count nor the defect cutoff supplies one.

## 5. Reproduction and remaining target

Run:

```sh
python3 tools/check_near_ones_channels.py
```

The checker verifies the defect embedding, its union identity, (3),
the defect support relation for nonzero coefficients, and the finite
table under two primes. The asymptotic proof of (5) is the ordinary
union-count argument above.
The recorded output is [near-ones-channels.json](../evidence/near-ones-channels.json).

The constructive target left by this analysis is a label family with
polynomial energy and a proved small channel factorization after all
annihilation terms are included. An analysis based only on leading
monomial unions is insufficient. No such family or improved factoring
rate is proved here.

The underlying operators and crossing identity are prior work from
the pinned Anthropic development. This analysis and checker used
substantial OpenAI assistant support.
