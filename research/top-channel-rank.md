# A sharp rank bound for the highest partition channel

There is an asymptotic result for one factor of the central piece.
It does not assert a rank bound for the complete piece.

Let `U_h` send a basis pair `(xi,eta)` of partitions of `h` to
`(t^xi)[h] t^eta` in `V_(2h)`. Put
`kappa=pi sqrt(2/3)/ln(2)`. Then

```math
\log_2\operatorname{rank}_{\mathbb Q} U_h
=\kappa\sqrt{2h}+o(\sqrt h).
\tag{1}
```

The domain here contains **all** partitions of `h`, at both input positions.
This is the full-domain `i=2h` row channel map, up to nonzero row signs.
The source algorithm instead restricts each label to its first `2^lambda`
partitions at the selected weight. The separate
[binary-prefix extension](binary-prefix-top-channel.md) proves that the
same leading rank bound survives this restriction. With the ordinary
Hardy--Ramanujan asymptotic, its quantitative error is `O(log lambda)`.

For the full domain, the identification follows directly from the formula.
Indeed, `binom(w+h,2h)` vanishes for `w<h`.
The homogeneous part of `a^xi` at weight `h` is
`(-1)^(h-length(xi)) t^xi`.


The channel map comes from the Appendix B construction in
[`AppHalfPiecesStatement.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean),
at commit `16a5ce41304fcdaea4e0be7c132e1c361dc26da5`.
Its integer operator recursion is
[`Boson/Defs.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Boson/Defs.lean).

## Proof: independent leading monomials

Suppose that a partition `pi` of `2h` splits as `xi union eta`,
where both parts have weight `h`. In the normal-order formula,
an annihilation reduces the number of output parts by two relative
to the term with creation only.

In that creation-only term, each creation index is at least its
corresponding part of `xi`. Their sums both equal `h`.
Every index therefore equals that part. The coefficient is one.
Thus

```math
(t^\xi)[h]t^\eta=t^\pi+
\text{terms with fewer than }\operatorname{length}(\pi)\text{ parts}.
\tag{2}
```

Choose one split for each splittable partition. Their images are
linearly independent: in a linear relation, the terms with the largest
number of parts have distinct leading monomials. Removing those terms
and repeating proves that every coefficient vanishes.

## Proof: enough partitions split

Set `t=floor(h^(3/4))` and `m=2h-t`. Take a partition of `m` whose
largest part is at most `t`, then append exactly `t` parts of size one.
The result has a subpartition of weight `h`.

To see this, greedily add original parts while the sum stays at most `h`.
At termination the gap is less than `t`, unless the sum already equals
`h`. The added ones fill that gap. The construction is injective,
because removing exactly `t` ones recovers the original partition.

The number of excluded partitions of `m` is at most
`m p(m-t-1)`: remove one largest part, which exceeds `t`.
The source's Hardy--Ramanujan bounds imply

```math
\log_2\frac{m p(m-t-1)}{p(m)}
\le -c h^{1/4}+O(\log h)\longrightarrow-\infty
```

for a fixed positive `c`. Thus at least `p(m)/2` partitions survive
for all sufficiently large `h`. Their logarithm is
`kappa sqrt(2h)-O(h^(1/4)+log h)`.
Equation (2) gives the same lower bound for `log2 rank(U_h)`.
The upper bound `rank(U_h)<=p(2h)` completes (1).

The partition estimates used here are the pinned source's
[`Part/Bounds.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Part/Bounds.lean).
The argument is an ordinary proof, without a new Lean formalization.

The finite checker gives the following values:

| `h` | Splittable partitions | Rank of `U_h` modulo `1000003` | `p(2h)` |
| ---: | ---: | ---: | ---: |
| 1 | 1 | 1 | 2 |
| 2 | 3 | 4 | 5 |
| 3 | 6 | 9 | 11 |
| 4 | 14 | 21 | 22 |
| 5 | 25 | 42 | 42 |
| 6 | 53 | 77 | 77 |

## The remaining product-rank problem

Equation (1) does not prove a rank lower bound for the complete piece.
The piece is a product `UV` through all intermediate weights up to `2h`.
Lower weights can create dependencies or cancellations in the product.
Separate lower bounds on `rank(U)` and `rank(V)` do not remove that problem.

In particular, the total channel count is of order
`sqrt(h) p(2h)` under the partition asymptotics.
Bounds of order `p(2h)` on the two factor ranks are insufficient for
Sylvester's inequality `rank(UV)>=rank(U)+rank(V)-R(h)`.
An asymptotic refactorization obstruction needs a further argument.
This note leaves that argument open.

The declared-count obstruction in the existing evaluator audit is a
different statement. It concerns the source's chosen widths, without
asserting that those widths equal the actual matrix ranks.

The binary-prefix extension also bounds the rank of the actual row
factor map. That stronger domain statement still leaves the product-rank
problem open.

## Reproduction and attribution

Run the supporting finite checks with:

```sh
python3 tools/check_partition_channels.py --max-weight 6
```

The [finite audit](half-piece-rank.md) explains the checker and its scope.
The output is [partition-channel-ranks.json](../evidence/partition-channel-ranks.json).
Finite computations support their stated cases. The ordinary proof above
supplies the asymptotic statement.

The partition operators are prior work from the pinned Anthropic development.
This channel-map argument and checker used substantial OpenAI assistant support.
The repository's Apache-2.0 license applies.
