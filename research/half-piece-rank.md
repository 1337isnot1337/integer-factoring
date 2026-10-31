# Actual ranks of the partition pieces: a finite audit

## Result and scope

The central partition piece has full row and column rank for every weight
`h=1,...,6` checked here. Its matrix has `p(h)^2` rows and columns.
The actual binary-label restrictions also have full rank for label widths
`lambda=1,2,3`.

These are finite results over the rational numbers. They give no
asymptotic lower bound for the pieces, the clique evaluator, or factoring.
They do not improve the current factoring exponent.

The construction comes from Appendix B of the pinned
[Anthropic development](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
Its integer definitions are
[`Boson/Defs.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Boson/Defs.lean)
and
[`Boson/Defs2.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Boson/Defs2.lean).
The public description of the pieces appears in
[`AppHalfPiecesStatement.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).

## What the checker computes

A partition is a nondecreasing tuple of positive integers. The checker
uses exact integer coefficients for the following operators.

For `a=t_m a'`, the coefficient recursion is

```math
\begin{aligned}
[t^\pi]a[r]t^\sigma
={}&\sum_{l\in\mathrm{supp}\pi}
\binom{l-1}{m-1}
[t^{\pi\setminus l}]a'[r-l]t^\sigma\\
&+\sum_{k\in\mathrm{supp}\sigma}
(-1)^{m-1}\binom{k+m-1}{m-1}k\,m_k(\sigma)
[t^\pi]a'[r+k]t^{\sigma\setminus k}.
\end{aligned}
```

The first sum vanishes when `l<m`. Creation has no output multiplicity.
Annihilation has the input multiplicity `m_k(sigma)`. This distinction
is necessary for the integer identity.

The checker expands

```math
a^\xi=\prod_{j\in\xi}(-1)^{j-1}
\left(\sum_{k=1}^{j}\binom jk t_k\right),\qquad
b^\nu=\frac{h!}{z_\nu}t^\nu.
```

It then constructs the central piece, and every piece needed for the
small identity checks. Rows use `(xi,eta)` and columns use `(nu,upsilon)`.
The second flattening exchanges `eta` and `upsilon`.

For weights at most three, two independent operator orders produce the
same entries. Their intermediate weights are `s+i` and `2h-s-i`.
Their sum over `-2h<=s<=2h` equals `h!` times the crossing matrix.
The computation checks 1,202 exchange entries and 98 summed entries
over the integers.

The rank computation uses Gaussian elimination over each of two prime
fields, with primes `65521` and `1000003`. Both computations give full
rank. A nonzero determinant modulo a prime implies a nonzero integer
determinant. Thus these computations also prove full rank over the
rationals. They do not use floating-point rank thresholds.

The checker also selects the first `2^lambda` partitions in the source's
specified order: descending part lists, sorted lexicographically.
This restriction matters because the source uses binary labels rather
than every partition at a weight.

## Recorded results

| Weight `h` | `p(h)` | Matrix dimension | Row rank | Column rank | Declared channels |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1 | 1 | 1 | 1 | 4 |
| 2 | 2 | 4 | 4 | 4 | 12 |
| 3 | 3 | 9 | 9 | 9 | 30 |
| 4 | 5 | 25 | 25 | 25 | 67 |
| 5 | 7 | 49 | 49 | 49 | 139 |
| 6 | 11 | 121 | 121 | 121 | 272 |

For the binary labels, `(lambda,h)` is `(1,2)`, `(2,4)`, and `(3,6)`.
Both ranks are respectively `4`, `16`, and `64`.

The declared central channel count is
`R(h)=sum_{j=0}^{2h} p(j)`. It first becomes smaller than the full
matrix dimension at `h=15`: `R(15)=28629<30976=p(15)^2`.
The small full-rank results therefore do not contradict the published
asymptotic compression. Directly constructing the weight-15 matrix
would require almost one billion entries.

## Relation to the asymptotic channel-map result

The separate [top-channel proof](top-channel-rank.md) establishes
`log2 rank(U_h)=kappa sqrt(2h)+o(sqrt h)` for the highest intermediate
channel map. Its small rank checks use this same checker.
That theorem concerns one factor, rather than the complete partition piece.
The [binary-prefix extension](binary-prefix-top-channel.md) proves the
corresponding rank bound on the source's actual binary labels.
Neither theorem removes the possible cancellations in the full piece.

## A failed positive-Gram shortcut

At weight two, the natural central matrix is

```math
A_0=\begin{pmatrix}
70&0&0&34\\
0&2&-6&0\\
0&-10&-6&0\\
10&0&0&10
\end{pmatrix}.
```

The basis order is `((11,11),(11,2),(2,11),(2,2))`.
The middle principal minor is `-72`, and the full determinant is
`360*(-72)=-25920`. The matrix is also nonsymmetric.
Thus this natural matrix is not a positive semidefinite Gram matrix.
No real change by the same congruence on both sides makes it one.

This counterexample does not exclude every separate row and column
transformation. Such a transformation needs an explicit identity.
The [top-channel note](top-channel-rank.md) states the remaining
product-rank problem and the limits of the current argument.

## Reproduction

Run the finite checker from the repository root:

```sh
python3 tools/check_partition_channels.py --max-weight 6
```

The recorded output is
[`evidence/partition-channel-ranks.json`](../evidence/partition-channel-ranks.json).
The checker uses only the Python standard library. It supports weights
through eight, with substantially larger time and memory requirements.
The recorded evidence uses six.

Two negative controls detect omitted alternating signs and omitted
annihilation multiplicities. Neither control changes the source theorem.
No Lean proof or implementation of the full large-input evaluator was
created for this experiment.

The partition-piece construction is prior work of the pinned Anthropic
development. The finite checker, channel-map proof, and audit used substantial OpenAI
assistant support and retain the repository's Apache-2.0 license.
