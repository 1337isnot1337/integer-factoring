# Complete ternary quadratic multipoint seeds with pure finite values

## Result and scope

Let `K(z;x,y,u,v)` have coefficients in `F_3`, degree at most two,
binary physical letters, and the strict endpoint

```math
K(1;x,y,u,v)=[x=v][y=u]=E(x,y,u,v).
```

Suppose both matrix cuts `(x,y)|(u,v)` and `(x,u)|(y,v)` have rank
at most one at each of `z=0` and `z=2`. For this entire family, the
[multipoint Hermite criterion](multipoint-hermite-crossings.md) using
the points `0,2,infinity` requires **`rho >= 3/4`**. This includes
all allocations of the degree budget and independent positive scales
in the two cut estimates. Therefore these seeds cannot improve the
existing `rho=1/sqrt(2)` crossing or the factoring score one half.

This is a finite-field construction obstruction. It is not a lower
bound for factoring, for the actual ranks of every tensor-power jet,
or for the separate two-place coefficient method. It does not cover
larger endpoint ranks, coefficients outside `F_3`, or interpolation at
new points of a coefficient extension. The same profiles and certificate
remain valid over extensions when the seed and three points are retained.

## 1. A small complete parameterization

A nonzero four-factor tensor of rank one in the first cut is
`A(x,y)B(u,v)`. Its second-cut matrix is a row/column permutation of
`A tensor B`, whose rank is `rank(A) rank(B)`. Thus rank one in both
cuts is equivalent to a pure tensor `a_x b_y c_u d_v`.

There are `2*4^4=512` nonzero pure tensors over `F_3`: each of the four
factors has one of four projective directions, with two total scalars.
Include the zero tensor to obtain 513 possible values at either point.
Every pair `A=K(0), B=K(2)` determines exactly one seed:

```math
K(z)=A+(B-E)z-(A+B+E)z^2.                         \tag{1}
```

Constant changes of physical basis

```math
G_x\otimes H_y\otimes H^{-T}_u\otimes G^{-T}_v
```

preserve `E` and both local Smith profiles. On a nonzero pure tensor,
the two pairings `a^T d` and `b^T c` are each zero or nonzero. Normalize
the `x,y` factors to `e_0`. A nonzero pairing lets the paired factor
become `e_0`; a zero pairing lets it become `e_1`. If either pairing is
zero, a diagonal stabilizer absorbs the remaining scalar. When both
are nonzero, their total contraction is the invariant scalar, either
one or two. There are consequently five nonzero orbits, plus zero.

The checker independently enumerates `GL_2(F_3)` and the resulting
pair-factor orbits. Their disjoint sizes are

```text
1, 144, 144, 96, 96, 32; total 513.
```

Fixing the first value to these six representatives and allowing all
513 second values therefore covers every seed in the stated family.
The **3,078** representatives have **23** ordered three-place profile
signatures, using **14** distinct local pairs of cut profiles. Both cuts
have generic rank four because the endpoint `E` is invertible.

At infinity the matrix is `t^2 K(1/t)`, including its full degree shift.
Using degree two for an actual lower-degree polynomial causes no
loophole: increasing the homogeneous degree adds the same integer to
every infinity order and to its allocated slope, leaving the local
pressure unchanged. Lower-degree criteria embed in this one.

## 2. An exact degree-budget certificate

For a local pair of ordered Smith profiles `nu,mu`, choose prefix lengths
`k,l` with `k*l >= 8`. Set

```math
b=\max\left(\frac1k\sum_{i=1}^k\nu_i,
             \frac1l\sum_{i=1}^l\mu_i\right).        \tag{2}
```

If the allocated slope `d >= b`, then for every `s,r >= 1`, the
arithmetic-geometric mean inequality gives

```math
\begin{aligned}
\left(\sum_i s^{-\nu_i}\right)s^d&\ge k,\\
\left(\sum_i r^{-\mu_i}\right)r^d&\ge l.
\end{aligned}                                        \tag{3}
```

Hence that point's two-cut pressure is at least `log(8)`.
For each of the 23 signatures, the saved certificate chooses such
prefixes at every point with

```math
b_0+b_2+b_\infty\le2.                              \tag{4}
```

The actual sums are `4/3`, `5/3`, or `2`. Every admissible degree
allocation has `d_0+d_2+d_infinity=2`. Let `J` be the nonempty set
of places with positive slope. If `d_j < b_j` at every `j` in `J`,
then `2=sum_J d_j < sum_J b_j <= sum_all b_j <= 2`, a contradiction.
Thus a positive-slope place has `d_j >= b_j`. Equations (3)--(4)
force the maximum pressure of the contributing places to
be at least `log(8)=4*(3/4)*log(2)`, as claimed. This argument uses
only rational means and the exact integer comparison `k*l>=8`.
It optimizes over every scale and allocation without a floating search,
including the omission of any zero-slope places.

## 3. Replay and independent controls

```sh
python3 tools/search_multipoint_rank_one_seeds.py \
  --output evidence/multipoint-rank-one-seeds.json
python3 tools/check_multipoint_rank_one_seeds.py \
  --output evidence/multipoint-rank-one-validation.json
```

The search computes Smith orders with minimum-valuation formal-series
pivots at precision nine, justified by the determinant degree bound eight.
The checker reruns the entire catalogue and compares every saved field.
For a witness of every signature it independently recovers both complete
profiles from full block-Toeplitz ranks: **138 cut profiles**. It also
reconstructs all **5,888 entries** of 23 two-copy crossings from jets at
the three places. These controls retain the strict endpoint, both cuts,
and the infinity shift. The source hashes are saved in the validation.

The next constructive question is whether larger local rank products,
additional coefficient-field points, or higher polynomial degree admit
a strict complete pressure witness below `1/sqrt(2)`. Finite favorable
profiles would still require the effective transfer's complete costs.

This argument and its exact computation used OpenAI Codex assistance.
