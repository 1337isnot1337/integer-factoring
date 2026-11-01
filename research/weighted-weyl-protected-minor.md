# A protected minor in the weighted Weyl symmetry span

## Result and scope

At energy `h=32`, the weighted Weyl crossing on `d=202` binary
partition labels has a common `1821` by `1821` protected minor in
both flattenings. It detects the crossing tensor. Consequently, any
decomposition of the crossing into arbitrary elements of the full
coefficient span, enlarged by every crossing-preserving permutation
of the four external labels, has one piece of both ranks at least
`1821`.

The strict integer comparison

```math
1821^{408}>202^{577},\qquad 577^2>2\cdot408^2
```

proves `1821>202^sqrt(2)`. Thus this entire finite symmetry span
fails the strict finite-seed rank-sum criterion. The result holds
over every field of characteristic two, with any number of pieces
and arbitrary coefficients in that field. It is stronger than
checking the 17 raw coefficients or enumerating their binary sums.
It does not give an all-energy obstruction or an improved factoring
bound.

## 1. The curve and its diagonal

Use the weighted Weyl curve from
[weighted-weyl-and-multivariate-extraction.md](weighted-weyl-and-multivariate-extraction.md),
with weights `w_j=2^j`, labels of weighted energy `h`, and coordinate
order `(alpha,y,u,beta)`:

```math
P(z)(\alpha,y,u,\beta)
=[\alpha+u=\beta+y]\,[t^\alpha]
\prod_j(1+t_j)^{\beta_j}
        (1+(1-z)^{w_j}t_j)^{y_j}.
\tag{1}
```

The equality in the indicator is equality of integer occupation
vectors. It is not reduced modulo the field characteristic.
The crossing is `E(alpha,y,u,beta)=[alpha=beta][y=u]`, and
`P(1)=E`.

On this crossing diagonal, Lucas' identity gives the exact formula

```math
P(z)(\alpha,y,y,\alpha)
=z^{\sum_jw_j(\alpha_j\mathbin{\&}y_j)}.
\tag{2}
```

Here `&` means bitwise intersection of nonnegative integers.
For one occupation coordinate, the coefficient being extracted is

```math
\sum_k\binom{\alpha_j}{k}\binom{y_j}{k}(1+z)^{w_jk}
=\prod_{b:\,2^b\text{ occurs in }\alpha_j\mathbin{\&}y_j}
  \bigl(1+(1+z)^{w_j2^b}\bigr)
=z^{w_j(\alpha_j\mathbin{\&}y_j)}.
```

The second equality uses Lucas' identity; the last uses the fact
that every `w_j2^b` is a power of two. Multiplying over coordinates
proves (2).

## 2. One minor works for every external symmetry

There are eight permutations of the four external positions that
preserve the crossing pairing. Under each permutation, the
occupation support in (1) becomes one of

```math
\alpha+u=\beta+y,\qquad
\alpha+y=\beta+u.
\tag{3}
```

The diagonal value (2) is unchanged, since it is symmetric in
`alpha,y`.

Choose a set `M` of ordered label pairs `(alpha,y)` such that:

1. `alpha_j & y_j=0` for every coordinate;
2. the integer sums `alpha+y` are all distinct;
3. the integer differences `alpha-y` are all distinct.

In the row cut, take rows `(alpha,y)` and matching columns
`(beta,u)=(alpha,y)` indexed by `M`. For two different elements of
`M`, the first support condition in (3) would equate their
differences; the second would equate their sums. Both are excluded.
Every diagonal entry is one by (2). Thus the restriction of every
symmetry image of the polynomial curve is exactly `I_|M|`.

In the column cut, use rows `(alpha,u)=(alpha,y)` and matching
columns `(beta,y)=(alpha,y)`. The same two uniqueness conditions
again kill all off-diagonal entries. The restriction is the same
identity matrix.

It follows that the coefficient at `z^0` of every symmetry image
restricts to the identity in both cuts, and every positive-degree
coefficient restricts to zero. Therefore there is a well-defined
linear functional `ell` on the full coefficient/symmetry span
whose restrictions in both cuts are

```math
T\big|_{M\times M}=\ell(T)I_{|M|},\qquad \ell(E)=1.
\tag{4}
```

The functional is intrinsic: it is any one of these diagonal
entries of the tensor. Its definition does not depend on choosing
a linearly independent coefficient representation.

If `E=sum_i T_i`, at least one `ell(T_i)` is nonzero, over any
characteristic-two field. That same piece has both ranks at least
`|M|`. Hence `sum_i sqrt(r_i c_i)>=|M|`, proving the claimed
finite-seed obstruction whenever `|M|>d^sqrt(2)`.

## 3. The exact energy-32 witness

Construct a bipartite graph. Left vertices are differences
`alpha-y`, right vertices are sums `alpha+y`, and there is an edge
for every ordered pair with coordinatewise bitwise intersection
zero. A matching gives exactly the set required above.

At energy 32, the exact counts are:

| Quantity | Count |
| --- | ---: |
| Binary partition labels | 202 |
| Eligible ordered label pairs | 9,882 |
| Distinct differences | 4,402 |
| Distinct sums | 1,821 |
| Matching size | 1,821 |

The supplied witness matches every right vertex. The checker
reconstructs the labels and graph, computes the matching, and
independently verifies all three conditions on the returned pairs.
It also checks the local diagonal binomial identity over the full
occupation ranges and verifies the eight symmetry support types.
No large tensor dump or numerical rank tolerance is used.

```sh
python3 tools/check_weighted_weyl_matching.py --output evidence/weighted-weyl-protected-minor.json
```

The evidence includes all label occupations and all matched label
pairs. The integer power comparison certifies the strict threshold;
the nearby decimals `1821` and `202^sqrt(2)≈1820.7871` are not used
as proof. A different energy or a family outside this coefficient
and external-symmetry span requires a new argument.
