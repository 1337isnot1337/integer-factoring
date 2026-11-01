# Height cutoffs and almost-linear extreme pieces

The [paired-path construction](paired-path-crossings.md) extends from
height three to any **fixed** height `H>=2`. This gives a quantitative
tradeoff between label entropy and the uniform extreme-piece rank.
For every fixed `epsilon>0`, a sufficiently large fixed cutoff gives
three extreme pieces with both ranks `O(L^(1+epsilon))`.
This strengthens individual pieces, not the complete factoring score.

## Construction and fixed-field cost

Use local splittings
`U_h:H_(h+1) direct-sum H_(h-1) -> H_h tensor F^2`
for `0<=h<=H`, with `dim H_h=h+1`. Retain paths in `0..H` ending
at zero, and use inverse duals as before. The proof of pairedness
is unchanged. A common binary initial segment preserves the endpoints.

There are two explicit choices for the weight-preserving blocks:

1. Classical blocks `[[1,-m],[1,h-m+1]]` over a fixed prime field
   of characteristic `p>H+1`. Their determinant is `h+1`.
2. The q-integer blocks of the paired-path note over a fixed binary
   extension containing an element of odd order `M>H+1`.
   Their determinants are `[h+1]_q`, which are nonzero because
   `[k]_q=0` precisely when `M` divides `k`.

For the second choice take a primitive element of `F_(2^d)` with
`2^d-1>H+1`. Both constructions have uniformly specified local maps.
The explicit implementation in `paired_path_crossing.py` remains
the height-three specialization; no arbitrary-height entry replay is
claimed for that tool. For any fixed `H`, the same dynamic program
has at most `(H+1)^4` magnetic states and polynomial-time entries.
The field and height are constants in the label-width parameter.

## Exact exponential rates

Write `L_(n,H)` for the number of retained label paths and `s_(n,H)`
for the number of supported physical words. A supported word has
prefix imbalance in `[-H,H]` and final imbalance zero, on primal
and inverse-dual legs alike. Let

```math
f_H=2\cos\frac\pi{H+2},\qquad
g_H=2\cos\frac\pi{2H+2}.
```

For even `n` tending to infinity,

```math
L_{n,H}=\Theta_H(f_H^n),\qquad
s_{n,H}=\Theta_H(g_H^n).                         \tag{1}
```

Here is a direct spectral proof. The adjacency matrix of the path
with `k` vertices has eigenvectors
`v_j(t)=sin(j*t*pi/(k+1))`, eigenvalues
`2*cos(j*pi/(k+1))`, and squared norm `(k+1)/2`.
Its `n`th power counts walks. The two extreme eigenvalues contribute
positively at even length. Apply this to the label path with `H+1`
vertices, starting at its first vertex, and to the support path with
`2H+1` vertices, starting at its middle vertex. The leading constants
are respectively
`4*sin(pi/(H+2))^2/(H+2)` and `2/(H+1)`.

The support-channel proof therefore bounds both ranks of every
extreme six-vertex piece by

```math
s_{n,H}=O_H(L_{n,H}^{\beta_H}),\qquad
\beta_H=\frac{\log g_H}{\log f_H}.               \tag{2}
```

For a least even binary-selected length, the full retained label
count is within a fixed `H`-dependent factor of the selected count.
Indeed, (1) bounds the ratio at consecutive even lengths for all
sufficiently large lengths; finitely many smaller ratios have a
finite maximum. Thus (2) holds for the exact binary alphabets too.

Since `f_H<g_H<2` and both tend to two, `beta_H>1` and
`beta_H -> 1`. Given any fixed `epsilon>0`, choose a fixed height
with `beta_H<1+epsilon`. Its field remains fixed as width grows.
This proves the almost-linear extreme-piece assertion without
using a growing field or a growing-height advice sequence.

The elementary uniform bound `beta_H<4/3` already holds for every
`H>=3`. At `H=3` it is the radical inequality in the
[bounded-support note](bounded-spin-support-ranks.md). For `H>=4`,
`f_H>=sqrt(3)` and `g_H<2<3^(2/3)<=f_H^(4/3)`;
the middle strict inequality is just `8<9`. At `H=2`, instead,
`beta_H=log(3)/log(2)>sqrt(2)`. Height three is the first cutoff
whose support estimate puts an extreme piece below the post-half
rank-cost threshold.

## What this does and does not resolve

For each fixed `H>=3`, the same channel construction gives

```math
r(T_{ABC}),c(T_{ABC})
 \le s_{n,H}\frac{n!}{(2A)!(2B)!(2C)!}.
```

Sublinear corner bands of width `o(n/log n)` still have total rank
`poly(n)g_H^n exp(o(n))`, whereas the crossing has rank of order
`f_H^(2n)` in each cut. Here `f_H^2>=phi^2>2>g_H` for `H>=3`.
So increasing the cutoff alone cannot make these bands a complete
decomposition. The almost-linear pieces also do not supply the
small-cut exponent below one required by the
[asymmetric rank obligation](paired-path-crossings.md#5-evidence-and-the-remaining-question).

The relevant open question remains a uniform compatible family of
asymmetric coefficient directions. Height is an additional parameter
for that problem, with explicit constant-field constructors and exact
entropy rates. The complete score remains `1/2`.

Run `make verify-paired-height` for exact walk recurrences and local
determinant controls at larger cutoffs. Those finite checks support
the implementation formulas; the all-height conclusions above follow
from the displayed constructions and spectral argument. This note
and checker used substantial Codex AI assistance.
