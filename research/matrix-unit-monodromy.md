# A matrix-unit monodromy crossing and its product-boundary restriction

## Exact construction outside the local parity span

Fix `d>=2`, let `D=d^2`, and index an auxiliary `D`-dimensional
space by ordered pairs of external letters. Put `p=(x,y)` and
`q=(v,u)`, and define

```math
K(z;p,q)=[p=q]I_D+(1+z)E_{q,p}
```

over a field of characteristic two. Here `E_(q,p)` is an ordinary
matrix unit. For `t` positions and the boundary functional
`ell(M)=M_(0,0)`, with boundary pair `0=(0,0)`, define

```math
P_t(z;\mathbf x,\mathbf y,\mathbf u,\mathbf v)
 =\ell\left(\prod_{j=1}^t K(z;p_j,q_j)\right).
```

The product retains its position order. Since `K(1;p,q)=[p=q]I`,
the degree-at-most-`t` polynomial satisfies `P_t(1)=E_(d^t)`.
Its coefficients therefore give an exact crossing decomposition
with at most `t+1` pieces. An entry is evaluated with fixed auxiliary
dimension; complete coefficient tables and their ordinary field rank
factorizations cost `d^O(t)` operations. No claim that matrix factor
halves commute is made or needed for this direct construction.

For binary external letters, this family is outside the local
even-parity tensor span. At two successive positions choose
`(p_1,q_1)=(01,00)` and `(p_2,q_2)=(00,01)`, interpreting each
pair consistently as the corresponding two external letters.
Both positions have odd four-leg parity, yet the scalar polynomial
entry is `(1+z)^2`, which is nonzero. Padding by the boundary pair
at later positions preserves a nonzero polynomial witness.

## Exact auxiliary walk and useful block structure

Start the auxiliary state at the boundary letter. If `p=q`, the
state remains unchanged; the scalar weight is `z` when that common
letter is the current state and is one otherwise. If `p!=q`, the
entry vanishes unless `q` is the current state. In the remaining
case it changes the state to `p` and contributes `1+z`.
The final state must again be the boundary letter.

Thus every nonzero entry is exactly `z^a(1+z)^b`, where `a` counts
the first kind of weighted stay and `b` counts actual state changes.
This gives a particularly short independent entry algorithm.

Equivalently, expanding by selected matrix-unit positions
`i_1<...<i_k` gives the conditions

```math
q_{i_1}=0,\qquad p_{i_k}=0,\qquad
p_{i_j}=q_{i_{j+1}}\ (1\le j<k),
```

with `p_i=q_i` at unselected positions. Therefore the row cut
preserves the full alphabet histogram of `p` and `q`. For binary
external letters, the column cut preserves separately the counts
of ones in `x,v` and in `y,u`. These are exact direct-sum blocks,
not approximate rank heuristics.

## A full coefficient-span restriction for the product boundary

**Proposition.** For every decomposition of the crossing tensor
into pieces in the coefficient span of `P_t`, the rank sum is at
least `(d^2-1)^t`. In particular this construction with this boundary
cannot satisfy the improving finite-seed criterion.

First restrict the row matrix to words over the pair alphabet with
the boundary letter removed. A first nonmatching pair `p_i!=q_i`
cannot change the initial state, because its `q_i` is not the
boundary. Thus the entire restricted polynomial matrix is the
constant identity of order `(D-1)^t`.

For the column cut use row words `a_i=(x_i,u_i)` and column words
`b_i=(v_i,y_i)`, both from the same nonboundary pair alphabet.
Order pairs lexicographically by their first coordinate, then their
second, and order words by their first differing position. On a
diagonal entry the auxiliary state stays zero and every local pair
is nonzero, so the diagonal polynomial is one.

At the first position where the row and column words differ, all
earlier positions have kept the state zero. A nonzero off-diagonal
entry requires `q=(v,u)=0`. Since both column and row pair labels
are nonzero, this forces `x!=0` and `y!=0`. Hence the row pair is
`(x,0)` and the column pair is `(0,y)`, and the row is later in the
chosen order. The restricted column matrix is therefore lower
triangular with constant diagonal one.

For any coefficient functional `L` with `L(1)!=0`, the scalar tensor
`L(P_t)` consequently has rank at least `(D-1)^t` in both cuts.
If coefficient-span pieces sum to the crossing, applying their
functionals to a protected diagonal entry gives `sum_j L_j(1)=1`.
At least one piece has `L_j(1)!=0`; its contribution alone gives
the claimed rank-sum lower bound. This proof permits arbitrary
recombinations and characteristic-two field extensions.

Finally `d^2-1>d^(sqrt(2))` for every integer `d>=2`. For example,
`d^2-1>=3d^2/4`; the ratio `(3/4)d^(2-sqrt(2))` increases with
`d`, and its value at two exceeds one because `3>2^(sqrt(2))`.
The binary case gives the leading restriction
`rho>=log_2(3)/2`, approximately `0.79248`.

The later [exterior occupation-graph argument](exterior-monodromy-and-occupation-graphs.md)
strengthens this conclusion for the displayed one-particle boundary:
a singleton occupation is a rectangle, so it proves
`S_t>=d^(2t)/(t+1)` and hence `rho>=1` for the full coefficient span.
The protected restriction above remains a shorter independent check.

The proposition is specific to the displayed product boundary.
A general trace-one auxiliary functional is a different coefficient
family. The finite checks below include the rank-three diagonal
functional `M_11+M_22+M_33`; it has worse tested profiles, but the
proposition is not asserted for it.

## Exact computation

`matrix_unit_monodromy.cpp` uses the histogram blocks and packed binary
elimination to compute actual coefficient ranks. The product boundary
at seven positions, with 128 external labels, has row ranks

`[12288,12273,12537,11757,9951,8091,6555,4096]`

and column ranks

`[10923,12046,13415,12513,11172,9368,7208,4096]`.

Its rank sum is approximately `79061.776`, far above the improving
threshold. The rank-three boundary at six positions has rank sum
approximately `24896.711`, compared with `16490.162` for the product
boundary at that length. These are finite profile comparisons only.

`check_matrix_unit_monodromy.py` independently multiplies the local
polynomial matrices, checks the auxiliary walk, forms both full cut
matrices through three positions, and checks the protected triangular
restriction through four positions. The evidence is recorded in
`matrix-unit-monodromy-verification.json`. The all-length proposition
is proved above, independently of these finite calculations.
