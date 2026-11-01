# Root-four bit circuits for the semiprime transform

## Scope and source

This ordinary arithmetic analysis uses the pinned
[`RiemannGaussian` source](https://github.com/dbsanfte/RiemannGaussian/tree/f5fd3ddb407972f9ca65298fbe9be3225fe06635).
It studies the exact word operations before claiming any global circuit-width
or factoring result. The source defines the transform ring as
`Z/(2^T+1)`, uses root four, and stores canonical residues in `T+1` bits:
[`SemiprimeMacroFermatNetwork.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroFermatNetwork.lean#L21-L49),
[`SemiprimeMacroBitNetwork.lean`](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroBitNetwork.lean#L21-L45).
The source's explicit twiddle stage multiplies by a power of four
([`MacroBitNetwork`, lines 52–78](https://github.com/dbsanfte/RiemannGaussian/blob/f5fd3ddb407972f9ca65298fbe9be3225fe06635/RiemannGaussian/SemiprimeMacroBitNetwork.lean#L52-L78)).

## Twiddles are signed cyclic rotations

Let `M=2^T+1`. Then `2^T=-1` in `Z/MZ`, so `4^T=1`. Reduce a twiddle
exponent to `0<=a<T` and write `2a=epsilon*T+s`, where
`epsilon` is 0 or 1 and `0<=s<T`. A canonical residue has the form

```math
x=h2^T+\sum_{i=0}^{T-1}b_i2^i,
\qquad h\in\{0,1\},\quad h=1\Longrightarrow b_i=0\text{ for all }i.
```

For `4^a x=2^{2a}x`, the low bits move by the cyclic shift `i -> i+s mod T`.
Each moved bit receives sign
`(-1)^(epsilon+floor((i+s)/T))`; the top bit contributes
`(-1)^(1+epsilon) h 2^s`. At destination `s`, the two contributions combine
to `(-1)^epsilon (b_0-h)2^s`. Canonicality makes `b_0-h` one of `-1,0,1`.
Thus the result before canonical reduction is a signed binary sum

```math
S=\sum_{j=0}^{T-1}c_j2^j,\qquad c_j\in\{-1,0,1\},\qquad
|S|\le 2^T-1,
```

congruent to `4^a x` modulo `M`. Its canonical representative is `S` when
`S>=0`, and `S+M` when `S<0`.

This gives a fixed `O(T)`-gate circuit. Read the signed digits from low to high
with carry state `k_j`:

```math
c_j+k_j=y_j+2k_{j+1},\qquad y_j\in\{0,1\},\quad k_j\in\{-1,0,1\},\quad k_0=0.
```

The terminal carry is `-1` exactly when `S<0`. If it is `-1`, increment the
low `T` output bits once; the increment carry is the canonical top bit. This
computes `S+2^T+1` and handles `S=-1` by producing the valid top-only word
`2^T`. Every position uses a constant-size signed-carry and increment cell.
The rotation wires are a permutation, not a generic constant multiplier.

## Modular additions and pointwise products

For canonical `x,y` in `[0,M-1]`, addition computes `x+y<2M` and conditionally
subtracts `M` once. Subtraction computes the signed difference in
`[-(M-1),M-1]` and conditionally adds `M` when negative. Ripple carry, one
sign/comparison chain, and per-bit selection use `O(T)` Boolean gates for
either operation. Therefore one radix-2 butterfly—two modular add/subtract
outputs and its root-four twiddle—has an `O(T)`-gate implementation.
The `T+1`-bit canonical-word condition is `h=1` only when every low bit is
zero; checking it takes `O(T)` one-bit gates per word.

For a pointwise product, schoolbook multiplication computes an integer
`P=x*y` with `O(T^2)` gates. Write `P=L+H2^T` with `0<=L<2^T` and
`0<=H<=2^T`. Since `2^T=-1`, the residue is `L-H mod M`; one signed
correction by `M` returns its canonical `T+1`-bit word. The usual row-wise
schoolbook circuit keeps one `O(T)`-bit accumulator and the multiplier word,
so an individual product gadget has a path decomposition of width `O(T)`:
process multiplier rows in order, keeping both operand words and the current
and next accumulator rows in each bag. Each partial-product row and its add
carry then meets only adjacent row bags.

## Redundant signed base-four digits defer transform carries

There is a second exact representation that removes canonical carry
normalization from the transform itself. Put `q=T/2`, so
`M=2^T+1=4^q+1`, and define

```math
\Phi(d)=\sum_{i=0}^{q-1}d_i4^i\pmod M,
\qquad d\in\mathbb Z^q.
```

Every residue has a representative with digits in `{0,1,2,3}`, except the
top residue `4^q`, represented by `(-1,0,...,0)`. Since `4^q=-1 mod M`,
the signed rotation

```math
R(d)_0=-d_{q-1},\qquad R(d)_i=d_{i-1}\quad(1\le i<q)
```

satisfies `Phi(R(d))=4 Phi(d)` for every integer digit vector, with no carry
or size growth. Thus every root-four twiddle is a wire permutation and sign
change in this representation, even after digits have become redundant.
Addition is coordinatewise. Multiplication is the negacyclic convolution

```math
c_k=\sum_{i+j=k}d_i e_j-\sum_{i+j=k+q}d_i e_j,
\qquad 0\le k<q,
```

because each wrap by `q` powers of four contributes `4^q=-1`. If
`|d_i|<=A` and `|e_i|<=B`, then `|c_k|<=qAB`.

For example, a `t`-term transform of canonical input words has digit bound
`3t`: rotations preserve the input bound three and summing `t` terms adds
their digit vectors. A following pointwise product has digit bound at most
`9q t^2`, and an `r`-term transform then has bound at most `9rq t^2`.
Any final scalar ring multiplier can be charged as one more negacyclic
convolution. These bounds are polynomial in `T`, so every signed digit needs
only `O(log T)` bits. Canonicalization can be postponed to the interface
with the pointwise multiplier and to final output; those conversions cost
polynomially many bit operations. This gives a valid ordinary arithmetic
representation, not a faster evaluator by itself.

The representation changes the separator question rather than removing it.
Using the direct negacyclic-convolution implementation creates all `q^2`
digit products. Its factor graph contains `K_{q,q}` between the two sets of
digit inputs; contracting `q-1` edges of a matching gives a `K_{q+1}` minor,
so each product gadget has treewidth at least `q`. Converting to canonical
binary words before multiplication instead retains the existing `O(T)`-width
schoolbook product gadget. Either choice leaves at least order-`T` local
multiplier width at each of the `T` frequency positions. This local bound is
not a lower bound for the whole network and does not determine its treewidth.

There is also an exact gauge obstruction in the top-level twiddle array.
Its signed-rotation offset is `a(j,k)=jk mod T`; on every in-range unit
plaquette,

```math
a(j+1,k+1)-a(j+1,k)-a(j,k+1)+a(j,k)=1\pmod T.
```

Adding independent row and column offsets `u(j)+v(k)` leaves this mixed
difference unchanged. Therefore such row/column reindexings cannot remove
all cell rotations at once. A general per-cell reindexing can move the
rotations onto neighboring transform wires, so this is not a lower bound on
all circuit representations or on their separators.

With `T` transform positions and `log_2 T` butterfly stages, canonical
`T+1`-bit words give `O(T^2 log T)` transform gates and `O(T^3)` gates for
all `T` pointwise products. Keeping `q` redundant digits of `O(log T)` bits
instead gives `O(T^2 (log T)^2)` transform gates. Direct schoolbook
negacyclic convolution uses `O(T^3 (log T)^2)` gates for all pointwise
products; canonicalizing at the product interface retains the `O(T^3)`
binary multiplication count and adds polynomial conversion work. Each
variant is polynomial in the original input length because
`T^2=Theta(ell)`. These gate-count bounds do not bound the treewidth of the
full composition.

## The global separator is still the question

The previous
[dense-table analysis](sanftenberg-transform-treewidth-floor.md) gives a
`Q_d` clique-blow-up minor and a width floor
`Omega(ell/sqrt(log ell))` for extensional four-word butterfly tables.
The circuits above change that local representation: add/subtract use carry
chains and the root-four twiddle is a signed rotation. But the complete
arithmetic network also connects `T` variable pointwise products, both forward
transforms, the inverse transform, and the coefficient/carry recovery. A
local `O(T)`-width product gadget does not by itself give a sub-`ell/sqrt(log
ell)` separator for this connected network.

### The ideal twisted-fiber model has a sharp half-rate floor

There is a useful conditional model for the missing global argument. Let
`H(B,q,pi)` have a `q`-vertex clique `F_v` for every vertex `v` of a simple
base graph `B`. For each base edge `uv`, join `F_u` to `F_v` by an arbitrary
perfect matching; `pi` denotes the collection of these permutations. Define

```math
\mu_{1/6}(B)=
\min_{\substack{U\subseteq V(B)\\
|U|,\,|V(B)\setminus U|\ge \lfloor n/6\rfloor}}
\nu(E_B(U,V(B)\setminus U)),\qquad n=|V(B)|.
```

Here `nu` is maximum matching size. For `|V(B)|>=6`,

```math
\operatorname{tw}(H(B,q,\pi)) \ge q\mu_{1/6}(B)-1,
```

independently of every matching permutation. To see this, take any balanced
separator `S` of `H`. If `|S|>=q|V(B)|/6`, use
`mu_{1/6}(B)<=floor(|V(B)|/6)`. Otherwise the components of `H-S`, each of
size at most half the graph, can be grouped into two unions of size at least
`q|V(B)|/6`: take one component of that size if it exists, and otherwise
accumulate components until the union first reaches that size. Since each
fiber is a clique, its surviving vertices lie in one group. Assign a fully
deleted fiber arbitrarily. The resulting base cut has both sides of size at
least `|V(B)|/6`. Every matching edge across a crossing base edge must meet
`S`; choosing a maximum matching of crossing base edges makes all `q`-edge
fiber matchings vertex-disjoint, so `|S|>=q mu_{1/6}(B)`. The balanced
separator characterization of treewidth gives the claim.

For `B=Q_d`, this parameter is sharp:

```math
\mu_{1/6}(Q_d)=\Theta(2^d/\sqrt d).
```

The upper bound uses a median Hamming-ball cut; its crossing edges meet a
central layer of size `O(2^d/sqrt(d))`. For the lower bound, the crossing-edge
graph of any cut is bipartite. By Koenig's theorem its maximum matching has
the size of a minimum vertex cover `C`. If `|C|` is already at least a fixed
fraction of `2^d`, the desired bound follows. Otherwise deleting `C` leaves
at least a fixed positive fraction of the cube on each side and no edge
between them. Thus the external vertex boundary of either side is contained
in `C`. Harper's vertex-isoperimetric theorem says this boundary is minimized
by an initial segment of the simplicial order. For a set whose size is between
fixed positive fractions of `2^d`, this prefix contains all layers below some
`k=d/2+O(sqrt(d))` and a subset `A_k` of layer `k`. Its external boundary
contains both the unused part of layer `k` and the upper shadow of `A_k` in
layer `k+1`. Double-counting incidences gives
`|upper(A_k)| >= |A_k|(d-k)/(k+1) >= |A_k|/2` for sufficiently large `d`.
Thus the boundary is at least a fixed fraction of `binom(d,k)`, which is
`Theta(2^d/sqrt(d))` by Stirling's estimate. Consequently
`|C|=Omega(2^d/sqrt(d))`.

There is also a matching upper bound on the bundle treewidth: lift every bag
of a tree decomposition of `B` by replacing each base vertex with its whole
fiber. Therefore

```math
\operatorname{tw}(H(B,q,\pi)) \le q(\operatorname{tw}(B)+1)-1.
```

The hypercube treewidth estimate from the
[dense-table analysis](sanftenberg-transform-treewidth-floor.md) gives
`tw(Q_d)=Theta(2^d/sqrt(d))`. Thus the ideal twisted-fiber graph has
`tw(H(Q_d,q,pi))=Theta(q2^d/sqrt(d))` for every choice of twists. With
`q=T/2`, `T=2^d`, and `T^2=Theta(ell)`, this is
`Theta(ell/sqrt(log ell))`: twists alone cannot take this graph-decomposition
route below the current half-rate scale. The vertex-isoperimetric input is
Harper, [*Optimal numberings and isoperimetric problems on graphs*](https://doi.org/10.1016/S0021-9800(66)80059-5).

This settles the matching-cut order for the abstract model, not its occurrence
in the actual bit circuit. The pointwise-product clique minors and the
butterfly connections must still be realized by pairwise-disjoint branch sets
in the same graph. Until that embedding is proved, the result is not a
treewidth lower bound for the source circuit. Even with the embedding, it
would obstruct dense tree-decomposition evaluation of this representation,
not factoring algorithms in general.

The next proof obligation is now specific: either construct the circuit-to-
bundle minor with disjoint branch sets, which would close this representation
at the half-rate width floor, or identify exactly which carry/circuit feature
prevents that embedding and could permit a smaller decomposition. Any
candidate improvement still needs the entire bit circuit, coefficient/carry
recovery, the source's semiprime promise, and complete primitive-machine cost
accounted for. No stronger factoring exponent or Lean result is claimed here.

No upstream code was copied. Substantial OpenAI assistant support was used in
the source and circuit analysis.
