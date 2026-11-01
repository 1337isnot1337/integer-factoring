# A noncommuting matrix auxiliary with a linear-rank leading channel

The [F4 braid path quotient](quantum-braid-path-quotient.md) admits a fixed
two-state auxiliary deformation. It preserves both strict crossing endpoints
and has degree `3n^2/2`. Its highest coefficient has **both cut ranks exactly
L** on every common L-letter path alphabet containing the zigzag path.
This is an all-length factorization of one coefficient. The rank sum for the
whole polynomial remains open; the complete factoring score stays `1/2`.

## 1. The auxiliary and the strict endpoints

Use `alpha^2=alpha+1`, the height-three local generators `e_i`, and

```math
Q=\begin{pmatrix}0&1\\1&\alpha\end{pmatrix},\qquad
N=\begin{pmatrix}0&1\\0&0\end{pmatrix}.
```

Thus `Q^2=I+alpha Q`, `Q^5=I`, `N^2=0`, and `QN != NQ`.
Replace each of the `D=n^2` gates of the block braid by

```math
G_i(t)=\bigl(I+t(Q-I)\bigr)\otimes I
       +\bigl(tI+t(1-t)N\bigr)\otimes e_i.        \tag{1}
```

All signs coincide in characteristic two. Start the auxiliary in `e_0`.
At zero the full product is identity. At one, each gate is
`Q tensor I+I tensor e_i`: the regular two-dimensional F4 representation
of the original F16 braid gate `qI+e_i`. The block endpoint is therefore
`Q^(D/2) tensor SWAP`, by the all-length braid endpoint already proved.

Put `k=-D/2 mod5`, `tau=(0,1,0,alpha,alpha+1)[k]`, and use the final
auxiliary row

```math
\ell=(1,\tau)Q^k.
```

The explicit five powers of Q give `ell e_0=1`; also
`ell Q^(D/2)e_0=1`. Therefore projecting the completed product onto
middle-height-zero path boundaries gives exactly `I` at zero and
`SWAP` at one. This applies to every common restriction of the path
alphabet. No local projection or ordinary invariant-subspace assumption
is used: the negligible-wall quotient justifies the full TL words.

For binary width lambda take the least even n with at least `2^lambda`
height-three return paths, and retain the last `2^lambda` lexicographic
paths. The zigzag is the last path, so it is retained. The existing path
count/unranking algorithm, with the offset `total-2^lambda`, remains
polynomial time. Here `n=O(lambda)`. Constructing all intermediate paths
of length 2n, or just propagating a single input column, costs
`2^O(lambda)` operations with polynomial-degree coefficient lists.
This supplies the allowed effective entry construction.

## 2. Degree and the unique highest auxiliary word

Every degree-two local term is `N tensor e_i`. Between any two such
terms a nonzero product requires a degree-one factor containing Q:
`N I N=N N N=0`. A factor before the first degree-two term is also
needed, since `Ne_0=0`. Consequently k degree-two terms require at
least k other gates. Thus `k<=D/2` and the degree is at most `3D/2`.

At that degree all remaining gates have degree one. There is exactly
one surviving pattern: degree-one, degree-two, repeated D/2 times
in application order. The degree-one term must be `(Q-I) tensor I`:
the other possibility `(I+N) tensor e_i` cannot reset N. Since
`N(Q-I)e_0=e_0`, the auxiliary factor is one, including the final row.
Hence the top coefficient is the boundary block of

```math
E_n=e_{i_D}e_{i_{D-2}}\cdots e_{i_2},           \tag{2}
```

where `i_1,...,i_D` is the block-braid word in application order.

## 3. An explicit all-length diagram and factorization

Write n=2m. Consecutive rows of (2) apply the commuting products

```math
e_j e_{j+2}\cdots e_{j+n-2},\qquad j=0,1,\ldots,n-1.
```

The following boundary invariant proves its diagram without a large
path calculation. After row j, input block P is capped in adjacent
pairs. For `0<=a<j`, input `Q_a` connects to output position a; for
`j<=a<n`, input `Q_a` still connects to output position `n+a`.
Output positions `j,...,j+n-1` are paired consecutively.
The initial row gives the invariant directly. Shifting the next
matching one position right creates a single alternating chain: it
transfers `Q_j` to output j and leaves the shifted matching. It creates
no closed loop. This proves the invariant by induction.

At the final row the connections are:

* adjacent caps on P;
* `Q_a` through to `R_a` for `a<n-1`;
* `Q_(n-1)` to `S_(n-1)`, and `R_(n-1)` to `S_0`;
* adjacent cups on `S_1,...,S_(n-2)`.

Let Z be the zigzag return path, and let `H(S)` be the maximum height
of S. The exact boundary coefficient is

```math
\boxed{\quad T_{3D/2}(P,Q,R,S)
       =[P=Z]\,[Q=R]\,[H(S)\le2].\quad}        \tag{3}
```

To see the factorization, the P caps select Z by orthogonality of the
fusion frame. Contracting S against its dual leaves a planar TL
endomorphism of one strand. Every such diagram is a scalar times
identity; after this contraction the Q-to-R diagram is identity.
The S factor lies in the common cup image of
`e_1,e_3,...,e_(n-3)` on the vacuum path space. This image has dimension
one: remove these m-1 cups, leaving the one-dimensional two-strand
vacuum space.

Its coordinate vector is proportional to `[H(S)<=2]`. Indeed apply
those commuting generators to Z. At every affected site the prefix
height is one and the input pair is DU. The local path formula sends
it to `alpha^2(DU+UD)`. Independent choices of those m-1 pairs produce
exactly all height-at-most-two return paths, with equal coefficient
`alpha^(2m-2)`. Finally the `(Z,Z,Z,Z)` matrix entry is one: closing
the displayed diagram with adjacent vacuum cups and caps gives n
loops, canceling the norm factor `alpha^n` of the 2n-strand zigzag.
This fixes the scalar and proves (3).

On a common alphabet containing Z, (3) is a product of a nonzero
P-vector, an identity Q/R channel, and a nonzero S-vector. Both actual
cut ranks are therefore L. It also proves that the degree bound is
attained. This direction survives the exact binary restriction above.
It is symmetric between the two rank bounds; a useful decomposition
of the remaining coefficients is still required.

## 4. Finite evidence and a design restriction

The optional NumPy implementation is
[quantum_braid_matrix_aux.py](../tools/quantum_braid_matrix_aux.py).
The [replay](../tools/check_quantum_braid_matrix_aux.py) checks complete
endpoints and top coefficients at n=4,6,8, all boundary connectivities
at even lengths through128, and independent scalar propagations over
F256. Sixty-four distinct values exceed the small polynomial degrees,
so they give complete polynomial comparisons for the sampled entries.
The source-pinned [receipt](../evidence/quantum-braid-matrix-auxiliary.json)
records exact counts and coefficient-span dimensions.

```sh
make verify-braid-matrix-aux
```

There is a uniform restriction on simpler auxiliary templates. If
every gate is `A(t) tensor I+b(t)I tensor e_i`, with the same matrix
A(t) at every site and scalar b(t), the two whole gate weights commute.
Expanding by the number j of selected e terms gives

```math
\sum_{j=0}^D A(t)^{D-j}b(t)^j\otimes A_j,
```

where `A_j` are the original unprojected braid coefficients. Every
scalar boundary projection stays in their original span. Noncommuting
individual coefficients of A(t) do not evade this identity. The two
whole weights in (1) do not commute, so that reduction does not apply.

The span enlargement and the linear-rank top coefficient are separate
facts. Neither establishes a favorable whole-family rank sum. No new
Lean claim or score above `1/2` is made. This construction builds on
the earlier ambient braid, local path quotient, and norm-dual frames;
their notes retain the underlying source and contributor attribution.
