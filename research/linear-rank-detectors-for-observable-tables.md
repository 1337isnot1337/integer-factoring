# Linear rank detectors for observable transition tables

An observable matrix table can have a crossing quotient without
a common crossing eigenvector. However, a linear map on the scalar
entry space can still force a lower bound on its word-rank spectral
transfer. For binary labels, the three standard pairings have a
particularly small detector: a `2 by 2` matrix of linear forms.

The conclusion concerns the spectral upper bound obtained from
exact individual-word ranks and transition counts. It is not a
lower bound on the true ranks of arbitrary grouped pieces after
additional cancellations.

## 1. A general detector lemma

Let `S` be a finite-dimensional scalar tensor space containing the
crossing `E`, and let `sigma(t)=sqrt(r(t)c(t))` be the geometric
mean of its two cut ranks. Suppose a field-linear map

```math
L:S\longrightarrow\operatorname{Mat}_{a,b}(F)
```

satisfies

```math
\operatorname{rank}L(E)=D>0,
\qquad \sigma(t)\ge\gamma\operatorname{rank}L(t)
\quad\text{for every }t\in S,                  \tag{1}
```

where `gamma>0` is real.

Let an `M by M` matrix table `T(q)` have entries in `S` and an
all-word crossing endpoint with fixed boundaries `c,d_0`:

```math
c^TT(q_1)\cdots T(q_k)d_0=\prod_{j=1}^kE(q_j).
                                                               \tag{2}
```

Put `N_ij=rank L(T_ij)` and `G_ij=sigma(T_ij)`. Let `s,t` be the
zero-one support vectors of `c,d_0`.

**Lemma.** The live word-rank spectral transfer satisfies

```math
\boxed{\operatorname{spr}(G_{\rm live})\ge\gamma D.} \tag{3}
```

**Proof.** Apply `L` independently to each of the `k` disjoint
input blocks in (2). The right side becomes `L(E) tensor ... tensor
L(E)`, of rank `D^k`. Each state path on the left gives a scalar
multiple of a Kronecker product of the matrices `L(T_ij)`. Its rank
is the product of their ranks. Rank subadditivity gives

```math
D^k\le s^TN^kt.                                \tag{4}
```

Thus the matrix `N` restricted to paths between these boundaries
has spectral radius at least `D`; otherwise its powers grow more
slowly than `D^k`. Entrywise, `G>=gamma N`, so monotonicity of the
spectral radius for nonnegative matrices proves (3). Removing
vertices that lie on no boundary-to-boundary path does not affect
(4). `QED`.

No matrix-entry commutation, common eigenvector, or triangularity
assumption occurs in this argument. The fixed-boundary all-word
endpoint is the essential input.

## 2. Binary three-pairing detector in characteristic two

Over any characteristic-two field, put

```math
E=[x=v][y=u],\qquad R=[x=y][u=v],\qquad
C=[x=u][y=v],
```

with binary external labels. These tensors are independent, so
`t=aE+bR+cC` has unique coefficients. Define

```math
\boxed{L(aE+bR+cC)=
\begin{pmatrix}a&a+b\\a+c&0\end{pmatrix}.}       \tag{5}
```

This map is injective and `rank L(E)=2`.

We compute both cut ranks explicitly. In the row cut, let `K`
be the binary swap matrix, and put
`r=(1,0,0,1)^T`, `s=(0,1,1,0)^T`. Then

```math
t_{\rm row}=aK+brr^T+cI.
```

If `a+c` is nonzero, `aK+cI` is invertible: its square is
`(a+c)^2 I`. Since `Kr=r` and `r^Tr=0`, its rank-one update by
`brr^T` remains invertible. If `a+c=0`, the matrix is
`a ss^T+b rr^T`, whose two blocks have disjoint supports. Hence

```math
r(t)=
\begin{cases}
4,&a+c\ne0,\\
[a\ne0]+[b\ne0],&a+c=0.
\end{cases}                                     \tag{6}
```

Exchanging the two cuts exchanges `b,c`, giving

```math
c(t)=
\begin{cases}
4,&a+b\ne0,\\
[a\ne0]+[c\ne0],&a+b=0.
\end{cases}                                     \tag{7}
```

The determinant of (5) is `(a+b)(a+c)`. When it is nonzero,
both ranks in (6)--(7) are four, and
`sigma(t)=4=2 rank L(t)`. When (5) has rank one, (6)--(7) give
`sigma(t)>=2`. The zero case is immediate. Therefore

```math
\boxed{\sigma(t)\ge2\operatorname{rank}L(t)}      \tag{8}
```

over every characteristic-two field. Over F2, equality holds for
all eight tensors; over larger fields, the additional cost
`sqrt(8)` can occur in a rank-one detector case.

Applying the general lemma with `gamma=D=2` proves:

**Corollary.** For every finite observable matrix table whose
binary scalar entries lie in `span(E,R,C)`, the live transition
cost matrix has spectral radius at least four. This holds over
every characteristic-two field, with arbitrary auxiliary size
and arbitrary finite observable-quotient endpoint.

Thus enlarging the auxiliary matrix, changing its basis, adding
nilpotent extensions, or adding other triangular diagonal
characters cannot produce a base below four while all entries
remain in this scalar space. In particular, no such table gives
an improving binary certificate through the word-rank spectral
bound. A later calculation of the actual grouped-piece ranks
would be a different method and is not excluded by this corollary.

## 3. An explicit scalar-space escape

Let `H=[x=y=u=v=0]`, a fully separable coordinate tensor. It is
independent of `E,R,C`, has both cut ranks one, and `E+H` has both
cut ranks three.

Any extension of (5) that kept the constant two in (8) would have
to send `H` to zero, since its cost is one. But then
`L(E+H)=L(E)` has rank two, while `sigma(E+H)=3<4`. Such an
extension is impossible. Adding this single coordinate tensor
therefore escapes this specific detector; it does not by itself
produce a better crossing construction.

## 4. Exact controls

`check_observable_rank_detector.py` checks all coefficients over
F2 and F4 against (6)--(8), the injectivity and crossing rank of
(5), the coordinate-tensor escape, and the Kronecker rank bounds
for the strict three-state observable example. The all-field and
all-state conclusions are the preceding proofs.
