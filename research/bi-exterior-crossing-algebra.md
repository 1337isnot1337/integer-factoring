# An exact exterior-algebra model for both crossing cuts

## Result

Over every field of characteristic two, a concrete four-dimensional
space of one-bit four-leg tensors gives an exact model of the
square-zero algebra on two generators. Tensoring this space produces
an algebra model in which **both cut ranks are multiplication ranks**.
This is an exact construction, not an associated-graded lower bound.

Let

```math
\mathcal A_\lambda
 =F[s_1,n_1,\ldots,s_\lambda,n_\lambda]/(s_i^2,n_i^2),
\qquad d=2^\lambda.
```

Write `M_f` for multiplication by `f` on this `d^2`-dimensional
algebra. Define the linear involution `H` independently at every
position by

```math
H(1)=sn,\qquad H(s)=s,\qquad H(n)=n,\qquad H(sn)=1.
\tag{1}
```

There is an explicit injective linear map `Phi` from `A_lambda`
to four-leg tensors on `d` labels such that

```math
r(\Phi(f))=\operatorname{rank}M_f,\qquad
c(\Phi(f))=\operatorname{rank}M_{Hf},
\qquad
\Phi\left(\prod_{i=1}^\lambda(1+s_in_i)\right)=E_d.
\tag{2}
```

Thus an identity decomposing the displayed algebra element yields an
actual crossing decomposition with these two exact ranks. No such
identity improving the current factoring exponent is claimed here.

## 1. The local matrices

Use the row index `(x,y)` and column index `(u,v)`, in binary order
`00,01,10,11`. Put

```math
L=\begin{pmatrix}0&1\\0&0\end{pmatrix},\quad U=L^T,
\quad A=L\otimes I_2+I_2\otimes L,
\quad B=U\otimes I_2+I_2\otimes U.
```

Explicitly,

```math
A=\begin{pmatrix}0&1&1&0\\0&0&0&1\\0&0&0&1\\0&0&0&0\end{pmatrix},
\qquad
B=\begin{pmatrix}0&0&0&0\\1&0&0&0\\1&0&0&0\\0&1&1&0\end{pmatrix}.
\tag{3}
```

Let `E_2` swap the two tensor factors, and set `N=E_2+I_4`.
Direct multiplication in characteristic two gives

```math
A^2=B^2=0,\qquad AB=BA=N.
\tag{4}
```

Define `Phi(1)=I_4`, `Phi(s)=A`, `Phi(n)=B`, `Phi(sn)=N`.
The vector `w=e_01` is cyclic: the ordered vectors

```math
(w,Aw,Bw,Nw)=(e_{01},e_{00},e_{11},e_{01}+e_{10})
\tag{5}
```

form a basis. In that basis, the operators `A,B` are exactly
multiplication by `s,n` on `F[s,n]/(s^2,n^2)`.
This proves the first rank equality in (2) locally.

## 2. The other cut

Let `C` reshuffle a row matrix into the column cut:

```math
C(T)_{(x,u),(y,v)}=T_{(x,y),(u,v)}.
```

Put `R=C(I_4)` and `J=C(N)`. In four-bit row-mask notation
(least-significant bit first in each row), the matrices are

| Matrix | Row masks |
|---|---|
| `R` | `(9,0,0,9)` |
| `J` | `(8,4,2,1)` |
| `A'=J C(A)` | `(2,0,9,2)` |
| `B'=J C(B)` | `(4,9,0,4)` |

Direct multiplication gives

```math
J^2=I_4,\quad JR=R,\quad (A')^2=(B')^2=0,
\quad A'B'=B'A'=R.
\tag{6}
```

The vectors `(e_00,A'e_00,B'e_00,Re_00)` form a basis:
their binary column masks are `(1,4,2,9)`.
Thus `A',B'` are another regular square-zero representation.
For `f=a+bs+cn+e sn`,

```math
J C(\Phi(f))=aR+bA'+cB'+eI_4,
```

which is multiplication by `Hf` in this second regular basis.
Multiplication by the invertible `J` preserves rank. This proves
the second equality in (2) locally.

## 3. Tensor powers and effectiveness

Tensor the local map across the positions. Reordering the bits in
each pair of indices turns the row representation into the regular
representation of `A_lambda`. The column reshuffle also factors
position by position, and multiplication by `J^tensor(lambda)` is
invertible. This proves both rank equalities for every polynomial
`f`, including arbitrary mixtures of degrees.

Finally, `Phi(1+sn)=I_4+N=E_2`. Tensoring gives the endpoint in
(2), exactly and with scale one.

There are `4^lambda=2^O(lambda)` algebra coordinates and each cut
matrix has this many rows and columns. Explicit coefficients,
entry evaluation, Gaussian elimination, and rank factorizations
therefore use `2^O(lambda)` scalar operations. Over `F_2`, or any
fixed effective finite field, this also gives that bit-time bound.
For growing fields, their representation, construction, and scalar
arithmetic costs must additionally satisfy the separate effective
field interface; matrix dimensions alone do not establish it. A useful
asymptotic factoring construction must additionally supply the
required small number of pieces and a better rank profile.

## 4. Every invertible linear change of the generators preserves both ranks

The map `H` is full monomial complementation followed by swapping
`s_i` with `n_i` in every pair. Let `P` be that pair-swap matrix.
For an invertible linear substitution `G` on the `2lambda`
generators, the top-degree pairing gives

```math
H\circ\Lambda G
 =\det(G)\,\Lambda(PG^{-T}P)\circ H.
\tag{7}
```

Here `Lambda G` denotes the induced algebra automorphism, and
`G^{-T}` is the inverse transpose. The usual signs in exterior
duality disappear in characteristic two. Formula (7) can also be
checked by pairing complementary monomials and using that `G`
acts on top degree by its determinant.

Algebra automorphisms conjugate multiplication operators. The
determinant is nonzero, so it does not affect rank. Consequently
both ranks in (2) are invariant under **every** invertible linear
change of all `2lambda` generators, including transformations that
mix the two generator sets.

## 5. A cheap actual direction that an associated-graded model missed

At `lambda=4`, consider

```math
F=\sum_{A\subseteq\{1,2,3,4\},\ |A|=2}s_A n_{A^c}.
\tag{8}
```

Its coefficient matrix between two-element subsets is the
complement permutation and has rank six. Also `HF=F`.
Exact multiplication has rank `36`, so (2) gives an **actual**
four-leg tensor with cut ranks `(36,36)`.

For comparison, inserting the same polynomial into a model with
Boolean selectors `s_i^2=s_i` gives ranks `(46,46)`.
The two constructions are different tensors. The square-zero
model above avoids that lift's rank increase because it supplies
a different compatible local representation of both cuts.

This example refutes the suggested inequality
`rank(M_F) rank(M_HF)>=d^2 rank(C_F)`:
`36^2=1296<256*6=1536`. It is not itself a crossing decomposition.
In particular, summing all coefficients of
`product_i(s_i+z n_i)` gives `product_i(s_i+n_i)`, not the crossing
element `product_i(1+s_i n_i)`.

Position-invariant families inside this four-letter span remain
subject to the separately established
[all-field binary-word obstruction, pinned at `3bdff157`](https://github.com/1337isnot1337/integer-factoring/blob/3bdff15730f9eb0a316cce03c842cb1962093aa5/research/all-field-binary-word-barrier.md).
Its theorem applies to every piece invariant under simultaneous
position permutations, with no restriction on the local tensors or
on the number of pieces. Therefore its scope includes the full
four-letter local span used here; some piece must have both ranks
at least `2^(3lambda/2-O(log(lambda+1)^2))`.
The exact model makes arbitrary position-dependent and mixed-degree
constructions accessible; it does not remove that obstruction.

The separate [degree-budget theorem](bi-exterior-degree-budget.md) now
excludes the strict rank-sum seed criterion throughout this entire
algebra through `256` labels, including arbitrary inhomogeneous pieces
and supplied formal limits. A binary seed in this model must start at
`512` labels or above. This does not limit crossing tensors outside the
model, and does not assert that a larger successful seed exists.

## Replay

```sh
python tools/check_bi_exterior_crossing_algebra.py \
  --output evidence/bi-exterior-crossing-algebra.json
```

The checker verifies both local regular representations, the endpoint,
independent dense four-tensor flattenings for random small algebra
elements, the actual `(36,36)` example, the independent Boolean-selector
`(46,46)` control, and rank invariance under invertible generator
changes. The Boolean control expands `n_i=1+g_i` in the binary group
algebra and computes every row and column convolution-block rank.
The all-size result is the tensor-product
argument above, not an inference from the finite controls.

The construction, proof, and exact computations used substantial
OpenAI assistant support. They are ordinary mathematics and code,
with no new Lean formalization.
