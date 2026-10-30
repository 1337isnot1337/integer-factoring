# General-characteristic trace-shear parameter and its two cut ranks

The binary construction and its sharper transition-count word bound were
already proved in the pinned [reduction-map Fierz note](https://github.com/1337isnot1337/integer-factoring/blob/bb3bf763ed054f6fdaf882649c037ad810e16b97/research/reduction-map-fierz-crossing.md).
This note extends its scalar normalization and opposite-cut ranks to
arbitrary characteristic. It also records an endpoint-specific
scalarization obstruction. No improved factoring exponent is claimed.

Let `F` be any field, let `q>=2`, and put `N(X)=Tr(X)I_q` on the auxiliary
space `V=Mat_q(F)`. Thus `N²=qN`, where the integer `q` is interpreted in
`F`. Choose `t!=0` with `kappa=1+tq!=0`. Then

```math
W=I+tN,\qquad W^{-1}=I-\frac{t}{\kappa}N,\qquad
W\eta=\kappa\eta,\quad\eta=I_q.
```

Define

```math
A_1=-\frac1t(E_{xy}\otimes E_{vu}),\qquad
A_2=\frac\kappa t W(E_{xu}\otimes E_{vy})W^{-1}.
```

The non-scalar boundary outputs cancel and the trace correction gives
`(A_1+A_2)eta=delta_(xv)delta_(yu)eta`. The original tensor factors and
those conjugated by `W` supply genuinely commuting halves for the two cuts.
As in the binary construction, `r_1=c_2=1`.

The two opposite ranks must be kept separate. Conjugating `A_1` into the
column factor structure gives the factors

```math
\left(I-\frac t\kappa |I\rangle\langle I|\right)
 \otimes (\operatorname{SWAP}+tI),
```

up to permutations and a nonzero scalar. The first factor is invertible.
The row flattening of `A_2` instead gives

```math
(I+t|I\rangle\langle I|)
 \otimes\left(\operatorname{SWAP}-\frac t\kappa I\right).
```

Therefore the exact profile is

```math
c_1=q^2\operatorname{rank}(\operatorname{SWAP}+tI),\qquad
r_2=q^2\operatorname{rank}\left(\operatorname{SWAP}-\frac t\kappa I\right).
```

In characteristic two with even `q` and `t=1`, both ranks are
`q³(q-1)/2`, recovering the previous note. In odd characteristic, when the
characteristic divides `q` and `t=-1`, they are **different**:

```math
c_1=\frac{q^3(q-1)}2,\qquad
r_2=\frac{q^3(q+1)}2.
```

The signs of the trace correction reverse between the two conjugations.
For example, over `F_3` at `q=3`, the exact profile is `(1,27),(54,1)`.
Direct construction of both full opposite flattenings confirms `27` and
`54`.

There is also an odd-characteristic parameter with equal antisymmetric
ranks: if `q=2` in the field and `t=-1`, then `kappa=-1` and
`-t/kappa=-1`. Both ranks become `q³(q-1)/2` again. For example this applies
at `q=5,8` over `F_3`. Here `W` is an involution with `W eta=-eta`; the
explicit scalar coefficients in the definition of the pieces account for
that boundary eigenvalue.

These constructions remain noncompetitive. They provide exact odd-field
catalytic endpoints and a warning against importing the symmetric binary
rank profile without checking both signs.

## A common product density exists but does not recover the endpoint

For `i!=j`, both `W` and its transpose fix the product vector `E_ij`.
Hence the rank-one product density
`D=|E_ij><E_ij|` is normalized and satisfies `W^(-1)DW=D`. Its trace functional
factorizes across both full pairs. This does not scalarize the catalytic
endpoint: for example, with `i=0,j=1`, the functional gives zero at the
crossing input `(0,0,0,0)`.

More strongly, for `q>=3` there is no functional which both factorizes
across the supplied full factor pairs and sends `A_1+A_2` to the crossing.
Such a functional would produce a scalar identity

```math
\delta_{xv}\delta_{yu}
 =a_{xy}b_{uv}+c_{xu}d_{yv}.                  \tag{1}
```

Here is a direct field-independent obstruction to (1). Under the row cut,
it says `P-ab^T=C tensor D`, with `P` the swap permutation matrix. The
right side has rank at least `q²-1`. For `q>=3`, a singular `q` by `q`
factor would give rank at most `q(q-1)<q²-1`, so `C,D` are both invertible.
Fix a nonzero vector `u`. The equation `b^T(u tensor v)=0` defines a
subspace of `v` of dimension at least `q-1`. For every vector in this
subspace, (1) gives

```math
Cu\otimes Dv=v\otimes u.
```

For nonzero `v`, equality of nonzero pure tensors forces `v` to be
proportional to the single vector `Cu`. This contradicts a subspace of
dimension at least two. Thus (1) is impossible.

The distinction is specific to the common-boundary endpoint. The
normalized common-product functional fixes the auxiliary identity, but
need not reproduce the desired scalar values on the endpoint family.

## Exact verification

The standard-library checker verifies all local endpoints for the stated
prime-field examples and reconstructs both opposite flattenings directly
at `p=q=3`. Run:

```sh
python3 tools/check_trace_shear_general.py \
  --output evidence/trace-shear-general-audit.json
```

The construction and exact checks used substantial OpenAI assistant support.
