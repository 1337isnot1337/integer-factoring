# A support barrier for Boolean-incidence crossing kernels

## Result

Let `G=F_2^lambda` be the label group, let `d=2^lambda`, and let
the coefficient field be arbitrary. For a scalar kernel `K:G×G→F`,
define the four-leg tensor

```math
T_K(x,y,u,v)=[x+y=u+v]K(x+y,x+u).
\tag{1}
```

The crossing tensor is `T_delta`, where `delta(S,w)=[S=w]`.
Suppose every piece in a crossing decomposition satisfies

```math
K_i(S,w)=0\quad\text{unless}\quad w\subseteq S,
\tag{2}
```

where inclusion refers to the supports of the binary labels. Then

```math
\sum_i\sqrt{r(T_{K_i})c(T_{K_i})}\ \ge\ d^{3/2}.
\tag{3}
```

There is no restriction on the number of pieces or their coefficients.
The ambient kernel space has dimension `3^lambda`. It includes the
unsymmetrized diagonal-word span and arbitrary row-dependent scalar
multiples and recombinations of its elementary coefficient tensors.
For polynomially many pieces, (3) forces contraction coefficient
`rho>=3/4` to leading order. It cannot improve the current factoring
exponent `1/2`.

The theorem is a support obstruction for the specified kernels. It is
not a lower bound for arbitrary crossing decompositions or factoring.

## 1. The two convolution cuts

For a function `f:G→F`, let `C_f` be the matrix
`C_f(a,b)=f(a+b)`. In the row cut, fix `S=x+y=u+v` and use `x,u`
as free coordinates. In the column cut, fix `w=x+u=y+v` and use
`x,y` as free coordinates. Consequently,

```math
r(T_K)=\sum_{S\in G}\operatorname{rank}C_{K(S,\cdot)},
\qquad
c(T_K)=\sum_{w\in G}\operatorname{rank}C_{K(\cdot,w)}.
\tag{4}
```

These block decompositions involve only permutations of rows and
columns and hold over every coefficient field.

**Support lemma.** If `f` is nonzero and has `k` nonzero coefficients,
then `rank C_f>=d/k`.

Indeed, every row is a translate of `f` and has support size `k`.
Choose `r=rank C_f` rows spanning the row space. Their supports must
cover all `d` columns: otherwise every row would vanish in an uncovered
column, whereas every column is a nonzero translate of `f`. The union
has size at most `rk`, proving `rk>=d`. This argument does not use
Fourier diagonalization or a characteristic assumption.

## 2. Diagonal support pays at least square-root dimension per label

For a kernel satisfying (2), put

```math
D_K=\{S:K(S,S)\ne0\}.
```

For `S in D_K`, the row `K(S,·)` is nonzero and has support size
at most `2^|S|`. The column `K(·,S)` is nonzero and has support
size at most `2^(lambda-|S|)`. By (4) and the support lemma,

```math
r(T_K)\ge\sum_{S\in D_K}2^{\lambda-|S|},
\qquad
c(T_K)\ge\sum_{S\in D_K}2^{|S|}.
\tag{5}
```

Cauchy--Schwarz in the form
`(sum a_S)(sum b_S)>=(sum sqrt(a_S b_S))^2` gives

```math
r(T_K)c(T_K)\ge d\,|D_K|^2.
\tag{6}
```

If `sum_i K_i=delta`, every one of the `d` diagonal positions belongs
to at least one `D_(K_i)`. Therefore `sum_i |D_(K_i)|>=d`.
Taking square roots in (6) and summing proves (3).

A useful extension allows each piece to satisfy its own incidence
condition after an invertible linear change `S→L_i S,w→L_i w`.
Such a change preserves the two convolution ranks and permutes the
diagonal positions, so (6) is unchanged. This does not assert that an
arbitrary linear combination of tensors from different incidence
spaces satisfies (6).

## 3. Which selector searches this closes

The elementary row polynomial is

```math
P_S(z)=\prod_{j:S_j=1}(1+zN_j),\qquad N_j=1+g_{e_j}
\quad\text{in }F[G].
\tag{7}
```

Every coefficient is supported on `w subset S`. Multiplying any
coefficient by an arbitrary function of `S` preserves this support.
Thus unrestricted row controls, arbitrary coefficient recombinations,
and arbitrary choices of final pieces inside their span all obey (3).
The controls need not have bounded Boolean degree or any symmetry.

Likewise, a relation correction that only uses group directions
supported inside `S` preserves (2). A correction using a direction
outside `S` can escape this theorem. Nonlinear common-label
permutations that do not preserve the group convolution blocks are
also outside its assumptions.

The support condition is essential. For the constant kernel `K=1`,
all `d` diagonal entries are nonzero, but every convolution block
has rank one. Hence `r=c=d`, violating (6) when `d>1`.

## 4. A separate bound for arbitrary conservation kernels

Even without (2), every nonzero `T_K` satisfies

```math
\sqrt{r(T_K)c(T_K)}\ge d.
\tag{8}
```

Let `u,v` be the numbers of nonzero rows and columns of `K`, and
let `m` be its number of nonzero entries. The support lemma and
Cauchy--Schwarz give

```math
r(T_K)\ge d\sum_{\text{nonzero rows}}\frac1{k_S}
       \ge\frac{d u^2}{m},
\qquad
c(T_K)\ge\frac{d v^2}{m}.
\tag{9}
```

Since `m<=uv`, (8) follows. Equality occurs for the constant kernel.
This fact bounds the number of nonzero pieces in a proposed small
rank-sum seed within the full conservation space. For example, a
passing `d=128` seed there could use at most seven nonzero pieces.
It does not provide the stronger `d^(3/2)` decomposition bound
without the incidence hypothesis.

## Replay

```sh
python tools/check_incidence_kernel_rank_sum.py \
  --output evidence/incidence-kernel-rank-sum.json
```

The checker exhausts all binary kernels at `d=4` for (8), exhausts
the binary incidence space through `lambda=2`, checks larger random
incidence kernels over `F_2` and `F_3`, and compares the convolution
formulas with independent dense four-tensor cuts. It includes the
constant off-incidence negative control. The all-size and all-field
claims follow from the support proof above, not the finite controls.

The proof and computations used substantial OpenAI assistant support.
They are ordinary mathematics and code without new Lean formalization.
