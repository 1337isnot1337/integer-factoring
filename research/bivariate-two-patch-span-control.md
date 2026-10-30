# An exact bivariate trace-span control at eight labels

Keeping a tagged trace coordinate as a separate parameter gives the
[bivariate endpoint](tagged-two-patch-span-controls.md) a coefficient
space containing the previous two scalar spaces. Its dimension is exactly
`4t-1` for every `t>=1`. A separate complete finite calculation closes
one small instance: at `t=3`, tag `(0,0)`,
the exact minimum binary-span rank cost is **64**, equal to `d^2`.

The rank-cost conclusion concerns one fully specified instance. No other bivariate tag, larger
width, different trace basis, or extension-field coefficient combination
is excluded by this result.

## 1. Exact coefficient dimension at every width

For `L=F_(2^t)`, set `e=S+w`, `B_0=Se`, and `B_1=(S+1)e`.
The bivariate kernel uses the factor `1+u tau_j(B_i)` at the tagged
coordinate and factors `1+z tau_j(B_i)` at all other coordinates.
Let `tau:L->F_2^t` be the invertible trace-coordinate map, and fix one
tagged coordinate in either patch. Let `b` be that bit and let `h` be the
Hamming weight of the other `2t-1` trace bits. The attainable pairs are
exactly

```math
(b,h)=(0,h),\quad 0\le h\le2t-1,
\qquad
(b,h)=(1,h),\quad 0\le h\le2t-2.
```

The pair `(0,0)` occurs when both patches are zero, for example `w=S`.
For `b=0,h>=1`, set the corresponding coordinate in the other patch to
one and choose `h-1` of the remaining `2t-2` coordinates to be one.
For `b=1`, set that corresponding coordinate to zero and choose `h`
of the remaining coordinates to be one. In either construction the two
trace vectors differ. Their field values therefore differ, and the
assignments `e=B_0+B_1`, `S=B_0/e`, `w=S+e` realize the chosen pair.
The omitted pair `(1,2t-1)` would make both trace vectors all ones.
Then `B_0=B_1` forces `e=0` and both patch values zero, a contradiction.

At a point of type `(b,h)`, the bivariate entry is

```math
H(z,u)=(1+ub)(1+z)^h.
```

Its coefficient functions are `A_k=binom(h,k) mod2` and
`C_k=b binom(h,k) mod2`. The function `C_(2t-1)` is zero. The remaining
`4t-1` functions are independent: restricting a relation to the points
`b=0,h=0,...,2t-1` gives a unit lower-triangular Pascal matrix, so every
`A` coefficient in the relation is zero. Restriction to
`b=1,h=0,...,2t-2` gives the same argument for the remaining `C`
coefficients. This proves the exact dimension `4t-1`. The triangular
evaluation minor has determinant one, so the dimension remains the same
over every characteristic-two coefficient extension.

Kernel independence transfers to the conservation tensors because every
`(S,w)` is realized at external labels `(0,S,w,S+w)`, where the
conservation indicator equals one. The argument includes `t=1`, whose
three attainable types are `(0,0),(0,1),(1,0)`.

This is an all-width dimension theorem. It gives no all-width optimum
for the rank cost.

## 2. Exact eight-label instance and coefficient basis

Use `L=F_2[X]/(X^3+X+1)`, `alpha=X`, and
`tau_j(a)=Tr(alpha^j a)` for `j=0,1,2`. Set

```math
e=S+w,\qquad B_0=Se,\qquad B_1=(S+1)e,
```

and define

```math
H(z,u;S,w)=
 (1+u\tau_0(B_0))
 \prod_{(i,j)\ne(0,0)}(1+z\tau_j(B_i)).
```

Multiply this kernel by the four-leg conservation condition
`[x+y=u_label+v]`, substituting `S=x+y,w=x+u_label`.
Here `u_label` denotes the external tensor label; `u` in `H` is the
independent polynomial parameter. At `(z,u)=(1,1)`, the tensor is exactly
the crossing `E`.

Write its coefficients as `A_0,...,A_5,B_0,...,B_5`, corresponding to
`sum_j z^j A_j + u sum_j z^j B_j`. The last coefficient `B_5` is zero,
because both trace vectors cannot be all ones. Exact binary elimination
finds eleven independent coefficient tensors, so the full binary span
has exactly 2,048 elements. The checker obtains this basis by elimination
and records the expression of every original coefficient in the basis.
It does not assume that the nominal coefficients are independent.

The coefficient symbols in this paragraph are tensors; the earlier
`B_0,B_1` are field-valued patch variables in the defining kernel.

## 3. Complete rank and potential certificate

Every span tensor receives two independent constructions of its cut
matrices: the additive-convolution block construction and the literal
full four-leg flattenings. All 4,096 row/column comparisons agree.

For a nonzero basis mask `a`, assign the exact integer lower cost

```math
L(a)=\left\lfloor\sqrt{r(Q_a)c(Q_a)}\right\rfloor.
```

The complete XOR-potential certificate has value zero at the zero mask
and value 64 at the target mask 2047. Its 4,192,256 inequalities verify

```math
\ell(x\mathbin{\oplus}a)\le\ell(x)+L(a)
```

for every state `x` and every nonzero span mask `a`. Thus any finite
decomposition of the crossing into tensors from this binary span has
total rank cost at least 64. The crossing itself has both cut ranks 64,
so its single-piece cost is exactly 64. This supplies the matching upper
bound and proves the exact optimum.

The argument includes every number of pieces, every binary recombination,
and every grouping within this particular span. No rational refinement
of the lower costs and no floating-point comparison is needed.

## 4. Reproduction and scope

[Exact checker](../tools/check_bivariate_two_patch_span.py) and
[complete finite certificate](../evidence/bivariate-two-patch-span-verification.json).

```sh
python3 tools/check_bivariate_two_patch_span.py \
  --output /tmp/bivariate-two-patch-span-verification.json
```

The compact evidence stores the independent basis indices, all coefficient
expressions, all rank pairs, all lower costs, and the complete potential.
The checker requires status `passed` only after both rank constructions
cover all 2,048 tensors and all potential inequalities have been checked.
It atomically preserves partial progress and exits with status 2 if its
default 110-second computation limit is reached. The saved result is
complete; no timeout or unfinished wider scan contributes to the claim.

The bivariate endpoint and linear coefficient-count bound remain useful
interfaces. This one eight-label binary span does not improve the trivial
one-piece cost, and supplies no improvement to the factoring exponent.
