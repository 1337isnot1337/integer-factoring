# A primary-algebra rank formula for the eight-site coherent crossing

## Result and scope

The [bounded-spin construction](bounded-spin-crossing.md) has thirteen
return-to-zero paths at length eight. Its nine coefficient tensors span a
six-dimensional space over `F_5`. In both cuts, after invertible column
operations, this is a commutative algebra with one double-root primary
component and four simple-root components. This gives an exact rank formula
for every coefficient combination over **every field extension of `F_5`**.

For every decomposition of the crossing into pieces from that coefficient
span, at least one piece satisfies

```math
\boxed{r(A)c(A)\ge2304=48^2>13^3.}                 \tag{1}
```

The minimum possible maximum rank product is exactly 2304; a six-piece
base-field decomposition attains it. Consequently every such decomposition
has square-root rank sum at least 48, exceeding both `13^(3/2)` and the
post-half target `13^sqrt(2)`.

This is a complete finite-block obstruction, including coefficient
recombination and extension fields. It is **not** a theorem about arbitrary
lengths, joint recombinations across tensor powers, smaller selected path
sets, a different local kernel, characteristic zero, or all coherent codes.
The uniform rank problem for the construction remains open. No factoring
score changes.

## 1. Independent tensor construction

Write `C_j=[z^j]G_8(z)` using the normalized coherent crossing from the
construction note. To avoid constructing a `2^32` ambient tensor, group
coordinate pairs `(x,y)` by `(x&y,x xor y)`. If a row and a crossed column
have the same signature, put

```math
s=|x\mathbin{\mathrm{xor}}y|,\qquad
b=|x\mathbin{\mathrm{xor}}u|.
```

Their elementary tensor entry is `z^b(1-z)^(s-b)`; otherwise it is zero.
Contract each signature block with the explicit path frames. The replay
computes this formula modulo five and compares forty complete entry
polynomials with the independent magnetic-state constructor: all sixteen
four-site entries and twenty-four deterministic eight-site samples.

Let `T_j` denote the tensor with primal encoders on all four legs, so the
dual normalization has not yet been applied. Permute columns in the row
cut `xy|uv` into `vu` order, and in the column cut `xu|yv` into `vy` order.
Divide column `(i,j)` by the product of the two path norms. Both endpoint
matrices are now `I_169`. Their individual coefficient ranks agree with
the corresponding cuts of `C_j`: all changes are invertible column
permutations and separate diagonal row or column scalings.

The raw coefficient rank pairs are

```text
(169,1), (169,1), (169,21), (112,56), (33,33),
(56,112), (21,169), (1,169), (1,169).
```

Use the common tensor basis

```math
(B_0,B_1,B_2,B_3,B_4,B_5)
       =(C_0,C_2,C_3,C_5,C_4,C_6).
```

The crossing coordinates in this basis are

```math
e=(4,3,2,0,1,4).                                   \tag{2}
```

The replay checks basis independence, the full span dimension, and (2)
directly on both complete cut matrices.

## 2. Two monogenic algebras and exact ranks

For the row cut choose generator coordinates `(4,2,1,1,2,3)` in the
displayed basis. For the column cut choose `(3,4,4,2,0,2)`.
The exact minimal polynomials, with coefficients in `F_5`, are

```math
\begin{aligned}
\mu_r(X)&=X^6+3X^5+4X^2+2X
          =X(X+1)(X+2)(X+3)^2(X+4),\\
\mu_c(X)&=X^6+2X^5+4X^2+3X
          =X(X+1)(X+2)^2(X+3)(X+4).
\end{aligned}                                      \tag{3}
```

The [model certificate](../evidence/bounded-spin-rank-model.json) records
each `B_i` as a degree-at-most-five polynomial in its cut's generator.
These are explicit matrix identities, checked on all `169^2` entries.
The first six generator powers are independent, and (3) annihilates the
generator. Thus each coefficient space is the full algebra
`F_5[X]/(mu)`. No simultaneous diagonalizability is assumed.

The primary dimensions and first kernel dimensions are:

| Factor | Row primary dimension | Row first kernel | Column primary dimension | Column first kernel |
| --- | ---: | ---: | ---: | ---: |
| `X` | 20 | 20 | 20 | 20 |
| `X+1` | 1 | 1 | 35 | 35 |
| `X+2` | 43 | 43 | 70 | 57 |
| `X+3` | 70 | 57 | 1 | 1 |
| `X+4` | 35 | 35 | 43 | 43 |

Each column of primary dimensions sums to 169. The double-root component
has dimension 70 and nilpotent rank `70-57=13`, with square zero.

For coordinates `w`, form `p_w(X)=sum_i w_i p_i(X)` using the certified
basis polynomials. A simple-root component contributes its whole dimension
to rank exactly when `p_w` is nonzero at its root. The double-root
component contributes:

- 70 if `p_w(root)` is nonzero;
- 13 if `p_w(root)=0` and `p'_w(root)` is nonzero;
- zero if both vanish.

This follows by factoring `p_w` into a unit times a power of `X-root`
on that primary component. It remains valid after every field extension:
the same primary factors, matrix identities and kernel dimensions persist.
There are twelve evaluation/derivative linear forms in six coordinates,
six for each cut.

## 3. A separator valid over every extension

Enumerate the common kernels of all subsets of the twelve forms, merging
identical subspaces by reduced row echelon form. There are exactly 759
distinct kernels. On the complement of its remaining evaluation
hyperplanes, each kernel has a specified rank pair given by section 2.
For every kernel whose specified rank product is below 2304, the entire
kernel lies in

```math
\ell(w)=3w_2+3w_3=0.                               \tag{4}
```

The span of these kernels has dimension five and is exactly (4).
The checker verifies containment by applying `ell` to a basis of every
such kernel. This is a linear-algebra certificate, rather than a screen
of points over one finite field.

To see why it covers every extension `K/F_5`, take any `w in K^6` and
collect the forms vanishing at `w`. Their common kernel is one of the
759 computed spaces, extended to `K`. A further form vanishes on this
space exactly when it belongs to their base-field row span. Every form
outside that span is nonzero at `w`, since otherwise it was already
included in the collected forms. Thus the specified rank pair for this
kernel is the actual pair at `w`. If the product is below 2304, (4)
holds at `w`.

But `ell(e)=1` by (2). Linearity therefore prevents any sum of pieces
of rank product below 2304 from equaling the crossing, proving (1).

For sharpness, the model certificate supplies six vectors summing to `e`.
Their directly checked rank pairs are recorded in the
[replay receipt](../evidence/bounded-spin-rank-algebra.json); their largest
product is 2304. The same decomposition works over every extension.
The replay also exhausts all `5^6=15,625` base-field vectors through the
rank formula, and compares sixteen additional deterministic directions
with direct full-matrix elimination.

## 4. Consequence for the next construction step

The visible nilpotent cancellation is real: a seventy-dimensional sector
can have rank thirteen. At this length, however, the simultaneous cheap
directions miss one coordinate of the crossing. Optimizing the coefficient
basis further, extending the field, or assuming semisimplicity cannot
resolve that defect.

The next structural question is whether the sequence of path-compressed
cut algebras admits an effective all-length primary decomposition, and
whether the missing-crossing separator persists or disappears. A uniform
statement must prove the algebra identities and primary multiplicities;
the eight-site pattern is not such a statement. Altering the local pencil
or the selected coherent subspace is also outside this certificate.

Run `make verify-bounded-spin-algebra` (NumPy required). The standard-library
constructor controls remain `make verify-bounded-spin`. The receipt pins
the constructor, checker, and algebra certificate. This result and its
implementation used substantial Codex AI assistance; independent review
is welcome.
