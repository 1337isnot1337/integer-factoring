# Capacity of all-subset orthogonal coefficient systems

## Result and scope

Let `M_q(r)` be the largest possible number of nonzero columns in a rank-r
matrix C over `F_q`, where `q=2^f`. Remove redundant rows, so C has r rows.
The matrix is subject to

```math
\mathrm{rank}(C_SC_S^{\mathsf T})
=\mathrm{rank}(C_S)\quad\text{for every column subset }S. \tag{1}
```

This is the coefficient condition for the
[orthogonal paired differential construction](orthogonal-differential-pair-codes.md).
Zero columns do not count, because any number of them can be appended.
This note proves the general bound

```math
fr\le M_{2^f}(r)
\le\left\lfloor fr\frac{q^r-1}{q^r-q^{r-1}}\right\rfloor
=\left\lfloor fr\left(1+q^{-1}+\cdots+q^{1-r}\right)\right\rfloor. \tag{2}
```

It also determines two infinite families exactly:

```math
M_2(r)=\left\lfloor\frac{3r}{2}\right\rfloor,\qquad
M_{2^f}(2)=
\begin{cases}
3,&f=1,\\
5,&f=2,\\
2f,&f\ge3.
\end{cases}                                                     \tag{3}
```

Every binary example has a structural classification: after invertible row
operations and a permutation of columns, it is a direct sum of odd circuits
and single independent columns. The construction attaining the binary
maximum uses disjoint three-column circuits.

The bounds concern the coefficient interface, not arbitrary physical
retractions. They neither improve the repository's factoring exponent nor
supply a crossing-pressure certificate. Their use is to delimit the
orthogonal search and to identify the field size needed to extend it.
These are ordinary proofs; no Lean formalization is asserted.

Equivalently, every puncturing of the row code must be Euclidean LCD.
The standard LCD condition is nondegeneracy of the row space for the
ordinary dot product. Requiring it for all column restrictions is stronger
than requiring it only for the full generator matrix. The
[scaling note](tensor-orthogonal-scaling.md) gives the connection to the
usual LCD equivalence theorem and records its primary sources. This note
does not make a literature-wide priority claim for the elementary
matroid and subspace-counting arguments.

## 1. The outer-product span remembers every column subspace

Write the nonzero columns of C as `c_1,...,c_m`, and let

```math
W=\mathrm{span}_{\mathbb F_2}
\{c_i c_i^{\mathsf T}:1\le i\le m\}
\subseteq\mathrm{Sym}_r(\mathbb F_q).
```

The outer products are independent over F2. A nonempty binary dependence
would give a zero Gram matrix for a subset containing a nonzero column,
contradicting (1). In particular, `dim_F2 W=m`.

For an F_q-linear subspace U of F_q^r, define

```math
\mathrm{Sym}(U)
=\{A\in\mathrm{Sym}_r(\mathbb F_q):\mathrm{im}(A)\subseteq U\}.
```

The key identity is

```math
W\cap\mathrm{Sym}(U)
=\mathrm{span}_{\mathbb F_2}
\{c_i c_i^{\mathsf T}:c_i\in U\}.                         \tag{4}
```

Indeed, every element of W is the Gram matrix of a uniquely specified
column subset S. Its image lies inside `span_Fq {c_i:i in S}`. Equality
of ranks in (1) makes these two subspaces equal. Its image is therefore
contained in U exactly when every selected column lies in U. This proves
both inclusions in (4). Hence

```math
\dim_{\mathbb F_2}(W\cap\mathrm{Sym}(U))
=m_U:=\#\{i:c_i\in U\}.                                  \tag{5}
```

If `dim_Fq U=k`, then `dim_Fq Sym(U)=k(k+1)/2`. This remains true when U
is degenerate for the standard bilinear form. To see it, choose a full
column rank matrix A spanning U and a left inverse L of A. Every symmetric
M with image in U satisfies

```math
M=A(LML^{\mathsf T})A^{\mathsf T}.
```

Conversely, every symmetric k-by-k K gives such a matrix `AKA^T`.
The map from K is injective because multiplication by L and L transposed
recovers K. This proves the dimension statement without a nondegeneracy
assumption on U.

## 2. Counting subspaces gives a linear rank bound

Grassmann's dimension inequality in the binary space of all symmetric
r-by-r matrices, combined with (5), gives

```math
m_U\ge m-d_k,\qquad
d_k=\frac f2\bigl(r(r+1)-k(k+1)\bigr).                    \tag{6}
```

Sum (6) over all k-dimensional F_q subspaces, with `1<=k<r`. Every
nonzero column belongs to the same fraction `(q^k-1)/(q^r-1)` of these
subspaces. Dividing by their total number yields

```math
m\frac{q^k-1}{q^r-1}\ge m-d_k.
```

Consequently every k gives the bound

```math
m\le\left\lfloor
\frac f2\bigl(r(r+1)-k(k+1)\bigr)
\frac{q^r-1}{q^r-q^k}
\right\rfloor.                                           \tag{7}
```

For hyperplanes, `k=r-1` and `d_k=fr`. This gives (2). The same argument
works directly at rank one: the only hyperplane is zero, so (6) becomes
`0>=m-f`.

For completeness, the hyperplane count can be written without Gaussian
binomial coefficients. There are `(q^r-1)/(q-1)` hyperplanes, and a fixed
nonzero column lies in `(q^{r-1}-1)/(q-1)` of them. Counting pairs of a
column and a containing hyperplane gives

```math
m\frac{q^{r-1}-1}{q-1}
\ge\frac{q^r-1}{q-1}(m-fr),
```

which rearranges to (2).

This improves the earlier quadratic bound
`m<=fr(r+1)/2` and its strict variant from the
[F4 rank-two note](tensor-f4-rank-two-limit.md).

## 3. Attainment by axis bases and exact rank-two capacity

Choose an F2 basis `b_1,...,b_f` of F_q. On each of the r coordinate
axes take the f columns

```math
b_1e_j,\ldots,b_fe_j.                                   \tag{8}
```

For a column subset, the Gram matrix is diagonal. Its j-th diagonal
entry is the square of the sum of the selected basis elements on that
axis. Binary independence makes this entry nonzero exactly when the
subset contains at least one column on that axis. The Gram rank and
column rank are therefore equal. This proves the lower bound `M_q(r)>=fr`.

The upper bound equals fr whenever

```math
fr\sum_{i=1}^{r-1}q^{-i}<1.                              \tag{9}
```

A convenient sufficient condition is `q-1>=fr`; the finite sum is
strictly less than `1/(q-1)` when r is finite and positive. Rank one
always has capacity f.

At rank two, (2) is

```math
M_{2^f}(2)\le 2f+\left\lfloor\frac{2f}{2^f}\right\rfloor.
```

For `f>=3`, `2f<2^f`, and the axis construction attains 2f. At f=1,
the binary triangle `[e_1,e_2,e_1+e_2]` attains three. At f=2, the
five-column construction from the
[complete F4 catalogue](tensor-f4-rank-two-limit.md) attains five.
This proves the second formula in (3).

As a non-axis example over F8, use `a^3=a+1`, with binary polynomial
encoding. The matrix

```math
\begin{pmatrix}
1&1&2&2&4&4\\
0&1&4&6&6&2
\end{pmatrix}                                           \tag{10}
```

has six distinct projective directions and satisfies all sixty-four
subset conditions. The checker verifies this construction directly.
It realizes the uniform rank-two matroid on six columns; its full
twelve-mode paired differential kernel has dimension 1168 by the
existing matroid Hilbert formula. This dimension alone is not a
crossing-pressure bound.

Direct sums preserve (1), since both the selected column ranks and
the Gram ranks add across blocks. Thus the F4 five-column rank-two
block, together with the rank-one two-column block, also determines:

| Field | Rank | Exact capacity | Attaining blocks |
| --- | ---: | ---: | --- |
| F4 | 1 | 2 | One axis basis |
| F4 | 2 | 5 | One five-column block |
| F4 | 3 | 7 | Ranks two and one |
| F4 | 4 | 10 | Two rank-two blocks |
| F4 | 6 | 15 | Three rank-two blocks |

In each row, the block construction matches the integer upper bound
in (2). No assertion of exact F4 capacity at the omitted ranks is made.

## 4. Binary classification and exact capacity in every rank

Suppose now that q=2. A circuit means a minimally dependent set of
columns. Its kernel is the line spanned by the all-ones vector, because
every nonzero coefficient in a binary minimal dependence equals one.
For a circuit with s columns, its row space is the orthogonal complement
of that all-ones vector.

For any matrix B, the rank drop from `rank(B)` to `rank(BB^T)` is the
dimension of the intersection of its row space with its kernel. In the
circuit case this intersection is nonzero exactly when the all-ones
vector is orthogonal to itself, which occurs exactly when s is even.
Condition (1) therefore requires every circuit to have odd size.

Choose a column basis, apply an invertible row transformation, and
permute columns so that this basis is the identity matrix I_r. These
operations preserve (1). Every additional column x has even, nonzero
Hamming weight: its fundamental circuit consists of x and the basis
columns in its support, and must have odd size. In particular, every
additional column uses at least two basis coordinates.

The supports of two additional columns x and y must be disjoint.
Repeated columns would already give a forbidden two-column circuit.
If distinct x and y had a common support coordinate, consider

```math
\{x,y\}\ \cup\ \{e_j:j\in\mathrm{supp}(x+y)\}.       \tag{11}
```

The sum of the columns in (11) is zero. This dependence is minimal.
Basis columns alone are independent. A dependence using only x, or
only y, would require a basis column at a common support coordinate,
and those coordinates are absent from (11). A dependence using both
x and y must include every indicated basis column. But the size of
(11) is `wt(x+y)+2`, which is even because x and y both have even
weight. This contradicts the odd-circuit requirement.

Thus the additional columns have disjoint even supports. Each support
and its additional column form one odd circuit, while unused basis
coordinates form independent singleton blocks. This proves the claimed
classification.

Conversely, a direct sum of odd circuits and singleton blocks satisfies
(1). Every proper subset of a circuit is independent; for independent
columns the Gram rank equals the column rank. For a complete odd
circuit, the all-ones kernel vector has nonzero self-pairing, so the
Gram rank again equals the column rank. Direct sums preserve equality.

There are at most `floor(r/2)` additional columns, hence

```math
m\le r+\lfloor r/2\rfloor=\lfloor3r/2\rfloor.
```

Take `floor(r/2)` disjoint triangle blocks, and one singleton if r is
odd. This attains equality and proves the first formula in (3).

## 5. Exact replay and limits

Run

```sh
python3 tools/check_tensor_orthogonal_capacity.py
```

The [checker](../tools/check_tensor_orthogonal_capacity.py) records an
[exact receipt](../evidence/tensor-orthogonal-capacity.json). It checks
the binary classification against all coefficient sets in F2^r through
rank four, using a separate subset-Gram calculation. It also checks
the finite-field witnesses, the projective hyperplane incidence counts,
and the intersection identity (4) on explicit examples.

The mathematical upper bounds do not depend on those finite checks.
The checker distinguishes directly checked constructions from larger
constructions whose validity follows by direct sum of checked blocks.

The [deterministic extension-field scaling theorem](tensor-orthogonal-scaling.md)
still allows any starting coefficient matroid after a sufficiently large
field extension. The present capacity bounds explain why extension can
be necessary. They do not obstruct that theorem and do not address
nonorthogonal retractions.
