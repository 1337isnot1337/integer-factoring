# Deterministic column scalings for every subset Gram matrix

## Result and scope

The [orthogonal paired-code construction](orthogonal-differential-pair-codes.md)
requires every column restriction of a coefficient matrix to have a
nondegenerate row space after scaling its columns. Its existence proof uses
a product over all column subsets and allows exhaustive search for a good
scaling. The following construction finds such a scaling in polynomial
time over an explicitly represented finite characteristic-two field.

The new step is a single curve of scalings, `lambda_i = theta^(w_i)`.
Give one column basis weight zero, and the remaining columns distinct
positive weights. A unique minimum-weight matroid basis proves that
every required Gram polynomial stays nonzero on this curve. The degree
budget is quadratic in the corank, even if the coefficient rank is large.
The
[existing finite-field specialization principle](cube-generic-smith-specialization.md)
then preserves all these polynomials simultaneously. A trace-one quadratic
tower constructs the extension and the element theta without a search for
irreducible polynomials.

In coding terminology, the requirement is that every puncturing of the
scaled row code is Euclidean LCD, meaning its intersection with its
orthogonal dual is zero. The standard theorem that every linear code
over F_q with q>3 is monomially equivalent to an LCD code addresses the
whole code. It does not supply simultaneous LCD puncturings. See Carlet,
Mesnager, Tang, Qi and Pellikaan,
[Linear codes over F_q are equivalent to LCD codes for q>3](https://doi.org/10.1109/TIT.2018.2789347),
with the [original construction preprint](https://arxiv.org/abs/1703.04346)
and [second-proof preprint](https://arxiv.org/abs/1707.08856).
The rank-one lower bound below shows why the simultaneous requirement
can force a field extension even when q>3.

This strengthens the effective construction of physical codes. It supplies
no new crossing-pressure bound and no improvement to the factoring
exponent. The complete factoring score remains one half. The field degree
here can be larger than in the preceding existential bound; the gain is
the deterministic polynomial construction.

## 1. The scaling theorem

Let F be a finite field of characteristic two, supplied with a valid
polynomial-time arithmetic representation and a computable basis over F2.
Let C be an r-by-m matrix over
F, with m positive, and let R be its rank. Choose a column basis J.
Give its R columns weight zero, and give the remaining k=m-R columns
weights 1 through k, in any fixed order. Put a=min(R,k) and set

```math
D=\sum_{j=k-a+1}^{k}j
 =\frac{a(2k-a+1)}2
 \le\frac{k(k+1)}2.                                  \tag{0}
```

Choose an extension E/F containing a nonzero element theta whose relative
degree d is strictly greater than D. If D is zero, theta=1 in F suffices.
Scale column i of C by

```math
\lambda_i=\theta^{w_i},\qquad 0\le i<m,
```

and call the resulting matrix B. Then for every subset S of columns,

```math
\boxed{\mathrm{rank}(B_SB_S^{\mathsf T})
       =\mathrm{rank}(B_S)=\mathrm{rank}(C_S).}       \tag{1}
```

The transpose in (1) is the ordinary transpose, without field conjugation.
The construction does not enumerate S.

### Unique minimum-weight bases

For a restriction to S, the zero-weight columns J intersect S are
independent. Every minimum-weight basis contains all of them. Otherwise,
adding an omitted zero-weight column creates a circuit containing a
positive-weight column: a circuit cannot lie entirely in the independent
set J. Exchanging that positive column for the omitted zero column
strictly decreases the total weight, a contradiction.

Contract the zero-weight columns in S. All remaining elements have
distinct weights, so greedy selection gives a unique minimum-weight
basis in the contracted matroid. For completeness, compare its greedy
basis G with a different basis H, and take the lowest-weight element e
where they differ. It belongs to G. The circuit formed by adding e to H
contains a heavier element f, since e is independent of the previously
retained lighter elements. Replacing f by e strictly improves H.
Consequently H cannot be minimum.

Restoring the forced zero-weight columns gives a unique minimum basis
of the original restriction. This argument permits ties among the
zero-weight columns because those columns form an independent set.

### The Gram polynomial

Fix S, put s=rank(C_S), and choose any s independent rows of C_S as a
row basis A. Define the maximal-minor polynomial

```math
Q_S(X)=\sum_{\substack{I\subseteq S\\|I|=s}}
       \det(A_I)X^{\sum_{i\in I}w_i}.              \tag{2}
```

For s=0 use Q_S=1. For positive s, the unique minimum-weight basis in the
previous paragraph supplies exactly one nonzero term at the smallest
exponent nu_S. Its coefficient cannot cancel. Therefore Q_S is nonzero.

Every exponent in (2) is a sum of at most min(s,k) distinct positive
weights from 1 through k; zero-weight columns contribute nothing.
Since s<=R, this sum is at most the sum of the largest a=min(R,k)
positive weights, which is D. Consequently

```math
Q_S(X)=X^{\nu_S}P_S(X),\qquad P_S(0)\ne0,
\qquad \deg Q_S\le D,\qquad\deg P_S\le D.          \tag{3}
```

Cauchy--Binet and characteristic two give the exact identity

```math
\det\!\left(A\,\mathrm{diag}_{i\in S}(X^{2w_i})
                 A^{\mathsf T}\right)
=\sum_I\det(A_I)^2X^{2\sum_{i\in I}w_i}
=Q_S(X)^2.                                         \tag{4}
```

Because theta has relative degree greater than D, the nonzero polynomial
P_S cannot vanish at theta. Theta is nonzero, so (3)--(4) make the scaled
row-basis Gram matrix invertible.

Write C_S=U A, where U has full column rank s. Its scaled Gram matrix
is `U G U transposed`, with G the preceding invertible s-by-s Gram matrix.
The left map U is injective and U transposed is surjective, so this product
has rank s. Column scaling also preserves rank. This proves (1), including
redundant rows, zero columns and the empty subset.

The proof applies the degree test separately to every S. It never forms
the product of their polynomials, so the number 2^m does not enter D.

### Improvement over distinct index weights

Giving every column its distinct index also works, with degree budget

```math
D_{\rm index}=\max_{0\le s\le R}s(m-s).
```

The basis-weight budget (0) never exceeds this bound. When R<=k it is
`Rk-R(R-1)/2`, compared with `D_index=Rk`. When R>k, it is `k(k+1)/2`,
while `D_index=floor(m^2/4)`. In particular:

| Rank R | Columns m | Index budget | Basis-weight budget |
| ---: | ---: | ---: | ---: |
| 2 | 4 | 4 | 3 |
| 3 | 6 | 9 | 6 |
| m-1 | m | floor(m squared / 4) | 1 |
| m | m | floor(m squared / 4) | 0 |

Every corank-one coefficient system therefore admits the construction
over a quadratic extension, regardless of its number of columns. A
full-column-rank matrix needs no extension and no scaling. At rank one,
the budget remains m-1, as required by the field-size lower bound below.

## 2. An explicit extension of controlled degree

Choose the least power of two d strictly greater than D. Then

```math
D<d\le2(D+1).
```

Write the absolute degree of F over F2 as f. The following deterministic
construction represents E with fd bits per field element.

First find a_0 in F with absolute trace one. Testing a basis of F over F2
is sufficient: the trace is a nonzero F2-linear map. Its nonzeroness also
follows from its polynomial
`X + X^2 + ... + X^(2^(f-1))`, whose degree is below the field size and
which is not the zero polynomial. Each test uses f successive squarings.

Starting with E_0=F, define

```math
E_j=E_{j-1}[z_j]/(z_j^2+z_j+a_{j-1}),\qquad
a_j=a_{j-1}z_j.                                    \tag{5}
```

If a has absolute trace one in a finite characteristic-two field, then
`X^2+X+a` is irreducible. Indeed, the image of x maps to x squared plus x
is exactly the trace-zero subspace: its kernel has two elements, its
image has half the field size, and every image has trace zero.

The relative trace of z_j is one, so the relative trace of a_j is
a_(j-1). Transitivity makes the absolute trace of every a_j equal to one.
Every quadratic in (5) is therefore irreducible. After k=log2(d) steps,
the relative extension degree is d.

For k positive choose theta=z_k. It is outside E_(k-1). Since d is a
power of two, every proper intermediate finite field over F is contained
in E_(k-1). Thus theta has full relative degree d. No multiplicative
primitivity is needed. For d=1 choose theta=1.

### Charged construction cost

The input coefficient rank and a column basis use ordinary Gaussian
elimination over F.
The trace-one search costs polynomially many operations in f. The tower
has O(log(d)) levels, and all of its defining coefficients occupy O(fd)
bits in total. Direct multiplication and inversion in the tower have
polynomial bit cost in f and d. The positive-weight scales are the first
k powers of theta, and the other scales equal one. They need only k
additional multiplications.

Since d=O(D+1)=O((m-R)^2+1), this is a deterministic polynomial-bit-time construction
in r, m and f, under the stated base-field representation. It produces
all scales explicitly and performs no column-subset search. The field
elements occupy O(f(D+1)) bits. The earlier O(f(Rm+1)) upper bound also
holds, since the new degree budget never exceeds the index budget.

This does not make construction of the entire physical kernel polynomial
in the number of modes. That kernel can have exponentially many occupation
coordinates. Its basis, retraction and complete crossing profiles still
incur their full output and computation costs.

## 3. A field-size lower bound already at rank one

Take the one-row matrix C=(1,...,1) with m columns. Any scaling satisfying
(1) over a characteristic-two field E must obey

```math
\sum_{i\in S}\lambda_i^2
=\left(\sum_{i\in S}\lambda_i\right)^2\ne0
\quad\hbox{for every nonempty }S.
```

Hence the m scalars are linearly independent over F2. In particular,

```math
[E:\mathbb F_2]\ge m,\qquad |E|\ge2^m.              \tag{6}
```

An extension of fixed cardinality cannot solve this interface for growing
m, even at rank one. Over a fixed base of absolute degree f, the preceding
construction has absolute degree less than 2fm in this rank-one case,
which matches (6) within a factor depending only on the fixed base.

## 4. Use in paired differential codes

Apply (1) to the coefficient matrix of the disjoint paired derivatives
in [the preceding construction](orthogonal-differential-pair-codes.md).
The resulting kernel has nondegenerate occupation pairing at every height
truncation. Given its explicitly constructed basis P, the retraction

```math
L=(P^{\mathsf T}P)^{-1}P^{\mathsf T}
```

is homogeneous and deletion closed, and LP=I. The complete natural zero
and infinity profiles coincide by the prior constant-Gram change-of-basis
proof. Their actual values and pressure remain separate obligations.

For a fixed finite code the enlarged coefficient field is fixed before
word-power amplification. If the code grows with the factoring input,
the fd-bit coefficient size and all field arithmetic must be included
in the evaluator cost. This note does not infer an improved exponent
from existence of these maps.

## 5. Exact implementation and replay

The [implementation](../tools/tensor_orthogonal_scaling.py) represents the
explicit quadratic tower over F2, including F4 and further tower bases.
The theorem permits any represented finite characteristic-two base; the
current implementation exposes this specific tower family. It constructs
the scales without enumerating subsets.

The [checker](../tools/check_tensor_orthogonal_scaling.py) independently
enumerates maximal minors and forms the formal Gram determinant. Its
[receipt](../evidence/tensor-orthogonal-scaling-controls.json) records:

- All 256 binary two-by-four coefficient matrices, including zero columns
  and dependent rows.
- Thirty F4 projective-column cases and twenty additional F4 matrices,
  for 306 coefficient matrices in total.
- 5,176 complete column-subset Gram checks, independent formal
  Cauchy--Binet identities and unique-minimum-basis checks.
- 596 comparisons with a separate four-product tower multiplier, plus
  trace-one, inverse and Frobenius controls through absolute degree 16.
- Rank-one controls through nine columns, checking independence over F2.
- Twenty-one corank-one F4 matrices through rank seven, using only a
  quadratic extension, with 1,524 additional direct subset-Gram checks.
- Seven full-column-rank controls, using no extension and identity scales.

The largest extension used by the matrix checks has absolute degree 16.
All arithmetic in these checks is exact. Run

```sh
python3 tools/check_tensor_orthogonal_scaling.py
```

The code and theorem are based on repository revision
`941cc1d8f56e6e80ec45938f1d660532de2f5046`. The general specialization
principle and the physical-code criterion are credited above. This note
does not claim a new general matroid theorem or independent mathematical
acceptance. Substantial OpenAI assistant assistance was used for the
derivation, implementation and text.
