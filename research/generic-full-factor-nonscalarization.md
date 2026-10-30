# Generic full factor pairs need not admit a common scalarization

## Precise scope

The `Mat_4` result in `commuting-halves-matrix-seed-transfer.md`
finds a normalized scalar functional that factorizes across two
full `2` by `2` tensor-factor structures. The larger proposed
auxiliary spaces behave differently.

**Proposition.** Over an algebraically closed field of characteristic
two, a nonempty Zariski-open set of changes of basis between two
full `3` by `4` tensor-factor structures on a twelve-dimensional
space admits no linear functional `ell` satisfying both
`ell(I)=1` and rank-one factorization across both factor pairs.
The same statement holds for two full `4` by `4` structures on a
sixteen-dimensional space.

The conclusion allows arbitrary linear functionals, not only
functionals of the form `phi T psi`. In particular, passing to a
further field extension cannot rescue the proposed scalarization
for a change of basis in this open set. Some such changes of basis
exist over finite extensions of `F_2`.

This proposition establishes that the common-scalar-functional
method does not collapse the entire larger catalytic interface.
It supplies no low-channel crossing identity. A particular finite
seed may still admit a scalarization on the smaller subspaces
spanned by its own factor entries.

## 1. Scalar functionals are product density matrices

The trace pairing on `Mat_b` is nondegenerate in every
characteristic: matrix units pair with their transposes. Thus
every linear functional has a unique representation

```math
\ell(T)=\operatorname{Tr}(DT).
```

Let `b=ac`. Rank-one factorization of
`ell((U tensor I)(I tensor V))` for all `U in Mat_a` and
`V in Mat_c` holds precisely when

```math
D=A\otimes B.
```

Indeed the two factors of the bilinear form are themselves trace
functionals, and pure tensors span the full matrix space. The
normalization is

```math
\operatorname{Tr}(A)\operatorname{Tr}(B)=1.        \tag{1}
```

For a second factor structure conjugated by `W`, the additional
condition is that `W^(-1)DW` also be a product density matrix.
Consequently the bad changes of basis are the projection of the
incidence set

```math
\mathcal I=\{(W,D,D'):\ W\in GL_b,\ D,D'\in\mathcal S,
                              \ WD'=DW\},         \tag{2}
```

where `S` is the variety of product matrices satisfying (1).

## 2. Dimension of the product strata

Suppose `rank A=u` and `rank B=v`, so `rank D=r=uv`.
The rank-`u` stratum of `a` by `a` matrices has dimension
`u(2a-u)`: a full-rank factorization `A=PQ` has `2au`
parameters, with a free `GL_u` change of basis of dimension `u^2`.
The corresponding dimension for `B` is `v(2c-v)`.

The nonzero product `A tensor B` determines its two factors up to
one common reciprocal scalar, and the nonzero trace normalization
cuts one further dimension. Thus the relevant product-density
stratum has dimension

```math
s_{u,v}=u(2a-u)+v(2c-v)-2.                         \tag{3}
```

For a fixed total rank `r`, put `s_r=max_(uv=r) s_(u,v)` over
the possible factor ranks. Allowing different factorizations of
the same total rank at the two ends of (2) gives dimension at most
`2s_r` for the pair of density matrices.

## 3. An elementary centralizer bound

For fixed similar matrices `D,D'`, the set of intertwiners `W`
in (2) is a coset of the invertible centralizer of `D`. Its dimension
is the dimension of the full linear centralizer.

A rank-`r` matrix of order `b` has

```math
\dim\operatorname{Cent}(D)\le (b-r)^2+r^2.         \tag{4}
```

To see this directly from Jordan form, the first column of the
conjugate zero-eigenvalue Jordan partition has `dim ker D=b-r`
boxes. All other columns, including every column for every
nonzero eigenvalue, have total size `r`. The centralizer dimension
is the sum of their squared column sizes: the space of maps between
two Jordan blocks of the same eigenvalue has dimension the smaller
block size, and blocks of distinct eigenvalues have no intertwiners.
The sum of the remaining squares is at most the square of their
sum, which proves (4).

Consequently a nonfull-rank part of the incidence (2) has dimension
at most

```math
2s_r+(b-r)^2+r^2.                                 \tag{5}
```

## 4. The twelve-dimensional case

For `a=3,c=4,b=12`, all possible nonfull ranks and bounds are:

| `r` | `s_r` | Centralizer bound | Incidence bound (5) |
|---:|---:|---:|---:|
| 1 | 10 | 122 | 142 |
| 2 | 15 | 104 | 134 |
| 3 | 18 | 90 | 126 |
| 4 | 19 | 80 | 118 |
| 6 | 21 | 72 | 114 |
| 8 | 22 | 80 | 124 |
| 9 | 22 | 90 | 134 |

For full rank, (1) implies `Tr B!=0`. A four-dimensional matrix
over characteristic two with only one eigenvalue has trace zero,
so every eigenvalue of `B` has algebraic multiplicity at most three.
Both factors are invertible. For any nonzero `lambda`, the
multiplicity of `lambda` as an eigenvalue of `A tensor B` is

```math
\sum_{\alpha\in\operatorname{spec}(A)}
 m_A(\alpha)m_B(\lambda/\alpha)\le 3\sum_\alpha m_A(\alpha)=9.
```

The centralizer is contained in the algebra of maps preserving all
generalized eigenspaces. Their dimensions sum to twelve, each is
at most nine, and hence its dimension is at most
`9^2+3^2=90`. Since `s_12=23`, the full-rank incidence has dimension
at most `2*23+90=136`.

Every stratum of (2) therefore has dimension at most 142, strictly
below `dim GL_12=144`. Its projection has closure of dimension at
most 142. The complement of that closure is the asserted nonempty
open set of changes of basis with no normalized common scalarization.

## 5. The sixteen-dimensional case

For `a=c=4,b=16`, the same calculation gives:

| `r` | `s_r` | Centralizer bound | Incidence bound (5) |
|---:|---:|---:|---:|
| 1 | 12 | 226 | 250 |
| 2 | 17 | 200 | 234 |
| 3 | 20 | 178 | 218 |
| 4 | 22 | 160 | 204 |
| 6 | 25 | 136 | 186 |
| 8 | 26 | 128 | 180 |
| 9 | 28 | 130 | 186 |
| 12 | 29 | 160 | 218 |

In full rank, each factor has nonzero trace and therefore has no
eigenvalue of multiplicity four. The product eigenvalue multiplicity
is at most twelve. The centralizer bound is
`12^2+4^2=160`; with `s_16=30`, the full-rank incidence dimension
is at most `220`. Every incidence stratum has dimension at most
250, below `dim GL_16=256`.

Finally, over the algebraic closure of `F_2`, the entries of any
chosen matrix in either nonempty open set lie in some finite
extension. The argument is an existence theorem for coefficient
structures. Producing an improving crossing tensor within such a
structure remains a separate construction problem.

## Exact arithmetic check

`check_generic_factor_nonscalarization.py` recomputes both dimension
tables from the possible factor ranks. It also checks the block and
column formulas for Jordan centralizers and their rank bound through
matrix order sixteen. The evidence is recorded in
`generic-factor-nonscalarization.json`. These finite checks verify the
displayed arithmetic; the generic-existence assertion is the incidence
argument above.
