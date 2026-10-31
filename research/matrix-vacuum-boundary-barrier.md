# Arbitrary matrix vacua: determinant and boundary-rank filters

The [quantum-Pascal trace construction](quantum-pascal-trace-crossings.md)
extends to two arbitrary polynomial matrix vacua. The square operator basis
still forces a first-cut determinant. A double boundary difference bounds
the zero-cut rank in terms of the two initial vacuum ranks. These are
ordinary proofs about the bare two-place Smith criterion, with exact finite
controls. They do not establish a stronger factoring exponent.

## 1. The full matrix-vacuum construction

Work over any field. Let the internal and physical dimensions both be `d`.
Use the divided lowering operators `D_t^[y]` and raising shift `X` from
the linked note. Let `P(t),Q(t)` be arbitrary polynomial `d` by `d` matrices
with `P(1)=Q(1)=E_00`, and define

```math
K(t;x,y,u,v)=\operatorname{Tr}
 \bigl(X^xP(t)D_t^{[y]}X^uQ(t)D_t^{[v]}\bigr).             \tag{1}
```

At one, `X^x E_00 D_1^[y]=E_xy`, so the endpoint is exactly
`[x=v][y=u]`. No diagonal-vacuum or charge-conservation assumption is used.

The coefficient of the matrix unit `E_rk` in `X^x P D^[y]` is

```math
[r\ge x][k\ge y]P_{r-x,k-y}(t)G_{k,y}(t).                \tag{2}
```

In lexicographic order the square basis matrix is triangular, with all
`d^2` diagonal entries equal to `P_00`, since `G_(y,y)=1`.
The trace pairing of matrix units is a permutation matrix. Hence

```math
\det K_{xy\mid uv}=\pm(P_{00}Q_{00})^{d^2}.               \tag{3}
```

Both factors are nonzero polynomials, since their endpoint values are one.
If `alpha_i` are the raw first-cut infinity orders and `n=d^2`, their sum
is `-n deg(P_00 Q_00)`. Thus, for every real `s>=1`,

```math
U(s)=\sum_i s^{-\alpha_i}
 \ge n s^{\deg(P_{00}Q_{00})}\ge n.                       \tag{4}
```

## 2. The zero-cut boundary identity

Put `A=P(0), B=Q(0)`. Since `G_(k,y)(0)=1`, the second cut is the matrix

```math
C_{(x,u),(y,v)}=
 \sum_{r=\max(x,v)}^{d-1}\sum_{k=\max(y,u)}^{d-1}
 A_{r-x,k-y}B_{k-u,r-v}.                                 \tag{5}
```

Extend entries with an index outside `0,...,d-1` by zero. Shifting both
`x,v` up by one removes the `r=d-1` boundary after reindexing. Shifting
`u,y` similarly removes the `k=d-1` boundary. Inclusion-exclusion gives

```math
C_{x,u;y,v}-C_{x+1,u;y,v+1}-C_{x,u+1;y+1,v}
 +C_{x+1,u+1;y+1,v+1}
 =A_{d-1-x,d-1-y}B_{d-1-u,d-1-v}.                         \tag{6}
```

The right side is a Kronecker product of reversed copies of `A` and `B`,
with rank `rank(A)rank(B)`. Each of the four terms on the left is obtained
from `C` by row and column shifts, so has rank at most `rank(C)`. Therefore

```math
\boxed{4\operatorname{rank}C\ge
             \operatorname{rank}A\operatorname{rank}B.}   \tag{7}
```

This holds in every characteristic, with no invertibility hypothesis.

There is a complementary exact boundary minor. Set `u=y=d-1` in (5).
Writing `a_i=A_(i,0)` and `b_i=B_(0,i)`, it becomes
`M_(x,v)=sum_r a_(r-x)b_(r-v)`. If `ord` is the first nonzero position,
or `d` for the zero vector, then

```math
\operatorname{rank}M=d-\max(\operatorname{ord}a,
                               \operatorname{ord}b).    \tag{8}
```

Indeed `M=U_a U_b^T`, where `U_a,U_b` are upper triangular Toeplitz
matrices. Factoring their first powers of the shift leaves invertible
outer factors and a middle matrix `Y^i X^j` of rank `d-max(i,j)`.
In particular, if `A_00 B_00!=0`, then `rank(C)>=d`.

## 3. Consequences and the complete affine catalogue

There are exactly `rank(C)` zero-order channels in the second cut, so
`V(s)>=rank(C)` at every `s>=1`. Combining this with (4) gives

```math
B(s)=U(s)V(s)\ge d^2\operatorname{rank}C
 \ge\frac{d^2\operatorname{rank}A\operatorname{rank}B}{4}. \tag{9}
```

For invertible `A,B` this implies `rho>=1-log(4)/(4 log(d))`, at least
`3/4` when `d>=4`. If `A_00 B_00!=0`, (8) gives `rho>=3/4` for every `d`.
These are profile restrictions, not lower bounds on factoring.

A complete binary catalogue tests a structured low-rank alternative at
`d=8`. Let `F` have rows `F(s)=(1,s_0,s_1,s_2)` for `s` in `F_2^3`.
Take `A=B=F L F^T`, for every rank-two `4` by `4` binary matrix `L`.
There are `35*35*6=7,350` such matrices: independently select two-dimensional
column and row spaces and an invertible two-dimensional linking map.
The evaluation matrix has full column rank, so every resulting `A` has rank
two, and the parametrization is injective.

The minimum exact rank of (5) is **24**, attained at coefficient mask
`L=19456` and packed vacuum mask `17361601741069221888` (row-major,
least significant bit first). Thus every polynomial path with one of
these initial pairs has

```math
B(s)\ge64\cdot24=1536,\qquad
\rho\ge\frac{\log1536}{4\log8}>\frac78,                  \tag{10}
```

where the strict comparison is the integer inequality `1536^2>8^7`.
This conclusion covers arbitrary polynomial interpolation paths with
those zero and cold endpoints, not just the linear interpolation.

## 4. Replay and next construction

[The checker](../tools/check_matrix_vacuum_barrier.py) uses prime-field
elimination, literal matrix traces, polynomial determinants, and direct
boundary differences. The optional catalogue uses NumPy 2 to accelerate
bit counts, with an independently formed literal minimum-rank matrix and
two elimination orders. [The evidence](../evidence/matrix-vacuum-barrier.json)
records the complete rank histogram and source hash.

```sh
make verify-matrix-vacuum
python3 tools/check_matrix_vacuum_barrier.py --affine-catalogue \
  --output evidence/matrix-vacuum-barrier.json
make verify
```

The standard-library target omits the optional complete catalogue.
Physical gauges, different paired types, mixed-place criteria, and arbitrary
coefficient-span recombinations are not excluded by (9). The square-basis
determinant argument also requires internal dimension equal to the physical
alphabet. A concrete next construction uses fewer physical ladder powers
in a larger internal operator space. Its endpoint remains the crossing,
while its selected operators no longer form the square triangular basis
used in (3). It still needs favorable complete profiles and a uniform cost
account. No new Lean work was used.
