# Nonabelian group shears give polynomial swap paths

Three group multiplication/division shears and a local inversion implement
the swap over every group. For a finite `p`-group in characteristic `p`,
interpolating the shears gives a determinant-one path up to a fixed sign.
This extends the construction beyond the involutory Steiner setting.
Six exact binary controls fail the bare criterion; no all-group obstruction
or improved factoring exponent is claimed.

## 1. A swap identity without commutativity

Let `G` be a finite group of size `q`. Define permutations of `G^2` by

```math
A(x,y)=(xy,y),\qquad B(a,y)=(a,ya^{-1}),\qquad
C(a,b)=(ba,b),\qquad D(c,b)=(c,b^{-1}).                    \tag{1}
```

The identity

```math
(x,y)\xrightarrow A(xy,y)\xrightarrow B(xy,x^{-1})
 \xrightarrow C(y,x^{-1})\xrightarrow D(y,x)               \tag{2}
```

uses associativity and inverses, with no commutativity assumption. Thus
`DCBA=S`, the swap. The final permutation `D` acts on just one coordinate.

Now suppose `q=p^h` and work over any field of characteristic `p`. Each
of `A,B,C` is a controlled group translation. Its order divides `q`, so
`(P-I)^q=P^q-I=0`. Consequently `W_P(t)=I+t(P-I)` is invertible over the
polynomial ring, with a finite geometric inverse, and has determinant one.
The cubic path

```math
K(t)=D W_C(t)W_B(t)W_A(t)                                 \tag{3}
```

has constant nonzero determinant `det D`, starts at the local permutation
`D`, and ends at the crossing `S`. Its second cut at zero has rank one.
In general the inverse path has higher degree; determinant constancy does
not imply a symmetric infinity profile or a degree-three inverse.

The [polynomial crossing transfer](polynomial-crossing-interpolation.md)
applies to every fixed finite group table in this family. At tensor length
`m` there are at most `3m+1` coefficient positions. The finite seed table,
polynomial coefficients, and rank factorizations retain the existing
`2^O(lambda)` construction cost. A favorable complete rank or profile
bound would still be needed for a better factoring machine.

## 2. Small exact certificates exclude the displayed binary controls

For characteristic two, let `r` be the first-cut rank of the degree-three
coefficient, and let `b_2` be the number of second-cut Smith orders equal
to two at zero. Since the complete polynomial degree is at most three,
there are exactly `r` raw infinity channels of order `-3`. Thus for every
real scale `s>=1`,

```math
U(s)\ge r s^3,\qquad V(s)\ge b_2s^{-2},\qquad
U(s)V(s)\ge r b_2.                                      \tag{4}
```

The number `b_2` is recovered from just three exact Toeplitz ranks. If
`R_j=sum_i max(0,j-nu_i)`, then

```math
b_0=R_1,\qquad b_1=R_2-2R_1,\qquad
b_2=R_3-2R_2+R_1.                                      \tag{5}
```

The following are complete first-jet certificates for the displayed seeds,
not selected numerical scales:

| Group | `q` | `r` | `b_2` | All-scale floor `r b_2` | `q^3` |
|---|---:|---:|---:|---:|---:|
| Cyclic `C_8` | 8 | 34 | 46 | 1,564 | 512 |
| Dihedral `D_8` | 8 | 24 | 39 | 936 | 512 |
| Quaternion `Q_8` | 8 | 26 | 40 | 1,040 | 512 |
| Cyclic `C_16` | 16 | 174 | 200 | 34,800 | 4,096 |
| Dihedral `D_16` | 16 | 122 | 191 | 23,302 | 4,096 |
| Dihedral `D_32` | 32 | 582 | 855 | 497,610 | 32,768 |

Here the subscript on the dihedral group denotes its order. In each row
`r b_2>q^3`, giving bare parameter strictly greater than `3/4`. This exact
comparison is sufficient to reject these seeds without optimizing a real
scale or computing their complete Smith profiles.

## 3. Replay and remaining scope

[The checker](../tools/check_group_shear_crossings.py) constructs literal
group laws, checks associativity, verifies the swap identity and the
`p`-power orders of the shears, expands every polynomial endpoint entry,
and computes the leading and Toeplitz ranks with two elimination orders.
[The evidence](../evidence/group-shear-crossings.json) includes the three
Toeplitz ranks for every row and a source hash. No optional numerical or
native package is needed for these certificates.

```sh
make verify-group-shears
make verify
```

Unlike the [Steiner cubic](steiner-cubic-crossings.md), nonabelian group
shears need not generate six-dimensional permutation orbits, and the
finite certificates above are not an all-size theorem. They also do not
exclude different polynomial shear schedules, physical word types,
independent local gauges, mixed-place transfers, or collective coefficient
recombinations. The current complete factoring exponent remains `1/2`;
no new Lean work was used.

The next constructive question is to change the shear schedule or the
coefficient representation so that many second-cut channels appear after
the leading first-cut pole order. A single zero channel is insufficient:
the order-two bulk in the table already rules out these cubics.
