# A uniform physical-matrix criterion for coupled Laurent poles

## Constructive theorem

Let `M>=4` be a power of two and let `q=2^(M/2)` be the physical
alphabet of the [reflection crossing](reflection-dual-rail-compression.md).
There is an explicit polynomial `q` by `q` matrix `B(w)`, with
`B(0)=I`, which decides the admissibility of every coupled square-zero
frame with pole order at most `M/2`.

More precisely, let `G` be any `q` by `q` matrix over a characteristic-two
field, with `G^2=0`, and let `1<=k<=M/2`. In Bell coordinates apply
`A=I+w^-k G` to legs `x,u` and `A^-T` to legs `y,v`. The complete
four-leg tensor is polynomial **if and only if**

```math
[G,B_{e+k}]+G B_{e+2k}G=0
\quad\text{for every }-2k\le e<0,                 \tag{1}
```

where `B_j=[w^j]B(w)`, `B_j=0` for negative `j`, and
`[G,C]=GC+CG`. Thus the decision uses only matrices of size `q`,
and only the first `2k` coefficients of `B`. It does not require
constructing the `q^4` entries of the crossing.

This is an all-interval result. In particular it extends the finite
[pole-two theorem](square-zero-pole-two-criterion.md): for every such
interval, pole two is admissible exactly when `GHG=0`, where `H` is
the normalized physical parity matrix defined below. The earlier
finite checkpoint remains an independent complete-coefficient check.

Every accepted frame has the original degree and strict crossing
endpoint. The result enlarges the effective construction family. It
does not give a weighted-pressure certificate or an improved
factoring exponent.

## 1. Literal block and leg conventions

Put `m=M/2`, `t=w^2`, and `u=1+w`. Define

```math
(G_M)_{ij}=t^{i\mathbin{\&}j},\qquad
(H_M)_{ij}=t^{M-1-(i\mathbin{|}j)},\qquad
(J_M)_{ij}=[i+j=M-1].
```

The one-particle source is conjugate, without changing its retained
minors, to

```math
D'(w)=
\begin{pmatrix}
uG_M&w^M I\\
w^M I&uH_M
\end{pmatrix}.                                  \tag{2}
```

Every retained physical letter selects one member of each of the
`m` reflection pairs. Denote by `rho(C)` the matrix of selected
`m` by `m` minors of an `M` by `M` matrix `C`.

The literal implementation stores tensor entries in `x,y,u,v` order.
Its selected rows in (2) are `(v,y)`, and its selected columns are
`(u,x)`. Consequently define

```math
Q(w)=\rho(uH_M),\qquad R(w)=\rho(uG_M).             \tag{3}
```

The selected-minor matrices are symmetric because `G_M,H_M` are
symmetric. The factor on the physical indices `(x,y)` is therefore
`Q`, while the factor on `(u,v)` is `R`. This order matters when
comparing the physical jets with the exact tensor implementation.

Write `Q_0=Q(0)` and define the normalized matrix

```math
B(w)=Q_0^{-1}Q(w).                               \tag{4}
```

Section 2 proves that `Q_0` is invertible and that these matrices
have the required local meaning.

## 2. Selected-minor multiplicativity over a truncated ring

Work first over

```math
\mathcal R=F[w]/(w^{2M}),
```

where `F` is any field of characteristic two. Put `U=uG_M` and
`V=uH_M`. The exact one-particle identities are

```math
UV=VU=(1+w^{2M})I,
\qquad U^T J_M U=V^T J_M V=(1+w^{2M})J_M.         \tag{5}
```

Thus `U,V` are inverse symplectic matrices over `R`.

The selected-minor representation is valid over this ring, including
its nilpotent elements. To see this without a field-only homology
argument, form the exterior complex on the `m` symplectic pairs with
differential

```math
d(a)=\left(\sum_{i=1}^{m}e_i\wedge f_i\right)\wedge a.
```

In characteristic two, `d^2=0`. For each pair the vacuum-to-double
occupation map is the identity between two free rank-one summands;
that summand is explicitly contractible. The two single-occupation
vectors form a free homology summand with zero differential.
Tensoring these decompositions over any characteristic-two
coefficient ring gives a free homology module with the `q` one-per-pair
basis vectors. Every other summand has an explicit contraction.
There is no dimension or flatness assumption here.

A symplectic matrix preserves the differential and acts on this
homology. Boundaries have zero coefficient on every admissible
one-per-pair vector. The action matrix is exactly its selected
minor matrix. It is therefore multiplicative over `R`, giving

```math
Q(w)R(w)=R(w)Q(w)=I\pmod{w^{2M}}.                \tag{6}
```

At zero, (6) proves `R(0)=Q_0^-1`. In particular, (4) satisfies
`B(0)=I`.

## 3. The third cut is an inner automorphism up to order `2M`

In a retained `M` by `M` minor of (2), there are `m` rows and
columns from each block. If a determinant permutation uses `r`
top-right entries, it must use exactly `r` bottom-left entries.
For `r>=1`, those entries contribute a factor `w^(2Mr)`.
Therefore, modulo `w^(2M)`, only block-diagonal permutations remain:

```math
T(w;x,y,u,v)=Q(w)_{xy}R(w)_{uv}\pmod{w^{2M}}.
                                                               \tag{7}
```

Apply the constant Bell normalization `Q_0^-1` on `x` and `Q_0^T`
on `v`. In the third flattening, rows are `(x,v)` and columns are
`(y,u)`, so a tensor acts on a matrix `X` by summing against `X_yu`.
Equation (7) becomes

```math
\mathcal T(w)=\operatorname{Ad}_{B(w)}
              +w^{2M}\mathcal E(w),
\qquad
\operatorname{Ad}_B(X)=BXB^{-1},                 \tag{8}
```

where `E(w)` is a matrix of formal power series. Indeed,
`R(w)Q_0=B(w)^-1 mod w^(2M)` by (6). The inverse in (8) exists as a
formal power series because `B(0)=I`. Only its finite truncation is
needed to obtain (7).

This is stronger than a first-jet approximation. Its error begins at
degree `2M`, uniformly in the interval size.

## 4. Reduce polynomiality to one physical matrix

For `A=I+w^-kG`, the hypothesis `G^2=0` gives `A^-1=A`.
The coupled physical action conjugates the third-cut operator by
`Ad_A`. Thus (8) transforms to

```math
\mathcal T'(w)=\operatorname{Ad}_{C(w)}
 +w^{2M}\operatorname{Ad}_A\mathcal E(w)
                    \operatorname{Ad}_A^{-1},
\qquad C=ABA^{-1}.                              \tag{9}
```

There are four physical factors of `A` or its inverse in the error
term. Its lowest possible exponent is `2M-4k`, which is nonnegative
when `k<=M/2`. Consequently `T'` is regular at zero if and only if
`Ad_C` is regular there.

The following elementary local fact makes that condition exact.
Let `O=F[[w]]` and suppose `C` is invertible over `F((w))`, with
`det C` a unit of `O`. If

```math
C\operatorname{End}(O^q)C^{-1}
\subseteq\operatorname{End}(O^q),                \tag{10}
```

then `C^-1 O^q` is stable under all matrix units. The diagonal units
split this lattice into fractional ideals on the coordinate axes;
the off-diagonal units force these ideals to be equal. Since `O` is
a discrete valuation ring, `C=w^h V` for an integral invertible
matrix `V`. But `val(det C)=q h=0`, so `h=0`. Hence

```math
\operatorname{Ad}_C\text{ is integral}
\quad\Longleftrightarrow\quad C\in GL_q(O)        \tag{11}
```

under the unit-determinant hypothesis.

In (9), `det C=det B`, which is a unit because `B(0)=I`. Thus
regularity is equivalent to `C` itself having no negative coefficient.
Its inverse is then automatically integral. Expanding gives

```math
C=B+w^{-k}(GB+BG)+w^{-2k}GBG.                    \tag{12}
```

The negative coefficient at exponent `e` is exactly (1).
This proves both directions of the criterion.

The original transformed tensor is a finite Laurent polynomial whose
only possible denominator is a power of `w`. Regularity at zero
therefore proves that it is a polynomial, not merely a formal series.
The identity part is the only contributor at the original highest
degree, and the paired physical maps preserve the crossing at one.

## 5. The first four pole orders

Each entry of `Q` is an `m` by `m` minor of `(1+w)H_M(w^2)`.
Its scalar factor is

```math
(1+w)^m=1+w^m.
```

Because `m` is an even power of two, all entries of `Q`, and hence
`B`, are even polynomials in `w`. Write

```math
B=I+w^2B_2+w^4B_4+w^6B_6+\cdots .
```

Let `P` be diagonal physical-letter parity and put
`H=Q_0^-1 P Q_0`. The marked-edge argument in the
[uniform first-jet theorem](square-zero-first-jet-frames.md) gives
`ad_(B_2)=ad_H`. Their difference is central, so

```math
B_2=H+cI
```

for some scalar `c`. It is unnecessary to choose this scalar:
`G(cI)G=cG^2=0`, and `[G,cI]=0`.

Substitution into (1) yields:

| Pole | Equivalent conditions | Interval range covered |
|---:|---|---|
| 1 | no additional condition beyond `G^2=0` | every `M>=4` |
| 2 | `GHG=0` | every `M>=4` |
| 3 | `[G,B_2]=0` and `G B_4 G=0` | every `M>=8` |
| 4 | `G B_2 G=0`, `G B_4 G=0`, and `[G,B_2]+G B_6 G=0` | every `M>=8` |

For pole three the commutation condition implies
`G B_2 G=B_2 G^2=0`, so the omitted quadratic condition is redundant.
The pole-four formula retains the possible cancellation between its
linear and quadratic terms. Testing them separately would reject
valid frames.

These formulas use the full `B_j` coefficients. They do not assume
that `B_4` commutes with `H`, or identify `B_4,B_6` with separately
extracted tensor jets without checking that normalization.

## 6. Effective construction and bit costs

The implementation constructs `rho(uG_M)` by a subset dynamic
program for each selected minor. A matrix entry of `G_M` is a monomial,
so multiplying it into a truncated polynomial is a shift followed
by a mask. The determinant signs coincide in characteristic two.
The factor `u^m=1+w^m` is added once per minor.

Since `H_M=J_M G_M J_M`, its selected-minor matrix is obtained by
complementing both logical bit labels. Constant binary inversion
then gives `Q_0^-1`, and (4) gives `B`. The code keeps the literal
`Q=rho(uH_M)` convention throughout.

For a cutoff of `2M`, this direct method uses
`O(q^3 M^2)` bit operations and `O(q^2 M)` bits of storage over fixed
`F2`, including elementary constant normalization. These conservative
bounds use `2^m=q` states per minor and linear bit cost for each
truncated shift or addition.

Given `G` over fixed `F4`, precomputing
`GB_j`, `B_jG`, and `GB_jG` for `j<M` costs `O(M q^3)` field
operations by ordinary matrix multiplication. Checking every
`k<=M/2` then costs `O(M^2 q^2)` further field additions. Sparse or
rank-one generators permit smaller costs, but no such improvement
is needed for effectiveness.

These costs decide admissibility. They do not construct every entry
or compute every weighted Smith profile of the transformed crossing
at unit cost. Any pressure or factoring transfer must still account
for those operations.

## 7. Independent exact controls

The checker constructs the physical matrices from selected minors
and compares their products with the independent Gaussian formula
for every tensor entry at intervals four and eight. It checks both
the original and Bell-normalized conventions, the two truncated
matrix inverses, evenness, and the scalar ambiguity in `B_2`.

Selected coordinate, nonbinary, rank-two, and rank-eight generators
are then tested in two ways: (1), using only the physical jets, and
the full four-leg Laurent expansion. Accepted and rejected poles
must agree through `M/2`.

An additional interval-16 control constructs physical matrices of
size 256, with 65,536 selected-minor entries, without constructing a
tensor of size `256^4`. It checks normalization and sampled inverse
product entries. This is a larger finite control of the uniform
constructor; the all-size theorem follows from Sections 2--4.

```sh
python3 tools/check_reflection_low_jet_frames.py \
  --include-interval16 \
  --output evidence/reflection-low-jet-frames.json
```

Files:

- [Physical jet constructor](../tools/reflection_low_jet_frames.py)
- [Independent finite checker](../tools/check_reflection_low_jet_frames.py)
- [Exact evidence](../evidence/reflection-low-jet-frames.json)

The theorem and all checks use ordinary algebra and exact finite
coefficient arithmetic. No new Lean formalization is asserted.
