# Polynomial joint-energy dressings reduce to constant conjugates

Replacing a constant joint-pair conjugation by a polynomial unimodular
matrix preserves one complete Smith profile and can change the other.
For the five-letter reference, a natural nilpotent sector admits an
all-degree reduction: if a dressing retains the original worst infinity
pole and rank-one leading coefficient, its resulting seed is already a
constant conjugate. A dressing equal to the identity at zero does not
change the seed at all.

The proof works over every coefficient field of characteristic two.
It uses a square-zero matrix sector, 180 explicit quadratic identities,
and two leading minors equal to the two remaining parameters. It is not
a stronger factoring theorem; the score remains `1/2`.

Section 6 removes the sector and rank assumptions on a constant
square-zero generator: every scalar-polynomial energy dressing
`I+f(t)P` retaining those pole/rank conditions is trivial. Its exact
quadratic elimination uses Frobenius roots over arbitrary
characteristic-two fields, rather than a finite-field census.

## 1. Polynomial unimodular conjugation

Use the notation and exact five-letter reference of
[joint energy conjugations](joint-energy-conjugations.md). In particular,
`A(t)=T(t)S`, `A(1)=I`, and the weighted infinity operator is

```math
B(z)=A(z^{-1})D_H(z)=\sum_{k=-8}^9 z^k B_k.
```

The ten-dimensional leading centralizer `L` has a common upper triangular
form on each integer energy block. Let `R` be the kernel of its six
independent diagonal characters. It has dimension four. Direct matrix
multiplication gives

```math
\mathcal R^2=0.                                        \tag{1}
```

Thus `R` is exactly the nilpotent sector of `L`, and every element is
square zero. Let `Z=R intersect Z_full`, where `Z_full` is the
four-dimensional commutant of every coefficient of the reference.
This intersection has dimension two, and its elements commute with
`A(t)` and all `B_k`.

Fix any field `K` of characteristic two. Allow a polynomial of arbitrary
degree

```math
P(t)=\sum_{j\ge0}t^j P_j,\qquad P_j\in\mathcal R\otimes K,
\qquad U(t)=I+P(t).
```

Equation (1) gives `U(t)^-1=U(t)` over `K[t]`. Define

```math
A_U(t)=U(t)A(t)U(t)^{-1},\qquad T_U(t)=A_U(t)S.           \tag{2}
```

The strict crossing endpoint survives because `A(1)=I`. Every coefficient
of `P(t)` preserves the integer energy blocks, so `U(t)` commutes with
`D_H(t)`. Both `U` and its inverse are polynomial units at zero.
Therefore (2) preserves the **entire zero-place weighted Smith profile**.
Finite products of dressings from this sector remain of the same form,
because all cross products in (1) vanish.

## 2. The quadratic terms vanish identically

For every weighted coefficient `B_k` and every pair of binary sector
basis matrices `R_i,R_j`, exact multiplication verifies

```math
R_i B_k R_i=0,\qquad
R_i B_k R_j+R_j B_k R_i=0\quad(i<j).                    \tag{3}
```

There are 180 identities. They imply `P(t)B(z)P(t)=0` over the polynomial
ring in the independent parameters `t,z` for any polynomial `P` with
coefficients in `R`. Since `P` commutes with `D_H`, the same statement
holds with `A(t)` in place of `B(z)`. Consequently the complete
transformation is linear in the polynomial dressing:

```math
A_U=A+[P,A],\qquad
B_U(z)=B(z)+[P(z^{-1}),B(z)].                           \tag{4}
```

No degree truncation or specialization to one finite field is used here.

Choose binary matrices `P_star,Q_star` completing a basis of `Z` to one
of `R`. The commutator map

```math
\mathcal R/\mathcal Z\longrightarrow\operatorname{Mat}_{25}(K),
\qquad X\longmapsto[X,B_{-7}]
```

is injective: its two displayed binary basis images are independent.
This is also a consequence of the two-coefficient stabilizer computation.

## 3. The worst pole bounds the effective degree

Coefficients in `Z` commute with the complete reference and contribute
nothing in (4), so discard all such polynomial terms for this argument.
If the remaining polynomial has degree `d>=2`, its first changed infinity
coefficient is

```math
z^{-7-d}[P_d,B_{-7}].                                  \tag{5}
```

It is nonzero by the injectivity above. The coefficient `B_-8` commutes
with all of `R`, so it creates no competing term, and all later
coefficients have larger powers. Equation (5) is strictly below `-8`.
Thus a dressing that introduces **no worse infinity pole** must have
effective degree at most one.

This pole condition is exactly that the smallest infinity Smith order
of the realigned matrix is at least the reference value `-8`: the
smallest order equals the minimum entry valuation, and realignment
only rearranges entries.

Write the remaining linear coefficient modulo `Z` as
`P_1=xi P_star+eta Q_star`. The new coefficient at `z^-8` is

```math
B_{-8}+[P_1,B_{-7}].                                  \tag{6}
```

In its realignment, two literal two-square minors are

```math
\det C[(0,3),(0,9)]=\xi,\qquad
\det C[(0,1),(0,23)]=\eta.                             \tag{7}
```

Indices use the reference cut order `x+5u` and `y+5v`. The checker builds
the symbolic entries as binary polynomials in `xi,eta` and verifies
these identities directly. Therefore leading rank at most one forces
`xi=eta=0` over **every** characteristic-two field. In fact the leading
coefficient then equals the original rank-one coefficient, so it cannot
disappear either.

It follows that every positive-degree coefficient of `P(t)` lies in
`Z`. Equation (4) now gives

```math
T_U(t)=T_{I+P_0}(t).                                   \tag{8}
```

This proves the reduction. If `U(0)=I`, then `P_0=0` and `T_U=T`.
More generally, (8) lies in the constant leading-coefficient-preserving
family from the predecessor note. Over `F4` that family has already been
checked against the reference at every scale at its natural heights.

## 4. Check and limits

The [checker](../tools/check_joint_energy_dressings.py) recomputes the
sector dimensions, all 16 square-zero products, the 180 identities (3),
the two-dimensional effective quotient and both symbolic minors (7).
Its [certificate](../evidence/joint-energy-dressings.json) retains the
quotient matrices, explicit minor indices and source hashes.

```sh
make verify-joint-energy-dressings
python3 tools/check_joint_energy_dressings.py --generate
```

The all-degree and all-field conclusions follow from (1)--(7), not from
enumerating sample dressings. This argument covers polynomial loops in
this one square-zero sector. Section 6 additionally covers arbitrary
square-zero constant energy generators with scalar polynomial multipliers.
Products leaving these classes, different heights, and dressings that permit a
worse pole or a higher leading rank remain outside it. A higher leading
rank is not itself a factoring impossibility result. No improved
coefficient-rank bound or new Lean work is claimed. The ordinary proof
and checker used substantial OpenAI assistant support and retain the
Gaussian reference's source attribution.

## 5. Every rank-one energy generator, even outside the sector

For a rank-one generator the centralizer assumption can be removed.
Let `P=v w^T` preserve the integer energy blocks, with `w^T v=0`, so
`P^2=0`. Let `f(t)` be any **nonconstant** scalar polynomial, and set
`U(t)=I+f(t)P`. Suppose its transformed infinity matrix has no entry
below order `-8` and its leading realignment has rank one. Then the
transformed seed is exactly the original reference, over every
characteristic-two field.

The additional fact about this reference is

```math
B_{-8}^4=I,\qquad (B_{-8}-I)^4=0.                      \tag{9}
```

The checker verifies the matrix identity directly. In particular its
only eigenvalue over any coefficient field is one.

Write `d=deg f>=1` and let `a_d!=0` be its leading coefficient. Before
using the sector identities, the general conjugation formula is

```math
B_U=B+f(z^{-1})[P,B]+f(z^{-1})^2PBP.
```

At power `z^(-8-d)`, the linear term is exactly
`a_d[P,B_-8]`; lower coefficients of `f` would require a nonexistent
coefficient of `B` below `-8`. Every quadratic contribution is a scalar
multiple of `P`, since `P B_k P=(w^T B_k v)P`. The no-new-pole condition
therefore forces `[P,B_-8]=beta P` for some scalar `beta`.

For a rank-one `P`, that relation implies
`B_-8 v=lambda v` and `w^T B_-8=(lambda+beta)w^T`.
Equation (9) forces both eigenvalues to be one, hence `beta=0`.
Thus `P` belongs to the leading centralizer after all. Its nilpotence
places it in the four-dimensional sector `R`.

The all-degree reduction now applies to `f(t)P`. If `P` were nonzero
modulo `Z`, a positive-degree coefficient would remain nonzero modulo
`Z`, contradicting sections 3--4. Hence `P` commutes with the entire
reference and the conjugation is trivial for every `t`.

This closes rank-one scalar-polynomial energy dressings under the stated
pole and leading-rank conditions. The next argument removes the rank-one
assumption. Transformations that relax the conditions remain open.

## 6. Arbitrary-rank scalar-polynomial generators

**Theorem.** Let `K` be any field of characteristic two. Let `P` be any
constant `25`-square matrix preserving the exact integer energy blocks
of the reference, with `P^2=0`, and let `f(t)` be any nonconstant scalar
polynomial. If the dressing `U(t)=I+f(t)P` introduces no infinity entry
below order `-8`, then `P` lies in `R tensor K`. If its leading
realignment also has rank at most one, then `P` lies in `Z tensor K`
and the transformed seed equals the reference.

The first conclusion uses only the pole condition. Neither the rank
of `P` nor its membership in the leading centralizer is assumed.

### 6.1 The top pole equations do not depend on lower coefficients

Absorb the nonzero leading coefficient of `f` into `P`, so `f` is monic
of degree `d`. Write `C_j=B_(-8+j)` for `0<=j<=17`, and `C_j=0` otherwise.
The general transformation, before any centralizer reduction, is

```math
B_U=B+f(z^{-1})[P,B]+f(z^{-1})^2PBP.                    \tag{10}
```

At powers `-8-2d+s`, `0<=s<d`, there is no linear contribution.
Characteristic two makes `f(z^-1)^2` contain only even powers.
Induction on `s` therefore forces

```math
P C_0 P=P C_1 P=\cdots=P C_{d-1}P=0.                  \tag{11}
```

Indeed the leading coefficient of the quadratic term is `P C_s P`;
every other term involves an already zero `P C_(s-2r) P`, `r>=1`.
At the next power, `-8-d`, the only linear term is `[P,C_0]` and all
lower quadratic contributions vanish by (11). Hence

```math
[P,C_0]+P C_d P=0.                                    \tag{12}
```

For `d>=18`, `C_d=0`. Thus `P` commutes with `C_0`; its square zero
places it in `R`, since the leading algebra has the common triangular
form and diagonal characters of section 1. This covers every larger
degree without a search or a bound on the coefficient field.

### 6.2 Exact quadratic elimination for degrees 1 through 17

For each remaining integer `d`, impose `P^2=0`, (11), and (12) on the
full 69-dimensional energy-block matrix space. These are equations of
degree at most two with coefficients in `F2`. For a binary basis
`E_1,...,E_s`, write `P=sum alpha_i E_i` and represent each matrix-entry
equation in the **distinct** formal monomials

```math
\alpha_i,\qquad \alpha_i\alpha_j\quad(i\le j).          \tag{13}
```

The following iteration preserves every solution over every
characteristic-two field:

1. Row-reduce the binary coefficient vectors of all equations in (13).
2. Extract their row-space intersections with the linear monomials and
   with the pure squares `alpha_i^2`.
3. Keep the linear consequences. Replace each pure-square consequence
   `sum b_i alpha_i^2=0` by its valid Frobenius root
   `sum b_i alpha_i=0`, where `b_i` are binary.
4. Parametrize the common linear kernel by a binary basis, substitute
   it into the original matrix equations, and repeat.

Binary row operations and kernel parametrizations commute with field
extension. Taking the root is valid because
`sum b_i alpha_i^2=(sum b_i alpha_i)^2`, and a field has no nonzero
square-zero scalars. **No identity `alpha_i^2=alpha_i` is imposed.**
The mixed equation `alpha+alpha^2=0` alone yields no linear constraint
in this procedure, as an explicit checker control confirms.

The complete dimension chains are:

| Degree `d` | Remaining dimensions, including initial and terminal spaces |
| --- | --- |
| 1 | `69 -> 8 -> 4` |
| 2 | `69 -> 15 -> 4` |
| 3 | `69 -> 14 -> 6 -> 4` |
| 4 | `69 -> 12 -> 4` |
| 7 | `69 -> 10 -> 4` |
| 5, 6, 8 through 17 | `69 -> 7 -> 4` |

In every case the terminal binary space is exactly `R`, and substitution
makes every original equation vanish identically. Thus this is an exact
solution-space classification for these necessary top-pole equations,
not merely a dimension count or a census of rational points.
Equations (11)--(12), the finite range `1<=d<=17`, and the separate
`d>=18` argument prove that the no-new-pole condition forces `P` into `R`
for every degree and every coefficient field in scope.

### 6.3 The leading rank closes the remaining quotient

Apply sections 2--3 to `f(t)P` once `P` lies in `R`. If `P` were nonzero
modulo `Z`, its nonconstant scalar multiplier would have a nonzero
positive-degree coefficient modulo `Z`. A degree at least two violates
the pole condition; degree one violates one of the two leading minors
(7). Hence `P` belongs to `Z`, and (10) is the original `B`.

The [Riccati checker](../tools/check_joint_energy_riccati.py) regenerates
all polynomial row spaces, their linear and pure-square intersections,
every kernel parametrization, and the final equality with `R`. The
[certificate](../evidence/joint-energy-riccati.json) records all stage
constraints, terminal bases, and source hashes. Literal dense `F4`
matrix multiplication independently checks the equation assembly at
four degrees. Those finite controls check implementation; the all-field
claim follows from the formal row-space and Frobenius argument above.
Only the Python standard library is required.

```sh
make verify-joint-energy-riccati
python3 tools/check_joint_energy_riccati.py --generate
```

For `d=1`, (12) is the invariant-graph equation for
`[[C_0,C_1],[0,C_0]]`. That viewpoint motivated the reduction, but the
certificate above starts from all 69 variables and needs no guessed
invariant subspace or external solver. The initial finite-field SAT
experiments are not dependencies of the proof or checker.

The remaining open classes include polynomial matrices with several
independent coefficient directions outside `R`, products of generators
that leave the classified sector, and gauges allowing a worse leading
pole or a larger leading rank. Closing this scalar-polynomial class
does not prove a lower bound for those classes or for factoring.
