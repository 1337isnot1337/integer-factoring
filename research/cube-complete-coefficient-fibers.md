# Complete coefficient fibers around the binary Gaussian cube

## Results and scope

The complete even-parity polynomial coefficient space of the degree-three
binary cube has 12,672 free coefficients after its literal endpoint is
fixed. In the three fixed Gaussian row frames, lowering every reference
row floor by at most three still leaves only the original Gaussian tensor.
Lowering every floor by four leaves a 22-dimensional affine fiber.

The first statement holds over every characteristic-two field: the
coefficient matrix is binary and has full rank. The complete second fiber
was enumerated over `F2`. All `2^22 = 4,194,304` tensors have each generic
cut profile weakly dominated by the Gaussian baseline in the direction
that makes its pressure no better. This is an all-scale comparison,
but its exhaustive field is only `F2`.

These are **fixed-frame coefficient fibers**, not all tensors with nearby
Smith profiles. A different frame can satisfy different row constraints.
The broader [homogeneous binary profile floor](cube-homogeneous-profile-budget.md)
already excludes the entire unweighted common-type binary criterion.
The useful additional results here are local rigidity, a complete
coefficient-space parameterization, and exact polynomial row certificates
for generic profiles. No stronger factoring score follows.

## 1. The complete coefficient space

For a physical occupation set `U` among the eight cube vertices, require

```math
T_U\in F[x,y,z]_{3|U|},\qquad
T_U=0\quad\hbox{if }|U|\hbox{ is odd},\qquad
T_U(1,1,1)=[U=\bar U].
```

The tensor has 128 potentially nonzero physical entries. Their total
number of homogeneous monomials is

```math
1+28\binom82+70\binom{14}2+28\binom{20}2+\binom{26}2=12800.
```

Each of the 128 endpoint equations is independent. Choose
`x^|U| y^|U| z^|U|` as its anchor. For every other monomial in that entry,
the sum of the monomial and its anchor is an endpoint-zero basis vector.
Adding an arbitrary linear combination of these vectors to the
[Gaussian seed](cube-homogeneous-gaussian.md) gives every allowed tensor
exactly once. The dimension is `12800-128=12672`.

The [parameterization](../tools/cube_coefficient_fibers.py) also supports
axis and translation symmetries, but the retained rigidity computation
imposes no such symmetry. All three physical cuts are constrained
separately after permuting their variables and vertices.

## 2. Fixed row floors and exact ranks

Let `L_a,R_a` be the explicit regular Gaussian frames for cut `a`.
They include the within-side Gaussian factors and the compounds of the
one-particle monomial Smith frames. The seed becomes diagonal. Its row
orders in the actual frame order are

```text
0,2,4,6,8,10,12,14,10,12,14,16,18,20,22,24.
```

For a loss `ell`, impose on every row `i` of every cut

```math
v_a((L_aT_aR_a)_{ij})\ge\max(0,s_i-\ell)
\quad\hbox{for every column }j.
```

These conditions are linear in the coefficient variables. The program
equates every offending coefficient in the two independent remaining
coordinates, without specializing either coordinate to a field point.
Binary Gaussian elimination gives:

| Uniform row-floor loss | Coefficient rank | Fiber dimension |
| --- | ---: | ---: |
| 3 | 12,672 | 0 |
| 4 | 12,650 | 22 |

The equations are homogeneous in the perturbation from the seed. Full
rank over `F2` remains full rank over every characteristic-two extension,
which proves the uniqueness result. Likewise the loss-four solution
space has dimension 22 over every such field, although the enumeration
below covers only its binary points.

The [retained input](../evidence/cube-coefficient-floor4-input.json)
contains the 22 free coordinates, their complete lifted coefficient
assignments, and their univariate images on all three cuts. The
specialization image also has dimension 22. Reconstruction checks every
coefficient equation and every literal endpoint.

## 3. Complete binary enumeration and an all-scale comparison

For each binary assignment, set the two other cut variables to one.
Even physical parity splits its `16`-by-`16` flattening into two
`8`-by-`8` blocks. Every block has invertible endpoint at the remaining
variable one. Its determinant degree is at most 96, so arithmetic modulo
the 128th power of that variable determines all eight Smith orders.

The [native census](../tools/cube_binary_fiber_census.cpp) traverses all
assignments in Gray-code order and computes all three complete profiles.
For every profile `p`, its sorted prefixes satisfy

```math
\sum_{i\le j}p_i\le\sum_{i\le j}b_i\quad(1\le j\le16),
```

where

```text
b = 0,2,4,6,8,10,10,12,12,14,14,16,18,20,22,24.
```

Specialization can only increase determinantal valuations. Thus the
generic orders satisfy the same prefix inequalities. For `0<t<1`,
`s -> t^s` is convex and decreasing; weak majorization gives
`Z_generic(t)>=Z_b(t)` at every scale. The entry degrees and physical
binary entropy have not changed, so none of these tensors improves the
baseline pressure under this criterion.

The [complete census receipt](../evidence/cube-coefficient-floor4-census.json)
records all 4,194,304 cases, zero exceptions, and 95 complete profile
controls. A separate Python series-elimination implementation replays
all 4,560 orders in those controls. The default gate replays the basis
and retained controls; the full option repeats every native case.

## 4. Generic profiles with explicit polynomial certificates

The [bivariate engine](../tools/cube_bivariate_smith.py) supplies a
different exact interface from the existing
[extension-field specialization engine](cube-generic-smith-specialization.md).
It retains polynomial row frames as well as the invariant orders.

Homogeneity permits setting one of the other variables to one without
losing a nonzero minor coefficient. Encode a monomial `x^i y^j` by the
binary polynomial bit at `i+256j`. Each parity-block minor has degree at
most 96. Fraction-free elimination forms products of degree at most
192, so the stride 256 cannot mix the two coordinates. All divisions
are exact polynomial divisions; no finite-field evaluation occurs.

The native engine outputs a polynomial row matrix `F` and candidate
orders `s_i`. An independent Python checker multiplies the original
matrix `A` and verifies

```math
FA=\operatorname{diag}(x^{s_i})H,\qquad
\det F(0,y)\ne0,\qquad\det H(0,y)\ne0.
```

Thus both `F` and `H` are units over `F2(y)[[x]]`, proving the complete
generic profile without trusting the native pivot sequence.
Its 15 retained generic profiles cover 240 orders, 1,920 polynomial
product entries, and 60 nonzero constant-term determinants.

After restoring the homogeneous third variable, these fixed rational
frames are compatible with the
[rational-frame coefficient-transfer lemma](fixed-rational-frame-transfer.md).
The frame interface is intended for future constructions with different
alphabets, weights, or constraints; the present unweighted binary family
has the broader pressure obstruction stated above.

## 5. Reproduction

The default replay uses Python's standard library and a C++17 compiler:

```sh
python3 tools/check_cube_coefficient_fibers.py
```

To recompute every binary case as well:

```sh
python3 tools/check_cube_coefficient_fibers.py --recompute-census
```

The compiler selects carryless multiplication when available and retains
a portable polynomial-arithmetic path. Source and input hashes are
recorded in the [control receipt](../evidence/cube-coefficient-fibers-controls.json).
Elapsed census time is diagnostic, not an asymptotic complexity claim.
The proof and controls use no new Lean formalization.
