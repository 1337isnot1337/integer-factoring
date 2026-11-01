# Exact generic polynomial frames over the four-element field

## Interface and scope

The [F4 profile engine](../tools/cube_f4_bivariate_smith.py) computes
complete generic Smith orders and explicit polynomial row frames for
matrices of order at most sixteen with bivariate entry degree at most
six. It also accepts eight-leg binary-physical tensors whose entries
are homogeneous degree-six polynomials in three variables. Physical
coefficients may have either occupation parity.

This extends the earlier [binary parity-block engine](cube-complete-coefficient-fibers.md)
to full sixteen-by-sixteen cuts and nonbinary coefficients. It is a
construction and verification tool, not an improved factoring bound.
The [flat parity floor](cube-flat-parity-floor.md) makes the full physical
matrix relevant: a search confined to even coefficients cannot cross
that criterion's threshold.

Write `F4=F2[alpha]/(alpha^2+alpha+1)`. Coefficients `0,1,2,3` mean
`0,1,alpha,1+alpha`. A bivariate polynomial is represented by two binary
polynomials, its constant and alpha parts. The binary bit at `i+256j`
represents the monomial `x^i y^j`.

The public functions are:

- `profile_matrix(matrix)`: generic orders and polynomial row frames.
- `verify_matrix(matrix, record)`: an independent exact certificate check.
- `profile_tensor(tensor, axis)`: the full physical cut on the given axis.
- `tensor_matrix(tensor, axis)`: its bivariate polynomial matrix.
- `terms(value)`: the nonzero triples `(x_degree,y_degree,F4_coefficient)`.

A tensor entry is a dictionary from exponent triples to field coefficients.
Binary monomial sets are also accepted. The tensor wrapper checks every
entry degree and the literal antipodal endpoint by default.

## 1. Why the generic computation is exact

Every `k`-minor of a degree-six homogeneous tensor is homogeneous of
degree `6k`. For one cut variable `x`, write a minor as

```math
\Delta(x,y,z)=\sum_i x^i f_i(y,z).
```

Each nonzero `f_i` is homogeneous of a known degree. Setting `z=1`
is injective on that homogeneous space, so the leading coefficient
in `x` cannot disappear. This preserves every determinantal valuation
and hence every generic Smith order.

The [native implementation](../tools/cube_f4_bivariate_smith.cpp) uses
fraction-free elimination on the resulting polynomials. Each minor has
total degree at most `16*6=96`; the products before exact division have
degree at most 192. The packing stride 256 therefore separates both
coordinates throughout the calculation. All divisions are checked for
zero remainder over `F4`; no coordinate is evaluated at a finite-field
point.

The computation returns a polynomial row matrix `F` and ascending
orders `s_i`. The Python verifier multiplies the original matrix `A`
independently and checks

```math
FA=\operatorname{diag}(x^{s_i})H,\qquad
\det F(0,y)\ne0,\qquad\det H(0,y)\ne0.                \tag{1}
```

Both matrices are therefore invertible over the local ring
`F4(y)[[x]]`. Equation (1) proves the complete profile without trusting
the native pivot decisions. The verifier also checks dimensions,
degree bounds, sorting, pivot permutations, and every reported prefix.

The right inverse of `H` is a fixed rational matrix regular in `x`.
Consequently these certificates fit the
[fixed rational-frame transfer](fixed-rational-frame-transfer.md).
Different physical cuts may use different coefficient expansion cones,
as specified in that transfer theorem.

## 2. Use in coefficient fibers

A flat degree-six tensor has 28 monomials per physical word. Prescribing
one endpoint equation for each of the 256 words leaves 6,912 free
`F4` coefficients, or 13,824 binary scalar coordinates. For a fixed row
frame, divisibility conditions on `F A` are linear in those coefficients.
Each `F4` equation splits exactly into two binary equations.

For a binary reference and binary frames, extending the field alone
does not change feasibility of a prescribed row-floor target: the
inhomogeneous binary system and its homogeneous companion separate.
Nonbinary samples can nevertheless yield new generic frames for a
subsequent search. The engine provides those frames and verifies their
regularity. No successful target crossing is claimed here.

## 3. Independent controls

The [checker](../tools/check_cube_f4_bivariate_smith.py) covers sixteen
complete profiles, 2,913 independently multiplied entries of `F A`,
32 nonzero constant-term determinants, and five rejection controls.
It includes all three cuts of the binary flat seed, compared with the
earlier binary engine; nonbinary constant-unit transformations; six
independently prescribed diagonal profiles; and literal endpoint-preserving
parity-breaking tensors over `F4`.

One generic polynomial matrix is singular at every specialization
`y in F4` while remaining nonsingular over `F4(y)`. This control shows
why enumerating the base field cannot replace the generic certificate.
Other controls reject altered orders, corrupted frames, and invalid
input bounds.

```sh
python3 tools/check_cube_f4_bivariate_smith.py \
  --output evidence/cube-f4-bivariate-smith.json
```

The [receipt](../evidence/cube-f4-bivariate-smith.json) records the exact
checks and source hashes. A C++17 compiler and Python's standard library
are sufficient. The native implementation retains both portable and
carryless-multiplication paths. No new Lean formalization is used.
