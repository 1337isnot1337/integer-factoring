# Complete coefficient fibers for flat cube tensors

## Result and scope

The [coefficient engine](../tools/cube_flat_coefficient_fibers.py) includes
every homogeneous degree-six coefficient of every one of the 256 physical
words, subject only to the literal antipodal equality endpoint. It allows
odd occupation, arbitrary axis asymmetry, and coefficients in F2 or F4.
The space has **6,912 field coordinates**. The F4 implementation represents
these as 13,824 binary coordinates.

Three exact generic polynomial row frames turn proposed divisibility floors
into affine linear equations. The solver returns a complete affine fiber.
A separate polynomial multiplication checks retained assignments. Complete
generic Smith profiles are then recomputed and independently certified.

The [replay](../tools/check_cube_flat_coefficient_fibers.py) and its
[evidence](../evidence/cube-flat-coefficient-fibers.json) establish these
fixed-frame dimensions. A loss of `k` means that row `i` has floor
`max(0,s_i-k)`, where `s_i` is the corresponding reference order.

| Reference and coefficient field | Binary variables | Dimensions at losses 0, 1, 2, 3, 4 |
| --- | ---: | --- |
| Rank-two-face reference over F2 | 6,912 | 0, 0, 0, 72, 284 |
| Odd edge-face reference over F2 | 6,912 | 0, 0, 0, 0, 60 |
| Rank-two-face reference over F4 | 13,824 | 0, 0, 0, 144, 568 |

The dimensions in the last column count binary coordinates. The F4
dimensions over F4 are half those numbers. The F2 systems also establish
the corresponding ranks over every characteristic-two field, because
their coefficient and augmented matrices have entries in F2.

These are exact local search tools and scoped rigidity results. They do
not exclude tensors outside the declared row floors and do not improve
the achieved factoring score, which remains **1/2**.

## 1. A complete endpoint parameterization

A homogeneous polynomial of degree six in three variables has 28
monomials. Fix the anchor `q=x^2 y^2 z^2`. For each physical word `w` and
each of the other 27 monomials `m`, use the direction

```math
e_w(m+q).
```

Every direction vanishes at `(1,1,1)`. Conversely, evaluation at that point
is the sum of the coefficients, so these directions form a basis of the
endpoint-zero space. Starting from a valid reference therefore gives all
`256*27=6912` allowed field coefficients, with a unique coordinate vector.
No parity condition or physical support restriction is imposed.

For F4, write a coefficient as `a+b alpha`, where
`alpha^2=alpha+1`. Multiplication by `c+d alpha` sends it to

```math
(ca+db)+(da+(c+d)b)\alpha.
```

The engine expands this identity into exact binary equations. It does not
enumerate the field or specialize any polynomial coordinate.

## 2. Generic frames and affine equations

At each cut let `A(x,y)` be the full sixteen-by-sixteen physical matrix.
The third polynomial variable is set to one only after homogeneity fixes
its exponent. This retains each polynomial identity and the generic
coordinate valuation over `F4(y)`.

The [generic frame engine](cube-f4-generic-frames.md) supplies a polynomial
matrix `F` with invertible constant term and an identity

```math
F A_0=\operatorname{diag}(x^{s_i})H,
```

where `H` also has invertible constant term. Its two determinant checks
certify every reference Smith order. For a proposed tensor `A`, the
condition that row `i` of `FA` is divisible by `x^{t_i}` is affine linear
in the complete coefficient vector. The engine equates every forbidden
coefficient of every column polynomial, separately in the two field
components. Binary elimination proves consistency and gives the entire
solution space.

The fixed matrix `F` is regular and invertible at the selected coordinate,
so any solution supplies valid row-divisibility bounds. The actual Smith
profile can be better than these bounds. The search always recomputes
the actual generic profile before retaining an improvement.

## 3. Independent verification and pressure

The control replay reconstructs three complete assignments, multiplies
their three physical matrices by the reference frames, and checks all
2,304 resulting polynomial entries. It then verifies nine complete generic
profiles independently of the linear-system assembly. It also checks
serialization, endpoint validity, field coordinates, and source hashes.

For complete profiles `s_(a,i)`, the numerical search minimizes

```math
\frac{\frac13\sum_a\log\sum_i e^{-u s_{a,i}}+2u}{4\log 2},\qquad u\ge0.
```

Convexity permits a monotone derivative bisection. This numerical value is
only a search guide. The engine also exposes a separate exact rational
certificate. At `t=a/b`, set

```math
Q=\prod_{j=1}^3\sum_i a^{s_{j,i}}b^{M_j-s_{j,i}},\qquad
M_j=\max_i s_{j,i},\qquad M=\sum_j M_j.
```

The strict integer inequality

```math
(Qb^6)^B<2^{12A}(b^M a^6)^B
```

proves pressure below `A/B` for this base. The controls certify the known
rank-two reference below `679/1000`, reject an invalid `577/1000` upper
claim, and explicitly record that the accepted bound is above the cube
target. No decimal optimizer output is promoted to a theorem.

## 4. Bounded searches and reproduction

The [search driver](../tools/search_cube_flat_coefficient_fibers.py) builds
the complete equations, raises individual row floors when consistent,
samples free coordinates, and retains every strict improvement over its
reference. It can also test the fifteen positive triples of second-row
floors whose sum is seven. These are useful in conjunction with the
[low-order necessary conditions](cube-flat-low-order-floors.md).

Each record includes the complete reference tensor, source hashes,
coefficient assignments, exact profiles, budgets, and the number of
completed attempts. The replay checks retained tensors and all recorded
feasibility systems. It does not replay the random trajectory or claim
that the affine fibers were exhaustively searched.

```sh
python3 tools/check_cube_flat_coefficient_fibers.py \
  --compare evidence/cube-flat-coefficient-fibers.json
python3 tools/search_cube_flat_coefficient_fibers.py \
  --field 2 --reference rank-two --second-order --seconds 600 \
  --output /tmp/flat-cube-search.json
python3 tools/check_cube_flat_coefficient_fibers.py \
  --search /tmp/flat-cube-search.json
```

A user-supplied F4 reference is accepted through `--reference-file`.
The two deterministic [seed constructors](../tools/cube_flat_seeds.py)
provide valid starting points without temporary-file dependencies.
