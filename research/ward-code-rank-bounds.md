# Computable Ward bounds and actual overlapping codes

The [graded Ward theorem](noncommutative-crossing-ward-rank-descent.md)
factors restricted coefficients through full Pfaffian matrices with smaller
row grades. The implementation here evaluates those rectangular matrices
using the constant stable decomposition, without constructing their full
ambient cut spaces.

The recorded upper bounds do not establish a better factoring exponent.
They are upper bounds for particular coefficient decompositions; exceeding
the target is not an obstruction for the underlying codes. Separately,
complete binary searches exclude the two finite overlapping-form codes
recorded below.

## 1. Rectangular stable blocks

Write a full rectangular cut as

\[
M_P(R_A,R_Y;C_V,C_U)[a,b]
\]

at `p=1+x,q=1+y,r=1`. The stable identity with a constant `G` ancilla
decomposes it into the same Frobenius blocks used for equal grades. The
four ancilla states have grades

\[
(1,0),\quad(0,1),\quad(2,1),\quad(1,2).
\]

For an ancilla row of grades `(a_r,b_r)`, a column of grades `(a_c,b_c)`,
and a stable summand with Frobenius weights `(m,n)` and determinant twist
`d`, the selected weights are

\[
\begin{aligned}
u_r&=(R_Y+b_r-R_A-a_r+m)/2,\\
u_c&=(C_U+b_c-C_V-a_c+m)/2,\\
v_r&=(R_A+R_Y+a_r+b_r-2P-2+n)/2,\\
v_c&=(C_V+C_U+a_c+b_c-2P-2+n)/2.
\end{aligned}
\]

A nonintegral or unavailable weight gives the zero block. For selected
binary masks with weighted intersections `w,z`, the matrix entry is

\[
x^d(1+x)^{u_c-w}(1+y)^{v_r-z}(1+x+y)^z.
\]

The code extracts its coefficients exactly in characteristic two, computes
each explicit block rank, and sums with the exact stable multiplicities.
It minimizes over all six constant entries of the ancilla. It also caps
the rank by the literal rectangular matrix dimensions. Equal-grade
blocks are included as a special case.

The Ward row formula needs no channel multiplicity. For the column
formula, an `ell`-mode channel is nonzero only if the corresponding
weight columns are linearly independent. Thus the generic factor
`2^ell binom(P,ell)` may be replaced by `2^ell` times the number of
independent `ell`-subsets of weight columns. For `k` disjoint forms on
blocks of `b` pairs each, this is exactly

\[
H_\ell=2^\ell\binom{k}{\ell}b^\ell.
\]

The final rank bound takes the minimum of the Ward bound, the original
full-grade bound, and `d^2` in each cut. The reported cost rounds every
square root upward to an integer. This makes the cost a certified upper
bound, without a floating-point comparison in its construction.

## 2. Label dimensions and recorded bounds

For a block of `b` pairs, the homogeneous radical dimensions are

\[
R_b(L)=\sum_s\binom bs2^s
       \binom{b-s-2}{(L-s)/2-1},
\]

where a summand is included only when `b-s` is positive and even and
all binomial indices are valid integers. The product-image dimensions
for disjoint blocks are coefficients of the product of these generating
polynomials. The orientation-half construction removes the zero-single
syndrome and divides the remaining dimension by two. The full code can
also be used; it need not cancel the second endpoint.

The reproducible screen includes both choices. For the full single-form
code the recorded examples include:

| Pairs `P` | Grade `L` | Labels | Computed `rho` |
|---:|---:|---:|---:|
| 16 | 10 | 8,811,817 | 0.933660 |
| 24 | 14 | 57,671,029,749 | 0.881961 |
| 32 | 18 | 402,083,973,823,997 | 0.854429 |

Here `rho=log(S)/(2 log(d))` for the recorded integer upper bound.
The target is `rho<1/sqrt(2)`. The largest explicit matrix in the last
calculation has only 80 rows or columns. These examples illustrate the
size reduction, but none certifies the target. The code also records
orientation halves and disjoint two-form and three-form examples.

The evidence covers exactly the cases listed in the two JSONL files.
It is not a complete search over grades, form matrices, or output maps.

## 3. Constructing overlapping forms literally

For binary weight columns `a_i`, the constructor works syndrome by
syndrome. On the even positions it checks independence and zero Gram
matrix. It then solves the simultaneous creation and contraction kernel
equations at the requested homogeneous grade. Their dimension is checked
against the all-size Ward dimension formula.

The constructor expands the resulting vectors in the original exterior
basis and chooses a coordinate left inverse by row reduction. It keeps
the orientation in which the least single pair uses its low mode. Thus
`Pi I=Id`, `I^T I=0`, and `Pi J I=0`. Both pure input vertices vanish,
while the remaining vertex is the crossing tensor.

These are actual tensor constructions, with every contracted entry and
both complete cut ranks computed. They are not dimension-only models.

| Weight-column types | Grade | Labels | Coefficient dimension | Common factor | Binary cost lower bound | Target |
|---|---:|---:|---:|---|---:|---:|
| `1,1,2,2,3,3` | 6 | 6 | 2 | `(q+r)^4(p+r)` | 36 | 12.603 |
| `1,1,1,2,2,3,3` | 7 | 34 | 13 | `(q+r)^2` | 1000 | 146.501 |

The types `1,2,3` encode the columns `(1,0),(0,1),(1,1)`. Both examples
therefore have two genuinely overlapping forms. The first has stronger
cancellation than the uniform Ward theorem requires. The second retains
13 independent coefficient directions.

All `3` and `8191` nonzero binary combinations are included. The rigorous
lower bounds come from shortest paths with each edge cost rounded down.
The second example has a concrete decomposition of cost approximately
`1001.268423`, using coefficient masks `494` and `495`. Their XOR is the
crossing mask. The complete rank tables, the centered jet ranks, and
the Ward bounds for those same jets are saved in the evidence.

These finite exclusions apply over the binary field to these two output
maps. They do not exclude other output maps or extension fields.

## 4. Independent checks and reproduction

The rectangle audit compares 18 complete small literal rectangle profiles
with the stable upper bounds. The Ward descent checker independently
constructs four restricted-code controls. A separate checker compares
every nonzero centered jet of each of the two overlapping-form examples
against both computed Ward rank bounds.

From the repository root:

```sh
python tools/check_noncommutative_crossing_ward_rank_descent.py
python tools/check_noncommutative_crossing_ward_rectangles.py
python tools/check_overlapping_symplectic_profiles.py --large \
  --output /tmp/overlapping-symplectic-code-profiles.json
python tools/screen_ward_code_bounds.py \
  --output /tmp/ward-code-rank-bounds.jsonl
python tools/screen_ward_code_bounds.py --full-code \
  --case 16:1:10 --case 24:1:14 --case 32:1:18 \
  --output /tmp/ward-full-code-rank-bounds.jsonl
```

The `--large` reconstruction includes the full 34-label contraction and
binary search. The default overlapping-code checker rebuilds the smaller
example, which is also part of the regular verification target.

Sources and evidence:

- [Rectangular bound engine](../tools/ward_rectangular_bounds.py).
- [Literal rectangle audit](../tools/check_noncommutative_crossing_ward_rectangles.py)
  and [audit evidence](../evidence/noncommutative-crossing-ward-rectangles.json).
- [Code-bound screen](../tools/screen_ward_code_bounds.py),
  [orientation-half bounds](../evidence/ward-code-rank-bounds.jsonl), and
  [full-code bounds](../evidence/ward-full-code-rank-bounds.jsonl).
- [Overlapping-code constructor](../tools/explore_overlapping_symplectic_codes.py),
  [jet and rank checker](../tools/check_overlapping_symplectic_profiles.py), and
  [complete finite evidence](../evidence/overlapping-symplectic-code-profiles.json).
