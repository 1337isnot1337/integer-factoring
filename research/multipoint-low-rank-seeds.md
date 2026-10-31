# Quadratic multipoint seeds with two-cut rank product at most two

## Result and scope

The [pure-value catalogue](multipoint-rank-one-seeds.md) extends to every
degree-at-most-two polynomial `K` over `F_3` on four binary letters with
the strict crossing endpoint `K(1)=E` and

```math
\operatorname{rank}_{xy|uv}K(a)\,
\operatorname{rank}_{xu|yv}K(a)\le2,\qquad a\in\{0,2\}.
```

For this whole family, the multipoint Hermite pressure criterion on
`0,2,infinity` requires **`rho >= 3/4`**, uniformly over all degree
allocations and both independent cut scales. The exact catalogue has
**32,265 representatives** and **321** ordered profile signatures.
Its **83** local profile pairs all receive rational prefix certificates.
This rules out an improvement on the existing `1/sqrt(2)` crossing
within the stated criterion and family. The factoring exponent remains
one half. No Lean proof is added.

The subsequent [common-cofactor proof](common-cofactor-barrier.md) also
closes this whole ternary family for actual tensor-power coefficient ranks,
with leading floor `log(3)/(2log(2))`. Its all-length conclusion uses
representation multiplicities and a complete projected-minor catalogue.

The physical alphabet is binary; the coefficient field is ternary.
The result does not cover arbitrary quadratic tensors, additional
extension-field places, larger alphabets, or higher polynomial degree.
It is a lower bound for the displayed sufficient pressure criterion,
rather than for actual tensor-power ranks or factoring.

## 1. Complete local values and nine representatives

If the first cut has rank one, write
`T(x,y,u,v)=A(x,y)B(u,v)`. The other cut has rank
`rank(A)rank(B)`. A nonzero tensor with product of cut ranks at most
two therefore has one of these forms, possibly after swapping the cuts:

- Both `A,B` have rank one: a pure four-factor tensor.
- One of `A,B` has rank one and the other has rank two.

The possibility of swapping cuts uses `y <-> u`, which preserves `E`.
Exchanging `(x,y)` with `(u,v)` also preserves `E`, and permits us to
put the rank-one matrix in the first pair. Constant physical basis
changes preserving `E` then normalize the two factors of that matrix
to `e_0(x)e_0(y)`.

For the remaining invertible two-by-two matrix `B`, the residual basis
changes give three normal forms over `F_3`:

```math
\begin{pmatrix}1&0\\0&1\end{pmatrix},\qquad
\begin{pmatrix}2&0\\0&1\end{pmatrix},\qquad
\begin{pmatrix}0&1\\1&0\end{pmatrix}.               \tag{1}
```

To see this, the residual upper triangular changes on the `x,y`
factors act by inverse transposes on `v,u`. After accounting for the
scalar introduced on `e_0(x)e_0(y)`, the entry `B_00` is invariant.
If it is nonzero, triangular shears kill the two off-diagonal entries;
the nonzero Schur complement can then be scaled to one. If `B_00=0`,
invertibility forces both off-diagonal entries to be nonzero. Independent
diagonal scalings make them one, and a remaining shear kills `B_11`.
This proves (1).

Together with the six pure/zero representatives, there are nine first
values. The second value ranges over every tensor satisfying the rank
condition: **3,585** values. The enumerator constructs these by all
projective first-cut rank-one factors and all nonzero second factors
with the required matrix-rank product, then includes the cut swap.

For an independent finite coverage check, the verifier acts on the nine
representatives with every pair from `GL_2(F_3)` and all eight physical
permutations preserving `E`. Their union equals the full 3,585-value
set exactly. Dense matrices independently check the two-cut rank product
of every value. Constant invertible leg changes preserve Smith orders;
the physical permutations can interchange the two cuts. Since the
pressure criterion treats both cuts symmetrically, this preserves its
optimized value. The 321 signatures describe the normalized catalogue;
globally swapping their cut labels need not give an additional record.

## 2. Exact all-scale exclusion and independent replay

Every pair of endpoint values defines a unique quadratic by the
Lagrange formula in the pure-value note. For each of the 321 signatures,
the same [prefix-mean certificate](multipoint-rank-one-seeds.md) gives
thresholds at the three places whose sum is at most two. At each
threshold, the chosen prefix lengths have product at least eight.
The degree-budget argument forces a contributing place to have
pressure at least `log(8)`. Consequently `rho >= 3/4` at every
allocation and scale, including schemes that omit a place.

Run:

```sh
python3 tools/search_multipoint_rank_two_seeds.py \
  --output evidence/multipoint-rank-two-seeds.json
python3 tools/check_multipoint_rank_two_seeds.py \
  --output evidence/multipoint-rank-two-validation.json
```

The checker freshly replays every representative and compares the
complete saved catalogue. It then independently checks **1,926**
block-Toeplitz cut profiles and reconstructs all **82,176 entries**
of **321 two-copy crossings**. The endpoint, infinity degree shift,
both cut orders, finite symmetry coverage, and every rational budget
are retained. No numerical minimization is used in the exclusion.

A cut-rank product of three cannot occur here: if one cut has rank one,
the other is a product of two ranks from `{1,2}`, so is `1,2,4`.
Thus this result also covers a local rank-product cutoff of three.
The [joint determinant theorem](joint-determinant-barrier.md) separately
closes pure-value cubic families over every field. Product four or higher
and new coefficient-field points remain outside this note's finite family.

This work used substantial OpenAI Codex assistance.
