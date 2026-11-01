# Five columns are the exact F4 rank-two orthogonal limit

## Result and scope

The all-subset Gram criterion from
[orthogonal paired differential codes](orthogonal-differential-pair-codes.md)
admits five nonzero columns at rank two over F4. Five is the maximum.
An explicit five-column matrix gives an orthogonal physical kernel of
dimension 320 on ten paired modes, with smaller homogeneous truncations
available to the complete crossing constructor.

The earlier four-column example therefore did not exhaust this coefficient
interface. The new code is a construction, not an improved crossing-pressure
certificate or factoring exponent. Both crossing profiles remain separate
computational obligations. No Lean formalization is asserted.

## 1. An outer-product dimension bound

Let E be a characteristic-two field of absolute degree f, and let C have
rank r. Discard zero columns and redundant rows. Suppose every column
subset S satisfies

```math
\mathrm{rank}(C_SC_S^{\mathsf T})=\mathrm{rank}(C_S). \tag{1}
```

The symmetric outer products c_i c_i transposed are linearly independent
over F2. A nonempty binary dependence would give

```math
\sum_{i\in S}c_ic_i^{\mathsf T}=0
```

for a nonempty subset of nonzero columns. Its Gram rank would be zero
while its column rank is positive, contradicting (1). Since the symmetric
r-by-r matrices have E-dimension r(r+1)/2,

```math
m\le f\,r(r+1)/2.                                  \tag{2}
```

The condition is unchanged after choosing a row basis: writing C=UA with
U of full column rank gives Gram matrices `U G U transposed`, which have
the same rank as G. Thus r in (2) is the actual coefficient rank.

Equation (2) generalizes the
[rank-one field-size requirement](tensor-orthogonal-scaling.md#3-a-field-size-lower-bound-already-at-rank-one).
The later [capacity theorem](tensor-orthogonal-capacity.md) strengthens it
to a bound linear in r and determines rank-two capacity over every
finite characteristic-two field. It also classifies every binary example.
Zero columns are excluded from m because adding them does not affect
(1). Repeated nonzero columns are forbidden: their two-element subset
already has zero Gram matrix.

## 2. Six nonzero columns are impossible over F4 at rank two

For F4 and rank two, (2) gives m at most six. Suppose m=6. The six outer
products would be an F2 basis of the entire space of symmetric two-by-two
F4 matrices.

There are exactly fifteen nonzero symmetric rank-one matrices in this
space. They are precisely

```math
vv^{\mathsf T},\qquad v\in\mathbb F_4^2\setminus\{0\}. \tag{3}
```

The map in (3) is injective because its diagonal entries determine both
coordinates by the inverse of squaring. It covers every nonzero symmetric
rank-one matrix: a rank-one symmetric matrix is a scalar multiple of one
outer product, and every F4 scalar has a square root.

Every matrix in (3) is a subset sum of the six basis outer products.
By (1), the corresponding subset of columns must have rank one. All its
columns therefore lie in one F4 projective line.

Let t_L be the number of selected columns on a projective line L. Their
outer products lie in a two-dimensional F2 space, so independence gives
t_L at most two. Nonempty subsets of those columns can produce at most
`2^t_L - 1` rank-one matrices. Summing over projective lines gives at most

```math
\sum_L(2^{t_L}-1)\le\frac32\sum_Lt_L=9.             \tag{4}
```

Rank-one matrices on distinct projective lines have distinct ranges.
Thus (4) can account for at most nine of the fifteen matrices required
by (3), a contradiction. Consequently m is at most five.

This excludes all rank-two F4 systems with six or more nonzero columns,
including repeated projective directions. It is a bound on this coefficient
interface, not on general physical maps or factoring algorithms.

### A strict dimension bound in higher rank

The same counting argument shows that equality in (2) is impossible
whenever

```math
\frac{q^r-1}{q-1}>\frac{r(r+1)}2,\qquad q=2^f.       \tag{4a}
```

If equality in (2) held, the outer products would span every symmetric
r-by-r matrix. There are q^r-1 nonzero symmetric rank-one matrices,
again in bijection with nonzero vectors v through v v transposed.
On any projective line there can be at most f chosen columns, and their
nonempty subsets give at most 2^t-1 rank-one matrices. For 0<=t<=f,

```math
2^t-1\le\frac{2^f-1}{f}\,t.
```

Summing over lines would give at most `r(r+1)(q-1)/2` rank-one matrices,
contradicting (4a). Therefore

```math
m\le f\,r(r+1)/2-1                                  \tag{4b}
```

for every r at least three, and also for r=2 when q is at least four.
For r at least three, (4a) follows already at q=2 from
`2^r-1>r(r+1)/2`; induction proves this starting at r=3.
The rank-two binary case is a real exception: its three nonzero columns
attain the original bound m=3 and satisfy every subset condition.

## 3. An explicit five-column construction

Use omega squared equal to omega plus one, encoded by omega=2 and
omega+1=3. Take

```math
C=\begin{pmatrix}
1&1&\omega&\omega&0\\
0&1&1+\omega&1&1+\omega
\end{pmatrix}.                                     \tag{5}
```

All five projective directions occur, so any two columns are independent.
All thirty-two subset Gram conditions hold. The exact checker verifies
them directly. Together with the upper bound, this proves that five is
the exact maximum.

The complete finite catalogue enumerates all 9,948 subsets of sizes one
through six from the fifteen nonzero F4 vectors. It gives:

| Number of columns | Valid unordered sets |
| --- | ---: |
| 1 | 15 |
| 2 | 105 |
| 3 | 390 |
| 4 | 450 |
| 5 | 18 |
| 6 | 0 |

All eighteen sets of size five have one column in each projective
direction. They form one orbit under GL2(F4). The checker enumerates all
180 invertible row transformations and verifies that their orbit on (5)
is exactly those eighteen sets.

This orbit statement concerns unordered coefficient sets. Assigning
columns to different physical mode pairs can change the crossing profiles;
the unordered orbit does not identify those physical assignments. Reducing
every ordering to row echelon form gives twelve distinct ordered kernels.
Frobenius conjugation pairs them into six classes, with no fixed kernel.
The checker retains one coefficient matrix from each class. Conjugation
preserves every minor valuation, so the six classes suffice for complete
natural-profile searches at any fixed height cutoff.

## 4. The associated differential code

Use disjoint pairs

```math
(0,9),(1,8),(2,7),(3,6),(4,5),
```

each of physical height twenty. Apply the two paired differential
operators with coefficient matrix (5). The subset rank is
`r(S)=min(2,|S|)`. The
[matroid Hilbert formula](f4-differential-pair-codes.md) gives

```math
\dim V
=2^5\left(1+\frac52+\frac{\binom52+\binom53+\binom54+\binom55}{4}\right)
=320.                                               \tag{6}
```

The same formula gives every height and occupation cell. Every homogeneous
height truncation has a nondegenerate pairing, by the preceding orthogonal
criterion. Its explicit retraction is the inverse-Gram retraction on an
explicit kernel basis. Unused modes may be appended.

For example, on fourteen modes the height-twenty truncation has 54 letters,
and the height-twenty-four truncation has 84 letters. The complete-profile
driver constructs their actual maps and computes the two cuts independently.
It does not infer a pressure from the dimensions in (6).

These two cutoffs retain every coordinate state below height sixteen.
They are consequently already covered by the
[protected-core pressure floor through height twenty-five](coherent-protected-height25.md).
Their finite profiles can compare constructions, but cannot pass the
post-half natural-pressure threshold. Higher cutoffs or changes to the
protected lower states require a separate analysis.

### Completed profiles and an exact ordering comparison

Both cuts were computed separately for the explicit matrix (5). The
completed natural-profile fits are:

| Height cutoff | Physical letters | Numerical rho |
| --- | ---: | ---: |
| 20 | 54 | 0.7109029748080219 |
| 24 | 84 | 0.7103546284991576 |

The [complete record](../evidence/tensor-f4-five-pair-orthogonal-profiles.json)
retains all 19,944 directly computed Smith orders across these four cuts.
The two profiles agree at each cutoff, as required by the orthogonal
identity. These values do not improve the existing best crossing bound.

At cutoff twenty, all six Frobenius classes of ordered coefficient
kernels were also computed. Their infinity profiles were derived from
the orthogonal identity, rather than recomputed. The exact best class
at every finite natural scale is represented by

```math
\begin{pmatrix}
1&0&1&\omega&\omega\\
0&1&\omega&1&\omega
\end{pmatrix}.                                         \tag{7}
```

All classes have the same physical height multiset. Write Z_j(r) for
the zero-profile counting polynomial of class j, with (7) indexed by
zero. For each of the other five classes, the exact checker proves

```math
Z_j(r)-Z_0(r)=r^{u_j}(1-r)Q_j(r)>0\qquad(0<r<1),       \tag{8}
```

where Q_j has a nonzero, nonnegative Bernstein expansion. The exponents
u_j are respectively 60, 64, 60, 60, 60; the quotient degrees are
39, 23, 39, 31, 39. Since both cuts have the same profile, (8) also
orders their two-cut products. With equal physical height polynomials,
it gives strict natural-pressure dominance at every finite scale and
every fixed physical type.

The [comparison checker](../tools/check_tensor_f4_five_pair_dominance.py)
replays these integer identities from the
[six complete profiles](../evidence/tensor-f4-five-pair-height20-catalogue.json)
and retains the [Bernstein certificates](../evidence/tensor-f4-five-pair-dominance.json).
It does not recompute the Smith orders. The source profile run computed
17,496 zero orders directly and supplied 17,496 infinity orders through
the proved orthogonal identity.

## 5. Exact replay

Run

```sh
python3 tools/check_tensor_f4_rank_two.py
```

The [checker](../tools/check_tensor_f4_rank_two.py) and its
[receipt](../evidence/tensor-f4-rank-two-catalogue.json) retain the entire
size-five catalogue, the GL2 orbit, the weighted Hilbert cells and the
source hashes. The rank decision in the complete catalogue uses the
explicit two-by-two determinant formula. A separate general Gaussian
elimination routine rechecks all 576 subset Grams across the eighteen
largest sets. The symmetric rank-one matrices are independently enumerated.

The [profile driver](../tools/tensor_f4_five_pair_profiles.py) exposes
orthogonal and coordinate retractions and records every completed pair of
profiles. The [six-class driver](../tools/tensor_f4_five_pair_catalogue_profiles.py)
computes the zero profile of every ordered kernel class and supplies its
infinity profile by the proved orthogonal identity. It explicitly records
that the latter profile was derived rather than independently recomputed.
The numerical fits are discovery data. This construction note does not
promote them to a stronger uniform factoring bound.

This note extends the paired differential and orthogonal-retraction work
at revision `941cc1d8f56e6e80ec45938f1d660532de2f5046`. Substantial OpenAI
assistant assistance was used for the proof, implementation and text.
