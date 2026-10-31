# A constant partial-transpose intertwiner rules out tensor-power savings

## Result and scope

A local algebraic symmetry can make a favorable finite jet rank incapable
of surviving tensor powers. The symmetry is a constant invertible matrix
`P` intertwining transpose on one matched pair of physical legs. It forces
the **two cut ranks to agree for every coefficient and every collective
coefficient recombination at every tensor length**. Polynomially many such
pieces summing to a crossing must then have leading contraction
`rho >= 1`.

This statement concerns actual ranks, rather than the sufficient Smith
bound in the [multipoint floor](multipoint-smith-floor.md). It holds over
every field and physical alphabet whenever the stated constant intertwiner
exists, including one over a coefficient extension. It also excludes a
union of such piece families, since each individual piece still has equal
ranks. Linear recombinations across spaces with different intertwiners
or new physical relabelings can break that property and remain outside
this argument. It is not a factoring lower bound.

An exact application excludes **32,454,657 of 135,005,697** binary-alphabet
ternary cubic seeds with pure values at `0,2,infinity` and strict crossing
at one. These counts describe actual coefficient tensors, before quotienting
by the crossing stabilizer. Among the 200 retained Smith-profile witnesses,
66 have this certificate, including both best two-copy witnesses.
The complete factoring score remains one half. No Lean work is added.

## 1. Local transpose becomes an invertible change on the two cuts

For each fixed `x,v`, view the local polynomial tensor as a matrix

```math
F_{xv}(z)_{y,u}=K(z;x,y,u,v).
```

Suppose one constant invertible matrix `P` satisfies, coefficientwise,

```math
F_{xv}(z)^T P=P F_{xv}(z)\qquad\text{for every }x,v. \tag{1}
```

There is no symmetry assumption on `P`. Equation (1) gives
`F^T=P F P^-1`. Swapping the `y,u` tensor legs therefore equals applying
`P` to leg `y` and `P^-T` to leg `u`. In the first flattening these
are invertible row and column changes. The swapped tensor's first
flattening is the original tensor's second flattening. Thus

```math
\operatorname{rank}_{xy|uv}K
=\operatorname{rank}_{xu|yv}K.                        \tag{2}
```

The same identity holds over polynomial coefficients and under linear
combinations. Tensoring `m` positions replaces `P` by `P^(tensor m)`;
coefficient extraction, Taylor shifts, infinity reversal, and arbitrary
collective coefficient recombinations all preserve (1). Consequently
every tensor in the full coefficient span of `K(z)^(tensor m)` has the
rank equality (2). An intertwiner on the `x,v` pair gives the same
conclusion, with a transpose and fixed row/column permutations in the
cut comparison.

Ranks do not change under a field extension, so the argument also applies
when an invertible solution `P` exists only after extending the field.

## 2. Polynomially many equal-rank pieces cannot give contraction

Let `q` be the physical alphabet size, and suppose `M(m)` tensors from
this coefficient span sum to a nonzero scalar multiple of the crossing
on `q^m` labels. The crossing's first-cut rank is `q^(2m)`. Rank
subadditivity implies that some piece has first-cut rank at least
`q^(2m)/M(m)`. Its other rank is the same. Hence

```math
\max_i r_i c_i\ge\frac{q^{4m}}{M(m)^2}.             \tag{3}
```

For `M(m)=poly(m)`, the normalized log rank product has lower limit
at least one:

```math
\frac{\log\max_i(r_i c_i)}{4m\log q}
\ge1-\frac{\log M(m)}{2m\log q}=1-o(1).             \tag{4}
```

Thus no fixed `rho<1` can meet the complete per-piece width condition
through these spans. Small finite rank products can be real while the
polynomial piece count still pays the entire leading width asymptotically.

### A necessary asymmetry for every stronger crossing decomposition

The same rank-subadditivity argument gives a design condition without
assuming an intertwiner. For any `M` pieces summing to the crossing,
write `Q=max_i r_i c_i`. Some row-large piece satisfies

```math
\frac{r_i}{c_i}\ge\frac{q^{4m}}{M^2Q}.             \tag{5}
```

Some column-large piece satisfies the reciprocal inequality. Indeed
`r_i>=q^(2m)/M`, while `c_i<=Q/r_i`; substitute these two bounds.
The column argument is identical. Thus a prospective bound
`Q <= poly(m) q^(4 rho m)` requires pieces with exponential rank
imbalance in **both** directions:

```math
\frac{r_i}{c_i},\quad\frac{c_j}{r_j}
\ge\frac{q^{4(1-\rho)m}}{\operatorname{poly}(m)}.   \tag{6}
```

At label width `lambda=m log_2 q`, beating one half needs a directional
width difference exceeding
`4(1-1/sqrt(2))*lambda = 1.171572...*lambda`, up to logarithmic costs.
This does not assert that asymmetry is sufficient. It shows why a family
whose individual pieces have only subexponential rank imbalance cannot
improve the leading contraction, even if many small finite ranks occur.
The constant intertwiner is an exact, computable instance of that failure.

## 3. Exact cubic coverage and extension-field completeness

For the completed ternary cubic family, the four values
`K(0)=A,K(1)=E,K(2)=B,K_infinity(0)=C` determine every coefficient.
The crossing `E` satisfies (1) for every `P`; its `F` matrices are
scalar identity matrices. It is therefore enough to intersect the
intertwiner solution spaces of the three pure values `A,B,C`.

There are **24** projective invertible two-by-two matrices over `F_3`.
The checker tests all of them on each matched pair. This decides the
existence of an invertible intertwiner over **every extension of `F_3`**:
the linear solution space is defined over `F_3`, and the determinant of
a linear combination of its basis matrices is a polynomial of degree at
most two. A nonzero polynomial of degree at most two in each variable
cannot vanish on the whole ternary grid, by the elementary univariate
root bound and induction on the number of variables. Thus an invertible
extension-field solution implies an invertible base-field solution.

Normalize `A` using the [six pure-value stabilizer representatives](multipoint-rank-one-seeds.md).
Allow every one of the 513 possible pure/zero tensors for `B` and `C`.
The numbers excluded for the six first values are

```text
180225, 54833, 54833, 71409, 71409, 86625.
```

Their sum is **519,334** out of **1,579,014** representatives.
The first-value orbit sizes are `1,144,144,96,96,32`. Weighting the
counts by these orbit sizes gives the **32,454,657** actual coefficient
seeds excluded above. The property is preserved by the crossing
stabilizer: its change on the matched legs transforms an intertwiner
by inverse congruence, preserving its invertibility and field of definition.
The checker obtains the same counts by direct pair enumeration and by
grouping the complete bit masks of admissible intertwiners.

The remaining seeds have no such constant intertwiner. That is a necessary
escape from this obstruction, not a favorable rank bound.

## 4. Actual finite jets and independent controls

The saved actual-rank search checks one witness for each of the 200 cubic
Smith signatures, at tensor lengths one, two and three. It does **not**
exhaust all coefficient tensors within a signature. At each length it
computes both actual ranks of every jet at `0,2,infinity`, then tries
every multiplicity allocation and its exact Hermite coefficients.
Zero reconstruction coefficients are omitted correctly. The winning
allocation reconstructs every crossing entry exactly.

The smallest maximum rank products over these 200 witnesses are:

| Tensor copies | Smallest maximum rank product | Finite normalized log product |
| --- | ---: | ---: |
| 1 | 4 | 0.500000 |
| 2 | 36 | 0.646241 |
| 3 | 729 | 0.792481 |

These are finite optima over the stated witnesses and schedules, not
asymptotic rates. The two-copy values alone would be misleading.
All ranks of every certified intertwiner witness agree across the two
cuts, as the theorem predicts.

```sh
python3 tools/search_multipoint_actual_jets.py --copies 1 2 3 \
  --output /tmp/multipoint-actual-jets-replay.json
python3 tools/check_transpose_intertwiner_barrier.py \
  --output /tmp/transpose-intertwiner-replay.json
```

The independent oracle shifts each local polynomial first, forms full
polynomial Kronecker matrices, and reduces their coefficient columns.
The search instead tensors physical entry tables first, shifts the
whole tensor, and uses row reduction. Their checks cover **1,512 jet
cut-rank values** across twelve selected witnesses, including the best
finite witnesses and seeds with no intertwiner. The search also retains
all **873,600 reconstructed crossing entries** at the three lengths.

A constructive route through these coefficient spaces needs an asymmetry
that survives tensor powers and escapes (1). A small finite rank should
be evaluated only after identifying that algebraic mechanism and charging
the complete reconstruction.

This argument and its computations used substantial OpenAI Codex assistance.
