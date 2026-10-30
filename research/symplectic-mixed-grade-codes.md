# Mixed exterior grades: endpoint contacts and exact repairs

Taking a direct sum of homogeneous radical codes increases the label
dimension, but it also reintroduces Pfaffian contact terms. The coordinate
identity `ΠI=Id` alone does not make the unrestricted mixed-grade endpoint
equal to the crossing. This note gives an explicit counterexample and
checks two repaired constructions.

The test uses four complementary mode pairs and the direct sum of the
radical spaces in grades 2, 4, and 6. Their dimensions are 1, 26, and 1, so the total
label dimension is 28. These dimensions follow from the
[symplectic radical formula](noncommutative-crossing-symplectic-radicals.md).
Odd grades have zero radical for this parity of the number of pairs.

## 1. The raw endpoint fails

Let `D_i=e_i f_i` be a filled mode pair. The first two code labels have
input vectors and coordinate output selectors
\[
 I_0=D_0+D_1+D_2+D_3,\qquad \Pi_0=[D_3],
\]
\[
 I_1=D_0D_1+D_0D_3+D_1D_2+D_2D_3,
 \qquad \Pi_1=[D_2D_3].
\]
At `p=q=0,r=1`, the literal Pfaffian contraction gives
\[
 G_{0,0,1,1}=1.
\]
The crossing at the same entry is
`δ_{0,1}δ_{0,1}=0`. Here the grades of `x,y` are 2 and those of `u,v`
are 4. The unused modes can participate in additional contractions.
The fixed-grade containment argument that forced every matching to the
crossing no longer applies.

The checker evaluates all 77,056 parity-compatible entries of the
28-label tensor. Exactly 604 constant-term entries differ from the
crossing. Every parity-incompatible entry vanishes by the same exact
local syndrome conservation used in the homogeneous checker.

Thus neither isotropy of the input code nor membership in
`ker Ω ∩ ker Λ` removes all mixed-grade endpoint contacts.

## 2. Two honest crossing constructions

Throughout the computation, set `r=1` and retain the complete polynomial
in `p,q`. Different output grades can have different total degrees.

### 2.1 Constant row-matrix normalization

Let `M(p,q)` be the row-cut matrix of the raw mixed-grade tensor, and
let `M₀=M(0,0)`. In this finite example, `M₀` is invertible in every
syndrome block. Define
\[
 \widetilde M(p,q)=M_0^{-1}M(p,q).
\]
Then `\widetilde M(0,0)=Id`, so the normalized tensor has the exact
crossing endpoint. The checker computes the inverse by binary
elimination and verifies this identity entry by entry.

This transformation preserves each row-cut coefficient rank. It can
change the other cut substantially. The checker therefore reshuffles
every normalized coefficient and recomputes all column-cut ranks. No
column-rank invariance is assumed.

The normalized family has 15 independent polynomial coefficients. The
complete binary screen tests all 32,767 nonzero combinations. Its 93
individually affordable vectors span only 11 dimensions and do not
contain the crossing. Every affordable vector has both ranks checked
again by the Python implementation.

### 2.2 Balancing the grades in each crossed pair

Alternatively, define
\[
 G^{\rm bal}_{xyuv}
 =\mathbf1_{|x|=|v|}\,\mathbf1_{|y|=|u|}\,G_{xyuv}.
\]
The filter retains arbitrary operator and state grades, but requires
equality within each crossed pair. At the `r` endpoint, the original
containment argument now applies separately to those pairs. The checker
also verifies the resulting crossing directly.

The balanced family has 9 independent coefficients. Its complete binary
screen has 15 individually affordable vectors, spanning 6 dimensions and
missing the crossing. A complete binary shortest-path calculation gives
a decomposition of cost 524.64671285 and a rigorous integer lower bound
of 523. The finite seed threshold is `28^{√2}`, which is below 112.

The complete four-element-field screen tests 87,381 projective vectors.
Its 37 affordable vectors again span only 6 dimensions and miss the
crossing. A separate scalar implementation checks all 37 vectors and
forty further deterministic controls.

| Construction | Field | Coefficient dimension | Vectors tested | Affordable vectors | Affordable span dimension |
|---|---|---:|---:|---:|---:|
| Row-normalized | `F₂` | 15 | 32,767 | 93 | 11 |
| Grade-balanced | `F₂` | 9 | 511 | 15 | 6 |
| Grade-balanced | `F₄` | 9 | 87,381 | 37 | 6 |

All three span tests exclude a winning decomposition in their stated
finite coefficient spaces. No four-element-field exclusion is claimed
for the 15-dimensional row-normalized family.

## 3. Reproduction and scope

```sh
c++ -O3 -std=c++20 tools/screen_pfaffian_span_f4.cpp -o /tmp/omega-span-screen
python3 tools/check_symplectic_mixed_grade_codes.py --native /tmp/omega-span-screen --field4 --export-dir /tmp/omega-mixed-grades --output /tmp/omega-mixed-grades.json
```

The checker has separate `raw`, `normalize`, and `balanced` modes. It
retains the input code, coordinate output choices, endpoint failures,
complete coefficient profiles, independent scalar coefficient
polynomials, native screen hashes, and all eligible-vector audits.

The files are:

- [Literal construction and checker](../tools/check_symplectic_mixed_grade_codes.py).
- [Complete report](../evidence/symplectic-mixed-grade-codes.json).
- [Three native completion logs](../evidence/symplectic-mixed-grade-screens).

The result is specific to the direct sum of the grade 2, 4, and 6 radicals on
four complementary pairs and the two repairs defined above. It does
not rule out other mixed-grade maps, different output projections, or
larger coefficient constructions.
