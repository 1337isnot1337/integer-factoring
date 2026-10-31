# Structured Smith searches for pure tensor corrections

## Scope

The [native structured engine](../tools/native_rank_one_smith_f4.cpp)
computes exploratory integer Smith orders for a diagonal polynomial matrix
plus a rank-one correction over `F4`. It alternates between that form and
an arrowhead matrix, using exact characteristic-two Schur identities.
Every pivot is chosen at the smallest valuation among all represented
entries. The implementation uses finite relative series precision.

The algebraic identities below hold over the rational function field.
The current implementation does **not** provide a uniform proof that
every discarded intermediate series tail is harmless. Its outputs are
search predictions until the original complete matrix is replayed.
Nothing in this note changes the factoring bound.

## 1. Two structured states

Write the first state as

```math
A=D+uv^{\mathsf T},\qquad D=\operatorname{diag}(d_i).
```

In characteristic two, eliminating the diagonal pivot
`p=d_i+u_i v_i` leaves

```math
D_{\widehat i}+\frac{d_i}{p}
 u_{\widehat i}v_{\widehat i}^{\mathsf T}.
```

If the pivot instead is `p=u_i v_j`, with `i != j`, the remaining
matrix has an arrowhead. Its distinguished row is the old row `j`,
and its distinguished column is the old column `i`. The center is

```math
a=\frac{d_i d_j}{u_i v_j}
  +\frac{d_j v_i}{v_j}+\frac{d_i u_j}{u_i}.
```

For every remaining leaf `k`, its diagonal remains `d_k`, its
distinguished-column entry is `d_i u_k/u_i`, and its
distinguished-row entry is `d_j v_k/v_j`. All other off-diagonal
leaf entries vanish. These formulas follow by substituting the
pivot into the Schur complement; they use no nonsingularity
assumption about `D`.

For the arrowhead state

```math
\begin{pmatrix}a&v^{\mathsf T}\\u&D\end{pmatrix},
```

there are four pivot types:

| Pivot | Remaining state |
| --- | --- |
| Center `a` | `D+uv^T/a` |
| Leaf diagonal `d_i` | Arrowhead with center `a+u_i v_i/d_i` |
| Leaf column `u_i` | Arrowhead with center `v_i+a d_i/u_i` and remaining column multiplied by `d_i/u_i` |
| Leaf row `v_i` | The transpose of the preceding formula |

At a global minimum-valuation pivot, every actual elimination
multiplier is regular. Thus exact versions of these operations
preserve the remaining Smith orders and expose the pivot order.
Some compact generator formulas can nevertheless involve intermediate
Laurent series. This is why finite relative-precision implementation
controls are kept separate from the algebraic proof.

The engine finds the smallest off-diagonal `u_i v_j` without forming
the full matrix: for each `i`, only the two smallest nonzero
valuations among the `v_j` are needed. The storage is linear in the
matrix size, and the number of series operations is quadratic.

## 2. Polynomial frames for two copies

For each cut of a fixed one-copy tensor `B`, the existing regular
frame routine constructs polynomial matrices `L,R` with unit
determinants in `F4[[z]]` and checks the full identity

```math
L B R=D\pmod {z^P}.
```

Consequently the two-copy matrix has frame
`(L tensor L),(R tensor R)` and diagonal `D tensor D` to the same
absolute precision. A physical correction pure at both cuts becomes
a rank-one correction in these frames. All physical monomial weights
are included before the transform.

The [pure two-copy adapter](../tools/pure_two_copy_frames.py) checks
the literal Cartesian support at each cut, including the row/column
index ordering. It factors the weighted support, applies the regular
frames, and forms the two Kronecker vectors. It retains the full
endpoint-zero polynomial, reversed at infinity inside the complete
degree bound. Subtracting the known clearing scalar recovers the
original weighted integer orders.

For the retained five-letter degree-18 reference, these are 625 by
625 matrices at each cut. The one-copy frames have size 25 by 25.
The construction is an exact finite-jet identity. Only the subsequent
structured elimination has the implementation limitation stated above.

## 3. Independent controls

The [engine controls](../evidence/rank-one-smith-f4-controls.json)
compare 192 random diagonal-plus-rank-one matrices with the complete
matrix oracle, using sizes 1, 2, 3, 4, 5, 8, 12, and 16. Every case
is compared at two structured work precisions. There are also 144
independent binary block-Toeplitz rank checks. All six pivot branches
occur, including both nonsymmetric arrowhead branches. Singular and
out-of-window negative controls return unresolved.

The [frame controls](../evidence/pure-two-copy-frame-controls.json)
replay 24 small physical examples and six full five-letter two-copy
examples against their original complete matrices, for 60 complete
profile comparisons. The small cases include negative physical
heights before scalar clearing. The six larger cases use different
pure supports, polynomial shifts, endpoint powers, and all nonzero
field coefficients.

In that retained run the structured two-profile queries took roughly
0.04--0.19 seconds on the six larger controls, while complete-matrix
queries took roughly 3.8--4.4 seconds. These timings describe this
finite run, not a general speed guarantee.

```sh
python3 tools/check_rank_one_smith_f4.py
OPENBLAS_NUM_THREADS=1 python3 tools/check_pure_two_copy_frames.py --large
```

The checks establish agreement on the stated finite cases. New
search records must distinguish structured queries from complete
matrix replays, retain unresolved cases, and replay every candidate
used to claim an improvement. A favorable floating pressure would
also need a rigorous pressure certificate before affecting the
factoring theorem.
