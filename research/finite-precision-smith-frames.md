# Exact regular Smith frames with bounded polynomial representatives

## Purpose and result

Complete-entry synthesis needs only finitely many coefficients of each
local regular frame. Expanding a fraction-free reduction as an unrestricted
polynomial can produce large irrelevant tails, especially after the source
has become non-Gaussian. This note gives an exact bounded computation for
the same local row conditions.

For a square binary polynomial matrix A and a positive precision P, the
algorithm returns polynomial matrices L,R with invertible constant matrices,
and a diagonal matrix D, such that

\[
LAR\equiv D\pmod{z^P}.
\tag{1}
\]

Every returned polynomial has degree below P. If all diagonal orders are
less than P, they are the complete Smith orders of A. Every row-divisibility
target at most P can be imposed using these representatives exactly. This
is not a claim that LAR is diagonal as an unrestricted polynomial matrix.

## 1. Operations and unit determinants

At each step choose a remaining matrix entry of smallest order r, and move
it to the pivot position. Write the pivot as z^r u and a same-column entry
as z^r v. The latter divisibility holds by the pivot's minimality. Apply
the row operation

\[
\text{row}_i\longmapsto u\,\text{row}_i-v\,\text{row}_k.
\tag{2}
\]

Its determinant multiplier is u, whose constant term is nonzero. The
operation is therefore invertible over the power-series ring. It clears
the selected column entry. Apply the corresponding operation on columns
to clear the pivot row. Record the row operations in L and the column
operations in R.

After every operation, reduce all polynomial entries modulo z^P. These
reductions preserve the product congruence (1). They also preserve the
constant matrices of L and R, so their polynomial representatives remain
actual regular units. There is no loss of precision proportional to r:
the factors u and v are extracted from known polynomial representatives
and then used as polynomial multipliers. No matrix is divided by z^r.

If the remaining matrix vanishes before all pivots are found, the routine
rejects the precision as unresolved. It does not infer a smaller generic
rank. The driver uses a separate complete native profile to choose
P=max(order)+allowance, then verifies that the bounded reduction recovers
those orders and the full congruence.

## 2. Why the finite diagonal determines the complete profile

Suppose (1) holds and every diagonal entry has order r_i<P. In the formal
power-series ring, the exact product has the form

\[
LAR=\operatorname{diag}(z^{r_i})\,U,
\]

where U has a diagonal unit constant matrix: the off-diagonal remainder
in row i has order at least P-r_i>0 after division by z^r_i. Thus U is
invertible. Since L and R are also invertible, the sorted r_i are exactly
the Smith orders of A. A later hidden channel is impossible in this case.

For a different tensor T in the same global coefficient family, use the
actual polynomial representatives L,R. Their coefficients below P determine
all coefficients of LTR below P, since the cleared source entries have
nonnegative powers. Consequently row conditions
ord((LTR)_ij)>=tau_i, with tau_i<=P, are complete affine coefficient
equations. No assumption about the uncomputed higher frame coefficients is
needed. Their high-degree tails can be zero in the chosen representative.

This construction complements the
[classification of regular row-frame cells](formal-row-flag-cells.md),
which describes all independent local frames for fixed targets. The bounded
Smith routine instead constructs one useful frame efficiently from a given
reference tensor. Neither claims that searching one frame covers every
crossing with the desired Smith orders.

## 3. Exact controls and current references

The independent checker constructs matrices with known complete profiles
using polynomial triangular units. Across sizes two through seven it checks
108 finite-frame cases, compares 108 full untruncated products with the
claimed congruence, and verifies invertibility of both constant frames.
Twelve small cases also agree with unrestricted fraction-free reduction.
Two unresolved precision or rank cases are rejected.

Two non-Gaussian four-symbol reference tensors supply a larger control.
Each comes from a complete-entry coefficient search followed by independent
optimization of all physical-symbol heights. Its rounded integer gauges
give exactly the same complete profile for its binary low coefficient plane
and its original F4 tensor. The replayed rates are approximately
0.7450485724 and 0.7449315952. These remain above the required threshold.

For each reference, the checker constructs the bounded frames, compares
the complete profiles with an independent restriction-of-scalars Smith
calculation, and verifies two literal synthesized F4 tensors against the
full polynomial row conditions. There are 316 independent block-Toeplitz
rank checks in the combined record.

```sh
make verify-regular-frame-cells verify-finite-smith-frames
```

Implementation and evidence:

- [`finite_precision_smith_frames.py`](../tools/finite_precision_smith_frames.py)
- [`check_finite_precision_smith_frames.py`](../tools/check_finite_precision_smith_frames.py)
- [`finite-precision-smith-frames.json`](../evidence/finite-precision-smith-frames.json)
- [`affine-degenerate-witness-0-integer-chamber.json`](../evidence/affine-degenerate-witness-0-integer-chamber.json)
- [`affine-degenerate-witness-1-integer-chamber.json`](../evidence/affine-degenerate-witness-1-integer-chamber.json)

The new reference-walk driver alternates this choice of frame and physical
heights with further complete-entry coefficient searches. Each walk keeps
its own frame and coefficient window fixed; changes between walks are
explicit. The finite searches are separate from the exact frame theorem.
