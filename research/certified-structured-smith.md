# Structured Smith elimination with explicit precision certificates

## Result

The [tracked-error engine](../tools/native_certified_rank_one_smith_f4.cpp)
removes the unresolved precision assumption from the earlier
[rank-one search accelerator](rank-one-smith-search.md). It returns a
complete profile only when every pivot is certified by explicit
absolute series-error bounds. Insufficient precision returns an error.

This is an ordinary algorithmic proof over `F4((z))`, accompanied by
independent finite controls. It is not a new Lean formalization and
does not establish a better factoring exponent.

## 1. Representation invariant

Each stored scalar represents

```math
a=z^v p(z)+O(z^E),\qquad \deg p<R.
```

Here `R` is the fixed working precision. If the stored polynomial is
nonzero, it is normalized to have nonzero constant coefficient, and
`v<E`. Its valuation `v` is therefore known exactly. If all stored
coefficients vanish, the object remains `O(z^E)` with valuation lower
bound `E`. It becomes exact zero only when the error bound is infinite.

The input vectors are exact finite polynomials, so their initial error
bounds are infinite. Every discarded coefficient creates a finite
bound. Renormalizing after cancellation preserves the absolute bound;
it never fills an unknown tail with asserted zero coefficients.

## 2. Propagating error bounds

Suppose the nonzero leading terms of `a,b` are known, with valuations
`v_a,v_b` and absolute error exponents `E_a,E_b`. Before the additional
truncation imposed by storage, valid output error bounds are:

| Operation | Absolute error exponent |
| --- | --- |
| `a+b` | `min(E_a,E_b)` |
| `ab` | `min(E_a+v_b,E_b+v_a)` |
| `a/b` | `min(E_a-v_b,E_b+v_a-2v_b)` |

The product's omitted cross-error has exponent `E_a+E_b`, which is
no smaller than the two displayed terms. The quotient formula follows
by expanding the inverse of the denominator's known unit part.
Division is forbidden unless the denominator's leading coefficient
has been certified nonzero.

If a computation is truncated at a window of `R` coefficients starting
at exponent `w`, the error bound is also capped at `w+R`. An exact
polynomial addition or multiplication whose complete degree fits in
the stored window need not incur this additional loss. Division by
an exact monomial is likewise exact on the stored polynomial part.
For a nonconstant unit denominator, Newton inversion is performed
modulo `z^R`, and the corresponding truncation bound is retained.

For an unknown zero `a=O(z^E_a)`, multiplication by a known-leading
`b` gives `O(z^(E_a+v_b))`, and division gives
`O(z^(E_a-v_b))`. Two unknown-zero factors give
`O(z^(E_a+E_b))`. These rules are implemented without inventing a
leading coefficient. Bounds that leave the supported integer exponent
range cause failure rather than being relabeled infinite.

The formulas preserve the representation invariant by induction over
every scalar operation. They allow intermediate Laurent valuations;
regularity is required for the actual matrix operations, not for
every compact generator appearing in a Schur formula.

## 3. Certifying each matrix pivot

The engine uses the exact diagonal-plus-rank-one and arrowhead identities
proved in [the earlier note](rank-one-smith-search.md). It evaluates a
certified lower bound for every active matrix entry. A candidate pivot
is accepted only if:

1. Its leading coefficient is known and nonzero.
2. Its valuation is no greater than every other entry's certified
   valuation lower bound.
3. Its valuation lies in the requested nonnegative output window.

The first two conditions prove that it is a true minimum-valuation
entry of the exact active matrix. Every row and column elimination
multiplier is then regular: its numerator valuation is at least the
pivot valuation. The elementary operations are invertible over
`F4[[z]]`, and they split off the pivot from its exact Schur complement.
Induction therefore gives the Smith orders, provided the routine
successfully returns all of them.

For `D+uv^T`, every diagonal entry is evaluated explicitly. Among the
off-diagonal products, a row needs only the two smallest column-factor
valuation lower bounds, since its matching diagonal column is excluded.
The same minima also expose every possible smaller unknown entry.
If such uncertainty prevents certifying the chosen pivot, the routine
fails safely. It may conservatively return unresolved even when a
higher-precision run would succeed.

For an arrowhead, the center, leaf diagonals, and two border vectors
list all possibly nonzero entries. The remaining leaf off-diagonals
are structurally zero. No unexamined dense matrix tail is assumed zero.

## 4. Transfer through finite polynomial frames

The complete one-copy frame identity is checked exactly modulo `z^P`:

```math
L M R=D\pmod {z^P},
```

with `L,R` regular and invertible. Tensor-product frames and transformed
rank-one vectors give the corresponding structured matrix to that same
absolute precision. Matrices equal modulo `z^P` have the same truncated
module maps modulo `z^k` for every `k<=P`.

If the tracked-error engine certifies all `n` orders below `P`, these
truncated maps determine every order. Indeed, the rank of multiplication
on truncated series is `sum_i max(k-a_i,0)`; successive differences
recover the number of orders below each `k`. Thus the complete original
matrix has exactly the returned profile. The independently checked
regular frames and the tracked-error pivot proof together justify this
transfer. The condition that all returned orders lie below `P` is
essential.

## 5. Independent controls

The [small controls](../evidence/certified-rank-one-smith-f4-controls.json)
repeat 192 complete-matrix comparisons at two work precisions, plus
144 independent binary block-Toeplitz rank checks. Every Schur pivot
branch occurs. Singular and out-of-window cases return unresolved.

The [large frame controls](../evidence/certified-two-copy-frame-controls.json)
revisit ten separately retained 625-by-625 physical examples, using
work precisions 256 and 512. This gives 40 complete profile comparisons.
Six examples use products of one-copy supports, and four use physical
vectors with correlations between copies. Every comparison agrees
with its previously computed complete original matrix. The record
retains the minimum certified pivot error margin and the source hashes.

```sh
python3 tools/check_certified_rank_one_smith_f4.py
OPENBLAS_NUM_THREADS=1 python3 tools/check_certified_two_copy_frames.py
```

Earlier search records retain the engine and hashes actually used at
their creation. This new proof does not silently relabel those records;
a separate replay through the tracked-error engine can certify their
recorded structured profiles.
