# Exact F4 row frames inside a compressed coefficient fiber

## Result and scope

The complete-entry coefficient solver now permits genuine F4 polynomial
row transformations at either local place. A common binary affine fiber
can be computed once and reused across many such frames. All subsequent
equations, including their nonbinary coefficients and affine constants,
are solved over F4. This enlarges the explicitly searched frame family.
It does not assert that a finite collection of frames is complete, or that
an improving crossing exists.

The construction applies to the
[typed rational coefficient lattice](typed-charge-lattice-synthesis.md),
including its strict polynomial specialization. It retains every entry
coefficient, the endpoint conditions, both complete local profiles, and all
physical-symbol costs. The current four-symbol strict polynomial instance
has 6,144 F4 coefficient variables; the stated common floors reduce the
working system to 744 F4 variables.

## 1. A common floor preserved by elementary row operations

At one local place, let the reference row orders be s_i. Fix integers
L>=0 and R>0, and put

\[
 f_i=\max(0,s_i-L),\qquad c_i=s_i+R.
 \tag{1}
\]

The coefficient equations are constructed completely through exponent
c_i-1 in row i, retaining the column index of every equation. Consider
an elementary transformation

\[
 \operatorname{row}_i\longleftarrow
 \operatorname{row}_i+\gamma z^p\operatorname{row}_j,
 \qquad i\ne j,\quad\gamma\in\mathbb F_4^\times,
 \quad p\ge\max(0,s_i-s_j).
 \tag{2}
\]

It is a polynomial unit matrix, with determinant one. Its inverse is the
same elementary operation because the characteristic is two. It preserves
the lattice in which row i is divisible by z^s_i. It also preserves the
common floor lattice: the map s to max(0,s-L) is increasing and
one-Lipschitz, so

\[
 p\ge f_i-f_j.
 \tag{3}
\]

Both the operation and its inverse preserve that lattice. Consequently,
the set of coefficient assignments satisfying all common floor equations
is exactly unchanged after any finite sequence of the permitted operations.
This is an equality of affine fibers, not merely one containment.

Coefficient truncation is exact for every requested target t_i<=c_i.
Indeed,

\[
 c_i-p\le c_j.
 \tag{4}
\]

Thus each source coefficient needed to update a target coefficient below
c_i was already constructed. Repeating operations preserves this property.
Coefficients outside the computed window are never assumed to vanish.

## 2. Reduce once over the binary field, then retain all F4 coefficients

The original affine equations have binary coefficients and binary right
sides. Their binary elimination remains valid over F4. Choose a basis of
free variables for the common floor system. Express each original
coefficient equation in those variables, preserving its exponent, column,
and affine right side.

Apply (2) directly to these reduced coefficient tables. Multiplication by
omega, with omega squared equal to omega plus one, acts on the two packed
binary coefficient planes by

\[
 (a,b)\longmapsto(b,a+b).
 \tag{5}
\]

The affine right side is multiplied by the same scalar. After these
transformations, the target equations form a genuine affine F4 system.
Packed Gaussian elimination gives exact feasibility and its full solution
dimension. Sampling free F4 variables, solving for pivots, and lifting
through the original binary affine equations produces a literal complete
tensor. The high coefficient plane uses the homogeneous binary lift; only
the low plane receives the original affine constant.

The implementation is
[`colored_row_fibers.py`](../tools/colored_row_fibers.py). The finite search
driver is [`search_colored_row_fibers.py`](../tools/search_colored_row_fibers.py).
Its records include every frame operation, exact row targets, source
hashes, complete tensor witnesses, and complete profiles. Numerical root
optimization is explicitly distinguished from the exact coefficient
construction.

## 3. Independent verification

The checker first enumerates all assignments in 160 small affine systems
over F4, using a separate literal field multiplication routine. It compares
existence and exact solution counts with elimination, then checks 1,136
sampled assignments directly.

Eight larger controls construct four-symbol tensors with nonzero high
coefficient planes. For each tensor, the checker verifies entrywise
polynomial divisibility independently of the coefficient tables. It forms
both full weighted matrices, multiplies by the original left and right
polynomial frames without truncation, and applies every recorded F4 row
operation literally. Each target coefficient is then checked for zero.
The original common floor conditions are checked independently as well.
Finally, complete native Smith profiles satisfy every asserted prefix bound.

```sh
make verify-colored-row-fibers
```

The finite control record is
[`colored-row-fiber-controls.json`](../evidence/colored-row-fiber-controls.json).
Search records can be replayed with the checker's `--search-records` option.
An interrupted search remains explicitly marked incomplete; its valid
saved tensors can still be replayed, but it is not a completed catalogue.

The operations in (2) preserve a selected reference lattice. They do not
range over every constant flag, every formal frame cell, or every possible
choice of coefficient degree. The
[general row-frame completeness theorem](formal-row-flag-cells.md)
describes the broader space; this implementation supplies a practical,
exact subset of that space.
