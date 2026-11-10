# Primitive central idempotents of the twelve-site coefficient algebra

## Result

The twelve-site coefficient algebra from the
[exact ten-dimensional table](bounded-spin-twelve-coefficient-table.md) has
exactly eight primitive central idempotents.  Six blocks are scalars and
two blocks are dual numbers:

```math
A_{12}\cong \mathbb F_5^{\,6}\times
  \mathbb F_5[\varepsilon_1]/(\varepsilon_1^2)\times
  \mathbb F_5[\varepsilon_2]/(\varepsilon_2^2).
```

The regular-representation ranks are therefore

```text
1, 1, 1, 1, 1, 1, 2, 2,
```

summing to the algebra dimension `10`.  The trace form has rank `8` and a
two-dimensional radical.  The radical is spanned by the two nilpotents
`\varepsilon_1,\varepsilon_2`, and it is square-zero.  The full idempotent
vectors, block bases, trace forms, and radical basis are recorded in the
[idempotent certificate](../evidence/bounded-spin-twelve-coefficient-idempotents.json).

This computes the **coefficient algebra's primary decomposition**.  It does
not compute the ranks of those idempotents on the `89^2=7921`-state cut
modules, the two-cut pairing, or any factoring exponent.  The complete
score remains `1/2`.

## Method

The [table builder](../tools/build_bounded_spin_twelve_coefficient_table.py)
supplies multiplication in the basis

```text
0, 2, 3, 4, 5, 6, 7, 8, 9, 10
```

with unit `C_0`.  The decomposition is purely algebraic and does not build
a cut operator:

1. Choose a generic element `a` in the current block.
2. Compute its minimal polynomial over `F_5`.
3. Factor the polynomial into powers of linear factors.  Every minimal
   polynomial in this algebra is split; the builder asserts this rather
   than assuming it.
4. Use the Chinese remainder theorem to form the orthogonal idempotents
   for the distinct linear factors.  A repeated factor remains in one
   local block; distinct factors split the block.
5. Repeat on every block with more than one scalar dimension.

The random choices use a fixed seed, and the primitive vectors are sorted
canonically.  Each resulting block is checked independently:

- scalar block: the trace form is `[[1]]`;
- dual block: the trace form has rank one, and an explicit nonzero
  nilpotent `\varepsilon` satisfies `\varepsilon^2=0` and
  `e\varepsilon=\varepsilon`.

Thus every block is local and every listed idempotent is primitive.  The
eight idempotents are pairwise orthogonal and sum to the unit.  The full
trace-form radical is computed separately from the kernel of the `10 x 10`
trace form; its dimension is `2` and its span equals the span of the two
block nilpotents.

The [checker](../tools/check_bounded_spin_twelve_coefficient_idempotents.py)
recomputes the decomposition from the table and compares the complete
mathematical certificate, excluding only the timing field.

## What changed

Before this computation the table gave the algebra's multiplication but not
its primary structure.  The ten-site algebra was `F_5^6 x dual`, with seven
primitive central idempotents; the twelve-site algebra was known to have
dimension `10` but its block structure was open.  The result here is the
exact twelve-site block structure, including the correct regular ranks.

The two extra dimensions beyond the eight-dimensional `R,T` orbit model
are not two further non-scalar blocks.  The full algebra has six scalar
blocks and two dual blocks; an eight-dimensional orbit model can only see
part of this structure.  The orbit algebra is identified abstractly as
`F5^4 x dual x dual` in
[the orbit-algebra note](bounded-spin-twelve-RT-orbit-algebra.md); its actual
embedding in the ten-dimensional algebra remains a separate check.

## What does not follow

The [diagonal-trace note](bounded-spin-diagonal-trace.md) gives normalized
trace forms modulo five.  Their values on these idempotents are a checksum,
not an exact rank: a rational lift or a direct exact rank computation is
still required.  Moreover, the
[image-floor note](bounded-spin-image-floor.md) records that the central
idempotent decomposition of the ten-site algebra has square-root rank sum
far above the target; the central idempotents are not automatically the
noncentral decomposition needed by the crossing.

The next object is the exact row and column ranks of these eight
idempotents on the twelve-site cut modules, followed by their two-cut
pairing and a comparison with the `89^sqrt(2)=571.286...` target.

## Reproduction

```sh
python3 tools/build_bounded_spin_twelve_coefficient_idempotents.py
python3 tools/check_bounded_spin_twelve_coefficient_idempotents.py
make verify-bounded-spin-twelve-coefficient-idempotents
```

The builder takes well under a second; the table build is separate and
takes about four minutes.  The certificate pins the SHA-256 of the table
receipt.  This note continues the bounded-spin lane and the
[Discussion #51 exchange](https://github.com/1337isnot1337/integer-factoring/discussions/51).
