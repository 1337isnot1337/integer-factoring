# The twelve-site coefficient algebra: exact ten-dimensional table

## Result

The twelve-site bounded-spin coefficient span now has an **exact
multiplication table**. In the coefficient basis

```text
0, 2, 3, 4, 5, 6, 7, 8, 9, 10
```

the table is commutative, associative, and unital. The unit is the first
basis vector `C_0`; the ten-site unit needed a combination. The
[table receipt](../evidence/bounded-spin-twelve-coefficient-table.json)
also records the regular representation, the rank-ten entry witness, and
two independent holdouts. The [build tool](../tools/build_bounded_spin_twelve_coefficient_table.py)
recomputes the table from the path constructor in about four minutes; the
[checker](../tools/check_bounded_spin_twelve_coefficient_table.py) replays
all algebra identities, the unit, the regular representation, and the
stored holdouts.

The ten-site control was reproduced before this run: the same
entry-functional method recovers the published eight-dimensional regular
representation exactly. This is evidence that the twelve-site table uses
the same cut normalization and multiplication orientation.

## Method

Use the ten rank-ten path quadruples from the
[dimension certificate](bounded-spin-twelve-coefficient-dimension.md) as
ten entry functionals. For a selected entry `((x,y),(u,v))`, the row
operator is indexed by `(u',v')` and the column operator by `(x',y')`.
The product entry is

```text
sum_s R_i[((x,y),s)] R_j[(s,(u,v))].
```

For each of the ten selected entries, the build stores the ten full row
vectors and ten full column vectors and solves the ten-by-ten linear system
over `F_5`. Two random entries outside the witness are then checked
directly against the table. Both holdouts agree in all one hundred basis
products.

## What this does and does not give

This removes the multiplication table from the twelve-site primary
calculation. The central idempotents, the primary dimensions, the cut
module ranks, and the rank formula are still open. In particular, the
table alone does not determine the two-cut rank product or a factoring
exponent. The complete factoring score remains `1/2`.

The primitive central idempotents are now computed in
[the primary-idempotent note](bounded-spin-twelve-primary-idempotents.md):
the algebra has eight primaries with regular block dimensions
`1,2,1,1,1,1,1,2`.  Their ranks on the twelve-site cut modules remain
open.
