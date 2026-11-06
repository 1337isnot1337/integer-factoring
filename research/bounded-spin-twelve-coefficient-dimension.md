# The twelve-site coefficient algebra has dimension ten

## Result

The twelve-site bounded-spin coefficient span is exactly ten-dimensional.
The upper bound is the uniform `n-2` coefficient identity; the lower
bound is an explicit rank-ten minor.

This corrects a natural but false guess that the twelve-site algebra
would remain eight-dimensional like the published ten-site model. The
ten-site `8 x 8` regular representation and its seven idempotents cannot
be transplanted to twelve sites. The next object is a ten-dimensional
twelve-site multiplication table, not a larger instance of the ten-site
one.

No primary decomposition, cut-rank bound, or factoring score changes.

## Upper bound

[The three uniform jet identities](bounded-spin-three-jets.md) give, for
even `n >= 6`,

```math
2C_1+nC_0=0,\qquad
2\sum_{j=0}^n jC_j-n\sum_{j=0}^n C_j=0,\qquad
2C_{n-1}+nC_n=0.
```

The three constraints are independent over every odd characteristic.
Hence the coefficient span has dimension at most `n-2`. At `n=12` this
is `10`.

## Lower bound

Evaluate the normalized original-pencil coefficient vector on ten
explicit path quadruples. The checker records the quadruples and verifies
that the resulting `10 x 13` matrix has rank ten; its `10 x 10` minor is

```text
coefficient indices: 0, 2, 3, 4, 5, 6, 7, 8, 9, 10.
```

The same evaluation at ten sites reproduces the published eight-dimensional
basis

```text
0, 2, 3, 4, 5, 6, 7, 10,
```

which is an independent layout control: the constructor and coefficient
ordering recover the known ten-site model before the twelve-site minor is
read.

Thus the twelve-site dimension is at least ten and at most ten, hence
exactly ten.

## What is now reduced

The twelve-site primary calculation needs:

1. the ten-dimensional coefficient multiplication table or regular
   representation at length twelve;
2. the central idempotents in that ten-dimensional algebra;
3. the exact rank of each idempotent, or an exact trace form that
   determines it.

The dimension is no longer an unknown. The ten-site idempotents and
their seven primary dimensions are not a valid short cut, because the
algebra grows by two dimensions.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_coefficient_dimension.py
make verify-bounded-spin-twelve-coefficient-dimension
```

The checker samples 50,000 quadruples to find the witness, then verifies
the stored ten-row minor exactly. It writes
`evidence/bounded-spin-twelve-coefficient-dimension.json`.

The bounded-spin construction, path frames, normalization, and the
three-jet upper bound are credited to the bounded-spin lane. This note
continues the [Discussion #51 exchange](https://github.com/1337isnot1337/integer-factoring/discussions/51).
