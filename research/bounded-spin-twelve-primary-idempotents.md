# Primitive central idempotents of the twelve-site coefficient algebra

## Result

The ten-dimensional twelve-site coefficient algebra over `F_5` has
exactly **eight primitive central idempotents**.  In the coefficient
basis `0,2,3,4,5,6,7,8,9,10` their vectors are

```text
e0 = 4 C8
e1 = C0+4C2+C3+3C4+2C6+2C7+3C8+4C9
e2 = C0+2C3+2C4+C5+2C6+4C7+4C8+3C9+4C10
e3 = 2C0+3C2+C3+4C4+2C5+2C6
e4 = 4C0+C2+2C3+4C5+4C6+C9
e5 = 4C0+C2+2C3+3C4+4C5+C7+2C8+3C9+4C10
e6 = 2C0+3C2+C3+4C4+2C5+4C8+4C9+2C10
e7 = 2C0+3C2+C3+4C4+2C5+3C7+2C8+C9.
```

They are pairwise orthogonal, each squares to itself, and their sum is
the table unit `C0`.  Their regular block dimensions are

```text
1, 2, 1, 1, 1, 1, 1, 2,
```

so the ten-dimensional algebra splits as a direct product of eight
local algebras.  The ten-site eight-dimensional algebra has seven
primaries; the two extra twelve-site dimensions do not simply add one
new block in the ten-site pattern.

This is an exact finite calculation over `F_5`.  It does not compute the
idempotent ranks on the `89^2=7921` path-pair cut modules, does not
supply an integer lift, and does not change the complete factoring
score.  The next target is the cut-module ranks of these eight
idempotents (or an exact trace form that determines them).

## Method

The [ten-dimensional table](bounded-spin-twelve-coefficient-table.md)
supplies the commutative associative structure constants.  For the
generic element

```text
x = C0+2C3+C4+2C5+2C6+4C7+C8+4C9
```

the minimal polynomial is

```text
m(x) = x(x-1)^2(x-2)(x-3)(x-4)^2
     = x + 4x^3 + 4x^5 + x^7       over F_5.
```

Its five distinct roots give five spectral projectors.  Their regular
dimensions are `1,4,1,1,3`.  The two repeated-root spectral blocks split
further.  Enumerating every idempotent in each generalized eigenspace
(after choosing an explicit image basis) gives the eight atoms above.
For each atom the same enumeration finds only the zero idempotent and
the atom itself, so each is primitive.  The checker also verifies all
products, the unit sum, the regular block dimensions, and the
ten-dimensional sum.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_primary_idempotents.py \
  --output evidence/bounded-spin-twelve-primary-idempotents.json
make verify-bounded-spin-twelve-primary-idempotents
```

The checker reads the stored exact table and does not rebuild the
length-twelve tensor.  The table's source hash is recorded in the
evidence receipt.
