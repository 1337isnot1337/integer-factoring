# The twelve-site separator direction completes the coefficient algebra

## Result

The twelve-site coefficient algebra is the ten-dimensional algebra with
basis `C_0,C_2,...,C_10` from
[the exact table](bounded-spin-twelve-coefficient-table.md).  In that
basis the images of the replica swap `R` and the compressed transposition
class sum `T` have coordinates

```text
R = [1,0,0,0,0,0,0,0,0,0],
T = [3,3,1,2,1,3,3,1,2,1].
```

The table's unit is `C_0`, so at twelve sites `R` is the unit in the
coefficient algebra; the replica swap is still nontrivial on the
`7921`-dimensional module.  Consequently

```text
dim <R,T>                 = 5,
dim <R,T,b_4>             = 7,
dim <R,T,b_2>             = 7,
dim <R,T,b_2+b_4>         = 10.
```

Here `b_i` denotes coordinate index `i` in the published basis
`[0,2,3,4,5,6,7,8,9,10]`; in particular `b_2` has polynomial coefficient
`3` and `b_4` has coefficient `5`.  No **single basis direction**
completes the algebra.  The scan of all pairwise sums finds two completing
single directions:

```text
b_2 + b_4   (coefficients 3 and 5),
b_1 + b_5   (coefficients 2 and 6).
```

These are pair sums in the basis coordinates; scalar multiples complete
as well, and higher multi-term combinations were not enumerated.

The independent pairs `[b_1,b_5]`, `[b_2,b_4]`, and `[b_2,b_5]` also
complete it.

## Ten-site control and the role reversal

The same entry-functional method, applied to eight independent ten-site
path quadruples, recovers the published coordinates

```text
R = [1,0,0,0,0,0,0,0],
T = [0,2,0,1,0,2,2,0].
```

With the published ten-site regular representation and unit, the checker
also recovers the ten-site dimensions

```text
dim <R,T>           = 6,
dim <R,T,b_4>       = 8,
dim <R,T,b_2+b_4>   = 7,
dim <R,T,b_2,b_4>   = 8.
```

Thus the ten-site result used the **basis direction** `b_4` and the
separator was the sum `b_2+b_4`.  At twelve sites the roles reverse: the
basis direction `b_4` does **not** complete, while the separator sum
`b_2+b_4` does complete.  The ten-site structural mechanism therefore
does transfer, but the direction to carry over is the separator
combination, not the basis coordinate alone.

The exact dimension `5` is a statement about the coefficient algebra's
regular representation.  It is not the dimension of the orbit of a
generic module vector under `R,T`; the latter was bounded below by `8` in
[the orbit note](bounded-spin-twelve-coefficient-dimension.md).  The two
representations need not have the same dimension because `R` collapses
to the algebra unit here.

## Character values on the primitive central idempotents

Using the trace-left form from the exact table and the eight primitive
central idempotents from
[the idempotent certificate](../evidence/bounded-spin-twelve-coefficient-idempotents.json),
the scalar character values of `T` and the separator `S=b_2+b_4` are

```text
index  rank  kind   T  S
   0     1   scalar 0  2
   1     1   scalar 2  0
   2     2   dual   3  2
   3     1   scalar 3  3
   4     2   dual   0  0
   5     1   scalar 0  3
   6     1   scalar 1  4
   7     1   scalar 3  1
```

The pair `(T,S)` takes eight distinct values.  Thus the two generators
separate all eight primitive central blocks, even though `T` alone has
only four distinct scalar values and `S` alone has five.  The replica
swap is the algebra unit and has value `1` on every block, as expected.
This is a statement about the coefficient algebra's characters; it does
not compute the action of the idempotents on the twelve-site cut modules.

## Consequence for the twelve-site primary target

The missing algebra direction at twelve sites is now a **single explicit
coefficient combination**, `b_2+b_4` in the published basis.  In
particular, adding one vector is enough to generate the full
ten-dimensional coefficient algebra; no two-dimensional complement is
needed.  This removes the algebra-generation part of the twelve-site
primary calculation.

It does **not** by itself determine the idempotent cut-module ranks, the
two-cut pairing, or the factoring exponent.  The primitive central
idempotents are now computed in
[the idempotent note](bounded-spin-twelve-coefficient-idempotents.md):
the algebra is `F5^6` times two dual-number blocks, with eight primitive
central idempotents of regular ranks `1,1,1,1,1,1,2,2`.  What remains is
to act with the idempotents on the twelve-site cut modules and to compare
the separator direction with the coarse-block characters.  The complete
bound is still `1/2`.

## Reproduction

```sh
python3 tools/check_bounded_spin_twelve_missing_direction.py
```

The checker reads the exact twelve-site table and the ten-site model,
computes the `R,T` coordinates from the stored rank-ten entry witness,
recovers the ten-site control, performs the basis and linear-combination
scans over `F_5`, and writes
[the receipt](../evidence/bounded-spin-twelve-missing-direction.json).
It runs in about one second.

The coefficient table, dimension certificate, orbit bound, and ten-site
model are credited to the bounded-spin lane.  This note continues the
[Discussion #51 exchange](https://github.com/1337isnot1337/integer-factoring/discussions/51).
