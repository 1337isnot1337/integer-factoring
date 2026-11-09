# One coefficient direction completes the ten-site algebra

## Result

The replica swap `R` and the compressed transposition class sum `T`
generate a six-dimensional subalgebra of the eight-dimensional
ten-site coefficient algebra. The
[transposition subalgebra note](bounded-spin-transposition-subalgebra.md)
showed that `R,T` are insufficient and that the coarse `335` block
must split. This note identifies a single missing direction.

In the basis of
[the ten-site algebra model](../evidence/bounded-spin-ten-algebra-model.json),
the coordinates are

```text
R = [1,0,0,0,0,0,0,0],
T = [0,2,0,1,0,2,2,0].
```

Adding the fourth basis coordinate (polynomial coefficient position `5`
in the model's `basis_indices` `[0,2,3,4,5,6,7,10]`) raises the generated
algebra from dimension `6` to the full dimension `8`:

```text
dim <R,T>            = 6
dim <R,T,b_4>        = 8
dim <R,T,b_2>        = 7
dim <R,T,b_2+b_4>    = 7.
```

Thus one further coefficient direction is enough; a two-dimensional
complement is not needed to generate the full ten-site algebra.

## Splitting the coarse `335` block

The two quotient directions individually do not separate all three
primaries `34,1,300` inside the `R,T`-coarse `335` block:

```text
b_2       -> (0,0,2)  on (34,1,300),
b_4       -> (0,1,0)  on (34,1,300),
b_2 + b_4 -> (0,1,2)  on (34,1,300).
```

So `b_2+b_4` is an explicit word-length-one separator for the coarse
block, while `b_4` is the single direction that completes the algebra.
The character values are read from the published regular representation;
no new primary decomposition is claimed.

## Consequence for the twelve-site target

The twelve-site primary calculation should include polynomial coefficient
position `5` (or its twelve-site analogue) as the first extra direction,
rather than searching a two-dimensional complement. Whether the
twelve-site analogue still completes the algebra, and whether its words
split the twelve-site `R,T`-coarse blocks, is open. This note does not
construct the twelve-site model, prove a rank bound, or change the
factoring exponent; the complete bound remains `1/2`.

## Update: the twelve-site analogue

The twelve-site case is now checked in
[the separator-direction note](bounded-spin-twelve-missing-direction.md).
The ten-site completing basis direction `b_4` alone gives dimension `7`
there, but the ten-site separator sum `b_2+b_4` completes the
ten-dimensional twelve-site coefficient algebra.  The roles of the
completing direction and the separator therefore reverse; the exact
twelve-site primary decomposition and cut ranks remain open.

## Reproduction

```sh
python3 tools/check_bounded_spin_missing_direction.py
```

The checker rebuilds the ten-site coefficient tensor in about four
minutes, solves for the `R,T` coordinates, and verifies the generated
dimensions and the three separator values. It writes
`evidence/bounded-spin-missing-direction.json`.

The bounded-spin construction, path frames, normalization, and ten-site
algebra certificate are credited to the bounded-spin lane. This note
continues the [Discussion #51 exchange](https://github.com/1337isnot1337/integer-factoring/discussions/51).
