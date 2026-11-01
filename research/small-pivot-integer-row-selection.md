# Small internal pivots retain the same original relations

## Result

The [small-pivot selector](../tools/small_pivot_integer_row_basis.py)
changes the internal pivot order of the
[primitive integer row basis](primitive-integer-row-selection.md).
It chooses the smallest absolute nonzero coefficient, breaking ties
by column index, and eliminates in insertion order. It returns exactly
the same original-row subsequence as rational Gaussian rank selection.

Three paired measurements reduce rank-selection CPU time from 0.137
to 0.068 seconds on the retained 96-bit core, and from 16.918 to 7.002
seconds on the 112-bit core. This is an implementation improvement to
ordinary exact Gaussian selection, with no new factoring exponent.

## Exactness with arbitrary pivot columns

For stored residuals `v_1,...,v_r` and distinct pivots `p_1,...,p_r`,
maintain

$$v_i[p_i]\ne0,\qquad v_i[p_j]=0\quad(j<i).$$

The columns need not be numerically ordered. Process them in insertion
order. At step i, replace the current residual u by

$$v_i[p_i]u-u[p_i]v_i,$$

then divide a nonzero result by its positive content. Sign normalization
is also allowed. This zeros the current pivot. Later stored rows vanish
in every earlier pivot, so later operations do not restore one.

Each operation changes the coset modulo the stored rational row span
only by a nonzero scalar. Thus an input outside that span cannot reduce
to zero. An input inside the span reduces to a vector in that span
which is zero in every pivot column. The stored pivot submatrix is
triangular with nonzero diagonal, so this remaining vector is zero.

A nonzero residual may therefore choose any nonzero remaining column.
The implementation chooses minimum absolute value to reduce future
multipliers, and makes that pivot positive. Induction gives the same
acceptance decision for each input row as every exact rational rank
test. The accepted original-row subsequence is consequently unchanged.

Only temporary rank-test residuals are divided by their content.
On acceptance, `.rows` receives the untouched original input. The
relations `(2,0)` and `(0,6)` remain those relations; internal unit
pivots are not emitted as unsupported divided class relations.
This does not claim that the selected subset spans the full integer
relation lattice. The existing saturation and torsion-image argument
continues to apply to the same original subset.

## Costs and exact controls

The prior primitive-row bit bound remains valid. Schur residuals are
represented by minors using the chosen pivot columns in insertion
order. Hadamard bounds primitive residuals and temporary products by
`O(d(b+log(d+1)))` bits for d columns and b-bit input entries. Every
content gcd and exact division is charged. The algorithm uses
`O(m d^2)` integer operations for m rows when `d>=1`, together with
the content operations at those bounded lengths. Pivot choice adds
at most d comparisons per accepted row. This is not a stronger
worst-case matrix-rank complexity theorem.

The [checker](../tools/check_small_pivot_integer_row_basis.py) compares
680 matrices in dimensions zero through sixteen with exact rational
selection. All 17000 decisions and complete original-row subsequences
agree. It checks the triangular invariant after every input, plus
positive primitive pivots, nonmonotone pivot columns, the relation-content
trap, and three invalid inputs. The
[control receipt](../evidence/small-pivot-integer-row-controls.json)
retains those checks and source hashes.

The [paired fixture receipt](../evidence/small-pivot-integer-row-basis.json)
uses final matrices from the completely replayed 96-bit and 112-bit
results. Each fixture has three paired repetitions with alternating
method order. Every repetition checks the complete selected-row list.

| Core | Leftmost primitive CPU median | Small-pivot CPU median | Elimination steps, old/new |
| --- | ---: | ---: | ---: |
| 321 rows, 89 columns | 0.136763 s | 0.068277 s | 4045 / 3014 |
| 619 rows, 289 columns | 16.917654 s | 7.002244 s | 60174 / 36192 |

Wall medians are 0.142/0.069 and 17.526/7.568 seconds. Other research
processes were active. These measurements cover rank selection only;
complete search, repeated checkpoints, fixture loading, and independent
factor-certificate replay have separate costs.

```sh
python3 tools/check_small_pivot_integer_row_basis.py \
  --reference evidence/optimized-oriented-power-96.json.xz \
  --reference evidence/optimized-oriented-power-112.json.xz \
  --repetitions 3 --output /tmp/small-pivot-row-check.json
```

This implementation and analysis used substantial OpenAI assistant support.
