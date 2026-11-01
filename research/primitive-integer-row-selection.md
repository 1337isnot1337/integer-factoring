# Primitive temporary rows for exact relation selection

## Result

The [integer row selector](../tools/primitive_integer_row_basis.py) returns
exactly the same subsequence of original rows as the existing rational
Gaussian selector, for every ordered integer matrix. It performs its
rank tests with primitive integer residuals instead of a `Fraction` at
every matrix entry.

Only temporary rank-test vectors are divided by their content gcd. Every
accepted relation is copied from the original input unchanged. This
distinction preserves the integral relations and their torsion content.

On the actual 321-by-89 core from the certified 96-bit prefix, five
alternating matched runs give median CPU times of 0.6057 seconds for the
rational selector and 0.1407 seconds for the integer selector. The same
89 original rows are selected in the same order. This measures the
selection component, not a complete factoring run or a new uniform
factoring bound.

## 1. Exact original-row preservation

For an integer residual vector v, remove zero entries, divide its nonzero
entries by their positive content gcd, and choose the sign that makes
the first nonzero entry positive. Denote the resulting primitive vector
by `prim(v)`. A zero residual stays zero.

Store one such pivot vector V for each pivot column p. Its pivot entry
`a=V_p` is a positive integer. To eliminate the p-entry `b=v_p` of a new
temporary residual, use

$$v\longmapsto\mathrm{prim}(a v-b V). \tag{1}$$

All operations in (1) are exact integer operations. After processing the
existing pivots in increasing column order, a nonzero residual creates
the next pivot. The selector then retains the untouched original row.

Compare this with normalized rational Gaussian elimination. Suppose its
pivot is P with `P_p=1`, and its current residual is R. Inductively, the
integer pivot and residual satisfy `V=aP` and `v=cR` for some nonzero
rational c. Before taking primitive content, (1) gives

$$a(cR)-(cR_p)(aP)=ac(R-R_pP). \tag{2}$$

Primitive normalization changes only a nonzero rational scale. Thus the
two procedures have the same zero residuals and the same first nonzero
pivot columns at every step. Induction over the ordered input rows proves
that they accept exactly the same original-row subsequence.

For example, the genuine relations `(2,0)` and `(0,6)` may have temporary
pivots `(1,0)` and `(0,1)`. The output remains `(2,0),(0,6)`. Emitting the
primitive pivots instead would alter the integer relation lattice and
could destroy the order-two information needed for factoring. The checker
includes this explicit control.

## 2. Coefficient sizes and charged cost

Let the matrix have d columns, and let every original entry have magnitude
at most `2^b`. After eliminating k independent pivot rows, each remaining
rational coordinate is a Schur-complement entry. By the determinant
formula, it is a `(k+1)`-row minor of the original matrix divided by the
same nonzero k-row pivot minor.

The vector of those numerator minors is integral. Its primitive integer
representative is the stored residual, up to sign. Content division can
only reduce the numerator magnitudes. Hadamard's inequality therefore
bounds every normalized residual entry by

$$\left(\sqrt{k+1}\,2^b\right)^{k+1}. \tag{3}$$

The temporary cross-products before content normalization multiply two
such bounded vectors, so they also have
`O(d (b+log(d+1)))` bits. This prevents the coefficient doubling that an
unchecked repeated cross-multiplication could produce.

For m supplied rows and `d>=1`, the direct implementation uses `O(m d^2)` integer
arithmetic operations: cross-products, subtractions, content gcd steps,
and exact content divisions. Their operand sizes are bounded as above.
All of those gcds and divisions are charged. Internal pivot storage uses
at most d sparse vectors; a dense upper bound is
`O(d^3 (b+log(d+1)))` bits, in addition to the input representation.
The returned original rows contribute at most `O(d^2 b)` bits.

These are ordinary polynomial component bounds. The implementation does
not claim a new asymptotic matrix-rank algorithm. Its benefit here is the
measured reduction in rational-object and per-entry normalization work.

## 3. Exact controls and actual-core measurement

The [checker](../tools/check_primitive_integer_row_basis.py) compares 17000
individual acceptance decisions across 680 ordered integer matrices.
Widths range from zero through sixteen. The matrices include zero rows,
signed entries, nonprimitive rows, deliberately dependent combinations,
and changing pivot-column patterns. Every accepted original-row list
agrees after every insertion. Every stored nonzero pivot has content one.

The larger comparison uses all 321 rows and 89 columns of the retained
[96-bit core](../evidence/optimized-oriented-power-96.json.xz). It repeats
each method five times with alternating order:

| Measurement | Rational selector | Primitive integer selector |
| --- | ---: | ---: |
| Median wall time | 1.366993 s | 0.310237 s |
| Median CPU time | 0.605675 s | 0.140686 s |
| Selected original rows | 89 | 89 |

The integer run performs 4045 elimination steps and 3843 nontrivial
content normalizations. Its largest observed temporary entry has 1109
bits. The selected original rows have canonical SHA-256
`b250bce2a302bc5ec32c32930e8d8cbc729139b9fa7628c893c85c5599af77b3`.

The [receipt](../evidence/primitive-integer-row-basis.json) retains all
timings, the input certificate hash, selection hash, and source guards.

```sh
python3 tools/check_primitive_integer_row_basis.py \
  --output evidence/primitive-integer-row-basis.json
```

This implementation and analysis used substantial OpenAI assistant support.
