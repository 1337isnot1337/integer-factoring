# Diagonal evaluation weights of the bounded-spin coefficient span

## Result and scope

The [six-vertex coefficient span](bounded-spin-six-vertex.md) acts on
the product cut states `(a,b)` by a scalar diagonal evaluation when the
input and output path pairs agree. Compressing those evaluations gives
a small list of projective linear forms in the coefficient coordinates.
This note records the exact finite lists at lengths eight, ten and
twelve. It does not identify a cut rank with a form multiplicity and
does not prove a bound.

For each even length `n`, let `L_n` be the number of return paths and
let `m_n=(n/2+1)(n/2+2)/2` be the number of six-vertex monomials. For a
pair `(a,b)` of paths, evaluate every monomial at `(a,b,a,b)` and
projectivize the resulting vector in `F_5^{m_n}`. The distinct
projective vectors and their pair multiplicities are:

| `n` | `L_n` | distinct forms | square-root multiplicity sum | target `L_n^sqrt(2)` |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 13 | 13 | 42.998532 | 37.614370 |
| 10 | 34 | 35 | 175.529811 | 146.500213 |
| 12 | 89 | 124 | 808.241885 | 571.286269 |

The square-root sum is shown only because it is the same functional
that appears in a two-cut rank sum. No theorem here equates it with a
cut-rank lower bound. It is a compact diagnostic for the next primary
decomposition: the forms are the diagonal weights of the coefficient
algebra on the product-state module, and their multiplicities are the
sizes of the corresponding evaluation classes.

Restricting the coefficient space to the original pencil span can only
merge projective form classes. Since `sqrt(a)+sqrt(b) <= sqrt(a+b)`,
the full six-vertex square-root sums above are lower bounds for the
same diagnostic on any subspace, including the original pencil span.

## What the numbers say

The number of distinct forms grows from `13` to `35` to `124`, roughly
`L_n^1.4`, while the target grows as `L_n^sqrt(2)`. The square-root sum
is above the target at all three lengths and its ratio grows:
`1.14`, `1.20`, `1.41`.

A direct rank interpretation is nevertheless false. At length eight
the actual cut ranks can exceed or fall below the number of nonzero
diagonal evaluations. Deterministic examples include rank pair
`(169,169)` with diagonal count `161`, rank pair `(133,169)` with count
`165`, and rank pair `(126,92)` with count `169`. Thus the diagonal
count is neither an upper nor a lower rank bound. The table is only a
reproducible weight-distribution diagnostic; a rank theorem would need
additional structure.

The diagonal forms are also a candidate replacement for the full
ten-site tensor construction: there are at most `L_n^2` of them, and
each is obtained from the published monomial-entry constructor in
polynomial time. The twelve-site run takes about two and a half minutes
in this checkout.

## Reproduction

```sh
python3 tools/check_bounded_spin_diagonal_weights.py --write-evidence
```

The script uses only the standard library and the published path and
monomial-entry constructors. The rank control at length eight is

```sh
python3 tools/check_bounded_spin_diagonal_rank_control.py
```

It builds the actual cut operators and asserts that both inequalities
fail in a deterministic sample. It writes
[`evidence/bounded-spin-diagonal-weights.json`](../evidence/bounded-spin-diagonal-weights.json)
with the lengths, monomial order, distinct-form counts, multiplicities
and square-root sums. It does not build a cut operator or a tensor.
No factoring exponent changes; the complete bound remains `1/2`.
