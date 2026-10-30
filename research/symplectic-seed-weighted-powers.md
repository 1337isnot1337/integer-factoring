# Weighted radical-seed powers and exact small modules

## A larger common-parameter family

For the 26-label radical seed, remove its known common scalar factor
\(u=q+r\). The resulting weighted cubic is

\[
K(u,v,s)=u^3E+u^2vA+uv^2B+v^3C+usD+vsF,
\qquad v=p+q,\quad s=qr.
\]

Give \(u,v\) weight 1 and \(s\) weight 2. The crossing is
\(K(1,0,0)=E\). This family includes the earlier curve \(R(t,0)\), but also
keeps the two nilpotent directions \(D,F\).

The complete six-generator action admits constant direct-sum decompositions
into **11 row types of size at most 28 and 15 column types of size at most 8**.
The delivered basis certificate verifies all six generators, including
\(D,F\). It contains 8,112 exact intertwining-column identities and full-rank
basis checks in every original block. The source is the frozen adapted
representation in `evidence/symplectic-radical-seed-algebra.json`.

| Summand dimension | Row multiplicity | Column multiplicity |
|---:|---:|---:|
| 1 | 24 | 26 |
| 2 | 96 | 153 |
| 3 | 32 | 0 |
| 4 | 0 | 33 |
| 6 | 8 | 34 |
| 8 | 0 | 1 |
| 12 | 24 | 0 |
| 28 | 1 | 0 |

Equal dimensions can represent distinct module types. The certificate retains
their exact matrices separately. Each cut has total dimension 676.

The proof data and evaluator are:

- `evidence/symplectic-seed-weighted-modules-certificate.json`;
- `evidence/symplectic-seed-weighted-modules.json`;
- `tools/symplectic_seed_weighted_modules.py`.

## Complete weighted coefficient spaces

For \(n\) tensor copies with common parameters, the possible monomials are

\[
u^{3n-a-2b}v^a s^b,
\qquad 0\leq b\leq n,\quad 0\leq a\leq 3n-2b.
\]

The evaluator records each monomial by its pair of exponents \((a,b)\). There are

\[
\sum_{s=0}^{n}(3n-2s+1)=(n+1)(2n+1)
\]

coordinates: 15 for two copies and 28 for three copies. The constant crossing
coordinate is \((0,0)\), namely \(u^{3n}\).

The evaluator groups identical small-module tensor factors by type counts.
It multiplies their polynomial matrices with common parameters, merges equal
monomials, and only then computes rank. The integer multiplicities are the
usual multinomial multiplicities of a direct-sum tensor power.

For two copies, all 66 row type blocks and 120 column type blocks are retained.
Their maximum dimensions are 784 and 64. For three copies, there are 286 and 680
type blocks. The checked cube bound retains the 127 row and 519 column blocks
of size at most 128. Omitted blocks are replaced by explicit zero summands of
the same weighted dimension. Therefore the resulting ranks are rigorous
**lower bounds** on the original coefficient ranks.

The omitted weighted dimensions in the cube are 171,968,320 on the row cut
and 27,325,952 on the column cut. The full cut dimension is
\(17576^2=308915776\). Omitting these summands weakens the exclusion test; it
cannot create a false exclusion.

## Exact finite gates

For label dimension \(d\), choose the integer

\[
g=\left\lceil d^{99/70}\right\rceil.
\]

The integer inequality \(99^2>2\cdot 70^2\) proves that
\(g\geq d^{\sqrt 2}\). A coefficient combination is eligible if its rank
product is less than \(g^2\).

If a decomposition had total cost below \(d^{\sqrt 2}\), every one of its
pieces would be eligible. Consequently, if the eligible vectors do not span
the crossing coordinate, no such decomposition exists. The same argument
works with lower bounds on the ranks: the computed eligible set is then a
superset of the truly affordable set.

The complete binary results are:

| Copies | Labels | Coordinates | Vectors tested | Gate | Eligible vectors | Eligible span dimension | Crossing in span |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| 2 | 676 | 15 | 32,767 | 10,055 | 14 | 7 | No |
| 3 | 17,576 | 28 | 268,435,455 | 1,008,161 | 2,073 | 17 | No |

The square uses exact ranks. All 14 eligible vectors and 40 additional vectors
were checked by independent scalar elimination. The cube uses the retained
direct-summand rank lower bounds. All 2,073 eligible vectors and 40 additional
vectors were independently checked. Its 32 consecutive scan chunks cover the
entire nonzero 28-bit coordinate space without gaps or overlaps.

Thus neither the square nor the cube of this particular 26-label seed passes
the binary finite-seed criterion anywhere in its complete common-parameter
weighted coefficient space. This is a larger exclusion than checking a
single chosen curve. It is not a statement about different label codes,
extension-field weighted spaces, or all tensor powers.

## Why reparameterization does not enlarge a fixed coefficient span

Let \(M(t)=\sum_{i=0}^m M_i t^i\). For a nonconstant rational function
\(a(t)/b(t)\), clearing denominators gives

\[
b(t)^m M(a(t)/b(t))
=\sum_{i=0}^m M_i a(t)^i b(t)^{m-i}.
\]

The scalar polynomials on the right are linearly independent over the ground
field: otherwise a nonzero constant-coefficient polynomial would annihilate
the nonconstant rational function. Therefore the new matrix coefficients
span precisely the original coefficient space. Multiplying or removing a
nonzero common scalar polynomial also preserves that space.

More generally, substituting any common rational functions into
\(K(u,v,s)^{\otimes n}\) and clearing denominators produces coefficients in
its complete weighted coefficient span. Over the binary field, the square
and cube gates above therefore also cover these common scalar-weighted and
rational curves. A change of coefficient basis cannot evade a complete-span
test.

## Different local curves

Different curves in different tensor factors need not lie in the common
weighted coefficient space. This is a genuine additional construction.
The checked cubic family is

\[
K(1+at,t,bt+ct^2),\qquad a,b,c\in\mathbb F_2.
\]

Each of its eight curves evaluates to the crossing at \(t=0\). The checker
tests all 28 distinct unordered pairs of curves. For each pair it computes
both exact rank tables for all 127 nonzero merged coefficient combinations
and solves the complete XOR shortest-path problem to the crossing.

None passes. The best checked cost lies in the integer interval
\([225868,225869]\), for the pair of parameter triples
\((0,0,1),(0,1,0)\), and also for its reversed-coordinate counterpart. Its
numerical cost is 225868.75793635. Here the label dimension remains
\(d=676\). Every one of the 28 integer lower bounds satisfies the stronger
comparison \(S_{\rm lower}^2>676^3\).

The complete mixed-curve evidence is
`evidence/symplectic-seed-mixed-curves.json`, reproduced by
`tools/check_symplectic_seed_mixed_curves.py`. These finite results do not
exclude higher-degree mixed curves or different local label spaces.

## Reproduction and native optimization

The full module certificate can be replayed directly:

```sh
python tools/symplectic_seed_weighted_modules.py \
  --output /tmp/symplectic-seed-weighted-modules.json
```

Generate the exact square input and replay its gate:

```sh
python tools/symplectic_seed_weighted_modules.py \
  --native /tmp/weighted-square.in --power 2
g++ -O3 -DNDEBUG -std=c++20 tools/screen_weighted_pfaffian_span.cpp \
  -o /tmp/weighted-screen
/tmp/weighted-screen --field2 < /tmp/weighted-square.in > /tmp/weighted-square.jsonl
python tools/check_weighted_pfaffian_screen.py \
  --input /tmp/weighted-square.in --screen /tmp/weighted-square.jsonl \
  --output /tmp/weighted-square-audit.json
```

Generate the capped cube input with `--power 3 --maximum-block-size 128`.
The separate `tools/screen_weighted_pfaffian_wide.cpp` accepts the 28 binary
coordinates. It leaves the frozen field-screen source unchanged. Its exact
fast path precomputes ranks of all binary matrices of size at most 4 and uses
bytewise coefficient maps for those blocks. Larger blocks use the original
elimination. The coefficient bit spreading is also performed by a fixed
integer operation sequence.

The wide tool accepts `--start` and `--count`. Consecutive chunk logs can be
concatenated and audited by `tools/check_weighted_pfaffian_chunks.py`. The
delivered cube log is
`evidence/symplectic-weighted-cube-cap128-F2.jsonl`; its deterministic export
metadata and independent audit are the corresponding `-export.json` and
`-audit.json` files. The audit checks source hashes, full interval coverage,
all eligible ranks, the integer gate, and the exact endpoint span remainder.
