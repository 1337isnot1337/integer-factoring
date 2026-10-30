# Exact field screens for shared-parameter seed powers

The 26-label radical seed has normalized polynomial
\[
 R(t)=E+tA+t^2B+t^3C.
\]
The [constant module decomposition](symplectic-seed-tensor-modules.md)
reduces its row cut to ten block types of size at most 22 and its column
cut to nine types of size at most 6. This note uses those exact
decompositions to screen the entire coefficient spaces of
`R(t)^{⊗n}` for `n=2,3,4` over both `F₂` and `F₄`.

All six spaces fail the finite seed criterion. For the square, the ranks
are complete. For the cube and fourth power, the screens use certified
rank lower bounds from smaller direct summands. The affordable
supersets already fail to span the crossing, so no ranks of omitted
large blocks are needed.

These are finite exclusions for shared-parameter powers of this curve.
They do not exclude other codes, nonzero second-parameter curves, or
different curves at different tensor factors.

## 1. Exact weighted blocks

In one cut, write the certified constant decomposition as
\[
 R(t)\cong\bigoplus_i R_i(t)^{\oplus m_i}.
\]
In the `n`th tensor power, a multiset of block types with counts `c_i`
contributes the polynomial matrix
\[
 \bigotimes_i R_i(t)^{\otimes c_i}
\]
with integer multiplicity
\[
 \frac{n!}{\prod_i c_i!}\prod_i m_i^{c_i}.
 \tag{1}
\]
Permuting factors only permutes matrix rows and columns, so one block
per multiset suffices. The multiplicities in (1) count direct summands;
they are never reduced modulo two.

The two cuts are decomposed independently. Their block sizes and
multiplicities can differ. Each rank is the sum of its weighted block
ranks for the same polynomial coefficient functional.

For `n=2`, the evaluator retains every block. The row and column cuts
have 55 and 45 weighted blocks, with largest sizes 484 and 36. All 127
binary combination rank pairs exactly match the independent Python
tensor-module checker.

For `n=3,4`, retain only blocks of matrix size at most 64. The remaining
blocks are replaced by zero matrices in the native input. This retains
the full dimension bookkeeping while explicitly discarding positive
rank contributions. For every coefficient vector `a`,
\[
 \underline r(a)\le r(a),\qquad
 \underline c(a)\le c(a).
 \tag{2}
\]
The retained dimensions are:

| Power | Labels `d=26ⁿ` | Full cut dimension `d²` | Retained row dimension | Retained column dimension |
|---:|---:|---:|---:|---:|
| 2 | 676 | 456,976 | 456,976 | 456,976 |
| 3 | 17,576 | 308,915,776 | 243,847,616 | 297,674,560 |
| 4 | 456,976 | 208,827,064,576 | 115,283,166,032 | 152,578,333,184 |

The source module certificate contains explicit constant basis columns.
Before building the tensor blocks, the evaluator checks all 5408
intertwining column identities and the source hashes.

## 2. Exhaustive coefficient-space exclusions

The degree is `3n`, so the coefficient coordinate count is `3n+1`.
The crossing is the coefficient at degree zero.

For an exact integer gate `g`, the evaluator verifies
\[
 g^{70}>d^{99},\qquad 99^2>2\cdot70^2.
\]
Thus `g>d^{√2}`. Every piece in a winning decomposition must have
cost less than `g`. By (2), its coefficient vector must appear in the
screen's list satisfying
\[
 \underline r(a)\,\underline c(a)<g^2.
 \tag{3}
\]
If the span of that list misses the crossing, no winning decomposition
exists. The same proof applies when the retained blocks give strict
lower bounds.

The coefficient vectors need not be assumed independent for this
argument. Every true linear relation gives a zero tensor and therefore
zero retained ranks, so its coefficient vector is included among the
eligible vectors. Such relations cannot invalidate the span test.

| Power | Field | Coordinates | Projective vectors tested | Eligible vectors | Eligible span dimension | Integer gate |
|---:|---|---:|---:|---:|---:|---:|
| 2 | `F₂` | 7 | 127 | 3 | 2 | 10,055 |
| 2 | `F₄` | 7 | 5,461 | 5 | 2 | 10,055 |
| 3 | `F₂` | 10 | 1,023 | 8 | 4 | 1,008,161 |
| 3 | `F₄` | 10 | 349,525 | 22 | 4 | 1,008,161 |
| 4 | `F₂` | 13 | 8,191 | 32 | 6 | 101,089,575 |
| 4 | `F₄` | 13 | 22,369,621 | 342 | 6 | 101,089,575 |

Every row of the table has a nonzero crossing remainder. In the
four-element field, projective enumeration divides the nonzero
coefficient vectors by the three possible nonzero scalars. Matrix ranks
are computed over `F₄=F₂[ω]/(ω²+ω+1)`.

## 3. Independent rank audits

The native evaluator uses high-pivot elimination in two binary bit
planes. A separate Python implementation uses low pivots and direct
field arithmetic in those planes. It checks every eligible vector and
forty additional deterministic samples from each field screen.

| Power | Binary vectors independently checked | Four-element-field vectors independently checked |
|---:|---:|---:|
| 2 | 43, plus all 127 ranks compared with the module checker | 45 |
| 3 | 48 | 62 |
| 4 | 72 | 382 |

For powers three and four, these independent calculations verify the
retained-block ranks used in (3). They do not claim that the omitted
blocks vanish. The complete native logs and every eligible coefficient
vector are preserved.

## 4. Scalar weighting and rational changes of parameter

The full coefficient-space exclusions also cover scalar weighting and
nonconstant rational changes of parameter of the same curve.

**Lemma.** Let `F(t)=Σ_{k=0}^m T_k t^k` be a polynomial with coefficients
in a finite-dimensional vector space over a field. Multiplication by
a nonzero scalar polynomial preserves the span of its coefficients.
If `a(x)/b(x)` is nonconstant, then clearing denominators in
`F(a(x)/b(x))` also preserves that span.

**Proof.** Every new coefficient is a linear combination of the old
ones. For scalar multiplication, division by the first nonzero scalar
coefficient recovers the old coefficients recursively. For rational
substitution, the scalar polynomials
\[
 a(x)^k b(x)^{m-k},\qquad 0\le k\le m,
\]
are linearly independent: a relation would make the nonconstant
rational function `a/b` algebraic over the constant field. Applying
any linear functional to the tensor coefficients then shows that the
old and new spans have the same annihilator. Their spans agree.
\(\square\)

The lemma applies directly to `F(t)=R(t)^{⊗n}`. A new coefficient space
therefore requires a construction change, such as activating the
second parameter in the full family or assigning different curves to
different tensor factors.

## 5. Reproduction and artifacts

Compile the weighted native evaluator:

```sh
c++ -O3 -std=c++20 tools/screen_weighted_pfaffian_span.cpp -o /tmp/weighted-pfaffian-screen
```

Replay all six screens and their independent audits:

```sh
python3 tools/check_symplectic_seed_tensor_fields.py --types evidence/symplectic-seed-tensor-modules-certificate.json --power 2 --native /tmp/weighted-pfaffian-screen --compare-binary evidence/symplectic-seed-tensor-modules.json --export-dir /tmp/omega-power-2 --output /tmp/omega-power-2-audit.json
python3 tools/check_symplectic_seed_tensor_fields.py --types evidence/symplectic-seed-tensor-modules-certificate.json --power 3 --maximum-block 64 --native /tmp/weighted-pfaffian-screen --export-dir /tmp/omega-power-3 --output /tmp/omega-power-3-audit.json
python3 tools/check_symplectic_seed_tensor_fields.py --types evidence/symplectic-seed-tensor-modules-certificate.json --power 4 --maximum-block 64 --native /tmp/weighted-pfaffian-screen --export-dir /tmp/omega-power-4 --output /tmp/omega-power-4-audit.json
```

The fourth-power extension-field enumeration checks more than 22 million
projective vectors. It is an exhaustive computation, even though the
small blocks reject most vectors early.

Durable files:

- [Weighted native evaluator](../tools/screen_weighted_pfaffian_span.cpp).
- [Tensor-block exporter and independent field auditor](../tools/check_symplectic_seed_tensor_fields.py).
- [Square evidence](../evidence/symplectic-seed-power-2-field-audit.json).
- [Cube evidence](../evidence/symplectic-seed-power-3-field-audit.json).
- [Fourth-power evidence](../evidence/symplectic-seed-power-4-field-audit.json).
- [Six complete native logs](../evidence/symplectic-seed-field-screens).

Every report records exact source, module-certificate, native-input, and
screen hashes. The native evaluator reuses the audited scalar field
elimination from `screen_pfaffian_span_f4.cpp`; its weighted layout and
independent two-cut block lists are the new parts.
