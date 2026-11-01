# Odd relation pivots with exact integer identities

## Result and scope

A fresh generator can be eliminated when its coefficient in a genuine
relation is any odd integer. The elimination preserves the image of the
presented two-primary torsion. It also produces genuine integer relations
for the existing [determinant-free extractor](determinant-free-ambiguous-extraction.md).
No class-group order, annihilator, or division in the group is needed.

The essential implementation detail is an unreduced common odd multiplier.
Cancelling that multiplier as one would simplify a rational expression can
turn a genuine relation into a false one in the odd-order part of the group.
The algorithm below keeps exact identities in the full group throughout.

This is a presentation compression theorem, not a bound on how quickly
useful relations can be selected. It leaves the uniform closure and
short-power-window questions in
[dynamic leading-norm closure](dynamic-leading-norm-closure.md) open.

## 1. An exact common multiplier

Write the abelian group additively. The active labels name actual elements
of the group. For every known label `p`, maintain an integer word `V_p`
in the active labels and one common positive odd integer `D`, with

\[
 Dg_p=\phi(V_p).
 \tag{1}
\]

Initially `D=1` and `V_p` is the corresponding unit vector. When a new
label is retained as active, its expression is `D` times its unit vector.
Stored relation rows are integer vectors that evaluate to zero under
`phi`; they do not need a common numerical scale.

Suppose a genuine raw row first introduces a label `q` with odd nonzero
coefficient `c`:

\[
 c g_q+\sum_{p\ne q}r_pg_p=0.
 \tag{2}
\]

All other fresh labels in this row first become active. Using the current
expressions, form

\[
 W=\sum_{p\ne q}r_pV_p.
 \tag{3}
\]

Then `c D g_q=-phi(W)` exactly. Set

\[
 D'=|c|D,\qquad
 V'_p=|c|V_p\quad(p\ne q),\qquad
 V'_q=-\mathrm{sign}(c)W.
 \tag{4}
\]

Delete the defining row and the fresh column `q` from the retained
presentation. Equation (1) remains exact for every known label. Existing
stored rows remain genuine without rescaling: the active group elements
themselves have not changed.

For a later genuine raw relation `r` with no chosen fresh odd pivot, store

\[
 R=\sum_p r_pV_p.
 \tag{5}
\]

Its group value is exactly `D` times the value of the raw relation, hence
zero. Therefore the relation-halving extractor can verify and use this
integer row without an odd-order error term. The implementation chooses
the largest fresh integer label with odd coefficient; this choice is a
heuristic and is not needed for correctness.

## 2. Why the common multiplier must not be cancelled

Take `G=Z/6Z`, an active element `g=1`, and a fresh element `q=4`.
The two independent raw rows

\[
 3q=0,\qquad 2g+q=0
 \tag{6}
\]

are genuine. Eliminating `q` gives `D=3`, `V_q=0`, and `V_g=3g`.
The second row then gives the genuine relation `6g=0`. Its half is the
actual involution `3g`.

Replacing `3q=0` by the rational-looking assignment `q=0` would instead
produce `2g=0`, which is false. Its proposed half `g` has order six.
In particular, even the expression `V_q=0` does not permit cancelling
the common odd multiplier. The engine never performs this cancellation.

## 3. The preserved object

Let `L` be the lattice spanned by the currently available raw relations
on all introduced labels, and let

\[
 A=\mathbb Z^n/L.
 \tag{7}
\]

The raw group values define a homomorphism `A -> G`. Let `A'` be the
corresponding quotient for the compressed active presentation. Put
`R=Z_(2)`, the rationals with odd denominator. Every pivot `c` and every
stored-row multiplier `D` is a unit of `R`. Thus (4) is ordinary pivot
elimination over `R`, and (5) changes a row only by a unit after
substitution. Induction on the stream gives an isomorphism

\[
 A\otimes R\;\simeq\;A'\otimes R
 \tag{8}
\]

that commutes with the maps to `G tensor R`. For finite `G`, the latter
is its two-primary subgroup. More generally, the same assertion applies
to the two-primary torsion of the finitely generated image involved here.

The torsion part of a finitely generated abelian group after tensoring
with `R` is exactly its two-primary torsion. Therefore (8) preserves the
two-primary part of the image of the presented torsion in `G`, and hence
the useful involutions in that image. In a rectangular presentation this
is only the torsion image: the free part can still map onto additional
two-torsion, just as in the existing partial-relation theorem.

All rows supplied to the extractor are genuine integer relations by
(5). Thus every evaluated half-relation is an actual involution in the
full group, not merely after localization. The localized proof establishes
preservation of available information; it is not a license to accept a
group equality only modulo odd-order elements.

If the raw rows are rationally independent, an elimination removes one
fresh pivot column and exactly one independent row. Other substitutions
and nonzero odd row scalings preserve independence. The compressed rows
therefore satisfy the rank promise of the generic extractor.

When a full finite presentation is available, the active subgroup also
contains the full two-primary part of the generated group: (1) says that
the quotient by the active subgroup is annihilated by the odd integer
`D`. This observation alone does not replace the partial-presentation
torsion-image statement above.

## 4. Size and the leading-norm application

Suppose `s` raw rows have integer `l1` norm at most `B`, with `B>=2`.
At every step, the joint maximum of 1, the common multiplier, the largest
expression norm, and the largest stored-row norm increases by at most a
factor `B`. Consequently

\[
 \log_2 D\le s\log_2B,\qquad
 \max_p\log_2(1+\|V_p\|_1),\;
 \max_R\log_2(1+\|R\|_1)=O(s\log B).
 \tag{9}
\]

The number of active labels is at most the number of written labels in
the input stream. Sparse integer substitution and the existing extractor
therefore have polynomial cost in the written stream size, group-operation
cost, and coefficient bit lengths. The claim does not require writing a
group power in unary, and does not treat a potentially large exponent as
a unit-cost integer.

For a reduced leading norm `A<=H`, any fresh prime larger than `sqrt(H)`
has exponent one and is the only prime of that size in `A`. It is
therefore the largest fresh odd pivot and is eliminated by this rule.
Every newly retained active prime is at most `sqrt(H)`, as with the
existing unit-pivot rule. Odd multiplicities such as a fresh `q^3` can
now also be eliminated. This can reduce the active presentation, but
does not prove any upper bound on the number of selected norms or on
the remaining deficiency.

## 5. Implementation and exact controls

The engine is
[odd_relation_presentation.py](../tools/odd_relation_presentation.py).
Its caller supplies genuine, rationally independent raw rows; the engine
does not infer group equality or check rational rank. It exposes the
active labels, exact common multiplier, integer expressions, retained
matrix, and elimination history. It does not change the stable dynamic
producer or any previous norm consumer.

The [checker](../tools/check_odd_pivot_two_primary_elimination.py) verifies
the full-group invariant after each update. Small product groups have an
independently constructed saturated relation lattice, so their complete
torsion-image involution sets can be compared before and after compression.
It also checks negative and higher odd pivots, deep multiplier growth,
the order-six counterexample, and every stored dynamic history.

```sh
python3 tools/check_odd_pivot_two_primary_elimination.py \
  --output evidence/odd-pivot-two-primary-elimination.json
```

The [evidence](../evidence/odd-pivot-two-primary-elimination.json) records
the exact control counts and source hashes. The existing dynamic histories
use only unit fresh odd pivots; dedicated synthetic controls exercise
the new odd coefficients. No finite control is used as a relation-selection
or factoring-exponent theorem.
