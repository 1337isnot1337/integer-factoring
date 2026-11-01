# Integer graph elimination preserves residual torsion

## Result and scope

The [oriented coprime refinement](oriented-coprime-basis.md) can turn one
composite residual ideal into several signed atom generators. Most words
in the retained finite controls have total absolute exponent at most two.
Those rows admit a graph reduction over the integers.

The reduction eliminates tree definitions and anchored components. It
keeps one double-root relation for every unanchored unbalanced component.
Larger residual words can be carried through the same substitution without
discarding them. The full result is an isomorphic integer presentation.
Removing only unused free direct summands then preserves its evaluated
torsion image exactly.

On the retained 40-bit prefix, the complete 320-row system on 252
generators becomes five rows on five generators. Two additional unit
pivots reduce it further to three rows on the original three base
generators. Both relation-halving extractions return 1000003. This is a
finite matrix reduction and a
supplied-witness theorem. It gives no uniform search bound and does not
improve the complete factoring score of one half.

Signed incidence matrices, switching, and the balanced/unbalanced
distinction are classical. See Thomas Zaslavsky,
[*Signed graphs*](https://doi.org/10.1016/0166-218X(82)90033-6),
Discrete Applied Mathematics 4 (1982), 47--74, and his
[*Matrices in the Theory of Signed Simple Graphs*](https://arxiv.org/abs/1303.3083).
The argument below is self-contained. The implementation contribution is
the exact presentation maps, their relation certificates, retention of
useful residual torsion, and the connection to unfactored ideal atoms.
No literature-wide novelty claim is made for graph elimination.

## 1. The integer presentation

Let the fixed-base generators be `x_1,...,x_b` and the residual generators
be `y_1,...,y_v`. All evaluate in a common abelian group G. The supplied
rows are genuine relations, and L is their integer row lattice.

A graph row has the form

```math
\beta_e+\sigma_{e,u}y_u+\sigma_{e,v}y_v=0,
\qquad \sigma_{e,u},\sigma_{e,v}\in\{-1,1\},          \tag{1}
```

where beta is an integer word in the fixed base. Also permit rows with
one residual coefficient `+/-1`, one coefficient `+/-2`, or no residual
coefficient. Equivalently, the residual absolute weight is at most two.
The first kind is an anchor; the second is an unbalanced loop. Zero
residual words are already base relations.

The theorem permits arbitrary additional rows. They are transported
through the graph reduction and remain charged in its input and output.
They need not satisfy the graph contract.

Write

```math
P=\mathbb Z^{b+v}/L.                              \tag{2}
```

We construct an isomorphic presentation P' using fewer generators, with
explicit mutually inverse maps. An image in G is called useful only
when the existing exact group extractor verifies it and obtains a proper
factor. A graph cycle by itself is not a factoring certificate.

## 2. Eliminate a spanning forest using unit coefficients

Use the two-vertex graph rows to choose a spanning forest. Choose one
root r in each component. For every vertex, compute

```math
y_v=\epsilon_v y_r+b_v \pmod L,
\qquad \epsilon_v\in\{-1,1\},                    \tag{3}
```

where b_v is a fixed-base word. Set epsilon_r=1 and b_r=0. If a tree row
joins parent u to child v, solving its unit coefficient gives

```math
\epsilon_v=-\sigma_{e,v}\sigma_{e,u}\epsilon_u,
\qquad
b_v=-\sigma_{e,v}\beta_e
    -\sigma_{e,v}\sigma_{e,u}b_u.                 \tag{4}
```

No division by a nonunit occurs. Each eliminated generator has an exact
integer certificate

```math
y_v-\epsilon_vy_r-b_v=\sum_e d_{v,e}R_e,          \tag{5}
```

where R_e is an original row. The implementation retains these sparse
source combinations. Thus (3) is not inferred from a rational nullspace.

After substitution, every graph row in that component becomes

```math
q_e+t_ey_r=0,\qquad t_e\in\{-2,-1,0,1,2\}.        \tag{6}
```

Only anchors have `t_e=+/-1`. A tree row becomes zero. A two-vertex row
or an unbalanced loop has residual coefficient zero or `+/-2`.

## 3. Anchored, balanced, and unbalanced components

If an anchor exists, its coefficient in (6) is a unit. Solve it for y_r,
substitute it in every row, and remove that root. The entire component
is now expressed in the fixed base. Multiple anchors give additional
base relations after the first is used.

Suppose no anchor exists. If all graph coefficients t_e vanish, call
the component balanced. Its graph rows give only base relations. If any
t_e is nonzero, choose one such row, multiply it by `+/-1`, and retain

```math
2y_r+q_r=0.                                       \tag{7}
```

Subtract an integer multiple of (7) from each other graph row to remove
its coefficient `+/-2`. This gives base relations and exactly one
retained root relation for that unbalanced component.

The root in (7) must remain. Dividing the relation by two can change the
presentation and destroy its useful order-two image. Even when q_r is
even, the element `y_r+q_r/2` may be a useful involution. When q_r has an
odd coefficient, (7) need not create torsion at all; the theorem does not
assign an independent factor of order two to every component in the
presence of arbitrary base words.

With all base words zero and no additional rows, the classification is
exact: every balanced unanchored component contributes a free Z, every
unbalanced unanchored component contributes Z/2, and an anchored component
contributes zero. This special case is independently checked by counting
characters into Z/3 and Z/4.

### Carry larger residual words without loss

Apply (3) and every anchor substitution to each additional row. It can
involve several roots and arbitrary integer coefficients. For each root
relation (7), write its coefficient `c=2a+r`, with `r` equal to zero or
one. Subtract `a` copies of (7), retaining the remainder coefficient.
This is an integer row operation even when c is odd or negative.

These extra rows can couple previously balanced roots. Such a root is
not free merely because its graph component was balanced. Only root
columns that are zero in every final row are removed as free summands.
All other roots remain available to the exact torsion extractor.

### Additional unit pivots

The [unit-elimination pass](../tools/signed_residual_unit_elimination.py)
can remove any remaining residual column having a coefficient `+/-1`
in a current row. Solve that row for its residual generator and substitute
it in every other row and every generator definition. Delete the defining
row and generator. The pass retains the fixed-base columns throughout.

Every step is an integer presentation isomorphism. If a double-root
relation also uses the eliminated generator, substitution transfers that
relation to the remaining generators; it is not discarded or divided by
two. An isolated relation `2y=0` therefore survives this pass.

The implementation chooses a unit in a row of minimum current support,
with deterministic index tie breaking. This is a finite sparsity heuristic,
not a claim that the resulting presentation has minimum dimension. Every
source combination and generator map is updated, so the four identities
in the next section replay the complete result directly.

## 4. Four identities certify the complete presentation

Use row-vector coordinates. Let R contain the old relation rows and S
the new ones. Let W send every old generator to its word in the new
generators. Every retained new generator is one old base or root
generator; let V be that inclusion map.

The certificate supplies integer matrices A, B and D with

```math
RW=AS,\qquad SV=BR,\qquad
I-WV=DR,\qquad VW=I.                             \tag{8}
```

Dimensions are determined by the old and new generator and row counts.
The first two equations say that W and V descend to the presented
quotients. The last two show that the descended maps are inverse.
Consequently P and P' are isomorphic over the integers, not merely over
the rationals. The direct replay evaluates (8) by coefficient addition
without traversing the graph or repeating the elimination choices.

Give each new generator the group value of the old generator selected
by V. The third equation and genuineness of R show that evaluation of
every old generator agrees with evaluation of its W word. The
presentation isomorphism therefore commutes with evaluation in G.

If f final root columns are zero in every row, then

```math
P'\simeq P_{\rm kept}\oplus\mathbb Z^f.             \tag{9}
```

Torsion elements have zero coordinates in the free summand, so deleting
those columns preserves the evaluated torsion image. Their values in G
may themselves have finite order; this does not make the free summand
part of the presentation's torsion.

The resulting rows can still be dependent. If the existing extractor
requires independent rows, retain a maximal rationally independent
subset afterward. The equal-saturation argument in
[oriented refinement](oriented-coprime-basis.md#5-torsion-image-preservation)
preserves the evaluated torsion image in that later step. The graph
reduction itself preserves the full integer presentation before (9).

## 5. Bit cost and the rough-norm interface

The spanning forest has at most v edges. Every constant b_v is a signed
sum of at most v base words, so its entries gain at most `O(log(v+1))`
bits. Anchor substitution has the same bound up to an absolute factor.
The expanded path certificates have at most v entries per vertex.

Let H bound the written coefficient bit lengths in all input rows.
Transporting an additional row involves at most v products of an input
coefficient with a tree constant. Its resulting entries have
`O(H+log(v+1))` bits. Subtracting double-root rows has the same polynomial
bound: their root columns are disjoint, so no repeated elimination chain
multiplies this bound across components. All stored words, source
combinations and direct certificate checks have polynomial bit cost in
the written input. No near-linear implementation bound is asserted.

The implementation accepts sparse integer-column dictionaries as well
as dense rows. This avoids allocating an input matrix full of zeros.
Its expanded path and relation certificates can still have quadratic
written size; that cost remains charged. The generator list, including
its v residual generators, is part of the input convention.

The additional unit pass also has polynomial bit cost. Its current
entries are Schur-complement minors of the initial reduced matrix.
The eliminated pivot block has determinant `+/-1`, since every pivot is
a unit. Hadamard's bound therefore controls the entry bit lengths by a
polynomial in the matrix dimension and the initial coefficient lengths.
Applying the same argument to the augmented matrix with an identity
block bounds the retained source combinations. There are at most v
residual-column eliminations.

The graph contract arises without factoring residual norms. Suppose
each supplied residual R is Y-rough, where `Y>=2`, meaning every prime
factor exceeds Y, and `R<=Y^3`. Pairwise coprime atom refinement gives

```math
R=\prod_a d_a^{|w_a|},\qquad d_a>Y.               \tag{10}
```

Three atom occurrences would make `R>Y^3`. Therefore
`sum_a |w_a|<=2`, including repeated powers. This proves the graph
interface from roughness and height alone. The algorithm only reads the
actual signed atom words; it needs neither atom primality tests nor
norm factorization. Witnesses outside that height range can be retained
as the additional rows of Section 3.

The costs of generating witnesses, refining their norm atoms, and
performing the final group operations remain separate. The theorem
does not bound any of those search costs by this matrix-reduction cost.

## 6. Exact controls and finite factoring examples

The [checker](../tools/check_signed_residual_graph.py) verifies:

- All 1728 three-vertex signed graph systems with optional unit anchors
  and double loops, with 3456 independent character counts.
- 420 mixed systems with arbitrary base coefficients and additional
  larger residual words, with 140 further character counts, 420 exact
  dense/sparse certificate comparisons, and 137 later unit pivots.
- Both directions of every presentation map, including JSON round trips.
- Nine invalid contracts and four corrupt certificates.

The four-cycle in [oriented refinement](oriented-coprime-basis.md) reduces
four rows on seven generators to the single base relation
`(-2,22,64)`. Its half extracts 1009 from 1022117.

A separate control shows why (7) matters. At discriminant -84, the
unramified ideals of norms 5, 17 and 41 used by the checker all have
reduced class `(5,4,5)`, of order two. The three relations

```math
y_1+y_2=0,\qquad y_2+y_3=0,\qquad y_1+y_3=0         \tag{11}
```

form an unbalanced triangle. There is no base generator or nonzero
base-only relation. The reduction retains `2y_1=0`; its image extracts
3 from 21. Keeping only balanced cycle relations would lose this
particular supplied torsion image.

For the [retained ordinary 40-bit search prefix](../evidence/oriented-cofactor-prefixes.json), `N=1000036000099`,
`Y=31`, and 320 norm queries produce 249 oriented atom generators in
addition to three base generators. Of its 320 rows, 317 satisfy the
graph contract and three are additional rows. The forest eliminates
105 tree definitions and 139 anchored roots. Three surviving roots are
unused free columns. The complete useful output is

```math
\begin{pmatrix}
0&0&29&-1&-2\\
-83&-237&276&0&1\\
-78&52&-11&0&0\\
-37&23&52&0&0\\
-51&169&-53&0&0
\end{pmatrix}
```

The five remaining generators are the three base generators and two
original atom generators. The extractor verifies all five relation
identities and finds `(500009,15,500009)`, yielding the factor 1000003.
The additional unit pass removes those two residual columns and leaves
the three base rows `(-78,52,-11)`, `(-37,23,52)`, and `(-51,169,-53)`.
It again returns the same ambiguous form and factor. Both full sets of
integer presentation maps are retained in the receipt.
The same graph rows alone, with the three additional rows discarded,
would leave only two base relations and would not extract a factor in
this control. Retaining the larger words is therefore material.

The [implementation](../tools/signed_residual_graph.py) retains every
coefficient certificate. Run the complete retained check with

```sh
python3 tools/check_signed_residual_graph.py \
  --prefix evidence/oriented-cofactor-prefixes.json
```

The [saved receipt](../evidence/signed-residual-graph.json) records the
exact input digest and source hashes. The optional prefix argument also
accepts a standalone complete prefix record. These are finite
results, not a claim about the probability or uniform cost of obtaining
such prefixes.

This implementation and analysis used substantial OpenAI assistant support.
