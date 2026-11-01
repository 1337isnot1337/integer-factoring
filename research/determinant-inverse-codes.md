# Physical crossing codes from determinant derivatives

## Mechanism and current status

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338858)
changes the physical code throughout its occupations. It replaces the
protected coordinate core with the derivative space of a determinant.
Same-row and same-column rook pairs disappear; other occupation states
become coherent minors. This differs from another boundary retraction of
the retained 96-letter code.

The physical maps below are acquired exactly. Their deletion closure and
split identity supply the existing coherent crossing construction. Complete
profiles for the prescribed finite panel are a separate diagnostic; no
favorable ranks, pressure or improved factoring exponent follow from the
dimension formula. The complete score remains **1/2**.

The protected-core floor requires literal retention of every low coordinate
state and therefore does not apply to these maps. The equal-block Gaussian
kernel floor does not automatically apply after this non-Gaussian physical
compression. This is not a claim to escape every invariant. No new Lean
formalization or subagents are used.

## 1. All minors form the derivative space, in characteristic two

On k^2 binary occupation modes, indexed `(i,j)` with `0<=i,j<k`, let

~~~
f=det(x_ij),             h_ij=2(ki+j)+1.
~~~

Work over a characteristic-two field, so determinant signs disappear.
Bit deletion D_ij deletes x_ij from a monomial containing it and kills a
monomial not containing it. If a set of deleted variables repeats a row
or column, its derivative of f is zero. Otherwise the deleted variables
form a partial matching. The surviving terms are exactly the determinant
of the complementary submatrix.

Conversely, every minor `p_(R,C)=det(X_R,C)` is such a derivative: delete
any matching between the complementary rows and columns. The minors of
order r have disjoint monomial supports for different pairs `(R,C)`,
because a monomial determines its row set and column set. They are all
nonzero, hence linearly independent. Thus the complete derivative space is

~~~
U=span{p_(R,C): |R|=|C|},
dim U_r=binomial(k,r)^2,       dim U=binomial(2k,k).
~~~

This dimension/minor-space fact is classical; see Shafiei,
[*Apolarity for determinants and permanents of generic matrices*](https://arxiv.org/pdf/1212.0515),
Lemma 1.3 and Section 2.1. That paper explicitly assumes characteristic
zero or characteristic greater than two. The elementary squarefree
deletion/disjoint-support proof above supplies the characteristic-two
statement used here, without importing its apolar-ideal theorem.

Each minor has one occupation r and one height

~~~
h(R,C)=2k*sum R+2*sum C+r.
~~~

The top height is k^3. A single deletion sends a minor to a smaller minor
or zero, so U is deletion closed. These facts hold for arbitrary k; they
are not inferred from the finite panel.

## 2. Two matching-field coordinate retractions

For each R,C in increasing order, choose its diagonal matching, or choose
its antidiagonal matching. Let L_(R,C) read that single occupation word.
Different minors have disjoint supports and coefficient one at their
chosen matching, so `LP=I` for the embedding P with columns p_(R,C).

Deleting any edge of a chosen matching leaves the same diagonal, or
antidiagonal, rule on the remaining row and column sets. The coordinate
row space of L is therefore deletion closed. Both maps preserve height
and occupation, meeting all hypotheses of the
[coherent compound interface](coherent-hot-compound-retractions.md).

## 3. A dense Cauchy dual at every occupation

Choose distinct field elements u_0,...,u_(k-1),v_0,...,v_(k-1), and put
`A_ij=1/(u_i+v_j)`. Every square minor of A is nonzero. For completeness,
if a linear dependence of an r-column submatrix vanished at its r chosen
u_i, the rational function `sum_j c_j/(X+v_j)` would have numerator
degree at most r-1 with r distinct roots. It is zero identically.
Evaluating the numerator at each v_j forces every c_j to vanish. This
proves the needed Cauchy property over the current field.

Set `q_(R,C)=det((A_ij*x_ij)_(R,C))`. These weighted minors are also
deletion closed: deleting x_ij gives A_ij times the smaller weighted
minor. Their ordinary bilinear pairing with p is diagonal, with

~~~
<q_(R,C),p_(R,C)>=det(A_R,C) != 0.
~~~

Consequently the rows
`L_(R,C)=q_(R,C)^T/det(A_R,C)` give `LP=I`. Their row space is deletion
closed and both gradings are preserved. The pairing uses no conjugation.
This is the existing
[bilinear dual-space interface](coherent-inverse-f4-pairings.md) applied
to a different explicitly acquired family.

The finite F8 choice is `F2[z]/(z^3+z+1)`, with u_i encoded by i and
v_j encoded by k+j. It works for k<=4. The all-k construction needs a
field with at least 2k elements and must charge its representation and
arithmetic. No field enlargement is part of the registered finite search.

## 4. A polynomial-size height character, separate from cut profiles

Writing e_r for the elementary symmetric polynomial, the full physical
charge/height inventory is

~~~
Q_k(y,t)=sum_(r=0)^k y^r*t^r
          *e_r(1,t^(2k),...,t^(2k(k-1)))
          *e_r(1,t^2,...,t^(2(k-1))).
~~~

Two subset-count dynamic programs and polynomial convolution compute
this character in polynomially many operations on polynomial-size integer
coefficients, since the height degree is k^3. Complementing both R and C
sends `(r,h)` to `(k-r,k^3-h)`, so the mean height is k^3/2. This inventory
does not determine cut ranks or supply a fast growing-k crossing evaluator.

The [constructor](../tools/determinant_inverse_codes.py) explicitly
enumerates finite sparse maps only for k<=4, paying for all permutation
monomials. Treating a growing k as free would ignore potentially large
basis-support and compound costs. For a fixed finite seed, the existing
uniform typed-word construction can instead charge its entire finite data.

## 5. Acquired controls and the bounded next diagnostic

The [map checker](../tools/check_determinant_inverse_codes.py) constructs
all three retractions for k=2,3. The
[receipt](../evidence/determinant-inverse-code-controls.json) contains all
six map pairs, 612 embedding deletions, 612 retraction deletions, 1,308
split-pairing entries, 26 Cauchy minors checked by independent scalar
elimination, and all 1,056 squarefree derivatives of the two top forms.
The dynamic-program character is checked through k=12 for dimension,
charge and complementary-height identities; these are character controls,
not a larger pressure search.

Run `make verify-determinant-inverse-codes`. The prescribed profile panel
is k=2 controls and k=3 with diagonal, antidiagonal and the fixed Cauchy
dual. Reconstruct both complete natural profiles and check selected
entries independently from the original one-particle matrix. A strict
numerical pressure witness needs an exact rational type/log certificate
and the complete uniform transfer before any score claim. Continue beyond
this panel only for such a witness or a new analytic uniform rank
mechanism. Above-critical finite fits alone do not justify a k, field or
retraction walk.

The [profile driver](../tools/probe_determinant_inverse_codes.py) requires
a clean frozen checkout and an external `--output` path. Precision starts
at cutoff 127 and increases only if an explicit tail remains. The endpoint
and entry degree bound give the finite determinant-order cap
`2*sum_(x,y in charge block)(h_x+h_y)`; no individual missing order is
inferred from this aggregate bound. The numerical fit uses an 81-point
log-scale grid on `[-9,log 4]`, followed by 64 golden-section steps at
each grid-local minimum. It claims no global scale optimum. The k=2
controls compare every supported charge-block polynomial entry and full
profile against original-matrix minors and binary field restriction.
The k=3 controls compare selected entries and their smaller complete
profiles; these are not a second reconstruction of the full k=3 matrices.

Credit: classical determinant derivative/minor spaces; existing coherent
inverse-system, dual-space and F8 compound constructors. This application,
maps and diagnostic were developed with substantial OpenAI assistant support.
