# Conjugate the physical code before compressing a crossing

## Registered structural question

Base: `07faf95`. Can a local basis change of an error-detecting code
retain its rate and one-site detection while breaking the exact Frobenius
identity that closed the earlier [CSS source](collective-css-pauli-compression.md)?
That identity identifies the full tensor with a two-cell source at `z^2`,
even at arbitrary physical heights. Changing the physical image before
contraction is the proposed mechanism; no improvement is presumed.

Use the same nine-dimensional code in four ternary coordinates over F4.
On its first coordinate apply the explicit orthogonal matrix

```
v=(1,omega,omega^2),  w=(1,1,0),
O=(I+v v^T)(I+w w^T),
P_new=(O tensor I tensor I tensor I) P.
```

Both vectors have zero dot square in characteristic two, so each
transvection is orthogonal. Verify `O^5=I`, `O!=I`, `P_new^T P_new=I`,
and explicitly check that conjugating a Weyl generator leaves the Weyl
list. This supplies a precise changed algebraic assumption. The projector
onto the code still detects every traceless one-site operator: local
conjugation preserves trace and the one-site matrix algebra. Check all
32 nonidentity one-site Weyl operators literally for the new split.

Construct the complete degree-twelve nine-letter crossing. Equivalently,
conjugate the first three-letter cell on its four matched physical legs,
then use the original sparse cosets; this avoids discarding new amplitudes
or treating an auxiliary matrix as a scalar. Preserve every coefficient
and verify the full endpoint and independent literal contraction controls.

The bare nine-letter criterion is already universally excluded by the
[box theorem](bare-seven-to-twelve-letters.md), so no bare optimization is
registered. Continue to weighted physical types only if the new tensor
has a nonzero odd-degree coefficient, which cannot arise by constant logical
basis changes of the former `z^2` family. A repeated exact Frobenius identity
stops this test before fitting. Nonzero odd coefficients alone are not a
rank or pressure improvement.

If that structural gate passes, run exactly three bounded fits with all
64 nonreference physical heights free: a near-zero balanced start, a
two-logical-digit graded balanced start, and a near-zero free-angle start.
Use fixed RNG seed `20261101`, scale `.1`, log-tangent zero initially,
height range `[-128,128]`, log-scale `[-9,log4]`, log-tangent `[-2,2]`,
at most 100 iterations and 20 seconds per fit. Retain query counts,
failures, complete best/last profiles and all starts. A support below
`1-1e-6` requires independent exact weighted orders and outward pressure
before any stronger claim. Otherwise archive the finite result and stop
this code. No rotation, field, code or parameter sweep follows.

The complete factoring score remains `1/2`; no actual improvement, Lean
addition or independent mathematical review is claimed. Substantial OpenAI
assistant support.

## Completed fixed-code diagnostic

The explicit matrix is

```
O = [[omega,0,omega^2],
     [omega,omega,1],
     [1,omega^2,omega^2]].
```

It is orthogonal of order five. Its conjugate of the shift Weyl matrix
is outside all 27 scalar Weyl matrices. The new split passes all 32
one-site nonidentity Weyl annihilations and the full crossing endpoint.
The tensor has nonzero odd-degree coefficients: the numbers of nonzero
entries in degrees `1,3,5,7,9,11` are `36,481,571,542,533,312`.
Thus the stated structural gate passes, and the former exact `z^2`
coefficient identity does not persist for this code.

The [independent checker](../tools/check_nonclifford_css_crossing.py)
compares all 81 transformed cell entries by a literal Kronecker basis
change, then all 6,561 compressed polynomials by scalar packed arithmetic,
examining every one of the 531,441 coset summands. It also recomputes 162
complete weighted orders with untruncated fraction-free arithmetic and
compares both native cut profiles. The
[receipt](../evidence/nonclifford-css-controls.json) retains the full
integer gauges, complete orders, split data and source hashes. These
finite checks passed in 14.67 seconds.

All three registered fits returned optimizer convergence: 12, 12 and 15
queries, respectively, with no rejected query. The complete fits took
0.145 seconds. Their best support values were `1.414385533173`,
`1.414385587608` and `1.414385533202`, well above the required one. The
[complete fit record](../evidence/nonclifford-css-fit.json) retains all
starts, best/last coordinates, both complete profiles, statuses and query
counts. These local outcomes do not establish a global optimum or exclude
other weighted gauges.

The continuation gate fails. The specific code stops here; no neighboring
rotation or larger code is searched. It supplies an explicit departure
from the old Frobenius identity, but no useful pressure improvement and
no stronger complete factoring bound. Run `make verify-nonclifford-css`
for exact construction/profile checks; the numerical fit archive is a
historical finite experiment rather than a deterministic verification target.
