# Height-homogeneous retractions with mixed occupation numbers

## Result and scope

The coherent reflection construction does not require a physical label to
have one particle number. A label can mix occupation words of the same
height, provided both physical map spaces remain closed under bit deletion.
The complete two-cut evaluator now supports this larger class. Its safe
block split uses parity.

This extends the 32-letter height-16 boundary search from 121 binary
occupation-preserving maps to all 496 binary height-preserving rank-four
boundary retractions, up to changes of basis within the boundary image.
No improved factoring exponent follows from the extension alone.

## 1. Why the complete reduction still applies

The [deletion-closed space theorem](coherent-hot-compound-retractions.md)
requires exact `LP=I`, a common height for each paired embedding column
and retraction row, and deletion closure of both map spaces. The regular
within-cut Gaussian factors descend because they are polynomials in
deletions. Their induced matrices are triangular with identity diagonal
because nonidentity deletions strictly lower height.

Particle number is not used in this argument. It was used only for the
finer block decomposition of the cross-compound matrix. Since each mode
height `2a+1` is odd, every occupation word satisfies

\[
|A|\equiv h(A)\pmod2.
\]

Thus a homogeneous physical label has a well-defined parity even when its
particle number varies. Both complete compressed cuts split into the two
parities of the sum of the physical label heights. The evaluator includes
both sectors.

For each sparse row/column term, the cross-compound entry is zero unless
the two selected occupation masks have equal size. On every surviving
term, the exact normalization

\[
\det B[R,S]
 =t^{(h(R)+h(S))/2}\det K[R,S]
\]

still applies. Its exponent is an integer and is common to every term of
one compressed entry. The complete coefficient degree bound remains
`h(R)+h(S)`. The native constructor sums all surviving terms and checks
every polynomial endpoint entry. It resolves the full generic rank or
retries at a larger series window.

Implementation:

- [Complete mixed-occupation interface](../tools/coherent_hot_mixed.py)
- [Native determinant and parity-block constructor](../tools/coherent_hot_mixed.cpp)

## 2. A literal five-mode validation

On five hot modes, the words `15={0,1,2,3}` and `24={3,4}` both have
height 16. They have four and two particles respectively. Retain the
fifteen coordinate states with masks `0,...,14`, the single mode with
mask 16, and the coherent column

\[
p=e_{15}+e_{24}.
\]

Use coordinate evaluation at word 24 as its paired retraction row. This
gives a 17-letter map with `LP=I`. Every derivative of either term in `p`
lies in the retained coordinate span, so both map spaces satisfy the
theorem. Transposing the embedding and retraction gives a second valid
orientation.

The [independent checker](../tools/check_coherent_hot_mixed.py) restricts
the original interval-16 alternating kernel to these five hot modes on
each physical leg. It then evaluates literal Pfaffians, compresses the
full tensors, restores infinity by polynomial reversal and the prescribed
physical heights, and computes both Smith profiles by the separate Python
polynomial elimination. These literal profiles match the parity evaluator
in both map orientations.

The [retained validation](../evidence/coherent-hot-mixed-controls.json)
records:

| Control | Count |
|---|---:|
| Selected literal Pfaffians evaluated | 104,329 |
| Complete mixed-occupation literal cut comparisons | 4 |
| Mixed-occupation Smith orders compared | 1,156 |
| Literal endpoint entries | 167,042 |
| Literal infinity-restoration entries | 167,042 |
| Previous 32-letter complete cut comparisons | 4 |

This exercises the new mixed-particle case directly. It does not infer
correctness only from occupation-preserving inputs.

```sh
python3 tools/check_coherent_hot_mixed.py \
  --output evidence/coherent-hot-mixed-controls.json
```

## 3. The complete binary height-16 boundary class

There are 28 coordinate words of height below 16. The height-16 cell has
the ordered basis

\[
(e_{24},e_{36},e_{66},e_{129},e_{15}).
\]

Choose a rank-four idempotent on this five-dimensional boundary and keep
all 28 lower states. Every such idempotent over `F2` has the unique form

\[
\Pi=I-w\ell,\qquad w\ne0,\quad\ell(w)=1.
\]

There are 31 choices for `w` and 16 choices for `ell`, hence 496 maps.
To produce an explicit embedding and retraction, choose a nonzero
coordinate `p` of `ell`. For every `i != p`, use

\[
P_i=e_i+\ell_i e_p,\qquad
L_i=e_i^{\mathsf T}+w_i\ell.
\]

These satisfy `LP=I_4` and `PL=Pi`. All boundary derivatives have lower
height and land in the retained coordinate part. No additional deletion
condition needs to be guessed.

Changing the basis of the four-dimensional image changes the cut matrices
by constant invertible row and column maps. Since all four boundary labels
have the same height, their pressure and physical height multiset are
unchanged. Thus the enumeration covers every binary retraction in this
declared class up to this irrelevant basis choice.

The 121 occupation-preserving cases consist of 120 rank-three idempotents
on the four two-particle words while retaining word 15, and the single
choice retaining all four two-particle words while omitting word 15.
The remaining 375 cases are new to the search. The
[driver](../tools/search_coherent_mixed_boundary.py) reuses the completed
exact profiles of the former cases and computes both complete parity
profiles of every new case. Its artifact records whether the full list
finished and retains any unresolved input without treating it as a gain.

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/search_coherent_mixed_boundary.py \
  --seconds 7200 --output evidence/coherent-q32-mixed-boundary.json
```

Numerical pressure roots guide this finite search. A result below the
factoring threshold would still need a strict exact pressure witness and
the complete effective crossing transfer.
