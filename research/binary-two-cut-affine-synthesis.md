# Two-cut affine synthesis for complete binary crossing tensors

## Result and scope

This construction assigns an independent polynomial to every entry of a
binary crossing tensor. It fixes the literal permutation endpoint, then
solves simultaneous row-lattice conditions at zero and infinity by exact
finite-field elimination. Constant row flags over either F2 or F4 can be
included independently at the two cuts. All source coefficients and the
full polynomial degree remain explicit.

The construction does not require occupation preservation or a Gaussian
factorization of the unknown tensor. The height ray and the two reference
polynomial frames are fixed. Row-lattice conditions are sufficient Smith
prefix bounds; they do not describe every matrix with those Smith bounds.
The reported finite searches supply no improved factoring exponent.

## 1. Every entry is free subject to the strict endpoint

Let the four physical labels be binary, in the order `(x,y,u,v)`, and let

\[
E_{xyuv}=\mathbf1_{x=v}\mathbf1_{y=u}.
\tag{1}
\]

For definiteness use the explicit polynomial tensor B whose sixteen
binary coefficient bitsets, in index order `x+2*y+4*u+8*v`, are

```text
(14, 30, 30, 60, 15, 15, 31, 30, 15, 31, 15, 30, 15, 15, 15, 14).
```

Bit j is the coefficient of the parameter raised to power j. Evaluating
each bitset at one verifies B(1)=E. The maximum degree is five.
For any nonnegative integer clearing exponent L and D at least 5+L,
every degree-at-most-D tensor with endpoint E is uniquely expressible as

\[
T(w)=w^L B(w)+(1+w)Q(w),\qquad \deg Q_{xyuv}<D.
\tag{2}
\]

This follows entry by entry from divisibility by w+1 in characteristic
two. There are exactly 16D free field coefficients. Equation (2) is an
affine parameterization of the whole stated endpoint fiber, not a small
perturbation neighborhood.

The parameter L makes Laurent searches possible without a new endpoint
argument. The tensor w^{-L}T has support from -L through D-L. Multiplying
it by w^L clears its negative powers once; the zero and infinity profiles
shift oppositely, so their product cost is unchanged. The physical seed
is the cleared polynomial T, and its complete degree D is retained.

## 2. Exact affine equations at both local places

The infinity flattening is `(x,u)|(y,v)`, and the zero flattening is
`(x,y)|(u,v)`. The recorded fixed height ray is

\[
h_\infty=(0,0,-1,-1),\qquad
h_0=(-1,-1,0,0).
\tag{3}
\]

To build complete polynomial matrices, the zero cut substitutes w=z and
multiplies by z^2; the infinity cut substitutes w=z^{-1} and multiplies
by z^{D+2}. Each entry also receives its four-label height from (3).
Both matrices are polynomial in z, including every coefficient of Q.
These common scalar shifts are subtracted again in the reported profiles.

The reference tensor w^L B supplies exact polynomial matrices U and V
with determinants that are units at the local place. They diagonalize
the reference matrix by the
[regular fraction-free elimination](weyl-lattice-filter-synthesis.md).
Its row and column operations multiply only by polynomials with nonzero
constant term, so the matrices and their inverses are regular.

For an invertible constant row basis W and an integer target vector t,
require row i of W U A(T) V to be divisible by z^{t_i}. Each coefficient
below that row target is an affine linear equation in the same 16D
unknowns from (2). Combine the equations from both places and solve them
over F4. When W has binary entries, one binary inhomogeneous solve and
one binary homogeneous solve give the two F4 coefficient planes.

Every reconstructed witness is checked against all the affine equations.
Its complete source polynomials are then passed to the independent Smith
profiler. There is no truncated-series feasibility decision.

Let t_{(1)},...,t_{(4)} be the row targets in increasing order. Every
k-by-k minor is divisible by at least the sum of the k smallest targets.
Regular multiplication by W, U, and V preserves the determinantal ideals.
Consequently the actual Smith orders satisfy

\[
\sum_{i=1}^k \sigma_i\ge\sum_{i=1}^k t_{(i)},
\qquad 1\le k\le4,
\tag{4}
\]

before subtracting the explicitly recorded common scalar shift. Full
profiles, both paired physical partitions, and subsequent numerical
scale optimization are kept separate from these exact feasibility claims.

## 3. Complete constant flags and finite target families

When targets are in increasing order, a row basis matters only through
its chain of tail spaces. Adding a later row to an earlier row preserves
the row-divisibility conditions. Thus complete flags represent every
constant row-basis choice. Their counts in dimension four are

\[
15\cdot7\cdot3=315\quad\text{over F2},\qquad
85\cdot21\cdot5=8925\quad\text{over F4}.
\tag{5}
\]

The binary enumerator deduplicates all invertible binary frames by their
tail spaces. The F4 enumerator chooses a projective vector in each
successive quotient. An independent check reconstructs all three tail
spaces of each F4 frame and verifies that the 8,925 chains are distinct.

At D=5 and L=0, the unshifted reference profiles are

\[
\sigma_\infty=(-5,-2,0,3),\qquad
\sigma_0=(-1,2,4,7).
\tag{6}
\]

The [complete binary-flag record](../evidence/binary-two-cut-complete-f2-flags-degree5.json)
tests all 315-by-315 flag pairs for the reference target and each of its
eight single-row increases by one. Only one flag pair is feasible at the
reference target; none of the eight increased targets is feasible. This
is exact infeasibility for those row-lattice targets, not an obstruction
to every crossing tensor or every height ray.

The search driver over F4 also tests balanced transfers of one, two, or
three units between two target rows, retaining the sorted targets whose
sampled critical cost is no larger than the reference target's. There
are 51 targets including the reference. The floating selection of that
finite target list is explicit; feasibility of each selected target is
decided exactly.

The [complete Laurent F4 flag record](../evidence/binary-two-cut-laurent-complete-f4-flags.json)
uses D=12,L=4 and all 8,925 flags independently at each cut. It completed
all 51 selected targets. The reference target has one feasible flag pair;
all fifty nonreference targets are infeasible. The reference solution is
one affine line over F4, and all four points have the profiles (6) after
removing the common clearing. This last statement is independently
checked by literal minors, not inferred from sampled optimizer values.

Before considering a flag, the solver imposes the necessary minimum-row
condition at both cuts in the original regular frames. Every invertible
constant row map preserves that common entry floor. These constraints
therefore reject whole groups of flag pairs without losing a feasible
solution. The logical coverage of 79,655,625 pairs per target counts
these grouped rejections; it is much larger than the number of individual
pair intersections or Smith calculations. Independent controls compare
pruned and unpruned affine systems and their solution dimensions.
Only 748 candidate pairs survived the single-cut constraints for direct
affine-system intersection across the complete 51-target run.

The [complete Laurent binary-flag run](../evidence/binary-two-cut-laurent-complete-f2-flags.json)
also checks the reference target and every single-row increase by one.
As in the polynomial run, all eight increased targets are infeasible.
The smaller exploratory records at
[129 F4 flags without Laurent clearing](../evidence/binary-two-cut-f4-flags-degree12.json)
and [129 flags with clearing](../evidence/binary-two-cut-laurent-f4-flags.json)
are retained separately. Both completed their stated 51 targets; the
full Laurent flag record supersedes their limited flag coverage for that
target family.

## 4. A finite Laurent window suffices for those targets

Every target in the preceding 51-target family has an infinity minimum
at least -8 and a zero minimum at least -4 after the common scalars are
removed. This alone gives an exact support bound, independent of the
chosen constant flags or the initially proposed degree.

Indeed, regular invertible local changes cannot turn a matrix whose
entries all have order at least a into one with an entry of order below a.
Write \(\widetilde T=w^{-L}T\) for the uncleared Laurent tensor
from Section 1. Applying the inverse changes therefore shows, entry by entry,

\[
-\deg \widetilde T_{xyuv}-u-v\ge-8,
\qquad
\operatorname{ord}\widetilde T_{xyuv}-x-y\ge-4.
\tag{7}
\]

Hence every nonzero Laurent entry is supported inside [-4,8]. Clearing
by w^4 produces a polynomial of degree at most twelve. Thus D=12,L=4
already contains **every Laurent tensor meeting these specified targets
in the fixed local frames**. Arbitrarily high coefficient degrees cannot
rescue an infeasible target in this list. Other targets, different frames,
and different physical heights are outside this conclusion.

## 5. Independent controls and reproduction

The [independent checker](../tools/check_binary_two_cut_lattices.py)
expands every minor of both four-by-four weighted matrices over F4.
It compares these literal determinantal profiles with the native engine,
checks every endpoint coefficient, reconstructs colored witnesses under
nonbinary flags, and tests common Laurent clearings of one, four, and
seven. It also replays all four F4 points on each one-dimensional
reference affine solution space at degrees five, eight, and twelve.
Its [saved evidence](../evidence/binary-two-cut-lattices.json) records the
exact counts and source hashes.

```sh
python3 tools/check_binary_two_cut_lattices.py --output /tmp/binary-lattice-replay.json
python3 tools/search_binary_flag_lattices.py --degree 5 --output /tmp/binary-flags.json
python3 tools/search_binary_f4_flags.py --degree 12 --baseline-shift 4 --flags 8925 --seconds 3600 --output /tmp/laurent-flags.json
```

The implementation consists of the
[affine lattice solver](../tools/search_binary_two_cut_lattices.py),
[binary flag enumerator](../tools/search_binary_flag_lattices.py), and
[F4 flag enumerator](../tools/search_binary_f4_flags.py).
The separate [affine-pencil](../tools/search_binary_affine_crossings.py)
and [polynomial-neighborhood](../tools/search_binary_polynomial_neighborhood.py)
drivers supply reusable endpoint constructors and pressure helpers.
Their presence does not assert completion of any unrecorded search run.

## 6. A complete finite schedule of simultaneous row-target changes

A further search allows all eight row targets to change together. Starting
from the reference vector, take every integer delta vector in
`[-2,2]^8` with sum zero and sorted targets at each cut, and include the
original single-row transfers of one, two, or three units. This gives
23,099 distinct raw sorted targets. A numerical critical-pressure
comparison at scale `0.29277969091981376`, with tolerance `1e-12`,
selects 7,079 targets. This floating filter defines a finite schedule;
it is not an all-scale sign certificate or target-space optimality claim.

The [complete record](../evidence/binary-two-cut-zero-sum-target-box.json)
uses degree twelve, common clearing four, and **129 selected constant F4
flags per cut**, consisting of 128 deterministic pseudorandom flags and
the identity. This is a larger target family with fewer flags than the
separate complete-flag experiment above. Every selected target and flag
pair was decided. Necessary common constraints reject 1,430 targets
independently of the flags. Across the selected flags, only five pairs
survive the separate-cut constraints for direct affine intersection.

There are exactly five feasible targets. They add a common integer
`g` to all four infinity targets and subtract it from all four zero
targets, for `g=-2,-1,0,1,2`. Each feasible system has affine dimension
one over F4. The [independent checker](../tools/check_binary_target_box.py)
reconstructs all four points of each line, checks literal endpoints and
every minor of both cuts, and multiplies the coefficients by the matching
monomial. All five normalized four-point families are exactly the
original degree-five reference line after its common clearing. Thus the
survivors are **five monomial shifts of the same four-point affine line**,
with the same complete cost, approximately `rho=0.804022810651`.

Every target still has the minimum bounds used in Section 4. Hence the
same finite support argument covers all Laurent coefficients for these
specified row targets and selected flags. It does not extend the 129-flag
coverage to all 8,925 flags.

The [replay record](../evidence/binary-target-box-replay.json) independently
reconstructs the 23,099/7,079 target sets, checks the flag manifest, and
re-eliminates the 17 distinct common-minimum systems underlying all 1,430
common contradictions. It verifies 50 literal weighted cut profiles and
400 endpoint coefficients across all fiber points and saved representatives.
The remaining flag-pair exclusions are checked against the completed
manifest, without independently rerunning the whole search.

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/search_binary_f4_flags.py --degree 12 --baseline-shift 4 --target-box 2 --flags 128 --seconds 750 --output /tmp/binary-target-box-search.json
OPENBLAS_NUM_THREADS=1 python3 tools/check_binary_target_box.py --output /tmp/binary-target-box-replay.json
```

A separate attempt to enlarge single-row transfers through six produced
no additional targets after its stated pressure filter. Its redundant
higher-degree replay was stopped at 22 of the same 51 targets, already
covered completely by the degree-twelve Laurent result. That interrupted
prefix is not presented as another completed search.
