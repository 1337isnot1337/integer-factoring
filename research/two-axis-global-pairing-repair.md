# A global-pairing replacement for the rejected two-axis regrouping

## Status and relation to the audit

The [audit](two-axis-cube-regrouping.md) correctly rejects the original
one-cut proposal in
[the retained proposal note](two-axis-cube-regrouping-proposal-rejected.md).
The rejected step chose the smallest cut independently for each cube.  That
is not one global plane orientation, a hidden-axis cut is not a parent cut,
and a permutation applied to one summand need not preserve the endpoint.

This note supplies a different replacement satisfying the reopening
conditions in Section 5 of the audit.  It uses **one common internal
coordinate for the whole piece assignment**, not independent per-cube
choices, and it asserts only the total parent-axis width used by the
multiaxis transfer, not a uniform per-piece fixed-cut width.  The result is
the same ordinary-proof claim,

\[
\alpha=\log_2(1/\rho)>\log_2(8/5)=0.6780719\ldots>\frac12,
\tag{1}
\]

with `rho<5/8` (or the sharper `0.6247...`), but the proof is the
total-width interface, not the rejected fixed-cut one.

This is an interface proof.  The inherited Gaussian factorization, its
three-cut profiles and rational certificate, and the two-axis evaluator are
not re-proved here.  The full transfer has not received independent
mathematical acceptance, and the project scoreboard remains `1/2` until
that review occurs.

## 1. The global-pairing lemma

Let `E_8` be the eight-leg antipodal endpoint.  The Gaussian family has an
exact decomposition `E_8=sum_i T_i`; each **same** eight-leg piece `T_i`
has exact factorizations at the three coordinate cuts, with widths

\[
0\le r_i^a\le 2\lambda',\qquad
S_i:=r_i^1+r_i^2+r_i^3\le 6\rho\lambda'+O(\log(\lambda'+2)),
\tag{2}
\]

where `lambda'=2 lambda`.  The family and (2) are inherited from
[cube-homogeneous-gaussian.md](cube-homogeneous-gaussian.md); putting
`lambda'=2 lambda` only rewrites `12 rho lambda` as `6 rho lambda'`.

> **Lemma (global internal coordinate).**  Fix one complete piece
> assignment to the `B=binom(g,2)^2` squares.  Choose a coordinate `c`
> maximizing `sum_i r_i^c` over the three coordinates, and let `a,b` be
> the other two.  Then:
>
> 1. a single coordinate permutation `P_c`, applied to **every** piece,
>    gives an exact identity `E_4=sum_i P_c(T_i)` on the fixed four-label
>    parent grid;
> 2. the pieces `P_c(T_i)` have parent factorizations at the cuts `a,b`
>    with the inherited widths `r_i^a,r_i^b`;
> 3. the total parent-pair width satisfies
>    `sum_i(r_i^a+r_i^b)<=4 rho lambda' B+O(B log(lambda'+2))`, and the
>    smaller of the two parent-axis totals is at most half of this;
> 4. the separate cap is `r_i^a<=2 lambda'`.

**Proof.**  The four super-labels are the two values of the internal
coordinate inside each of the four quadrants.  For a fixed internal
coordinate `c`, the map from the eight physical legs to the four
super-labels is a fixed bijection of the eight leg positions; call it
`P_c`.  The antipodal endpoint is invariant under coordinate
permutations, and the two antipodal quadrants are read in the opposite
internal order; the checker verifies `E_8=regrouped_endpoint` for all
`2^8` assignments and all three choices of `c`.  Hence
`P_c(E_8)=E_8`.  Applying the common map to the exact identity
`E_8=sum_i T_i` gives

\[
E_8=P_c(E_8)=\sum_i P_c(T_i).
\tag{3}
\]

The relabeled endpoint is exactly the four-label endpoint
`E_4(x,y,u,v)=[x=v][y=u]`; this is the same computation as the pairing
identity in the audit, now used with one common map rather than a
different map per summand.  The parent factorizations of `P_c(T_i)` at
`a,b` are the inherited coordinate-cut factorizations, with the same
widths.

For the widths, the maximum is at least the average:

\[
\sum_i(r_i^a+r_i^b)
 =\sum_i S_i-\sum_i r_i^c
 \le \frac23\sum_i S_i
 \le 4\rho\lambda'B+O\!\left(B\log(\lambda'+2)\right).
\tag{4}
\]

The smaller of `sum_i r_i^a` and `sum_i r_i^b` is at most half of (4).
The cap is inherited from (2).  `square`

The choice of `c` depends on the complete assignment, exactly as the
axis choice in the multiaxis transfer depends on the assignment.  It is
one global orientation for that assignment, not one choice per cube.

## 2. Why the audit's controls do not block this interface

The audit's four exact controls remain valid.  They rule out the rejected
per-cube construction; they do not rule out the common-map total-width
construction above.

1. **Independent minima.**  The audit's widths `(0,3,4)` and `(3,0,4)`
   give independent-minimum sum `0` and global-minimum axis sum `3`.
   The repair chooses the internal coordinate `c=2` (sum `8`); the parent
   pair is `(0,3)+(3,0)=6`, and the smaller parent total is `3`.  This is
   the global minimum, not the independent-minimum sum.  The checker
   verifies the equality with the audit's `global_best`.

2. **Hidden-axis cut.**  The repair never inserts the hidden cut as a
   factor between parent planes.  The parent cuts `a,b` are the other two
   coordinate cuts; they partition whole super-labels.  The hidden cut is
   charged only through `S_i` and the choice of `c`.  The audit's
   incidence check confirms that axes `0,1` cut whole parent labels while
   axis `2` splits all four super-labels; the repair uses the former as
   parents and the latter only as an internal bit.

3. **Piece-specific permutation.**  The audit's counterexample applies a
   coordinate swap to `T_1` only.  The repair applies the **same** `P_c`
   to every summand.  Linearity gives
   `sum_i P_c(T_i)=P_c(sum_i T_i)=P_c(E_8)=E_8`.  The audit itself states
   the common-permutation case in Section 3; the repair is that case.  A
   checker positive control applies the same map to the audit's
   `T_0=E_8+A`, `T_1=A` pair and verifies the sum for all `2^8`
   assignments, while the piece-specific control still fails.

4. **Fixed-cut rank.**  The audit's rank obstruction assumes every piece
   at the fixed cut has width at most `2 rho lambda'+C' log`.  The repair
   does **not** assert that uniform per-piece bound.  It asserts the
   total selected-axis width `W<=2 rho lambda' B+O(B log)` and the
   separate cap `r_i^a<=2 lambda'`.  Rank subadditivity at the selected
   parent cut requires `2^(2 lambda')<=sum_i 2^(r_i^selected)`, which is
   feasible with one piece at the cap `2 lambda'` and the rest small,
   provided `B>=1/rho`.  The
   audit's Section 4 explicitly leaves open the "single orientation
   chosen from total assignment widths" case; the repair is precisely
   that case.

The finite rank contrast is checked explicitly: at `lambda'=32`, `C=4`,
`rho=5/8`, the audit's uniform bound gives `33^4*2^40<2^64`, while a
total-width assignment with one `2 lambda'=64`-bit piece and the rest zero
has `sum_i 2^(r_i)>=2^64` and total width `64<=2 rho lambda' B`.  The
audit's rank inequality is not violated by the total-width interface; it
only excludes the rejected uniform fixed-cut claim.

## 3. The two-axis transfer

The common map gives the two-axis endpoint and the parent factorizations
required by [the multiaxis transfer](multiaxis-cube-transfer.md) with
`d=2`, under the weaker total-width hypothesis.  Its proof uses:

- the selected parent-axis total `W<=2 rho lambda' B+O(B log)`, not a
  per-piece pair bound;
- the separate cap `2 lambda'` for the maximum child-edge weight;
- the inherited exact pieces, table queries and regrouping costs.

Equations (5)--(6) of the multiaxis note then give

\[
\frac UV\le \rho+\frac1g
 +O\!\left(\frac{\log(\lambda'+2)}{\lambda'}\right),
\tag{5}
\]

and the `d=2` recursion schedule gives `alpha=log_2(1/rho)`.  The
smaller parent-axis total is the one selected for the merge; the other
parent cut remains part of the endpoint but is not the merge orientation.
The common `P_c` fixes the actual parent labels for the whole assignment,
so the child network is consistently wired.

Substituting the inherited certificate `rho<5/8` gives (1); the sharper
`rho<312359221/500000000` gives `alpha>0.6786...`.  Both are above `1/2`.

## 4. Remaining obligations

- The Gaussian piece decomposition, its three-cut profiles and the
  rational pressure certificate are inherited.  A failure there removes
  the improvement.
- The two-axis evaluator, equitable grouping, table queries, witness
  recovery and terminal contraction are inherited and charged in the
  cited transfer.  This note checks the regrouping interface, not a new
  implementation.
- The repair is an ordinary proof and an exact interface checker.  It is
  not Lean-formalized and has not received independent mathematical
  acceptance.  The scoreboard therefore remains `1/2` until that review.

## 5. Checks

The [checker](../tools/check_two_axis_global_pairing.py) verifies:

1. the common-map endpoint identity `E_8=sum_i P_c(T_i)` for all three
   choices of `c` and all `2^8` assignments, including the audit's
   `T_0,T_1` pair under a common map;
2. the audit's piece-specific counterexample still fails (128 mismatches),
   so the positive common-map control is not vacuous;
3. the global width selection: `c` maximizes `sum_i r_i^c`, the parent
   pair is at most `(2/3) sum_i S_i`, the smaller parent total is at most
   half the pair, and the separate cap holds;
4. the audit's `(0,3,4),(3,0,4)` control: the selected parent total is
   `3`, equal to the established global minimum, not the independent-min
   sum `0`;
5. the hidden-axis incidence: the parent cuts partition whole super-labels
   and the hidden cut splits all four;
6. the finite rank contrast: the uniform fixed-cut bound fails while the
   total-width assignment meets the endpoint rank;
7. the inherited one-particle profiles `(2,4,8,10)` and the rational
   pressure certificate, replayed from
   [the audit checker](../tools/check_two_axis_cube_regrouping.py).

The audit's exact controls and rejected proposal remain in place and are
not overwritten.  This note is a replacement interface, not a claim that
the rejected fixed-cut proof was correct.
