# Regrouping the three-axis cube as a two-axis crossing

## Result and scope

The degree-six homogeneous Gaussian cube family already has a certified
three-cut pressure

\[
\rho<\frac{312359221}{500000000}=0.624718442
\]

from [the optimized rational type](cube-optimized-rational-type.md), and
the coarse bound `rho<5/8` from
[the explicit construction](cube-homogeneous-gaussian.md).  The
three-axis transfer gives the exponent `log_3(1/rho)`, which is below
`1/2`.

This note observes that the **same** eight-leg family also satisfies a
two-axis grid interface.  Pair the eight legs into four super-legs of
length `lambda'=2 lambda`, using any two of the three cube coordinates.
For every piece assignment, use the third coordinate as the internal
super-leg bit and map the smallest of the piece's three cut widths to
the parent row cut.  The resulting two-axis transfer has pressure
`rho'=rho` and exponent

\[
\alpha=\log_2(1/\rho)>\log_2(8/5)=0.6780719\ldots>\frac12.
\tag{1}
\]

The improvement is not a new tensor family and not a smaller constant in
an intermediate profile.  It is a change of the grid growth base from
three to two, made possible by a regrouping that charges only one of the
three available cut factorizations of each piece.

The proof is ordinary and complete at the level of the interface.  The
Lean development and the established `1/2` theorem are unchanged.  The
new claim has not been independently reviewed; the scope is stated
precisely below.  The two-axis profile floor `3/4` of
[the multiaxis profile note](multiaxis-gaussian-profile-floor.md) does
**not** apply: each super-leg here carries two Gaussian legs, so this is
a weighted physical mode outside that floor's hypotheses.

## 1. The relaxed two-axis transfer

We use the exact two-axis grid of
[the multiaxis transfer](multiaxis-cube-transfer.md), with `d=2`.
Recall its geometry.  The parent is a clique on `[g]^2`.  A generic edge
determines a square with four vertices.  The four vertices are grouped
into four labels `X,Y,U,V`; the endpoint is

\[
E_4(x,y,u,v)=[x=v]\,[y=u].
\tag{2}
\]

There are `B=binom(g,2)^2` squares.  The two coordinate cuts are
`{X,Y}|{U,V}` and `{X,U}|{Y,V}`.  Merging the `g` vertices in each plane
along one coordinate produces a child clique on `g` vertices.  The exact
edge and volume estimates are those of the cited note; only the width
hypothesis changes.

The cited two-axis theorem assumes, for every piece `i`, factorizations
at **both** cuts with

\[
r_i^1+r_i^2\le 4\rho'\lambda'
\tag{3}
\]

and separate caps `r_i^a<=2 lambda'`.  The proof uses (3) only through
the total over all cubes and axes; the separate cap is used for the
maximum child-edge weight.  We therefore state the weaker hypothesis
needed here.

> **Lemma (one-cut regrouping).**  Fix `lambda'=2 lambda`.  Suppose the
> four-label endpoint (2) has an exact piece decomposition
> `E_4=sum_i T_i`, and each **same** eight-leg tensor `T_i` has exact
> factorizations at the three coordinate cuts of an eight-vertex cube,
> with widths satisfying
>
> \[
> 0\le r_i^a\le 2\lambda',\qquad
> r_i^1+r_i^2+r_i^3\le 6\rho\lambda'.
> \tag{4}
> \]
>
> Then for every piece assignment to the `B` squares, one can choose for
> each cube one of its three cuts, and map that cut to the parent row
> cut, so that the total row-cut width is at most `2 rho lambda' B`.
> The two-axis grid step then produces a child clique of size `g` with
> the volume contraction
>
> \[
> \frac UV\le \rho+\frac1g+O\!\left(\frac{\log(\lambda'+2)}{\lambda'}\right),
> \tag{5}
> \]
>
> and hence the uniform exponent `alpha=log_2(1/rho)`.

**Proof.**  For a cube whose assigned piece is `T_i`, choose a cut `a`
minimizing `r_i^a`.  Since the three widths are nonnegative,

\[
\min_a r_i^a\le \frac{r_i^1+r_i^2+r_i^3}{3}
\le 2\rho\lambda'.
\tag{6}
\]

Summing over the `B` cubes gives the total row width `2 rho lambda' B`.
The separate cap in (4) is exactly the two-axis cap `2^{d-1}lambda'`
with `d=2`.  The endpoint of the eight-leg cube is the product of its
four antipodal equalities.  Under any two-coordinate pairing it becomes
(2), with the antipodal quadrants read in the opposite internal order.
The eight legs are the same variables; changing the pairing only changes
which two legs form each of the four labels.  Therefore the local
identity `E_8=sum_i T_i` holds in every regrouping, and the piece `T_i`
is not replaced by a different tensor.  A coordinate permutation of the
cube moves the chosen cut to the parent row coordinate and preserves the
endpoint; it only relabels the eight legs.  Thus the local tensor is a
valid two-axis piece with the required row factorization.  The child-edge
estimate uses only the row cap and the total row width, so the proof of
the cited two-axis theorem applies verbatim.  This proves (5).  The
recursion schedule of the cited note with `d=2` gives
`alpha=log_2(1/rho)`.  `square`

Two points about the hypothesis are important.

- The lemma does **not** assert a per-piece bound
  `r_i^a+r_i^b<=4 rho lambda'` for a fixed pair of cuts.  Such a bound
  can fail: the three-axis family only controls the sum of all three
  widths.  The choice is made per cube, and only the chosen cut is
  charged.
- The other parent cut is not used in the contraction.  The two-axis
  grid endpoint contains both diagonal equalities, but the row
  contraction needs only the row factorization.  This is why a one-cut
  factorization per cube is sufficient.  The cited two-axis theorem
  assumed both cuts so that it could choose either orientation after
  seeing the assignment; here the parent row orientation is fixed and
  the cube coordinate realizing it is chosen locally.

## 2. Application to the Gaussian cube

Take the existing eight-leg family of
[cube-homogeneous-gaussian.md](cube-homogeneous-gaussian.md).  It has
the exact endpoint at (1,1,1), a polynomial number of pieces, exact
factorizations at all three coordinate cuts, separate widths at most
`4 lambda`, and total three-cut width

\[
r_i^x+r_i^y+r_i^z
 \le 12\rho\lambda+O(\log(\lambda+2)).
\tag{7}
\]

Put `lambda'=2 lambda`.  Then (7) is

\[
r_i^x+r_i^y+r_i^z
 \le 6\rho\lambda'+O(\log(\lambda'+2)),
\tag{8}
\]

and each separate width is at most `2 lambda'`.  The lemma applies with
the same `rho`.  The pairing of the eight legs into four super-legs is
explicit: fix two cube coordinates, pair the two values of the remaining
coordinate inside each quadrant.  The four quadrants are `X,Y,U,V`; the
antipodal label pairs are the two diagonals, so the endpoint is (2).
The two antipodal quadrants `U,V` are read in the opposite internal
order, because the antipode reverses the remaining coordinate.  This is
a fixed basis permutation inside the column labels; it does not change
any cut rank.  The three possible coordinate pairings are the three
choices of the parent axes.

The cited construction's pieces and factor entries are uniformly
computable in `2^{O(lambda)}=2^{O(lambda')}` bit operations.  The piece
count is polynomial in `lambda'`.  The top constructor, equitable
grouping, table queries, witness recovery, and terminal contraction are
those already charged in the two-axis transfer.  No new oracle or
nonconstructive existence statement is introduced.

Using the coarse certified bound `rho<5/8` gives

\[
\alpha>\log_2(8/5)=0.6780719\ldots.
\tag{9}
\]

Using the sharper exact certificate
`rho<312359221/500000000` gives

\[
\alpha>\log_2\!\left(\frac{500000000}{312359221}\right)
 =0.6786\ldots.
\tag{10}
\]

Both are strictly above `1/2`.

## 3. What is and is not claimed

- The complete deterministic factoring time is
  `exp(O(n/(log n)^alpha))` with `alpha` as in (9) or (10), with all
  costs of the existing two-axis evaluator included.  The derivation is
  an ordinary proof; it is not Lean-formalized.
- The result is a **stronger complete bound** only in the sense of the
  new regrouping lemma.  The local Gaussian family, its exact
  factorizations, and its rational pressure certificate are inherited
  from the cited notes.  A failure of those inherited inputs would
  remove the improvement.
- The per-piece two-cut pressure of the family can be as large as
  `6 rho lambda'` and is not below `4 rho' lambda'` for a fixed pair.
  The improvement uses the assignment-dependent choice of one cut per
  cube; it is not a claim that the family satisfies the stronger
  two-cut interface of the earlier note.
- The two-axis profile floor `3/4` is for a Gaussian kernel with one
  binary mode per physical leg.  Here each physical super-leg is a pair
  of Gaussian legs.  The floor's weighted-mode exclusion applies.
- The existing three-axis exponent `log_3(1/rho)` and the established
  `1/2` theorem are unchanged.

## 4. Checks

The [checker](../tools/check_two_axis_cube_regrouping.py) verifies

1. for `lambda=1`, all three two-coordinate pairings turn the eight-leg
   antipodal endpoint into (2), and the two parent cuts are the
   corresponding coordinate cuts;
2. the arithmetic `min(r_1,r_2,r_3)<=(r_1+r_2+r_3)/3`, the total
   `2 rho lambda' B`, the separate cap `2 lambda'`, and the two-axis
   edge and volume formulas for finite `g`;
3. the exact comparisons `2(5/8)^2<1` and
   `2(312359221/500000000)^2<1`, so both exponents exceed `1/2`;
4. the existing optimized-type certificate and the three exact
   one-particle profiles `(2,4,8,10)` are still replayed.

The checker does not re-prove the Gaussian tensor factorization or the
two-axis evaluator; those are inherited from the cited notes and their
checkers.  It also does not claim independent review of the new
one-cut regrouping lemma.
