# Cell-set minor caps for the gauged two-place relaxation

## Result and scope

This note adds a new family of necessary prefix inequalities to the
[shared-minor profile relaxation](shared-minor-profile-relaxation.md).
For every set of `k` cells `(x,y)`, the cell-set minor with rows in that
set and columns at the swapped cells is nonzero at the endpoint, and its
degree is at most `k d` plus the cell energies of the set. Each normalized
Smith profile therefore satisfies

\[
r_1+\dots+r_k\le kd+\mathrm{Cell}_k,\qquad
s_1+\dots+s_k\le kd+\mathrm{Cell}_k\qquad(1\le k\le q^2),
\]

where `Cell_k` is the sum of the `k` smallest of the `q^2` cell energies
`c_X(x)+c_Y(y)`. These cell-set caps contain the rectangular caps (5) of
the predecessor note and are strictly stronger. Within the class of bounds
that use one minor, an endpoint value and the support bounds alone, they
are the strongest one-cut prefix bounds. The shared lines of the
predecessor note are the only minors common to both cuts.

The hypotheses are those of the predecessor's Section 1: a scalar
Laurent-polynomial strict crossing over any field, the cut `(x,y)|(u,v)`
at zero and `(x,u)|(y,v)` at infinity (or the reverse), and independent
integer heights on all four legs at both places. Rational heights follow
by ramification. Real heights are covered when the profiles are defined
by determinantal prefixes, as in the
[universal binary floor](universal-binary-smith-floor.md), Section 6.
Auxiliary dimension `D>1`, multipoint criteria and type-projected
endpoints are outside the hypotheses.

Consequences for the relaxation, all decided by exact outward dyadic
bounds:

| Point (same paired energies on both sides) | Old caps | With cell caps | Status |
|---|---|---|---|
| q=4, d=1, energies (0,0,2,8/3) | -0.0231521 | **+0.0430257** | predecessor target **excluded** |
| q=5, d=1, energies (0,0,2,8/3,3) | -0.0976380 | -0.0042091 | still negative |
| q=5, d=7/6, energies (0,0,19/9,23/9,37/9) | -0.0889196 | **-0.0190509** | sharper negative target |
| q=6, d=8/7, energies (0,0,15/8,27/10,4,4) | -0.1155333 | **-0.0436118** | sharper negative target |

On the exploration side, labelled floating-point searches found the
following. At `q=4` no negative value was found over 2,429 local
searches. The infimum is `0`, approached but not attained as some
energies diverge (the one-letter limit gives `G=d`), and the smallest
interior local minimum found is about `+0.00509`. At `q=3` every
value found is positive, and values approach zero only as the energies
diverge. At `q=5` and `q=6` the best values found are about `-0.01907` and
`-0.04367`; negativity there is certified by the exact targets above.
The search minima are not theorems. The exact ternary statement proved here
is partial: every point with all paired energies in `[0,1]` has
`G_relaxed>=1/100`, at every `d>=0` (Section 7).

Relaxation points are not tensors. No crossing realizing a negative
target is known. Nothing here is a crossing certificate or changes the
factoring score, which remains `1/2`.

This builds directly on the
[shared-minor profile relaxation](shared-minor-profile-relaxation.md)
(upstream PR #39, branch `research/coupled-frontier-c480f802` at
commit `a3c5bf1`). It uses that note's normalization, convex-minorant
reduction and outward arithmetic, as well as the shared-minor mechanism
of the [universal binary floor](universal-binary-smith-floor.md).
The cell-set inequality and a floating-point probe came from the
coordinating research session (contributed by orda-aska with Claude);
this note proves the inequality and checks it exactly. Two independent
referee passes re-derived Theorem 1 and recomputed every certified
pressure; their corrections to the floating-point descriptions are
included. The work was developed with substantial Anthropic Claude
assistance.

## 1. Setting

We use the predecessor's equations (1)–(3) verbatim. In summary, after
normalization the zero-place matrix has polynomial entries `R_a(w)`,
indexed by `a=(x,y,u,v)`. Its endpoint is `R(1)=[x=v][y=u]`, and

\[
\deg R_{xyuv}\le d+c_x(x)+c_y(y)+c_u(u)+c_v(v),
\qquad c_X=c_x+c_v,\quad c_Y=c_y+c_u .
\tag{1}
\]

The infinity matrix has entries `z^(d+sum_l c_l(a_l)) R_a(1/z)`. The
zero profile `r` is the sorted list of Smith orders at `w=0` of one cut
of `R`. The infinity profile `s` is the sorted list of Smith orders at
`z=0` of the other cut of the infinity matrix. Both lists start at zero.
The paired energies satisfy `c_X,c_Y>=0`, with `c_X(0)=c_Y(0)=0`. The
individual leg energies `c_l` may be negative. We write `X_k` and `Y_k`
for the sums of the `k` smallest paired energies. The relaxation
variable is `gamma=(r+s)/2`, with prefix sums `S_k`.

## 2. The cell-set minors

**Theorem 1.** Let `S` be any set of `k` cells `(x,y)`, and put
`Cell(S)=sum_{(x,y) in S}(c_X(x)+c_Y(y))`.

1. In the cut `(x,y)|(u,v)`, take rows `(x,y)` in `S` and columns
   `(u,v)=(y,x)`. The resulting `k`-square minor `m_S` of `R` equals `1`
   at `w=1`, and its degree is at most `kd+Cell(S)`.
2. In the cut `(x,u)|(y,v)`, take rows `(x,u)` in `S` and columns
   `(y,v)=(u,x)`. The resulting minor `m'_S` has the same two properties.
3. In either cut, the corresponding minor of the infinity matrix equals
   `z^(kd+Cell(S)) m(1/z)`.

Consequently, whichever cut is placed at zero,

\[
\sum_{i\le k}r_i\le kd+\mathrm{Cell}_k,\qquad
\sum_{i\le k}s_i\le kd+\mathrm{Cell}_k,\qquad
S_k\le kd+\mathrm{Cell}_k .
\tag{2}
\]

*Proof.* Order the columns of `m_S` by the bijection
`(x,y) -> (y,x)` from rows to columns. At `w=1` the entry in row `(x,y)`
and column `(y',x')` is `[x=x'][y=y']`, so `m_S(1)` is the determinant of
an identity matrix.

Now consider a determinant term, the product over the rows `(x,y)` of
`R_{x,y,pi(x,y)}` for a bijection `pi` onto the column set. By (1), its
degree is at most the following sum. Each row contributes `d+c_x(x)+c_y(y)`.
Each column `(u,v)=(y',x')` contributes `c_u(y')+c_v(x')`, once over the
column set. The total is

\[
kd+\sum_{(x,y)\in S}\bigl(c_x(x)+c_v(x)+c_y(y)+c_u(y)\bigr)
=kd+\mathrm{Cell}(S).
\]

The same bound holds for every term, independently of `pi`. The column
cut is identical after exchanging the roles of `y` and `u`. Row
`(x,u)` contributes `d+c_x(x)+c_u(u)`, and column `(u,x)` contributes
`c_y(u)+c_v(x)`.

Because each term carries the same total shift `kd+Cell(S)`, the
infinity minor is the reciprocal polynomial with that shift. A nonzero
polynomial has zero order at most its degree. The infinity order is
therefore `kd+Cell(S)-deg m`, which is at most `kd+Cell(S)`.

The sum of the first `k` Smith orders at a place is the minimum
valuation of the `k`-square minors there. The zero cut contains `m_S`,
so its `k`-prefix is at most `ord_0 m_S<=deg m_S<=kd+Cell(S)`. The
infinity cut contains the reciprocal of `m'_S`, so its `k`-prefix is at
most `kd+Cell(S)-deg m'_S<=kd+Cell(S)`. Minimizing over `S` gives the
`k` smallest cells, which proves (2). Averaging the two bounds proves
the bound on `S_k`. QED

**Remarks.**

- *Field and scale.* Only valuations and degrees were used, so (2)
  holds over every field. A nonzero endpoint scale `lambda E` changes
  `m_S(1)` to `lambda^k`. The pressure parameter `t` multiplies `d`, `c`,
  `r` and `s` together, and (2) is positively homogeneous, so it is
  unaffected.
- *Rational and real heights.* Ramification multiplies every quantity
  in (2) by the same integer. For real heights, the proof applies
  verbatim to the determinantal-prefix definition of the weighted
  profiles. Each weighted minor valuation is an affine function of the
  heights, and the support bounds hold entrywise with real exponents.
- *Containment.* A product set `S=A x B` has
  `Cell(S)=|B| X_A+|A| Y_B`, which recovers the rectangular caps (5) of
  the predecessor note. Cell sets of non-product shape are new. A
  line `A x {0}` gives only `kd+X_k`, which is twice the shared-line
  bound. The shared lines (4) therefore remain separate constraints.
- *Sharpness.* Take only matched monomial entries, with every infinity
  order zero. Then `d=0`, `s=0`, and `r` is the sorted list of cell
  energies. Every inequality in (2) on `r` then holds with equality.
  The checker builds 24 such tensors from random heights. A single cell
  cap therefore cannot be lowered, although the two profiles need not
  reach their caps together.

## 3. Which minors are shared, and what single minors can give

**Proposition 2.**

1. *Endpoint minors.* Over the strict endpoint, the `k`-square minors
   of a cut whose endpoint value is nonzero are exactly the cell-set
   minors of Theorem 1, up to sign. The support bounds (1) bound every
   determinant term of `m_S` by the same number, `kd+Cell(S)`. Hence no
   better degree bound for a fixed minor can be derived from (1) alone.
   Within the class of arguments of the form
   *prefix ≤ order of one endpoint-nonzero minor ≤ its support degree*,
   the per-profile caps (2) are best possible.
2. *Literally shared minors.* Let a `k`-square array of entries
   `R(x_i,y_i,u_j,v_j)` with `k>=2` be a submatrix of both cuts, allowing
   transposition. Then either `y_i` and `u_j` are constant, or `x_i` and
   `v_j` are constant. If the minor is nonzero at the endpoint, the two
   constants coincide. These are exactly the shared lines of the
   predecessor's (4), and the cheapest choice of constant is the
   baseline symbol.

*Proof.* (1) In each cut the endpoint is the permutation matrix of
`(x,y) -> (y,x)` or `(x,u) -> (u,x)`. A minor of a permutation matrix is
nonzero exactly when its column set is the image of its row set; it is
then `±1`. The uniform term bound was shown in Theorem 1.

(2) Suppose the array is a submatrix of the column cut with the same
row index `i`. The row label `(x_i,u_j)` must then depend only on `i`.
For `k>=2` this forces `u_j` to be constant. Likewise the column label
`(y_i,v_j)` must depend only on `j`, which forces `y_i` to be constant.
If instead the array is transposed, the same reasoning forces `x_i` and
`v_j` to be constant. At the endpoint the entries are
`[x_i=v_j][y_i=u_j]`. With `y` and `u` fixed, a nonzero minor needs
`y=u`, and symmetrically in the other case. A line at symbol `b` costs
`kd+X_A+k c_Y(b)`, which is minimized at `b=0`. QED

**Lossless averaging.** Every constraint used here is symmetric in `r`
and `s`, or bounds `r+s`. Any `gamma` satisfying the averaged caps
therefore gives the feasible pair `r=s=gamma`. Conversely, every feasible
pair satisfies `log Z(r)+log Z(s)>=2 log Z(gamma)`. The averaged problem
of Section 4 consequently has the same minimum as the two-profile
problem with the same constraints. Relations among *different* minors
are not used. These include Plücker relations, Laplace expansions, and
the partial-transpose relation between the two cut matrices on the
same entry set. Such relations could exclude further points and remain
open.

## 4. The strengthened relaxation

We take the predecessor's caps, namely `S_0=S_1=0`, the shared lines (4)
and the rectangles (5), and add (2) for every `1<=k<=q^2`. As in the
predecessor's Section 3, `C(k)` denotes the greatest convex minorant of
the cap points, and `gamma*` its slopes. The argument there,
`sum_i (u-gamma_i)_+ = max_k(ku-S_k)` integrated against `e^-u`, gives
`Z(gamma)>=Z(gamma*)`. Hence the complete pressure (3) is at least

\[
G_{\rm relaxed}=d+2\log\sum_i e^{-\gamma_i^*}
-\sqrt2\Bigl(\log\sum_x e^{-c_X(x)/\sqrt2}+\log\sum_y e^{-c_Y(y)/\sqrt2}\Bigr).
\tag{3}
\]

The caps are nondecreasing in `d` and in every paired energy. Larger
caps give a pointwise larger minorant and so a smaller `Z(gamma*)`. This
monotonicity justifies the corner bound of Section 7.

## 5. Exact verdicts

The [checker](../tools/check_cell_minor_caps.py) builds every cap with
`fractions.Fraction`, forms the exact hull, and bounds (3) with 120-bit
outward dyadic arithmetic. Exponentials use Taylor remainders. Logarithms
use binary argument reduction with an `atanh` series and remainder.
`sqrt(2)` uses integer square roots. A sign is accepted only when the
whole interval lies on one side of zero. The intervals have width below
`10^-30`.

The q=4 predecessor target has the hull vertices
`(0,0),(1,0),(2,1),(4,4),(8,16),(12,92/3),(13,107/3),(15,47),(16,160/3)`.
The vertices at `k=8,12,13,15` are cell caps of non-product cell sets.
The sixteen cell energies are `0` and `2` and `8/3` four times each,
`4` once, `14/3` twice and `16/3` once.
The average orders become
`0;1;3/2^2;3^4;11/3^4;5;17/3^2;19/3`, and `G_relaxed` rises to
`+0.0430257`. **The q=4 target of the predecessor note is excluded.**

The q=5 predecessor target remains negative, at `-0.0042091`.

The two sharper targets in the table were located by the labelled
search of Section 6 and are certified exactly:

| q | d | Paired energies (both sides) | Average orders `gamma*` | `G_relaxed` |
|---|---|---|---|---|
| 5 | 7/6 | 0,0,19/9,23/9,37/9 | 0; 7/6; 59/36; 67/36; 95/36; 377/108 ×3; 67/18 ×4; 95/18 ×4; 97/18; 35/6 ×2; 113/18; 133/18 ×2; 47/6 ×2; 169/18 | -0.01905094 |
| 6 | 8/7 | 0,0,15/8,27/10,4,4 | listed in the evidence record | -0.04361177 |

The binding families at these points, and at the q=4 and q=5
predecessor points, are the shared lines and the cell sets. A rectangle
cap equals the cell cap of its product set, so rectangles are redundant
once cell caps are present; the checker labels ties as cell sets. The common denominators of the data and
orders are 108 and 1120. After multiplying every height by that integer,
all target orders and energies are integral, with compensating scale
`t=1/108` or `1/1120`.

These are construction targets: a realizing tensor needs simultaneous
profiles `r=s=gamma*` or better, with the stated energies. The
[complete-entry affine solver](full-alphabet-two-cut-affine-synthesis.md)
is the existing interface for finite frames and windows.

## 6. Exploration (floating point; not a certificate)

The scripts minimize (3) with double precision. They use a coordinate
pattern search over `d` and the nonbaseline energies, with random and
two-zero-symbol starts, both symmetric (`c_X=c_Y`) and asymmetric. The
recorded values are as follows.

| q | Best value found | Approximate minimizer | Searches |
|---|---|---|---|
| 3 | positive at every point found; infimum 0 | values below `10^-4` occur only with all energies near 50 and `d=0` (for example, `G≈4e^(-t/2)` along `(0,t,t)`); finite local minimum `≈0.0789` at `d=0`, `(0,1.65,4.95)` | 4,037 local searches |
| 4 | no negative value; infimum 0 in degenerate tails | smallest interior local minimum **+0.00509** at d≈1.198, both sides ≈(0,0,2.395,2.395); `G(d=0,(0,T,T,T))` falls to about `2.7e-4` at `T=20` | 2,429 local searches |
| 5 | -0.01907 | d≈1.165, both sides ≈(0,0,2.115,2.545,4.066) | over 1,700, including asymmetric starts |
| 6 | ≤ -0.04367 | d≈1.14, both sides ≈(0,0,1.86,2.70,4.03,4.03) | about 50; not converged globally |

Pinning some energies at `60`, which effectively deletes those letters,
gives the rectangular limits. Searches over the `4x1`, `4x2`, `4x3`,
`3x3`, `2x4`, `3x1`, `3x2`, `5x1` and `5x2` alphabets are all positive.
Their smallest values are either about `0.05` to `0.1`, or tend to zero
only as the remaining energies also diverge.

With the predecessor's caps alone, the research lead's corresponding
q=4 probe finds about `-0.0550`, at `d≈1.159` and `(0,0,2.240,2.397)`. If
`G_relaxed>0` everywhere at `q=4`, the strengthened inequalities would
exclude every gauged four-symbol seed under this criterion. That needs a
proof, which is not supplied here. Because sending one letter on each
side to infinity turns the q=4 relaxation into the q=3 one, such a proof
contains the whole ternary closure of Section 7, plus tail lemmas for
every rectangular sub-limit up to `3x4` and `4x3`, plus an interior bound
with margin about `0.005`. The smallest gauged alphabet where these
inequalities certainly leave room is `q=5`.

## 7. Ternary status

Every limit in which some energies diverge is a relaxation for a
smaller rectangular alphabet. The one-letter limit gives `G=d>=0`, so
the q=3 infimum is zero and is not attained. A closure proof must
therefore control nested tails. The template is the binary reduction of
the [universal binary floor](universal-binary-smith-floor.md),
Sections 4–6.

**Proved here (exact).** For `q=3`, every point with all paired
energies in `[0,1]` and any `d>=0` has `G_relaxed>=1/100`.

- For `d in [0,4]`, the checker covers the region with 15,000 dyadic
  boxes. On each box it uses the monotone corner bound: caps at the
  upper corner, and the deficit and physical term at the lower corner.
- Boxes are discarded only outside the fundamental domain of sorting
  and of the exchange `X<->Y`.
- An audit re-derives the split tree and checks that every saved leaf
  is reached exactly once. Each leaf is then re-evaluated at a second
  precision.
- For `d>=4`, `Z(gamma*)>=1` and each physical partition is at most
  `3`, which gives `G>=4-2sqrt2 log 3>0.89`.

**What remains for a ternary theorem.**

1. *All energies large.* With the line caps at `k=2,3`, the leading
   terms are `2(e^-gamma_2+e^-gamma_3)`, with
   `gamma_2<=d+min(a_1,b_1)/2`. The physical term is
   `sqrt2 sum e^(-c/sqrt2)`. Since `1/2<1/sqrt2`, this tail looks
   provable in the binary manner, but no proof is written.

   The mechanism is a norm comparison. If one side has a single letter
   and `d=0`, the caps give `gamma*=c/2` exactly. `G_relaxed` is then
   `log||x||_(1/2)-log||x||_(1/sqrt2)>=0` with `x_i=e^(-c_i)`, by
   monotonicity of `l^p` quasi-norms in `p`. A proof for general `d`
   and two-sided alphabets must control the cross terms between the
   two sides.
2. *Mixed tails.* Some energies are large and others are not. These
   reduce to the rectangular `2x3`, `1x3`, `2x2` and `1x2` relaxations,
   each of which needs a uniform floor.
3. *The compact interior.* The corner bound loses at first order in
   the box width. Exact runs needed 2,683 leaves at energy edge `1/2`
   and 20,162 at edge `1` before the `X<->Y` reduction (15,000 after
   it). An independent floating-point reimplementation of the same split
   rule and domain reduction finishes edge `2` with about 125,500 leaves,
   so the leaf count grows about eightfold per doubling of the edge. In
   any case no fixed-margin cover can reach the tails, where
   `G_relaxed` tends to zero. A different method is therefore needed: a two- or
   three-parameter analytic reduction as in the binary Section 5, or
   derivative-aware interval bounds. Pushing this cover outward is not
   a viable route.

## 8. Literal control

Part 1 of the checker is independent of the proof. It works over
`F_10007` with random integer heights in `[-2,2]` on all four legs at
both places. It builds random sparse Laurent strict crossings, matched
monomial crossings, and the tight family. For each tensor it:

- recomputes `mu_0`, `mu_infinity`, `d` and the paired energies from
  the raw heights;
- checks the degree bounds (1), and that the infinity matrix is the
  reciprocal normalization;
- computes all four complete profiles by exact power-series elimination
  (both cut assignments, at both places);
- checks every cell, shared-line and averaged-hull prefix;
- expands random cell-set minors by exact interpolation, checking their
  endpoint value, their degree bound, and the reciprocal identity at
  infinity.

The recorded run covers 104 tensors (q=2,3,4), 416 complete profiles and
352 minors. It records 502 tight cell prefixes and 131 tight
shared-line prefixes, and every tight-family tensor attains all its zero
cell caps. These counts include trivially tight `k=1` prefixes, and an
instrumented replay shows that almost all tight cell prefixes (494) come
from the purpose-built `d=0` tight family. The random tensors have `d`
between 2 and 12, and their caps are loose. The control therefore checks
the transcription of the caps and the profile computation; it does not
exercise sharpness at `d>0`. (A referee construction, matched `w^N`
entries plus one unmatched `1-w^N` entry, attains the full-determinant
cap `q^2 d` at `q=2,3`, so the `kd` term is sharp at `k=q^2`.) As
mutation controls, lowering every relaxation cell cap by `1/10`, or
every literal cell cap by `1` (detected through the tight family), makes
the checker exit with failure. The proof is Theorem 1.

## Reproduction

```sh
make verify-cell-minor-caps
```

This runs
[`check_cell_minor_caps.py`](../tools/check_cell_minor_caps.py) in about
30 seconds. The ternary cover takes most of that time; the literal
control and verdicts alone take under 2 seconds with `--skip-cover`.
The record is [`cell-minor-caps.json`](../evidence/cell-minor-caps.json).
The tool uses the standard library and
`tools/charged_gaussian_floor.py`; `numpy`, when present, only
accelerates series products, and both paths give the same record hash.
It does not import the predecessor checker. In an environment without
`scipy`, that checker's import chain fails before running, and so does
`make verify-shared-minor-targets`. The predecessor's two pressures are
recomputed here with the cell caps switched off.
