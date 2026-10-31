# The gauged binary Smith floor extends to auxiliary matrix kernels

## Result and scope

Let `F` be any field and `D>=1`. Let

\[
K(w;x,y,u,v)\in\operatorname{Mat}_D\bigl(F[w,w^{-1}]\bigr),
\qquad x,y,u,v\in\{0,1\},
\]

be a Laurent-polynomial auxiliary (matrix-valued) binary seed with a
nonzero constant vector `e` such that

\[
K(1;x,y,u,v)\,e=[x=v][y=u]\,e .
\tag{1}
\]

This is the endpoint of the [Smith-profile transfer](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/smith-profile-crossing-transfer.md)
of PR #30; a boundary covector `f` with `f^T e=1` then exists. Use the
cost cuts `(x,y)|(u,v)` at zero and `(x,u)|(y,v)` at infinity, or the
swapped assignment. In each cut both auxiliary indices sit on one
physical side; the side may differ between the two cuts. Allow:

* every real physical-letter height on every leg at both places;
* an independent auxiliary conjugation gauge at each place,
  `K -> G_0 K G_0^{-1}` at zero and `K -> G_inf K G_inf^{-1}` at
  infinity, with arbitrary `G_0, G_inf` in `GL_D(F(w))`. These are the
  "auxiliary rational gauges" permitted in
  [type-optimized pressure](type-optimized-crossing-pressure.md) section 2
  and [fixed blocks](block-smith-crossings.md) section 3.

Let `G_aux(t)` be the independent paired-type critical pressure of the
[universal binary floor](universal-binary-smith-floor.md), equation (4),
computed from the four sorted Smith orders of the two rectangular gauged
cut matrices (`4D^2`-by-`4` when the auxiliary pair sits with the row
side of the cut, `4`-by-`4D^2` when it sits with the column side). Then
for every `t>0`

\[
\boxed{G_{\rm aux}(t)\ \ge\ G_{\widetilde S}(t)\ >\ \frac{e^{-\max(a,b)t}}{28}\ >\ 0,}
\tag{2}
\]

where `a,b>=0` are the canonical paired energies of that note's (5), and
`S~` is a scalar binary Laurent-polynomial strict crossing built from `K`
by a twisted boundary contraction (Section 3). Consequently **no
auxiliary binary seed, of any fixed auxiliary dimension, over any field,
with any physical heights and any auxiliary rational gauges, can improve
the factoring exponent through this two-place Smith criterion.** The
constant `1/28` of the scalar theorem is kept; no weaker constant is
needed.

Without auxiliary gauges, or with one common conjugation
`G` in `GL_D(F[w,w^{-1}])` (a Laurent matrix with Laurent inverse) at
both places, the hypothesis (1) can be weakened to: some constant linear
functional `phi` on the auxiliary pair index, over `F` or an extension,
satisfies `phi(K(1))=[x=v][y=u]` (Section 3). A common conjugation by a
merely rational `G` is a special case of independent gauges and needs
(1).

The proof is a reduction, not a new estimate. A one-sided contraction
with weights regular at a place cannot lower any Smith order there,
coordinatewise ([fixed blocks](block-smith-crossings.md) section 2), so
the auxiliary partition function is at least that of the contraction.
The only point where the auxiliary case needs a new idea is that with
independent gauges no **constant** contraction need be regular at both
places at once. The explicit seed of Section 5 shows that every constant
contraction can fail the per-place domination, and that the boundary
contraction `f^T K e` of a legitimate boundary pair can have strictly
larger pressure than the auxiliary seed. A Laurent twist of the boundary
vectors, obtained from Grothendieck's splitting theorem (explicitly for
diagonal gauges), repairs this.

What is **not** covered:

* block certificates of several binary positions: a block is an auxiliary
  seed on alphabet `2^r>=4`, and the reduction only lands on scalar
  four-or-more-letter seeds, where the gauged question is open (see
  [shared-minor profile relaxation](shared-minor-profile-relaxation.md));
* cancellation-module or minimal-boundary-module word-rate bounds
  ([cancellation modules](cancellation-module-smith-automata.md),
  [minimal modules](minimal-boundary-module-completeness.md)), which are
  not local Smith profiles;
* auxiliary row weights that are not conjugations (they do not telescope
  in products, and no transfer is known for them), and auxiliary indices
  split across a cut;
* rational auxiliary kernels with poles at `w=1`; the
  [typed rational strictification](typed-rational-crossing-strictification.md)
  is scalar;
* multipoint criteria and actual coefficient ranks.

This is an obstruction for one sufficient certificate. It is not a
factoring lower bound. The established factoring exponent remains `1/2`.

Predecessors: the scalar theorem and its 421-leaf compact certificate
are PR #39's [universal binary floor](universal-binary-smith-floor.md)
(source revision
[a3c5bf15](https://github.com/1337isnot1337/integer-factoring/blob/a3c5bf15fcc80ea8ab9af7b4853a1a70ae6a88b0/research/universal-binary-smith-floor.md)).
The coordinatewise contraction inequality is section 2 of the in-tree
[fixed blocks](block-smith-crossings.md) note, an instance of the
classical fact that the invariant factors of `A` divide those of `CA`.
The Cauchy--Binet treatment of bare (ungauged) auxiliary cuts is from
[PR #40](https://github.com/1337isnot1337/integer-factoring/pull/40) at
revision e3dadb07: section 1 of
[shared-minor-smith-floor.md](https://github.com/1337isnot1337/integer-factoring/blob/e3dadb07ea81424d20ff933eaa3d2e8462c5bf8a/research/shared-minor-smith-floor.md)
and of
[all-degree-small-alphabet-smith-floor.md](https://github.com/1337isnot1337/integer-factoring/blob/e3dadb07ea81424d20ff933eaa3d2e8462c5bf8a/research/all-degree-small-alphabet-smith-floor.md).
The auxiliary interface is PR #30's. Two independent referees of this
work supplied the Section 5 seed and the conjugated controls of
Section 6.
This note was developed with substantial Anthropic Claude assistance. It
is an ordinary written proof with exact finite controls; there is no new
Lean formalization.

## 1. The auxiliary two-place pressure

Write `lambda=(x,y,u,v)` and let `alpha` range over the `D^2`
auxiliary index pairs. For a place `P` in `{0,infinity}`, with local parameter `w` at zero
and `tau=1/w` at infinity, the gauged rectangular cut matrix has rows
`(physical row pair, alpha)` and four physical columns (the column-side
case is its analogue with columns `(physical column pair, alpha)`; it is
handled by transposition throughout), with entries

\[
M_P[(r,\alpha),c]=\pi_P^{\sum_l h_{P,l}(\lambda_l)}
\bigl(G_PKG_P^{-1}\bigr)_\alpha(\lambda),
\qquad \pi_0=w,\ \pi_\infty=\tau,
\tag{3}
\]

with `K(lambda)` read at `tau^{-1}` at infinity. The row physical pair is
`(x,y)` in the zero cut and `(x,u)` in the infinity cut. The physical
weight is a product of a row factor and a column factor and does not
depend on `alpha`.

Over the discrete valuation ring of the place, `M_P` has generic rank at
most four. Let `delta_k(M_P)` be the minimum valuation of its nonzero
`k`-by-`k` minors. The Smith orders are
`sigma_{P,k}=delta_k-delta_{k-1}`, nondecreasing in `k`. For real
heights, `delta_k` is the minimum over minors of an affine function of
the heights, as in the scalar note section 6. It is continuous, so the
sortedness proved for rational heights by ramification persists. The
pressure `B(t)` and `G_aux(t)` are the scalar note's (3)--(4) with these
four orders. Section 3 shows that the rank is four.

## 2. Contraction monotonicity

**Lemma 1 (regular contraction).** Let `M` be a matrix over `F((pi))`
whose rows are indexed by pairs `(r,alpha)`. Let
`Phi[r',(r,alpha)]=[r'=r] phi_alpha` with every `phi_alpha` in
`F[[pi]]`, and suppose `Phi M` and `M` have the same generic rank `n`
(in the application `n=4`, the number of physical columns).
Then their sorted Smith orders satisfy

\[
\sigma_i(\Phi M)\ \ge\ \sigma_i(M)\qquad(1\le i\le n),
\tag{3a}
\]

and hence `delta_k(Phi M)>=delta_k(M)` for every `k`. The same holds
after multiplying the rows of `M` and of `Phi M` by `pi^{rho(r)}`
(depending only on the physical row) and all columns by `pi^{kappa(c)}`,
with real weights. It also holds for column contractions `M Phi^T`, with
the roles of rows and columns exchanged.

*Proof.* Because `Phi` acts trivially on the physical row label,
multiplying by the physical row weights commutes with `Phi`; column
weights commute with any left multiplication. For integer weights,
multiply by a common power of `pi` so that `M` is regular; then (3a) is
equation (3) of [fixed blocks](block-smith-crossings.md) section 2:
multiplication by `Phi` maps the image of `M` modulo `pi^l` onto the
image of `Phi M` modulo `pi^l`, so `Phi M` has at most as many orders
below `l` as `M`, for every `l`. For rational weights, replace `pi` by
`varpi^N`; `F[[varpi]]` is again a discrete valuation ring and `Phi`
stays regular. For real weights, each `delta_k` is the minimum, over the
weight-independent set of nonzero `k`-minors, of affine functions of the
weights, so every `sigma_i=delta_i-delta_{i-1}` is continuous in the
weights and (3a) passes from rational to real weights. Column
contractions follow by transposition. `QED`

**Corollary.** Suppose a scalar tensor `S` has, at each place, a cut
matrix of the form `Phi_P M_P` (or `M_P Phi_P^T`) as in Lemma 1, with
the same physical heights, and has rank four. Then the auxiliary
matrices also have rank four (by Cauchy--Binet a nonzero `4`-minor of
`Phi_P M_P` is a combination of `4`-minors of `M_P`, which has only four
physical columns, or rows in the column-side case), and each auxiliary order is at most the
corresponding order of `S`, so `B_aux(t)>=B_S(t)` term by term. The
physical term of the critical pressure depends only on the heights.
Therefore `G_aux(t)>=G_S(t)` for all `t>0`.

**Remark (prefix version).** Cauchy--Binet gives only the weaker prefix
inequalities `delta_k(Phi M)>=delta_k(M)`. These already suffice for the
pressure comparison, by the Tomic--Weyl inequality: if
`x_1<=...<=x_n`, `y_1<=...<=y_n` and
`sum_{i<=k} x_i <= sum_{i<=k} y_i` for every `k`, then
`sum_i e^{-t x_i} >= sum_i e^{-t y_i}` for `t>0`. (Let `c_i` be the
divided difference of `g(z)=e^{-tz}` at `x_i,y_i`; then `c_i<=0` is
nondecreasing by convexity, and Abel summation with
`X_k=sum_{i<=k}(x_i-y_i)<=0` gives
`sum_i(g(x_i)-g(y_i))=c_nX_n+sum_{i<n}(c_i-c_{i+1})X_i>=0`.) The prefix
form is the one available for criteria that only record determinantal
prefixes; the main proof uses (3a).

## 3. A scalar contraction regular at both places

The contraction must be one scalar tensor `S` for **both** places. This
keeps the shared-minor coupling of the scalar proof. It must also be
regular in Lemma 1's sense at each place for that place's gauge.

**No auxiliary gauge.** Take `S=f^T K e`, or more generally `phi(K)`
with `phi(K(1))=E`. The weights are constants, and `S` is a scalar
Laurent seed with `S(1)=E`. If `phi` is defined over an extension
`L/F`, nothing changes: determinantal valuations do not change under
field extension, and the scalar theorem holds over every field.

**One common Laurent conjugation.** Let `G` in `GL_D(F[w,w^{-1}])` act at
both places. Then `K'=GKG^{-1}` is again a Laurent seed, and the
auxiliary profiles at both places are those of `K'` without gauge. The
determinant of `G` is a unit of `F[w,w^{-1}]`, a monomial `c w^k` with
`c != 0`, so `G(1)` is invertible. Put
`phi'(X)=phi(G(1)^{-1} X G(1))`. Then
`phi'(K'(1))=phi(K(1))=E`, and the ungauged case applies to `K'` and
`phi'`. If `G` is only rational, `K'` can have poles away from
`0,infinity` (for example at `w=1`), and this argument does not apply;
such a common gauge is a special case of the independent gauges below,
which need (1).

**Diagonal monomial gauges.** Let `G_0=diag(w^{gamma_c})` and
`G_inf=diag(w^{eta_c})`, both written in `w`. The gauged entry `(c,b)`
is `w^{gamma_c-gamma_b}K_{cb}` at zero and `w^{eta_c-eta_b}K_{cb}` at
infinity. Choose `c^*` in the support of `e` that maximizes
`eta_b-gamma_b` over that support, and put

\[
\widetilde S=\sum_b \frac{e_b}{e_{c^*}}\,w^{\gamma_{c^*}-\gamma_b}\,K_{c^*b}.
\tag{4}
\]

At zero the weight of the gauged entry `(c^*,b)` is `w^0`. At infinity it
is `w^{(eta_b-gamma_b)-(eta_{c^*}-gamma_{c^*})}`, with a nonpositive
exponent, so it is regular in `tau`. By (1), `S~(1)` is
`(K(1)e)_{c^*}/e_{c^*}=[x=v][y=u]`.

**General rational gauges.** Let `A=F[w,w^{-1}]`. For an integer `n`,
let `Gamma(n)` be the set of `v` in `A^D` with `G_0 v` regular at zero
and `w^{-n}G_inf v` regular at infinity. Let `Gamma*(n)` be the set of
`phi` in `A^D` with `G_0^{-T}phi` regular at zero and
`w^{n}G_inf^{-T}phi` regular at infinity. These are the global sections
of a rank-`D` vector bundle `V(n)` on the projective line and of its dual
`V^*(-n)`. The bundle is trivial on `w != 0, infinity`, so its fiber at
`w=1` is `F^D`. By Grothendieck's splitting theorem, valid over every
field (for example Hazewinkel and Martin, *J. Pure Appl. Algebra* 25
(1982)), `V` is isomorphic to `O(m_1)+...+O(m_D)`. Let `L_i` be the
fiber lines at `w=1` and put `Phi_{>=k}` equal to the sum of the `L_i`
with `m_i>=k`. Evaluation at one maps `Gamma(n)` onto `Phi_{>=-n}`, and
maps `Gamma*(n)` onto the annihilator of `Phi_{>=-n+1}`.

Since `e != 0`, there is a largest `k` with `e` in `Phi_{>=k}`; take
`n=-k`. Then some `e~` in `Gamma(n)` has `e~(1)=e`. Also `e` is not in
`Phi_{>=k+1}`, so some `f~` in `Gamma*(n)` satisfies `f~(1)^T e=1`. Put

\[
\widetilde S=\widetilde f^{\,T}K\,\widetilde e
=(G_P^{-T}\widetilde f)^T\,(G_PKG_P^{-1})\,(G_P\widetilde e)
\quad(P=0,\infty).
\tag{5}
\]

At zero, both outer vectors are regular. At infinity, insert
`w^{n}w^{-n}`: both twisted factors are regular. So the products
`f-hat_c e-hat_b` are regular weights at each place. `S~` is a scalar
binary Laurent-polynomial tensor. By (1),
`S~(1)=f~(1)^T K(1) e = [x=v][y=u]`.

`S~` therefore satisfies every hypothesis of the scalar theorem. Its
cut matrices have rank four, and by Lemma 1 so do the auxiliary ones.
The case `K(1)e=(E+N)e`, with `N` of strictly positive paired charge,
gives `S~(1)=E+N`. The scalar note's section 7 then applies unchanged.

## 4. Proof of (2) and audit of the scalar steps

Section 3 and the Corollary give `G_aux>=G_S~`. The
[universal binary floor](universal-binary-smith-floor.md), (29), applies
to `S~` with the same heights, so `G_S~(t)>e^{-max(a,b)t}/28`. Rational
and real heights, and the swapped cut assignment by `w -> 1/w`, are
handled as in that note. This proves (2).

The task of checking each scalar step separately leads to the same
conclusion. It also shows exactly where scalarity entered.

| Scalar step (universal note) | Auxiliary status |
| --- | --- |
| (5) pair flips, letter-0 shifts | relabel physical rows and columns only; unchanged |
| (6)--(8) normalization by the minimum orders `mu_0, mu_inf` and entry supports | under independent gauges the gauged entries at the two places are different Laurent polynomials, so (8) has no entrywise auxiliary meaning (without gauges it does); it is used only for `S~`. The auxiliary minima are `mu_P^aux=sigma_{P,1}^aux<=mu_P^{S~}` by (3a) with `i=1`, so `d_aux>=d_{S~}`, and `S~` obeys (8) also with the auxiliary normalization |
| `d>=0` from `R_0000(1)=1` | `d_aux>=d_{S~}>=0`, using the matched entry of `S~` |
| (11)--(12) determinant budgets | `delta_4^aux(0)<=delta_4^{S~}(0)<=-4(a+b)-4mu_inf^{S~}<=-4(a+b)-4mu_inf^aux`, by (3a), then (12) for `S~`, then (3a) with `i=1` at infinity; this is `sum r_i^aux<=C_aux=4d_aux+4(a+b)`. The infinity budget is symmetric |
| (13)--(15) shared `2x2` line minor | the shared minor is a minor of `S~`, not of the auxiliary matrices. By (3a), `delta_2^aux(P)<=delta_2^{S~}(P)` at both places, and (14)--(15) for `S~` give `delta_2^{S~}(0)+delta_2^{S~}(inf)<=2 sum_l H_l+2min(a,b)`; subtracting `2mu_0^aux+2mu_inf^aux` gives `r_2+s_2<=2d_aux+2min(a,b)` |
| endpoint nonvanishing | used only through `S~(1)` |
| (17)--(29) convexity, balancing, compact cover | profile-only; unchanged |

Here `r,s,d,C` in the auxiliary column are the scalar note's quantities
(6), (10), (11) computed from the auxiliary minimum orders. The one
hidden scalar assumption is that a single contraction is regular at both
places. Section 5 shows that no constant contraction need have this
property under independent gauges. The direct audit, with `d` defined by
the auxiliary minimum orders, is checked on every control below.

## 5. Constant contractions can fail under independent gauges

Over any field `F` take `D=2`, `e=(1,1)`, all physical heights zero and

\[
K=\begin{pmatrix} 0 & w^5E\\ E & 0\end{pmatrix},
\qquad G_0=\operatorname{diag}(1,w^3),\qquad G_\infty=\operatorname{diag}(1,w).
\]

Then `K(1)e=Ee`. In either cut, the four-by-four cut matrix of `E` is a
permutation matrix, so a scalar seed `c(w)E` has four equal orders
`ord_P c` at each place. The gauged entries are `w^2E` and `w^3E` at zero
and `w^4E` and `wE` at infinity. Hence, for both cut assignments and
with the auxiliary pair on either side, the auxiliary orders are
`(2,2,2,2)` at zero and `(-4,-4,-4,-4)` at infinity, and
`B_aux(t)=16e^{2t}`.

*Every constant contraction fails the per-place domination.* Let `phi`
be any constant functional on the auxiliary pair index, over any
extension of `F`, with `phi(K(1))=E`. Then
`phi(K)=(phi_12 w^5+phi_21)E` and `phi_12+phi_21=1`. If
`phi_21 != 0`, the zero orders of `phi(K)` are `0<2`. Otherwise
`phi(K)=w^5E` has infinity orders `-5<-4`. Either way the conclusion of
Lemma 1, even in prefix form, fails at one place. So with independent
gauges a constant boundary contraction cannot replace the twist in the
proof.

*A boundary contraction with strictly larger pressure.* If `F` has at
least three elements, take `f=(lambda,1-lambda)` with
`lambda` not in `{0,1}`, for example `f=(2,2)` over `F_3`. Then
`f^T e=1`, and `f^TKe=(lambda w^5+1-lambda)E` has orders `(0,0,0,0)`
and `(-5,-5,-5,-5)`, so `B(t)=16e^{5t}`. The physical terms agree, so

\[
G_{f^TKe}(t)-G_{\rm aux}(t)=3t>0\qquad(t>0).
\]

The boundary contraction of a legitimate boundary pair therefore does
not bound the auxiliary pressure from below.

*The twist.* Formula (4) has `c^*=1` and gives `S~=w^2E`, with orders
`(2,2,2,2)` and `(-2,-2,-2,-2)`. These dominate the auxiliary orders
coordinatewise, and `G_aux(t)=G_{S~}(t)+2t`.

*What this does not show.* The two pure functionals `phi=e_12` and
`phi=e_21` give `w^5E` and `E`, whose pressure equals `G_{S~}`. So in this
seed some constant contraction still bounds `G_aux` from below at the
pressure level, although none dominates per place. Whether some seed
forces a nonconstant twist at the pressure level, against every constant
functional, is not settled here. The note's theorem does not depend on
this question.

The checker recomputes all of these profiles. Over `F_2` and `F_3` it
enumerates every constant `phi` with `phi(K(1))=E` (8 and 27 of them)
and confirms that each fails the per-place domination for every cut
assignment and auxiliary side. It confirms the exact gap `3t` for
`f=(2,2)` over `F_3`, and that both formula (4) and the general solver
of Section 3 give a dominating `S~`.

## 6. Exact controls

The [checker](../tools/check_auxiliary_binary_gauge_floor.py) works over
`F_2`, `F_3`, and `F_5` with `D=1,2,3`. It uses random Laurent seeds
satisfying (1) literally and random integer heights in `[-3,3]` for both
letters on every leg at both places. Gauges are independent at the two
places and are of four kinds:

* absent;
* diagonal monomial;
* framed: `diag(w^g)(I+w^jN)P`, with `N` nilpotent and `P` an invertible
  constant matrix;
* conjugated (`D=2`, over `F_2` and `F_3`, heights in `[-2,2]`):
  `K=V[[A,0],[E,0]]V^{-1}` with `A` a random scalar strict crossing,
  `e=V(1)(1,1)`, and `G_0=U_0 diag(1,w^g) V^{-1}`,
  `G_inf=U_inf diag(1,w^h) V^{-1}`, where `U_0` and `U_inf` are unipotent
  and regular at their places, `g>=4` and `h<=-4`. This family follows a
  referee's construction. It makes the constant contraction fail under
  non-diagonal gauges, so that the general twist is exercised where it
  is needed.

There are 230 seeds: 210 random ones and 20 conjugated ones. For each
seed, both cut assignments, and (for `D>=2`) each of the four choices of
auxiliary side at zero and at infinity, the checker:

* computes the rectangular Smith orders exactly by local elimination over
  `F_p[[w]]` with adaptive precision. On 96 rectangular control profiles
  (the twelve `D=2` seeds over `F_2` and `F_3`), and on the 32 profiles of
  the Section 5 seed, it also recomputes them from the definition, by
  minimum valuations of all minors;
* finds `e~,f~` by exact linear algebra over a Laurent window, and checks
  their values at one, the regularity of every weight at both places,
  the identity (5), and that `S~` is a strict crossing. All 920 `S~`
  profiles are computed both by elimination and by all minors;
* checks the coordinatewise domination (3a), hence the prefix
  domination, and compares the partition sums exactly at
  `e^{-t}=1/2, 9/10, 99/100`;
* checks the normalized inequalities `d>=0`, `sum r,sum s<=C`, and
  `r_2+s_2<=2d+2min(a,b)` directly on the auxiliary profile;
* checks the outward dyadic bound `G_aux(t)>=Phi(dt,at,bt)>0` of the
  scalar (18), and `G_S~(t)>0`, at `t=1/10,1/2,1,2`. The stored values
  of `Phi` are the certified lower endpoints rounded down to a multiple
  of `10^-6`;
* compares the constant contraction `f^TKe` with the auxiliary seed, per
  place and exactly at the pressure level (the physical terms cancel, so
  this is an exact rational comparison of `B(t)` at the three values of
  `e^{-t}`).

It replays the stored 421-leaf certificate
[binary-fixed-ray-floor.json](../evidence/binary-fixed-ray-floor.json):
it checks that the certificate's domain is `[0,3/2]^2` and its margin
`1/16`, evaluates every leaf afresh with outward arithmetic, and rejects
a deleted leaf, an altered bound, and an enlarged domain. The hash of
that file is recorded. It also recomputes the Section 5 seed.

Recorded outcomes are in the
[evidence record](../evidence/auxiliary-binary-gauge-floor.json):

* all 1660 (seed, assignment, side) combinations pass every check;
* in 1449 of them, all with `D>=2`, the auxiliary orders are strictly
  below those of `S~` at some place;
* in the 120 random diagonal or framed gauged seeds the constant
  contraction happened to dominate per place every time;
* in the conjugated family it fails per place in 104 of 160
  combinations, and has strictly larger pressure at some tested `t` in
  68. In 9 of the 20 conjugated seeds this pressure failure occurs with
  non-diagonal gauges and a nonconstant twisted pair; the twisted `S~`
  dominates in every case.

Scope of the controls. They use integer heights, Laurent gauges with
Laurent inverses, and the hypothesis `K(1)e=Ee`. Real heights, merely
rational gauges, the weakened hypothesis for common conjugations, other
fields, and larger `D` are covered by the written proof only. Random
seeds lie far above the floor. These controls check the transcription
and the constructions, not the general theorem, which is the argument
above.

## Reproduction

```sh
make verify-auxiliary-binary-floor
python3 tools/check_auxiliary_binary_gauge_floor.py --output /tmp/auxiliary-binary-replay.json
```

The run takes about 45 seconds and needs only the Python standard
library. Because of its running time it is not part of the default
`make verify`. PR #39's
`make verify-binary-smith-floor` imports SciPy through a search module.
This checker instead replays that target's compact certificate directly.
