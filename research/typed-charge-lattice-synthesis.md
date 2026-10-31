# Complete-entry rational synthesis with paired-type crossing endpoints

## Result and scope

This construction frees every coefficient entry of the reflection tensor,
subject to explicit divisibility conditions at chi=1+w. Local tensors may
have poles at w=1. Every paired full-letter-type restriction of every tensor
power is nevertheless polynomial there and has the crossing endpoint.
The construction includes arbitrary non-Gaussian coefficient changes.

At four physical symbols, one family has 6,144 independent F4 coefficients,
numerator degree 25, and effective degree 24 after the exact typed division.
A larger family has 6,237 coefficients, numerator degree 27, and the same
effective degree 24. These counts describe complete affine coefficient
spaces under the stated local conditions. They are not counts of inequivalent
tensors or a classification of all possible typed endpoints.

The endpoint and degree theorems hold for every power-of-two reflection
interval M at least 4. The finite affine solver currently implements M=4.
No pressure crossing or improved factoring bound is claimed.

## 1. A family of admissible charges

Put m=M/2. Choose rational mode weights a_i,b_i satisfying

\[
 A=\max_i|a_i|,\qquad B=\max_i|b_i|,\qquad A+B\le2.
 \tag{1}
\]

For a complete physical occupation entry lambda=(x,y,u,v), define

\[
 E(\lambda)=\sum_i a_i(x_i-v_i)+\sum_i b_i(y_i-u_i).
 \tag{2}
\]

Paired full-letter word types force the total E to be zero. The X,V type
can differ from the Y,U type. The bookkeeping weights in (2) are separate
from the physical heights used in the zero and infinity cost matrices.

Let T be the original polynomial reflection tensor and v the chi-adic
valuation. Its Gaussian kernel gives

\[
 v(T_\lambda)+E(\lambda)\ge0
 \quad\text{for every nonzero entry.}
 \tag{3}
\]

Indeed, a constant matched edge has charge zero and valuation zero. Every
external nonconstant edge has valuation two and charge at least -(A+B).
An internal reversal edge has valuation 2M and charge at least -2A or -2B.
The latter magnitudes are at most four, while 2M is at least eight.
Thus every edge satisfies (3). Charges and valuations add in each Gaussian
matching term, and cancellation can only increase the valuation of the sum.
The common scalar is a power of w and is a chi-adic unit.

The indexed charge used in
[the quartic construction](reflection-charged-quartic-tiles.md) is recovered
by a_i=b_i=i/(m-1). Other useful choices are

\[
 a_i=b_i=1,
 \qquad\text{or}\qquad a_i=2,\ b_i=0.
 \tag{4}
\]

Independent signs on these weights are permitted when (1) holds. For the
first charge in (4), a creation on X or Y and an annihilation on U or V
each increase E by one. A one-per-leg CC-AA quartet therefore admits poles
of order at most three. For the second charge, the quartet
C_X(i) C_X(j) A_V(i) A_V(j) has increment eight and admits poles of order
at most seven. These statements follow from the same strict charge-budget
proof as the indexed quartic family.

## 2. The typed endpoint theorem

Perturb each entry independently by a rational function delta T_lambda
whose only possible finite pole is at chi=0. Require either

\[
 v(\delta T_\lambda)+E(\lambda)>0,                         \tag{5}
\]

or the following larger, explicitly stated condition:

\[
 v(\delta T_\lambda)+E(\lambda)\ge0,
 \quad\text{with equality allowed only when }E(\lambda)<0.
 \tag{6}
\]

Zero perturbation entries satisfy either condition by convention.

**Theorem.** Under (3) and either (5) or (6), every paired full-letter-type
restriction of every tensor power of T+delta T is regular at w=1 and has
the crossing endpoint. If all coefficients are Laurent polynomials in w
and chi^{-1}, the retained coefficients are Laurent polynomials in w alone.

For (5), expand a retained coefficient into products of baseline entries
and perturbation entries. The total charge is zero. Every factor has
nonnegative valuation plus charge, and a nonempty perturbation selection
has at least one strictly positive factor. Its total valuation is therefore
positive, so every pole cancels and the term vanishes at w=1.

For (6), suppose a nonempty perturbation selection had total valuation
zero. Every factor would then have valuation plus charge zero. Any baseline
factor with this property has E at most zero: a noncrossing baseline entry
vanishes at w=1 and hence has positive valuation, while a matched entry has
E=0. Every selected perturbation factor with equality has E strictly below
zero by (6). The total charge would be negative, contradicting the type
restriction. Thus nonempty perturbation selections again have positive
valuation. The empty selection specializes to the original crossing.

For integer charges, the larger condition has a direct polynomial
interpretation. The entrywise gauge G_lambda=chi^{E(lambda)}T_lambda is
polynomial, and G(1) is crossing plus entries of strictly negative charge.
One may freely change those negative-charge endpoint entries while keeping
the zero-charge crossing entries fixed. Paired word types kill every
nonempty negative-charge endpoint term. The chi gauge cancels identically
on each retained word. The valuation proof above also covers rational
charges without requiring fractional powers of chi.

Equality at positive charge is not allowed by this argument. At M=4 with
the count charge, take lambda=(1,1,0,0) and mu=(0,0,1,1), using the first
logical occupation bit. Their charges are two and minus two, and T_mu has
chi-order two. Adding chi^{-2} to T_lambda creates a nonzero endpoint in
their two-copy paired-type product, although that word is not a crossing
word. The exact checker retains this negative control.

## 3. The complete affine coefficient space

For the strict condition (5), set

\[
 k_\lambda=\lfloor-E(\lambda)\rfloor+1.
 \tag{7}
\]

For the larger condition (6), replace k_lambda by -E(lambda) exactly when
E(lambda) is a negative integer. In all other cases keep (7). Choose

\[
 S=\max(0,-\min_\lambda k_\lambda),
 \qquad D\ge D_{\rm base}+S,
 \qquad D_{\rm base}=3M^2/2.
 \tag{8}
\]

Every allowed degree-D common numerator in this family has the form

\[
 N_\lambda
 =\chi^S T_\lambda+\chi^{S+k_\lambda}P_\lambda(w),
 \qquad
 \deg P_\lambda\le D-S-k_\lambda.                          \tag{9}
\]

The P_lambda are independent polynomials over a characteristic-two field.
A negative degree cap means P_lambda=0. Conversely, every entry satisfying
the stated congruence and degree bound has this form, by exact polynomial
division. Thus (9) is the complete affine space for these conditions.

The maximum possible charge is sum |a_i|+sum |b_i|, at most 2m. Therefore
S is at most 2m-1. With D=D_base+S+e, the exact number of coefficients is

\[
 \sum_\lambda\max(0,D_{\rm base}+e-k_\lambda+1).          \tag{10}
\]

For integer weights and the strict condition, this becomes
q^4(D_base+e), where q=2^m. Indeed, each entry has D_base+e+E(lambda)
coefficients and the sum of E over the complete alphabet is zero.

At M=4 and e=0 the counts are:

| Charge weights (a_0,a_1;b_0,b_1) | Boundary condition | S | Numerator degree | Coefficients |
|---|---|---:|---:|---:|
| (0,1;0,1) | Strict | 1 | 25 | 6,144 |
| (1,1;1,1) | Strict | 3 | 27 | 6,144 |
| (2,2;0,0) | Strict | 3 | 27 | 6,144 |
| (1,1;1,1) | Negative boundary allowed | 3 | 27 | 6,237 |
| (2,2;0,0) | Negative boundary allowed | 3 | 27 | 6,224 |

The negative-boundary extension adds one coefficient at each negative-charge
entry. For count charge there are 93 such entries among the 256 entries.
All rows in the table have effective degree 24.

## 4. Degree, rank, and two-place costs

The retained numerator of an n-fold tensor power of (9) is divisible by
chi^{Sn}. After that exact division it is polynomial of degree at most
n(D-S). This is the effective degree used in the interpolation construction.
The implementation takes polynomial P_lambda and charges every numerator
degree. If Laurent P_lambda are used in a separate construction, their
common w-clearing must also be included; a maximum exponent alone is not
a support-width bound.

Every local scalar cut has full generic rank by the
[paired-type endpoint rank theorem](paired-type-endpoint-full-rank.md).
This remains true after selecting any fixed paired physical subalphabets,
including unequal alphabet sizes. Local regularity at w=1 is unnecessary;
the all-length regular type restrictions supply the rank argument.

Since the local tensor is N/chi^S, its original weighted Smith orders are

\[
 \sigma_0=\sigma_0(N),\qquad
 \sigma_\infty=\sigma_\infty(N)+S.                         \tag{11}
\]

The infinity shift is included before evaluating either cost partition.
The physical word types and their full symbol entropies remain part of the
criterion. A relative gain at fixed heights, a small absolute pressure, or
an omitted common denominator does not certify a factoring improvement.

For O(M^2) effective degree, the complete local table and its coefficients
have size 2^O(M). All coefficient operations and the exact retained division
remain within that construction scale. The optional negative boundary and
the common pole of order O(M) are explicitly included in these bounds.

## 5. Exact two-cut coefficient synthesis

[`typed_charge_lattice_synthesis.py`](../tools/typed_charge_lattice_synthesis.py)
constructs the complete M=4 affine space over F4. Fixed regular left and
right Smith frames turn desired row orders in both cuts into binary affine
linear equations on the low coefficient plane and homogeneous equations
on the high plane. Every prescribed chi divisibility is imposed before
these row constraints. The solver can reuse the exact common-fiber
compression of the [strict affine interface](full-alphabet-two-cut-affine-synthesis.md).

The default frames come from the Gaussian source. A binary numerator from
the same affine family can instead supply new frames for coefficient/frame
alternation. Those frames may be computed to a finite precision that covers
every requested row target; their polynomial representatives are regular
units, and the required product identity is exact through that precision.
The coefficient space is unchanged by this choice of frames.

Row divisibility gives sufficient lower bounds on all Smith prefixes.
It does not describe every tensor with those prefixes. Retained candidates
are evaluated with their actual complete profiles in both cuts. A failed
walk in a fixed frame or degree window is not an obstruction theorem.

## Verification

[`check_typed_charge_lattice_synthesis.py`](../tools/check_typed_charge_lattice_synthesis.py)
checks the generalized Gaussian edge budget, original coefficient charges,
both affine coefficient tables against literal polynomial matrix products,
entrywise chi divisibility, paired two-copy endpoint cancellation, the
complete restored profiles, and independent rational DVR profiles. It also
tests the positive-boundary failure and rejects an invalid charge budget.

```sh
python tools/check_typed_charge_lattice_synthesis.py \
  --output evidence/typed-charge-lattice-synthesis.json
```

The theorem and exact interface are distinct from the exploratory target
walks. A score change still requires a strict full-pressure witness and its
complete recursive cost accounting.
