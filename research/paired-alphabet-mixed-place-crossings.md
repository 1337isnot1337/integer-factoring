# Mixed-place pressure for unequal paired alphabets

The weighted mixed-place construction applies to a tensor with physical
dimensions `(q_x,q_y,q_y,q_x)`. The two paired types can have different
lengths. A closed irreducible factor can also contribute some of its
conjugate places to each cut, provided the selected conjugates are
disjoint. These two features are useful for explicit subalphabet
restrictions whose new determinant factors have degree greater than two.

The evaluator
[`paired_mixed_place_pressure.py`](../tools/paired_mixed_place_pressure.py)
implements both features. It retains the complete degree, local profile,
physical mean, and common paired-type costs. The formula is an extension
of the existing
[mixed-place Hermite criterion](mixed-place-hermite-crossings.md), not
an improved factoring bound by itself.

## 1. The complete cost

Let K have degree D and endpoint

\[
 K(1;x,y,u,v)=\mathbf1_{x=v}\mathbf1_{y=u},
 \quad x,v\in X,\ y,u\in Y.
\]

Use one physical type p_X on x,v and one type p_Y on y,u at every place.
Put

\[
 H=H(p_X)+H(p_Y).
\]

Each cut has N=q_x q_y generic Smith channels. A selected group j
contains k_j conjugate places with the same profile and gauge. Give
the group allocation alpha_j, so each root gets alpha_j/k_j. The
allocations sum to one separately on the two cuts. For a common
logarithmic outer scale t>0, the local logarithmic scale is

\[
 s_j=\frac{t k_j}{\alpha_j}.
\]

The physical height vector has leg lengths `(q_x,q_y,q_y,q_x)`.
Write nu_i(h_j,s_j) for the complete weighted local orders when the
unweighted valuation is multiplied by s_j and the physical weights
are h_j. The mean shift is

\[
 \mu_j=\sum_xp_X(x)(h_{j,x}(x)+h_{j,v}(x))
       +\sum_yp_Y(y)(h_{j,y}(y)+h_{j,u}(y)).
\]

Define

\[
 C_j=\log\sum_{i=1}^{N}\exp(-\nu_i(h_j,s_j))+\mu_j,
 \qquad L_R=\max_{j\in R}C_j,\quad L_C=\max_{j\in C}C_j.
\]

The resulting balanced parameter is

\[
 \boxed{\rho=\frac{D t+L_R+L_C}{2H}}.
\]

For equal alphabets and uniform types this denominator is 4 log q,
as in the existing criterion. Both types, rather than a separate type
for each place, are essential. The implementation uses one common pair
and includes its entropy derivative in the optimizer.

For a rational certificate, the ratios h_j/s_j are rational physical
gauge exponents. A fixed local ramification clears their denominators.
The local coefficient estimates and opposite Hermite flags are those
of the linked criterion. Floating optimizer parameters alone are not
such a certificate.

## 2. Reduction to equal external alphabets

The denominator can also be checked by an explicit tensor construction.
Take K together with the copy obtained by exchanging its two paired
species:

\[
 \widetilde K_{(x,x'),(y,y'),(u,u'),(v,v')}
   =K_{x,y,u,v}K_{y',x',v',u'}.
\]

Each external alphabet now has size q_x q_y. The endpoint is the literal
crossing after the fixed pair-coordinate ordering. At length m, its
external label sets are products of the two type classes, with common
size

\[
 |\mathcal T_{p_X}(m)|\,|\mathcal T_{p_Y}(m)|
    =\exp(mH+O(\log m)).
\]

This product label set need not be a single joint type on X times Y.
Every product gauge used here is additive in its two species, so its
diagonal weight is constant on the product of the two type classes.
That is the condition needed by the weighted projection and Hermite proof.

The degree doubles. On each cut the complete local orders are all
pairwise sums of the original orders; the second factor is transposed
on the second cut, which preserves its Smith profile. Its logarithmic
channel cost and type mean therefore double. The two external paired
entropies are both H. Applying the equal-alphabet criterion gives

\[
 \frac{2D t+2L_R+2L_C}{4H}
   =\frac{D t+L_R+L_C}{2H}.
\]

Thus the unequal-alphabet evaluator does not obtain an entropy benefit
from an omitted tensor copy. The usual rational-type padding and
all-width construction apply to these explicit product label sets.

## 3. Closed places and their conjugates

Let f be monic irreducible over F4 with residue degree e. Its e roots
are distinct in a fixed finite splitting extension. The
[closed-place oracle](closed-place-weighted-oracles.md) computes the
same complete profile at each conjugate, because the original
coefficients are in F4 and the physical heights are shared.

The native matrix uses restriction of scalars from F4[alpha]/f. Each
true Smith order occurs e times in that expanded matrix. Subtracting
log e from its logarithmic channel sum recovers the N-channel cost.
This does not account for the separate Hermite capacities of the
selected roots. If k of those roots are assigned to a cut, the local
scale still contains k, as displayed in Section 1. In particular,
k can be smaller than e when the conjugates are split between cuts.

The evaluator identifies a conjugate by its monic factor and an index
from zero through e-1. It rejects repeated conjugates, overlap between
the two cuts, the strict endpoint, and aliases that would count a
rational point twice. The same factor may occur on both cuts only
through disjoint conjugate indices.

For unweighted profile screening, conjugates of the same factor are
interchangeable. The driver
[`screen_restricted_closed_partitions.py`](../tools/screen_restricted_closed_partitions.py)
therefore enumerates the number assigned to the row cut, rather than
all permutations of equivalent assignments. It retains every distinct
unweighted partition. Its interpolated adaptive rates are numerical
candidate selectors.

## 4. Independent controls

The checker
[`check_paired_mixed_place_pressure.py`](../tools/check_paired_mixed_place_pressure.py)
compares the new evaluator with the established equal-alphabet formula,
checks analytic gradients for both 2-by-3 and 3-by-2 paired alphabets,
and constructs the pair-swapped tensor in Section 2 from full polynomial
coefficients. Its replay records:

- 12 equal-alphabet value and gradient identities, with zero discrepancy;
- 360 finite-difference gradient coordinates for unequal alphabets;
- 1,296 literal endpoint entries of the pair-swapped tensor;
- four complete product profiles, agreeing with pairwise order sums;
- zero discrepancy in the resulting paired pressure identity;
- five invalid-place or alias rejections.

The largest gradient discrepancy was below 4.7e-10, and the largest
profile-convolution discrepancy was below 3.6e-15. Exact field
arithmetic supplies all underlying coefficients and local matrices;
the reported pressure and derivatives are floating evaluations.

```sh
python tools/check_paired_mixed_place_pressure.py \
  --output evidence/paired-mixed-place-pressure-controls.json
```

The [control record](../evidence/paired-mixed-place-pressure-controls.json)
binds these checks to the full input and source hashes. It supplies no
claim of global optimization or a new factoring exponent.
