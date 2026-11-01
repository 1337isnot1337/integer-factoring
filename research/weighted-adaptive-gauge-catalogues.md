# Weighted adaptive flags with a catalogue of local gauges

## Result and scope

The adaptive opposite-flag construction extends to fixed physical types,
rational local gauges, and a separate choice of gauge for every jet order.
The actual polynomial tensor, physical basis, and paired types remain
fixed. A gauge is used only to bound the rank of a specified original jet.
Consequently different gauges can certify different jets without changing
the reconstruction functional or its endpoint.

The exact coefficient bound below also removes a loose clearing-power
factor from earlier finite checks. For the existing ternary basis witness
at tensor length three, the finite row jet at place two has certified rank
bound one, equal to its actual rank; its earlier recorded bound was seven.
This is a finite improvement, not a new asymptotic crossing parameter.

The implementation supplies exact dynamic programming and complete small
typed reconstructions. No improved factoring score is claimed.

## 1. Fix one tensor and one pair of physical types

Let `K(z)` be a degree-`D` polynomial crossing on an alphabet of size `q`.
Any constant physical basis change is applied once to this entire tensor.
Fix rational types `pi_X` on the paired `x,v` legs and `pi_Y` on the paired
`y,u` legs. At every admissible tensor length `m`, all type counts are
integers. Both cost cuts have the same exact typed dimension

\[
 d_m=\binom{m}{m\pi_X}\binom{m}{m\pi_Y}.
 \tag{1}
\]

Its exponential rate is `H=H(pi_X)+H(pi_Y)`.

Choose disjoint row and column place sets, avoiding the crossing endpoint.
All selected places are rational over one fixed coefficient extension,
as in the mixed-place theorem.
At each place use the ordinary local parameter, with the degree-`Dm`
homogeneous reversal at infinity. Every reconstruction functional acts on
the original polynomial space of degree at most `Dm`.

For one local gauge, choose an integer ramification `L>=1` that makes
all four physical-letter heights integral. If the original parameter is
`t`, put `t=tau^L`. Denote the resulting integral letter heights by
`h_x,h_y,h_u,h_v`, and their common paired-type mean by

\[
 A=\sum_a\pi_X(a)(h_x(a)+h_v(a))
   +\sum_a\pi_Y(a)(h_y(a)+h_u(a)).
 \tag{2}
\]

Let `nu_1,...,nu_(q^2)` be the complete local Smith orders of the ramified,
gauged seed in the chosen cut. They may be negative. The input types make
`mA` an integer, since every retained entry has total gauge height `mA`.

After restriction to those types, the original order-`j` Hasse jet is
therefore the coefficient of order

\[
 Q=Lj+mA
 \tag{3}
\]

in the gauged tensor power. The ramification changes this analysis of the
jet. It does not multiply the global polynomial degree budget.

## 2. An exact bound for one jet

For a channel word `w=(w_1,...,w_m)`, put
`nu(w)=nu_(w_1)+...+nu_(w_m)`. The rank of the specified typed jet is at
most

\[
 c_g(j)=\sum_{w\in\{1,\ldots,q^2\}^m}
              \max\{Lj+mA-\nu(w)+1,0\},
 \tag{4}
\]

and at most `d_m`.

To prove (4), write the local Smith decomposition as analytic invertible
left and right matrices around a diagonal matrix of powers `tau^nu_i`.
Tensoring gives one diagonal channel for every word `w`. Its coefficient
of order `Q` is a sum of at most `Q-nu(w)+1` outer products when that
number is positive, and is zero otherwise. Restricting rows and columns
to the physical type classes cannot increase rank. This proof works
directly with negative Smith orders; no extra scalar clearing factor is
needed in the final bound.

The formula is invariant under translating all orders and the mean shift
by the same constant. Omitting the mean shift would violate (3) and can
give a false zero-rank conclusion.

For fixed seed and gauge, (4) is computable in polynomial time in `m`.
Subtract the minimum channel order, convolve the resulting nonnegative
order multiplicities `m` times, take cumulative channel counts, and then
take one more cumulative sum. Evaluating the result at

\[
 Lj+m(A-\min_i\nu_i)
\]

gives (4). The support range and the bit length of every count are
polynomial in `m`. The implementation uses the independently checked
integer-packing convolution from the adaptive flag code.

## 3. A different gauge may certify each jet

At a fixed place, let `G` be a finite catalogue of valid gauges, each
with its own ramification, complete orders, and mean shift. All these
gauges refer to the same physical basis and the same pair of types.
Define

\[
 c(j)=\min_{g\in G}c_g(j).
 \tag{5}
\]

Every function `c_g(j)` is nondecreasing, by (4). Their pointwise minimum
is also nondecreasing. Therefore the same greedy merge used in the
[adaptive opposite-flag theorem](adaptive-mixed-place-crossings.md) applies
to the sequences (5): repeatedly select the least-cost next available jet
among the places on a given side. Only the next jet of each place is
available, so every prefix remains an initial segment at every place.

The selected gauge is recorded with each jet. It certifies that jet's
rank; it does not alter the jet vector in the polynomial dual space.
The row and column flags still have complementary intersections because
their place supports are disjoint. Thus their intersections give the
same common reconstruction pieces as before.

If `r_k` is the sum of the first `k+1` row-jet costs and `s_l` is the sum
of the first `l+1` column-jet costs, set

\[
 R_k=\min\{d_m,r_k\},\qquad C_l=\min\{d_m,s_l\}.
 \tag{6}
\]

The piece indexed by `k=0,...,Dm` has rank product at most

\[
 R_k C_{Dm-k}.
 \tag{7}
\]

The type restriction and the physical basis are common to every piece.
Choosing a different physical basis or a different type for each jet
would not be justified by this argument.

The greedy prefix minimizes the sum of the available costs among all
admissible initial-segment allocations of the same size. Hence adding a
valid gauge to a catalogue can only improve the resulting exact bound.
Zero-cost jets, negative local orders, and changing optimal gauges are
all handled by the same finite algorithm.

## 4. Asymptotic use of the catalogue

For a fixed gauge, a Chernoff bound on (4) gives the local exponential
rate upper bound

\[
 \Phi_g(d)=\inf_{x\ge0}
 \left\{\log\sum_i e^{-x\nu_i}+x(Ld+A)\right\}
 \tag{8}
\]

at jet order `j=dm+O(1)`, with a polynomial factor in `m`. The exact type
dimension supplies the additional rate bound `H`. At each place the
catalogue therefore gives

\[
 \Phi(d)=\min\left\{H,\min_{g\in G}\Phi_g(d)\right\}.
 \tag{9}
\]

This is the lower envelope of the available affine Chernoff bounds in
`d`; it is nondecreasing and concave wherever finite. The number of
cheap jets available below a proposed exponential cost is the sum of
the capacities supplied by the separate places. Inverting those summed
capacities gives the adaptive prefix rates, as in the unweighted
construction. At the end of a zero-cost prefix, the inverse must use
the right limit: one additional jet is needed for the flag dimension.
The exact algorithm (4)--(7) treats that boundary without an interpolation
convention.

Consequently a finite collection of rational gauges and rational
supporting bounds can certify different portions of the jet range.
Their comparison allocations need not be nested, because the one
actual greedy flag is no more expensive than each admissible comparison
prefix. This is a way to turn numerical local optimizers into a
reviewable adaptive certificate while retaining one common typed tensor.

No asymptotic improvement is inferred merely from a finite run of the
dynamic program. A new uniform `rho` claim still needs a bound on (7)
for every sufficiently large admissible `m`, or an equivalent exact
interval certificate with strict margin, followed by the existing
paired-type transfer.
That final all-width transfer requires positive type entropy `H>0`.
The finite rank bounds themselves also apply to zero-entropy types.

## 5. Exact controls

```sh
python3 tools/check_weighted_adaptive_catalogues.py \
  --output evidence/weighted-adaptive-gauge-catalogues.json
```

The checker uses the already certified ternary inverse-Pascal basis.
It computes 12 complete ramified Smith profiles from independent block
Toeplitz ranks. The three locations each have four gauges, including
zero heights and negative-order profiles.

Two common paired-type cases are then built in full: tensor length two
with unequal types, and tensor length three with uniform types. The
checks include 36 exact typed jet ranks, 24 ranks of actual opposite-flag
pieces, and 1,312 endpoint entries. Every piece is reconstructed from the
actual intersection functional, and the pieces sum to the original
typed crossing endpoint.

Independent word enumeration verifies 208 coefficient bounds and their
208 translated-gauge counterparts. An abstract cost example checks all
nine greedy prefix minima and changes the optimal gauge along a place.
The missing-mean negative control falsely predicts rank zero at a jet
whose independently computed rank is one, and is rejected.

The code and evidence are
[`weighted_adaptive_mixed_place.py`](../tools/weighted_adaptive_mixed_place.py),
[`check_weighted_adaptive_catalogues.py`](../tools/check_weighted_adaptive_catalogues.py),
and [`weighted-adaptive-gauge-catalogues.json`](../evidence/weighted-adaptive-gauge-catalogues.json).
