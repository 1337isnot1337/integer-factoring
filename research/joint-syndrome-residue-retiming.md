# Arbitrary coefficient fields and ramification reduce to binary templates and convex offsets

## The reduction

Fix a binary syndrome quotient of rank `c`, with `B=2^c` coefficient
polynomials. Every valid joint polynomial filter is weakly dominated, in
both cost factors and for every physical measure, by a sum of at most `B`
binary polynomial vectors with independent colors and fractional exponent
offsets. The colors can lie in a field of degree at most `B` over `F_2`.

For each fixed collection of binary vectors and each fixed charge, scale,
and physical mode assignment, optimizing all offsets is a compact convex
problem. Offsets may change order or coincide. With fixed
integer normalized degree caps, there are finitely many binary templates,
and the resulting finite collection of convex problems has exactly the
same infimum as the original family with arbitrary coefficient fields and
arbitrary allowed ramification.

The finite collection can be very large. This is an exact structural
reduction and a constructive search interface, not an exhaustive search
of all those templates or a factoring improvement.

## 1. The binary polynomial transform

Use the notation of
[joint syndrome synthesis](joint-syndrome-polynomial-filters.md). Let `m=nP`,
`X=x^m`, and write the coefficient vector as `f=(f_b)`, indexed by the `B`
syndromes. Its group-state transform is

\[
 \delta_z=\sum_b X^cW_{z,b}(X^{-1})f_b.
 \tag{1}
\]

Every matrix entry in (1) lies in `F_2[X]`. At `X=1`,

\[
 \delta_z(1)=f_{\sigma(z)}(1),
 \tag{2}
\]

where `sigma(z)` is the syndrome of the group parity state. All coefficient
endpoints are nonzero. With normalized zero orders `b_z` and infinity
degrees `d_b`, the full relative log cost is

\[
 L=\log\sum_z q_z e^{-2\tau(b_z-c)}
      +\log\sum_b p_b e^{2\tau d_b}.
 \tag{3}
\]

The physical and charge-shifted probabilities `p` and `q` remain separate.
The result applies to arbitrary nonnegative normalized measures in (3).

## 2. Residue and coefficient-plane decomposition

Take an `F_2` basis `xi_1,...,xi_h` for the span of the finitely many original
coefficients. The coefficient field itself need not be finite. Euclidean
division of every exponent by `m` gives

\[
 f_b(x)=\sum_{0\le r<m}\sum_{\ell=1}^h
      \xi_\ell x^r F_{r,\ell,b}(X),
 \qquad F_{r,\ell,b}\in\mathbf F_2[X].
 \tag{4}
\]

Apply (1) to each binary vector separately, obtaining `Delta_(r,ell)`.
Different residues in (4) have disjoint exponent supports modulo `m`.
At a fixed residue, linear independence of the `xi_ell` prevents
cancellation between coefficient planes. Therefore, with a zero
polynomial omitted from the corresponding minimum or maximum,

\[
 \frac{\operatorname{ord}\delta_z}{m}
   =\min_{r,\ell}\left(\operatorname{ord}\Delta_{r,\ell,z}+\frac rm\right),
 \qquad
 \frac{\deg f_b}{m}
   =\max_{r,\ell}\left(\deg F_{r,\ell,b}+\frac rm\right).
 \tag{5}
\]

This independence argument is needed before assigning new colors. It
shows that recoloring does not silently remove a cancellation that was
present in the original zero or infinity profile.

## 3. Keep only a small endpoint cover

The endpoint support of a binary component `F_j` is

\[
 E_j=\{b:F_{j,b}(1)=1\}.
 \tag{6}
\]

Since every original `f_b(1)` is nonzero, these supports cover all `B`
syndromes. Choose one component covering each still-uncovered syndrome.
After removing duplicates, at most `B` components remain.

Give the selected components `eta_1,...,eta_l` independent colors over
`F_2`, with `l<=B`, and retain their original offsets. One may use an
`F_2` basis in `F_(2^l)`, or use `l` independent elements of `F_(2^B)`.
For every coefficient endpoint, at least one selected binary component
contributes 1. Independent colors make its sum nonzero. Equation (2)
gives the transformed endpoints as well.

Removing components increases each minimum in (5) and decreases each
maximum. Thus it weakly decreases both factors of (3), simultaneously for
every charge, scale, and physical mode assignment. There is no numerical
optimization in this dominance statement.

The cover need not be unique. Choosing one cover gives a dominating
filter; it need not give the best cover. An individually admissible
component, whose endpoint support is already all `B` syndromes, dominates
the entire independently colored sum by itself. More generally, every
redundant component can be deleted whenever the remaining endpoint
supports still cover all syndromes.

## 4. Independent colors permit arbitrary offset retiming

Fix selected binary component vectors `F_j` and their transforms `Delta_j`.
For rational offsets `theta_j` in `[0,1]` whose common denominator divides
an allowed `m=nP` with `P` odd, choose such an `m` and form

\[
 \widetilde f_b(x)=\sum_{j=1}^l
       \eta_j x^{m\theta_j}F_{j,b}(x^m).
 \tag{7}
\]

Independent colors prevent cancellation between components even when
offsets coincide, or when an offset 1 meets a neighboring integer exponent
of an offset-0 component. The exact normalized profiles are

\[
 b_z(\theta)=\min_j(b_{zj}+\theta_j),\qquad
 d_b(\theta)=\max_j(d_{bj}+\theta_j),
 \tag{8}
\]

where `b_zj=ord Delta_j,z` and `d_bj=deg F_j,b`. Zero component entries are
omitted. Endpoint support is independent of the offsets. Thus (7) stays
valid under arbitrary reorderings and collisions.

A fixed normalized degree cap `D_b` gives the explicit upper bounds

\[
 \theta_j\le U_j:=\min\left(1,\min_{b:F_{j,b}\ne0}(D_b-d_{bj})\right).
 \tag{9}
\]

If normalized transformed order targets `T_z` are also imposed, they give

\[
 \theta_j\ge A_j:=\max\left(0,\max_{z:\Delta_{j,z}\ne0}(T_z-b_{zj})\right).
 \tag{10}
\]

The feasible offset set is the closed box `A_j<=theta_j<=U_j`. An empty
interval rejects the template. These constraints charge the actual
infinity degrees; offset changes are not free degree changes.

Adding the same constant to every offset multiplies all coefficients by
a common monomial. The two factors of (3) change reciprocally, so their
product is invariant. This symmetry does not require imposing the
nonconvex constraint `min theta=0` during optimization.

## 5. The offset optimization is convex

Substitute (8) into (3). The negative of each minimum of affine functions
is convex, and each maximum of affine functions is convex. The
`logsumexp` function is convex and increasing in every coordinate.
Consequently `L(theta)` is convex on the entire offset box, including
across changes of the controlling component.

An explicit smooth epigraph formulation uses auxiliary variables `u_z`
and `v_b` and linear constraints

\[
 u_z\le b_{zj}+\theta_j,\qquad
 v_b\ge d_{bj}+\theta_j
 \tag{11}
\]

for every finite component entry. Minimize

\[
 \log\sum_zq_z e^{-2\tau(u_z-c)}
  +\log\sum_bp_b e^{2\tau v_b}
 \tag{12}
\]

over (9)--(11). The objective decreases as `u_z` increases and as `v_b`
decreases, so its optimum realizes the exact profiles in (8).

Real offsets define the continuous optimization problem. At any strictly
negative full-pressure point, sufficiently close rational offsets preserve
the strict sign. Common denominators with odd ramification are dense, so
this approximation produces an actual polynomial filter when no
incompatible exact equality constraint is imposed on an offset.

With integer normalized degree caps and no extra order targets, every
nonempty interval in (9) is `[0,1]` or the singleton `{0}`. Allowed rational
offsets are dense in those boxes, including their boundaries. The real
minimum is therefore exactly the infimum over arbitrary ramifications.
For arbitrary real degree/order caps, the continuous minimum always gives
a rigorous lower bound; equality with a discrete infimum additionally
requires density of the allowed rational points in the prescribed box.

## 6. A finite complete reduction at fixed normalized degree

Suppose the degree caps `D_b` are nonnegative integers. There are only

\[
 2^{\sum_b(D_b+1)}
 \tag{13}
\]

binary component vectors within these caps, including the zero vector.
Select at most `B` of them with endpoint supports covering every syndrome,
assign independent colors, and solve the corresponding convex offset box.

Every original admissible filter with these normalized degree caps is
dominated by one of these templates, by (4)--(6). Conversely, every
rational-offset template is an actual admissible filter within the caps.
Density from Section 5 identifies the two infima. This proves the stated
finite reduction over all coefficient fields, polynomial coefficients,
and allowed ramifications. It does not enumerate the enormous collection
in (13), nor does it impose an unproved degree bound on the full research
problem. No joint convexity in scale, charge, or the final exponent root
is asserted.

## 7. A transport bound for pruning offset boxes

The following bound applies to the full cube `[0,1]^l`; it also lower
bounds any smaller box. Evaluate the cost when all offsets coincide.
For each zero state let `S_z` be the components with smallest integer
order `b_zj`. For each infinity state let `T_b` be the components with
largest integer degree `d_bj`. Let `alpha_z` and `gamma_b` be their
normalized thermal weights at this common offset. Define

\[
 d=\max_{U\subseteq\{1,\ldots,l\}}
 \left(\sum_{z:S_z\subseteq U}\alpha_z
       +\sum_{b:T_b\subseteq U^c}\gamma_b-1\right).
 \tag{14}
\]

The value is nonnegative because `U` empty or full gives zero. Then

\[
 \boxed{L(\theta)\ge L(0)-2\tau d
        \quad(\theta\in[0,1]^l).}
 \tag{15}
\]

To prove this, put a bipartite edge between zero state `z` and infinity
state `b` exactly when `S_z` and `T_b` intersect. Supplies are `alpha_z`
and demands are `gamma_b`. The minimum unmatched mass is the Hall deficit
(14): for a set of zero states, use the union `U` of their component sets;
including all other `S_z` contained in `U` can only increase its deficit.
The finite max-flow/min-cut theorem therefore gives a matching of mass
`1-d`.

For each matched piece choose a common component in the intersection.
Assign the corresponding negative zero subgradient and positive infinity
subgradient to that component; they cancel. The unmatched pieces have
total mass `d` on each side. This constructs a subgradient of `L` at the
center of the cube with l1 norm at most `4 tau d`. Convexity and the
distance at most `1/2` in every coordinate prove (15). If `d=0`, coincident
offsets are a global minimizer. This can reject a complete continuous
offset family without running an optimizer.

## Implementation and verification

[`joint_syndrome_residue_retiming.py`](../tools/joint_syndrome_residue_retiming.py)
implements residue and binary-plane decomposition, endpoint-cover
selection, exact color lifts at rational offsets, and the convex offset
program. Numerical optimizer output is a proposal and is recomputed from
the actual min/max profiles before use. Its transport bound is evaluated
numerically during searches; the checker below verifies rational transport
certificates independently.

The exact controls include original-profile preservation before dropping
components, pointwise dominance after the endpoint cover, retimings with
collisions and reversed residue order, explicit nonzero endpoints, degree
cap boxes, and a negative control that fails if independent colors are
discarded. A strict factoring claim would still require a complete
pressure certificate and the construction's recursive cost accounting.

```sh
python tools/check_joint_syndrome_residue_retiming.py \
  --output evidence/joint-syndrome-residue-retiming.json
```

The exact replay passed 48 original-profile controls, 48 endpoint covers,
240 rational retimings with 2,920 literal physical-state comparisons,
240 degree-cap checks, 64 exact fractional max-flow comparisons with the
Hall deficit, and 256 rational instances of the transport inequality.
