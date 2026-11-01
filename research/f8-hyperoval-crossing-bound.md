# A ten-point F8 code improves the 96-letter crossing

## Result and scope

An explicit 96-letter coherent crossing over F8 has the exact strict bound

```math
\rho < \frac{142001167}{200000000}=0.710005835.
```

This strengthens the [retained F4 bound](f4-hyperoval-crossing-bound.md)
`0.710005861` by `0.000000026` without increasing the alphabet size.
Moreover, **both natural channel profiles improve at every positive
scale**, not just at one numerical optimum.

This remains above `1/sqrt(2)` and does **not** improve the complete
factoring exponent `1/2`. It is a stronger bound for this finite
96-letter construction, not a global best over arbitrary alphabet sizes.
Larger full Gaussian alphabets already approach the threshold from above.
No new Lean formalization or independent mathematical acceptance is claimed.

## 1. The physical maps

Use twelve hot modes and `h(w)=sum(2i+1)` over occupied positions.
Retain all 92 coordinate words of height below 25, and the five-particle
word 31 of height 25. These give 93 labels.

The remaining three labels use all ten three-particle words of height 25:

```text
52,76,82,97,138,145,262,265,517,1027.
```

Work over `F8=F2[x]/(x^3+x+1)`, with elements encoded by integers zero
through seven. Use the ten columns

```math
(1,t,t^2)^{\mathsf T}\quad(t\in\mathrm{F8}),
\qquad (0,1,0)^{\mathsf T},\quad(0,0,1)^{\mathsf T},
```

in the displayed field-element order, followed by the two final columns.
The rows of this three-by-ten matrix give the three embedding vectors.
Retraction reads columns 0, 8 and 9, which form the identity matrix.
Thus `LP=I` exactly.

All 120 three-column minors are nonzero: the finite-column case is
Vandermonde, and the cases with one or both final columns reduce to
nonzero differences, their squares, or one. This is the standard even-field
hyperoval construction; the code itself is not claimed to be new.
Every deletion of a boundary word lies in the retained coordinate core,
so both physical map spaces are deletion closed and homogeneous.
The [F8 coherent construction](coherent-hot-f8-retractions.md) gives the
strict paired endpoint and both natural interfaces.

## 2. Complete profiles and independent entry controls

The [full-profile table](../evidence/f8-hyperoval-q96-full-profile.json)
contains 9,216 orders in each cut. The
[clean construction replay](../evidence/f8-hyperoval-q96-clean-replay.json)
reconstructs both lists with cutoff 127 and requires **zero unresolved
channels**. The highest orders are 116 and 108.
Historical prototype metadata in the profile table is provenance only;
the operative source validation is the clean replay and its current hashes.

The separate [entry controls](../evidence/f8-hyperoval-q96-entry-controls.json)
compare 76 complete polynomial entries against determinants of the
original binary one-particle matrix, using 164 independent minors.
Fourteen entries have a nonzero third F8 plane. All eleven occupation
levels and both cuts are covered. This is not a second reconstruction
of every matrix entry.

Run the source checks and selected-entry controls with:

```sh
python3 tools/replay_f8_hyperoval_q96.py
```

Rebuild both full profiles with `--full`. The
[F8 engine controls](coherent-hot-f8-retractions.md) additionally compare
70 complete block profiles with a separate binary-restriction method.

## 3. Exact pressure and a rational physical type

Let `A=142001167`, `B=200000000` and
`z=557123173/557123174`. From the two complete order lists and physical
heights, form

```math
Z_0=\sum_j z^{A\nu_{0,j}},\qquad
Z_\infty=\sum_j z^{A\nu_{\infty,j}},\qquad
Q=\sum_i z^{Bh_i}.
```

The [arithmetic checker](../tools/certify_f8_hyperoval_q96.py) proves

```math
4A\log Q-B(\log Z_0+\log Z_\infty)>0.
```

At 192 bits and 96 logarithm terms, the positive lower numerator is

```text
3447922615016541404695491904464166232117665971367324507388
```

in units of `2^-192`. A separate 224-bit, 112-term calculation passes.
A nearby test of `0.710005834`, at an explicitly adjusted rational base,
has strictly negative gap. That negative control is not a global lower
bound over other scales or constructions.

The [pressure receipt](../evidence/f8-hyperoval-q96-pressure.json) also
specifies a rational physical type with denominator `2^36`. Its smallest
count is 2,813,418, and its own outward entropy inequality is strictly
positive. Thus the construction uses an explicit rational type, not only
an optimized real-valued distribution. The fixed-seed uniform word and
cost construction in the F4 note carries over; F8 arithmetic adds only
a constant cost factor.

## 4. Strict improvement at all natural scales

The two constructions have exactly the same physical height inventory.
For each cut, let `mu` be the old sorted order list and `nu` the new one.
The checker verifies all 9,216 exact prefix inequalities

```math
\sum_{i=1}^k\nu_i\ge\sum_{i=1}^k\mu_i.
```

Both first become strict at prefix 3,320. Each cut has 5,897 strict
prefixes. The total order gains are 528 at zero and 310 at infinity.

For every real `u`, the hinge sum is

```math
\sum_i(u-\nu_i)^+=\max_k\left(ku-\sum_{i=1}^k\nu_i\right).
```

Thus the new hinge sum is never larger. Integrating against
`a^2 exp(-a u)` for `a>0` proves that the new channel partition sum is
smaller at every positive natural scale. Strictness follows because the
lists differ. This comparison does not extend to arbitrary physical
weights or asymmetric criteria.

The broader [protected-core barrier](height36-protected-core-floor.md)
still applies. A factoring-exponent improvement needs an escape from
its stated hypotheses, rather than further tuning only this boundary.
