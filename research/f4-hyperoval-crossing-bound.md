# A slightly stronger 96-letter crossing from an F4 hyperoval

## Result and limitation

An explicit F4 physical retraction gives a complete 96-letter crossing
with the exact strict typed pressure bound

```math
\boxed{\rho<\frac{710005861}{10^9}=0.710005861.}
```

This is a small improvement over the pinned
[coordinate-code certificate](exact-natural-pressure-intervals.md),
which gives rho < 0.710005874881191. The improvement in these displayed
upper bounds is 0.000000013881191.

This does **not** improve the complete integer-factoring exponent.
The new bound remains above 1/sqrt(2), and the established complete
exponent is still 1/2. The result is a stronger explicit crossing
input, using ordinary proofs and exact finite computations, without
new Lean formalization or independent mathematical acceptance.

## 1. The explicit physical maps

Use twelve hot modes with word height

```math
h(w)=\sum_{a\in w}(2a+1).
```

There are 92 occupation words with height below 25. Retain all of them,
and retain the five-particle word 31, whose height is 25. These give
93 coordinate directions. The ten three-particle words of height 25 are

```text
52, 76, 82, 97, 138, 145, 262, 265, 517, 1027.
```

On the first six of these words, define a three-dimensional embedding
by the following six columns, in that order:

```math
(1,t,t^2)^{\mathsf T}\quad(t=0,1,\omega,1+\omega),
\qquad (0,1,0)^{\mathsf T},\quad(0,0,1)^{\mathsf T},
\qquad \omega^2=\omega+1.
```

The three embedding vectors are the rows of this 3-by-6 matrix, viewed
as vectors in the six-dimensional physical occupation space. Retract
by coordinate evaluation on columns 0, 4, 5. Those columns form the
identity, so LP=I exactly. Together with the 93 retained directions,
this gives 96 paired labels.

Every deletion of a selected height-25 word lands among the retained
lower-height words. Both map spaces are therefore deletion closed.
Each paired label has one height and occupation. The
[coherent constructor](coherent-hot-f4-retractions.md) consequently
supplies a strict crossing and both complete natural-cut interfaces.
An ambient reflection interval M=32 suffices: the maximum physical
energy is 4 times 25, strictly below 4M, and all twelve modes fit.

For completeness, all twenty three-column minors of the six-column
matrix are nonzero. Three finite columns give a Vandermonde product.
Two finite columns and the final coordinate column give t+s; using
the other coordinate column gives (t+s) squared. One finite column
and both coordinate columns give determinant one. Distinct field
elements make these expressions nonzero. This is the standard F4
hyperoval code used as a physical embedding, not a claim that the code
itself is new.

## 2. Complete profiles and source controls

The [complete construction receipt](../evidence/coherent-hot-f4-hyperoval-q96.json)
contains both sets of 9,216 Smith orders in the deflated variable t=w
squared. The two profiles differ; neither is substituted for the other.
The native construction checks 35,973,996 occupation-block endpoint
entries across the two cuts. Charge-incompatible entries vanish by
construction. The exact complete-profile run took about 1,551 seconds
in the recorded environment; this timing is not a mathematical bound.

The entry degree bound is checked before any series truncation.
The strict endpoint proves full generic rank in every block, and every
one of those pivots is resolved. Source fingerprints match the completed
receipt. No partial profile or numerical rank decision is used.

A separate [entry checker](../tools/check_coherent_hot_f4_hyperoval_entries.py)
computes determinants of the original B(t) by subset dynamic programming,
rather than using the native height-factored K(t) determinant method.
Its [receipt](../evidence/coherent-hot-f4-hyperoval-entry-controls.json)
compares 76 complete F4 polynomial entries from all eleven occupation
levels and both cuts. This uses 107 independent determinants, includes
three nonbinary entries, and checks their endpoints and degree bounds.
It is not a second construction of the full 96-letter matrices.

## 3. An exact strict pressure inequality

Let nu_0 and nu_infinity be the two complete order lists, and let h_i
be the 96 physical heights. Set

```math
A=710005861,\qquad B=10^9,\qquad
z=\frac{2785578335}{2785578336},\qquad s=A\log(1/z).
```

Define positive rational sums

```math
Z_0=\sum_j z^{A\nu_{0,j}},\qquad
Z_\infty=\sum_j z^{A\nu_{\infty,j}},\qquad
Q=\sum_i z^{Bh_i}.
```

The exact certificate proves

```math
B(\log Z_0+\log Z_\infty)<4A\log Q.                 \tag{1}
```

It uses the outward integer power and logarithm intervals from the
[prior pressure certificate](exact-natural-pressure-intervals.md).
At 192 bits and 96 logarithm terms, the lower integer numerator of
4A log Q minus B(log Z_0+log Z_infinity), in units of 2 to the power
minus 192, is

```text
20761173205562230853374803669111538315310322338358836517868
```

It is strictly positive. A separate 224-bit, 112-term replay also
passes. The same-base test with the slightly smaller proposed bound
0.710005860 fails, providing a negative control. That failure is not
a lower bound on all constructions or all scales.

To see the implication of (1), use the rational paired physical type
p_i=z^(Bh_i)/Q on both sides. Its entropy is

```math
H(p)=\log Q+\frac{s}{A/B}\sum_i p_i h_i.
```

The full two-cut typed cost is

```math
\frac{\log Z_0+\log Z_\infty+4s\sum_i p_i h_i}{4H(p)}.
```

Substitution of (1) makes this strictly less than A/B. Both complete
channel sums and the full physical means occur in this expression.
No field-element count is interpreted as an integer count.

## 4. An explicit small-denominator type

The certificate also gives 96 positive integer type counts summing to
2^36 = 68,719,476,736, with minimum count 2,813,121. The counts are
constructed using integer lower enclosures of the Gibbs probabilities,
then distributing the remaining integer count. No floating-point
rounding decision enters their verification.

For these counts k_i and L=2^36, the checker bounds directly

```math
4\frac AB\left(\log L-\frac1L\sum_i k_i\log k_i\right)
-\frac{4s}{L}\sum_i k_i h_i-\log Z_0-\log Z_\infty>0.
```

All entropy, mean and logarithm terms use outward rational intervals.
Thus the displayed strict bound holds for this explicit finite type,
not merely an unrecorded real-valued optimizer. The standard
[finite-block amplification](finite-block-crossing-amplification.md)
applies over the fixed field F4; its charged reconstruction and entry
costs are unchanged. The strict gap absorbs the polynomial factors.

## 5. Replay and next obstruction

Run `python3 tools/check_coherent_hot_f4_hyperoval_entries.py` and
`python3 tools/certify_coherent_hot_f4_hyperoval.py`. The
[pressure receipt](../evidence/coherent-hot-f4-hyperoval-pressure.json)
includes both precision replays, the negative control, all 96 type
counts, source hashes and the strict integer inequalities. These quick
replays do not claim to reconstruct the large profile. The separate
probe script reproduces that longer computation.

The new counting-polynomial product is not better than the coordinate
product at every scale. Exact rational evaluations show a gain at some
scales and a loss at others; only the displayed strict certificate is
claimed. Its first product-polynomial difference is at order 63, which
helps explain why this boundary modification gives a small gain. A
larger change to the low-order profile is still needed to reach the
post-half threshold.

Substantial OpenAI assistant support was used for the construction,
computation, proof, verification, and text.
