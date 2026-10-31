# Polynomial row flags with exact common-floor reuse

## Result and scope

The [equal-order row-flag construction](f4-equal-order-smith-flags.md)
extends to polynomial row transvections across unequal reference orders.
The exponent in each transvection preserves the common-floor module and
ensures that every coefficient needed by a later row condition is already
known. The original floor system is eliminated once.

The implementation is
[`f4_polynomial_smith_flags.py`](../tools/f4_polynomial_smith_flags.py), with
the driver
[`search_f4_polynomial_smith_flags.py`](../tools/search_f4_polynomial_smith_flags.py).
It changes the regular row frames used by a finite coefficient search.
It does not establish a new factoring bound.

## 1. A transvection preserving all three row requirements

At one place, let \(\sigma_i\) be the reference row orders,
\(f_i\) the common floors, and \(c_i\) the exclusive coefficient caps.
Consider the elementary row operation

\[
F_i\longleftarrow F_i+\lambda z^eF_j,\qquad i\ne j,
\quad\lambda\in\mathbf F_4^\times.
\tag{1}
\]

Require

\[
e\ge\max(0,\sigma_i-\sigma_j,f_i-f_j,c_i-c_j).
\tag{2}
\]

The elementary matrix in (1) is a polynomial unit. In characteristic two
it is its own inverse. Condition (2) implies three facts.

First, a reference row of order at least \(\sigma_j\) contributes to
row \(i\) only at order at least \(\sigma_i\). Thus every original
reference row target remains feasible after the operation.

Second, (1) preserves the row-floor module
\(\mathrm{diag}(z^{f_i})\mathrm{Mat}_n(\mathbf F_4[[z]])\).
Its inverse has the same property, so the common-floor tensor fiber is
unchanged, not merely included in a larger fiber.

Third, for any requested coefficient power \(k<c_i\), the source power
needed from row \(j\) satisfies

\[
k-e<c_j.
\]

If that source power is below \(f_j\), it vanishes on the entire common
floor fiber. Otherwise it is among the explicitly stored coefficients.
No unknown coefficient above a cap is treated as zero.

For the usual windows

\[
f_i=\max(0,\sigma_i-a),\qquad c_i=\sigma_i+b,
\]

the condition \(e\ge\max(0,\sigma_i-\sigma_j)\) already implies (2).
The search can add a further nonnegative exponent. Every operation and
its coefficient are retained in the evidence.

## 2. Applying the operation in free coefficient coordinates

The column-resolved affine expressions from the equal-order construction
are enough. For each stored coefficient of recipient row \(i\), the
implementation adds the scalar multiple of source row \(j\) at power
\(k-e\). It performs this directly in the already compressed free
coordinates. Products of transvections therefore reuse the same affine
parameterization throughout.

In parallel, the implementation multiplies the explicit row-frame change
\(U(z)\). It keeps its polynomial jet modulo \(z^P\), where
\(P=\max_i c_i>\max_i\sigma_i\). This does not change any coefficient
used by the search. Its constant determinant remains nonzero, so the
displayed polynomial representative is still a regular unit.

The same statement holds after conjugation by the floor or reference
diagonal. Indeed, the omitted error has conjugated order at least
\(P-f_i+f_j>0\), or at least \(P-\sigma_i+\sigma_j>0\).
Thus the constant terms of the conjugated matrices remain invertible.
The independent replay verifies these unit conditions by binary
restriction of scalars.

## 3. Exact controls and a completed finite pilot

Run the construction controls with:

```sh
python3 tools/check_f4_polynomial_smith_flags.py \
  --output /tmp/f4-polynomial-smith-flag-controls.json
```

The [record](../evidence/f4-polynomial-smith-flag-controls.json) contains
fifteen polynomial flags, 340 row transvections, 6,360 literal coefficient
comparisons, 84 feasible tensor substitutions, three identity controls,
and nine rejected invalid operations. Literal comparisons multiply the
full polynomial matrices and changed frames without truncation.

A [five-symbol pilot](../evidence/f4-polynomial-smith-flags-five-symbol-prefix.json)
completed three flags and 120 target walks. The first flag is the unchanged
frame. The other two contain positive-degree, nonbinary transvections.
Each produced one target vector and one complete-profile class. Its nine
profile queries left the reference ratio at approximately
\(0.8335422390\), above the required threshold.

The pilot uses the degree-18 reflection restriction, the paired energies
\((0,12,6,0,8)\), row floors four below the reference orders, and caps
five above them. Its exact common-floor fiber has 379 free F4 coordinates.
This is a finite search at one physical gauge and one common-floor window;
it does not exclude other row flags, coefficient fibers, or seeds.
