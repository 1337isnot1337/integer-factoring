# Retained coefficient precision and Smith-profile reuse

## Result

The Fock frontends now distinguish the coefficient precision used during
construction from the precision retained in their exported arrays.
`builder_known_below` records the former. `known_below` is a conservative
common absolute cutoff for the arrays that a caller actually receives.

This is a metadata hardening change. No incorrect profile was reproduced
in the existing standard path, which uses the same window for export and
Smith elimination. That path also has a native weighted window cap.
The distinction matters when a caller aligns, combines, pads, or otherwise
reuses the exported arrays with a larger elimination window.

## 1. The two cutoffs

Suppose construction determines every Laurent coefficient below absolute
order `K`. An exported nonzero entry has valuation `v_ij`; the frontend
stores only `W` coefficients after this valuation. The exported entry is
therefore known below

\[
\min(K,v_{ij}+W).
\]

Using `v_min=min v_ij` gives the safe common bound

\[
\boxed{K_{\rm retained}=\min(K,v_{\min}+W).}
\tag{1}
\]

An entry that vanishes throughout construction is known below `K` and
does not weaken (1). If the tracked degree bound and the complete output
degree prove that every polynomial coefficient is present, the existing
`exact_polynomial_matrix` flag remains sufficient and no truncation guard
is required.

For positive scale `ell` and row and column weights `r_i,c_j`, a common
weighted cutoff is

\[
K_{\rm retained}\ell+\min_i r_i+\min_j c_j.
\tag{2}
\]

Returned Smith orders must stay strictly below this bound when the
matrix is not known to be exact. Increasing the elimination window
does not recover coefficients that were discarded during export.

## 2. An adversarial example

Over the binary polynomial ring, consider

\[
A(z)=
\begin{pmatrix}
1&1&0&0\\
1&1+z^{150}&z^{200}&0\\
0&z^{200}&z^{300}&0\\
0&0&0&1
\end{pmatrix}.
\tag{3}
\]

Its exact Smith orders are

\[
(0,0,150,250).
\tag{4}
\]

After the first two unit pivots, the relevant two-by-two block has
smallest entry order 150 and determinant of order 400, giving (4).

Construct (3) with 1024 coefficients, then export only 128 coefficients
after each entry's valuation. The entry `1+z^150` becomes `1`; the
monomials at orders 200 and 300 retain their separate Laurent shifts.
If the exported arrays are padded and interpreted as exact polynomials,
their Smith orders are

\[
(0,0,200,200).
\tag{5}
\]

The construction cutoff 1024 would incorrectly certify (5) if used as
a common retained cutoff. Formula (1) gives 128 and rejects it.
A fresh export at window 512 retains the missing term and recovers (4).

With the original 128-coefficient Smith window, the native engine rejects
this example because it cannot resolve every promised pivot. Its cap is
the smallest weighted input-entry valuation plus 128 times the scale.
Thus the example demonstrates unsafe reuse with an enlarged window,
not an incorrect return from the existing same-window path.

## 3. Implementation and replay

The shared `retained_known_below` function is used by both
[rectangular_defect_profiles.py](../tools/rectangular_defect_profiles.py)
and
[rectangular_defect_coordinates.py](../tools/rectangular_defect_coordinates.py).
The original construction cutoff is retained separately for diagnosis.

Run

```sh
python3 tools/check_retained_profile_precision.py \
  --output evidence/retained-profile-precision.json
```

The checker constructs (3) through the actual Fock program backend. It
checks the same-window native rejection, the wrong padded-array profile,
the corrected metadata rejection, a fresh exact replay, and a weighted
scaled control. The [evidence](../evidence/retained-profile-precision.json)
records each result. The original exact-matrix path is unchanged.
