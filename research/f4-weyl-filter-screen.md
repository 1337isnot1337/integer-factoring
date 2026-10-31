# Exact F4 profiles for nonlinear three-direction filters

## Result and scope

The binary occupation-filter search now evaluates genuinely correlated
three-direction filters over `F4`. The exact
[Weyl reduction](identity-operator-weyl-reduction.md) replaces the
`64` by `64` active zero matrix by an `8` by `8` matrix, with explicit
auxiliary shifts. Its occupation blocks have sizes `1,3,3,1`.

A complete finite search at interval length `n=32` includes 169 exponent
tables and every one of the 2,187 normalized nonzero coefficient tables.
The 369,603 cases produce 3,850 distinct exact profile classes. None
improves the Gaussian cost on the recorded charge/scale grid or in the
recorded local refinements. No strict critical-pressure witness or new
factoring exponent is claimed. This finite observation does not exclude
arbitrary polynomial filters, other exponent tables, or other gauges.

## 1. Admissible filters and the coefficient field

Work in characteristic two, and encode
`F4=F2[omega]/(omega^2+omega+1)` by the integers `0,1,2,3` corresponding
to `0,1,omega,1+omega`. Occupation states remain binary. For three
selected identity directions, use

\[
\alpha_d(x)=c_d x^{-e_d},\qquad
d\in\{0,1\}^3,\quad c_d\in\mathbb F_4^*.
\tag{1}
\]

The state order is lexicographic, matching Python's
`product(range(2), repeat=3)`. Common scalar multiplication and a common
Laurent monomial do not change the two-cut cost, so the search fixes
`c_000=1` and `e_000=0`.

Every endpoint value in (1) is nonzero. The physical-`Y` diagonal repair
from [the endpoint theorem](repaired-endpoint-state-filters.md) restores
the original crossing coefficient by coefficient, without changing
either weighted Smith profile or polynomial degree. A common monomial
clears all Laurent powers; its opposite shifts at the two cuts cancel.
The fixed field extension has degree two.

For the recorded initial-direction search, the selected indices are
`0,1,2`, the retained `XV` indices are `29,30,31`, and the ramification
is `P=1`. The active infinity orders are exactly

\[
-63d_0-61d_1-59d_2+2e_d.
\tag{2}
\]

All other infinity modes are Gaussian. The active zero orders and their
charge blocks come from the exact reduced matrices, not from a rank
estimate or an asymptotic approximation.

## 2. Compile the full polynomial matrices once per exponent table

The reduced geometric coefficients lie in `F2`. For fixed exponents
`e_d`, precompute each state's contribution to every reduced matrix
entry, including the complete Laurent shift and every polynomial
coefficient. Multiplication by `c_d` only determines which of the two
`F4` bit planes receives that contribution. Changing the color table
therefore uses bitwise exclusive-or operations on complete polynomials.

The native Smith routine receives both complete planes. It enlarges
the working series window if a promised pivot has not yet appeared.
It rejects a source polynomial beyond the supported window instead of
discarding its high coefficients. The common denominators in the Weyl
construction have nonzero constant term, so clearing them does not
change a zero order. They are auxiliary basis changes; the original
polynomial filter and its ordinary degree remain those in (1).

After computing the reduced orders `lambda` in occupation `d`, the
auxiliary state `s` contributes

\[
\lambda-2P\sum_j j s_j,
\qquad\hbox{in charge }d-|s|.
\tag{3}
\]

Here the selected indices are the actual physical indices. Formula (3)
is at aggregate gauge `J=1`. Changing `J` adds
`P(1-J)(d-|s|)` to each order. This permits exact profiles to be reused
throughout the charge search.

The implementation is
[the F4 profile wrapper](../tools/f4_repaired_identity_filters.py),
specifically `compiled_monomial_profiles`. Its
[independent compiler checker](../tools/check_f4_weyl_filter_compiler.py)
and [exact evidence](../evidence/f4-weyl-filter-compiler.json) compare
30 complete unreduced active matrices with the compiled profiles.
The cases cover two and three directions, initial, terminal, and mixed
geometries, and ramification one and three. A separate source with
exponent 600 verifies complete high-degree packing.

The same checker expands the full profiles and physical partition
directly in 48 small cases. It agrees with both search pressure
assemblers to within `5.33e-15`; the comparison tolerance is `1e-10`.
The exact algebraic checks are separate from these floating arithmetic
controls.

## 3. The finite search and its portable record

The exponent catalog in
[the search driver](../tools/search_f4_repaired_identity_filters.py)
contains single-state shifts, functions of total occupation, and
independent pair and triple interactions. Its construction is explicit
in `candidate_exponents`; after normalization and deduplication it has
169 tables. For each table, all `3^7=2187` normalized nonzero `F4`
coefficient vectors are tested.

The
[complete compact record](../evidence/f4-repaired-three-mode-initial-n32-screen.json)
retains every distinct exact zero profile, its exponent table, a color
representative, and the number of coefficient tables giving that
profile. The multiplicities sum to 369,603. It also retains the full
charge and scale grids and every reported local refinement.

The smallest refined critical pressure is approximately
`0.0330679625084`, attained by the Gaussian profile. The smallest
sampled relative log-cost is approximately `-8.9e-16`, within floating
roundoff. Both the same-gauge-gain flag and the negative-critical-pressure
flag are false at tolerance `10^-9`.

To avoid repeated metadata, each sorted charge profile is run-length
encoded. The record is 1,558,561 bytes, and its exact catalogue digest is

```text
33e96cc8a5504772c8483417c1f602a145ec955dc7664cb0f6c8a0e2a906ba13
```

[The compaction tool](../tools/compact_f4_filter_screen.py) checks a
round trip of every exact order, representative, multiplicity, and
stored floating observation. Its `expand_rows` function also reads the
ordinary uncompressed record format. The digest covers the exact
catalogue data, independently of elapsed time and floating observations.

The [catalogue checker](../tools/check_f4_filter_catalogue.py) also
recomputes every one of the 3,850 representative profiles from its
recorded exponent and coefficient vectors. Its
[replay record](../evidence/f4-three-mode-catalogue-replay.json) confirms
the exact catalogue digest and the total coefficient multiplicity.

Reproduce the checks and the complete catalog from the repository root:

```bash
python tools/check_f4_weyl_filter_compiler.py \
  --output evidence/f4-weyl-filter-compiler.json
python tools/search_f4_repaired_identity_filters.py \
  --n 32 --modes 3 --geometry initial \
  --output /tmp/f4-three-mode-full.json
python tools/compact_f4_filter_screen.py \
  --input /tmp/f4-three-mode-full.json \
  --output evidence/f4-repaired-three-mode-initial-n32-screen.json
```
