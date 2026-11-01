# Coherent physical retractions over F8

This extends the [F4 coherent constructor](coherent-hot-f4-retractions.md)
to `F8=F2[x]/(x^3+x+1)`. It supplies exact low-order profiles with explicit
unknown tails. It does not, by itself, establish a better factoring bound.

## Construction and representation

Sparse vectors map occupation words to field elements encoded by integers
zero through seven. The three bits are coefficients of `1,x,x^2`.
Pairings are bilinear, not Hermitian. The code validates `LP=I`, common
height and occupation for each paired label, and deletion closure of both
map spaces.

The field-extension argument is the same as in the F4 construction.
All original Gaussian identities are binary. Extension to F8 preserves
them. The natural cuts use rows `L tensor L` and columns `P tensor P`,
or rows `L tensor P` and columns `P tensor L`, with the paired swap.
Their endpoint identities follow from `LP=I`.

Each hot-kernel minor remains a binary polynomial. The native engine
computes it with the existing height-based degree bound, multiplies its
coefficient by the map scalars, and accumulates **three binary polynomial
planes**. The field generator is not substituted for the series variable.
F8 multiplication reduces `x^3=x+1` and `x^4=x^2+x`.

The coefficient engine and exact truncated DVR elimination are in
[`coherent_hot_compounds_f8.cpp`](../tools/coherent_hot_compounds_f8.cpp).
The Python interface is
[`coherent_hot_compounds_f8.py`](../tools/coherent_hot_compounds_f8.py),
with map validation in [`coherent_f8_maps.py`](../tools/coherent_f8_maps.py).

## Known orders and unknown tails

Given cutoff `K`, the engine returns every Smith order at most `K`,
plus the number of unresolved channels. It does not declare unresolved
channels zero or infer rank deficiency from truncation.

For `0<t<1`, if `Z_known(t)` sums the resolved orders and `r` channels
remain, then

```math
Z_{\rm known}(t)\le Z_{\rm full}(t)
 \le Z_{\rm known}(t)+r t^{K+1}.
```

The endpoint supplies the required full generic rank in the validated
physical-code construction. If the tail count is zero, all orders are
resolved. The returned `all_orders_resolved` flag states exactly this.
Entry-order pruning is used only where it is invisible modulo `t^(K+1)`;
it is not a substitute for the Python endpoint and closure checks.

## Independent controls

Run:

```sh
python3 tools/check_coherent_hot_f8.py
```

The [receipt](../evidence/coherent-hot-f8-controls.json) records:

- All 64 scalar products against a separate binary-extension implementation.
- 16,884 complete polynomial entries checked by subset dynamic programming
  on the original binary one-particle matrix, using 1,676 distinct minors.
- 70 complete block profiles checked by
  [restriction to binary matrices](binary-extension-smith-restriction.md).
- Exact retained-order and tail comparisons at cutoffs 8, 16 and 24.
- 70 profile comparisons under field Frobenius, and nine rejected invalid
  map or block inputs.

These controls use all seven nonzero F8 coefficients in an eight-label
mixed-pair code. They validate the engine, not a universal optimization
claim. Larger code experiments require their own construction receipts,
complete-cost accounting, and exact pressure certificates.
