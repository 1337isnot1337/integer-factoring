# Pure physical corrections with correlations between copies

## Construction

Let `B` be the retained strict five-letter crossing of degree 18.
Group two copies into a 25-letter physical alphabet and choose four
arbitrary nonzero vectors `a,b,c,d` in `F4^25`. Define

```math
T_{XYUV}(z)=(B\otimes B)_{XYUV}(z)
             +h(z)a_Xb_Yc_Ud_V,
\qquad h(1)=0,\quad\deg h\le36.
```

The endpoint is exactly `[X=V][Y=U]`. Each physical vector is indexed
by a pair of original letters. It can be a full-rank five-by-five
matrix under that indexing. Thus this family includes correlations
between the two copies inside each physical label, beyond the
product-of-one-copy-support family in
[two-copy couplings](two-copy-crossing-couplings.md).

The update has rank one at each of the two required cuts. That fact
does not require any physical vector to factor across copies. The
complete matrices still have size 625 by 625. Their physical weights
and paired energies are exactly those inherited from the two copies.

## Applying the tensor-product regular frames

Let `L,R` be the 25-by-25 regular polynomial frames of one weighted
cut. Arrange a two-copy correction column as a 25-by-25 matrix `U`
in copy order. The transformed column is the vectorization of
`L U L^T`. For a correction row arranged as `V`, its transform is
the vectorization of `R^T V R`.

These are polynomial matrix sandwiches modulo the chosen absolute
jet precision. The [native helper](../tools/native_tensor_frame_f4.cpp)
uses ordinary truncated polynomial multiplication, with no division
or valuation approximation. The
[adapter](../tools/entangled_two_copy_frames.py) includes all four
physical weight factors before those transforms. It then calls the
[structured rank-one engine](rank-one-smith-search.md).

The helper's finite-jet identities are exact. The older structured
engine used by this initial search has the intermediate-precision
limitation stated in its note, so the search records count original
complete-matrix replays separately from structured predictions.
Every improving prediction is replayed before it is retained as an
improvement. Floating scale optimization still needs a rigorous
inequality certificate before changing the factoring theorem.

## A safe profile-cost pruning criterion

Let the two complete integer profiles be `p` and `q`. Their product
of exponential sums is the exponential sum of the multiset

```math
C=\{p_i+q_j:1\le i,j\le625\}.
```

Suppose `C` and the unchanged reference multiset `C0` have the same
cardinality, and every ascending prefix sum of `C` is at most the
corresponding prefix sum of `C0`. Then for every `s>0`,

```math
\sum_{x\in C}e^{-sx}\ge\sum_{x\in C_0}e^{-sx}.
```

Here is a direct proof. For any real `t`, the sum of positive parts
`sum(t-x)_+` equals the maximum over `k` of `kt` minus the ascending
`k`-term prefix sum, including `k=0`. The prefix assumption therefore
makes this positive-part sum no smaller for `C`. Finally,

```math
e^{-sx}=s^2\int_x^\infty e^{-st}(t-x)\,dt.
```

Integrate the positive-part inequality. Since the physical energy
lists remain unchanged, the resulting profile pair cannot improve
the pressure at any positive scale or proposed exponent. This is
an exact sufficient obstruction, not a necessary one.

The implementation convolves profile histograms and compares prefix
sums at their joint multiplicity boundaries. Between boundaries the
prefix difference is affine, so checking the endpoints checks every
integer prefix. The structured search labels this as *predicted*
domination unless the underlying profiles have been independently
replayed.

## Finite controls and reproduction

The [control record](../evidence/entangled-two-copy-frame-controls.json)
checks 206 full-precision polynomial matrix entries against the
independent Python multiplication routine, at precisions 128 and
256. Twelve physical examples give 24 original-matrix profile
comparisons and 70,936 independently reconstructed weighted entries.
The four large examples have genuinely correlated physical vectors:
their five-by-five ranks range from four to five.

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/check_entangled_two_copy_frames.py
OPENBLAS_NUM_THREADS=1 python3 tools/search_entangled_two_copy_updates.py \
  --cases 12000 --seconds 5400 \
  --output /tmp/entangled-two-copy-pure-search.json
```

The finite driver selects physical support windows by their actual
energies, uses the full binomial coefficients of `(1+z)^r`, retains
unresolved cases, and states whether the requested schedule finished.
This construction and its controls establish no stronger factoring bound.
