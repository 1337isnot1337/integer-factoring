# Nonconstant regular flags for complete binary crossing synthesis

## Result and scope

The two-cut affine solver now accepts arbitrary polynomial row frames
over F4 whose constant matrices are invertible. This enlarges the earlier
[constant-flag search](binary-two-cut-affine-synthesis.md) without imposing
occupation preservation or a Gaussian form on the unknown tensor.

The first completed catalogue uses every frame `I+zN` with strictly lower
triangular N at each cut. All 59 specified targets were decided exactly.
Only the reference target is feasible, and its complete solution fiber
is the same four-point affine line already present with constant frames.
This finite search supplies no stronger factoring exponent. Higher jets,
other constant flags, other targets, and other physical height rays are
outside this completed catalogue.

## 1. Polynomial frames keep the problem affine

Use the complete degree-at-most-12 endpoint fiber

```math
T(w)=w^4B(w)+(1+w)Q(w),\qquad \deg Q_{xyuv}<12,
```

where B and the two fixed physical height vectors are specified in the
earlier note. There are 192 independent coefficients over F4. The
literal endpoint is the crossing permutation on all sixteen entries.

At each place, let A(T) be the complete weighted polynomial flattening,
and let U,V be the fixed regular left and right reference matrices.
Choose a polynomial matrix W(z) with invertible W(0). For each row i,
impose

```math
\operatorname{ord}_z\bigl(W(z)U(z)A(T;z)V(z)\bigr)_{ij}
\ge \tau_i\qquad\text{for every column }j.                 \tag{1}
```

Every coefficient below its row target is affine in the 192 unknowns.
Multiplying by a prescribed nonconstant W changes the coefficient
convolution but does not change linearity. The solver forms the complete
convolutions before extracting these equations. It rejects a singular
constant matrix rather than treating a nonunit row operation as free.

If the sorted row targets are t_1,...,t_4, every k-by-k minor in (1) has
order at least t_1+...+t_k. Since W,U,V and their inverses are regular,
the same lower bound holds for the original determinantal ideal. Thus
the actual Smith profile satisfies every corresponding prefix bound.
The explicitly recorded common scalar shift is removed afterward.

The equations at zero and infinity share the same unknown coefficients.
Their intersection constructs an actual tensor, rather than two
unrelated local matrices. Each retained solution is reconstructed from
its full polynomials and checked with the independent Smith profiler.

## 2. Finite jets suffice for a fixed row-lattice target

Suppose every entry of U A(T) V has order at least m, and every target
is at most T_max. For a regular formal frame W, retain only its terms
of degree less than T_max-m. The discarded part has order at least
T_max-m, so its product with U A(T) V has entry order at least T_max.
It cannot affect any equation in (1). The truncated frame has the same
invertible constant matrix as W. If T_max=m, the common entry floor
already supplies all the target inequalities.

This proves finite-jet sufficiency for the prescribed targets. It does
not say that lower triangular jets represent every regular frame.

The common floor can be imposed before choosing W: a regular matrix and
its regular inverse preserve the minimum entry order. The Laurent
support lemma from the earlier note therefore continues to hold for
nonconstant frames. In particular, the current target catalogue needs
only support from -4 through 8 after removing the common clearing.
Degree 12 with clearing exponent 4 contains the complete coefficient
fiber for those targets. This coefficient completeness is distinct from
the finite scope of the chosen frame catalogue.

## 3. Completed first-jet catalogue

At each cut choose

```math
W(z)=I+zN,\qquad N\text{ strictly lower triangular over }\mathbb F_4.
```

There are 4^6=4,096 frames per cut. Every determinant is exactly one.
The reference target is `(5,8,10,13)` at each place. The other targets
are its eight single-row increases and the sorted balanced transfers
of up to three units retained by the driver's explicit numerical
critical-pressure screen. The saved manifest records the resulting
59-element list. The screen defines the finite schedule; it supplies
no assertion about omitted targets.

The complete record is
[binary-polynomial-flags-first-jet.json](../evidence/binary-polynomial-flags-first-jet.json).
Exact common-floor and single-cut elimination decide groups of frame
pairs at once. The search covers 989,855,744 logical target/frame-pair
combinations and performs 945 direct intersections after pruning.
These counts describe different operations and should not be conflated.

Only one frame pair is feasible, at the reference target. Its affine
dimension is one over F4. All four points have restored full profiles

```text
infinity: (-9,-6,-4,-1)
zero:     ( 3, 6, 8,11).
```

The common clearing shifts these oppositely relative to the unshifted
reference tensor, leaving the product cost unchanged. The numerical
fixed-ray rate is approximately 0.8040228107, well above the required
1/sqrt(2). The rate is reported as a numerical search result, not as a
new certified constant.

## 4. Independent controls and reproduction

The checker uses dictionary polynomials and literal matrix multiplication
to compare 480 transformed rows from arbitrary invertible constant
frames with polynomial terms through degree four. It also checks 256
constant-frame identities, 1,024 literal endpoint coefficients, and
34 complete all-minor cut profiles. It enumerates all four points of
the surviving affine line, checks irrelevant high jets, and rejects a
singular constant frame. The saved
[control record](../evidence/binary-polynomial-flag-controls.json)
includes source hashes.

```sh
make verify-polynomial-flags
OPENBLAS_NUM_THREADS=1 python3 tools/search_binary_polynomial_flags.py \
  --degree 12 --baseline-shift 4 --field-order 4 --jet-power 1 \
  --transfers 3 --seconds 1800 \
  --output /tmp/binary-polynomial-flags-first-jet.json
```

The complete source is
[search_binary_polynomial_flags.py](../tools/search_binary_polynomial_flags.py).
Neither this constructor nor its finite catalogue requires new Lean
formalization.
