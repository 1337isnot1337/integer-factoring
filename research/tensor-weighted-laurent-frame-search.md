# Literal Laurent frames and completed coefficient searches

This checkpoint extends the literal Pfaffian constructor to finite
parity-preserving Laurent physical maps. It also retains the completed
large-kernel and finite-field searches listed below. None improves the
current factoring exponent or the compact pressure benchmark. The
numerical pressure fits are separate from the exact capped Smith data.

## 1. A split Laurent frame

Order the ambient occupation words as retained words followed by one
anchor word. For Laurent column vectors \(a,b\), set

\[
 L=[I\;a],\qquad
 P=\begin{bmatrix}I+ab^T\\b^T\end{bmatrix}.
\]

In characteristic two, \(LP=I\) identically. The implementation checks
this identity by sparse Laurent-polynomial multiplication for every
constructed frame. This is stronger than checking only the parameter-one
endpoint.

The retained occupation words have their natural heights
\(h(S)=\sum_{j\in S}(2j+1)\), with parameter \(t=w^2\). If the anchor
has height \(h_*\), the tested flow assigns

\[
 a_i=\alpha_i t^{(h_*-h_i)/2},\qquad
 b_i=\beta_i t^{-(h_*-h_i)/2}.
\]

Every supported word has the anchor's occupation parity. Thus the
exponents are integral, and the two parity blocks can be computed
separately. The native constructor currently requires this parity
condition. It does not support arbitrary parity-mixing finite maps.

The zero-place contraction includes the physical height gauges on its
two output legs. The infinity-place contraction reverses every Laurent
exponent in every physical map, applies the prescribed opposite-leg
gauges, and contracts the literal reversed kernel. The code clears each
cut explicitly and subtracts that cut's clearing power from every
returned Smith order. It does not assume that independently chosen local
clearing powers cancel.

The relevant implementations are
[`tensor_literal_weighted_maps_f16.py`](../tools/tensor_literal_weighted_maps_f16.py),
its [native constructor](../tools/tensor_literal_weighted_maps_f16.cpp), and
[`tensor_literal_moving_frames.py`](../tools/tensor_literal_moving_frames.py).
They extend the constant-map literal constructor described in
[`tensor-literal-low-cell-retractions.md`](tensor-literal-low-cell-retractions.md).

## 2. Independent finite controls

The [checker](../tools/check_tensor_literal_moving_frames.py) contracts
the original raw Laurent tensor in
[`reflection_dual_rail_gaussian.py`](../tools/reflection_dual_rail_gaussian.py).
It does not obtain its reference coefficients from the new weighted
native constructor. The retained
[receipt](../evidence/tensor-literal-moving-frame-controls.json) contains:

- 207,360 complete reference coefficient comparisons;
- 1,101 independent binary Toeplitz-rank checks;
- 4,320 capped Smith-order comparisons;
- 103,680 endpoint checks and 103,680 infinity-reciprocity checks.

Replay without replacing the retained receipt:

```sh
python tools/check_tensor_literal_moving_frames.py \
  --output /tmp/tensor-literal-moving-frame-controls.json
```

## 3. Completed moving-frame census

The ambient alphabet contains all occupation words of height at most 14.
The anchor is word 9, of height 8. The four supported retained words are
\(0,3,5,6\). Deleting the anchor leaves 23 physical letters. The complete
binary census tests all 16 support masks for each of \(a\) and \(b\),
giving 256 frames. Both cuts are computed directly with actual-order
cap 48. All 270,848 returned capped orders are retained in the
[complete census](../evidence/tensor-literal-moving-frames-height14-flow1-binary-census.json).

The best observed fitted upper-counting root is
\(0.7161709647129381\), attained at masks \((0,8)\). This winning frame
only mixes the two equal-height words 6 and 9. The broader Laurent flows
in this census gave no better observed value. This is a finite search
result, not an exclusion theorem for arbitrary moving frames.

The earlier [25-case binary pilot](../evidence/tensor-literal-moving-frames-height14-flow1-profiles.json)
and [25-case F4 pilot](../evidence/tensor-literal-moving-frames-height14-flow1-f4-pilot.json)
are also retained. The F4 pilot uses coefficient 2 on each selected
support; it is not an exhaustive F4 census. Both pilots have best observed
upper-counting fit \(0.716170986060009\).

## 4. Completed F16 low-cell census

The [F16 receipt](../evidence/tensor-literal-f16-rank-one-cells-height14-profiles.json)
extends the constant two-word rank-one census in
[`tensor-literal-low-cell-retractions.md`](tensor-literal-low-cell-retractions.md).
At each of heights 8 and 9, it covers every normalized pair
\(P,L\) with \(LP=1\) in a two-dimensional height cell. There are 272
normalized F16 pairs per cell and 76 representatives under Frobenius.
Thus the 152 completed representatives cover 544 normalized pairs.
All other words of height at most 14 remain unchanged, leaving 23 letters.
Both cuts are computed directly at order cap 64.

| Changed height | Best observed upper-counting fit | Winning pair in the ordered height cell |
|---|---:|---|
| 8 | 0.7161709624583944 | \(P=(1,0)^T,\ L=(1,1)\) |
| 9 | 0.7149501056884414 | \(P=(0,1)^T,\ L=(0,1)\) |

These are the same best observed values already present over F4. The
exhaustiveness refers to the declared normalized pairs, not to the
one-variable numerical minimization used for each pressure fit.

## 5. Completed larger profiles

All four previously running profiles below finished. Every profile uses
order cap 64. The first two and the orthogonal height-9 profile compute
one cut and derive the other through the documented orthogonal identity.
The coordinate height-9 profile computes both cuts directly.

| Family | Physical letters | Fitted lower-counting root | Fitted upper-counting root | Receipt |
|---|---:|---:|---:|---|
| F16 four-pair kernel, height at most 27 | 115 | 0.7105996597892456 | 0.7106153481942494 | [profile](../evidence/tensor-f16-fast-four-pair-height27-profiles.json) |
| F4 mixed cubic graph, height at most 25 | 91 | 0.7110023549600893 | 0.7110060452527294 | [profile](../evidence/tensor-f4-cubic-graph-height25-profiles.json) |
| F4 orthogonal height-9 replacement, height at most 20 | 55 | 0.7138884069117858 | 0.7138884274589950 | [profile](../evidence/tensor-literal-height9-f4-orthogonal-cap20-profiles.json) |
| Binary coordinate height-9 retraction, height at most 20 | 55 | 0.7138884075004080 | 0.7138884280475133 | [profile](../evidence/tensor-literal-height9-binary-coordinate-cap20-profiles.json) |

The kernel constructions and checks are documented in
[`tensor-f16-capped-search.md`](tensor-f16-capped-search.md),
[`tensor-cubic-graph-orthogonality.md`](tensor-cubic-graph-orthogonality.md),
and the literal low-cell note. Every retained receipt records source
hashes and its completed status.

For clarity about capped tails, if the true orders are \(\delta_i\) and
the cap is \(T\), then, for \(0<z<1\),

\[
 \sum_{\delta_i<T}z^{\delta_i}
 \;\leq\; Z(z)
 \;\leq\;\sum_i z^{\min(\delta_i,T)}.
\]

The lower and upper columns fit the corresponding counting functions.
The exact inequalities do not make the floating-point searches into
certified global minima. An order recorded at the cap is an unresolved
tail entry, not an assertion that its actual order equals the cap.

## 6. Scope of the checkpoint

The positive implementation result is a literal constructor for
parity-preserving Laurent physical frames with independently checked
reciprocity and explicit clearing costs. The completed searches preserve
their negative numerical results and exact capped data. No stronger
factoring exponent or certified pressure improvement follows from this
checkpoint.
