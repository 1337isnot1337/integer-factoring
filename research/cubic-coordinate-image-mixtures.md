# Coordinate-image output mixtures: an all-power barrier

## Result

Keep the cubic input map and cyclic output selector fixed.
Allow every complementary output map \(C=e_iv^T\), where the changed output row \(i\) and covector \(v\) may both vary.
Form the full combined coefficient space of all these maps, including linear combinations made before either cut rank is evaluated.

Over every characteristic-two field and at every tensor power \(n\), each decomposition of the crossing has cost
\[
 \boxed{\quad S_n\geq392^{n/2}>7^{\sqrt2 n}.\quad}
\]
This closes the last configuration in [the mixed F8 search](cubic-mixed-F8-spaces.md): the pair \(zE_{01},z^2E_{23}\).
It also permits mixtures across all seven changed-row choices, with arbitrary covectors.

The [earlier untouched-output gate](cubic-output-perturbations.md) applies to a fixed changed-row set.
The present result allows different row choices within one combined coefficient space, even when no row remains untouched across the entire mixture.

## Larger space and exact restrictions

Let \(\mathcal B\) denote the baseline coefficient space at \(C=0\), and write \(e(B)\) for its crossing coefficient.
The physical parent identities proved in [the fixed-image theorem](cubic-fixed-image-rank-one-barrier.md) show that the current family lies in
\[
 \mathcal W=
 \left\{B+U:\ B\in\mathcal B,\quad
 U(x,y,v,w)=\delta_{xw}K_x(y,v)\right\}.
\]
Here all seven matrices \(K_x\) may be arbitrary. Symmetry is not needed for this enlarged-space proof.

Use the fourteen unordered nonedges
\[
 \mathcal P=\{(a,b):0\leq a<b\leq6,\ b-a\notin\{1,6\}\}
\]
and the twenty-eight oriented nonedges \(\mathcal O=\mathcal P\cup\mathcal P^T\).
The actual baseline matrices satisfy
\[
 R(B)|_{\mathcal P\times\mathcal P}=e(B)I_{14},
 \qquad
 C(B)|_{\mathcal O\times\mathcal O}=e(B)I_{28}.
\]

In the column restriction, every row index \((x,w)\in\mathcal O\) has \(x\ne w\).
Thus every update vanishes on the entire restriction, leaving the scalar matrix \(e(T)I_{28}\).

In the row restriction, a nonzero update from column \((v,w)\) to row \((x,y)\) requires
\[
 v<w=x<y.
\]
Ordering \(\mathcal P\) by decreasing first label makes every update strictly upper triangular.
The full restricted row matrix is therefore upper triangular with all fourteen diagonal entries equal to \(e(T)\).
There are ten possible strict upper positions for the entire enlarged update class.

These identities also make \(e\) a well-defined linear functional on \(\mathcal W\), regardless of how a tensor is represented as baseline plus update.

## Tensor powers

At tensor power \(n\), the row restriction remains upper triangular and the column restriction remains scalar.
Every diagonal entry in both cuts equals \(e^{\otimes n}(T)\), even for arbitrary recombinations of tensor-product directions.
Hence
\[
 e^{\otimes n}(T)\ne0
 \quad\Longrightarrow\quad
 \operatorname{rank}R(T)\geq14^n,
 \qquad\operatorname{rank}C(T)\geq28^n.
\]
Because \(e(E)=1\), at least one piece in any crossing decomposition has nonzero value of this functional.
The resulting cost lower bound is \(\sqrt{14^n28^n}\).
The comparisons \(392>7^3\) and \(3^2>8\) establish the strict seed barrier.

## Reproduction and limits

The [checker](../tools/check_cubic_coordinate_image_mixtures.py) verifies the baseline restrictions and all 343 matrix-unit updates in the enlarged class.
It also checks the complete seven-direction space from the recorded two-row-choice configuration directly.
The [certificate](../evidence/cubic-coordinate-image-mixtures.json) records the row ordering, column coordinates, ten permitted update positions, and source hashes.

```bash
python tools/check_cubic_coordinate_image_mixtures.py \
  --output /tmp/cubic-coordinate-image-mixtures.json
```

The theorem covers recombinations between maps that each change one complementary-output row.
A single complementary matrix with several nonzero rows also creates cross terms between those rows.
Those cross terms need not be diagonal in the two output indices and are not covered by this theorem.
