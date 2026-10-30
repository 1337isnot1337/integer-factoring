# Varying rank-one images: an affordable-span escape

## Positive result and remaining gate

Allowing both vectors in \(C=uv^T\) to vary across output maps defeats the fixed-image affordable-span obstruction.
A concrete set of **154 binary rank-one maps** has a full mixed coefficient space in which individually affordable tensors span the crossing.
An explicit decomposition uses 22 affordable pieces with total cost 49.
The same space has a three-piece decomposition of total cost **36 over F2**.

Neither displayed decomposition meets the strict seven-label budget \(7^{\sqrt2}\), which is below 16.
The result passes a necessary span condition and identifies a larger search space; it does not establish a strict seed, a minimum cost for this space, or a factoring exponent improvement.

This scope differs from [the fixed-image barrier](cubic-fixed-image-rank-one-barrier.md): the image vector now varies within one local coefficient space.
Non-coordinate image vectors are essential to this construction, so it also lies outside [the coordinate-image mixture theorem](cubic-coordinate-image-mixtures.md).

## A finite family

For every \(0\leq i<j\leq6\), define
\[
 U_{ij}=\{e_i,e_j,e_i+e_j\}.
\]
Include the output maps \(C=uv^T\) for all \(u,v\in U_{ij}\), and deduplicate across pairs.
There are 154 distinct nonzero rank-one binary matrices: the nine maps per pair share the seven diagonal maps across pairs.
The construction uses the fixed cubic input and the standard cyclic output selector.

The raw \(q^3\) coefficient of a complementary output map is
\[
 H_C(x,y,v,w)=C_{xy}C_{wv}.
\]
When \(C=uv^T\), its column flattening is the outer product
\[
 (u\otimes u)(v\otimes v)^T.
\]
Both flattenings have rank one, so every one of these 154 coefficient tensors has cost one.

## Explicit affordable pieces

Write
\[
 a_{ij}=e_i\otimes e_j+e_j\otimes e_i.
\]
In characteristic two,
\[
 a_{ij}=(e_i+e_j)^{\otimes2}+e_i^{\otimes2}+e_j^{\otimes2}.
\]
Consequently, the tensor whose column flattening is \(a_{ij}a_{ij}^T\) is the sum of the nine \(H_C\) tensors from that pair.
Call this tensor \(A_{ij}\).
Its column rank is one, and its row rank is four: the corresponding two-by-two alternating matrix has rank two in each of the two row factors.
Thus \(c(A_{ij})=2\).

Let \(H\) be the noncrossing pairing
\[
 H(x,y,v,w)=\delta_{xy}\delta_{wv}.
\]
Its row and column ranks are \((1,49)\), so its cost is seven.
In the column cut, \(H\) is the swap operator and \(E\) is the identity.
Therefore
\[
 E+H=\sum_{i<j}A_{ij}.
\]
The crossing \(E\) is a raw coefficient of every allowed output map, so this identity also proves that \(H\) belongs to the full mixed space.
In an expanded membership witness, 147 distinct \(q^3\) tensors remain after the seven shared diagonal-map terms cancel.

We have the decomposition
\[
 E=H+\sum_{i<j}A_{ij},\qquad
 S=7+21\cdot2=49.
\]
All 22 pieces have cost strictly below \(7^{\sqrt2}\).
Thus the affordable-piece span really contains the crossing, despite the excessive total cost of this particular decomposition.

## A cost-36 binary decomposition

Let the other noncrossing pairing be
\[
 N(x,y,v,w)=\delta_{xw}\delta_{yv}.
\]
The finite family includes \(C=e_0(e_0+e_1)^T\), which satisfies \(CC^T=0\).
The general complementary-output identity gives
\[
 N=A(C)+D(C)+H_C,
\]
where \(A(C)\) is the raw \(pr^2\) coefficient and \(D(C)\) is the raw \(qr^2\) coefficient.
Thus \(N\) also belongs to the full mixed space.

The three tensors
\[
 H,\qquad N,\qquad E+H+N
\]
sum to \(E\), with rank pairs
\[
 (1,49),\qquad(49,1),\qquad(22,22).
\]
Their total cost is \(7+7+22=36\).
The earlier [scalar F4 output construction](cubic-complementary-output-maps.md) attained the same cost within a single output-map space.
The present witness attains it over F2 by combining rank-one output-map spaces.
No claim is made that 36 is optimal in the enlarged mixed space.

## Exact verification

The [checker](../tools/check_cubic_varying_image_escape.py) regenerates all five physical raw coefficient matrices for every one of the 154 maps.
It checks every \(q^3\) entry against the direct outer-product formula in both cuts, totaling 739,508 entry comparisons.
It verifies the nine-term edge identities, the crossing identity, the memberships of \(H\) and \(N\), and every displayed rank using two pivot orders.
The [certificate](../evidence/cubic-varying-image-escape.json) records every output map, every edge-piece witness, and all decomposition data.

```bash
python tools/check_cubic_varying_image_escape.py --output /tmp/cubic-varying-image-escape.json
```

The remaining task in this larger space is a total-cost improvement, not another affordable-span test on the same family.
