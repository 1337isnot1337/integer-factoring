# Scalar reweighting does not repair the single-block favorable cube family

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338747)
relaxes the unit-coefficient rule in the rejected
[global-pairing construction](two-axis-global-pairing-repair.md). Could
new scalar coefficients restore the complete endpoint while retaining only
orientations that hide a maximum-width axis?

For the fixed actual Gaussian coefficient family used in that audit, the
answer is **no**, even allowing every tied favorable orientation and all
scalar coefficients over every characteristic-two extension field. This is
a finite single-block statement; it is not an asymptotic obstruction or an
exclusion of coupled blocks. The complete factoring score stays **1/2**.

## Exact family and certificate

Use exactly the original m=5, weight-one alphabet `{00001,00010}`. If G
is the eight-mode Gaussian state, its tensor entries are the polynomials
`G_mask*G_(mask xor 255)`. There are 22 nonzero coefficient tensors A_i,
each with 256 binary entries. For each, compute its three exact cut ranks
and widths `ceil(log2 rank)`, as in the predecessor audit. An orientation
is favorable if its hidden axis has maximum width, including ties.

With the fixed four-parent labeling, mask 240 means
`X=Y=(0,0), U=V=(1,1)`; mask 255 means all four labels equal `(1,1)`.
The linear functional

~~~
L(T)=T[240]+T[255]
~~~

vanishes on **all 30 favorable reoriented coefficient tensors**.
The endpoint has values zero and one at these two assignments, so
`L(E4)=1`. Therefore E4 is outside their linear span. The coefficients
and functional are binary, so the same separating equation works over
every characteristic-two extension, without enumerating possible scalar
weights. This permits several favorable orientations of the same piece,
and is more permissive than choosing exactly one.

The [checker](../tools/check_cube_favorable_reweighting.py) reconstructs
the actual family, all reoriented truth tables and exact cut widths.
Independent binary elimination gives dimension eight for the 22
first-choice generators and again eight for all 30 favorable generators.
The unrestricted 66-generator span has dimension fourteen and contains
E4; an explicit solution is retained as a positive control. The
[receipt](../evidence/cube-favorable-reweighting.json) records every
generator and all three span outcomes.

Run `make verify-cube-favorable-reweighting`. The explicit two-entry
functional checks the negative conclusion directly, independently of
the rank of the span. No pressure optimization or new Lean is involved.

## Exact limitation: the functional does not tensorize under global selection

The favorable condition above is per block. For products, the global
width rule hides an axis maximizing the **sum** of block widths. That
does not require it to maximize every block separately.

For example, two original coefficient pieces have widths `(4,0,0)` and
`(0,4,0)`. Their joint widths are `(4,4,0)`, so hiding axis zero is
globally favorable. The corresponding two-block tensor need not satisfy
the single-block two-entry relation. Consequently this certificate cannot
be promoted to a coupled-block exclusion by tensoring the statement.
Whether a global coupled span can restore the endpoint requires its own
identity test and, if successful, a uniform construction and full rank
and recursion accounting.
