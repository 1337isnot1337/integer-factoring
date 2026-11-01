# Mixing the recorded cubic output spaces

Different output maps have the same logical crossing tensor, so their
coefficient pieces can be mixed in one decomposition. Checking each
space separately does not exclude that possibility.

For the 680 maps in
[the multirow cubic experiment](cubic-arbitrary-output-identity.md), a
single coordinate excludes all mixtures of the individually affordable
binary pieces. The statement is limited to those recorded pieces; new
linear combinations across spaces remain an open search direction.

## A one-entry certificate

Collect every nonzero binary coefficient tensor in each recorded space
whose cut ranks satisfy

\[
r(T)c(T)<256.
\]

After removing duplicates there are 659 tensors. Their binary span has
dimension 341. Every one of them satisfies

\[
T_{0,3,0,3}=0,
\]

whereas the crossing has `E_(0,3,0,3)=1`.

Consequently no linear combination of these tensors equals `E`. This
coordinate argument also allows arbitrary characteristic-two scalars
on the fixed binary tensors. It does not add the new tensor directions
that arise from arbitrary extension-field coefficient vectors.

The integer gate is sufficient because

\[
16^{70}>7^{99},\qquad99^2>2\cdot70^2.
\]

Any piece of a decomposition with total cost below `7^sqrt(2)` must
have individual cost below 16 and therefore occurs in the collected
pool if it is a binary coefficient tensor from one of the 680 spaces.
Thus mixing those spaces' individually affordable binary pieces cannot
give a strict seed.

## What remains open

The union of all variable directions `A_i,B_i,N`, including expensive
ones, has dimension 819 and **does contain the crossing**. The separating
coordinate does not annihilate this larger space. A linear combination
of expensive tensors from different spaces could have smaller cut
ranks than its summands; such newly formed tensors were not enumerated
in the 659-piece pool.

The [complementary-output construction](cubic-complementary-output-maps.md)
also changes the fixed tensor directions by allowing `Pi Omega f` to be
nonzero. It is a separate escape from this recorded-pool certificate.

## Reproduction

```sh
python tools/check_cubic_affordable_pool.py \
  --input evidence/cubic-multirow-outputs.json \
  --output /tmp/cubic-affordable-pool.json
```

The checker regenerates every tensor, checks the separating entry,
deduplicates by the full tensor, and recomputes both span dimensions.
The underlying rank tables have the independent complete replay
described in the multirow note.

- [Pool checker](../tools/check_cubic_affordable_pool.py).
- [All pool representatives and the certificate](../evidence/cubic-affordable-pool.json).
