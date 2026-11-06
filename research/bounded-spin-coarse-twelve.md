# Twelve-site coarse transposition spectrum

## Result

The [transposition subalgebra note](bounded-spin-transposition-subalgebra.md)
left the twelve-site coarse computation open.  Run the same construction:
let `R` be the replica swap and let

```math
T=\sum_{i<j}T_{ij}\otimes T_{ij}
```

be the compressed diagonal transposition class sum on the `89^2=7921`
path-pair module over `F_5`.  The six nonzero joint generalized
eigenspaces of `(T,R)` have dimensions

```text
(T=0,R=+):1320   (T=0,R=-):165
(T=1,R=+):2015   (T=1,R=-):1266
(T=4,R=+):670    (T=4,R=-):2485
```

They sum to `7921`.  The coarse eigenvalues are `0,1,4`; the ten-site
coarse spectrum used `0,1,2`.  The smallest block, `165`, is below the
twelve-site label target `89^sqrt(2)=571.286...`, but this is a module
dimension, not a cut-rank pair.  The other five blocks are above the
target.  As at ten sites, the full coefficient algebra can refine a
coarse block: the ten-site coarse `335` splits as `34+1+300`.

## What this does and does not settle

This closes the explicit twelve-site coarse computation in the
transposition note.  It shows that `R` and `T` alone do not give a
favorable primary partition at twelve sites: the coarse blocks are not
the full central idempotents, and one cannot read a two-cut rank sum
from their dimensions.  The full twelve-site primary computation still
needs the additional coefficient directions and the exact rank pairs.

No score changes; the complete bound remains `1/2`.

## Reproduction

The full replay is optional and takes about two hours in this checkout:

```sh
python3 tools/check_bounded_spin_coarse_twelve.py --recompute \
  --output evidence/bounded-spin-coarse-twelve.json
```

The fast evidence check is

```sh
make verify-bounded-spin-coarse-twelve
```

It checks the six recorded dimensions, their sum, and the recorded
scope.  The replay uses the same frame, transposition, class-sum and
generalized-eigenspace code as the eight- and ten-site note.
