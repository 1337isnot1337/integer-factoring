# Constant F4 flags for the complete Pascal Clifford evaluator

## Exact family

The [weighted Pascal flag](weighted-pascal-flag-crossings.md) constructor
now permits every constant `C in GL(L,F4)`, alongside its existing binary
flags. This changes the polynomial kernel while retaining the full typed
endpoint. It is distinct from changing a derivative inverse-system
physical code. Both complete cut profiles and both physical costs are
still recomputed at every search query.

Write the original blocks as `A,B,P` and put `s=u^(2L)`. They satisfy
`AB=(1+u)^(2L)P` and `BPA=(1+u)^(2L)I`. For any constant invertible
`C`, use

```math
 D_C=\begin{pmatrix}AC^{-1}&sP\\I&CB\end{pmatrix},
 \qquad
 D_C^{-1}=\begin{pmatrix}CBP&sI\\P&PAC^{-1}\end{pmatrix}.
```

Direct block multiplication gives the identity, since
`(1+u)^(2L)+u^(2L)=1`. The off-diagonal terms cancel in characteristic
two. At `u=1` the endpoint remains `[[0,P],[I,0]]`, independent of `C`.
The same polynomial inverse and Clifford formulas consequently apply
over the fixed field `F4=F2[omega]/(omega^2+omega+1)`.

The [constructor](../tools/f4_pascal_flags.py) retains the two exact
binary coefficient planes. It never interprets a field coefficient as
an integer multiplicity. The existing native engine already computes
exact polynomial arithmetic over this field. Its weight decisions and
the optimization remain floating calculations.

## Search interface and scalar initialization

Use `--flag-field 4` with `--flag-shears` in the
[mode-gauge driver](../tools/search_weighted_pascal_flags.py). Elementary
row additions use a nonzero F4 multiplier, so every generated flag is
invertible. The actual flag matrix is retained for reproduction.

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/search_weighted_pascal_flags.py \
  --images 0 1 2 3 --flag-field 4 --flag-shears 6 \
  --restarts 2 --iterations 180 --seconds 30 --seed 604402 \
  --output /tmp/f4-pascal-flag-search.json
```

The mode driver now fits its inherited common scale from one complete
Clifford query: multiplying all heights and the valuation parameter by
`t>0` multiplies all orders by `t`. It retains the scalar Gaussian
normalization as well as the profile partition functions. This avoids
repeating an identical pivot problem during scalar initialization.

The original-Fock symbol backends explicitly reject an F4 mode-start
record for now; their selected-minor implementations currently use
binary coefficients. The complete F4 mode-Clifford path is implemented
and is separately validated below.

## Independent controls and finite boundary

The [checker](../tools/check_f4_pascal_flags.py) enumerates all 256 two-by-two
coefficient matrices: all 180 invertible flags have both exact kernel
inverse identities, and all 76 singular flags are rejected. Every binary
specialization agrees with the earlier kernel and both complete Clifford
matrices, giving twelve complete specialization comparisons.

Fourteen specified invertible flags also have a separate packed-F4
subset-determinant reconstruction of the entire original Fock tensor.
Independent fraction-free polynomial Smith elimination agrees with
**448 complete original-Fock orders** under nontrivial integer physical
heights. Three additional numerical gradient controls are kept separate
from those exact identities. The
[receipt](../evidence/f4-pascal-flag-controls.json) retains the actual flags,
heights, original profiles, and source hashes.

Two small CLI controls, one over each field, reduce their numerical
pressure to approximately `0.81128`; neither approaches the half-score
threshold. No optimum or contraction improvement is claimed. Larger
flags and independently retimed physical heights remain a constructive
search space. The complete factoring exponent remains **1/2**.

## A degree-normalized symbol-gradient probe

The first larger 24-letter symbol run stopped at its startup guard.
The ordinary scale-coordinate central difference crossed neighboring
valuation regions: at step `10^-7`, its derivative discrepancy was
approximately `5.86*10^-4`; at `10^-9`, it was below `6.6*10^-8`.
The other four tested coordinates agreed at their original step.

The symbol driver now divides the scale probe by `outer*mode_count`,
matching the earlier mode driver's treatment of polynomial orders. The
same fixed 24-letter start then passes the original `10^-4` guard with
maximum discrepancy below `4.3*10^-7`. The actual probe sizes and all
coordinate comparisons are saved. The
[regression packet](../evidence/pascal-q24-scale-gradient-control.json.xz)
retains the original failure, the decreasing-step diagnosis and the
successful bounded CLI replay. This fixes the numerical guard; it does
not establish an exact pressure inequality.

## Pinned default verification

`make verify` passed at clean source
`566b596782afef3db6aba97c39678849e6446566`, tree
`63bab50efffe27c2066130a62904813d96a34117`: 341 logged Python command
blocks and unit-test suites of 66 and five tests. The
[receipt](../evidence/verify-566b596.json) and its compressed stdout retain
the exact source and output hashes. This pass predates the whole-symbol,
occupation-sector and F4 additions. Those have the focused checks stated
in their notes; a later-source full-suite pass is not claimed. No new
Lean target was run.

This work reuses and credits the preceding Pascal, F4 arithmetic,
Clifford, original-Fock and native Smith implementations. Substantial
OpenAI assistant support was used for the proof, code, computation and
text. No stronger uniform factoring bound is claimed.
