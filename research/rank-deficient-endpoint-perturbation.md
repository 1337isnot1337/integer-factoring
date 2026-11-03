# Lowering the zero-scale rank does not improve the q4 natural pressure

## Result and scope

The q4 Gaussian identity endpoint has natural pressure
`0.7310809359736006` at scale `0.3290111256807374`. A rank-one
perturbation that subtracts a tensor at degree zero and restores it at
degree `d` preserves the endpoint at `w=1` while lowering the constant
term. In a deterministic 9,000-trial random screen over
`d=1,2,3,4,6,8`, no perturbation improved the natural pressure. The best
value in every degree was the unperturbed base value.

This is a finite numerical screen, not a global optimum or a pressure
certificate. It does not exclude a different rank-deficient seed, a
non-rank-one perturbation, or a nonlocal change of the low block. The
complete factoring exponent remains `1/2`.

## Mechanism

The q4 seed is the strict Gaussian tensor `tensor(4)` from
[`reflection_dual_rail_gaussian.py`](../tools/reflection_dual_rail_gaussian.py).
Its natural gauge uses the letter energies `[0,12,4,16]`; the base
orders give the recorded natural pressure. The zero-scale rank of a cut
is not by itself the pressure: the complete Smith orders and the height
gauge determine the support value.

For a rank-one `F4` tensor `R`, subtract `R` at degree zero and add `R`
at degree `d`. At `w=1` the two contributions cancel, so the endpoint
`[x=v][y=u]` is unchanged. The constant term changes from the identity
cut to `I-R_cut`, lowering the zero-scale rank by one on each cut. The
added degree-`d` term contributes positive orders. The screen asks
whether the rank drop can beat the added order cost at the natural
gauge.

The answer in the tested family is no. The full-rank identity endpoint
is locally optimal under this perturbation. The known rank-deficient
cyclic Pauli cell has unequal zero-scale cut ranks `(1,9)` but its
recorded support is `1.412832752034821`, far above the q4 base; unequal
rank alone is not a crossing mechanism.

## Reproduction

```sh
python3 tools/check_rank_deficient_endpoint.py \
  --output evidence/rank-deficient-endpoint-perturbation.json
make verify-rank-deficient-endpoint
```

The checker fixes the random seed, reconstructs the q4 base, and replays
the screen with the existing rectangular Smith oracle and natural
pressure fit. It records the best value and trial index for each
compensation degree and rejects any value below the base.

## Relation to the protected-core and q96 results

The [height-38 protected-core floor](jensen-tail-protected-core-floor.md)
and its [height-39 transition](protected-core-height39-transition.md)
apply to the retained height-below-24 core. The q4 seed here is a
different low block, so those floors are not asserted for it. The
[asymmetric zero-scale identity](asymmetric-degenerate-scaling-floor.md)
shows that rectangular scale and angle changes approach one from above;
this note tests the remaining rank-lowering perturbation of the same
identity endpoint. Neither route crosses. A height-39 support or a
genuinely different low block remains the structural escape, not another
scale, angle, or constant-term perturbation.
