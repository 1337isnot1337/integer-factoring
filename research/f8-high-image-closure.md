# Height-25 boundary-image freedom is already closed

## Question

The q96 F8 hyperoval map keeps the 92 coordinate words of height below
25 and adds four paired labels from the eleven height-25 words (the
height-25 word `31` and three hyperoval boundary embeddings). Can a
different rank-three boundary image on the same ten boundary triples
lower the pressure below `1/sqrt(2)` and cross the half exponent?

## Answer: no, for every height-25 boundary image

The [common-core natural-pressure floor](hot-boundary-natural-pressure-floor.md)
fixes the 92 height-below-25 coordinate labels and allows **four
arbitrary paired labels from the eleven-dimensional height-25 space**,
with any height-preserving embedding and retraction satisfying `LP=I`,
possibly different on the two physical sides. Over every
characteristic-two field it proves

```math
\rho \ge \frac{17}{24} > \frac1{\sqrt2}.
```

The boundary maps may be nonbinary, nonorthogonal, and may mix the
three- and five-particle words. The q96 map is exactly one instance:
its four height-25 labels are word `31` and the three hyperoval
embeddings. Consequently **every** replacement of the boundary image,
including a different rank-three triple basis, a different retraction,
or an arbitrary mixing over F8, is excluded by this theorem. No
boundary-only q96 variant can cross by the natural two-cut pressure
route.

This also explains two earlier screens:

- the 141 same-height/charge orthogonal `2 x 2` label mixings gave the
  identical cutoff-31 value `0.715030753621605`;
- 80 random rank-three images on the same pivots gave the identical
  cutoff-31 value.

Those are confirmations of the common-core plateau, not independent
evidence and not a new optimization lever. In particular, a
constant invertible change of the code-label basis preserves both cut
Smith profiles exactly (the cut matrices are multiplied by the same
invertible constant matrix on both sides), so basis changes cannot move
the pressure even without the floor theorem. The floor theorem is the
stronger statement, covering genuinely different images.

## Cutoff-47 control for two different images

As a finite control, the hyperoval image and one random rank-three image
on pivots `(7,8,9)` were resolved at cutoff 47 for q96. The 935 retained
orders per cut have the **same SHA-256 hash** for both cuts and both
images (`7ae6d94a...27eeb`), with 8,281 unresolved tails per cut.

| quantity | hyperoval | random image |
| --- | ---: | ---: |
| known-order lower `rho` | `0.7098760177608998` | `0.7098760177608998` |
| tail-augmented upper `rho` | `0.7107288703954673` | `0.7107288703954673` |

The known-order fit is a lower bound: unresolved orders only add
positive terms to `Z`, so the full pressure cannot fall below it. It is
already above `1/sqrt(2)` by about `2.77e-3`. The complete q96 profile
is `0.7100057872`, between the two bounds. Thus no image in this tested
class crosses; the universal exclusion remains the `17/24` theorem.

## What this does and does not close

The theorem covers height-25 boundary labels only. It does not cover
labels of height 26 or more, changes to the 92-word lower core, or a
different ambient tensor.

The tested escapes outside the theorem remain above the threshold:

- q97 with the height-26 word `46`: cutoff-47 fit `0.7098643`;
- q99 with the height-39 chain `160,192,224`: cutoff-47 fit
  `0.7098582`, complete numerical pressure `0.7099930431`;
- q94 dropping mode 11 (leaving the old protected core): cutoff-47 fit
  `0.7099289`;
- q94 plus the q99 chain: cutoff-47 fit `0.7099123`;
- the q45 minimal deletion closure: complete `0.7156911630`.

Supports through height 38 are still covered by the `29/41` floor.
Height 39 is the only remaining protected-core escape; the Jensen
condition admits the narrow near-threshold band around `0.7061--0.7071`,
but no split in that band is constructed. The next useful step is a
different protected low block, a different ambient tensor, or an actual
height-39 split -- not another boundary-image, retraction, or word
addition search on the fixed q96 core.

## Scope

The `17/24` theorem is the universal statement. The cutoff-47 order hash
and pressure values are a finite control: they reproduce the saved
order lists and fits, not a fresh complete profile. The pressure is a
numerical fit; no exact pressure certificate, uniform physical type,
transfer, or score change is claimed. The complete factoring exponent
remains `1/2`.

## Reproduction

```sh
make verify-hot-boundary-floor
make verify-f8-high-image-order-invariance
```

The second target replays the saved cutoff-47 order lists, source pins,
fits, and hashes. The profile run that produced them used the scripts in
`source_sha256` and took about seven minutes per image on the recorded
host; the replay does not repeat that computation.
