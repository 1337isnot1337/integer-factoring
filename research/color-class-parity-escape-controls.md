# Parity escapes in the color-class cube lose the endpoint or the spectrum

## Question and scope

The [sector obstruction](color-distinct-sector-obstruction.md) assumes that
an entry vanishes whenever an individual color occurs an odd number of
times. The one-hot color-class tensor has that property. A parity-breaking
variant would evade the hypothesis. This note records a finite `q=4`
control of the natural variants: support forgetfulness, grouping colors,
an endpoint-restoring antipodal delta, pair-weighted support, and a
color-averaged mixture of one-hot tensors.

The test is structural, not a new pressure theorem. The individual
antipodal endpoint and the finite Smith spectrum are checked at `q=4`
through order six. No uniform statement or score change is claimed.

## Variants and checks

All variants use the same eight-mode Gaussian cube, the same `q=4`
alphabet, and evaluation at `(y,z)=(2,3)` over `F_256`.

- **One-hot control.** `T(a)=prod_c G_{S_c}`; individual endpoint correct.
- **Support forgetful.** `T(a)=G_{S(a)}`, where `S(a)` is the set of
  nonvacuum vertices. Colors are forgotten.
- **Two groups.** `T(a)=G_{S_1}G_{S_2}` with colors `{1,2}` and `{3}`.
- **Support diagonal.** Multiply the support tensor by the antipodal
  delta; the cut matrix becomes diagonal.
- **Group diagonal.** The same delta applied to the two-group tensor.
- **Pair-weighted support.** Multiply the support tensor by
  `lambda` to the number of mismatched antipodal pairs, so the endpoint
  is unchanged; `lambda=2`.
- **Color average.** Sum the one-hot tensors over all translations of the
  `F_2^2` color labels, with weights normalized to sum one. The endpoint
  is preserved and individual parity is broken.

The checker records the low finite Smith orders through six, the number of
zero orders, the number of finite orders, and a grid pressure objective
using only those finite orders. For a fixed finite-order set, adding the
missing orders can only increase the objective, so the listed value is a
lower bound for the full finite spectrum. The one-hot control is evaluated
on its full committed `q=4` profile.

| variant | individual endpoint | finite orders through 6 | zero orders | lower-bound objective |
|---|---|---:|---:|---:|
| one-hot | yes | 20 | 1 | `0.541` (`0.667` full) |
| support forgetful | **no** (9744 spurious) | 4 | 1 | `0.252` |
| two groups | **no** (1040 spurious) | 10 | 1 | `0.417` |
| support diagonal | yes | 256 | 136 | `1.004` |
| group diagonal | yes | 256 | 72 | `1.003` |
| pair-weighted support, `lambda=2` | yes | 256 | 136 | `1.004` |
| color average | yes | 72 | 4 | `0.773` |

The spurious counts are full eight-vertex colorings that are
group-antipodal but not individually antipodal, out of `4^8=65536`.
The support-forgetful count is `9744`; the two-group count is `1040`.
Those variants are not endpoint tensors and their favorable low spectra
are not candidates.

Every endpoint-correct variant is worse under the finite-spectrum
objective. The diagonal restorations collapse the cut matrix and create
72 or 136 zero orders. The pair weighting keeps the endpoint but leaves
136 zero orders. The color average has only 72 finite low orders and a
lower-bound objective `0.773`, already above the one-hot full value
`0.667`.

## Conclusion

The natural parity escapes split into two failures. Grouping or forgetting
colors lowers the low orders but destroys the individual endpoint.
Restoring that endpoint with an antipodal delta or a pair weight collapses
the cut matrix to a diagonal or near-diagonal spectrum with many zero
orders. Averaging over color translations preserves the endpoint but adds
low orders without recovering the one-hot finite spectrum.

A viable parity-breaking construction must preserve both the individual
endpoint and a dense finite low-order spectrum. These controls do not
provide one. The one-hot profile remains the comparison point, and its
large-`q` obstruction is not removed by the tested escapes.

## Reproduction and credit

```sh
python3 tools/check_color_class_parity_escapes.py
make verify-color-class-parity-escapes
```

The checker uses the one-hot constructor and profile extraction from
[`check_color_class_cube_profile.py`](../tools/check_color_class_cube_profile.py),
the Gaussian kernel from the shared cube tools, and the committed
`q4_full_profile` receipt. The sector obstruction and its source credit
are in [color-distinct-sector-obstruction.md](color-distinct-sector-obstruction.md).
No Lean formalization or independent mathematical review is claimed.
