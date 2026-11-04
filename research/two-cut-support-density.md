# The Gaussian two-cut support is about one third full

## Result and scope

For the Gaussian reflection seed, the two cuts `(x,y)|(u,v)` and
`(x,u)|(y,v)` share the same four-index tensor. Their joint image is the
tensor support. For the checked seeds:

| Interval | Alphabet | Entries | Nonzero | Support fraction | Product/support |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `4` | `4` | `256` | `104` | `0.40625` | `2.4615` |
| `8` | `16` | `65,536` | `22,328` | `0.34070` | `2.9351` |

The support is not full, but the compression factor is small and slowly
increasing toward about `3`. It is a constant-factor saving, not an
exponent change. The Smith codec compresses each cut separately; the
two-cut product is still the product of the two codebooks. This note
does not claim a rate improvement.

## Relation to the local-solver route

The bundle and fan-out diagnostics show large local compressions when
the two codecs read the same input digits. The two-cut joint image of the
tensor is smaller than the raw product, but only by a factor below `3`
on these seeds. A local solver that exploits the raw support can
therefore save at most a constant factor here. A crossing needs the
Smith codec or a different source family, not the raw two-cut support.

The support fractions are consistent with a product structure: the
one-particle constraints remove a fixed fraction of the entries, and the
remaining fraction stabilizes. This is a finite observation, not a proof
of the limiting density.

## What remains

A codec-level joint image would be the relevant object if the Smith
transforms correlate the two cuts more strongly than the raw support.
Computing that image requires the explicit cut codecs, not just the
support. The next bounded diagnostic is to extract the q4 cut codecs
from the native Smith oracle and measure their joint image against the
product of the two code dimensions.

## Reproduction

```sh
python3 tools/check_two_cut_support_density.py \
  --output evidence/two-cut-support-density.json
make verify-two-cut-support-density
```

The checker counts the exact support of the interval-4 and interval-8
Gaussian seeds.
