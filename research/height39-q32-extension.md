# A height-39 extension of the height-16 core stays on the plateau

## Exact construction

Start from the binary q32 coherent winner, which retains the 28
coordinate words of height below 16 and a rank-four retract of the five
height-16 words.  Add the four one-hot coordinates

```math
96=(5,6),\quad 160=(5,7),\quad 192=(6,7),\quad 224=(5,6,7),
```

of heights `24,26,28,39`, to both the embedding and the retraction.
They are outside the old span, and the old maps are orthogonal to them,
so `LP=I` still holds.  The map is deletion closed: the derivatives of
`96,160,192` lie in the height-below-16 core, and deleting a mode from
`224` gives `96,160,192`, which are now retained.  The result is a valid
q36 map with maximum height 39.

Both complete natural cuts resolve all `36^2=1296` orders, with no
unresolved tail.  The numerical natural pressures are

```math
\rho_{32}=0.7121150125007152,\qquad
\rho_{36}=0.7121133074903934.
```

The height-39 extension gains `1.705e-6` and remains above
`1/sqrt2` by about `0.00502`.

## Interpretation

The minimal height-39 extension of the q96/F8 hyperoval profile gained
only `1.27e-5` and reached `0.709993...`; the present q36 extension of
the much smaller height-16 core gains only `1.7e-6` and remains at
`0.712113...`.  Adding high words to either a small or a large protected
core leaves the pressure near the core's plateau.  The bottleneck is the
protected low block and its profile, not the availability of support
heights 26--39.

This is a finite negative control.  It does not prove a floor, construct
a split, or change the factoring exponent.  The q32 and q99 maps differ
in field, core size, and boundary, so the numerical comparison is
illustrative; the common conclusion is that the tested height-39
boundary additions are far too small.

## Reproduction

```sh
python3 tools/check_height39_q32_extension.py \
  --output evidence/height39-q32-extension.json
make verify-height39-q32-extension
```

The [receipt](../evidence/height39-q32-extension.json) records both
complete profiles' pressure fits, the extension words and heights, and
the exact `q=36` and `1296` order counts.
