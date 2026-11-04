# The identity-core cutoff trend stops above the height-38 floor

## Question

The protected-core maps used so far keep a low core and add a deletion-closed
high boundary.  The simplest baseline keeps only the identity core
`C_H = {words of height < H}` and restricts the two cuts to it.  Does raising
`H` lower this baseline pressure enough to cross `1/sqrt(2)`?

## Construction

For each cutoff `H`, take the one-hot coordinates of all physical words with
`height(word) < H`, with the natural mode heights and the common balanced cut
scale.  The map is the identity on those coordinates, so `LP=I` and deletion
closure hold.  The two cut profiles are the restricted tensor profiles; the
pressure is the numerical natural two-cut pressure.

## Measurements

| `H` | `q` | pressure | `scale_t` | wall seconds |
|---:|---:|---:|---:|---:|
| 16 | 28 | 0.7124841635989507 | 0.3407915124178001 | 3.1 |
| 18 | 38 | 0.7116418233895179 | 0.31469232290926774 | 23.7 |
| 20 | 49 | 0.7110773177567727 | 0.29566496532486897 | 215.4 |
| 22 | 64 | 0.7105903583541862 | 0.2780334985986807 | 758.3 |
| 24 | 81 | 0.7102275227039717 | 0.26396351552252123 | 1559.1 |
| 25 | 92 | 0.7100552807213277 | 0.2569742280260311 | 3712.9 |

The trend decreases monotonically but flattens.  The successive drops are
about `8.42e-4`, `5.65e-4`, `4.87e-4`, `3.63e-4`, and `1.72e-4`.  The
q96 full map has pressure `0.7100058347`; it is `4.9e-5` below the `H=25`
identity core and `2.2e-4` below the `H=24` core, so the q96 boundary helps
the 81-word core but does not change the scale of the plateau.

The q94 mode-drop pilot reported in
[the minimal-core note](f8-hyperoval-minimal-core.md) has numerical fit
`0.709929`; it is below q96 and q99, and is the best recorded support-at-most-25
pilot.  It is not a complete profile.  The q99 height-39 map is
[`0.7099930431`](f8-hyperoval-height39-extension.md).

## Why the identity-core route cannot cross

The [height-38 protected-core floor](jensen-tail-protected-core-floor.md)
applies to every constant height-homogeneous paired map that retains the
height-below-24 core and is supported through height 38.  For every `H <= 38`,
the identity core `C_H` retains that low core and has support below 38, so

```math
\rho(C_H) \ge 29/41 = 0.7073170731707317.
```

The measured `H=25` value is `0.7100553`, about `2.74e-3` above the floor and
`2.95e-3` above `1/sqrt(2)`.  Adding core words down to height 38 therefore
cannot cross: the floor is a proved barrier, not a numerical extrapolation.
The q96 and q94 maps are also support-at-most-25, so the floor contains them.
The q99 map escapes the floor only because it uses a height-39 coordinate;
its pressure is still `0.7099930`.

The [q45 minimal-core note](f8-hyperoval-minimal-core.md) shows the opposite
failure mode: deleting the q96 core words unused by the boundary image raises
the pressure to `0.7156912`.  The core words are load-bearing, but the
height-38 floor still prevents the enlarged identity cores from crossing.

## Consequence

Larger identity cores and additional height-24/25 boundary words are not a
route below `1/sqrt(2)`.  The next change must alter the low block or its
gauge while keeping a height-39-supported boundary.  A natural bounded test is
the q94 mode-drop low block plus the q99 height-39 coordinates; it remains a
finite map optimization, not an escape from the low-core bottleneck.

## Reproduce

```sh
make verify-core-cutoff-trend          # fast replay: H = 16,18,20
python3 tools/measure_core_cutoff_trend.py   # full table; about 1.6 hours
```

The complete receipt is
[`evidence/core-cutoff-pressure-trend.json`](../evidence/core-cutoff-pressure-trend.json).
The pressures are numerical fits, not exact rational certificates; the
complete factoring exponent remains `1/2`.
