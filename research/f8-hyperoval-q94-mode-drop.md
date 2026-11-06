# Dropping the height-23 mode does not lower the complete q94 pressure

## Question

The q96 F8 hyperoval map retains the ten canonical height-25 boundary triples.
Mode 11 (height 23) is the only ambient mode that no boundary triple uses.
Removing it gives a q94 map on 11 modes.  A cutoff-47 pilot reported a
numerical fit `0.709929`, slightly below q99, but explicitly did not resolve
the tails.  Does the complete profile confirm the improvement?

## Construction

Start from the q96 map and delete every embedding and retraction coordinate
whose support uses mode 11.  The remaining map is the q94 map: it keeps all
q96 coordinates not using mode 11, including all ten canonical height-25
boundary triples, and is deletion closed.  The complete profile has `8836`
orders per cut at cutoff 135 and zero unresolved tails.

## Result

The complete natural pressure is

```math
\rho \approx 0.7100458192
```

with `scale_t = 0.2563463846`.  This is **above** the q96 value
`0.7100058347` by about `4.0e-5`, above the q99 value
`0.7099930431`, and above `1/sqrt(2)` by about `2.94e-3`.  The cutoff-47 pilot
`0.709929` did not survive the complete profile: its unresolved tails moved
the fit in the favorable direction.  The q94 map does not lower the pressure.

## Interpretation

The q94 map removes part of the full height-below-24 core, so the height-38
protected-core floor does not apply to it.  It nevertheless has a higher
pressure than q96 and q99: the removed height-23 words are load-bearing, in
the same direction as the [q45 minimal-core test](f8-hyperoval-minimal-core.md),
whose more aggressive deletion has pressure `0.7156912`.  Thus the
single-mode-drop low-block change is not a route.  The q96 and q99 low blocks
remain the best recorded maps; the q99 height-39 map is still the best
height-39 construction.

A q97 composition of the q94 low block with the q99 height-39 coordinates
would start from a worse base and is not pursued.  The next low-block change
must alter more than the highest unused mode, for example a different
supported-height set or ambient structure.

## Scope

This is a complete finite profile of the specified q94 map.  The pressure is
a numerical fit, not an exact rational certificate.  It is above `1/sqrt(2)`
and the complete factoring exponent remains `1/2`.  No uniform transfer or
global optimality is claimed.

## Reproduce

```sh
make verify-f8-q94-profile
make build-f8-q94-profile   # about 3.5 hours
```

The receipt is
[`evidence/f8-hyperoval-q94-profile.json`](../evidence/f8-hyperoval-q94-profile.json).
