# A complete height-39 extension of the F8 hyperoval profile

## Exact construction

Start from the complete q96 F8 hyperoval profile. Add the three
single-word coordinates

```math
160=(5,7),\qquad 192=(6,7),\qquad 224=(5,6,7),
```

of heights `26,28,39`, to both the embedding and the retraction. The
three words are outside the old span, and the old retraction and
embedding vectors are orthogonal to them, so the enlarged map still has
`LP=I`. The map is deletion closed: deleting a new pair gives a
single-mode core word, and deleting a mode from `(5,6,7)` gives the two
new pairs together with the height-24 word `(5,6)`, which is already in
the q96 core. The 81-word height-below-24 core is retained unchanged.
The result is a valid q99 map with maximum height 39.

The complete natural profile resolves at cutoff 191: all `9801` orders
in each cut are retained, with no unresolved tails. The numerical
natural pressure is

```math
\rho \approx 0.7099930431440509,
```

compared with `0.7100057871974` for the q96 profile. The extension gains
about `1.27e-5`, and it is still above `1/sqrt(2)` by about `0.002886`.
The height-39 barrier no longer excludes this map, but the map does not
cross.

## Interpretation

This is a complete protected-core construction with support height 39,
and the gain over the height-25 q96 profile is only in the last decimal. The low core still dominates the pressure. Adding a
minimal deletion-closed high boundary lowers the pressure slightly, but
it does not close the `0.0029` gap to the threshold. The next useful
change is a different protected low block, not another boundary word.

The construction also confirms that the height-39 transition is not an
artifact of the Jensen-tail test: a valid map at that support height
exists, and its pressure sits near the q96 plateau rather than near the
height-39 floor `70/99 = 0.707070...`.

## Scope

The profile is a complete finite computation for the specified q99 map.
The pressure is a numerical fit, not an exact rational certificate; the
complete factoring exponent remains `1/2`. The map has the fixed F8
hyperoval image, the retained height-below-24 core, and the common
balanced natural cut scale. No uniform physical type, transfer, or
global optimality is claimed. Other height-39 triples or larger
extensions are not screened.

The full build takes about `3.7` hours and is reproducible with
`tools/build_f8_hyperoval_height39_extension.py`. The fast checker
`tools/check_f8_hyperoval_height39_extension.py` replays the map
identity, deletion closure, recorded orders, and numerical pressure; it
does not regenerate the full profile.

```sh
make verify-f8-height39-extension
```

The [receipt](../evidence/f8-hyperoval-height39-extension.json) contains
the complete orders, heights, charges, base-profile hash, source pins,
and numerical pressure.
