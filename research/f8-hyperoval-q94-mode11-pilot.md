# A q94 mode-11 drop pilot is slightly below q99

## Construction

Start from the complete q96 F8 hyperoval profile. Delete the two letter
coordinates whose supports use mode 11, at indices `90` and `91`. The
remaining 94 letters are unchanged. Because the surviving map is a
principal submap of a map with `LP=I` and deletion closure, both
properties persist; the checker revalidates them directly on 11 modes.

The two deleted letters are the only parts of the map that use mode 11.
The resulting map has 94 letters, maximum height 25, and no height-39
boundary.

## Cutoff-47 pilot

A `partial_profiles` run at cutoff 47 retains 932 orders in each cut and
leaves 7,904 orders per cut unresolved. The numerical natural pressure
is

```math
\rho \approx 0.709928941871497,
```

at scale `0.247919041256729`. This reproduces the pilot quoted in
[the minimal-core note](f8-hyperoval-minimal-core.md) and is below the
complete q99 value `0.7099930431440509` by `6.41e-5`.

The value is still above `1/sqrt(2)` by about `0.002822`, and the pilot
leaves about 90 percent of the orders unresolved. It is a lower bound,
not a complete profile or a score claim: the map has a height-zero
letter and all orders are nonnegative, so adding the unresolved orders
can only increase the pressure (see the
[truncated-profile lower bound](truncated-profile-lower-bound.md)).
Therefore the complete q94 pressure is at least `0.709928941871497`,
and the q94 drop **cannot cross**; a full q94 build is unnecessary for
the crossing question.

The pilot being below q99 is an artifact of the omitted orders, not
evidence of a better complete low block. It does establish that the
11-mode drop is valid and gives a reproducible lower bound; it does not
establish a complete value, a certificate, or a score change. A
different protected low block remains the structural requirement for a
crossing.

## Reproduction

```sh
make verify-f8-q94-mode11-pilot
```

The checker rebuilds the 94-letter map from the q96 profile, validates
the map identity and deletion closure, and replays the recorded pilot
pressure and its defining equation. It does not regenerate the profile.
The [receipt](../evidence/f8-hyperoval-q94-mode11-pilot.json) contains
the recorded orders, heights, charges, tail counts, base-profile hash,
and source pins.

The full q94 profile is not computed here; at cutoff 47 the pilot took
about six minutes in this checkout.
