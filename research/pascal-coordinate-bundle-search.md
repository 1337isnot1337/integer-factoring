# Fixed-type bundle steps on complete coordinate crossings

The [new driver](../tools/search_pascal_coordinate_bundle.py) adapts the
existing [rectangular core bundle method](rectangular-convex-pressure.md)
to the [independently paired coordinate constructor](rectangular-pascal-coordinate-alphabets.md).
It retains the complete original-polynomial profiles. This is a different
numerical optimizer for the same crossing objective, not a new transfer
theorem or a larger feasible family.

For fixed paired types, `log B + <pi_B,G_B> + <pi_X,G_X>` is convex in
the raw physical heights and positive valuation scale. Each Smith prefix
is a minimum of affine minor valuations. For any decreasing probability
vector, its weighted sum of ordered Smith values is a nonnegative linear
combination of these prefixes. The entropy variational formula for
`log(sum exp(-order))` then makes `log B` convex. The prior note gives the
equivalent full minor-constraint formulation.

At a queried profile, the attaining nested minors provide a supporting
plane for `log B`, including at valuation ties. The driver recovers this
subgradient from the established Gibbs-root derivative by removing its
current physical-type contribution. It accumulates the planes, solves
a proximal quadratic master, and accepts only candidates whose complete
profile evaluation gives sufficient actual improvement. Gibbs type updates
and restarts remain exploratory; joint optimization over types and heights
is not asserted to be convex.

Every new plane is tested against every earlier query in both directions.
The report retains the actual planes, query heights, best complete profiles,
master failures, profile rejections and separate construction/search times.
Neither the finite plane checks nor a numerical master stopping condition
proves a global optimum. The search has explicit raw-height, normalized-height
and scale bounds; it makes no claim outside them.

The [adapter checker](../tools/check_pascal_coordinate_bundle.py) controls
48 derivatives and 132 held-out supporting-plane inequalities on twelve
nonbinary unequal-alphabet queries. It also checks the fixed-type objective
identity. The [receipt](../evidence/pascal-coordinate-bundle-controls.json)
retains the queried planes and source hashes. A separate
[small driver control](../evidence/pascal-coordinate-bundle-smoke.json.xz)
checks three resolved masters and their actual null-step evaluations.

```sh
python3 tools/check_pascal_coordinate_bundle.py
```

The [registered comparison](../evidence/pascal-coordinate-bundle-manifest.json)
fixes two starts before evaluation: the equal 32-letter control and the
selected 48-by-24 coordinate restriction. Each receives thirty minutes
of search time, with three restarts and four hundred iterations per restart
as ceilings. This complements the separate one-hour direct search on the
selected rectangle. No stronger contraction is asserted at registration.

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/search_pascal_coordinate_bundle.py \
  --manifest evidence/pascal-coordinate-bundle-manifest.json \
  --output /tmp/pascal-coordinate-bundle.json
```

No new Lean formalization or improved complete factoring score is included.
