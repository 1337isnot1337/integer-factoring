# Two-even-site codes with different output sections

This search changes the logical code and the output retract. It therefore
tests constructions beyond tensoring the fixed 24-label coefficient space.
It finds no strict finite seed. The conclusions below concern the recorded
finite coefficient spaces over the stated fields; they do not exclude
larger two-even-site codes or all tensor powers of the modified spaces.

For the specific five-site code with the upper output section, the
later [common-functional compression theorem](common-functional-rank-compressions.md)
does exclude every tensor power over every characteristic-two field.
The modified four-site sections below retain their stated finite scope.

## Construction

Write `D_i` for the occupation of both modes at physical site `i`.
Choose two distinct sites `a,b`. At every other site choose one of its two
single modes, and let `s` denote their product. Use the input vector

\[
I_{\{a,b\},s}=s(D_a+D_b).
\]

An output functional may select either of the two coordinates `s D_a`
and `s D_b`. Different labels have disjoint input support, so either
choice gives a dual retract. Each input has self-pairing zero in
characteristic two, and the inputs are mutually orthogonal. All labels
have degree `P`, where `P` is the number of sites.

Keeping every pair and every single-mode choice gives

\[
d=\binom P2 2^{P-2}.
\]

This is `d=24` at four sites and `d=80` at five sites. The literal Pfaffian
entry calculation checks the crossing endpoint and the zero first
noncrossing endpoint. The five-site construction is the part of the
middle-degree symplectic radical with exactly three single-occupied sites.

## Complete four-site graph and section search

For each of the six edges of the complete graph on four sites, the search
chooses one of three states:

1. Omit that pair of even sites.
2. Keep the pair and select its lower-numbered site in every output.
3. Keep the pair and select its higher-numbered site in every output.

All four single-mode orientations are retained for each chosen edge.
Thus there are exactly `3^6-1 = 728` nonempty cases. The search includes
every labeled graph, without assuming that the chosen output sections
respect graph isomorphisms. It does not include output choices that vary
with the single-mode orientation within an edge.

The program builds one literal rectangular tensor containing both output
choices, then restricts its two flattenings to each case. Three separate
literal constructions verify the restriction interface. The coefficient
spaces have the following dimensions:

| Coefficient dimension | Cases |
| ---: | ---: |
| 3 | 24 |
| 5 | 62 |
| 6 | 642 |

Across the 728 spaces, every one of the **42,536 nonzero binary coefficient
vectors** is checked. Exactly 3,942 vectors satisfy the deliberately
generous single-piece cost gate. In every case, these affordable vectors
span a subspace that does not contain the crossing tensor. Consequently,
every binary decomposition in every recorded space fails the strict
finite-seed criterion.

The independent Python replay computes both exact cut ranks for all
42,536 vectors, including every rejected vector. It also rechecks each
affordable-span calculation and the exact crossing endpoint. The native
screens are retained in a single JSON archive.

## Five-site extension-field search

The five-site code has 80 labels and coefficient dimension eight.
The complete finite searches give:

| Field | Nonzero vectors or projective points tested | Affordable points | Affordable span dimension | Crossing in span? |
| --- | ---: | ---: | ---: | --- |
| `F2` | 255 | 6 | 5 | No |
| `F4` | 21,845 | 8 | 5 | No |

Projective normalization is valid because multiplying a tensor by a
nonzero field scalar does not change either cut rank. The gate is 492:
it satisfies `492^70 > 80^99` and `99^2 > 2*70^2`, so every piece in a
decomposition with total cost below `80^sqrt(2)` must pass it. Because
the affordable span misses the crossing, neither field gives a strict
seed in this coefficient space.

The extension-field replay checks all eight affordable points and sixteen
deterministically selected additional points by independent scalar
`F4` elimination. The complete native enumeration provides the absence
claim for the other projective points. The binary affordable ranks are
included among these extension-field checks, where they agree with the
binary ranks.

## Reproduction

The archived checker regenerates the tensor bases from the construction;
no temporary basis file is an input dependency.

```sh
python tools/check_two_even_site_sections.py \
  --work-dir /tmp/two-even-site-replay \
  --output /tmp/two-even-site-replay.json
```

Use `--scope graphs` or `--scope five-site` to replay one part. The
graph explorer can regenerate the native screens with
`tools/explore_two_even_site_sections.py --native PATH --work-dir DIR
--output REPORT`. Its optional `--fixed-upper-only` mode selects the
63 nonempty labeled graphs with the upper output section.

- [Construction and graph explorer](../tools/explore_two_even_site_sections.py).
- [Independent replay](../tools/check_two_even_site_sections.py).
- [All graph and section records](../evidence/two-even-site-sections.json).
- [Archived native graph screens](../evidence/two-even-site-section-screens.json).
- [Five-site binary screen](../evidence/two-even-site-P5-F2.jsonl).
- [Five-site extension-field screen](../evidence/two-even-site-P5-F4.jsonl).
- [Combined replay certificate](../evidence/two-even-site-sections-audit.json).
