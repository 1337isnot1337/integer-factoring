# Principal paired alphabets on the retained five-letter pencil

The completed [free-height asymmetric search](two-place-free-gauge-asymmetric-search.md)
uses bounded normalized heights on all five physical letters. Sending
unwanted letters to infinite height is outside that finite search box.
An exact principal restriction instead removes those letters before any
optimization. This supplies another constructive boundary search, not a
global exclusion or a new factoring exponent.

Choose any two subsets `S_X,S_Y` of the five original letters, each of
size at least two. Restrict the four legs to `S_X,S_Y,S_Y,S_X`. Every
retained source has literal endpoint `[x=v][y=u]`, so its restriction has
the same strict rectangular endpoint. Both cut matrices are square of
order `|S_X| |S_Y|` and invertible at the endpoint. Their generic full
rank is therefore proved, independently of any finite series window.

The [constructor](../tools/pencil_principal_alphabets.py) keeps every
original selected polynomial entry. Natural heights are restricted on
each leg and shifted only by that leg's constant reference value. Both
restored physical energy vectors are recomputed. This operation does not
assume that the old natural-profile classes remain equivalent at new
heights or after restriction.

The [checker](../tools/check_pencil_principal_alphabets.py) verifies 209
literal selected entries and endpoints for three unequal restrictions,
including nonbinary sources. Fifty complete orders at new integer gauges
agree with separate untruncated, fraction-free Smith elimination. The
[receipt](../evidence/pencil-principal-alphabets-controls.json) retains
all tested heights, orders and source hashes.

```sh
python3 tools/check_pencil_principal_alphabets.py
```

## Registration and numerical boundary

There are 26 subsets of size at least two, with total cardinality 75.
The [manifest](../evidence/pencil-principal-alphabets-manifest.json) fixes
all `67*26^2=45,292` matched-pair restrictions, including the 67 full
controls, before evaluation. A completed screen must resolve
`2*67*75^2=753,750` individual channel orders. The 67 source points are
the previously selected representatives, not all `4^11` pencil members.

The [driver](../tools/screen_pencil_principal_alphabets.py) independently
queries every restricted tensor. It retains an append-only profile ledger,
source texts, literal-input hashes and scalar-fit records. Missing native
pivots leave it incomplete. Numerical scalar fits reuse only exactly equal
inherited curves. To check that equivalence, center each profile at its
minimum, each physical energy list at its minimum, and retain the exact
gap `min(G_X)+min(G_Y)-min(profile0)-min(profile1)`. These data determine
the entire scalar pressure equation. They say nothing about new gauges.

The scalar fitting procedure is the existing 25-point grid with local
refinement around its three best points. It does not certify a global
minimum. Its zero-entropy guard is inherited from that exploratory tool;
no numerical fit is an exact positive or negative pressure certificate.

The registered follow-on selects up to twenty distinct proper restricted
tensors by their inherited balanced rate, then tries the existing seven
free-height asymmetric starts on each. Trial and global ceilings are
explicit. Successful local optimizer returns do not close unsearched
heights, angles, source parameters or other tensors.

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/screen_pencil_principal_alphabets.py \
  --manifest evidence/pencil-principal-alphabets-manifest.json \
  --directory /tmp/pencil-principal-screen
```

No completed screen or improved score is asserted at registration.
