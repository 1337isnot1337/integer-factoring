# Nonbinary flags with complete coordinate Fock profiles

The selected-entry occupation-sector constructor now extends to F4-valued
Pascal kernels and constant `GL(L,F4)` flags. It allows arbitrary potentials
on retained occupation words, including correlations between modes. Every
selected original polynomial entry and every promised Smith order is retained.
This is a constructive search extension, with exact finite controls and a
registered larger experiment. No improving contraction is claimed.

## Construction and endpoint

The [constructor](../tools/f4_pascal_coordinate_sectors.py) computes each
selected original minor by exact F4 polynomial arithmetic. It retains the
bounded sparse minor cache of the binary predecessor. Occupation conservation
still gives the same exact block decomposition at both cuts: a nonzero
minor has equal row and column cardinality, independently of its coefficients.

The [constant F4 flag construction](f4-pascal-flags.md) changes the two
diagonal one-particle blocks by inverse flags. Those blocks vanish at the
endpoint, so the original endpoint remains
`T_1(b,x,v,y)=[b=y] det P[v,x]`. The principal physical coordinate restriction
therefore retains the same [typed crossing proof](pascal-symbol-potentials.md):
nontrivial Pascal transitions strictly increase the mode-weight statistic
and disappear after imposing matching types. This does not require deletion
closure or a Gaussian profile formula after restriction.

The finite implementation accepts up to sixteen modes and sixty-four
distinct selected words. It preserves both independent matching-pair costs.
Historical binary constructors and source-pinned receipts remain unchanged.

## Independent checks

Run `make verify-f4-pascal-coordinate-sectors`. The
[checker](../tools/check_f4_pascal_coordinate_sectors.py) verifies 4,977
selected original entries and endpoints with an independent subset dynamic
program using packed field polynomials. Its recursive constructor instead
uses sparse minor expansion and two-plane arithmetic. The full endpoint
comparison includes triangular Pascal transitions, rather than incorrectly
assuming a literal crossing before type projection.

For nonbinary flags, 210 complete Smith orders with arbitrary integer
symbol gauges agree with the independent fraction-free polynomial oracle.
Another 210 binary-specialization orders agree with the old sector backend.
Fifteen finite-difference checks control symbol and valuation-scale gradients.
The [receipt](../evidence/f4-pascal-coordinate-controls.json) stores every
input and complete nonbinary profile. These finite checks are distinct
from a uniform proof or an exact improving pressure certificate.

The [search driver](../tools/search_f4_pascal_coordinate_flags.py) uses the
same Gibbs-root objective on complete original-polynomial profiles. It fixes
the eight constant-height gauge directions and checks its gradient before
each run. It retains its source texts, rejected queries, physical costs,
best complete profiles and separate construction/search times. A
[small smoke run](../evidence/f4-coordinate-search-smoke.json.xz) exhausted
its two-restart ceiling with 110 successful queries; its rate is weak.

## Larger registered experiment

The [manifest](../evidence/f4-pascal-coordinate-screen-manifest.json) fixes
eight eight-mode, thirty-two-letter starts before their evaluation: the
identity flag, six single nonbinary transvections, and a double transvection.
The alphabet consists of the deterministic lowest-cost excitations around
the filled half interval, as in the existing
[binary coordinate starts](pascal-coordinate-sectors.md). Frobenius-conjugate
single coefficients are not duplicated.

The declared follow-on policy selects the lowest inherited rate among the
seven nonidentity flags and keeps the identity control. Each gets a one-hour
whole-symbol search budget, with recorded restart/iteration ceilings and
fixed seeds. Registration is not completed coverage; a time-budget stop
does not claim all restarts were exhausted. Any candidate below the critical
threshold still needs independent exact profiles and outward pressure before
it can support a stronger crossing or factoring claim.

Substantial OpenAI assistant support was used. No new Lean formalization or
improved complete factoring score is included; that score remains **1/2**.
