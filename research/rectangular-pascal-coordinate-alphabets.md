# Independent coordinate alphabets on the two physical pairs

## Construction and scope

Choose any coordinate sets `S_B,S_X` in the same ambient Pascal Fock
space, of sizes `A,B`. The original polynomial tensor is restricted to
`S_B,S_X,S_X,S_B` on its `(b,x,v,y)` legs. In particular, a fixed
occupation or parity restriction on one matched pair need not be imposed
on the other. Neither set needs deletion closure.

The [constructor](../tools/rectangular_pascal_sectors.py) evaluates the
original minors of `D`, over F4, rather than replacing them by one-particle
profiles. A bounded `2^20`-entry minor cache can be shared between coordinate
restrictions of the same source tensor. Cache eviction affects only cost.
The original equal-alphabet constructors are unchanged.

This extends the [coordinate Pascal implementation](pascal-coordinate-sectors.md)
and uses the existing [independent paired-type transfer](paired-type-crossing-pressure.md).
It supplies a constructive search space, not a better factoring exponent.

## Endpoint and complete occupation blocks

At `z=1`, the original tensor is

```math
 T_1(b,x,v,y)=[b=y]\det P[v,x].
```

Every nonzero off-diagonal entry of the Pascal compound strictly increases
the occupied-mode-index sum from `x` to `v`. On a tensor power, restricting
both `x,v` words to the same type makes their total sums equal. Hence all
off-diagonal transitions vanish. Independent matching types on `b,y`
preserve their delta factor. This argument applies to any two principal
coordinate restrictions and gives the required crossing endpoint.

Occupation conservation gives square blocks in both cuts:

```text
zero:     |b|-|v|=|y|-|x|,  row b+A*v, column x+B*y;
infinity: |b|+|x|=|v|+|y|,  row b+A*x, column v+B*y.
```

Indices in these formulas refer to the separate coordinate lists, not the
ambient word values. Both total matrix orders are `AB`. The endpoint
matrix in each block consists of copies of a principal compound Pascal
matrix, with a permutation of columns. Each principal submatrix is unit
triangular when ordered by the increasing mode-index sum. Thus every block
is generically invertible. This supplies the exact full-rank promise to
the native Smith engine; unresolved pivots raise an error.

All block orders contribute to each complete profile. Row and column
enumerations use different strides when `A != B`. The restored physical
energies are independent vectors of lengths `A` and `B`. For a candidate
rate `R=2 rho`, the numerical equation is

```math
 \log B=R\left(\log\sum_{b\in S_B}e^{-G_B(b)/R}
                   +\log\sum_{x\in S_X}e^{-G_X(x)/R}\right).
```

Both physical costs and the common valuation-scale derivative are retained.
The balanced threshold is `rho < 1/sqrt(2)`. The existing supporting-line
criterion can also be evaluated from these complete profiles, with target
value below one; numerical fits alone certify neither criterion.

## Independent validation

The [checker](../tools/check_rectangular_pascal_sectors.py) verifies 514
selected original entries and endpoints using a separate subset dynamic
programme. It compares 76 complete nonbinary orders with untruncated,
fraction-free Smith elimination, across three unequal-alphabet examples
and both cuts. Twenty-six finite differences test arbitrary symbol
derivatives and the valuation scale. Equal-alphabet specialization agrees
with the preceding F4 backend at 132 complete orders and every derivative
component. The [receipt](../evidence/rectangular-pascal-coordinate-controls.json)
includes source hashes and all tested alphabet lists.

```sh
python3 tools/check_rectangular_pascal_sectors.py
```

These finite controls do not reconstruct every entry of the larger screen.

## Registered larger screen

The [manifest](../evidence/rectangular-pascal-coordinate-manifest.json)
registers nineteen eight-mode cases before evaluating them. They include
unequal alphabet sizes, single-pair occupation restrictions, single-pair
parity restrictions, swapped pairs, and an equal-32 control. Filled-sea
excitation order is deterministic, with word-value tie breaks.

The [driver](../tools/screen_rectangular_pascal_sectors.py) retains every
complete profile, all numerical pressure fits, unresolved-precision
attempts and full as-run sources. The precision policy is 512, then 1024,
then 2048 coefficients; a wholly unresolved case leaves the screen
incomplete. After the screen, the best balanced case with distinct
paired alphabets is selected for the registered one-hour whole-symbol
follow-on. No result from that follow-on is asserted at registration.

```sh
python3 tools/screen_rectangular_pascal_sectors.py \
  --manifest evidence/rectangular-pascal-coordinate-manifest.json \
  --output /tmp/rectangular-coordinate-screen.json
```

The [whole-symbol driver](../tools/search_rectangular_pascal_sectors.py)
implements that registered selection directly from a completed screen.
It replays the chosen complete baseline, checks five new derivatives,
removes only the eight constant-height gauge directions, and records
successful and rejected queries separately. Its budget counts actual
search time; complete-profile and derivative initialization are recorded
separately. Twelve restarts and one thousand iterations per restart are
ceilings, not promises that every restart completes within the hour.

```sh
python3 tools/search_rectangular_pascal_sectors.py \
  --screen /tmp/rectangular-coordinate-screen.json \
  --output /tmp/rectangular-coordinate-search.json
```

A [separate two-by-three smoke receipt](../evidence/rectangular-pascal-driver-smoke.json.xz)
checks the screen/selection/search pipeline with an explicitly shortened
one-second budget and two-restart ceiling. It is a driver control, not part
of the nineteen-case screen or its one-hour follow-on.

## Completed nineteen-case screen

The frozen `b04ca4c` screen completed all nineteen cases in 466.92 seconds.
The [full archive](../evidence/rectangular-pascal-coordinate-screen-complete.json.xz)
retains **34,304 complete orders**, both physical type vectors, all fitted
parameters, stdout and as-run source texts. All 57 supporting-line fits
reported local success. These status flags are numerical, not proofs of
global optimality. The [coverage receipt](../evidence/rectangular-pascal-coordinate-screen-controls.json)
checks dimensions, registration agreement and source-text hashes; it does
not independently reconstruct all larger native profiles.

| Restriction, with both pair orientations tested | Best inherited numerical rho |
| --- | ---: |
| Equal 32-letter control | 0.7216933250053 |
| 24 by 48 lowest-excitation words | 0.7218724642968 |
| 16 by 64 lowest-excitation words | 0.7249020226556 |
| 8 by 64 lowest-excitation words | 0.7356440245148 |
| 32 even words by 32 odd words | 0.7867432427894 |
| 24 even words by 48 unrestricted words | 0.8186209971961 |
| 24 odd words by 48 unrestricted words | 0.8201486375873 |
| 32 occupation-four words by 32 unrestricted words | 0.8578975618474 |
| 16 occupation-four words by 32 unrestricted words | 0.8613599937479 |
| 16 occupation-three words by 32 unrestricted words | 0.9031303929582 |

The lowest supporting-line value is `1.0001056729915`, still above one.
The registered rule selects `full-24-by-48-swapped`, of sizes 48 and 24,
for the one-hour whole-symbol follow-on. Its last-bit advantage over the
unswapped case is numerical roundoff, not evidence of a structural
orientation advantage. The follow-on is running separately from this
completed screen. No exact stronger crossing or factoring bound follows.

Full default verification plus six named focused gates passed in a clean
checkout at exact `7dd1903`; the [receipt](../evidence/verify-7dd1903.json)
and [compressed log](../evidence/verify-7dd1903.stdout.txt.xz) retain its
source tree and command. That pin includes the rectangular constructor
and cache controls, but predates the completed-screen archive and bundle
adapter. Their focused checks and source pins are separate.
