# Adaptive regular LDU chambers for complete profile queries

A fixed rational-series diagonalization can answer many weighted Smith
queries on the same source matrix. A cache hit checks regularity of every
row and column operation, then sorts the complete diagonal. It also permits
resuming a certified prefix when a later operation needs to change.

The [implementation](../tools/native_smith_chambers_f4.py) is optional.
Its finite controls agree with independent untruncated profiles, but current
large moving-trajectory timings are mixed. Repeated-point speed does not
establish a faster search, and no improving crossing or factoring exponent
is supplied. Historical native sources are unchanged.

## 1. A regularity cone, rather than a pivot-order cone

Let `A(z)` be a fixed Laurent-polynomial matrix of **known exact generic
rank** `r`. Take a resolved rational-series Gaussian path on `A`. Its
remaining entries at a pivot are `a_ij`, with pivot order `d`. Normalize
the pivot row by the inverse of its unit, making the pivot exactly `z^d`.
The coefficient clearing row `i` has raw order
`ord(a_ip)-d`; the coefficient clearing column `j` has order
`ord(a_pj)-d` after the unit normalization.

For the weighted source `D_R(z) A(z^lambda) D_C(z)`, `lambda>0`, the two
operations are regular and regularly invertible whenever

```math
R_i-R_p+\lambda\{ord(a_{ip})-d\}\ge0,
\qquad
C_j-C_p+\lambda\{ord(a_{pj})-d\}\ge0.                  \tag{1}
```

Unit normalization and row/column permutations are always regular units.
Once row clearing has made the pivot column zero below the pivot, right
clearing leaves the remaining lower-right block unchanged. Thus ordinary
Gaussian residuals suffice to obtain the conceptual two-sided diagonalization.

If (1) holds at every step, the resulting diagonalization uses regular
invertible factors. Its complete finite profile is the **sorted multiset**

```math
\{\lambda d_k+R_{i_k}+C_{j_k}:1\le k\le r\}.            \tag{2}
```

This is valid even when the diagonal orders change order. Requiring the
same globally smallest pivot at every step would impose unnecessary
conditions on the remaining interior entries. Only the row and column
clearing operations need to remain regular. There are at most `n(n-1)`
such conditions for a full-rank `n`-square matrix.

Rational weights follow by ramification; every finite input double is an
exact binary rational. A supplied rank must be exact: asking for the first
few pivots of a higher-rank matrix is outside this cache contract. When
`r<n`, the remaining mathematical block vanishes by that rank promise,
not because a truncated calculation happened to return zero.

## 2. Explicit knowledge of every partial series

The new backend tracks a finite series together with an integer `P` meaning
that all coefficients below absolute exponent `P` are known. An exact
finite polynomial has `P=infinity`. A zero stored prefix at finite `P`
means only that the actual order is at least `P`.

Addition retains knowledge through `min(P_a,P_b)`. Multiplication retains
knowledge through

```math
\min\{P_a+v_b,\ P_b+v_a\},                            \tag{3}
```

where each `v` is a known order or a proved order lower bound. Word-window
truncation can further lower the recorded precision. Unit inversion uses
only the available coefficients and records the resulting precision.
Unknown zero prefixes participate in subsequent error propagation; only
exact zeros can be skipped.

Two algebraic identities retain stronger information than a generic
interval calculation: dividing a pivot by its **actual** unit makes that
pivot an exact monomial, and subtracting its exact elimination multiple
makes the cleared entry exactly zero. Other products retain the conservative
knowledge (3). Reference-weight pruning explicitly records every discarded
coefficient as unknown; it never treats discarded tails as zero.

If a multiplier has only an order lower bound `P`, using `P-d` in (1)
still gives a sufficient regularity test. A successful complete construction
therefore yields an actual conditional regular diagonalization. If a promised
pivot cannot be resolved, or the resulting cone does not certify the current
point, the wrapper invokes the existing complete-profile oracle. An unresolved
legacy query still raises an error. No unknown tail becomes a rank claim.

## 3. Exact comparisons and prefix reuse

Cone tests compare an integer multiple of the positive scale plus differences
of input doubles. A separated-sign long-double filter handles ordinary
cases; close cases use error-free `TwoSum` expansions and an `FMA` product
residual. This checks the sign for the actual binary rational inputs rather
than accepting a negative condition within an arbitrary epsilon. The checker
also compares the signs with independent Python rational arithmetic.

The reference path uses row/column local minima, found by alternating scans.
Each strict move lowers the selected weighted entry, so the finite procedure
terminates. A globally smallest remaining entry is unnecessary for (1).

Four chambers per place are retained. A bounded residual checkpoint stores
the matrix, row/column labels and coefficient knowledge after a prefix.
When every condition before that checkpoint still holds, those earlier
regular transformations can be reused and only the remaining block rebuilt.
The saved partial series represents the same source Schur complement;
its truth does not depend on the next query's weights. Any further pruning
again records its coefficient loss explicitly.

Prefix reuse does not mean its diagonal entries remain the smallest ones.
The final sort includes both reused and newly computed pivots. The checker
forces an example in which a reused leading pivot moves to the last profile
position. It independently reconstructs the resulting complete profile.

Returned row and column IDs label the affine weights in (2), which gives
the correct derivative inside a chamber. After diagonal reordering they
are **not advertised as source-minor certificates**. Use the untruncated
minor oracle when an explicit attaining source minor is needed.

## 4. Exact and numerical controls

Run `make verify-smith-chambers-f4`. The
[checker](../tools/check_smith_chambers_f4.py) verifies:

- 512 signs against exact binary-rational arithmetic, including close
  cancellation and widely varying scales;
- 1,344 complete integer/half-integer orders and 432 orders at actual
  real-valued input doubles against untruncated fraction-free elimination;
- 64 additional orders on a source whose determinant order is 420 while
  the series window is 256, and eight orders in the forced prefix-resume case;
- at least 384 forced repeat hits, diagonal reordering, recomputation outside
  a cone, a true rank-two source, rejection of an insufficient window and
  successful reconstruction at the larger window.

The [receipt](../evidence/smith-chambers-f4-controls.json.xz) retains the
actual matrices, queries, comparisons, cache diagnostics and source hashes.
The independent profile checks use exact field polynomials and rational
weights; floating output is compared after its final conversion.

An actual thirty-two-letter original-Fock pilot compares five nearby queries
with the legacy engine: **10,240 complete orders** and **1,285 gradient
components** agree. Its source uses the identity constant flag; the small
controls include genuinely nonbinary matrices. The
[pilot packet](../evidence/smith-chambers-q32-pilots.json.xz) retains each
measured implementation, its source text, full results and timings. The
moving trajectory and identical repeats are reported separately. These
comparisons do not replace an untruncated reconstruction of every larger
Fock profile.

In the latest run, five legacy trajectory queries took **15.67 seconds**;
five cached queries including cold frame construction took **19.67 seconds**.
Installing the cache, including compilation in that process, took a further
**3.47 seconds**. Twenty identical queries after construction took **0.201
seconds**. The moving portion made 90 block queries, with 39 hits, 51 builds
and 898 reused pivots. No fallback or unresolved tracked construction occurred.
All earlier slower measurements are retained. These observations justify
keeping the feature optional, rather than using repeat-only timings as a
claim of faster moving searches.

## 5. Search integration and limits

Use `--cached-profiles` with the
[F4 coordinate driver](../tools/search_f4_pascal_coordinate_flags.py) to
opt in. It records source texts, complete profiles, cache hits, rebuilds,
reused pivots, legacy fallbacks and periodic query counts. The existing
frozen one-hour searches continue with their original implementation.

No universal runtime speedup is asserted. Repeated points are cheap once
a cone has been built; moving points can require expensive rebuilds. The
larger pressure searches still require both physical costs and every promised
profile order. A numerical candidate still needs independent exact replay
and an outward pressure certificate before it supports a stronger crossing
or factoring claim. Substantial OpenAI assistant support was used. No new
Lean formalization or improved score is included.
