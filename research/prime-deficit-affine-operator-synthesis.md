# Affine operator synthesis over the prime field

The qutrit projector catalogues restrict the polynomial dependence to
one scalar function times a fixed idempotent. This note removes that
restriction within a prescribed Laurent interval and prescribed local
row lattices. The endpoint remains exactly the identity; no extension
field or additional endpoint repair is assumed.

For two qutrit directions, the occupation blocks have sizes
\(1,2,3,2,1\). There are therefore
\(1+4+9+4+1=19\) independent matrix coefficient polynomials.
The interval \(q^{-1},1,q\) gives 57 coefficient variables over
\(\mathbb F_3\). The source of each variable is the exact
[arbitrary-selection Weyl construction](arbitrary-one-deficit-operator-search.md).

## 1. A finite polynomial row frame is sufficient

Let \(M_G(q)\) be one Gaussian occupation block of the reduced
Weyl matrix. Choose a positive coefficient precision \(T\).
Gaussian elimination modulo \(q^T\) proceeds by a minimum-valuation
pivot, row and column permutations, and operations

\[
R_i\longmapsto uR_i-bR_k,
\qquad M_{kk}=q^t u,\qquad M_{ik}=q^t b,
\qquad u(0)\ne0.
\tag{1}
\]

All multipliers in (1) are polynomials. Recording the row operations
modulo \(q^T\) gives a finite polynomial matrix \(U_T\).
Its constant term is invertible: every operation in (1) has invertible
constant determinant, and truncation does not change that constant
term. Thus the chosen polynomial representative \(U_T\) is an
actual invertible matrix over \(\mathbb F_p[[q]]\).

Column permutations only reorder coordinates. Once a pivot column
has been cleared below its pivot, subsequent column operations that
would clear its pivot row do not change the remaining Schur block.
The recorded pivot orders are therefore the Smith orders whenever
all pivots are visible. The implementation checks these orders against
the independent complete-polynomial Smith routine, and verifies the
full congruence \(U_TM_G\) modulo \(q^T\).

Only row floors are needed afterward. A prescribed bound
\(\operatorname{ord}_q(U_TM)_{ij}\ge z_i\) is a sufficient
local lattice condition. It is not asserted to describe every matrix
with a given Smith profile.

## 2. Exact coefficient equations and precision

Write the operator as

\[
F_{a,b}(q)=\sum_{d=L}^{H}f_{a,b,d}q^d,
\qquad |a|=|b|.
\tag{2}
\]

One may instead prescribe a finite subset of the integer powers in
\([L,H]\), by omitting the other variables. The sparse parent-power
search below uses this option; it does not contain every power between
\(-27\) and \(27\).

The complete Weyl sources before multiplication by (2) are
polynomials of nonnegative degree. Their transformed coefficients
below \(q^T\) are exact. Multiplication by (2) can lower an
omitted source degree only to \(T+L\). Consequently the row
conditions are exact whenever

\[
\boxed{\max_i z_i\le T+L.}
\tag{3}
\]

The constructor rejects a query outside (3). It does not assume that
an omitted coefficient vanishes. Truncation determines which
particular regular polynomial frame is used; it does not weaken the
subsequent coefficient equations.

At infinity, use the original occupation rows and column energies

\[
\kappa(b)=P\sum_{j\in S}(2j+1-2n)b_j.
\tag{4}
\]

The coefficient \(f_{a,b,d}\) has exact infinity order
\(\kappa(b)-2d\). A target row floor \(h_a\) simply forces
that coefficient to zero when \(\kappa(b)-2d<h_a\).

More generally, a fixed invertible constant matrix \(B\), preserving
occupation, may act on the infinity rows. The corresponding equation
is \(\sum_a B_{i,a}f_{a,b,d}=0\) whenever
\(\kappa(b)-2d<h_i\). This is again linear. Since \(B\) is
a constant local unit, it preserves the infinity Smith profile.
It is never applied as an unverified global endpoint repair.

A concrete identity-endpoint starting kernel exists for such a flag.
For nonnegative integers \(g_i,H\), put

\[
F(q)=B^{-1}\operatorname{diag}(q^{-g_i})
       +q^{-H}(I-B^{-1}).
\tag{4a}
\]

Then \(F(1)=I\), and
\(BF=\operatorname{diag}(q^{-g_i})+q^{-H}(B-I)\).
It has infinity row floors \(\kappa_i+2g_i\) whenever
\(\kappa_j+2H\ge\kappa_i+2g_i\) for each nonzero entry
\((B-I)_{ij}\). The polynomial coefficient system can subsequently
seek better zero floors while preserving those infinity targets.

Finally impose the affine equations

\[
\sum_{d=L}^{H}f_{a,b,d}=\delta_{a,b}.
\tag{5}
\]

These enforce the full identity endpoint. All constraints are linear
over the original prime field. Reducing each zero-row prefix to a
basis changes neither its solution space nor the endpoint equations.
An exact solution supplies a concrete kernel; an inconsistent system
excludes only these fixed row frames, targets, and coefficient bounds.

## 3. Physical accounting and implementation

The output kernel is passed to the previously checked complete-profile
evaluator. Its formal dummy factor is removed exactly, its remaining
spectators are restored, and both physical paired-type partitions are
included in the pressure. One common Laurent clearing has the same
degree accounting as the existing operator construction. Neither the
local frame \(U_T\) nor its determinant is multiplied into the
global physical tensor.

The prototype is
[`prime_deficit_operator_synthesis.py`](../tools/prime_deficit_operator_synthesis.py).
`row_source` constructs the regular Gaussian row frame;
`coefficient_system` fixes the complete coefficient interval;
`synthesize` returns a checked identity-endpoint kernel or `None`.
Every returned affine solution is substituted back into all equations
before its kernel is accepted.

## 4. Independent verification

The
[`checker`](../tools/check_prime_deficit_operator_synthesis.py) and
[`saved record`](../evidence/prime-deficit-operator-synthesis.json)
contain these completed controls:

| Control | Completed count |
|---|---:|
| Literal transformed polynomial entry comparisons | 1,692 |
| Complete Smith-profile prefix comparisons | 84 |
| Affine target queries checked by exhaustive coefficient enumeration | 27 |
| Complete assignments in each exhaustive query | 19,683 |
| Nontrivial constant infinity-flag witnesses | 6 |
| Original physical cut profiles | 4 |
| Original endpoint coefficient equalities | 1,649 |
| Deliberately insufficient-precision query rejected | 1 |

The source comparisons include both ordinary powers \(-1,0,1\)
and sparse parent powers \(-27,-1,0,1,27\). The full physical
controls include qutrit matrices of order 81 and binary matrices of
order 64. Two other assistant workers also reviewed the row-frame, affine-equation,
and constant-flag arguments during this session.

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/check_prime_deficit_operator_synthesis.py
```

## 5. Completed finite synthesis paths

The
[`search driver`](../tools/search_prime_deficit_affine_operators.py)
starts from feasible row targets, raises individual zero or infinity
floors in two specified coordinate orders, and solves each affine
system exactly. It profiles the initial kernel and distinct terminal outputs of those
paths. It does not profile every intermediate feasible witness.
These paths explore selected row lattices;
they do not enumerate the entire polynomial coefficient space.

All four runs completed. Their
[`combined record`](../evidence/qutrit-deficit-affine-operator-searches.json)
contains every distinct returned profile and its exact kernel.

| Selection | Allowed powers | Infinity frame | Paths | Exact affine queries | Distinct returned profiles |
|---|---|---|---:|---:|---:|
| \((7,18)\) | \(-1,0,1\) | Identity | 72 | 523 | 3 |
| \((0,13)\) | \(-1,0,1\) | Identity | 72 | 523 | 3 |
| \((0,13)\) | \(-27,-1,0,1,27\) | Identity | 90 | 2,608 | 13 |
| \((0,13)\) | \(-27,-1,0,1,27\) | Pascal flag | 8 | 1,331 | 8 |

The nontrivial flag is the upper Pascal matrix
\(\left(\begin{smallmatrix}1&1&1\\0&1&2\\0&0&1\end{smallmatrix}\right)\)
on the occupation-two block and the identity elsewhere. Formula
(4a), with \(H=27\), proves that its starting targets are feasible.

The 242 paths made 4,985 exact affine queries. The table counts 27
profile records across the four runs, corresponding to 23 distinct
selection/profile pairs. Each distinct output
received the same 72-point full-cost pressure grid, with the best
profiles refined numerically. The best refined critical pressure was
again the Gaussian value \(0.07235825135226\); relative deviations
below \(1.7\times10^{-14}\) are numerical roundoff. Thus this
completed synthesis step supplies a broader exact search constructor,
but no demonstrated improvement in the factoring bound.

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/search_prime_deficit_affine_operators.py \
  --selected 0,13 --powers parent --infinity-flag pascal \
  --output /tmp/qutrit-affine-pascal.json
```

## 6. Fresh replay against the final source

The original runs did not retain historical source-hash manifests. A
separate [current-source replay](../evidence/qutrit-deficit-affine-current-source-replay.json)
checks all 27 saved kernels against the final published implementation.
All 16,721 affine coefficient equations, all 54 complete cut profiles,
every endpoint, every degree bound, and every saved pressure point
match. The record includes the source hashes used for that fresh replay.
It makes no claim about an unrecorded historical source hash.

All four complete search drivers were also rerun against the final source.
Every saved JSON field matched, apart from elapsed seconds. The
[full-run replay manifest](../evidence/qutrit-affine-complete-current-source-replay.json)
records all 242 paths, 4,985 affine queries, the per-run profile counts,
and the source and original-record hashes. This validates present-source
reproducibility of the entire prescribed schedule.
