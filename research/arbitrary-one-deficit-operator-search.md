# Exact operators on arbitrary one-deficit selections

## Construction and API

The [arbitrary-selection geometry theorem](arbitrary-one-deficit-identity-cores.md)
extends the prime-field operator evaluator from end selections to any
increasing tuple of actual indices at \(n=p^k-1\), provided
\(r+1\le n\). The reduction uses \(2r+1\) actual active
zero modes, a \(p^{r+1}\)-dimensional auxiliary Weyl matrix,
and exact removal of its single formal Gaussian factor.

The implementation is the existing
[`prime_deficit_identity_operators.py`](../tools/prime_deficit_identity_operators.py),
with an optional final `selected_indices` tuple. Existing positional
calls retain their previous meaning. With an explicit selection, the
`initial_modes` and `terminal_modes` arguments specify only its total
size; all physical energies, spectators, endpoint conditions, and
degree accounting use the supplied indices.

For example:

```python
from prime_deficit_identity_operators import zero_profiles, infinity_profiles

selected = (7, 18)
zero = zero_profiles(3, 27, 1, 0, 2, 1, kernel,
                     selected_indices=selected)
infinity = infinity_profiles(3, 27, 1, 0, 2, 1, kernel,
                             selected_indices=selected)
```

The same keyword is available in the literal active tensor, original
physical seed, complete profile, and pressure routines. Every compact
one-deficit contraction is cross-checked by exact polynomial products
against the independent adjugate contraction before use. Input
polynomials remain complete; no fixed source-degree window is used.

An admissible physical crossing still needs either the exact identity
endpoint or an invertible energy-filtered endpoint with the existing
constant physical repair. The arbitrary-index API does not relax
that condition. The local eigenbasis from the geometry theorem is
not applied to the physical endpoint.

## 1. Independent exact controls

The new
[`checker`](../tools/check_selected_deficit_operators.py) and
[`saved evidence`](../evidence/selected-deficit-identity-operators.json)
contain the following completed controls:

| Control | Completed scope |
|---|---:|
| Complete original physical cut profiles | 10 |
| Original endpoint coefficient equalities | 4,328 |
| Literal rectangular core versus Weyl profiles | 3 |
| Prior end-selection API comparisons, both cuts | 5 |
| Arbitrary-selection Gaussian recovery at \(n=26\) | 4 |
| Independent complete pressure comparisons | 36 |
| Invalid-selection rejections | 5 |

The rectangular comparisons include a nontrivial qutrit operator at
\(n=8\), \(S=(2,5)\): its literal zero matrix has order
\(3^5=243\), while the Weyl matrix has order \(3^3=27\).
The original physical controls use complete binary matrices of order
64 and qutrit matrices of order 81, both ramification values 1 and 3,
and identity or filtered endpoints. These tests include the initial
two-index selection that was absent from the earlier end-only API.

The independent pressure expansion restores all Gaussian spectators
and both physical paired-type partitions. Its largest observed
absolute discrepancy was below \(3.0\times10^{-13}\).
The previous prime-deficit checker also passes without changing its
default arguments or saved evidence.

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/check_selected_deficit_operators.py --extended
```

## 2. A separately specified parent-exponent family

The scalar part of the one-deficit contraction contains
\(1-\tau^N\), with parent length \(N=27\). This motivates
testing the parent exponents \(\pm27\), separately from the
previously completed actual-length exponents \(\pm26\).

For two qutrit directions, let \(P_s\) run over every rank-one
projector in one of the occupation blocks of dimensions \(2,3,2\).
The family is

\[
F(q)=I+\lambda(1-q^h)P_s,
\qquad \lambda\in\mathbb F_3^*,\qquad
h\in\{-27,-1,1,27\}.
\tag{1}
\]

The projector counts are \(12,117,12\), so (1) contains exactly
1,128 nonidentity kernels, together with the Gaussian identity.
Every kernel has endpoint exactly \(I\), including all
off-diagonal entries of its chosen projector. The one common Laurent
clearing has exponent at most 27. At \(n=26\), \(P=1\),
the sufficient complete physical polynomial degree bound is

\[
8n^2+8n+27=5643.
\tag{2}
\]

The [search driver](../tools/search_qutrit_deficit_operators.py)
records the exponent convention explicitly. Its default remains the
earlier actual-length family. An arbitrary-selection parent-family
run is reproducible with:

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/search_qutrit_deficit_operators.py \
  --parent 27 --selected 7,18 --exponent-length parent \
  --output /tmp/qutrit-parent-interior.json
```

These finite pressure screens compare each exact profile with the
Gaussian at the same aggregate charge and scale. A numerical absence
of improvement in such a screen is not an all-scale obstruction.

All four untwisted parent-exponent runs completed. The
[`saved catalogue`](../evidence/qutrit-parent-exponent-projectors-n26.json)
contains every exact profile class and a representative kernel.

| Actual selected indices | Retained auxiliary indices | Kernels | Exact profile classes |
|---|---|---:|---:|
| \((24,25)\) | \((0,1,2)\) | 1,129 | 101 |
| \((0,25)\) | \((0,1,25)\) | 1,129 | 207 |
| \((7,18)\) | \((0,8,19)\) | 1,129 | 101 |
| \((0,13)\) | \((0,13,25)\) | 1,129 | 207 |

Each run evaluated a 72-point charge/scale grid and refined its 30
best profile classes. The best refined critical pressure was
\(0.07235825135226\), attained by the Gaussian. The smallest
reported relative grid values differed from zero by less than
\(1.7\times10^{-14}\). None gives a demonstrated improvement.
