# Energy-conjugated operators in one-deficit cores

The arbitrary-index one-deficit reduction permits polynomial dependence
inside an occupation block as well as between its eigenvalues. A simple
additional exact family comes from conjugating a constant projector by
selected-index monomials. This changes the physical kernel; it is not
treated as a free basis change at either cut.

## 1. Construction and endpoint

For an occupation label \(a\), put

\[
E_S(a)=\sum_{j\in S}j a_j,
\qquad
D_\epsilon(q)_{a,a}=q^{\epsilon P E_S(a)},
\qquad \epsilon\in\{-1,1\}.
\tag{1}
\]

Starting from an occupation-preserving kernel \(F\), define

\[
F_\epsilon(q)=D_\epsilon(q)F(q)D_\epsilon(q)^{-1}.
\tag{2}
\]

In entries, (2) multiplies \(F_{a,b}\) by
\(q^{\epsilon P(E_S(a)-E_S(b))}\). This operation preserves
occupation and has the exact inverse \(\epsilon\mapsto-\epsilon\).
Since \(D_\epsilon(1)=I\), it preserves the full endpoint
matrix. In particular, an identity endpoint remains the identity.

Every coefficient is still a finite Laurent polynomial. A single
common monomial clearing gives the admissible physical polynomial
crossing. Both resulting Smith profiles are computed from these
changed coefficients, and both physical paired-type partitions remain
those of the original crossing alphabet. No intermediate diagonal
matrix in (1) is assumed regular for free.

The implementation is `energy_conjugate` in the
[`qutrit search driver`](../tools/search_qutrit_deficit_operators.py).
It is tested independently by
[`check_deficit_energy_conjugation.py`](../tools/check_deficit_energy_conjugation.py).
The [exact record](../evidence/deficit-energy-conjugation.json) contains
4,516 endpoint and inverse-conjugation controls, and four complete
original binary cut profiles of order 64. The latter include both
signs, ramification 1 and 3, and all 1,569 endpoint coefficients.

## 2. Complete degree accounting

For two qutrit indices \(j_0<j_1\), occupation preservation gives

\[
|E_S(a)-E_S(b)|\le2(j_1-j_0).
\tag{3}
\]

For a projector kernel
\(I+\lambda(1-q^h)P_s\), with \(|h|\le27\),
\(P=1\), the coefficient support after (2) has width at most

\[
27+4(j_1-j_0).
\tag{4}
\]

Each representative records its actual clearing and resulting degree;
(4) is only a uniform upper bound. With \(n=26\), the complete
physical polynomial degree is bounded by \(5616\) plus (4).
Thus \(S=(7,18)\) has sufficient bound 5687, while
\(S=(0,13)\) has sufficient bound 5695. This includes the one
common clearing; no degree or zero-profile penalty is omitted.

## 3. Specified finite family

The base family has all 1,129 kernels from the
[parent-exponent catalogue](arbitrary-one-deficit-operator-search.md):
all rank-one projectors in the three nontrivial occupation blocks,
both nonzero coefficients in \(\mathbb F_3\), and
\(h\in\{-27,-1,1,27\}\), together with the identity.
Conjugation is bijective, so each fixed selection and sign still
contains exactly 1,129 distinct kernels.

A run is reproduced by:

```bash
OPENBLAS_NUM_THREADS=1 python3 tools/search_qutrit_deficit_operators.py \
  --parent 27 --selected 7,18 --exponent-length parent \
  --energy-twist -1 --output /tmp/qutrit-energy-interior-minus.json
```

The numerical objective includes both exact active profiles, every
Gaussian spectator, and the complete physical paired-type costs.
Finite numerical screening does not certify an all-scale obstruction.

All four specified runs completed, with the full profile-class records
saved in
[`qutrit-energy-conjugated-projectors-n26.json`](../evidence/qutrit-energy-conjugated-projectors-n26.json).

| Selected indices | Conjugation sign | Kernels | Exact profile classes |
|---|---:|---:|---:|
| \((7,18)\) | \(-1\) | 1,129 | 101 |
| \((7,18)\) | \(+1\) | 1,129 | 313 |
| \((0,13)\) | \(-1\) | 1,129 | 207 |
| \((0,13)\) | \(+1\) | 1,129 | 317 |

Each run used the same 72-point charge/scale grid and refined its 30
best profile classes. The Gaussian again attained the best refined
critical pressure, \(0.07235825135226\). No relative improvement
exceeded numerical roundoff. These 4,516 kernel-selection-sign cases
are separate from the 4,516 untwisted parent-exponent cases.
