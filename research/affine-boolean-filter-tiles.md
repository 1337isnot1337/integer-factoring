# Affine Boolean relabelings of symmetric polynomial tiles

## Status and scope

This note gives an exact way to leave occupation symmetry while retaining all
coefficient polynomials, endpoint values, and infinity degrees. It also records
a completed finite search. The search did not improve the genuine two-cut
ratio: every recorded optimum chose the unmodified Gaussian tile on every
block. This is a finite search result, not an exclusion of arbitrary affine
relabelings or arbitrary polynomial filters.

The input family is the 517 four-mode profiles in
[`symmetric-tiled-filters-r4.json`](../evidence/symmetric-tiled-filters-r4.json),
constructed at source parameters \(m=4,12,20\). See
[`symmetric-multilevel-filter-synthesis.md`](symmetric-multilevel-filter-synthesis.md)
for the coefficient construction and exact disjoint-block cost formula.

## 1. Relabel the coefficients, then recompute the transform

Let \(f_0(x),\ldots,f_4(x)\) be a symmetric four-mode coefficient witness over
a finite extension of \(\mathbf F_2\). For an invertible affine Boolean map
\(a(k)=Ak+b\), define the new coefficient polynomials by

\[
 \widetilde f_k(x)=f_{|Ak+b|}(x),\qquad k\in\{0,1\}^4.
\]

Every endpoint remains nonzero. The complete zero-place transform is

\[
 \widetilde\beta_k(x)
   =\sum_{\ell\subseteq k}(1+x^m)^{|k|-|\ell|}
                              f_{|A\ell+b|}(x).                 \tag{1}
\]

In general, (1) is **not** obtained by permuting the old beta orders. The
Boolean-subset incidence relation in the transform depends on the original
physical occupation coordinates. For example, start with
\(f_0=1+x^8+x^{16}\) and \(f_1=\cdots=f_4=1\), at \(m=4\), and translate by
the first coordinate vector. The new coefficient at state zero is 1, so
\(\operatorname{ord}_0\widetilde\beta_0=0\). Permuting the old symmetric beta
orders would incorrectly give 4. The checker retains this negative control.

The implementation expands (1) literally in the additive binary planes of
the coefficient field. Multiplication by \((1+x^m)^j\) uses exactly the shifts
\(x^{mh}\) with \(h\) a binary submask of \(j\). No field multiplication is
needed in this transform, and the complete field-valued endpoint colors are
retained.

### Common monomial normalization

Let \(s=\min_k\operatorname{ord}_0\widetilde\beta_k\). The transform and its
inverse are polynomial matrices with determinant 1. Consequently \(x^s\)
divides every \(\widetilde f_k\), and it is legitimate to divide the entire
coefficient vector and its transformed vector by \(x^s\). After division,

\[
 B_k=\operatorname{ord}_0\widetilde\beta_k-s,
 \qquad D_k=\deg\widetilde f_k-s\ge0.                         \tag{2}
\]

Both endpoint colors and the product of the two relative costs are unchanged.
The two factors acquire reciprocal scalar monomial weights. Literal replay
checks divisibility of every coefficient polynomial, not just its degree.

If \(B_k\le m|k|\) for all states, the zero-place expectation is at least 1.
The infinity expectation is at least 1 by \(D_k\ge0\). Such a profile can be
discarded for every charge and every positive scale. This pruning criterion
is exact and makes no assumption about the physical occupation probabilities.

## 2. Complete finite orbit manifest

The tested linear maps are the identity, all 12 ordered one-bit XOR shears,
the four-bit Gray map, and its inverse. Each is combined with all 16
translations. Thus there are 240 specified affine maps. These 15 linear maps
form a chosen subfamily; they are not the full group \(\mathrm{GL}_4(2)\).
The exact matrix rows, translations, labels, and occupation weights are saved
in the evidence manifest.

Applying every map to every one of the 517 source witnesses gives:

| Quantity | Count |
|---|---:|
| Literal coefficient cases | 124,080 |
| Cases discarded by the all-state baseline criterion | 39,012 |
| Distinct retained full-state \((B,D)\) profiles | 65,475 |
| State coordinates per profile | 16 |
| Common degree denominator \(m\) | 60 |

Every source witness is regenerated and checked against the existing source
polynomial hash. Degree and order lists from \(m=4,12,20\) are transported to
\(m=60\) before profile deduplication. The manifest records a SHA-256 digest of
the complete deterministic ordered catalogue, together with the source
catalogue hash. The checker regenerates all 124,080 cases and the complete
digest. This avoids storing tens of megabytes of duplicate derived vectors
without weakening the finite-family specification. Seven literal
representatives, and every profile selected by a reported optimum, are stored
with all 16 coefficient and transformed polynomials.

## 3. Full physical costs and completed searches

For a physical block with its actual mode energies \(\lambda_i\), define
the normalized product measures

\[
 p(k)\ \propto\ \exp\!\left(-t\sum_i\lambda_i k_i\right),\qquad
 q(k)\ \propto\ \exp\!\left(-t\sum_i(\lambda_i+J)k_i\right).
\]

The complete relative log cost of a tile is

\[
 \log\sum_kq(k)\exp\bigl(2t(n|k|-B_k/P)\bigr)
 +\log\sum_kp(k)\exp(2tD_k/P).                             \tag{3}
\]

The implementation evaluates all 16 states in (3), using the actual mode
energies rather than a symmetric or average-bias approximation. Independent
blocks may select different profiles. The baseline is always an available
choice. Contiguous, interleaved, reflected hot/cold, and seeded shuffled
partitions are each tested at \(n=64\) and \(n=128\); every energy list is saved.

The optimizer minimizes the positive root \(\rho\) of the **full** pressure,
including the Gaussian term and (3), over the recorded scale and charge
domain. It does not treat a small absolute pressure at collapsing entropy as
an improved ratio. All four partitions at a fixed size returned the same
baseline result:

| Physical size | Best recorded genuine \(\rho\) | Full pressure at \(1/\sqrt2\) | Selected nonbaseline blocks |
|---|---:|---:|---:|
| 64 | 0.7075072433080506 | 0.0195069230211 | 0 |
| 128 | 0.7072351467479230 | 0.0109546875847 | 0 |

These numerical optimizer outputs are reproducible search evidence. They do
not certify a global lower bound over the family. In particular, the positive
critical pressures do not extend the proved homogeneous-filter exclusion to
this multilevel affine family.

## 4. Verification and reproduction

The independent checker regenerates the complete finite orbit, replays every
stored literal polynomial witness, and evaluates 64 additional source/map
pairs through an original subset-sum transform. Those comparisons cover 1,024
full-state transformed polynomials. It also checks the old-beta-permutation
negative control, validates every physical partition, and replays all eight
reported pressure roots with positive physical entropy.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/search_affine_symmetric_tiles.py \
  --output evidence/affine-symmetric-tiled-filters.json
OPENBLAS_NUM_THREADS=1 python3 tools/check_affine_symmetric_tiles.py
```

The implementation is split between
[`affine_symmetric_binary_filters.py`](../tools/affine_symmetric_binary_filters.py),
[`search_affine_symmetric_tiles.py`](../tools/search_affine_symmetric_tiles.py),
and [`check_affine_symmetric_tiles.py`](../tools/check_affine_symmetric_tiles.py).
The saved records are
[`affine-symmetric-tiled-filters.json`](../evidence/affine-symmetric-tiled-filters.json)
and [`affine-symmetric-tiled-filter-check.json`](../evidence/affine-symmetric-tiled-filter-check.json).
