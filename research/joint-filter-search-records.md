# Joint polynomial identity-filter search records

## What is established

The [residue synthesis theorem](joint-binary-filter-synthesis.md) gives
an exact degree/order feasibility oracle for joint diagonal identity
filters in every prime characteristic. The
[prime-power identity theorem](prime-power-joint-identity-filters.md)
then supplies both physical profiles for every selected mode subset.
The records below contain completed finite searches using these tools.

No negative critical pressure or improved factoring exponent was found
in the reported pressure grids. This is a finite search observation.
The exact theorems are the profile reduction, residue compression,
endpoint criterion, and literal witness construction.

## 1. Completed coefficient catalogues

All binary catalogues use \(n=32\). The qutrit catalogues use \(n=27\).
Here \(P\) denotes ramification, \(r\) the number of selected modes,
and the offset is added to every degree cap, including the full
occupation state. The full degree grids and row priority permutations
are stored in each record. Except for the stated individual-cap search,
degree caps depend on total occupation.

| Catalogue | Characteristic | \(r\) | \(P\) | Offset | Degree vectors | Exact feasibility queries | Distinct bound pairs | Exact ordered profiles |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| [Binary three modes](../evidence/joint-binary-polynomial-n32-r3.json) | 2 | 3 | 1 | 0 | 512 | 786,895 | 3,244 | 1,478 |
| [Binary four modes](../evidence/joint-binary-polynomial-n32-r4.json) | 2 | 4 | 1 | 0 | 256 | 415,147 | 789 | 302 |
| [Binary ramification three](../evidence/joint-binary-polynomial-n32-p3-r3.json) | 2 | 3 | 3 | 0 | 512 | 896,387 | 2,773 | 841 |
| [Binary offset 64](../evidence/joint-binary-polynomial-n32-r3-offset64.json) | 2 | 3 | 1 | 64 | 216 | 415,466 | 1,510 | 705 |
| [Four modes, offset 32](../evidence/joint-binary-polynomial-n32-r4-offset32.json) | 2 | 4 | 1 | 32 | 256 | 471,625 | 1,053 | 471 |
| [Individual caps, offset 32](../evidence/joint-binary-polynomial-n32-r3-individual-offset32.json) | 2 | 3 | 1 | 32 | 2,187 | 956,297 | 10,504 | 2,361 |
| [Qutrit two modes](../evidence/joint-qutrit-polynomial-n27-r2.json) | 3 | 2 | 1 | 0 | 256 | 149,068 | 762 | 350 |
| [Qutrit offset 27](../evidence/joint-qutrit-polynomial-n27-r2-offset27.json) | 3 | 2 | 1 | 27 | 256 | 182,544 | 873 | 391 |

The filenames containing `p3` in the binary records refer to
\(P=3\), not characteristic three. Binary three-mode witnesses use
\(\mathbf F_{16}\), binary four-mode witnesses use
\(\mathbf F_{32}\), and the joint qutrit witnesses use
\(\mathbf F_{27}\).

The [catalogue replay](../evidence/joint-filter-catalogue-replay.json)
reconstructs all **6,899 literal polynomial witnesses** and verifies
62,117 nonzero endpoint values. It checks every stored alpha degree,
beta order, endpoint image dimension, coefficient field order, witness
hash, catalogue hash, and bound-pair multiplicity. Repeated profiles in
different catalogues are retained as separate reproducibility records.

The greedy row-priority procedure returns exact feasible targets, but
does not enumerate the entire Pareto set. A catalogue is complete for
its stated generator and finite inputs, not for all polynomial filters.

## 2. Completed all-index pressure screens

An ordered coefficient profile can be reused at every subset of
physical indices. The expectation identity in the synthesis theorem
evaluates its exact profile cost. Floating point is used only for the
pressure comparison.

The main binary scale grid is \(t\in\{2,4,5,6,8\}/32\), with
\(J\in\{-1,0,1\}\). The qutrit grid uses the same numerator set
divided by 27 and
\(J\in\{-27,-9,-1,0,1,9,27\}\). Each record gives its full grids.

| Screen | Ordered profiles after stated mode assignments | Selected subsets | Scale/charge pairs | Profile/subset/grid evaluations |
|---|---:|---:|---:|---:|
| [Binary three modes](../evidence/all-index-three-mode-residue-polynomial-n32-screen.json) | 4,479 | 4,960 | 15 | 333,237,600 |
| [Binary four modes](../evidence/all-index-four-mode-residue-polynomial-n32-screen.json) | 2,970 | 35,960 | 15 | 1,602,018,000 |
| [Binary \(P=3\), stored assignment](../evidence/all-index-three-mode-residue-polynomial-n32-p3-screen.json) | 841 | 4,960 | 15 | 62,570,400 |
| [Binary offset 64, stored assignment](../evidence/all-index-three-mode-residue-polynomial-n32-offset64-screen.json) | 705 | 4,960 | 15 | 52,452,000 |
| [Four-mode offset 32, stored assignment](../evidence/all-index-four-mode-residue-polynomial-n32-offset32-screen.json) | 471 | 35,960 | 15 | 254,057,400 |
| [Individual caps, every mode assignment](../evidence/all-index-three-mode-residue-polynomial-n32-individual-offset32-screen.json) | 7,402 | 4,960 | 15 | 550,708,800 |
| [Qutrit two modes](../evidence/all-index-qutrit-polynomial-n27-r2-screen.json) | 657 | 351 | 35 | 8,071,245 |
| [Qutrit offset 27](../evidence/all-index-qutrit-polynomial-n27-r2-offset27-screen.json) | 757 | 351 | 35 | 9,299,745 |

Rows without the qualification “stored assignment” include all
permutations of the filter modes. Thus each physical subset is tested
with every assignment of the retained coefficient profiles to its
indices. The qualified rows evaluate the catalogue's stored ordering;
they make no additional permutation-completeness claim.

Three completed broader-charge records are also retained:

| Screen | Ordered profiles | Subsets | Scale/charge pairs | Evaluations |
|---|---:|---:|---:|---:|
| [Three-mode \(\mathbf F_4\) monomials](../evidence/all-index-three-mode-f4-n32-wide-screen.json) | 3,850 | 4,960 | 55 | 1,050,280,000 |
| [Polynomial offset 64, every mode assignment](../evidence/all-index-three-mode-residue-polynomial-n32-offset64-wide-screen.json) | 2,161 | 4,960 | 55 | 589,520,800 |
| [Four-mode offset 32, every mode assignment](../evidence/all-index-four-mode-residue-polynomial-n32-offset32-wide-screen.json) | 5,124 | 35,960 | 15 | 2,763,885,600 |

The first two use \(t\in\{2,4,6,8,16\}/32\) and
\(J\in\{-63,-32,-16,-8,-3,0,3,8,16,32,63\}\).
The last uses \(t\in\{4,6,8\}/32\) and
\(J\in\{-32,-16,0,16,32\}\).
The \(\mathbf F_4\) family was already closed under mode permutations;
its screen does not add a separate permutation expansion.

All these grids found no same-gauge cost gain below the detection
tolerance \(10^{-9}\), and no negative critical pressure. Tiny negative
relative logarithms, at most about \(8\cdot10^{-15}\) in magnitude,
are floating roundoff at Gaussian equality. The minimum critical
pressure in the binary grids was approximately 0.0333616450; in the
qutrit grids it was approximately 0.0703114174. These are grid minima,
not globally optimized lower bounds.

The [older transported-row screen](../evidence/all-index-three-mode-polynomial-n32-screen.json)
evaluates 2,044 ordered profiles derived from 4,270 fixed-basis row
certificates. Its scope is those sufficient certificates. It does not
replace the literal coefficient catalogues above.

## 3. Reproduction

Run the independent feasibility and complete catalogue gates from the
repository root:

```sh
python tools/check_joint_binary_filter_synthesis.py
python tools/check_joint_prime_filter_synthesis.py
python tools/check_joint_filter_catalogues.py
```

The first two gates compare against uncompressed coefficient systems:
288 binary comparisons with 32 literal replays, and 100 prime-state
comparisons with 23 independent tensor-transform replays. The catalogue
gate is a portable replay of every retained witness, not a rerun of
all feasibility queries from the search.

Representative catalogue and pressure commands are:

```sh
python tools/search_joint_binary_filter_synthesis.py --n 32 --modes 3 --output /tmp/joint-r3.json
python tools/search_joint_binary_filter_synthesis.py --n 32 --modes 4 --degree-grid 0,1,32,33 --priority-limit 16 --output /tmp/joint-r4.json
python tools/search_joint_prime_filter_synthesis.py --prime 3 --n 27 --modes 2 --priority-limit 12 --output /tmp/joint-qutrit.json
python tools/search_all_index_joint_filters.py --n 32 --modes 3 --joint-source evidence/joint-binary-polynomial-n32-r3.json --permute-modes --output /tmp/joint-r3-screen.json
python tools/search_all_prime_index_filters.py --source evidence/joint-qutrit-polynomial-n27-r2.json --output /tmp/joint-qutrit-screen.json
```

To reproduce another catalogue, pass its stored degree grid, offset,
ramification, degree family, and row-priority limit. Each complete
record stores the actual priority permutations, so the search domain
is inspectable independently of the command-line defaults.

## 4. Interrupted work

The execution environment became unavailable during additional
fractional-degree binary three-/four-mode and qutrit two-mode searches,
continuous index/charge optimization, and a near-prime-power
\(n=31\) flag screen. Their temporary outputs did not survive the
outage. They are excluded from all completed counts and pressure
claims in this note. No conclusion is drawn from those attempts.

The surviving individual-cap binary catalogue is complete, and its
all-subset discrete screen was rerun to completion after recovery.
Its separate continuous optimizer was interrupted. Continuous choices
of real or repeated indices would in any case require a fresh replay
at distinct physical indices before becoming an admissible witness.

Near-prime-power calculations that were rerun to completion after the
outage are documented separately in the
[flag and unit-reduction note](weyl-filter-flags-and-unit-reduction.md).
