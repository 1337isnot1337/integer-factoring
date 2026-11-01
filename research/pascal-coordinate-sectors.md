# Complete coordinate crossings without ambient Fock expansion

## Construction

The [whole-symbol evaluator](pascal-symbol-potentials.md) now has a
selected-entry occupation-sector backend. It computes arbitrary
coordinate subalphabets of larger Pascal mode sets without writing
the complete ambient Fock tensor. All nonempty sectors and every
promised Smith order remain in the pressure calculation. The backend
is not a low-order or unknown-tail estimate.

The [constructor](../tools/pascal_coordinate_sectors.py) evaluates each
required original minor of the polynomial matrix `D`. Cached exact
determinant expansion chooses a sparse row. Its cache is bounded at
`2^20` entries and cleared after construction; clearing changes only
recomputation cost. The current finite implementation permits up to
sixteen modes and sixty-four selected coordinate words. Individual
entries are binary polynomials, with no truncated determinant expansion.

The exact identity

```math
 T(b,x,v,y)\ne0\quad\Longrightarrow\quad
 |b|+|x|=|v|+|y|
```

holds because each entry is a minor with row set `v union (y+L)`
and column set `b union (x+L)`. Thus the first cut separates into
sectors `|b|-|v|=|y|-|x|`, and the second separates into sectors
`|b|+|x|=|v|+|y|`. Diagonal symbol potentials preserve these sectors.
Each sector is square and generically invertible by the triangular
endpoint argument in the preceding note. The complete profile is
the multiset union of its sector profiles.

This is an exact block decomposition at both places. For sixteen
physical letters, the full cut has order 256 while its largest
occupation block has order 70. The selected-entry constructor also
computes eight-mode examples with thirty-two physical letters; their
ambient alphabet would have 256 symbols. No stronger factoring score
follows from the construction alone.

## Common-scale initialization

Scaling both physical potentials and the common valuation parameter by
the same positive `t` scales every complete Smith order by `t`. One
complete profile query therefore supplies the entire scalar curve for
an inherited start. The scalar optimizer need not repeat native
elimination at each scale. The resulting minimization remains numerical;
its best value is not an exact pressure certificate.

Use `--backend sectors` in the existing
[driver](../tools/search_pascal_letter_gauges.py), with an explicit
`--alphabet` for larger mode sets. Both independent matching-pair
physical costs and arbitrary symbol potentials are unchanged.

## Exact finite validation

The [checker](../tools/check_pascal_coordinate_sectors.py) compares all
selected entries with the previous full-Fock constructor for three
small kernels, including a nonmonotone exponent permutation. For an
eight-mode nine-letter example, a separate subset dynamic programme
checks every one of its 6,561 original entries.

Nine complete integer-profile comparisons agree at **1,782 orders**
between the occupation blocks and the previous dense original-Fock
backend. Twelve nonadditive gradient comparisons are numerical. Three
positive-scale controls check the full profile scaling identity, and
the scalar-start control recovers the retained ordered reference.

```sh
OPENBLAS_NUM_THREADS=1 python3 tools/check_pascal_coordinate_sectors.py \
  --output /tmp/pascal-coordinate-sectors-controls.json
```

The [receipt](../evidence/pascal-coordinate-sectors-controls.json)
retains the actual parameters, sector sizes, comparison counts and
source hashes. These controls do not claim an independent reconstruction
of every larger example.

## Filled-sea coordinate starts

The first low-height thirty-two-word restriction at eight modes gave
a weak numerical pressure, approximately `0.85730`. Its inherited
physical energy favors a filled half interval, so low-height words
around the empty vacuum omit thermally important states.

For the ordered eight-mode kernel, the filled reference is word 15.
The excitation costs are `14,10,6,2,2,6,10,14`. Select the lowest-cost
words relative to this reference, with deterministic integer tie breaks.
The [registered packet](../evidence/pascal-fermi-coordinate-starts.json.xz)
contains the full candidate list and complete profiles for its first
three completed cases:

| Selected letters | Inherited numerical pressure |
|---:|---:|
| 16 | 0.728681631882 |
| 24 | 0.723799386492 |
| 32 | 0.721693325006 |

These values are exploratory common-scale fits. They remain above
`1/sqrt(2)` and do not improve the strongest existing crossing
certificate. The larger fourth case is explicitly pending in this
checkpoint. Arbitrary symbol-potential optimization, different flags,
and changes to the physical maps remain constructive next questions.

This work retains the preceding Pascal, original-Fock, paired-type,
and native Smith provenance. Substantial OpenAI assistant support was
used for the construction, proof, computations and text. No new Lean
formalization or improved uniform factoring exponent is claimed.
