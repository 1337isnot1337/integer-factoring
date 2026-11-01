# Distinct physical alphabets on the two Gaussian boundary pairs

The retained five-letter crossing uses the same four-particle alphabet
on all four boundaries. This checkpoint permits distinct alphabets
on the paired boundaries `X/V` and `Y/U`. It preserves the strict
endpoint and retains every polynomial coefficient and Smith order.

The completed searches did not improve the reference rate
`rho = 0.7263479477044683`. The one-symbol search covered all 325
replacements on `Y/U`; a separate low-grade search covered all 4,368
five-state choices from its sixteen-state shell. These are finite
screens at the stated natural physical heights, not exclusions for
other heights or all asymmetric alphabets.

## Exact construction

Let `C_X` and `C_Y` each contain five distinct four-element subsets
of eight modes. No reflection-admissibility condition is imposed.
The full strict Gaussian kernel has sixteen modes. For physical
labels `x,y,u,v`, the tensor entry is its minor with row mask

\[
 C_X(v)\;\cup\;(C_Y(y)+8)
\]

and column mask

\[
 C_Y(u)\;\cup\;(C_X(x)+8).
\]

At the strict endpoint, this is the paired crossing
`delta[x,v] delta[y,u]`. All entries have even exponents in the
original variable. The constructor compresses the square variable,
computes the exact common polynomial factor of all 625 entries,
and divides by it. That factor has value one at the endpoint.
The resulting normalized tensor is polynomial and still strict.

This is the same complete-minor construction as the
[balanced physical-code search](balanced-four-particle-codes.md),
with the two paired alphabets supplied separately. The
[aligned complement construction](aligned-gaussian-complement-perturbations.md)
gives the closed kernel and the even-exponent identity used here.

For a state `C`, define its grade by

\[
 g(C)=\sum_{i\in C}(2i+1).
\]

On each alphabet separately, subtract the minimum grade to obtain
energies `e_X` and `e_Y`. All energies are even. The natural local
heights used in this screen are

\[
 (0,0,e_Y/2,e_X/2)\quad\text{at infinity},\qquad
 (e_X/2,e_Y/2,0,0)\quad\text{at zero}.
\]

Thus the infinity matrix on `(x,u)|(y,v)` has row weight
`e_Y(u)/2` and column weight `e_X(v)/2`. The zero matrix on
`(x,y)|(u,v)` has row weight `(e_X(x)+e_Y(y))/2`.
The complete two-place profiles feed the paired-type cost with
the two distinct energy arrays.

## Complete finite schedules

The fixed `X/V` alphabet, in increasing grade order, is

\[
 C_X=(15,23,43,51,77),
\]

where integers encode eight-bit occupation masks.

| Screen | Declared schedule | Exact tensor queries | Best recorded rate |
|---|---:|---:|---:|
| One paired-symbol replacement | Reference plus all `5*65` replacements on `Y/U` | 326 | 0.7263479477044683 |
| Sixteen-state low-grade shell | Every five-element subset of the first sixteen states ordered by `(grade, mask)` | 4,312 new, 56 reused | 0.7263479477044683 |

The shell record contains every declared alphabet and identifies
the 56 cases reused from the first record. Its 4,368 covered
alphabets are distinct. The reuse is checked by the exact hash
and complete literal rows of the earlier record.

The best nonreflection alphabets in the one-symbol screen include
`C_Y=(15,23,43,51,216)` and `(15,23,43,51,228)`, with recorded
rate approximately `0.7415036109683197`. The latter is included
among the retained independent controls. The best genuinely new shell
case is `C_Y=(15,23,27,39,43)`, with rate approximately
`0.8433425419176845`. The reference remains best among the full
union of searched cases.

The decimal rates come from numerical pressure optimization.
No strict negative-pressure candidate was recorded. A finite
screen at these heights does not certify global optimality;
different heights remain a separate search direction.

## Independent controls

The [checker](../tools/check_asymmetric_four_particle_codes.py)
rebuilds both complete schedules independently and verifies their
coverage, reuse, source hashes, and retained records. For each
retained tensor, it expands all 625 determinants in the closed
sixteen-mode strict kernel using Python polynomial arithmetic.
It also checks the exact common-factor identity and every strict
endpoint entry.

Each retained complete Smith profile is then certified by regular
polynomial row and column frames on the literal weighted matrix.
This replay does not call the native selected-minor constructor
or native Smith oracle. Across the two records, it verifies
30,000 complete polynomial minor identities and 96 complete
profiles, including the best nonreflection controls. The baseline also
matches the retained five-letter source in all 625 entries after
its explicit label permutation.

The [control record](../evidence/asymmetric-four-particle-code-checks.json)
reports these exact identities separately from the numerical
cost checks.

## Reproduction

```sh
python3 tools/search_asymmetric_four_particle_codes.py --output evidence/asymmetric-four-particle-codes-q5.json
python3 tools/search_asymmetric_four_particle_shell.py --states 16 --seconds 1800 --reuse evidence/asymmetric-four-particle-codes-q5.json --output evidence/asymmetric-four-particle-shell16-q5.json
python3 tools/check_asymmetric_four_particle_codes.py
```

The construction is in
[asymmetric_four_particle_codes.py](../tools/asymmetric_four_particle_codes.py).
The complete records are the
[one-symbol schedule](../evidence/asymmetric-four-particle-codes-q5.json)
and [sixteen-state shell](../evidence/asymmetric-four-particle-shell16-q5.json).
The factoring exponent claim is unchanged.

This research and its implementation used substantial OpenAI assistant support.
