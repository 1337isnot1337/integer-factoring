# Arithmetic presentation search: a framework shift

## Objective and current status

The owner requested a change of framework after the crossing/evaluator
campaign reached a wall. [PR #40](https://github.com/1337isnot1337/integer-factoring/pull/40)
retains that campaign's proofs and exact controls. This branch instead
studies deterministic arithmetic relation generation in quadratic class
groups. The strongest complete factoring score remains 1/2.

The components proved here are a positive Fourier translate-hitting lemma
that produces bounded full-rank relation rows given a pair extractor, and
a deterministic black-box two-group basis algorithm using only an
annihilator. A powering seed plus published expander amplification also
supplies biased words without knowing the class-group decomposition.
The arithmetic primitives that would turn these into a faster
factoring machine are still open. This is not a claimed breakthrough bound.

## Primary-source audit

The baseline is H. W. Lenstra Jr. and C. Pomerance, *A rigorous time bound
for factoring integers*, JAMS 5 (1992),483-516. The
[publisher copy](https://community.ams.org/journals/jams/1992-05-03/S0894-0347-1992-1137100-0/S0894-0347-1992-1137100-0.pdf)
and [retrieved mirror](https://www.cs.umd.edu/~gasarch/TOPICS/factoring/rigorous.pdf)
identify the source. The retrieved 34-page PDF has SHA-256
`ebf5d89ed9c89f0e4c51abbc2b49be8e52d2c9326dc8d2b347f533ceaef51e14`.
Relevant scanned pages were inspected visually; OCR was used only as a
navigation aid.

The source's L_N[1/2,1+o(1)] expected time is exp(O(sqrt(b log b))) for
b=bits(N). Its first parameter 1/2 is a different scale from the speedrun's
score 1/2. A deterministic algorithm at that arithmetic scale would beat
every fixed logarithmic-saving score.

| Source location | What actually requires attention |
| --- | --- |
| Theorem 4.1, p.495; Algorithm 4.4, p.496 | Prime-form generators are sampled up to exp(c(log|Delta|)^2). Enumerating that interval is not a cheap derandomization. |
| Section 2.7, pp.488-489 | Constructing a prime form also uses modular square-root computation; the source uses a probabilistic finite-field routine. |
| Algorithm 10.1, p.506 | Relation exponent vectors are random. Smoothness density is only one part of their role. |
| Algorithm 10.1 and proof, pp.507-509 | Nontrivial ambiguous output uses uniformity among hidden lattice cosets modulo 2L; many valid relations alone do not prove it. |
| Sections 7 and 10 | Smooth recognition and sparse linear algebra use randomized subroutines. Bounded trial division and dense elimination can remove these at an L_N[1/2,O(1)] cost, without addressing the harder sampling issues. |

In particular, the paper does not prove that a subexponential factor base
generates the whole class group, or that deterministic smooth-relation
selection has the needed parity distribution.

## A different sufficient contract

For an odd composite non-prime-power N, it suffices to produce, for one
valid bounded multiplier d coprime to N:

1. Actual generators of the negative-discriminant group C_(-dN).
2. A square full-rank integer relation matrix on those generators.

Its determinant is an annihilator. The
[deterministic two-part algorithm](arithmetic-torsion-extraction.md)
then recovers a basis of the ambiguous subgroup, from which a proper
factor follows. No exact class number or saturated presentation is required.

For an asymptotic improvement the generator count, relation bit lengths,
construction cost and group operations all have to satisfy the requested
uniform bit bound. After trial division through 2^ceil(sqrt(b)), the
two-rank is O(sqrt(b)), so generic torsion extraction costs only
2^O(sqrt(b)) times its polynomial input costs. For example, if generator
and relation construction take 2^O(b/(log b)^beta) with beta>1/2, while
their encoded sizes are 2^O(sqrt(b log b)), the remaining extraction and
recursion preserve score beta. This is a conditional transfer, not a
supplied generator/relation construction.

## Concrete positive-kernel program

The [positive relation kernel](positive-relation-kernel.md) gives the
following deterministic construction from specified primitives:

- A set S of small smooth reduced forms with a lower density bound and
  short factor-base words.
- An explicit low-bias distribution on the class group, supplied as words
  in the generators; the [biased-word construction](arithmetic-biased-words.md)
  now supplies this from generators and a group-size upper bound.
- A procedure that extracts a smooth pair from the guaranteed translated
  correlation intersection.

For each generator, choose a sufficiently large power, extract a pair,
and form a relation row with a dominant diagonal coordinate. All rows are
independent by construction. This replaces random relation-span growth
with a checked integer determinant. The two-group algorithm separately
replaces random coverage of ambiguous forms.

The open arithmetic task is the efficient, uniform construction of the
generators and smooth pair witnesses. The finite kernel
checker explicitly uses exhaustive pair search and does not count that
as an efficient extraction algorithm. This separation is the central
research contract, rather than a hidden oracle assumption.

## Current evidence and next action

`python3 tools/check_arithmetic_two_group.py` verifies 703 product-group
controls, including mixed orders, odd components, redundant generators
and oversized annihilators. `python3 tools/check_positive_relation_kernel.py`
checks the exact Fourier lemma and a complete finite presentation/torsion
pipeline with an unsaturated relation lattice. Neither is evidence of a
new asymptotic factoring bound.
The [actual class-group controls](arithmetic-class-group-controls.md) also
check 40 discriminants, all associative triples, smooth prime-form words,
and proper-factor extraction for 18 targets. Those controls enumerate
classes and orders; they do not implement fast generator/relation search.
The [norm-sieve experiment](arithmetic-norm-sieve.md) goes further in a
finite run: it searches arithmetic norm witnesses for 88 dominant rows,
then factors a 60-bit input by determinant/torsion extraction without
class enumeration. Its fixed-box success has no uniform guarantee.

The [bounded smooth-norm interface](bounded-smooth-norm.md) now proves the
exact arithmetic characterization and implements conversion of supplied
witnesses. Its exact divisor step costs 2^O(b/log b); locating a witness
quickly remains open. The next work is deterministic arithmetic search
within that interface. The prime-form generator requirement remains a
separate obligation for complete two-torsion recovery. The
[density/index transfer](arithmetic-density-lift.md) gives an alternative:
relations on the possibly proper small factor-base subgroup suffice for a
global class-number multiple. That removes full generation from this
intermediate stage, but still needs deterministic splitting or coverage of
the missing ambiguous directions. Repeating a crossing seed search is not part of this
branch's active scope.

This work uses substantial OpenAI Codex assistance. It adds ordinary
proofs, exact arithmetic and an explicit research contract, with no new
Lean formalization.
