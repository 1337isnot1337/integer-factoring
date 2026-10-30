# Current crossing frontier and resolved shortcuts

The target remains a uniform leading width coefficient below `2sqrt(2)`
(or a different evaluator/reduction). The repo's strongest complete ordinary
factoring exponent is `1/2`. This ledger records which present results meet
which part of that target; it is not a list of new factoring bounds.

## Effective families with an open rank theorem

- [Bounded bit weights](hamming-pfaffian.md) and
  [fixed antipodal-pair sectors](hamming-sectors.md) have linear degree,
  complete-span optimization cost `2^O(lambda)`, and scale-one crossing.
  Their missing requirement is a uniform rank-product bound. A fixed number
  of sectors may be combined in that same budget.
- [Sparse-width interpolation](sparse-width-interpolation.md) permits an
  unbounded uniform sequence with polynomial gaps. Finitely many favorable
  examples cannot substitute for that unbounded sequence.
- [Rational-family transfer](rational-family-transfer.md) permits bounded
  rational parameters and effective specialization/descent. A symbolic
  low-rank identity still has to be constructed and proved.

The [characteristic-zero word barrier](word-symmetry-barrier.md) rules out
full position-invariant word constructions on fixed seeds of sizes 2, 3, and 4
even with collective recombinations, including products of `o(t)` position
classes. Its semisimple block argument is not a
characteristic-two theorem.

The [characteristic-two pairing-type barrier](pairing-word-barrier.md)
closes the full universal three-pairing type span for all seed dimensions
at least five, including collective recombinations. It uses explicit local
nilpotent blocks and covers extension fields. Its bound leaves dimensions
2, 3, and 4 open.

## Reliable finite checkpoints

- [Three-label source span](modular-recombination.md): exact optimum four
  over `F_2`, with a complete audit and an independent short lower proof.
- [One-pair family](hamming-pfaffian.md): full-span optima 4, 36, and 256
  at widths 1, 2, and 3.
- [Joint one/two-pair family](hamming-joint-optimum.md): exact width-two
  optimum 25 over `F_2`; its rounded width sum is still six.
- [Extension-field scope control](hamming-joint-optimum.md#a-strict-extension-field-scope-control):
  the separating row-rank lemma fails over `F_4`, with explicit ranks
  `(4,10)`. The product 40 is not an improvement.

The width-three joint space has 130 dimensions. Restricted single/pair
and old-low-span coset searches did not find a product below 256; a larger
rank-factorization probe timed out. These are exploratory limits reported
on [board #8](https://github.com/1337isnot1337/integer-factoring/issues/8),
not a lower-bound certificate for that whole space. Reopening requires a
new rank argument, a materially different search constraint, or a new
family; repeating the same bounded combinations supplies no new evidence.

## Why a finite frozen-mode sector is not invariant

A proposed shortcut is to keep a finite set of active modes and fix the
occupancy of every other retained mode. The actual homogeneous source
operators still create modes above any such finite core.

Work in the characteristic-two regular Fock quotient: distinct odd modes
have occupancy zero or one. Fix `h=2^a>=8`, and odd parts
`1<n<m=h-n`. Put `A=n-1` and `B=m-1`. Their non-low bits are disjoint,
nonempty complements below the leading `h` bit.

Let `F` be any finite occupied set and let the proposed active core also
be finite. Choose a power of two `R` larger than `h` and every retained
mode, and put

```math
p=n+R,\qquad q=m+2R.
```

Both new odd modes are outside the core and unoccupied. The pure-creation
coefficient of the source monomial operator `t_n t_m`, from `F` to
`F union {p,q}`, is the determinant over `F_2` of

```math
\begin{pmatrix}
[ A subset p-1 ] &[ A subset q-1 ]\\
[ B subset p-1 ] &[ B subset q-1 ]
\end{pmatrix}
=\begin{pmatrix}1&0\\0&1\end{pmatrix}.
```

The diagonal conditions hold because adding the high bits causes no carry
in the low positions. The off-diagonal conditions fail because `A,B` are
disjoint and nonempty. Hence this coefficient is one. An annihilation or
mixed term cannot reach this output: its mode count is larger by two, and
the parameter has only two parts. This is a nonzero source action of
weight shift `p+q` outside the proposed fixed-occupancy sector.

Therefore retaining all modes of this particular source operator is
incompatible with that finite active-core replacement. Filling finitely
many additional modes does not repair it: choose `R` larger again.
An infinitely filled tail is not a finite state of the source Fock module
and cannot be substituted without a new representation and cost proof.
This is a restricted operator-invariance obstruction, not a lower bound on
arbitrary deformed crossing families or factoring.

The creation coefficient is the source's binomial coefficient reduced by
Lucas's identity. The source recursion and the independent determinant
formula are already derived in [modular-crossing.md](modular-crossing.md),
from Claude (Anthropic)'s Appendix B at the pinned `16a5ce4` revision.
The argument here used substantial OpenAI assistant support, without new
Lean formalization. It complements the earlier
[projection counterexample](post-half-directions.md#projecting-away-high-modes-breaks-the-full-crossing-identity).

```sh
python3 tools/check_frozen_core_leak.py
make verify
```

The finite checker compares the source recursion with independent
normal-order determinants on explicit frozen states, including creation
above each cutoff. The all-size conclusion follows from the displayed
identity, not from extrapolating the checked cutoffs.
