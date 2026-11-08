# Twelve-site coefficient action: word-pair expansion obstruction

## Scope

The twelve-site coefficient algebra has dimension `10`
([proof](bounded-spin-twelve-coefficient-dimension.md)), and the exact
ten-dimensional multiplication table and primary decomposition are now
available ([table](bounded-spin-twelve-coefficient-table.md),
[idempotents](bounded-spin-twelve-coefficient-idempotents.md)). The generic
`R,T` orbit has dimension `8` with an explicit `8 x 8` orbit model; the
exact identification of the `R,T` subalgebra inside the ten-dimensional
algebra remains open. This note records a computational obstruction to
computing the **action of the coefficient directions and idempotents on the
cut modules** by the existing word-pair expansion. It is a negative
diagnostic, not a no-go theorem and not a bound.

## Exact finite data

At length twelve the bounded-spin cut has `q = 89` paths, so the cut space
has dimension `q^2 = 7921`. The word-pair interaction expansion uses all
`C(12,6)^2 = 924^2 = 853776` ordered word pairs; the cut matrix is obtained
by summing local blocks. The coefficient action on a vector is a sum over
these blocks.

A direct action for coefficient index `5` (the ten-site completing
direction) was timed on worker-1 in the repository's native GF5 arithmetic:

- weight `0`: `924` pairs, `924` blocks, `0.3 s`;
- weight `2`: `33264` pairs, `16632` blocks, `72.2 s`;
- the next weights exceeded a `900 s` budget before completing.

At weight `2` alone the intermediate matrix has `7921 x 33264` entries
(`263,316,144` bytes before temporaries). Later weights have more blocks and
larger local dimensions. The earlier attempt to build the full coefficient
matrix `coefficient_matrix(12,5)` exceeded a twenty-minute budget without
output. The ten-site control has `q = 34` and only `63504` word pairs, which
is why the same construction is routine there.

## What this rules out

- The full word-pair tensor, and the naive action obtained from it, is not a
  practical route to the twelve-site cut-module action on this machine.
- Trace forms do not distinguish the coefficient directions: the tested
  coefficient sets all have trace-form rank `6`, and indices `8,9` span only
  a rank-`2` trace image. Trace quotients cannot replace the action.

## What is still open

- The exact row and column ranks of the eight primitive central idempotents
  on the `89^2=7921`-state cut modules, and their two-cut pairing.
- The exact identification of the `R,T` subalgebra inside the ten-dimensional
  coefficient algebra, and of the eight-dimensional orbit model as a quotient.
- A transfer-matrix or matrix-product-operator contraction of the local
  six-vertex model is the natural replacement for the word-pair expansion.
  A quotient that does not expand all `924^2` word pairs would also suffice.
- The complete factoring bound is unchanged; no crossing or rank partition
  is claimed here.

## Reproduction

The timing diagnostic is the `/tmp/apply12.py` prototype used during the
campaign; it is not a repository checker. The established numbers are in
`evidence/bounded-spin-twelve-RT-orbit.json` and the companion
[note](bounded-spin-twelve-coefficient-dimension.md). The peer table and
idempotent certificates are reproducible with
`make verify-bounded-spin-twelve-coefficient-table` and
`make verify-bounded-spin-twelve-coefficient-idempotents`.
