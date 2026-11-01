# Polynomial shifts of endpoint-zero affine directions

The constant-companion reflection catalogue gives a 1,911-dimensional
space of polynomial affine differences at alphabet four. Its vectors
have degree at most 24, zero crossing endpoint, and zero paired traces.
Multiplying a difference by `w` preserves those zero identities whenever
the resulting degree is at most 24.

The [bounded shift-closure program](../tools/close_polynomial_shift_span.py)
enlarges this space to dimension 3,674. This is an exact subsidiary
search space, not a pressure improvement. The subsequently saturated
independent-paired space is larger and already contains these directions.

## Exact closure algorithm

First reorder polynomial coordinates by degree, then by tensor entry.
Use echelon pivots in decreasing coordinate order. The basis rows whose
leading degree is at most `d` span the entire part of the current space
of degree at most `d`: any nonzero combination has the largest pivot
of its selected rows, which cannot cancel against a smaller pivot.

Multiply every row of degree at most 23 by `w`, reduce it against the
current basis, and insert any nonzero remainder. New rows of degree at
most 23 join the same queue. At termination, the space is stable under
every degree-preserving shift of every one of its elements. The queue
adds at most 6,400 independent field coordinates, so it terminates.

Every step uses linear combinations and a legal polynomial shift.
Consequently the result is the smallest subspace containing the input
and closed under this bounded operation. This statement does not assert
membership in a larger module obtained by allowing intermediate degrees
above 24 and then canceling them.

The replay performs 3,456 independent fixed-point shift reductions.
It checks every input and output vector's endpoint and both complete
paired-trace polynomial contractions, and checks both packing directions.
The [record](../evidence/polynomial-shift-affine-span.json) binds the
input and output hashes. The [compressed basis](../evidence/polynomial-shift-affine-basis.json.gz)
uses the existing `entry*25+degree` packing, so the higher-Hasse-jet
target solver can consume it directly.

```sh
python3 tools/close_polynomial_shift_span.py \
  evidence/companion-affine-polynomial-span.json.gz \
  --output evidence/polynomial-shift-affine-span.json \
  --basis-output evidence/polynomial-shift-affine-basis.json.gz
```

The [uniform low-jet degree truncation](low-jet-reflection-degree-truncations.md)
is a separate direct construction. It does not depend on this finite
closure or on solving any of its affine constraints.
