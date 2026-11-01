# Exact quotient updates for the oriented atom product tree

## Result

The [new tree variant](../tools/quotient_update_oriented_basis.py) preserves
every product, overlap query, oriented atom, and signed usage word from
the [sparse incremental basis](incremental-oriented-coprime-basis.md).
It changes only the arithmetic used to update ancestor products after
a leaf insertion or removal.

On the retained 80-bit prefix, three matched construction-and-export runs
reduce the median component wall time from 17.864 seconds to 1.171 seconds,
about a factor of 15.2. All 9776 supplied ideal identities and the complete
mathematical output agree. This is a component timing improvement; the
whole factor search is measured separately.

## 1. Exact update identity

Suppose a leaf changes from value d to value e. Every ancestor stores a
product P containing that leaf exactly once. Its new value is

$$P'=\frac{P}{d}e. \tag{1}$$

The quotient is integral. No modular inverse, factorization, approximation,
or cancellation assumption is needed. The implementation uses exact
integer division and asserts a zero remainder.

An insertion replaces a unit vacancy, so `d=1` and the update is simply
`P'=Pe`. A removal sets `e=1`, so the update is `P'=P/d`. These are the only
two updates made by the current atom basis. Capacity doubling retains the
original balanced product-tree rebuild.

Induction from the leaf to the root proves that every stored product is
identical to the original child-product update. Therefore every later gcd
query visits exactly the same nodes and returns the same overlapping
atoms. Local refinement, root orientations, sparse usage propagation, and
the literal ideal identities are unchanged.

The variant subclasses the original tree and basis. The original modules
remain available for exact matched comparisons. It exports
`QuotientUpdateOrientedBasis` and the convenience function
`quotient_update_oriented_coprime_basis` with the same sparse output
contract and optional dense compatibility output.

## 2. Charged arithmetic cost

The earlier implementation recomputes each changed ancestor by multiplying
its two child products. Near the root, both children can be large even
when the changed leaf is small. Equation (1) instead divides or multiplies
the current large product by the leaf value.

For an L-bit ancestor and a b-bit leaf, the update costs one exact division
or one multiplication with operand sizes `(L,b)`, as applicable. Elementary
binary arithmetic bounds either operation by `O(Lb)` bit operations.
Fast unbalanced integer arithmetic can have a different cost; the bound
does not treat an update or a division as a unit-cost operation.

The height, number of leaf updates, gcd decisions, and storage remain the
same. Rebuild multiplications at capacity doubling are still charged.
The receipts add exact-division counts, division operand-bit totals, and
maximum update multiplier and divisor sizes. The method does not remove
the cost of repeated large-product gcds or establish an essentially-linear
bound for the entire incremental sequence.

## 3. Exact controls and retained measurement

The [checker](../tools/check_quotient_update_oriented_basis.py) compares
every tree node after 96 insertions and removals, including six capacity
doublings and reused slots. It compares 12192 overlap queries with the
original tree. Eight ideal batches in two input orders replay 332 literal
ideal identities and return identical atom words.

The larger fixture is the complete 9776-key basis from
[the certified 80-bit power prefix](sparse-oriented-power-prefixes.md).
Both methods are rebuilt from empty on the same ordered keys. Verification
occurs after the timed construction and sparse export. Run order alternates.

| Quantity | Child-product updates | Exact quotient updates |
| --- | ---: | ---: |
| Median wall time, three runs | 17.863553 s | 1.171483 s |
| Median CPU time | 17.681335 s | 1.168561 s |
| Gcd decisions, tree plus local | 172726 | 172726 |
| Tree leaf updates | 11587 | 11587 |
| Capacity doublings | 14 | 14 |
| Tree multiplications | 175907 | 166691 |
| Exact tree divisions | 0 | 9216 |
| Final atom count | 10079 | 10079 |
| Sparse usage entries | 15045 | 15045 |

The final root product has 271098 bits. Every divisor and update multiplier
in the quotient variant has at most 40 bits. The multiplication operand-bit
sum remains above 2.4 billion: the improvement comes from the unequal
operand sizes, not from pretending those large products disappear.

The full atom-and-word payload has SHA-256
`16dd0ed0f9a8a08cff09e06dea363b790cfc3a5b3cab821322c670ee45fecea9`
under the checker's canonical encoding. Independent HNF replay confirms
all 9776 literal ideal identities. The
[receipt](../evidence/quotient-update-oriented-basis.json) retains all
timings, input hashes, arithmetic counters, and source hashes.

```sh
python3 tools/check_quotient_update_oriented_basis.py \
  --output evidence/quotient-update-oriented-basis.json
```

This implementation and analysis used substantial OpenAI assistant support.
