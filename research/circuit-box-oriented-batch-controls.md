# Executable two-pass construction of supplied oriented ideals

The two-pass construction in
[circuit-box-orientation-batch.md](circuit-box-orientation-batch.md) now has an
executable implementation in
[circuit_box_oriented_batch.py](../tools/circuit_box_oriented_batch.py).
It calls the repository's classical coprime-basis engine twice. The preferred
path first takes one gcd for each nonzero orientation difference, replacing
that difference by a unitary divisor of its own parent. The second pass then
batches every orientation split. Subsequent sign extraction uses sparse
membership queries and performs no further orientation gcds. The original
raw-difference path remains an explicitly tested option.

This result concerns **supplied odd, unramified ideal keys**. It does not bound
the cost of finding those keys, obtaining enough relations, or factoring the
original input. The fast-arithmetic complexity in the proof remains distinct
from the timings of this Python implementation.

## Construction and representation

The input is a negative quadratic-order discriminant \(D\) and a list of
keys \((R_i,b_i)\). Each key satisfies

\[
R_i\ge1,\qquad R_i\text{ odd},\qquad
\gcd(R_i,D)=1,\qquad b_i^2\equiv D\pmod {4R_i}.
\]

The routine normalizes \(b_i\) modulo \(2R_i\). Distinct orientations remain
distinct input rows even when their norms agree. Unit and repeated keys are
accepted.

The first coprime-basis pass produces pairwise coprime parents \(A_j\) and
sparse norm words \(R_i=\prod_j A_j^{e_{ij}}\). A product/remainder tree
reduces each input root only at its incident parents. For one reference root
\(\beta_{0j}\) at each parent, form

\[
\delta_{ij}=
\frac{\beta_{ij}-\beta_{0j}}2\bmod A_j.
\]

For a nonzero difference, the preferred path replaces it by

\[
g_{ij}=\gcd(A_j,\delta_{ij}).
\]

The half-difference and half-sum have product zero modulo \(A_j\). No prime
dividing \(A_j\) divides both, because the root is unramified. Therefore
\(g_{ij}\) is a **unitary divisor**: \(\gcd(g_{ij},A_j/g_{ij})=1\).
It records precisely the same orientation choices as the original
difference. The second pass receives all parents and all these nonzero
preconditioned differences, including possible units.

Every second-pass input divides one parent as a unitary divisor. Its natural
coprime basis therefore consists entirely of unitary parent factors, with
no auxiliary components. To see this at each prime, write the input
valuation vector as \(e\) times a zero-one vector. The parent coordinate is
one, so the gcd of that vector's entries is exactly \(e\). The natural
basis retains the full prime power \(p^e\), grouping primes only when their
zero-one patterns agree. In particular, every parent-to-child exponent is
one. The implementation checks both conclusions.

Setting `unitary_preconditioning=False` instead sends every original
nonzero difference to the second pass, exactly as in the preceding proof
note. In this option, a component is retained exactly when it occurs in a
parent word. It then has a unique parent. Components that occur only in
differences are discarded.

If \(C_\ell^{u_{j\ell}}\) occurs in \(A_j\), the unramified root condition
ensures that a nonzero difference contains either no copy of this component
or at least \(u_{j\ell}\) copies. Membership in the difference word therefore
determines its orientation. A zero difference has the reference orientation
on every child. The final coefficient has magnitude
\(e_{ij}u_{j\ell}\); its sign also incorporates the canonical choice between
the reference root and its negative modulo \(2C_\ell\).

A positive signed coefficient denotes a positive power of the displayed
integral ideal. A negative coefficient denotes a positive power of its
**conjugate integral ideal**. The optional verifier checks the resulting
literal integral-ideal identity by Hermite normal form. It does not interpret
a negative coefficient as an integral inverse.

## Charged bounds

Let \(F=\sum_i\lfloor\log_2 R_i\rfloor\). The implementation checks

\[
\text{second-pass input bits}\le4F,
\qquad
\text{final sparse nonzeros}\le F.
\]

The unitary preconditioning also has total gcd operand length at most
\(4F\). Indeed, for incidence \((i,j)\), both operands have at most
\(\mathrm{bits}(A_j)\) bits, while
\(\sum_{j:e_{ij}>0}\lfloor\log_2 A_j\rfloor
\le\lfloor\log_2 R_i\rfloor\). Fast gcd arithmetic therefore gives
\(O(S\,\mathrm{polylog}(S))\) bit cost for all preconditioning, where
\(S=\sum_i\mathrm{bits}(R_i)\). There is no extra factor equal to the
largest parent length or the number of supplied keys.

The proof note gives the incidence argument behind these bounds. Every
product/remainder-tree multiplication and remainder has an arithmetic count
and an operand-bit count. The receipt also includes the complete counters
from both classical coprime-basis calls. In the raw-difference option,
auxiliary components are charged in the second pass even though they are
absent from the output.

The classical integer engine follows Bernstein's coprime-basis algorithms;
the relevant source is [Factoring into coprimes in essentially linear
time](https://cr.yp.to/lineartime/dcba-20040404.pdf), especially Theorems 18.2
and 21.3. With the input-size reduction proved here, those two passes have
essentially linear integer-data cost. Root validation, form reduction, sparse
assembly, and optional literal-ideal verification have polynomial overhead
in the discriminant and norm bit sizes. No factorization oracle or residual
primality test is used by the construction.

## Independent controls

Run the complete deterministic suite with:

```bash
python tools/check_circuit_box_oriented_batch.py \
  --output evidence/circuit-box-oriented-batch.json
```

The saved [receipt](../evidence/circuit-box-oriented-batch.json) records:

| Control | Result |
| --- | ---: |
| Ordered batches over five quadratic discriminants | 15 |
| Literal integral-ideal identities | 2,932 |
| Independent prime-expansion comparisons | 3,357 |
| Comparisons with the single-pass oriented constructor | 15 |
| Unitary/raw-difference variant comparisons | 24 |
| Invalid inputs rejected | 5 |
| Mixed prime-power orientations | 256 keys |
| Largest norm in the mixed family | 186 bits |

The small batches contain every valid odd unramified key of norm below 140
for \(D=-4,-23,-84,-260,-575\), in forward, reverse, and shuffled orders.
The shuffled order includes repeated rows. Independent trial factorization
is used only by these small control cases, not by the construction.

The mixed family uses all 256 sign patterns on eight explicitly supplied
primes, with norm exponents from one through five. The first pass keeps one
coarse parent; both variants produce the same eight retained atoms and
2,048 output nonzeros. Its 28,598 input norm bits have the following charged
second-pass costs:

| Quantity | Preferred unitary path | Raw differences |
| --- | ---: | ---: |
| Preconditioning gcd calls | 255 | 0 |
| Second-pass input bits | 4,890 | 9,327 |
| Second-pass gcd calls | 62,442 | 209,237 |
| Second-pass multiplication operand bits | 586,534 | 9,120,891 |
| Discarded auxiliary components | 0 | 212 |

These are arithmetic counts on one fixed control, not an asymptotic speedup
claim or a matched wall-clock benchmark. The raw path's discarded
components exercise the relevance filter and remain visible in its counters.

The suite also checks the **first 128 keys** of each retained 80-bit and
96-bit factoring fixture. Those source fixtures contain 9,776 and 19,705
keys, respectively; this receipt does not claim that their full key sets were
replayed. The two checked prefixes produce 137 and 129 atoms. Source hashes,
prefix lengths, and arithmetic counters are included in the receipt.

Additional cases cover empty input, repeated units, two orientations of 65,
unequal valuations in a split of 325, and repeated powers of 13. A raw
difference may have strictly greater valuation than its parent component;
this is accepted and tested. Hashes of every implementation dependency are checked
both before and after the suite to detect concurrent source changes.

The raw atom lists need not agree. The explicit input
\((29,24),(29,34),(65,16)\) at discriminant \(-4\) leaves atoms 29 and 65
in the preferred path. The raw difference 5 from parent 29 additionally
splits the unrelated parent 65, producing atoms 5, 13, and 29. Both outputs
are checked. Comparisons therefore use independent prime expansions and
literal ideal identities, which are the invariant outputs relevant to
relation processing.
