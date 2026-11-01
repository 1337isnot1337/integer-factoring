# Executable classical coprime construction for oriented ideals

## Result and scope

The [classical integer engine](../tools/classical_coprime_basis.py) now
constructs the norm basis and every sparse exponent word required by
the [linear supplied-data theorem](oriented-basis-linear-charges.md).
Its complete wrapper then attaches quadratic orientations and verifies
the resulting literal integral ideal identities. No prime factorization
or primality oracle is used by this construction.

The implementation follows D. J. Bernstein,
[Factoring into coprimes in essentially linear time](https://cr.yp.to/lineartime/dcba-20040404.pdf),
Journal of Algorithms 54 (2005), 1–30, particularly Algorithms 13.2,
18.1, and 21.2 and Theorems 18.2 and 21.3. The essentially linear
integer algorithm is classical. This checkpoint supplies an executable
reference and connects it to the repository's exact oriented output.

Complete retained 64-, 80-, and 96-bit fixtures reproduce exactly the
previously saved keys, atoms, and signed words. This is a construction
and correctness result. It neither supplies new factoring relations
nor improves the proved score exponent, which remains `1/2`.

The Python implementation is substantially slower on these fixtures
than the [prepared bucket implementation](batched-coprime-bucket-basis.md).
The active searches continue to use that practical implementation.
The classical fast-arithmetic complexity theorem and these Python
timings have different purposes and must not be identified.

## 1. Interfaces and exact output

`factor_into_coprimes(values, verify=True)` accepts a finite list of
positive integers and returns pairwise coprime integers `A_j > 1`
together with every original sparse factorization

$$v_i=\prod_j A_j^{e_{ij}},\qquad e_{ij}>0.$$

An atom `A_j` can be composite or a perfect power. The task is to
separate incompatible prime-exponent patterns across the supplied
integers, not to recover a prime factorization. For example, a common
composite factor with the same exponent pattern can remain one atom.

The implementation deduplicates equal input values for norm
construction, then restores a word for every original occurrence.
Units have empty words. Empty families are permitted. Invalid
nonpositive inputs and insufficient supplied bases are rejected.

`classical_oriented_coprime_basis(discriminant, keys, verify=True)`
accepts the same valid keys as the supplied-data orientation theorem:
negative quadratic-order discriminant `Delta`, positive odd norms
coprime to `Delta`, and roots satisfying

$$B_i^2\equiv\Delta\pmod{4R_i}.$$

It first constructs the norm factorizations of all `R_i`. It then calls
the already published [orientation lift](../tools/oriented_norm_basis_lift.py).
For canonical atom roots `C_j`, the output signed words mean

$$I(R_i,B_i)=\prod_{w_{ij}>0}I(A_j,C_j)^{w_{ij}}
\prod_{w_{ij}<0}\overline{I(A_j,C_j)}^{-w_{ij}}.$$

Negative signs select conjugate integral ideal factors in this
identity. They are not a claim that conjugation equals inversion as
an equality of integral ideals. That distinction remains necessary
when a downstream calculation asks for an actual ideal square.

## 2. Constructive integer stages

The code retains the paper's decomposition into small exact routines.
All arithmetic uses integer multiplication, exact division, gcd, and
comparison. Every division checks that its remainder is zero.

| Routine | Published algorithm | Role |
| --- | --- | --- |
| `support_split` | 11.3 | Separate the complete part supported on another integer's primes. |
| `comparison_split` | 11.4 | Separate primes according to which of two exponents is larger. |
| `pair_basis` | 13.2 | Construct the natural basis of two integers by shifted exponent ranges. |
| `product` | 14.1 | Multiply a family by a balanced product tree. |
| `split_over` | 15.3 | Associate each coprime basis member with its supported input part. |
| `extend_basis` | 16.2 | Refine a coprime family together with one additional integer. |
| `merge_bases` | 17.3 | Merge two coprime families using binary-index products. |
| `natural_basis` | 18.1 | Recursively merge the bases of balanced input halves. |
| `reduce_power` | 19.2 | Recover an exponent by recursive exponent doubling. |
| `factor_one` | 20.1 | Recover one sparse word over a supplied basis. |
| `factor_family` | 21.2 | Recover all words while restricting each recursive branch to its relevant basis. |

Two details matter for the published analysis. First, the ordered
recursive arguments in Algorithm 13.2 are preserved. Replacing them
by a repeated subtractive gcd refinement would be a different
algorithm with a different dependence on large exponents. Second,
the family factorization filters the relevant basis before dividing
the input family. It does not repeatedly scan an unrestricted global
basis for each input.

One optional shortcut precedes the ordinary merge. If the products
of the two already coprime families have gcd one, their concatenation
is already a valid merged basis. Otherwise the published merge runs.
The extra product and gcd work at balanced recursion levels is
essentially linear up to logarithmic factors in the supplied data
under the same fast-arithmetic model. This shortcut does not invoke
factoring or alter the required output.

The instrumented arithmetic object records operation counts, operand
bit totals, maximum operand length, and recursive pair depth. These
are diagnostic charges, not a substitute for a bit-complexity proof:
an operation on a large integer is not a unit-cost operation.

## 3. Complexity and its boundary

Let

$$S=\sum_i\mathrm{bitlength}(R_i)$$

count every supplied key, including repeated norms, and let `n` bound
the discriminant and individual norm bit lengths. The classical
construction and sparse factorization theorems, followed by the
orientation and validation bounds already proved in
[the linear-charge note](oriented-basis-linear-charges.md), give

$$S\,(n+\log(S+2))^{O(1)}$$

bit operations for the complete supplied-key construction. No dense
key-by-atom table is needed. The exponent incidences and the final
oriented output are sparse and have the stated input-length charges.

The wrapper verifies pairwise coprimality by product/remainder trees,
reconstructs all norm words, and independently checks the ideal HNF
identities. Norm reconstruction stops if an intermediate value exceeds
the claimed norm. The orientation step processes every original key,
even if its norm was deduplicated in the first stage.

This bound describes an algorithm with fast integer arithmetic.
Python's native integer operations and the deliberately transparent
instrumentation do not by themselves realize that asymptotic model
or a practical speed improvement. In particular, the timings below
are not evidence for an improved factoring exponent.

The unresolved factoring task remains the supply of enough useful
relations, with every generation, enumeration, and extraction cost
charged. Efficiently organizing a supplied list cannot prove that
such a list appears in fewer worst-case search steps.

## 4. Independent controls

The [checker](../tools/check_classical_coprime_basis.py) and
[complete receipt](../evidence/classical-coprime-basis.json) contain:

- All 10,000 ordered pairs with entries from 1 through 100, compared
  with an independent prime-exponent-direction oracle.
- 246 integer families and 3,233 complete sparse words, including
  duplicate values and malformed or insufficient inputs. Six invalid
  cases are rejected.
- Six exponent-ladder cases, reaching 21,175-bit supplied integers.
  The largest case has pair recursion depth 16 and 1,147 gcd calls.
- Five small quadratic ideal banks, totaling 239 literal ideal
  identities, with independent integral HNF verification.
- Complete comparisons against the retained 64-, 80-, and 96-bit
  factoring fixtures. Every key, atom, and signed word agrees with
  its saved reference; every reconstructed ideal is checked.

The small pair and family controls may use known prime factors to
construct their independent expected answers. Those factors are
never passed to the construction under test.

| Retained fixture | Supplied keys | Oriented atoms | CPU seconds | Wall seconds |
| --- | ---: | ---: | ---: | ---: |
| 64 bits | 1,806 | 1,854 | 1.981064 | 1.981143 |
| 80 bits | 9,776 | 10,079 | 27.208621 | 27.273541 |
| 96 bits | 19,705 | 20,402 | 100.607930 | 105.214634 |

These are single complete construction-and-verification observations
from the saved receipt. They are not matched end-to-end factoring
benchmarks. The prepared-bucket component's separate 96-bit matched
median is 1.732033 CPU seconds; the difference in scope and protocol
also prevents treating their ratio as a controlled speedup estimate.

The largest retained fixture uses 723,130-bit intermediate operands
in the norm engine. Those large products are explicitly instrumented.
The subsequent orientation lift performs 27,584 orientation gcds and
checks 19,705 complete ideal identities.

## 5. Reproduction and use

Run all finite controls and the default complete 64-bit fixture with:

```sh
python3 tools/check_classical_coprime_basis.py
```

To reproduce the larger saved construction comparisons:

```sh
python3 tools/check_classical_coprime_basis.py \
  --fixture-bits 64 80 96 --output /tmp/classical-coprime-replay.json
```

The receipt binds the construction, checker, orientation lift, and
underlying ideal arithmetic by SHA-256. The complete fixture hashes
are also retained. No new Lean formalization is claimed by this
checkpoint.

This implementation and analysis used substantial OpenAI assistant support.
