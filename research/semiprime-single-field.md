# One primitive field component for promised semiprimes

## Result and promise

For `N=p r` with distinct primes, one primitive field component at every
coprime prime `q>=7`, `q=3 mod 4`, determines the unordered signed factor
residues. No roughness, inverse-orbit or sparse-distance condition is needed.
Polynomially many supplied components recover the promised semiprime with
polynomial non-provider cost. The original-index coefficient engine remains
open. No arbitrary-input reduction or stronger complete exponent is claimed;
the complete score stays `1/2`.

This changes the arithmetic input promise of the
[faithful-character decoder](sparse-faithful-character.md). The next target
is the restricted coefficient engine, rather than further decoder searches.

## Quadratic recovery in one field

Put `u=(q-1)/2`. Choose a primitive generator g modulo q and any irreducible
factor f of `Phi_u` over `F_2`. In `F=F_2[T]/f`, omega=T has order u. Write

```math
\eta(a)=\omega^{\log_g(a)\bmod u},\qquad
C_N=\sum_{d\mid N}\eta(d),\qquad c=\eta(N).
```

Only `C_N` is supplied by the missing engine. The residue of N gives c by
a q-element discrete-log table. The four divisors give
`C_N=1+eta(p)+eta(r)+c`. Thus `s=C_N+1+c` and c are the sum and product
of the factor-character values. They are exactly the roots of

```math
X^2+sX+c.
```

Scan the u powers of omega. Keep a unique root twice: distinct integer
primes can have equal character values. Eta's kernel is exactly `{1,-1}`,
so a root exponent j determines a factor up to sign modulo q, and its
square uniquely as `g^(2j)`. The two squared residues retain multiplicity.
The complete cyclic divisor word is also determined: for root exponents
j,k it is `1+X^j+X^k+X^(j+k)` modulo `X^u-1`, with binary cancellations.

## Symmetric CRT and complete additional cost

At each level add the decoded squares to get `v_q=p^2+r^2 mod q`. If the
product M of distinct levels exceeds `N^2+1`, CRT uniquely lifts
`V=p^2+r^2`, because `N^2+1-p^2-r^2=(p^2-1)(r^2-1)>0`. Now
`(p+r)^2=V+2N`; integer square roots recover the sum and then the difference
from `(p+r)^2-4N`. There is no cross-level root matching. Both square
identities and the final proper-factor product are checked.

Primality, integer-square and level-gcd preprocessing handle exceptional
promised inputs. A square semiprime is split by its integer root. A level
dividing a composite input yields a proper factor unless the input equals
that prime, handled by primality. These steps do not certify a semiprime
promise for arbitrary N.

Trial factoring of q-1, generator search, cyclotomic construction, the
existing deterministic Berlekamp field acquisition, the u-element root scan,
and conversion cost `poly(q,log N)` bit operations. Field degree is at most u;
complete element words are charged. CRT and integer roots have polynomial
cost in the complete input and level encoding length.

The [existing polynomial level acquisition](sparse-faithful-character.md#42-polynomial-size-levels-by-a-classical-lower-sieve)
can be reused at constant sparsity parameter four. It supplies polynomial-size
distinct q=3 mod 4 levels. Its extra roughness is harmless and not needed by
this decoder. Taking `2b+O(1)` levels for a b-bit input ensures the product
threshold after removing exceptional levels; a level gcd can terminate sooner.
No safe-prime assumption or favorable-level search is introduced.

A uniform one-field coefficient engine costing `(b+q)^O(1) F(b)` would
therefore give a promised-semiprime algorithm costing `poly(b)(F(b)+1)`.
This preserves any fixed logarithmic-saving exponent supplied by the engine.
It neither constructs that engine nor extends its guarantee to all integers.
The parent's inversion-normalized weight-one target can supply the same
input by multiplying its field value by the known normalizing unit.

## Scope, controls and credit

The parent's q503 reciprocal sparse collision was not asserted to be actual
divisor data. The quadratic explains why it does not contradict this
semiprime-plus-known-product promise. At q7, q19 and q31 the new proof works
when the preceding sparse-distance certificate is unavailable. It does not
withdraw that certificate for arbitrary sparse inputs.

The readout is elementary quadratic and symmetric-moment algebra. Compare
the [Sanftenberg factor-moment source audit](riemann-gaussian-semiprime-source-audit.md#second-reuse-priority-factor-moments-and-honest-count-recovery)
and [existing moment recovery](divisor-moment-recovery.md). No novelty claim
is made for sum/product reconstruction. This specialization and implementation
were developed with OpenAI Codex assistance; no Lean or full independent
rederivation is claimed.

An earlier Eisenstein sum/product recovery argument appears in Peter Bruin,
[Computing coefficients of modular forms](https://pub.math.leidenuniv.nl/~pbruin/coefficients.pdf),
Remark 1.3. The coefficient algorithm in its Theorem 1.1 requires a
factored index, as recorded in the [source-contract audit](weight-one-divisor-coefficients.md#4-encoding-and-arithmetic-are-charged).

`make verify-semiprime-single-field` runs the
[decoder controls](../tools/check_semiprime_single_field.py): 486 literal-divisor
fixtures at q7,11,19,23,31,47,503, including 35 equal-root cases, three complete
CRT reconstructions, and input guards. The
[decoder](../tools/semiprime_single_field.py) never constructs these fixtures
or lists candidate factors of N. The [receipt](../evidence/semiprime-single-field.json)
pins all sources. Fixture construction is not the missing coefficient engine.

The [full default replay](../evidence/semiprime-single-field-full-verification.json)
passed at clean frozen `35b0d7492858f5d27d2d7793a7dfbed45c68f299` in
800.82 seconds. Complete compressed stdout, source/tree and log/archive
hashes are retained. Later prose/archive changes do not claim another replay.

Sol worker 2 later supplied a [scoped ordinary inspection](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19360111)
at the frozen source above, finding no gap in the quadratic/CRT and charged
cost steps without rederiving all imported field constructors. The unchanged
[independent replay](../evidence/semiprime-field-independent-review/replay.py)
passes 98 integer-cyclotomic controls across every q47/q503 field and two
generators, including 18 equal-root cases. Its
[provenance](../evidence/semiprime-field-independent-review/provenance.json)
preserves the original hashes. Run it from the repository root with
`python3 evidence/semiprime-field-independent-review/replay.py`. It writes
its result to the original temporary path and supplies no fast coefficient
engine. This later archival replay does not change the full-source pin.
