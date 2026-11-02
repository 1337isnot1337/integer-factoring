# One half-width trace and a modular-inverse interval test

## Research question and status

Base `af12893`, following the [two-input Witt-trace construction](witt-trace-native-word-circuits.md).
Can arithmetic determine the cofactor from the smaller candidate instead
of independently summing two native words? Yes: a division-free inverse
polynomial and a certified no-wrap range replace the second summation input
and product-equality test. All pointwise bit costs and exact divisions are
charged below.

This is a **complete conditional reduction**, not a faster contractor.
The represented query handles every odd input, and an even-input adapter
returns 2 immediately. Its uniform contraction remains open. The complete
proven factoring score stays **1/2**. The existing
[pointwise-selection budget](trace-selection-weil-budget.md),
[inverse-window audit](circuit-box-inverse-windows.md) and
[interval-weighted Kloosterman blocks](circuit-box-kloosterman-blocks.md)
are applicable constraints; no escape from them is asserted.

The mechanism changes the arithmetic representation, rather than a field
cutoff or search parameter. The immediate continuation question is whether
this particular unary modular-inverse predicate has an evaluator that
preserves its structure across the explicitly admitted exact divisions.
The present polynomial carry compiler does not provide one.

The later [trace-dual coordinate construction](trace-dual-native-counts.md)
keeps this inverse/no-wrap test and replaces the nonuniform word map by
a bijection onto odd candidates. It uses several coordinate traces of one
Teichmuller input, permits trial-certified smaller unweighted count
precision, and replaces inverse-bit extraction by polynomial Boolean gates
with no exact halvings. It supplies no fast circuit contractor and does
not change this note's original fibre-weighted query or evidence.

## 1. Determine a cofactor from one smaller candidate

Let N>=3 be odd, b=bits(N), and put

```
a=ceil(b/2),      k=a+b,      M=2^k.
```

For odd x with 0<=x<2^a define y to be the least nonnegative residue
of N*x^(-1) modulo M. Then

```
x>=2 and 2<=y<2^b   <=>   x>=2, y>=2 and xy=N.          (1)
```

Indeed the left side implies xy=N modulo M. Both xy and N lie in
[0,M), so they are equal as integers. Conversely, any proper factor x in
this candidate range is odd; its integer cofactor is between 2 and 2^b-1
and equals the specified inverse residue. The x=1 and y=1 exclusions
remove both trivial decompositions.

Every odd composite has a prime divisor at most sqrt(N)<2^a. Therefore
the candidate restriction preserves existence and the smallest prime
factor. It is harmless if the range also contains larger proper divisors.
For a prefix query, prefix bits describe this a-bit candidate word, most
significant bit first.

The modulus matters. With N=7, x=3 and modulus only 2^b=8, the inverse
residue is y=5. It passes the cofactor range while 3*5=15 differs from N.
Our k=a+b range excludes this wrapped-product false positive. An
inverse-interval test without a no-wrap proof is insufficient.

## 2. A polynomial inverse with no variable division

At available scalar precision h>=1, set l=ceil(h/3) and use

```
J_h(X)=X sum_(j=0)^(l-1) (1-X^2)^j.
X J_h(X)=1-(1-X^2)^l.                                  (2)
```

For every odd integer x, x^2=1 modulo 8. Hence J_h(x)=x^(-1)
modulo 2^h. The geometric sum has ordinary degree 2l-1 and O(h)
arithmetic gates, acquired uniformly by Horner evaluation; it uses no
division by X or by a nonunit. It is defined on even x too. Those inputs
are masked by the candidate's least significant bit, so no inverse claim
is made for them. For h=1,2 the same identity remains valid.

This elementary finite geometric identity is also compatible with the
existing [collected unit-inverse interface](precision-efficient-residue-indicators.md).
The new difficulty is the canonical word and range selection surrounding
the inverse, not pointwise unit inversion.

## 3. One certified trace input and the full precision invoice

Retain the preceding Witt-trace theorem at word width a. Put

```
d=4a+1,       K=d+1,       P=K+a+k=6a+b+2.              (3)
```

Work in W_P(F_(2^d)) and sum over **one** Teichmuller input U. The scalar
trace Tr_W(U), reduced modulo 2^a, attains every a-bit word. This follows
from the predecessor's Artin--Schreier--Witt twist proof with genus
g_a=(4^a-3*2^a+2)/6 and the strict certificate 2^d>4*g_a^2. No new
surjectivity or equal-fibre assertion is needed here.

Use its acquired parity polynomial H_j and known exact-division gate

```
e=H_j(t),       t'=(t-e)/2,       precision j -> j-1.
```

First apply a such gates to the trace. Reconstruct
x=sum_(i=0)^(a-1) 2^i e_i as a scalar at remaining precision
h=P-a=K+k. This reconstruction is essential: the untruncated trace
can have nonzero higher digits, which cannot be fed directly to the
inverse in (1).

Evaluate N*J_h(x), then extract its first k bits by k more known exact
divisions, leaving precision exactly K. Write those bits f_i. The product
of these flags is the pointwise predicate:

```
e_0 (1-prod_(i=1)^(a-1)(1-e_i))          # odd x>=2
* prod_(i=b)^(k-1)(1-f_i)               # y<2^b
* (1-prod_(i=1)^(b-1)(1-f_i))           # y>=2
* specified candidate-prefix flags.
```

Empty products are one. Thus the a=1 case correctly accepts no proper
odd candidate. Equation (1) proves that the predicate A_(N,prefix)(U)
is exactly the proper factor/prefix indicator, modulo 2^K. A total of
a+k=b+2a precision bits is spent on known exact divisions. All other
divisions are by fixed integers during polynomial coefficient acquisition,
as already charged in the predecessor; (2) adds none.

The coefficient ring has degree d=O(b) and precision P=O(b), hence
O(b^2) bits per represented ring element. Trace acquisition, all Horner
selectors, native-word reconstruction, polynomial inversion and mask
products have polynomial construction and pointwise bit cost. For example
the a+k dense parity Horner loops together have at most
P(P+1)-K(K+1)=O(b^2) scalar multiply/add iterations. School integer
arithmetic already suffices. The existing deterministic fixed-characteristic
field-model construction and unramified-ring acquisition are polynomial
in b, as cited in the predecessor; the finite checker supplies its small
model explicitly rather than implementing that general constructor.

Compared with its two b-bit inputs, this interface has one a-bit input,
field degree 4a+1 instead of 4b+1, and precision 6a+b+2 instead of
11b+4. These are resource invoices for different exact predicates. They
do not constitute a runtime comparison between fast domain evaluators,
since neither evaluator has been supplied.

## 4. Exact count, complete consumer and remaining evaluator

Ask for

```
C_(N,prefix)=sum_(U in T_d) A_(N,prefix)(U) mod 2^K.      (4)
```

There are 2^d summands, each zero or one, and K=d+1. The returned
residue is therefore the exact nonnegative count. Each valid native
candidate has a nonempty trace fibre, so C>0 is equivalent to a matching
proper divisor. Fibre multiplicities need not be equal or known.

On input N>=2, return 2 if N is even. Otherwise query the empty prefix.
If its count is zero, N is prime and is returned. If positive, make a
successive queries, choosing zero when the corresponding child prefix
has a positive count and choosing one otherwise. The result is the
least accepted proper divisor, hence the smallest prime factor. This
handles primes, prime powers and all other composites without a semiprime
promise or a separate primality oracle.

**Conditional transfer.** A uniform evaluator for this exact family (4)
at F(b)*poly(b) bit cost supplies an all-input smallest-prime-factor
algorithm at the same cost. There are at most a+1 oracle calls;
constructing the ring, circuit and prefixes and decoding the result
costs poly(b). If F(b)=exp(O(b/log^alpha b)), its score is alpha.
Every input cost is measured against the original b, not an unspecified
polynomial of the circuit description length. Even correct zero/nonzero
answers suffice for this consumer.

No evaluator is proved here. In particular the canonical reconstruction
and inverse range bits use known exact divisions *before* summation.
Replacing them by division of the final sum, using the unit-inverse
compiler alone, assuming uniform trace fibres, or dropping the outer
ratio-class cost in the Kloosterman algorithm is unjustified. The
ordinary polynomial expansion may have exponential degree. Fewer
summation inputs do not by themselves eliminate this obstruction.

## 5. Acquired controls, credit and continuation criterion

Run `make verify-one-trace-inverse-words`. The
[pointwise implementation](../tools/one_trace_inverse_words.py) uses the
actual predecessor's parity Horner and exact-division gates and the
polynomial (2), not an inverse oracle. Its conditional consumer takes an
explicit count callback and contains no domain-summation implementation.
The [independent checker](../tools/check_one_trace_inverse_words.py)
uses literal integer divisibility and an independently computed modular
inverse for its comparisons. Its [receipt](../evidence/one-trace-inverse-word-controls.json)
records 511 inverse-polynomial controls, 50,694 scalar pointwise queries,
3,423 prefix-count queries and 132 scalar minimum-factor recoveries.

One complete degree-nine, precision-18 ring model enumerates all 512
Teichmuller inputs. Its measured two-bit word multiplicities are
136,136,120,120; they are deliberately not replaced by uniform weights.
All proper-factor/prefix counts for six targets match the literal weighted
answer and recover the true minimum factor. Scalar checks also exercise
40-,60- and 128-bit supplied factor pairs without enumerating their fields.
The wrapped-product false positive, even-input zero-oracle-call adapter,
high trace digits and rejected uncertified field parameters are retained.
These are finite controls for the stated reduction, not a uniform count
algorithm or independent mathematical review.

The [full default receipt](../evidence/one-trace-inverse-full-verification.json)
and [compressed log](../evidence/one-trace-inverse-full-verification.stdout.txt.xz)
pin an exit-zero `make verify` at clean source
`2a09eb7a75761fc249cace3dd0c311b70085fc64`, in 961.27 seconds. This
includes the new focused target. Later source integration, including
PR90's connected-event work, is outside that full-replay pin. The unary
module, checker and imported source hashes are unchanged by integration.
No Lean target was run.

Credit: trace surjectivity, parity-polynomial acquisition, exact-division
semantics and the Galois-ring implementation retain the preceding authors'
and Davis--Wan--Xiao's credit via the linked proof. Modular inversion,
the smaller-factor restriction and binary prefix recovery are classical.
This unary composition and its controls were developed by Sol worker 2
with substantial OpenAI Codex assistance. No new Lean work is included.

Continue only with a contraction mechanism preserving the native inverse
and interval structure, or a different representation with a complete
cost invoice. Direct low-degree pointwise trace compilation remains
subject to the existing selector budget; another local sweep is not a
continuation of this question.
