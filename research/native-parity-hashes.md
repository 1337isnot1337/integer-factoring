# One-bit native queries suffice with direct hash masks

## Result and provider contract

The [trace-dual native circuit](trace-dual-native-counts.md) supports a
complete all-input factoring reduction at **coefficient precision one**.
Use its actual Boolean candidate bits to test hash buckets directly,
then apply the existing [bit-tagged isolation theorem](bit-tagged-witness-recovery.md#5-odd-characteristic-allows-a-smaller-hash-family-in-polynomial-space).
There is no Fourier projection or coefficient-algebra enlargement in this
composition. The provider must evaluate the specified **masked** queries;
one unmasked parity by itself does not establish existence of a factor.

This is an explicit restricted provider interface and implemented consumer.
No fast provider or stronger complete factoring exponent is supplied.
The established exponent remains **1/2**. The polynomial circuit may have
exponentially large expanded degree and support.

For odd N>=5 let `b=bits(N)`, `a=ceil(b/2)`, `d=a-1` and use the acquired
field `F_(2^d)`. Its trace-dual coordinate bits e_i(U) are the d bits of
`w=(x-1)/2`, where x ranges bijectively over odd candidates below `2^a`.
At precision one, their parity polynomial is just H_1(T)=T. Thus all
coordinate, bit-serial cofactor and mask operations are ordinary finite-field
arithmetic. The predecessor proves that its unprefixed Boolean mask A_N(U)
is exactly `[1<x<N and x|N]` on that candidate domain.

Supply binary-packed control constants `v_i in {0,1}^h`, a bucket
`z in {0,1}^h`, and either no tag or an index j in `[0,d)`. The new provider
returns the scalar bit

```math
\sum_{U\in\mathbb F_{2^d}} A_N(U)\,
 \prod_{r=0}^{h-1}\left(1+z_r+\sum_{i=0}^{d-1}(v_i)_r e_i(U)\right)
 \begin{cases}1,&\text{no tag},\\e_j(U),&\text{tag }j.\end{cases}       \tag{1}
```

Everything in (1) is evaluated in characteristic two. Each factor tests
one bit of the binary linear hash. Coefficients and input encodings have
polynomial size in b,h. This is a family with one native field input and
one output bit, with an explicit polynomial-size pointwise circuit.
The contract charges its uniform global evaluation against the original b.

## 1. Reuse polynomial hashing and verify every decoded candidate

After completed trial division through `2^t`, the true candidate-word set
S has size at most

```math
D=2^{\min\{d,\lceil b/t\rceil\}}.
```

Acquire a binary **control** field `H=F_(2^h)` at the least positive h with

```math
Q=2^h>(D-1)\lceil d/h\rceil.
```

This is the field and chunk-polynomial hash already used in the credited
bit-tagged construction. Here it supplies only binary control constants;
its field is not embedded into the native coefficient field.
Split a d-bit word into h-bit chunks, padding the last by zero, and set

```math
W_w(T)=\sum_{r=0}^{\lceil d/h\rceil-1}\operatorname{chunk}_r(w)T^r.
```

For each theta in H, the hash `W_w(theta)` is binary linear in w, with
acquired basis values v_i. For fixed w and each distinct w', the difference
polynomial is nonzero and has at most `ceil(d/h)-1` roots. The union of
these collision sets has size below Q. Thus every true witness is alone
in a bucket for at least one theta. This is the predecessor's ordinary
isolation argument, reused without a new general hashing claim.

For every theta and bucket request (1) without a tag. When it is one,
request each of its d bit tags and decode their bit vector as a word w.
Check the actual integer candidate `x=1+2w`: retain it only if `1<x<N`
and `x|N`. A multi-witness odd bucket can have misleading XOR tags;
verification is essential. Every accepted candidate is a true divisor,
and each true candidate is eventually isolated and decoded. Keep only
the smallest accepted candidate. If none is found, the input is prime,
because every odd composite has a prime factor below `2^a`.

The [implemented consumer](../tools/native_parity_hashes.py) performs the
trial pretest before using D. It returns two on even inputs and handles
N=3 separately. This is a complete all-input reduction, with no semiprime,
squarefreeness or independent-fibre promise.

## 2. Complete cost and a concrete coefficient target

The control size Q is `poly(b,D)`; for instance a control degree exceeding
`log_2(Dd)` suffices. The imported exhaustive small-field acquisition is
polynomial in Q and hence polynomial in b,D. The native degree-d field
retains the predecessor's polynomial-time fixed-characteristic construction.
There are at most `Q^2(d+1)` scalar queries, each with polynomial-size
hash coefficients and pointwise circuit. Recovery can stream its queries
and retain one candidate.

If (1) is uniformly evaluated in `F(b)poly(b,h)` bit operations, the
complete smallest-factor algorithm takes

```math
(2^t+F(b))\operatorname{poly}(b,2^{\lceil b/t\rceil}).                 \tag{2}
```

The provider's parameter is b, rather than an unspecified power of its
circuit description length. Taking `t=ceil(sqrt(b))` makes the additional
cost `exp(O(sqrt(b)))`, preserving every fixed positive logarithmic-saving
exponent of F. The existing exact-count interface is a separate sufficient
provider contract; this one explicitly charges all direct masks and tags.

There is also a precise univariate form of the missing operation. Let P(X)
be the ordinary division-free circuit in (1) before summation, and let R be
its degree-below-`2^d` remainder modulo `X^(2^d)-X`. Finite-field power sums give

```math
\sum_{u\in\mathbb F_{2^d}}P(u)=[X^{2^d-1}]R(X).                    \tag{3}
```

The constant term sums to zero; the only positive degree below `2^d`
whose power sum is nonzero is `2^d-1`. Its sum is one. This also holds
at d=1. Thus a provider may target that one coefficient of the reduced
shared circuit. Explicitly writing R or all field values costs exponential
size in b; (3) does not supply a fast reduction or coefficient algorithm.

## 3. Finite evidence and attribution

```sh
make verify-native-parity-hashes
```

The [checker](../tools/check_native_parity_hashes.py) constructs all native
field points for its fixed models and runs the actual precision-one trace
and cofactor circuits. Independent literal divisibility, chunk Horner
evaluation and bitwise hash arithmetic check the masks and recovery.
The [receipt](../evidence/native-parity-hash-controls.json) records 1172
native point controls, 504 direct hash-mask controls and 3208 scalar parity
queries. The test oracle caches its exhaustively evaluated field states;
it is not an implementation of the missing uniform fast sum.

In the compressed d=10 fixture `N=67*131*257`, t=6 gives D=16 and a
32-element control field. Its point-zero hash puts factors 67 and 131 in
the same even bucket, so that bucket's parity is zero. The full hash family
still recovers minimum factor 67. A deliberately zero mask outside that
acquired family also checks the false-singleton guard: candidates 3,7,11
for N=231 have scalar one but XOR tags decode to 15, rejected by division.

The [author's full default receipt](../evidence/native-parity-full-verification.json)
and [complete compressed log](../evidence/native-parity-full-verification.stdout.txt.xz)
record exit-zero `make verify` in 864.10 seconds at clean frozen
`d76f0eee995405be8a407be2f11df1f7bbffc998`, tree
`7314ca3abcce8ae473e1896a50b75e0584e73a73`. The original metadata records
the same source/tree after the run, clean tracked state, and a fresh
private-/tmp marker and different mount namespace, with read-only host
root, only that checkout writable and network unshared. Sol worker 2
checked and archived this existing author receipt and complete stdout,
adding raw/archive hashes; this is not an independent full replay.
Later integration at `aadc9f1` has separate native-parity, coordinate-word,
completion, spatial-fusion and docs1205 focused checks, outside that pin.

Sol worker 2's [scoped ordinary inspection](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19363530)
at that frozen source found no gap in the inspected direct-mask,
hash isolation, verified-candidate recovery and original-b cost steps.
It did not rederive every imported field/native constructor or fixture,
independently replay the full suite, or supply the missing evaluator.

Polynomial hashing, sparse isolation, bit tags and candidate verification
retain the existing bit-tagged recovery credit and implementation. Trace
duality, repeated-squaring traces, polynomial Boolean gates and the no-wrap
predicate retain the native-word construction's source credit. The change
here is direct masked-parity composition and its implemented consumer,
developed with OpenAI Codex assistance. No new Lean, independent mathematical
review or faster factoring bound is claimed.
