# Direct exact divisor counts in sublinear workspace

## The unconditional half-rate result

For an `n`-bit integer `N>=2` and a read-only binary prefix,
the exact number of divisors with that prefix can be obtained from
**one weighted evaluator query** with

\[
T(n)=\exp\!\left(O\!\left(\frac n{\sqrt{\log(n+2)}}\right)\right),
\qquad
S(n)=O\!\left(\frac n{\sqrt{\log(n+2)}}\right).                  \tag{1}
\]

The input is read-only and the result is streamed to a write-only
output. The prefix itself is read-only; it need not fit in the
workspace. In particular, the empty-prefix query gives `tau(N)`
directly. An interval with read-only endpoints uses at most `2n+1`
such queries and has the same bounds.

This combines the integral scale-one family from
[integral-counting-transfer.md](integral-counting-transfer.md)
with logarithmic type blocking and the short uniform divisor bound.
It improves the recovery interface for counts: no preliminary
factorization, divisor enumeration, or witness-isolation transform
is used to obtain the count.

The time exponent remains one half. This note does not assert that
the usual materialized `n`-bit prefix self-reduction for the least
factor fits in (1). Least-factor and full-factorization outputs
continue to use the short-descriptor construction in
[sublinear-space-half-rate.md](sublinear-space-half-rate.md).

## 1. A modulus longer than the count, shorter than the input

For `n>=4`, let

\[
L=\lfloor\log_2 n\rfloor,\qquad
q=\min\left(n,\left\lceil\frac{9n}{L}\right\rceil\right).
\]

Section 1 of the sublinear-space theorem proves uniformly that

\[
\tau(N)\le2^q,\qquad q=O(n/\log(n+2)).                           \tag{2}
\]

The proof splits small and large primes and uses the elementary
inequality `L(L+1)<=6*2^(L/2)`. It does not require knowing a
factorization of `N`. The finitely many smaller inputs are handled
directly.

Use the coefficient ring

\[
R=\mathbb Z/2^{q+1}\mathbb Z.                                  \tag{3}
\]

For a prefix with `k` unrestricted low bits and read-only high
part `x`, put

\[
C(N,k,x)=\#\{d\mid N:\lfloor d/2^k\rfloor=x\}.
\]

The multiplication-network contraction over the integers is exactly
this count, by unique extension of each factor pair. The integral
crossing identities preserve the contraction in the ring (3).
Consequently the canonical output residue equals the count itself:

\[
0\le C(N,k,x)\le\tau(N)\le2^q<2^{q+1}.                         \tag{4}
\]

No field property is required. Ring elements have `q+1` bits,
and arithmetic modulo (3) is ordinary integer arithmetic with
truncation. All exact divisions in the integral coefficient generator
finish before reduction.

## 2. Charge the large residues separately from the binary labels

The translated integral Schur family has polynomially many pieces,
uniform `2^O(t)` bit-time factor and metadata queries, separate
binary widths at most `2t`, and total width

\[
2\sqrt2\,t+O(\log(t+2)).
\]

These are the proved interfaces in Sections 1--2 and equation (7)
of the integral-counting note. Apply the logarithmic type-blocking
construction in
[simultaneous-time-space-transfer.md](simultaneous-time-space-transfer.md).
It uses `r=Theta(log(lambda))`-bit original blocks, with a fixed
small constant chosen from their runtime bound. Their raw integer
or rational scratch and outputs then occupy at most
`O(sqrt(lambda+2))` bits. They are reduced to (3) only after
finishing the exact integer calculation.

The new family retains the half-rate contraction, has
`polylog(lambda)` piece-index length, and has `O(lambda)`
binary label/channel workspace. Complete-piece and capped-factor
entries use explicit type-word enumeration, which is valid over
every commutative ring. Only one block oracle and a constant number
of residue accumulators are live in a local query. Thus its workspace
is

\[
O(\lambda+q).                                                \tag{5}
\]

Use the original top constructor and the `n^(3/4)` early cutoff.
The packed grid proof bounds all binary assignments, channel
counters, piece descriptors, and grouping data by
`O(n/sqrt(log(n+2)))` bits. There are `O(loglog(n+2))` active
structural/query frames, with a constant number of residues in
each. The extra scalar storage is therefore

\[
O(q\log\log(n+2))
=O\!\left(\frac{n\log\log(n+2)}{\log(n+2)}\right)
=o\!\left(\frac n{\sqrt{\log(n+2)}}\right).                     \tag{6}
\]

The local scratch in (5) is reused; it is not charged once per
block or once per crossing site. No array of residues is attached
to the binary label coordinates. This proves the space bound in
(1).

Each residue operation has polynomial bit cost in `n`. Hence
charging it changes only polynomial factors in the established
evaluator time. The single canonical residue in (4) is the exact
answer. It occupies only `q+1` working bits before being emitted.

## 3. Read-only intervals need only short prefix descriptors

For a read-only `n`-bit endpoint `X`, decompose `[0,X)` into
prefix blocks. For every position `j` at which the bit of `X`
is one, use the block whose higher bits equal those of `X`,
whose bit at `j` is zero, and whose lower `j` bits are free.
These disjoint blocks cover exactly `[0,X)`.

A query needs only the endpoint identity and the position `j`.
On a prefix-bit request, read the corresponding endpoint bit or
return the single changed zero. The prefix is never materialized.
The descriptor and bit-position counters use `O(log n)` bits.
There are at most `n` blocks, and their count accumulator is
bounded by `tau(N)`, so it uses at most `q+1` bits.

For read-only `0<=A<=B<2^n`, the interval count is

\[
\#\{d\mid N:A\le d\le B\}
=\#\{d\mid N:d<B\}+C(N,0,B)-\#\{d\mid N:d<A\}.
\]

This takes at most `2n+1` direct prefix queries. Each endpoint
remains read-only, and all additional workspace is absorbed in (1).

## 4. A rational-family version and the endpoint distinction

Suppose the rational family satisfies the hypotheses of
[rational-simultaneous-transfer.md](rational-simultaneous-transfer.md),
with `alpha=log_2(1/rho)`. Its small-block denominator product
`D` has `O(n^(1/32))` bits and can be constructed in polynomial
time. Replace its coefficient modulus by

\[
J=2^{q+1}D+1.                                                \tag{7}
\]

Every normalized small-block factor denominator divides `D` and is
coprime to `J`. Type blocking uses only addition and multiplication,
so all its factor entries are defined in `Z/JZ`. The scale-one
crossing identities therefore survive
in `Z/JZ`. Also `J>2^q`, so the output residue is again the
exact count. The coefficient ring has `O(q+n^(1/32))`-bit
elements; forming it and reducing one rational output need only
polynomial extra work.

The same direct-query argument gives time
`exp(O(n/(log(n+2))^alpha))` and workspace

\[
O\!\left(
\frac n{(\log(n+2))^\alpha}
+\frac{n\log\log(n+2)}{\log(n+2)}
\right).                                                     \tag{8}
\]

For every fixed `0<alpha<1`, the second term is smaller, so the
simultaneous saving is retained. At `alpha=1`, the current direct
evaluation stores `O(loglog n)` full scalar registers and proves
`O(n loglog(n+2)/log(n+2))` space. This extra factor is a limit
of the present single-query accounting, not a space lower bound.
Section 5 removes it with polynomially many direct queries.

Applying Section 8 of the rational transfer with `alpha=1/2`
gives a second unconditional derivation of (1). Its denominator
product can be the product of the positive source scales over
the required small block lengths, with polylogarithmic bit length.
One query over the odd modulus (7) then returns the exact count.
Thus the one-query theorem can use either the translated integral
Schur family or the normalized pinned Appendix B source family.
The latter derivation does not depend on the new integral-Schur
construction; it still uses the same packed global evaluator.

## 5. Small-modulus CRT removes the endpoint scalar overhead

The single-query accounting in (8) can be improved at `alpha=1`
by using polynomially many direct queries. The same method also
works at every `0<alpha<1`. It gives

\[
T(n)=\exp\!\left(O\!\left(\frac n{(\log(n+2))^\alpha}\right)\right),
\qquad
S(n)=O\!\left(\frac n{(\log(n+2))^\alpha}\right)                 \tag{9}
\]

for exact read-only prefix and interval counts under the rational
family hypotheses. It does not use witness isolation. All coefficient
moduli in an individual evaluator call have only `O(log n)` bits.

Keep the small-block denominator product `D` from Section 4;
for an integral family, take `D=1`. Set

\[
b=\operatorname{bitlength}(D),\qquad a=b+1,\qquad
H=4(q+a^2+4).
\]

Since `b=O(n^(1/32))`, we have `H=O(n/log(n+2))` for large
inputs, with family-dependent constants. Enumerate the integers
`3<=j<=H` satisfying `gcd(j,2D)=1`. Evaluate the same count
modulo each such `j`. Every factor denominator is a unit in
`Z/jZ`, so each returned residue is the reduction of the same
integer count. Composite and overlapping moduli are permitted.

### An elementary supply of coprime moduli

Let `L_H=lcm(1,...,H)`. For every `0<=k<H`,

\[
H{H-1\choose k}\ \mid\ L_H.                                  \tag{10}
\]

Indeed its valuation at a prime `p` is the sum over `r=p^e`
of

\[
\left\lfloor\frac H r\right\rfloor
-\left\lfloor\frac k r\right\rfloor
-\left\lfloor\frac{H-1-k}r\right\rfloor.
\]

Each term is zero or one: writing the two latter numerators as
quotients and remainders, the remaining numerator is the sum of
their two remainders plus one, between one and `2r-1`.
There are at most `floor(log_p H)` nonzero terms, exactly the
maximum prime-power exponent available in `L_H`. This proves
(10). Some binomial coefficient in that row is at least
`2^(H-1)/H`, so

\[
L_H\ge2^{H-1}.                                                \tag{11}
\]

Let `G_H` be the lcm of the eligible integers `j<=H`.
It is obtained from `L_H` by deleting the maximal prime powers
whose primes divide `2D`. There are at most `a` such primes,
and each deleted prime power is at most `H`. Consequently

\[
G_H\ge\frac{2^{H-1}}{H^a}.                                    \tag{12}
\]

By construction `H>=16` and `a<=sqrt H/2`. The elementary
inequality `log_2 H<=sqrt H` for `H>=16` gives
`a log_2 H<=H/2`. Equation (12) therefore implies

\[
\log_2G_H\ge H/2-1>q.                                        \tag{13}
\]

No prime enumeration is required. The lcm and prime-power reasoning
in this proof establish the supply; the algorithm tests eligibility
using gcd alone.

### Generalized CRT with a short retained state

Retain a canonical residue `z` and its modulus `M`, initially
`z=0,M=1`. After obtaining the count residue `c` modulo an
eligible `j`, put `g=gcd(M,j)` and `h=j/g`.
The residues are compatible because they come from the same count.
If `h=1`, nothing changes. Otherwise solve

\[
t\equiv
\frac{c-z}{g}\,(M/g)^{-1}\pmod h,\qquad 0\le t<h,
\]
\[
z\leftarrow z+Mt,\qquad M\leftarrow Mh.                       \tag{14}
\]

The inverse exists because `gcd(M/g,h)=1`. For computation,
first reduce `c-z` modulo `j`, then divide this compatible
remainder by `g`. The update in (14) is already canonical:
`0<=z<M` before the update and `0<=t<h` imply
`0<=z+Mt<Mh`.

Stop as soon as `M>2^q`. Equation (13) guarantees that this
happens before the eligible list ends. The update multiplies the
previous modulus by at most `H`, so the stopping modulus and
residue have at most

\[
q+O(\log H)
\]

bits. At stopping, the residue is the exact count by (2).
All CRT operations have polynomial bit cost and use
`O(q+log H)` space.

Each evaluator call now has only `O(log H)=O(log n)`-bit
scalar registers. Across all active frames these cost
`O(log n loglog n)` bits. The only long CRT residue and modulus
are held once, outside the evaluator, in `O(q)` bits.
The packed binary-label workspace remains
`O(n/(log(n+2))^alpha)`. Since `q=O(n/log n)` and
`alpha<=1`, all these terms fit (9).
There are at most `H` evaluator calls and polynomial extra
arithmetic, so the time bound also remains (9).
The interval construction of Section 3 adds another polynomial
number of calls and only a constant number of `O(q)`-bit
count accumulators.

### Fixed number fields

The fixed-number-field extension in Section 7 of the rational
transfer has the same consequence. Use a fixed integral primitive
element with monic minimal polynomial `f` of degree `e`, and
include all coordinate denominators in `D`. For each eligible
`j`, evaluate in `(Z/jZ)[X]/(f)`. The exact integer count
reduces to its constant coordinate, with all other coordinates
zero. Feed that constant residue into (14). Monicity gives the
unique degree-`<e` coefficient representation, and fixed `e`
makes scalar storage `O(log n)`. The bounds (9) follow unchanged.
This does not assert a result for growing algebraic degree.

## Verification

The checker tests actual integral Schur crossing entries over
power-of-two rings, exact prefix residues, the read-only interval
cover, and rational reductions under (7). Long power-of-two inputs
with known divisors exercise cases where `q<n`. It does not
implement or benchmark the full asymptotic grid evaluator; that
algorithm and its cost bounds are supplied by the linked proofs.
The endpoint checks also verify the binomial/lcm divisibilities,
the exact supply inequality, exhaustive CRT reconstructions with
overlapping composite moduli, and constant-coordinate extraction
in a fixed quadratic quotient ring.

```sh
python tools/check_direct_sublinear_counts.py \
  --output evidence/direct-sublinear-divisor-counts.json
```
