# Smooth leading norms on consecutive prime-form powers

## Status and finite result

A fixed-base consecutive-power search factors the 92-bit integer

\[
 4951760154835678088235319297
 =2147483647\cdot2305843009213693951.
\]

With prime cutoff 1021, it uses 85 prime-form generators and finds a
factor after 82 independent relations. It tests 218,738 reduced leading
norms using one modular smoothness test per candidate, followed by only
5,433 divisions by known base primes. Every tested norm
has at most 47 bits. The largest tested exponent is 9892.

The 124-bit control also factors:

\[
 21267647932558653899591465697288388633
 =2305843009213693951\cdot9223372036854775783.
\]

It uses cutoff 4093, 283 generators, and 274 independent relations.
The search tests 3,650,183 leading norms of at most 63 bits, followed by
56,571 known-prime divisions. One target exhausts its 100,000-step
window; the remaining partial presentation still supplies a factor.

This search uses neither a factoring oracle nor a supplied norm list.
It constructs the target powers itself, reduces them, and tests their
leading coefficients for smoothness. It introduces no new primes and
performs no coordinate or ellipse search.

The missing asymptotic statement is a uniform bound on the required
power window. The finite result does not establish a stronger complete
factoring exponent. The algorithm and the recorded search are exact.

Files:

- [leading_power_smooth_relations.py](../tools/leading_power_smooth_relations.py)
- [check_leading_power_smooth_relations.py](../tools/check_leading_power_smooth_relations.py)
- [leading-power-smooth-relations.json](../evidence/leading-power-smooth-relations.json)

## 1. A directly recognizable candidate family

Use the discriminant and height from
[dynamic leading-norm closure](dynamic-leading-norm-closure.md):

\[
 \Delta=\begin{cases}-N&N\equiv3\pmod4,\\-4N&N\equiv1\pmod4,
 \end{cases}
 \qquad H=\left\lfloor\sqrt{|\Delta|/3}\right\rfloor.
 \tag{1}
\]

Build a fixed set of unramified split prime forms \(g_p\) for primes
\(p\le Y\). If \(N\equiv1\pmod4\), include the norm-two ambiguous
class \(\tau\) as well. A ramified odd base prime supplies a proper
gcd and is handled immediately. Let \(m\) be the number of generators.

For each generator, set

\[
 K_p=\min\{k\ge1:p^k>H\}.
 \tag{2}
\]

Test the reduced forms

\[
 Q_{p,k}=\mathrm{red}(g_p^k),
 \qquad k=K_p,K_p+1,\ldots,K_p+W-1.
 \tag{3}
\]

After one binary exponentiation to construct the first form, successive
forms need one composition with \(g_p\). Their leading coefficients
\(A_{p,k}\) all lie between 1 and \(H\).

When the current leading coefficient is coprime to \(p\), the code uses
the ordinary coprime CRT composition directly and then reduces. At a
shared prime, it falls back to the general HNF composition. This changes
the implementation cost, not the target class or the candidate sequence.

Compute \(\gcd(A_{p,k},N)\). Let \(P\) be the product of the known
base primes, and put \(B=\mathrm{bits}(H)\). The exact test

\[
 P^B\equiv0\pmod {A_{p,k}}
 \tag{4}
\]

holds precisely when the leading norm is supported on the base. Indeed,
every prime exponent of a number at most \(H\) is at most \(B\), so
the test is equivalent to \(A_{p,k}\mid P^B\). The modular power is
computed without writing \(P^B\). Only candidates passing this test
are divided over the prime list. A unit residual gives their complete
factorization directly;
there is no call to a general integer-factoring routine. The middle
coefficient supplies every prime orientation. The ramified-two case is
decoded by the exact identity \([A,B,C]=\tau[A/2,B,2C]\) proved in the
dynamic note.

Accept the first smooth leading norm for each target prime. Decode its
signed word \(v\), and retain

\[
 r_p=k e_p-v.
 \tag{5}
\]

The coordinate witness is always \((1,0)\) in the reduced target form.
Thus this is a constructive search for both a represented norm and its
coordinates, with a variable target exponent.

## 2. Independence for arbitrarily wide windows

The exact decoder gives

\[
 \prod_q q^{|v_q|}=A_{p,k}\le H<p^k.
 \tag{6}
\]

Consequently the minor on the distinct accepted target primes is
strictly diagonally dominant after weighting column \(q\) by
\(\log q\). Every accepted prefix is independent over \(\mathbb Q\).
The exponent can be much larger than \(K_p\); independence still holds.

The implementation does not reuse the small-exponent product guard of
the height-only relation wrapper. It establishes (6) by checking
\(k\ge K_p\), verifying the factored norm, and using the reduced-form
height bound. Computing \(K_p\) requires only multiplying until the
product first exceeds \(H\). No integer of the size \(p^k\) is
materialized for a large selected exponent.

Every relation coefficient has bit length
\(O(\log(b+W))\), apart from the already polynomial size of the input
base. Here \(b\) denotes the bit length of \(N\). The
[generic determinant-free extractor](determinant-free-ambiguous-extraction.md)
has polynomial cost in the number of rows, \(b\), and this coefficient
bit length. It can use a partial presentation.

The implementation uses an incremental binary row reduction to decide
when to try extraction early. It also tries extraction on the final
accepted matrix. This final call is needed because identity half-relations
can reveal new dependencies that are absent from the original rows
modulo two.

## 3. Exact cost and the remaining theorem

There are at most \(mW\) candidate forms. Each needs one class-group
composition, a gcd, and a modular power whose base has at most
\(m\log_2Y\) bits and whose exponent is \(O(b)\). Trial division
is needed only for the accepted candidates. The complete bounded
attempt therefore costs

\[
 O\bigl(m^2W\,\mathrm{poly}(b,\log Y)\bigr)
 +\mathrm{poly}(m,b,\log W),
 \tag{7}
\]

in addition to constructing the base. All rejected candidates are
included. In particular, recorded prime-factor certificates do not stand
in for an uncharged factoring operation.

More explicitly, put \(L=\mathrm{bits}(P)\le m\log_2Y+1\).
With classical integer arithmetic, reducing \(P\) modulo one leading
norm costs \(O(Lb)\) bit operations. The subsequent modular exponent
has \(O(\log b)\) multiplications of \(O(b)\)-bit residues, for
\(O(b^2\log b)\) additional bit operations. Thus the complete modular
filter costs
\(O(mW(Lb+b^2\log b))\). Accepted candidates require at most
\(m(m+O(b))\) divisions in total: one pass over the base per accepted
norm, plus at most \(O(b)\) successful prime divisions per norm.

Each fast CRT update uses one modular inverse, a constant number of
integer products and divisions, and reduction of a form whose entries
have \(O(b+\log Y)\) bits. The shared-prime fallback is the already
proved polynomial-time constant-dimensional HNF composition. The first
power costs \(O(\log K_p)\) compositions per generator, and at most
\(mW\) further updates are made. These bounds account for the full
search, rather than only for verification of its accepted witnesses.
Memory is polynomial in \(m,b,\log Y,\log W\): rejected forms are discarded,
and the stored relation matrix has at most \(m^2\) entries.

If the starting base generates the class group and every target has a
smooth leading norm in its window, the accepted matrix is square and
nonsingular. For an odd composite with at least two distinct prime
factors, the ambiguous-class classification and extractor then give a
proper factor. Standard primality, perfect-power, and recursive splitting
steps complete factorization with polynomial additional overhead.

The
[polylogarithmic generator theorem](polylog-prime-form-generators.md)
provides a base of size polynomial in \(b\), under the pinned
all-number-field finite-order Hecke zero-free strip. Thus the following
would suffice for a score \(\alpha\):

> For that generating base, a deterministic choice of windows of length
> \(W(b)=\exp(O(b/(\log b)^\alpha))\) always supplies the required
> independent relations, or produces a useful partial presentation sooner.

No such uniform short-window theorem is proved here. A polynomial
window bound would give a polynomial-time algorithm through (7), which
shows how strong that possible statement would be.

Consecutive exponents are an implementation choice. The same proof
accepts any deterministic list of candidate exponents \(k\ge K_p\)
whose written bit lengths are polynomial in \(b\). Binary exponentiation
then costs a polynomial number of group operations per candidate.
Therefore a future selector can instead supply a sparse list of large
exponents; its generation and all rejected tests must still be charged.
A sufficient bound on that total list size would give the same factoring
score without requiring a short contiguous window.

The generating base may also be enlarged. For any fixed \(c\ge1\), a
cutoff \(Y=2^{O((\log b)^c)}\), chosen above the generator cutoff,
has setup and fixed-degree polynomial
processing costs negligible on the logarithmic scale
\(b/(\log b)^\alpha\), for every fixed \(\alpha\). In particular,
the [quasipolynomial specialization](quasipolynomial-base-norm-processing.md)
uses \(Y=2^{32\lceil\log_2 b\rceil^2}\). Its existing smooth-set
count gives a global fraction at least

\[
 2^{-(1/64+o(1))b/\log_2 b}
\]

of classes with an odd, base-smooth reduced leading coefficient below
the injective height cutoff. This statement uses the pinned analytic
input and holds for sufficiently large input lengths. The same fraction
is a lower bound for the larger smooth-leading-norm set recognized here.
It is a global class-group density, not a density or hitting-time bound
along each power orbit. A short-list theorem at this enlarged cutoff
would therefore be another sufficient route to score one.

There is a much weaker automatic termination bound. If \(o_p\) is the
order of \(g_p\), every interval of \(o_p\) consecutive exponents
contains a multiple of \(o_p\). Its reduced form is principal and has
leading coefficient 1. Hence a window of length the full class number
always succeeds for every generator. This is exponential in \(b\),
and does not improve the existing factoring bound.

The proved prime-word mixing theorem does not supply the missing
short-window estimate: it averages products of independently chosen
prime forms, while (3) remains in one cyclic orbit. Finite smoothness
frequencies are not used as a worst-case theorem.
The [powered smooth-norm existence theorem](powered-prime-smooth-norms.md)
supplies a represented smooth norm somewhere in each class; it does not
assert that the reduced leading coefficient is smooth. Reduction need
not preserve the smoothness of a represented norm, so that theorem does
not certify the windows used here.

## 4. Replay and interpretation

The checker reruns the complete bounded search, including rejected norm
candidates, its modular tests, and its division count. The plain-division
reference accepts the same candidates; for the 92-bit case it uses
18,850,094 divisions. It separately checks each selected
target power, the complete norm factorization, prime-form orientations,
the exact group relation, and the weighted independence hypotheses.
It verifies the returned factor by division.

The 92-bit control uses window 10,000; all accepted exponents before the
factor are at most 9892. Its final partial presentation has deficiency
three. The extractor evaluates two involutions and returns a proper
factor without an identity saturation or a determinant computation.
The small cutoff is a finite experimental choice. Its success does not
assert that cutoff 1021 satisfies the uniform analytic generation bound.

Run the complete replay with:

```sh
python3 tools/check_leading_power_smooth_relations.py
```

The same local arithmetic supports both discriminant branches and
repeated-prime inputs. A bounded attempt that lacks enough relations is
reported as an exhausted window. A complete presentation of a small
nongenerating subgroup can also return no factor; neither outcome is a
primality certificate.
