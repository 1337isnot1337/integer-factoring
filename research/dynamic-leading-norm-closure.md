# Dynamic leading norms and a smaller exact presentation

## Status

This checkpoint gives a deterministic adaptive relation constructor. It
reduces a target prime-form power, factors its leading coefficient, and
adds the resulting prime forms when they are new. Every integer submitted
to the factoring routine has at most half the input bits, up to an
absolute additive constant. Fresh prime forms with unit coefficients can
be eliminated immediately from the presentation.

The exact 60-bit control factors

\[
  1000000016000000063=1000000007\cdot1000000009
\]

using 2,208 factorizations of integers of at most 30 bits. Its 276 original
prime-form relations reduce to eight active generators and eight retained
relations. It performs no coordinate or ellipse search.

There is no proved subexponential bound on the number of prime forms
introduced. The finite control is therefore an adaptive algorithm and a
precise complexity interface, not an improved complete factoring exponent.
The present unconditional score remains unchanged.

The implementation is
[dynamic_leading_norm_relations.py](../tools/dynamic_leading_norm_relations.py).
The finite factoring helper
[certified_small_factor.py](../tools/certified_small_factor.py) emits
recursive Lucas certificates. Its Pollard schedule is only an accelerator;
complete trial division is its fallback. No subexponential bound is
claimed for that helper.

## 1. A discriminant that makes leading-norm decoding direct

Let \(N\ge3\) be odd. Set

\[
 \Delta=\begin{cases}-N&N\equiv3\pmod4,\\-4N&N\equiv1\pmod4.
 \end{cases}
 \qquad D=|\Delta|,\qquad H=\left\lfloor\sqrt{D/3}\right\rfloor.
 \tag{1}
\]

Every reduced primitive positive definite form \(Q=(A,B,C)\) of
discriminant \(\Delta\) has

\[
 1\le A\le H,
 \tag{2}
\]

because \(|B|\le A\le C\) gives \(D=4AC-B^2\ge3A^2\).
Thus \(A\) has at most \(\tfrac12\log_2N+O(1)\) bits.

First compute \(\gcd(A,N)\). A proper gcd already factors the input.
Otherwise \(\gcd(A,N)=1\), since \(A<N\). Suppose a complete prime
factorization of \(A\) is available.

For each odd \(p^e\Vert A\), the residue \(B\bmod2p\) selects an
orientation of the prime form above \(p\). Its whole exponent \(e\)
has that orientation. This follows from primitivity and the primitive
ideal of norm \(A\) associated to \(Q\), as in the
[direct norm converter](direct-discriminant-short-norms.md) and
[HNF arithmetic](../tools/quadratic_ideal_hnf.py).
No new modular square-root search is required for an introduced prime:
the middle coefficient already supplies its root.

### Odd discriminant

If \(N\equiv3\pmod4\), then \(\Delta=-N\) is odd. All prime divisors
of \(A\) are coprime to \(\Delta\). If 2 is inert, then
\(\Delta\equiv5\pmod8\), and \(B^2-4AC\equiv5\pmod8\) forces
\(AC\) odd. Thus 2 cannot divide \(A\). If 2 splits, the root of
\(\Delta\) modulo 8 is supplied by odd \(B\), and the ordinary
prime-form decoder applies.

### Even discriminant

If \(N\equiv1\pmod4\), explicitly include

\[
 \tau=\left[2,2,(N+1)/2\right]
 \tag{3}
\]

among the generators, with weight base 2. It has order at most two.
If \(A\) is even, then primitivity and \(B^2-4AC=-4N\) imply

\[
 A\equiv B\equiv2\pmod4,\qquad C\equiv1\pmod2.
 \tag{4}
\]

Indeed, writing \(B=2b\) gives \(AC=N+b^2\). Even \(A\) forces
odd \(b\), hence \(AC\equiv2\pmod4\). Primitivity then makes
\(C\) odd and \(A\equiv2\pmod4\).

The form

\[
 Q'=(A/2,B,2C)
 \tag{5}
\]

is primitive with the same discriminant and odd leading coefficient.
The usual coprime leading-coefficient composition with (3) gives

\[
 [Q]=\tau[Q'].
 \tag{6}
\]

For clarity, choose the middle coefficient in that composition to be
\(B\): it is congruent to 2 modulo 4, and is already the middle
coefficient of \(Q'\) modulo \(A\). The product then has leading
coefficient \(A\) and trailing coefficient \(C\). Decode the odd
factorization of \(A/2\), and append \(\tau\) once.

In both branches, decoding produces an exact signed prime word \(v\)
with

\[
 [Q]=\prod_p g_p^{v_p},\qquad
 \prod_p p^{|v_p|}=A\le H.
 \tag{7}
\]

The implementation checks the class identity by exact HNF composition.
Repeated prime factors of \(N\) cause no exception: the discriminant
can define a nonmaximal order throughout.

## 2. A finite closure construction with independent rows

Start with a finite set of unramified split prime forms. Include (3) in
the even-discriminant branch. Keep a queue of known primes that have not
yet supplied a relation.

When processing prime \(p\), choose an exponent

\[
 k\ge K_p,\qquad K_p=\min\{k\ge1:p^k>H\}.
 \tag{8}
\]

Reduce \(g_p^k\), factor its leading coefficient, and decode (7).
Add any new prime forms to the queue. Retain the relation

\[
 r_p=k e_p-v.
 \tag{9}
\]

Each prime is processed only once. The finite scheduler tests a prescribed
number of consecutive exponents beginning at \(K_p\). It selects the
candidate that introduces the fewest new free generators after the unit
elimination in the next section. The tie breaks are the number of new
primes, their total bit length, the leading norm, and the exponent.
All tested norms are charged and recorded, including rejected candidates.

**Independence.** Multiply column \(q\) of the original relation matrix
by the positive real number \(\log q\). In row \(p\),

\[
 k\log p>\log H\ge\sum_q |v_q|\log q.
 \tag{10}
\]

The square minor on the distinct processed pivot primes is strictly
row diagonally dominant. In particular, every prefix of the original
rows is independent over \(\mathbb Q\). The implementation compares
integer products, so it does not approximate logarithms.

**Finite closure.** Every introduced prime divides a leading norm at most
\(H\). There are at most \(\pi(H)\) such primes. Processing each prime
once therefore closes after at most the initial base size plus
\(\pi(H)\) steps. At closure there is one independent relation for
every original generator.

This is a termination statement. The bound is exponential in the input
bit length and does not improve the factoring score.

## 3. Eliminate fresh unit columns before building a large matrix

Suppose (9) introduces a prime \(q\) whose exponent in \(v\) is
\(\varepsilon\in\{-1,1\}\). No earlier relation contains \(q\).
After substituting already known expressions, (9) has the form

\[
 w-\varepsilon e_q=0.
 \tag{11}
\]

Store \(g_q\) as the word \(\varepsilon w\) in the current active
generators. Eliminate the defining relation and the column for \(q\).
Keep \(q\) in the queue: its own target power can still produce a later
relation. If several fresh primes have unit coefficients, the finite
scheduler eliminates the largest one.

This is an integral change of presentation, with a unit pivot. It changes
neither the presented abelian group nor its map to the form class group.
The other fresh primes remain as active generators. Repeating the process
preserves independence of the retained relations: an integral elimination
of \(u\) unit rows and columns decreases the rank of the original
independent relation matrix by exactly \(u\).

**Smaller active primes.** Every newly retained active prime is at most
\(\sqrt A\le\sqrt H=O(N^{1/4})\). A prime \(q>\sqrt A\) dividing
\(A\) is unique and has exponent one. If it is fresh, the rule that
eliminates the largest fresh unit prime must eliminate it. Thus the active
basis has at most the initial base size plus \(\pi(\sqrt H)\) primes.
Already known primes do not affect this argument. The dictionary of
eliminated primes can still grow throughout the larger range up to
\(H\); this bound on active dimension is not a bound on the number of
norm queries or on the pending queue.

The implementation stores all prime-form expressions in a fixed growing
set of active generators. It does not recursively expand expressions at
evaluation time. If \(m\) original rows have been processed and each step
tests at most \(\ell\) consecutive exponents, a new expression has
coefficient \(\ell_1\) norm at most \(O(b+\ell)\) times the previous
maximum, where \(b\) is the input bit length. Consequently its entries
have

\[
 O\bigl(m\log(b+\ell)\bigr)
 \tag{12}
\]

bits. Substitution, exact group evaluation, and retained relation
processing have polynomial cost in \(b,m,\ell\). Large numerical
exponents do not hide exponentially large written expressions.

After each retained row, call the
[determinant-free partial extractor](determinant-free-ambiguous-extraction.md).
It can find a factor before closure. If it does not, the failure concerns
the image of the presented quotient's torsion subgroup; it does not show
that the full class group lacks useful involutions.

## 4. Why a full generated closure factors a non-prime-power

For \(\Delta=-4N\), the
[direct-discriminant classification](direct-discriminant-short-norms.md)
shows that at most one nonidentity ambiguous class is useless for
factoring: \(\tau\). The determinant-free extractor has exactly that
allowance.

For \(\Delta=-N\), every nonidentity ambiguous class is useful.
A reduced ambiguous form has \(B=A\) or \(A=C\), because odd
discriminant excludes \(B=0\). In the first case
\(N=A(4C-A)\); \(A=1\) is principal, while \(A>1\) gives a proper
factor. In the second case
\(N=(2A-B)(2A+B)\); both factors are proper unless the form is the
principal form at \(N=3\).

A useful ambiguous class exists whenever \(N\) has at least two distinct
prime factors. This can be seen without computing a genus rank. Choose
an unknown coprime factorization \(N=uv\) with \(1<u\le v\).
For \(N\equiv1\pmod4\), the primitive diagonal form \((u,0,v)\)
supplies one. For \(N\equiv3\pmod4\), the factors have opposite
residues modulo 4, and

\[
 (u,u,(u+v)/4)
 \tag{13}
\]

is primitive of discriminant \(-N\) and is equivalent to its inverse.
It is nonprincipal: its value at \((x,y)\) equals
\(u(x+y/2)^2+vy^2/4\), which cannot be 1. For \(y=0\), the value
is either zero or at least \(u>1\); for \(y\ne0\), one has
\(v\ge5\) and hence the value exceeds 1.

Therefore a full independent presentation of a generating set must
yield a factor. The
[polylogarithmic prime generator theorem](polylog-prime-form-generators.md)
supplies a deterministic starting base under the pinned all-number-field
finite-order Hecke zero-free-strip input. The algebra and the finite
replays here do not use that analytic theorem. Small experimental bases
are not claimed to meet its uniform cutoff.

## 5. The precise complexity question

Assume a deterministic scheduler has a bound \(R(b)\) on all norm
factorization calls and the size of its presentation before it finds a
factor, for every preprocessed \(b\)-bit composite with at least two
distinct prime factors. Include all tested exponents and introduced
generators in this bound. Starting-base construction and all other
processing has cost \(C(b)\), polynomial in \(b\) and \(R(b)\).

Every norm to factor has at most \(b/2+O(1)\) bits. Complete recursive
factorization therefore obeys a bound of the form

\[
 F(b)\le R(b)\,\mathrm{poly}(b)\,
 F(\lceil b/2\rceil+O(1))+
 C(b).
 \tag{14}
\]

The polynomial factor includes completely splitting a norm, input
primality and perfect-power preprocessing, and recording all prime factors.
It also includes the at most \(O(b)\) nodes of the outer factor tree.
A proper factor can leave a cofactor of nearly \(b\) bits, so that tree
cannot itself be treated as halving at every split. Bounding all its nodes
by a monotone input-length envelope adds only this polynomial overhead.
With these polynomial overheads included in the effective call count,
this is exactly \(T(b)\le C(b)+R(b)T(b/2+O(1))\). Primality can be
handled by a deterministic polynomial-time test; the finite helper's
recursive certificate-generation cost is not assumed free or polynomial.
For a fixed \(\alpha\ge0\), a uniform bound

\[
 R(b)=\exp\left(O\left(\frac b{(\log b)^\alpha}\right)\right)
 \tag{15}
\]

would preserve the same \(\alpha\) in the complete factoring bound.
Indeed, taking logarithms of the multiplicative recurrence gives a sum
over input sizes \(b,b/2,b/4,\ldots\); its leading terms form a
geometrically decreasing sequence. The final small-size tail contributes
only a polynomial factor. In particular, a polynomial closure bound
would give \(T(b)=2^{O((\log b)^2)}\), a quasipolynomial bound in the
input bit length. Neither such closure bound is established here.

Equation (15) is the missing theorem. The present construction proves
that every newly introduced prime is bounded by \(H\), but not that the
largest prime, the number of pending primes, or another useful potential
decreases. Recursive factoring of the half-size norms alone does not
establish a fast closure bound.

## 6. Exact finite evidence

The checker verifies every supplied prime using recursive Lucas
certificates, every norm factorization by multiplication, every selected
prime word by HNF composition, and the final retained matrix by exact
rational rank. Certificate-only replay runs the same deterministic
scheduler on the recorded complete factorizations; it does not rerun
Pollard factoring.

For the 60-bit input above, the cutoff is 31 and the scheduler tests eight
consecutive exponents per processed prime. It starts with six generators,
processes 276 distinct primes, and performs 268 fresh unit eliminations.
The eight retained relations factor the input. The largest leading norm
has 30 bits.

The evidence also records both discriminant branches, repeated-prime
inputs, direct gcd and target-involution exits, a prime input that closes
without a factor, and a bounded 92-bit run that remains unresolved. An
unresolved bounded run is recorded as such; its partial presentation is
not treated as evidence of a complete factoring bound.

Run:

```sh
python3 tools/check_dynamic_leading_norm_relations.py
```

The default uses the stored certificates in
[dynamic-leading-norm-relations.json](../evidence/dynamic-leading-norm-relations.json).
Use `--regenerate` only to construct fresh finite factoring certificates.
