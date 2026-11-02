# Sparse divisor data from one faithful character

## Question and result

Base `252128d`, following the [full character decoder](cyclotomic-character-decoder.md).
That decoder requires every divisor-index character order for arbitrary
binary paired vectors. Actual divisor vectors are sparse. Does the
information restriction change under a certified sparsity bound?

Yes. Write q-1=2eu, with e a power of two and u odd. In the full
modulo-two power-basis algebra R=F_2[T]/(Phi_u(T)^e), a single faithful
observation P(T) determines every degree-below-eu binary polynomial P
of weight at most s whenever

```
u=1, or 2s < the least prime factor of u.                 (1)
```

There is a deterministic decoder with polynomial bit cost in q,s.
The kernel's minimum nonzero weight is exactly the least prime factor
of u when u>1, independently of e. The full nonreduced coordinates are
retained. This does not contradict the earlier necessity theorem, whose
quantifier is arbitrary binary vectors.

Classical lower-sieve prime acquisition then gives a **complete conditional
factoring reduction using only faithful-character coefficient vectors**,
preserving every fixed positive logarithmic-saving exponent. An additional
CRT/Linnik construction supplies a more expensive acquisition bound.
All level and arithmetic costs are charged below. A faster algorithm
for those coefficients at the original composite index is still missing.
The complete proven factoring score remains **1/2**.

This changes the coefficient engine's required input family. Reducing a
polynomial number of queries alone does not improve an exponent. The
continuation target is a uniform algorithm for this restricted faithful
coefficient family, rather than additional decoder or prime searches.

## 1. Exact kernel distance at odd order

Let u>1 be odd, ell its least prime factor, and phi=Phi_u modulo two.
Consider the binary cyclic code

```
C_u={P: degree(P)<u, phi divides P}.
```

Its minimum nonzero Hamming weight is ell. Here is a direct proof.

Choose any primitive u-th root omega in characteristic two. If P is in
C_u, then P(omega^r)=0 for every r coprime to u. Suppose P has w<ell
nonzero positions j_1,...,j_w. Every r=1,...,w is coprime to u, so

```
sum_(i=1)^w (omega^j_i)^r=0,       r=1,...,w.             (2)
```

The nodes z_i=omega^j_i are distinct and nonzero. The matrix in (2)
has determinant (product_i z_i) times a Vandermonde determinant, which
is nonzero. It cannot annihilate the nonzero all-one coefficient vector.
Thus w>=ell.

The bound is attained by

```
1+X^(u/ell)+...+X^((ell-1)u/ell).                         (3)
```

It vanishes at every primitive u-th root, because omega^(u/ell) is a
nontrivial ell-th root. Since u is odd, X^u-1 is squarefree; hence phi
divides (3). Its degree is below u and its weight is ell.

Only ordinary Vandermonde and cyclotomic facts are used. This is a
classical cyclic-code argument applied to the preceding character map,
not a new coding lower-bound claim.

## 2. The nilpotent code is an interleaving

For e=2^v, characteristic-two Frobenius gives

```
Phi_u(X)^e=Phi_u(X^e).
```

Write, uniquely,

```
P(X)=sum_(r=0)^(e-1) X^r P_r(X^e),       degree(P_r)<u.
```

Division by Phi_u(X^e) preserves exponent residues modulo e. Therefore
P is in the full kernel exactly when each P_r is in C_u. Its weight
is the sum of the block weights. For u>1 the minimum nonzero weight
remains ell, attained by substituting X^e in (3). For u=1 the modulus
has degree e, so the degree-below-e target has zero kernel.

Two vectors of weight at most s differ in weight at most 2s. The kernel
distance proves (1). It is a sufficient uniform certificate; failure
of (1) is not a claim of ambiguity for every particular input. A sharp
small collision is u=3,e=1: X^2 and 1+X have the same observation.

This interleaving also extracts the observations needed for decoding.
Split the canonical remainder M=P modulo Phi_u(X)^e into
M=sum_r X^r M_r(X^e). Then M_r=P_r modulo Phi_u. Reducing R to its
residue field *before* this split can destroy information. At e=4,u=1,
the distinct weight-one vectors 1 and X coincide modulo X+1, while the
full faithful observation distinguishes them.

## 3. A polynomial decoder with its field acquisition paid

Choose one irreducible factor f of phi and work in F=F_2[Y]/f, with
omega=Y modulo f. It has order exactly u. The field is acquired in
polynomial time in u by the classical Berlekamp construction, directly
implemented here: find the F_2-kernel of a -> a^2-a in F_2[Y]/phi.
Its elements take values zero or one on every irreducible factor.
Successive gcds with a kernel basis leave a single factor, since that
basis separates any two remaining field components. No finite-field
or class-group model is an uncharged oracle.

For each interleaved block, obtain the consecutive moments

```
m_j=M_r(omega^j),             1<=j<=2s.
```

Every such j is a unit modulo u by (1). Consequently phi(omega^j)=0,
so these are exactly the moments of P_r's sparse support. The **full**
remainder M_r is used in these substitutions, even though their answers
lie in one field component. Keeping only M_r modulo f at the start
would not justify evaluation at all omega^j.

This distinction persists when the full coefficient algebra is already
reduced. At q=47, u=23, one acquired field factor is
f=1+X+X^5+X^6+X^7+X^9+X^11. The sparse words
1+X+X^5 and X^6+X^7+X^9+X^11 have the same residue modulo f, but
different full faithful observations modulo Phi_23. Their weights are
three and four, within the certificate 2*4<23. This is a counterexample
for arbitrary sparse vectors to discarding input components; it is not
a lower bound for an alternative arithmetic-specific coefficient interface.
In fact, actual divisor reciprocity recovers the other q47 field from
one field and the known residue of N; see Section 3.1 below. The arbitrary
sparse collision must not be interpreted as an arithmetic obstruction.

If the block has t nonzero positions, its leading t-by-t Hankel matrix
H_ij=m_(i+j+1), 0<=i,j<t, is nonsingular: it factors as
V diag(z_1,...,z_t) V^T. Every larger such matrix is singular. Choose
the largest nonsingular size at most s and solve

```
sum_(j=0)^(t-1) c_j m_(i+j+1)=m_(i+t+1),  0<=i<t.
```

The monic polynomial X^t+sum_j c_j X^j is the support-node polynomial.
Evaluate it at omega^k for 0<=k<u to recover the binary support. The
zero observation decodes to zero by the distance certificate. Interleave
the recovered blocks, verify the sparsity bound and re-encode the
observation. Inconsistent supplied data is rejected.

The field degree is at most u. Kernel acquisition, binary polynomial
arithmetic, at most s dense field systems per block, and the u-node
scan all have polynomial cost and workspace in eu,s. The implementation
uses ordinary Gaussian elimination over the acquired field. The bound
does not enumerate the coefficient index N, divisors of N or an
exponential set of candidate support vectors.

For the stronger distinct-semiprime promise, [one primitive field always suffices](semiprime-single-field.md) at every coprime q>=7, q=3 mod 4, using a quadratic and the known product residue. This separate arithmetic decoder does not assert one-field recovery for arbitrary sparse words.

### 3.1 Actual divisor reciprocity and an inversion-fixed target

The arithmetic inputs have more structure than arbitrary sparse words.
Let q=3 modulo 4, q>=7, u=(q-1)/2 and q not divide N. In the full cyclic
space F_2[X]/(X^u-1), let

```
P_N(X)=sum_(d|N) X^(log_g(d) mod u),
c=log_g(N) mod u,       a=(u+1)/2.
```

All exponents and coefficients are reduced in the stated cyclic space.
The involution d -> N/d proves

```
P_N(X^-1)=X^-c P_N(X).
Q_N(X)=X^(-ac)P_N(X),       Q_N(X^-1)=Q_N(X).            (3a)
```

Here 2a=1 modulo u. Moreover, the coefficient of X^0 in the **cyclic**
word Q_N is the known square indicator of N: all noncentral divisor pairs
cancel, and the only fixed divisor is sqrt(N) when that is integral.
This is not a claim about the constant coordinate of its canonical
remainder modulo Phi_u, which can differ.

Projection of Q_N lies in the inversion-fixed subalgebra of
R=F_2[T]/Phi_u(T). It is a reduced algebra of dimension phi(u)/2, acquired
in polynomial time in q by binary linear algebra. Equivalently, it is
generated by T+T^-1: over an algebraic closure the primitive roots pair as
omega,omega^-1, with distinct sums for distinct pairs. The corresponding
real cyclotomic polynomial is separable. This root order is **u**, unlike
the phase-q nodal algebra in the [Tate/theta compiler](theta-coordinate-programs.md).

This changes the arithmetic target, not the exponent of any supplied
algorithm. Multiplication by the known unit T^(ac) recovers the original
faithful observation with polynomial overhead. Primitive-generator search,
the known-residue discrete logarithm by a q-element scan, fixed-space
acquisition and basis conversion all cost poly(b,q) bit operations. To
keep the complete reduction within this promise, acquire L+b+1 prime
levels and discard q3 and those dividing N. There are at most b dividing primes;
primality/eligibility acquisition and all gcds remain polynomial in b,D.

For an irreducible factor f of Phi_u, let f* be its monic reciprocal.
Equation (3a) determines the f* component from the f component and c.
Indeed the map T -> T^-1 interchanges those factors, and the missing raw
component is T^c P_f(T^-1) modulo f*. Thus one supplied field determines
its whole **inversion orbit**, not automatically every field component.
At q47, Phi_23 has exactly two factors, binary masks 2787 and 3189,
exchanged by inversion. Consequently **one field component suffices for
the full faithful observation on actual divisor inputs at this level**,
without a sparsity assumption. The earlier masks 35 and 2752 satisfy no
twisted-reciprocal promise at any center c.

The converse shortcut remains invalid in general even with reciprocity.
At q503, u251 is prime and Phi_251 has five degree-50 fields. One acquired
factor f has binary mask 1986479280909191 and is self-reciprocal. The word
X^-25 f(X)(X+X^-1) modulo X^251-1 is nonzero, inversion-fixed, has zero
central coefficient, and vanishes modulo f. Split its fourteen opposite
position pairs into two groups of seven. This gives two reciprocal words
of weight fourteen with the same f component and different full faithful
observations, despite 2*14<251. This is a counterexample for the
**reciprocal sparse promise**; neither word is asserted to arise from the
divisors of an actual N. Further arithmetic constraints remain available.
The rough-prime sieve supplies no theorem that its eligible levels have
only one inversion orbit. No stronger acquisition condition is assumed.

There is also an ordinary fixed-weight modular realization of the
normalized target. Let chi be the original complex odd character of
order 2u, and set alpha=chi^-a and beta=chi^(1-a). Both are nontrivial
primitive characters of conductor q. The series

```
E_(1,alpha,beta)(z)
 =sum_(n>=1) (sum_(d|n) alpha(n/d) beta(d)) z^n          (3b)
```

has weight one, level q^2 and quadratic character
alpha*beta=chi^-u=chi^u. Its constant term is zero and all coefficients
are cyclotomic integers. These statements use the coefficient formula,
constant term and parity theorem in
[Stein, Chapter 5, formula (4) and Theorem 5.8](https://wstein.org/books/modform/modform/eisenstein.html).
For q not divide N its Nth coefficient is exactly
chi(N)^-a sum_(d|N)chi(d), which reduces modulo two to Q_N(T).
The complex coefficient need not itself be real: conjugation multiplies
it by chi^u(N), a sign that disappears modulo two. The q3 degenerate case
is excluded from (3b); the existing interface already handles that level.
Raising the level from q to q^2 and converting the known normalization cost
only polynomial overhead in q. A fast original-index coefficient engine
for (3b) is still missing; neither reciprocal symmetry nor the smaller
coefficient algebra supplies one.

The focused checker independently reconstructs all 390 coprime targets
2<=N<400 at q47 from one supplied field and the known c, using literal
divisors for the supplied data and integer cyclotomic observations.
It checks the square indicator and the q503 reciprocal collision. This
bounded panel verifies the implementations of the stated identities,
not a uniform coefficient engine or an arithmetic impossibility result.

## 4. Eligible prime acquisition and its full cost

### 4.1. A direct CRT/Linnik construction

Fix positive integers s and L. Let

```
B=product_(odd primes p<=2s) p,        W=4B.
```

Choose a prime R>max(2s,L), with R=O(s+L), by Bertrand's theorem and
trial division. For j=1,...,L, CRT supplies a residue a_j modulo WR
satisfying

```
a_j=3 modulo 4,
a_j=2 modulo every odd prime p<=2s,
a_j=j modulo R.                                         (4)
```

Each residue is coprime to WR. Search its progression for the least
prime q_j. The levels are distinct by their different nonzero residues
modulo R. They satisfy q_j=3 modulo 4 and

```
q_j-1=2u_j,        u_j=1 or least_prime_factor(u_j)>2s.
```

Thus every acquired level meets (1), with e=1. No unproved statement
about safe primes or shifted-prime density is needed.

The only analytic input is the effective classical least-prime bound
P(a,m)<=C*m^A for coprime a,m. The primary source
[Xylouris, Acta Arithmetica 150 (2011), pp.65--66](https://www.impan.pl/shop/publication/transaction/download/product/83078)
states Linnik's effective theorem in the introduction and the stronger
exponent 5.18 in Theorem 1.1. Only an absolute effective A is used here;
the algorithm stops on finding a prime and does not require a numerical
value of C as advice.

The elementary bound log B=O(s) follows by bounding the primes in each
dyadic interval with the central binomial coefficient. Hence

```
q_j <= exp(O(s)) (L+s)^O(1).                             (5)
```

Computing the CRT data, searching all L progressions by trial primality
tests, and polynomial-in-q acquisition/decoding therefore cost
exp(O(s))*poly(L,s). All constants are absolute. The large numerical
levels are charged, including the full coefficient-vector output length.
This is an asymptotic construction theorem; the finite checker does not
numerically certify Linnik's constants.

### 4.2. Polynomial-size levels by a classical lower sieve

There are L eligible levels bounded by poly(s+L), and they can be found
deterministically in polynomial time in s+L. Here is the analytic input
and its exact specialization. This is an ordinary application of classical
sieve theorems, not a new distribution theorem for shifted primes.

Use [Kevin Ford's Sieve Methods notes (2023), Theorems 3.3 and 3.4,
pp.32,35](https://ford126.web.illinois.edu/sieve2023.pdf). The lower sieve
has a fixed positive main coefficient at dimension one and remainder
moduli through z^3.35. Bombieri--Vinogradov bounds the sum of progression
errors through x^(1/2)/log^B(x), for any requested logarithmic error saving.

Apply these to the integer sequence

```
A={(q-1)/2 : q<=x prime, q=3 modulo 4}.
```

Every member is odd, so g(2)=0. For odd squarefree d, divisibility by d
selects exactly the residue 2d+1 modulo 4d. Thus, with X=li(x)/2,

```
|A_d|=pi(x;4d,2d+1)=X/phi(d)+r_d,       g(d)=1/phi(d).
```

The residue is coprime to 4d and phi(4d)=2phi(d). The local sieve density
is dimension one, and V(z) is bounded above and below by constant
multiples of 1/log z, by Mertens' product estimates and convergence of
the odd-prime product of 1-1/(p-1)^2.

Choose z=x^(1/10). The required moduli satisfy d<=x^(67/200).
The fixed factor four still leaves them below the Bombieri--Vinogradov
level, since 67/200<1/2. Its remainder is O(x/log^10(x)), which is
smaller than the positive main term. Consequently there are at least
c*x/log^2(x) such primes with (q-1)/2 free of primes at most x^(1/10),
for all sufficiently large x and an absolute c>0.

Take x a sufficiently large absolute multiple of (s+L+2)^20. Then
x^(1/10)>=2s and the count exceeds L. This proves the polynomial bound.
The algorithm simply scans q=3,7,11,..., tests primality by trial
division, and rejects q for which an odd prime at most 2s divides q-1.
It stops after L successes. Every operation and its stopping point have
polynomial cost in s+L. Analytic constants and the sufficiently-large
threshold are not input advice: the search does not use them, and an
absolute constant absorbs the remaining finite parameter pairs.

The special progression 3 modulo 4 is fixed. No distribution theorem
with an uncharged growing modulus, safe-prime conjecture, or uniform
claim at a numerical small threshold is used. The published lower-sieve
and Bombieri--Vinogradov theorems supply the asymptotic bound; the finite
controls check only eligibility and the CRT/local-density identities.

## 5. Complete faithful-only conditional factoring invoice

After trial division through 2^t finds no factor, the existing
[trial-budget bound](sparse-divisor-residues.md#3-trial-division-makes-the-ring-degree-small)
gives D=2^ceil(b/t)>=tau(N). At a prime level q, the paired binary
polynomial of the [weight-one reduction](weight-one-divisor-coefficients.md)
has degree below (q-1)/2 and at most D nonzero positions. A single
faithful coefficient, the power-basis vector of

```
sum_(d|N) chi_1(d) modulo two,
```

is its observation P(T). This is one **whole coefficient vector**, not
one scalar bit. The character sends a supplied primitive generator
modulo q to the primitive (q-1)-st cyclotomic root. Divisors divisible
by q contribute zero, as in the predecessor.

Acquire L distinct levels by Section 4.2, with s=D and

```
S=2D,
L=(b+1)S(S-1)/2+S(b+2)+3.
```

Every queried level has q=3 modulo 4, hence e=1. The coefficient engine
therefore needs only the full **reduced** algebra F_2[T]/Phi_u(T).
Its primitive field components are all retained; the one-component
counterexample in Section 3 still applies to the generic sparse input.

Decode each faithful observation and copy its paired coefficients to
the opposite residue positions. This supplies the same signed Laurent
residues as the previous reduction. Its collision and CRT proofs use
only that the L prime moduli are distinct, not that they are the first
L primes. At most (b+1)S(S-1)/2 primes cause any collision between the
signed integer exponents. Maximum support identifies the separated
levels, which include enough primes for exact CRT. Squaring opposite
residue pairs reconstructs product_(d|N)(Y-d^2), with coefficient height
at most 2^(D(2b+1)). Remove Y-1, isolate the least remaining integer root
and take its integer square root. This is the smallest prime factor,
including for primes, prime powers and nonsquarefree inputs.

Suppose the faithful coefficient vector at any of these supplied
levels is obtained uniformly in F(b)*poly(b,q) bit operations, with
fixed polynomial degrees. The decoder, polynomial level acquisition and
earlier complete recovery give total bit cost

```
2^t poly(b) + F(b) poly(b,D),         D=2^ceil(b/t).       (6)
```

This includes parameter acquisition, coefficient-ring acquisition, output,
sparse decoding, signed support recovery, CRT and integer root isolation.
The algorithm's described searches and finite exceptions are deterministic.

One fixed computable schedule t=ceil(sqrt(b)) gives
D<=2^(sqrt(b)+1). For F(b)>=1, equation (6) is therefore

```
F(b) exp(O(sqrt(b))).                                    (7)
```

For every fixed beta>0, sqrt(b)=o(b/log^beta b). Thus an engine with
F(b)=exp(O(b/log^beta b)) retains that logarithmic-saving exponent under
the same schedule, including beta>1. The more economical trial schedule
in the predecessor can also be used when polynomial-size levels in b
are desired. This is not a faster factoring theorem without the engine.

Using only the CRT/Linnik alternative gives the separately valid invoice
2^t poly(b)+F(b)exp(O(D))poly(b). That version preserves exponents through
one with the schedules recorded in the original checkpoint; Section 4.2
is what removes that restriction. No engine or score is inferred from
the better level bound alone.

### The faithful coefficient family preserves both directions of the score

Conversely, a smallest-prime-factor algorithm at a fixed saving exponent
beta>0 computes the faithful coefficient vectors at the same exponent,
with polynomial overhead in b,q. Fully factor N with at most b calls on
integers no larger than N, then use the predecessor's multiplicative formula

```
M_1(N)=product_(p^a || N) (1+chi_1(p)+...+chi_1(p)^a).
```

Discrete logarithms modulo the auxiliary q are acquired by scanning all
q-1 powers of its primitive generator, so their cost is polynomial in q.
Use division-free powering and geometric-sum doubling in the explicit
cyclotomic algebra; no inverse of 1-chi_1(p) is assumed. For p=q the
character is zero and the local factor is one. All coefficient and ring
costs are polynomial in b,q, as in the
[baseline engine](weight-one-divisor-coefficients.md#4-encoding-and-arithmetic-are-charged).

Thus, for each fixed positive logarithmic-saving exponent, factoring and
computing these whole faithful vectors for all prime levels q=3 modulo 4
have equivalent upper-bound targets under the stated polynomial-in-q
invoice. The same conclusion holds if the coefficient engine accepts the
certified roughness promise from the trial phase. The forward direction
needs the eligible levels only. This does not establish either task at a
stronger score or prove an unconditional complexity lower bound.

## 6. Controls, credit and stopping point

Run `make verify-sparse-faithful-character`. The
[decoder/acquisition implementation](../tools/sparse_faithful_character.py),
[independent checker](../tools/check_sparse_faithful_character.py) and
[source-pinned receipt](../evidence/sparse-faithful-character.json) retain
exhaustive kernel distances, sparse balls and bounded deterministic samples,
mixed and nilpotent orders, actual supplied divisor vectors, the sharp
collision outside the certificate, reduced-quotient information loss and
three fixed CRT acquisition controls. Supplied character observations use
the predecessor's exact **integer** cyclotomic root powers independently
of the new binary reduction. Actual divisor controls explicitly enumerate
known divisors and exclude instances outside their declared sparsity bound;
they are not implementations of the missing coefficient engine.

The panel has 19,642 nonzero kernel codewords, 1,902 sparse roundtrips,
457 actual supplied divisor vectors and 303 actual vectors explicitly
outside their fixed sparse promises. Twelve CRT-acquired and twelve
directly sieved prime levels have their eligibility checked. Forty
fixed-four CRT density controls check the local progression formula.
Three complete supplied-data CRT factorizations
(9,11,15) recover their exact divisor-square polynomials and smallest
factors from six eligible faithful-character vectors each. Their collision
and height conditions are independently checked; these small certificates
do not replace the uniform L-level bound in Section 5.

The [full default receipt](../evidence/sparse-faithful-full-verification.json)
and [complete compressed log](../evidence/sparse-faithful-full-verification.stdout.txt.xz)
pin exit-zero `make verify` at clean frozen source
`a4e607f5132e56ee34f5f5fc9864e4beabf81414`, tree
`d2d40ac03967d1d4b9d846058c6907a9b0dfd975`, in 907.08 seconds.
The run used a verified private /tmp namespace, with the host root read-only
and only the frozen checkout writable; tracked state was clean before and
after. It includes the sparse decoder, both prime-acquisition algorithms
and their controls. The later one-field-component guard and the reverse
score argument have their focused checker/link replay separately; the
decoder implementation and imported arithmetic sources are unchanged.
Later main integrations are outside that full-replay pin. No Lean target
or independent mathematical proof review is included.

The [combined default receipt](../evidence/charge-reciprocity-full-verification.json)
and [complete log](../evidence/charge-reciprocity-full-verification.stdout.txt.xz)
also cover the later arithmetic-reciprocity checker and its updated receipt:
clean frozen source `a21c955`, tree `ab36f031d27fb3dcbe41eec969bfab9d47a8a4ed`,
exit zero in 888.40 seconds. That private-namespace run includes
PR96/97/98/99 and excludes PR100 and subsequent changes. No Lean target
was run; later documentation changes have a separate focused scope.

The faithful map, signed sparse-CRT consumer, trial-budget inequality and
source attribution retain the preceding contributors' credit. Berlekamp,
Prony/Vandermonde recovery, cyclotomic distance arguments, CRT, Linnik,
the lower sieve and Bombieri--Vinogradov are classical. The present
composition, charged restricted-engine transfer
and implementation were developed by Sol worker 2 with substantial OpenAI
Codex assistance. No new Lean or independent mathematical review is claimed.

The information and level-acquisition questions are closed at this scope.
Further progress requires a faster faithful coefficient engine, a different
complete arithmetic algorithm, or another genuine change of representation.
Additional support/prime sweeps do not advance that remaining obligation.


### A field-query contract for the complete coefficient target

For the acquired q=3 modulo4 levels, a full faithful coefficient vector
can instead be supplied by polynomially many original-index queries in ONE
fixed finite field, provided the engine is uniform in a supplied primitive
generator. This requires additional coefficient queries; it does not infer
all components from one fixed observation and does not assume a semiprime.

For u=(q-1)/2>1, acquire one irreducible f of Phi_u and put omega=T in
F=F_2[T]/f, of degree nu. Choose one representative r from each doubling
orbit on the units modulo u. Every orbit has size nu. For the canonical
generator g choose an odd unit k modulo2u with k^(-1)=r modulo u, and
query the engine using generator g^k, still sending that supplied generator
to the SAME omega. Its output is P_N(omega^r). Frobenius gives all values
in that orbit because P_N has binary coefficients. The phi(u)/nu queries
therefore give its values on every primitive u-th root.

On degree-below-phi(u) binary polynomials this is a square injective binary
linear map: a polynomial in its kernel vanishes at all primitive roots,
so Phi_u divides it, and it is zero. Explicit binary elimination acquires
the inverse in polynomial bit cost in q. Apply it to recover the complete
faithful remainder, then use the earlier sparse decoder and complete
factoring reduction. Field acquisition, the changed generators, all query
outputs and the inverse are charged. Conversely, a full coefficient vector
answers each such field query by evaluation. Thus the complete fixed-saving
complexity target in Section5 is equivalently a single-field coefficient
engine uniform in the supplied primitive generator, with polynomial-in-q
overhead. This is a standard finite-field/CRT encoding observation, not a
new coefficient algorithm or an exponent gain. It applies to these reduced
e=1 levels; it does not discard nilpotent coordinates at other levels.

The [semiprime specialization](semiprime-single-field.md) is stronger in a
different direction: it uses only ONE fixed component per level and a
quadratic, with an explicit semiprime promise. A provider restricted to
that promise does not satisfy the all-input field-query contract above.
