# Known polynomial product fibres and rough integer inputs

There is a classical alternative to enumerating every pair of digit
settings **when their complete product polynomial is supplied**. Factor
the primitive polynomial over Q, enumerate its polynomial divisors, and
filter their coefficient boxes. For a primitive degree-D product the
cost is `2^O(D) poly(input bits)`, rather than Cartesian enumeration in
all coefficient bits. After actual trial division below the radix,
every accepted digit product is primitive.

This gives a concrete local solver outside the literal-enumeration
hypothesis of [half-top-structure](half-top-structure.md). It does not
choose the unknown carries, reduce the Appendix B pad, or provide a
faster global contraction. The complete factoring exponent remains
**1/2**. The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19366648)
and [roughness refinement](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19366707)
distinguish those obligations.

## 1. Exact factor-box algorithm, including the content

Let `W in Z[X]` be nonzero of degree D, with coefficients of at most
H bits. Let A and B range over explicit finite integer coefficient
boxes; call their full written bit length S.
Trailing zero coefficients are permitted in each box.

Compute the positive content `c=gcd(coefficients(W))` by integer gcds,
and the sign epsilon of W's leading coefficient. Classical Gauss
factorization gives

```math
W=\varepsilon c\prod_{i=1}^s f_i^{e_i},
\qquad f_i\in\mathbb Z[X]\text{ primitive, irreducible, positive leading},
\qquad\sum_i e_i\deg f_i=D.
```

The deterministic polynomial-time acquisition of the primitive factors
is the theorem of [A. K. Lenstra, H. W. Lenstra and Lovasz (1982)](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Lovasz/LovaszLenstraLenstrafactor.pdf),
Theorem3.6 (PDF page18), with primitive input explicitly required.
Its bit cost is polynomial in D and `log ||W||_2=O(H+log(D+1))`.
The source PDF SHA256 is
`dabefb8bcfb5dbb8b36a43745081f8ab6f18aadc85b1e252fd7c85abefc5f704`;
the introduction and theorem were visually checked. This imported
theorem factors the polynomial over Q; it does not factor c over Z.

Every ordered integer polynomial factor pair is uniquely of the form

```math
A=a\prod_i f_i^{u_i},\qquad
B=\frac{\varepsilon c}{a}\prod_i f_i^{e_i-u_i},
\qquad a\in\mathbb Z\setminus\{0\},\quad |a|\mid c,
\quad0\le u_i\le e_i.
```

This follows by taking the positive-leading primitive part of A and
applying Gauss's lemma and unique factorization in Z[X]. Enumerate these
pairs, pad their coefficients with zeros, and test the supplied boxes.
Integer-evaluation prefix, range or proper-divisor filters can also be
tested on each candidate with polynomial bit cost.

There are at most

```math
2\tau(c)\prod_i(e_i+1)\le2\tau(c)2^D
```

candidates. The last inequality uses `e+1<=2^e` and `deg f_i>=1`.
The LLL output has polynomial bit length. Multiplying at most D of
those factors and an H-bit scalar, evaluating them at a written radix,
and testing their boxes all have polynomial bit cost. No unit-cost
large-integer arithmetic is used.

For `c=1` this proves the stated `2^O(D) poly(D,H,S)` algorithm. For
general c its integer factorization must be acquired and paid for.
The existing [complete half-rate machine](half-rate-factoring.md),
applied recursively to this at-most-H-bit content, costs
`exp(O(H/sqrt(log(H+2)))) poly(H)`, including full factorization by
polynomially many smallest-factor calls. Divisor enumeration costs
`tau(c) poly(H)` once that factorization is known.

An elementary uniform estimate suffices here:
`tau(c)<=exp(O(H/log(H+2)))`. For large H, split prime divisors at
`z=floor(H/(log H)^2)`. There are at most z small primes and each
exponent is O(H); their total contribution to `log tau` is
`O(z log H)=O(H/log H)`. For the large primes, the sum of their
exponents is `O(H/log z)=O(H/log H)`, and `log(e+1)<=e log2`.
The bounded small-H range is absorbed in an absolute constant.

Thus the fully acquired general nonzero fibre costs

```math
\exp\!\left(O\!\left(D+\frac{H}{\sqrt{\log(H+2)}}\right)\right)
\operatorname{poly}(D,H,S).
```

This bound inherits its content term from the existing factoring
theorem. The degree-zero case already contains ordinary integer
factorization, so it must not be described as a polynomial-time
general integer-content algorithm.

If W=0, the integral-domain property says A=0 or B=0. With box sizes
`a,b` and zero-membership indicators `a0,b0`, the count is exactly
`a*b0+b*a0-a0*b0`, computed directly. This handles all zero settings
without invoking a factorization of zero.

## 2. A known carry polynomial specifies one such fibre

Let the radix be `R=2^w`, use d digits for each nonnegative operand,
and T=2d positions for their product. Put `K=d(R-1)`. Write the
base-R target polynomial as `N_R(X)=sum_(j<T)n_j X^j`.
For an internal carry vector `r_1,...,r_(T-1)` in `[0,K]`, put

```math
C(X)=\sum_{j=0}^{T-2}r_{j+1}X^j,\qquad
W_C(X)=N_R(X)+(R-X)C(X),\qquad r_0=r_T=0.
```

The coefficient equations are exactly

```math
a_j+r_j=n_j+Rr_{j+1}.
```

Consequently every digit-polynomial factorization `AB=W_C` has
`A(R)B(R)=N`. Conversely every operand pair with product N has one
unique carry vector and hence one unique such fibre. The bound K and
uniqueness follow from the carry recurrence and
`a_j<=d(R-1)^2`; these are the same integer hypotheses as
[the adaptive constructor](adaptive-transform-constructor.md#3-carry-bounds-and-the-exact-factor-assignment-bijection)
and [the weighted carry reformulation](weighted-dyadic-carry-contraction.md).
Different carry vectors give different W because multiplication by
`R-X` is injective in Z[X].

Knowing the single value `W(R)=N` does not specify W. For example, at
R=4, N=12 both `W=3X` from operands3,4 and `W=4+2X` from operands2,6
are valid. Their carry polynomials are0 and1, respectively. Factoring
only the target digit polynomial would miss the second fibre.

## 3. Actual roughness removes the scalar-content call

Suppose N>0 has no prime divisor below R, as established by an actual
trial pretest. If A is a nonzero digit polynomial, its positive content
is at most R-1 and divides A(R). For an accepted factor pair it therefore
divides N. Any content greater than1 has a prime divisor below R, a
contradiction. Thus A and B are primitive, and Gauss's lemma gives
`content(W_C)=1` for every accepted fibre.

Nonprimitive W_C can therefore be rejected without factoring its
content under this promise. Every remaining known-W fibre has the
`2^O(d) poly(input bits)` primitive algorithm. This promise is not
obtained merely by supplying a trial cutoff parameter. The complete
pretest costs `2^w poly(b)` bit operations for a b-bit target.

Choosing `d=ceil(b/w)` includes every possible proper factor, with
filters excluding1 and out-of-range operands. The two local costs
`2^w` and `2^O(b/w)` would balance near w=sqrt(b). **They are not the
cost of a complete factoring algorithm:** the carries remain unknown.
The explicit carry domain has up to `(K+1)^(T-1)` vectors; listing it
is still charged. This is an enumeration bound, not a lower bound
against every possible collective carry algorithm.

## 4. Source placement and the remaining operation

The adaptive source's full inverse transform gives canonical coefficient
residues. Since `M>dR^2`, they are exact integer coefficients, so all
of them together can supply W to this algorithm. Partial Fourier
samples, residues without that bound, and the value at X=R do not.
Contracting a whole transform into this fibre function does not
automatically supply its dense table over all boundary coefficients.
It also does not change the leading pad in the existing evaluator.

### Generic factorization does not select the exceptional carries

This distinction can be made exact. Leave C0 free in

```math
F(X,\mathbf C)=N_R(X)+(R-X)\sum_{j=0}^{T-2}C_jX^j,
\qquad N_R(R)=N\ne0.
```

Then F is irreducible over `Q[C_0,...,C_(T-2),X]`, and hence over
`Q(C_0,...,C_(T-2))[X]`. Indeed regard it as a polynomial in C0 over
the UFD `A=Q[X,C_1,...]`. Its coefficient is R-X. The remaining
term reduces to the nonzero rational constant N modulo R-X, so the
two coefficients are coprime in A. Thus F is primitive and linear in
C0; it is irreducible over the fraction field of A, and Gauss's lemma
applies. Its positive X-degree permits localization to `Q(C)[X]`.
The same argument holds after fixing any other carry coefficients,
provided C0 remains free. This is an elementary Gauss-lemma application,
not a new general irreducibility theorem.

Nevertheless useful specializations exist. At R=4,N=35, fix the unused
top carry to zero and consider the degree-two family

```math
F(X,U,V)=3+2X^2+(4-X)(U+VX).
```

It is generically irreducible, but U=0,V=1 gives
`3+4X+X^2=(X+1)(X+3)`. Both factors have digits below4 and evaluate
to5,7 at X=4; N has no prime divisor below4. The specialization also
preserves degree two, so this phenomenon is not just a leading-degree
drop.

The [existing orbit-gcd source audit](interval-orbit-gcd.md#4-what-existing-algebraic-results-actually-supply)
already records the constant-depth, polynomial-total-degree Q-circuit
hypotheses of [the 2025 factorization theorem](https://arxiv.org/html/2504.08063v2#S5.SS1).
That result factors a supplied polynomial and preserves generic factors
under its variable reduction. It does not find exceptional integer
specializations of an irreducible family. No bit-cost or characteristic-two
transfer from that theorem is assumed here. This is the carry-family
instance of the generic-specialization distinction in the earlier note,
not a lower bound against an algorithm that actively searches for those
specializations.

At degree two the arithmetic selection problem has a useful exact form.
If `W=aX^2+bX+c`, a!=0, and W(R)=N, its discriminant satisfies

```math
b^2-4ac=(b+2aR)^2-4aN.
```

This follows by substituting c=N-aR^2-bR. A rational factorization
requires this integer discriminant to be a square; integrality and
digit bounds must still be checked. It is a difference-of-squares/norm
search across the carry choices, not a supplied generic factor. No
faster collective norm search is supplied by this identity.

The next relevant question is a **collective** carry selection or sum
using this primitive-fibre structure, with every query, intermediate
size and deterministic bit operation counted. A method making Q fibre
queries must pay all Q calls; cheap individual coefficients are not a
free evaluation of their generating polynomial. The present result
does not provide that operation or claim a new uniform score.

## 5. Finite evidence and credit

```sh
make verify-polynomial-digit-fibres
```

The [standard-library checker](../tools/check_polynomial_digit_fibres.py)
replays fixed full boxes `(R,d)=(2,4),(4,3),(8,2)`: 8448 literal ordered
products, 11331 full/prefix-box fibre checks, 8448 exact carry identities
and891 rough-product primitive controls. It also checks a signed
two-coefficient box, including zero and negative products.

The [compressed factor certificates](../evidence/polynomial-digit-fibres.factorizations.json.xz)
were acquired with SymPy1.14 in an isolated temporary venv. The replay
does not import SymPy: it verifies primitive contents, exact integer
products, and irreducibility of every supplied degree1/2/3 factor by the
rational-root criterion, then enumerates the proved complete factors.
The [receipt](../evidence/polynomial-digit-fibres.json) retains the exact
counts, source/data hashes and acquisition scope. This finite acquisition
is not an implementation or benchmark of the imported LLL theorem.

Polynomial factorization, Gauss's lemma, divisor enumeration, the
integer carry correspondence and the existing factoring bound retain
their classical and parent-repository credit. The contribution here
is the explicit factor-box invoice, roughness consequence and source
interface/certificate check, developed by Sol worker 2 with OpenAI
Codex assistance. No novelty claim for the classical algorithms, new
Lean or independent mathematical acceptance is made.
