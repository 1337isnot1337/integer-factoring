# Native integer words from Witt traces and known exact divisions

## Research question and scope

Base: `c9a4f1c`, after the precision carry and packed-system checkpoints.
Can a constant number of Teichmuller inputs represent a complete integer
factor/prefix query, if we explicitly admit known exact divisions by two?

The answer at the **pointwise circuit** level is yes. A Witt trace supplies
every native word in a field of linear extension degree. The acquired parity
polynomial supplies bit extraction, with one bit of precision charged for
each division. This gives a concrete two-input family whose nonzero-count
oracle would factor integers with polynomial overhead.

The [implemented polynomial carry backend](teichmuller-cyclic-carries.md)
does **not** contract these exact-division circuits. Their polynomial
expansion may have exponential degree. No faster counter, stronger complete
factoring bound, new Lean or independent mathematical review is supplied.
The complete score remains **1/2**.

This changes the input interface, rather than asserting that the earlier
binary-linear or direct elliptic word barriers have been bypassed within
their hypotheses. The [elliptic audit](elliptic-word-multiplication-cost.md)
concerns direct algebraic encodings; this construction charges precision
and introduces an additional operation whose summation cost is open.

## 1. Surjective native words from a trace

Let `b>=1`. Put

```
g_b=(4^b-3*2^b+2)/6.
```

Choose an odd `d` with `2^d>4*g_b^2`. The uniform choice `d=4b+1`
suffices. In `R_b=W_b(F_(2^d))`, let `T_d` be the Teichmuller set and
`sigma` its lifted Frobenius. Define

```
phi_b(U)=Tr(U)=sum_(r=0)^(d-1) sigma^r(U) in Z/(2^b).
```

**Claim.** Every element of `Z/(2^b)` occurs as `phi_b(U)`.

Here is an ordinary proof with the geometric input explicit. The
Artin--Schreier--Witt curve

```
F(Y)-Y=[X],        Y a Witt vector of length b,
```

is a smooth affine cover of degree `2^b`, with one smooth rational point
at infinity. Its smooth projective completion has genus `g_b`. These facts
follow by taking `p=2` and `f(x)=x` in Davis--Wan--Xiao,
[*Newton slopes for Artin--Schreier--Witt towers*, arXiv:1310.5311v3](https://arxiv.org/pdf/1310.5311v3),
Introduction, pages 1--3: the displayed genus formula gives
`2g_b=sum_(j=1)^b 2^(j-1)*(2^(j-1)-1)`.
Their purity statement gives the usual Hasse--Weil point bound.

For a target `a in Z/(2^b)`, take the constant `c=a/d` in the base Witt
ring; division by the odd integer `d` is a unit operation. Twist the curve to

```
F(Y)-Y=[X]-c.                                          (1)
```

Over the algebraic closure, translate `Y` by a solution of the constant
Artin--Schreier--Witt equation. This identifies (1) geometrically with the
original curve. Thus its genus is unchanged, and its unique geometric point
at infinity is again rational over `F_(2^d)`.

The trace onto the base Witt ring is surjective: on base constants it is
multiplication by the unit `d`. The kernel of `sigma-1` has size `2^b`,
and its image is contained in the trace kernel. Both the image and that
kernel have size `2^(b*(d-1))`, so they agree. Consequently an `X=x`
has rational points in (1) exactly when `Tr([x])=a`; in that case it has
`2^b` points in its fibre.

With `q=2^d`, Hasse--Weil gives at least

```
q-2*g_b*sqrt(q)>0
```

affine points on every twist. Hence every target trace occurs. The condition
is certified by the exact integer inequality `q>4*g_b^2`. For `d=4b+1`,
`2g_b<4^b<sqrt(q)`, so the choice is uniform in `b`.

The assertion concerns the trace modulo `2^b`. It does not assert that every
value occurs at arbitrary higher precision in the same extension degree,
that the encoding is injective, or that its fibre multiplicities are equal.
The existence reduction below only requires nonempty fibres.

## 2. Pointwise bit extraction with a charged operation

Use the preceding [parity polynomial](precision-efficient-residue-indicators.md),
written `H_k`, of degree `2k-1`. On base scalars modulo `2^k`, it is exactly
zero on even values and one on odd values. At available precision `k`, put

```
e=H_k(u),       u'=(u-e)/2.                             (2)
```

The numerator is known to be divisible by two for every input. The result
is well-defined modulo `2^(k-1)`, so (2) loses **one precision bit**. Repeat
it `b` times and retain the bits `e_0,...,e_(b-1)`. Then

```
X=sum_(j=0)^(b-1) 2^j e_j
```

is the canonical integer word in `[0,2^b)` congruent to the original trace.
The operation in (2) is an explicitly admitted known-exact-division gate;
two has no inverse in the coefficient ring. No such gate is silently
provided by the unit-inverse compiler.

Each `H_k` is acquired and evaluated by its integer coefficients and Horner
arithmetic, with polynomial bit cost in `k`. The whole pointwise procedure
uses polynomially many operations in the starting precision and `b`.
This is a bound on circuit construction and evaluation at a supplied point.
It is not a bound on a sum over all points.

## 3. A complete factor/prefix predicate

Let `N>=2` have bit length `b`, choose `d=4b+1`, and put

```
K=2d+1,        L=2b+1,        P=K+b+L=11b+4.            (3)
```

Work in `R=W_P(F_(2^d))`. Its Teichmuller set reduces compatibly to the
one in Section 1. There are just **two summation inputs**, `U,V in T_d`.
Compute their lifted traces and apply (2) to acquire native words `X,Y`.
After the depth-`b` extraction, retain precision `K+L`.

Compute `f=XY-N` and apply (2) for another `L` bits, leaving precision `K`.
The product of the complements of those bits is the indicator of
`f=0 mod 2^L`. Since `X,Y<2^b` and `N<2^b`,

```
abs(XY-N)<2^(2b)<2^L,
```

so this congruence is equivalent to the integer equation `XY=N`.
If `x_j,y_j` are the word bits, the flags

```
1-product_(j=1)^(b-1)(1-x_j),
1-product_(j=1)^(b-1)(1-y_j)
```

test `X>=2` and `Y>=2`. A specified most-significant-bit prefix of `X`
is tested by multiplying the corresponding bit or complemented-bit factors.
Multiply these flags and the product-zero indicator. The resulting
`A_(N,prefix)(U,V)` is exactly zero or one modulo `2^K` and tests the full
integer factor/prefix predicate.

The two operand extractions run in parallel. Their precision loss is `b`,
followed by `L` more for the product test. Thus (3) accounts for every
division depth; the number of gates is polynomial in `b`. None of the
prefix, comparison or integer multiplication semantics is delegated to an
unwritten encoding.

## 4. The exact missing evaluator and complete transfer

Request

```
C_(N,prefix)=sum_(U,V in T_d) A_(N,prefix)(U,V) mod 2^K. (4)
```

The actual count is at most `2^(2d)<2^K`. Thus its residue recovers the
nonnegative integer, and zero detects absence of a witness. Surjectivity
ensures that every native factor pair has at least one preimage; unequal
positive fibre weights cannot turn existence into absence.

An initial empty-prefix query decides whether a proper factor exists.
If it does, choose prefix bit zero whenever the corresponding count is
positive, and one otherwise. After `b` further queries this gives the least
proper divisor, which is prime. If the initial count is zero, output `N`.
This handles primes, squares and general composites without a semiprime
promise.

The field model is also charged. Deterministic construction of a prescribed
degree irreducible polynomial in fixed characteristic is polynomial time;
see Shoup, [*New algorithms for finding irreducible polynomials over finite
fields*](https://www.shoup.net/papers/detirred.pdf), Math. Comp. 54 (1990),
435--447. The finite controls below supply their small models explicitly
and do not implement that general constructor. Lift the model with the
existing unramified-ring construction; its arithmetic, Frobenius and
Teichmuller acquisition have polynomial bit cost in `d,P`. Here both are
linear in `b` and the encoded ring elements have polynomial bit length.

**Conditional transfer.** If the particular family (4), with all parameters
in (3), is evaluated uniformly in `F(b) poly(b)` bit operations, then
smallest-prime-factor computation takes `F(b) poly(b)` bit operations.
The cost parameter is the original input bit length `b`, not an unspecified
power of the circuit description length. In particular a bound
`F(b)=exp(O(b/log^alpha(b)))` would retain the saving exponent `alpha`.

The evaluator in this hypothesis remains open. Known exact division before
summation cannot be replaced by division of a final sum without a proof.
Returning the correct zero/nonzero answer for (4) is already sufficient;
the factoring reduction does not require the other count bits.
Expanding the composed selectors as ordinary polynomials can cost
exponential degree. The current norm, collected-unit-inverse and packed
residue-system algorithms do not prove the hypothesis. Further work on this
representation requires a contractor for the explicit trace/word family
with its complete cost below the desired bound.

## 5. Finite controls and attribution

The [pointwise module](../tools/witt_trace_word_circuits.py) implements trace,
Horner parity, known halving and the factor/prefix predicate. It has no
domain-summation routine. The
[checker](../tools/check_witt_trace_word_circuits.py) enumerates bounded fields
and then uses measured word fibres to check weighted prefix counts against
independent native integer predicates. It distinguishes those controls from
scalar word-pair checks for larger fields that are not enumerated. Its
[receipt](../evidence/witt-trace-word-controls.json) retains the model
certificates, fibre counts, query digest and source hashes. Run
`make verify-witt-trace-words`.

The fixed panel enumerates all 544 points of its two field models and
checks 1,056 decoded trace words across the certified word sizes. It checks
19,936 pointwise word-pair circuits, 74 prefix counts weighted by the
complete measured fibres, and 62 additional scalar-prefix controls. The
test oracle recovers the least factor on its eight finite targets using
22 prefix queries; it enumerates the small word domains and is not a fast
implementation of (4). The degree-one negative control omits words two
and three at word precision two, demonstrating why surjectivity needs
an actual certificate.

The curve input and genus calculation retain Davis--Wan--Xiao's credit;
uniform finite-field construction retains Shoup's credit. The ring,
Frobenius and Hermite parity implementation retains the preceding repository
contributors' credit and Rabinoff's Witt-vector reference. The contribution
here is the explicit surjective word interface, charged precision circuit
and conditional complete factoring transfer. Substantial OpenAI assistant
support was used. These are ordinary arguments and finite controls.

The later [pointwise-selection analysis](trace-selection-weil-budget.md)
bounds the degree/rank budget of an exact polynomial binary-trace
selector for this native query. Its full fibre rectangles also obstruct
direct low-degree ordinary equations without auxiliaries. These are
restrictions on those acquisition interfaces, not a lower bound for
the admitted exact-division circuit family or a structured contractor.
