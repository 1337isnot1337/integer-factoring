# Typed rational scalar crossings reduce to strict crossings without a pressure change

The later [approximate extension](approximate-rational-strictification.md)
also treats other finite poles when every paired-type word is globally
Laurent-regular. It uses weak potentials, Frobenius rounding and a bounded
charged correction; negative pressure witnesses survive, although finite-D
pressure is no longer identical. The exact theorem below retains its
original endpoint-pole hypotheses.

## Result and scope

Let K be a finite scalar rational crossing kernel over a field of positive
characteristic. Assume explicitly that its entries belong to
`k[w,w^(-1),(1-w)^(-1)]`, so its only possible finite poles are at zero
and one, including over the algebraic closure. Suppose every paired
full-letter-type restriction of every tensor power
has the required crossing endpoint at one. Then an explicit ramification
and paired physical-leg gauge turn K into a strict Laurent-polynomial
crossing. The complete two-place Smith profiles and physical type energies
scale by the same factor, so the full pressure is exactly unchanged.

Consequently the [universal binary Smith floor](universal-binary-smith-floor.md)
also excludes every scalar rational binary kernel satisfying these
all-length type hypotheses. Local endpoint poles do not escape that floor.
This strengthens the earlier triangular polynomial endpoint corollary.

For binary alphabets in characteristic two there is a uniform finite
version: checking the paired-type endpoint through tensor length four is
sufficient, and four is necessary in general. Twelve primitive charge
circuits give the complete test. The substitution w=z^8 always suffices
for the strictification, followed by an integer physical-leg gauge.

The statement does not cover auxiliary matrix kernels whose endpoint
appears only after a boundary contraction, a selected sequence of typical
types, arbitrary other finite poles, or all characteristics. Larger
alphabets remain open for a pressure improvement even after strictification.
No stronger factoring exponent is asserted here.

## 1. All paired-type endpoints are exact cycle inequalities

Let the paired alphabets be A for X,V and B for Y,U. For a local entry
lambda=(x,y,u,v), define its integer charge vector

\[
 a_\lambda=(e_x-e_v,\ e_y-e_u).
 \tag{1}
\]

One coordinate in each block may be omitted, leaving dimension
d=|A|+|B|-2. Call an entry matched if x=v and y=u. Write chi=1-w,
and let nu_lambda be the chi valuation of a nonzero rational entry.
Identically zero entries can be omitted from every support calculation.

A tensor-word coefficient is a single product of local scalar entries.
Its two paired word types agree precisely when the sum of the vectors
(1) is zero. The crossing coefficient at the endpoint is one if every
local entry is matched, and zero otherwise.

The one-copy endpoint condition therefore requires each matched entry to
be regular at one and to have value one. Such entries have valuation zero.
After removing matched factors from a word, the remaining endpoint
conditions are exactly

\[
 \sum_\lambda n_\lambda\nu_\lambda>0
 \quad\text{whenever}\quad
 n_\lambda\in\mathbb Z_{\ge0},\quad
 \sum n_\lambda a_\lambda=0,\quad
 \sum n_\lambda>0,
 \tag{2}
\]

where the sums use only nonzero offmatched entries. A negative total
valuation would be a pole; total valuation zero would give a nonzero
endpoint coefficient where crossing requires zero. Products of nonzero
scalars cannot cancel. Conversely, positive total valuation gives
regularity and a zero endpoint. Thus (2), together with the matched
one-copy conditions, is equivalent to every paired-type endpoint at every
length. This equivalence itself is valid in every coefficient field.

## 2. Strict separation supplies one common paired potential

**Lemma.** For finitely many integer charge vectors a_i and integer
valuations nu_i satisfying (2), there is a real vector h such that

\[
 \nu_i+h\mathbin{\cdot}a_i>0\qquad\text{for every }i.
 \tag{3}
\]

Here is a direct separation argument. Let C be the convex hull of the
finitely many points (a_i,nu_i), and let R be the ray {(0,t):t<=0}.
They are disjoint. Indeed, if they intersect, minimizing the last
coordinate over the rational polytope

\[
 \lambda_i\ge0,\quad\sum\lambda_i=1,
 \quad\sum\lambda_i a_i=0
\]

gives a rational minimizing vertex with nonpositive value. Clearing its
denominators contradicts (2).

The compact set C has positive distance from the closed ray R. A
hyperplane strictly separates them, giving a vector (h,beta) with

\[
 h\mathbin{\cdot}a_i+\beta\nu_i>0\quad\text{for every }i,
 \qquad\beta\ge0.
\]

The sign condition is forced by boundedness of that functional above
on R. If beta is zero, replacing it by a sufficiently small positive
number preserves all the finitely many strict inequalities. Divide by
beta to obtain (3). This also treats supports with no nonzero balanced
cycle at all.

Write the coordinates of h as anchored potentials f on A and g on B.
Then (3) is

\[
 \nu_{xyuv}+f(x)+g(y)-g(u)-f(v)>0.
 \tag{4}
\]

The finite system is open. For a field of characteristic p, approximate
f and g closely enough by numbers in Z[1/p] that every inequality stays
strict. Choose D=p^k clearing those denominators. Both Df and Dg are
integer-valued. No extension of the coefficient field is needed.

## 3. Frobenius ramification removes every endpoint pole

Substitute w=z^D and define

\[
 K'_{xyuv}(z)=
 (1-z)^{D[f(x)+g(y)-g(u)-f(v)]}K_{xyuv}(z^D).
 \tag{5}
\]

In characteristic p,

\[
 1-z^D=(1-z)^D.
 \tag{6}
\]

Each nonzero offmatched entry of K' therefore has positive valuation
at z=1 by (4). Each matched entry has gauge exponent zero and endpoint
one. Thus K'(1) is the literal crossing.

The original rational entries had no finite poles except zero and one.
After (6), these remain the only possible finite pole locations. The
preceding valuation argument removes all poles at one, leaving only
possible powers of z in the denominator. Hence K' is a Laurent-polynomial
tensor. A common monomial can clear its negative z powers if a polynomial
representative is wanted; that clearing has its usual exact two-place
cost cancellation.

For any retained paired-type word, the product of the gauge factors
in (5) is one. This is consistent with the original endpoint theorem,
but the new tensor also satisfies the crossing endpoint before any
type restriction.

## 4. The full two-place pressure is identical

Let h_(0,l) and h_(infinity,l) be arbitrary physical-letter height
functions for the original kernel. Define the four gauge potentials

\[
 e_x=f,\quad e_y=g,\quad e_u=-g,\quad e_v=-f.
 \tag{7}
\]

Give K' the heights

\[
 h'_{0,l}=D h_{0,l},\qquad
 h'_{\infty,l}=D h_{\infty,l}+D e_l.
 \tag{8}
\]

At zero, 1-z is a unit, so multiplication by its physical-leg powers
is multiplication by invertible row and column matrices over the local
ring. At infinity, its order is minus one. The added heights in (8)
cancel this order, leaving row and column units there as well.
Substitution z^D multiplies every original valuation by D. Consequently,
for either complete cost-cut profile,

\[
 \sigma'_0=D\sigma_0,\qquad
 \sigma'_\infty=D\sigma_\infty.
 \tag{9}
\]

For a matched X,V symbol, e_x+e_v=0; similarly e_y+e_u=0. The two
physical type energy lists therefore satisfy

\[
 E'_X=D E_X,\qquad E'_Y=D E_Y.
 \tag{10}
\]

Evaluating the new pressure at scale t/D gives exactly the old pressure
at t. This preserves the complete type-optimized criterion, not merely
one minimum order, one rank, or the difference of determinant degrees.
It also preserves a strict negative-pressure witness if one existed.

Apply the universal strict binary floor to K'. It has positive critical
pressure for all of the compensated heights. Equations (9)--(10) prove
the same positive pressure for K. Thus the entire stated scalar rational
binary class is excluded from a post-half witness by this criterion.

## 5. A finite endpoint test at every fixed alphabet

The normalized balanced distributions from Section 2 form a rational
polytope. Its vertices correspond to positive circuits of the charge
columns. A circuit whose columns have rank r has r+1 support entries.
Its primitive integer multiplicities are bounded by the absolute values
of r-by-r minors of those columns.

Every column of (1), even after one coordinate per block is removed,
has Euclidean norm at most two. The determinant bound from the product
of column norms therefore makes every primitive multiplicity at most
2^r. Since r<=d, each primitive circuit has total length at most

\[
 L_d=(d+1)2^d.
 \tag{11}
\]

If any balanced distribution has nonpositive valuation, a minimizing
vertex does too. Its primitive integer circuit gives a failing word of
length at most (11). Thus all paired-type endpoints can, in principle,
be checked through that finite length. The bound is deliberately coarse;
the potential certificate (4) is generally a much smaller proof object.

This finite test concerns scalar endpoint identities. It does not
classify the two-place profile cost or solve the factoring construction
problem at larger alphabets.

## 6. Binary endpoints: twelve circuits and a sharp length four

For binary paired alphabets, each nonzero offmatched charge is one of
the eight directions in {-1,0,1} squared, with the origin removed.
For each direction retain the smallest valuation among its nonzero
entries. Missing directions can be omitted.

There are exactly twelve primitive positive circuits:

| Circuit kind | Number | Representative directions | Multiplicities | Word length |
|---|---:|---|---|---:|
| Opposite pair | 4 | a,-a | 1,1 | 2 |
| Three-direction triangle | 4 | (1,0),(0,1),(-1,-1), and quarter-turns | 1,1,1 | 3 |
| Three-direction weighted triangle | 4 | (1,1),(-1,1),(0,-1), and quarter-turns | 1,1,2 | 4 |

Each applicable weighted valuation sum must be positive. These are
necessary by Section 1. They are sufficient because every normalized
balanced distribution is a convex combination of the circuit vertices.
Thus checking all paired-type tensor powers through length four is
equivalent to checking all lengths.

Four cannot be replaced by three. Put every matched entry equal to one,
set all other entries to zero except

\[
 K_{1,1,0,0}=\chi,\qquad
 K_{0,1,0,1}=\chi,\qquad
 K_{0,0,1,0}=\chi^{-1}.
 \tag{12}
\]

The available offmatched charges are (1,1),(-1,1),(0,-1). Their first
balanced word uses multiplicities 1,1,2. It has length four, and its
product is one, whereas its crossing coefficient should be zero. No
nonempty offmatched balanced word of length at most three exists. The
kernel therefore passes all the shorter endpoint checks and fails at
four exactly as claimed.

## 7. A uniform eighth-power construction in characteristic two

For a valid binary kernel, every primitive circuit has an integer
valuation sum at least one and word length at most four. Its normalized
mean valuation is at least one quarter. Convex combinations give the
same lower bound for every normalized balanced distribution.

The non-strict alternative theorem, or linear-programming duality,
therefore gives a feasible potential polygon

\[
 \nu_i+a_i\mathbin{\cdot}h\ge\frac14
 \qquad\text{for all present directions }i.
 \tag{13}
\]

If there is no balanced distribution, the strict separation argument
gives an arbitrary margin, so (13) is still feasible. Intersect the
nonempty polygon with a sufficiently large integer-coordinate box.
At a vertex of that bounded polygon, two independent active normals
have entries in {-1,0,1}, so their determinant is plus or minus one or
two. All right sides are quarter-integral. The resulting potential is
therefore eighth-integral.

Take D=8 in (5). Every offmatched order is then an integer at least two.
The implementation constructs such a potential by enumerating exact
intersection vertices in successively enlarged boxes. All arithmetic is
rational, and the circuit theorem proves termination. The factor eight
is a uniform sufficient choice; no optimality claim is made for it.

## Verification

```sh
make verify-binary-typed-endpoints
```

The exact constructor is
[`binary_typed_endpoint_gauges.py`](../tools/binary_typed_endpoint_gauges.py).
Its checker verifies all twelve circuits against the complete list of
balanced multisets through length four, all 6,561 complete valuations
in the stated finite test box, and 600 partial-support valuation tests.
The sharp four-copy counterexample is retained as a negative control.

Sixteen literal rational F4 tensors are then converted to strict tensors.
The checker verifies all 256 strict endpoint entries, 64 complete native
profiles, exact scaling of both physical energy lists, and independent
binary Smith calculations with 28 block-Toeplitz rank checks. The record
is [`binary-typed-endpoint-gauges.json`](../evidence/binary-typed-endpoint-gauges.json).

The proof closes the scalar binary endpoint-pole route for this two-place
criterion. It also gives a constructive normal form for scalar typed
kernels on larger alphabets, where a strict pressure improvement remains
an open construction problem.
