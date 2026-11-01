# Weil bounds and the native-word pointwise-selection budget

## Result and scope

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338238)
asks for a correct tradeoff between the **number** and **degree** of binary
polynomial trace predicates selecting a thin nonempty set. Integer Fourier
orthogonality, Artin--Schreier reduction, an exceptional-fibre polynomial
and the classical Weil bound give such a tradeoff.

Over F_q, q=2^w, let S be the nonempty solution set of r predicates
`Tr(f_i(U))=b_i` on the full v-variable field domain. Suppose their
ordinary polynomial degrees are at most D>=1. If

~~~
|S|/q^v <= 2^(-h),       0<=h<=w/2,
~~~

then

~~~
r + log_2(D(D+1)) >= h.                                (1)
~~~

For the [native Witt-word interface](witt-trace-native-word-circuits.md)
with word size b and w=4b+1, a nonempty proper-factor query has density
at most `8*2^(-3b/2)`. Consequently any exact polynomial-trace selector
for it obeys

~~~
r + log_2(D(D+1)) >= 3b/2-3.                            (2)
~~~

The same density argument allows extra full-F_q auxiliary variables if
every accepted tuple projects to a native factor pair and the accepted
set is nonempty. It is not a statement about arbitrary mixed domains or
arbitrary existential encodings without those hypotheses.

There is also an elementary degree obstruction for the **full pointwise
two-input predicate** described by ordinary field equations without
auxiliaries: some nonzero equation must have degree greater than 2^(3b).
Section 5 proves this from the Cartesian product of two native-word fibres.

These results concern acquiring the stipulated pointwise input for the
current trace or residue-equation interfaces. They are **not runtime lower
bounds** for their implementations or for factoring. Weighted cancellation
formulas, different word maps, non-polynomial exact-division circuits,
different auxiliary domains and structured global evaluators remain open.
The complete factoring score stays **1/2**.

## 1. Reduce formal Artin--Schreier classes, not finite-field functions

For a nonconstant monomial `c U^alpha`, divide all exponents by their
largest common power 2^k and apply inverse Frobenius to c. The binary
trace is unchanged. Every remaining exponent vector has at least one odd
coordinate; combine equal vectors by addition in F_q.

Remove constants separately and XOR their **trace bits** into the target.
Do not place a target bit under the trace, which would lose it at even w.
These operations are linear over F_2 and do not increase total degree.

Binary elimination on the nonconstant coefficient vectors yields R<=r
independent **formal AS classes**. Since S is nonempty, every relation
among those classes has a consistent target after constant traces are
removed. Thus R basis predicates give the same set S, still at degree
at most D. Every nonzero binary combination of the basis has a nonzero
formally reduced polynomial.

This is not necessarily the exact rank of functions on a small field.
For example, `XY^2+X^2Y` is formally reduced but its trace is identically
zero on F_4^2. The character estimate below remains valid in that case,
with an exceptional set large enough to make its bound uninformative.

The [sparse implementation](../tools/trace_artin_schreier.py) performs
this reduction and elimination in polynomial time in the written input,
coefficient-field degree and exponent bit lengths. It does not enumerate
the domain. It can prepare valid nonlinear trace constraints for a caller
of the existing compiler; it does not itself contract them or silently
change the existing API's output precision.

## 2. A nonzero polynomial detects exceptional fibres

Let g be any nonzero formally reduced combination, of degree at most D.
Choose a variable Y appearing in g, and write the remaining variables T:

~~~
g(T,Y)=a_0(T)+sum_(j>=1) a_j(T)Y^j.
~~~

Group j=m*2^k by its odd part m. Choose the largest occurring odd part M,
and let K be the largest k in that group. For a fixed T=t, univariate
AS reduction makes the coefficient of Y^M equal to

~~~
C(t)=sum_k a_(M*2^k)(t)^(2^(-k)).
~~~

Inverse Frobenius here is an operation on the coefficient **values** in
F_q, not an uncharged polynomial root. Raise to 2^K instead:

~~~
A(T)=sum_k a_(M*2^k)(T)^(2^(K-k)).                      (3)
~~~

This is an ordinary polynomial, of degree at most D*2^K<=D^2, and
`A(t)=C(t)^(2^K)`. It is nonzero. Indeed, cancellation of monomials from
two different k would require their full original exponent vectors,
including the Y coordinate, to be positive powers of two of one another.
Both vectors were already globally reduced, so this is impossible.
Within a coefficient polynomial Frobenius is injective as well.

Except where A(t)=0, the reduced fibre has positive odd degree exactly
M<=D. Schwarz--Zippel bounds the exceptional tuples by
`degree(A)*q^(v-2)` when v>=2. With v=1, A is a nonzero constant and
there are no exceptions. A nonzero formal A may vanish at every point
of a small finite field; the degree bound correctly permits that.

## 3. Weil, integer Fourier projection and density

For a univariate polynomial of positive odd degree M in characteristic
two, the curve `Z^2+Z=f(Y)` has genus (M-1)/2 and one rational point at
infinity. The Hasse--Weil bound therefore gives

~~~
abs(sum_(Y in F_q) (-1)^Tr(f(Y))) <= (M-1)*sqrt(q).
~~~

For primary-source context, the curve, unique infinity, genus formula
and purity are in Davis--Wan--Xiao,
[*Newton slopes for Artin--Schreier--Witt towers*, Introduction](https://arxiv.org/pdf/1310.5311v3),
at tower level one and p=2. A nonmonic leading coefficient can be normalized
over the algebraic closure without changing the geometric genus. The
projective point count is `q+1` plus the displayed binary character sum.

Apply this on good fibres, and use the trivial bound q on each bad fibre.
The triangle inequality gives the sharper witness-specific estimate

~~~
abs(sum_U (-1)^Tr(g(U)))
 <= (M-1)q^(v-1/2)+degree(A)q^(v-1)
 <= (D-1)q^(v-1/2)+D^2 q^(v-1) = E.                   (4)
~~~

Now use ordinary integer sign characters. Fourier orthogonality over
the R basis predicates gives

~~~
2^R |S| = q^v + sum_(lambda!=0) signed_character_sum(lambda).
~~~

Each nonzero combination has bound (4). It follows that

~~~
2^(-R) <= |S|/q^v +(D-1)/sqrt(q)+D^2/q.                (5)
~~~

If the density is at most 2^(-h), h<=w/2, both q^(-1/2) and q^(-1)
are at most 2^(-h). The right side is therefore at most
`D(D+1)*2^(-h)`, proving (1) first for R and then for r>=R.

This argument takes place over the integers before any reduction modulo
two. It does not use the
[withdrawn ANF degree assertion](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338124).
The latter omitted r and changed coefficient characteristics incorrectly.
Its F_8 counterexample, with two degree-one predicates selecting {3,6},
is preserved in the new receipt and is consistent with (1).

## 4. The native Witt factor set is thin

Retain the certified word map

~~~
phi_b(U)=Tr_W([U]) mod 2^b,       q=2^(4b+1).
~~~

Its ASW twist has genus `g_b=(4^b-3*2^b+2)/6`. Each word a has fibre
size M_a satisfying

~~~
(q-2g_b sqrt(q))/2^b <= M_a <= (q+2g_b sqrt(q))/2^b.
~~~

In particular the certified Hasse gap implies `0<M_a<2q/2^b`.
A proper-factor query for a b-bit integer N accepts at most tau(N)
ordered native pairs. Its full two-input density is at most

~~~
sum_(xy=N, x,y>=2) M_x M_y /q^2
 < 4*tau(N)/2^(2b)
 <= 8*2^(-3b/2),
~~~

using only the elementary divisor bound tau(N)<=2sqrt(N) and N<2^b.
Prefixes can only reduce this density. For a nonempty query take
`h=3b/2-3<=w/2` in (1), proving (2). An empty prefix is nonempty on
every composite because every native factor pair has positive fibres.

Extra full-field auxiliary variables cannot increase the density if
every accepted tuple projects to an accepted native pair: at most
q^(v-2) auxiliary assignments lie above each pair. Thus the same argument
applies to any nonempty such relation with v>=2. This does not cover
different domains or a representation with signed cancellations.

In the [Artin--Hasse compiler](artin-hasse-trace-projection.md), its degree
base K+2 is at least D+1 after legitimate preprocessing. Its declared
exponential parameter consequently satisfies

~~~
r+v log_2(K+2) >= r+2log_2(D+1) >= 3b/2-3.
~~~

The sufficient sublinear-parameter criterion fails for this pointwise
selection interface. An upper-bound parameter being large is **not** a
lower bound on every implementation's runtime. A new contraction estimate
or a different arithmetic representation remains possible.

## 5. Ordinary field equations must vanish on large fibre rectangles

For this uniform choice of q, `g_b<4^b/6` gives `q>16g_b^2`. Hence

~~~
M_a > q/2^(b+1) = 2^(3b)                               (6)
~~~

for every native word a. Fix any proper factor pair x*y=N. The full
pointwise native predicate accepts the entire rectangle
`phi_b^(-1)(x) times phi_b^(-1)(y)`.

A nonzero ordinary polynomial F(U,V) over F_q, or any coefficient-field
extension, cannot vanish on A times B if its total degree is smaller
than min(|A|,|B|). For each u in A, the root bound first makes F(u,V)
identically zero; the same root bound on all its U coefficient polynomials
then makes F identically zero, a contradiction.

If an ordinary polynomial-equation system without auxiliaries defines
the full native predicate, every equation vanishes on that rectangle.
At least one equation is nonzero, since the predicate is not the full
domain. Equations (6) and the root bound force its degree to exceed
2^(3b). This blocks direct low-degree acquisition of this **full predicate**
for the ordinary residue-system interface. It does not exclude a sparse
witness cross-section, another word map, auxiliary projection, weighted
cancellation or a non-polynomial division-circuit evaluator.

## 6. Exact controls and remaining work

The [checker](../tools/check_trace_selection_weil.py) and
[receipt](../evidence/trace-selection-weil-controls.json) retain 25 character
cases, 5,808 pointwise AS reductions, 436 fibres, 12 trace systems with
480 pointwise equivalences and eight integer Fourier identities. They
check the sharper form of (4) by squared integer comparisons, without
floating square roots. The F_4 exceptional polynomial X+X^4 vanishes at
every field point and is explicitly retained; no false function-rank
claim is made. Even-width constants, Frobenius coefficient roots,
subfields, nonlinear redundancy, contradictions and invalid inputs are
included. Three rectangle controls exhaust 4,672 subthreshold polynomials
and exhibit nonzero examples at the exact degree threshold.

Run make verify-trace-selection-weil. The finite controls verify the
implementation and conventions; the published curve bound and ordinary
slicing/root arguments supply the global statements. No native large
field or factor query is being exhaustively enumerated here.

The word-map geometry retains the preceding contributors' and
Davis--Wan--Xiao's credit. The earlier withdrawn mechanism and counterexample
are credited through the linked coordination report. AS reduction, Weil
and Schwarz--Zippel are standard tools; the application here was developed
with OpenAI assistant support, with no new Lean or independent mathematical
acceptance. The missing efficient native exact-division contractor and
representations outside these pointwise interfaces remain open.
