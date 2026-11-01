# Exact trace projection with Artin--Hasse characters

## Question and result

Base: `21f0790`. The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19337728)
asks whether classical splitting functions remove the growing precision
from the exponential dependence on the number of binary trace predicates.
They do, with an explicit precision cost for Fourier projection.

Use a supplied unramified ring R_P=GR(2^P,d), with w dividing d, t dividing
w, C=w/t, and v Teichmuller variables over F_(2^w), optionally punctured.
There are r explicit ordinary polynomial predicates

~~~
Tr_(F_(2^w)/F_2) f_i(U) = b_i,    b_i in {0,1}.
~~~

Only coefficient residues must lie in F_(2^w). Write D for the sum of the
total monomial degrees over the written nonzero residue terms of all f_i.
This is an input parameter, not the degree of a free arithmetic circuit.
Scalar synchronized numerator and total-unit-inverse weights have degree
totals S_A,S_B, as in the [canonical rational compiler](precision-efficient-residue-indicators.md).
Choose output precision H>=1 and **supply input precision P=H+r**. Put

~~~
B = S_A + ((2P+1)d-1)S_B + (P-1)wD.
~~~

The implemented deterministic cost for the weighted constrained sum is

~~~
2^r * (B+2)^O(v*w/t) * poly(input_length,P,d,w,v,r).       (1)
~~~

At t=w the full conjugacy orbits give the sharper bound

~~~
K = S_A + (2P+1)(d/w)S_B + (P-1)D,
cost = 2^r * (K+2)^O(v) * poly(input_length,P,d,w,v,r).    (2)
~~~

The output is modulo 2^H. With weight one and H>v*w it is the exact integer
count. The ring model, coefficient embeddings and their construction remain
charged. Nonlinear predicate lists use their written length r here; no
unproved rank reduction is used. This replaces the
[divided-power lift](divided-power-trace-lifting.md)'s `(2H)^O(r)` cost by
`2^r` queries at increased precision. Both implementations are retained.
No native integer digit/division evaluator or stronger factoring bound is
supplied. The complete factoring score remains **1/2**.

## 1. Standard splitting function, explicitly acquired

The classical Artin--Hasse exponential is

~~~
E_2(X) = exp(sum_(j>=0) X^(2^j)/2^j) = sum_(n>=0) a_n X^n.
~~~

Its integrality follows from the formal identity

~~~
E_2(X) = product_(m>=1, m odd) (1-X^m)^(-mu(m)/m).
~~~

Taking logarithms proves the identity coefficientwise by Mobius inversion.
Each exponent is in Z_2, whose binomial coefficients are integral: approximate
it by nonnegative integers at arbitrary 2-adic precision and use the fixed
finite denominator in the binomial polynomial. Thus a_n lies in Z_(2).

The existing [exact acquisition](../tools/witt_artin_hasse.py) is reused:

~~~
a_0=1,   n*a_n = sum_(2^j<=n) a_(n-2^j).
~~~

Only n<P is needed. Exact rational arithmetic and odd-denominator reduction
are polynomial in P. For completeness, the formal exponential shows that
each denominator divides n!, while positivity and evaluation at X=1/2 bound
a_n by a constant times 2^n. Numerators and denominators therefore have
O(P log(P+2)) bits. There is no enumeration of a field of size 2^w.

There is a unique pi in 2Z_2 with E_2(pi)=-1. Indeed E_2(X)+1 has constant
term two, linear coefficient one, and derivative odd at every even X;
Hensel lifting starts at zero modulo two. Reduction modulo four gives
pi=2 modulo four, so v_2(pi)=1. The implementation lifts this root one
bit at a time using the acquired polynomial. Define

~~~
theta(X) = E_2(pi X) = sum_(n>=0) a_n*pi^n X^n.
~~~

The nth coefficient is divisible by 2^n. Modulo 2^P, terms n>=P vanish
on integral inputs; acquire the explicit polynomial of degree at most P-1.
All precision and acquisition costs are included in (1).

This is standard Dwork splitting, specialized to characteristic two.
For primary-source context see Lin and Wan,
[*On inverted Kloosterman sums over finite fields*, Section 2.3.1](https://arxiv.org/pdf/2301.04287),
which defines E_p, pi, theta and the monomialwise Teichmuller product.
The following elementary argument fixes the exact binary normalization used
here; no Fredholm determinant or cohomological algorithm is assumed.

## 2. The binary character identity

Work first in the complete unramified integer ring O containing F_(2^w).
For a Teichmuller lift x, let sigma(x)=x^2. Then

~~~
N(x) = product_(s=0)^(w-1) theta(sigma^s(x))
     = (-1)^(Tr_(F_(2^w)/F_2)(x mod 2)).                 (3)
~~~

Here the product is fixed by sigma and belongs to 1+2Z_2. The convergent
2-adic logarithm gives

~~~
log N(x)
 = sum_(j>=0) (pi^(2^j)/2^j) sum_(s=0)^(w-1) x^(2^(s+j))
 = (sum_(s=0)^(w-1) x^(2^s)) * log E_2(pi)
 = 0.
~~~

The inner sum is unchanged by a cyclic shift because x^(2^w)=x.
All rearrangements converge. The kernel of log on 1+2Z_2 is {1,-1}:
squaring puts a kernel element into 1+8Z_2, and for a nonzero z in 4Z_2
the leading term z of log(1+z) has strictly smaller valuation than every
later term, so log is nonzero there. Hence N(x) is a sign. Modulo four,

~~~
N(x) = 1 + pi*sum_s sigma^s(x)
     = 1 + 2*Tr_(F_(2^w)/F_2)(x mod 2)    mod 4,
~~~

which distinguishes the two signs and proves (3). The identity includes
x=0. It also applies when w is a proper divisor of the coefficient degree
d: the Teichmuller subfield elements are fixed by sigma^w, so the same
product is fixed by sigma. Reduction modulo 2^P preserves the identity.

## 3. Split monomials, not arbitrary polynomial lifts

Write f=sum_mu c_mu U^mu in the residue field. Replace every coefficient
by the Teichmuller lift of its residue. Each monomial value

~~~
z_mu = [c_mu] * product_a U_a^(mu_a)
~~~

is Teichmuller. Apply (3) to each z_mu and multiply. Binary trace is
additive, so this computes the sign character of Tr f. Higher digits of
the original coefficients have no effect on the predicate.

Applying theta to the whole polynomial value is generally **wrong**.
For example, in GR(32,2), two constant Teichmuller monomials 1+1 have
residue trace zero and hence character one. With pi=10 modulo 32,
theta(2)^2=25 modulo 32. The monomialwise product gives one. This
counterexample is a mandatory exact checker control.

For a monomial mu and shift s=rho+t*j, acquire the synchronized factor

~~~
sum_(n=0)^(P-1) a_n*pi^n * sigma^s([c_mu])^n
                     * product_a Y_(a,j)^(n*mu_a).
~~~

It has at most P terms and degree at most (P-1)*|mu|. There are w such
factors per nonzero residue monomial, acquired with polynomially many
represented ring operations. Their product is the desired sign character.
The total formal degree over all predicates is at most (P-1)wD, without
extra summation variables or relative columns.

## 4. Fourier projection and exact precision loss

Let chi_i(U) be these sign characters. Pointwise in O,

~~~
sum_(lambda in {0,1}^r) (-1)^(sum_i lambda_i*b_i)
                       * product_i chi_i(U)^(lambda_i)
 = product_i (1+(-1)^b_i chi_i(U))
 = 2^r * product_i [Tr f_i(U)=b_i].                     (4)
~~~

For each lambda the canonical rational compiler evaluates the corresponding
scalar synchronized product, retaining the supplied scalar numerator and
total unit inverse denominator factors. All operations occur at precision
P=H+r. Its at most P+1 inverse queries and degree contribution
`((2P+1)d-1)S_B` are charged as before.

Sum these 2^r outputs with their indicated signs. By (4), the numerator
belongs to 2^r R_P even with arbitrary scalar weights. The unramified
polynomial basis is free over Z/(2^P), so every canonical coefficient is
divisible by 2^r as an ordinary integer. Divide those coefficients and
reduce modulo 2^H in the **same defining polynomial**. This quotient is
well defined modulo 2^H, independently of the numerator representatives.
No nonunit ring inverse is used. The interface requires the extra r input
precision bits; it does not claim a precision-free projection.

Equation (4), the exact character identity and the existing carry theorem
prove correctness and the conservative cost (1). Empty predicate lists,
contradictions and punctured domains all satisfy the same pointwise identity.

## 5. Full-row carry bound

For t=w there are v formal variables per row. Every character monomial
contributes one conjugate factor per row, so its row degree is at most
(P-1)|mu|. For any one variable, the sum of its maximal exponents over
all these factors, weighted by the row powers of two, is at most

~~~
(P-1)*mu_a*sum_(s=0)^(w-1) 2^s
 = (P-1)*mu_a*(2^w-1).
~~~

The unique-cycle carry bound divides this number by 2^w-1. Thus the
character contribution to every initial carry is at most (P-1)D.
This is a structural orbit bound, not a fitted finite maximum.

A denominator's d conjugates visit every row d/w times. The norm power
is at most 2P, and its adjugate is bounded by another full orbit. Their
row degrees and initial carries are therefore at most
`(2P+1)(d/w)S_B`. Numerators contribute at most S_A. Consequently every
row degree and initial carry is at most K from (2). Explicit row tables,
carry states and ring operations cost `(K+2)^O(v)` times polynomial
overhead, proving (2). Constant factors have degree zero but their written
length and acquisition cost remain in the polynomial input factor.

## 6. Implementation, controls and next question

The [compiler](../tools/artin_hasse_trace_carries.py) acquires the scalar
splitting polynomial, Teichmuller-lifts monomial coefficients, makes the
2^r canonical rational queries and checks coefficientwise divisibility
before returning a precision-H ring element. It never enumerates its
summation field. The supplied input weights are represented at precision P;
the output ring uses precision H and the same quotient polynomial.

The [checker](../tools/check_artin_hasse_trace_carries.py) and
[receipt](../evidence/artin-hasse-trace-controls.json) retain 120 full-field
character cases, four proper-subfield cases, the non-Teichmuller sum
counterexample, 14 contractions against 80 direct assignments, 116
pointwise polynomial characters and 14 comparisons with the independent
divided-power coefficient-algebra implementation. Full-row degree/carry
caps are checked directly from the acquired factors. Five invalid-input
controls cover insufficient precision, coefficient subfield, negative
exponents, target bits and incompatible dimensions. The panel includes
nonlinear/coupled traces, even-width constants, higher coefficient digits,
relative columns, punctures, scalar poles/nonunits and a nonscalar weighted
quotient. Run make verify-artin-hasse-traces.

For H,K polynomial in an input bit bound b, (2)'s exponential part is
`2^O(r+v log b)`. A complete arithmetic representation with
`r+v log b=O(b/log b)`, with all other costs at that scale, would meet a
score-one budget. This is a sufficient parameter criterion, **not** such
a representation. Native Witt digit extraction, arbitrary exact-division
circuits and high-conductor characters are outside this binary interface.

The splitting function, Artin--Hasse integrality and Fourier orthogonality
are classical. Existing ring carries, scalar inverses, trace packets and
Witt coefficient acquisition retain their source credits. This application
and its ordinary argument were developed with OpenAI assistant support;
finite controls are not independent mathematical acceptance. The next
research step must obtain or evaluate actual integer arithmetic within
the charged budget, rather than enlarge this generic character panel.
