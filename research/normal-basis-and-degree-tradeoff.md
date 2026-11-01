# Exact degree/column tradeoff for a direct packed Boolean AND

## Research question and result

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19337929)
tests an actual word gate: pack Boolean coordinates in a normal basis of
F_(2^w), then try to use low-degree synchronized polynomials for bitwise
AND. Normal bases make bit rotation an inexpensive Frobenius operation;
field multiplication is different from coordinatewise Boolean AND.

Fix **any** normal basis, w>=1, and t dividing w. Allow a direct
representation as a **sum** of synchronized row polynomials in X,Y, with
C=w/t relative columns and no auxiliary variables. The exact minimum
maximum formal row degree is

~~~
delta_min(t) = 2^floor(t/2)+1.                           (1)
~~~

Thus polynomial degree permits only O(log delta) rows. If a proposed
compiler is charged by the present generic parameter

~~~
(v*w/t) * log_2(S+2),      S>=delta_min(t),
~~~

that parameter is at least v*w/3. Trading degree for fewer relative columns
does not make this parameter sublinear in the packed bit count. This is
a restriction on this direct normal-basis gate representation and the
specified estimate. It is **not** a time lower bound for every algorithm,
a factoring lower bound, or a theorem about arbitrary bases, auxiliaries,
predicate-only encodings or special structured contractions.

The representation class even allows adding row polynomials. The current
carry API accepts synchronized **products**; no automatic conversion of a
row sum into that API is supplied. This permissive gate test fails the
degree/dimension criterion before that additional compilation problem.
The complete factoring score remains **1/2**.

The later [opposite-direction diagnostic](mixed-direction-frobenius-carries.md)
changes the aligned-row assumption used here. Its fixed F16 panel admits
degree four after reversing one word's direction and retiming its input;
the present aligned minimum remains five. No asymptotic compression or
complete arithmetic compiler follows from that finite improvement.

## 1. The packed gate and its unique bilinear polynomial

Let beta_j=beta^(2^j), for 0<=j<w, be a normal basis, and let
gamma_j=gamma^(2^j) be its trace-dual normal basis. Write

~~~
X=sum_j x_j beta_j,    Y=sum_j y_j beta_j,
H(X,Y)=sum_j x_j y_j beta_j,       x_j,y_j in F_2.
~~~

The dual coordinates are x_j=Tr(gamma_j X). Expanding both traces gives
the bilinear linearized representation

~~~
H(X,Y)=sum_(i,j=0)^(w-1) c_ij X^(2^i) Y^(2^j),
c_ij=Tr(beta * gamma^(2^i) * gamma^(2^j)) in F_2.        (2)
~~~

Indeed, the coefficient obtained directly from the coordinate definition
is the sum over k of
`beta^(2^k) gamma^(2^(k+i)) gamma^(2^(k+j))`, exactly the displayed trace.
The polynomial is unique as a function on F_(2^w)^2: each variable has
degree at most 2^w-1, and iterated univariate root bounds give uniqueness
in that full function basis. This also handles w=1.

**Every cyclic offset occurs.** For any k modulo w, not all
c_(i,i+k) vanish. Put a=gamma*gamma^(2^k), which is nonzero. Then

~~~
c_(i,i+k) = Tr(beta * a^(2^i))
          = Tr(beta^(2^(w-i)) * a).
~~~

As i varies, beta's conjugates span the entire field. If all these traces
vanished, nondegeneracy of the trace pairing would force a=0, a
contradiction. Consequently (2) contains at least one monomial in every
offset class i-j modulo w, and hence in every offset class modulo t.
This argument is uniform over all normal generators, not a finite-field
pattern extrapolation.

## 2. Subthreshold synchronized monomials miss an offset

A row rho, 0<=rho<t, has formal symbols

~~~
X_j=X^(2^(rho+t*j)),    Y_j=Y^(2^(rho+t*j)),    0<=j<C.
~~~

A formal monomial has coefficient in the field and counts alpha_j,beta_j,
with total degree sum(alpha_j+beta_j). Its actual exponents are

~~~
A=2^rho * sum_j alpha_j*2^(t*j),
B=2^rho * sum_j beta_j *2^(t*j).                         (3)
~~~

On the **full** field, exponent zero stays zero. A positive exponent e
reduces to `1+((e-1) mod (2^w-1))`, in {1,...,2^w-1}; it is not reduced
to exponent zero. This distinction retains the correct value at X=0.
Multiplying by 2^rho cyclically rotates a positive exponent's w binary
digits, with the all-ones exponent fixed.

Put h=floor(t/2), and suppose the total degree is below 2^h+1, hence at
most 2^h. Every individual count is below 2^t. The unshifted integer in
(3) therefore has separate t-bit blocks, with no overlap between blocks,
and lies in [0,2^w-1]. To reduce to a bilinear exponent 2^i it must
already be a power of two: precisely one alpha_j is nonzero, of value
2^a, 0<=a<t. Similarly the Y counts have a single value 2^b.

Thus a bilinear term at offset i-j modulo t must have

~~~
i-j = a-b mod t,          formal degree = 2^a+2^b.
~~~

The offset h modulo t needs at least one of a,b to be h or larger.
If both are below h, their difference cannot equal h modulo t; for
t=2h or t=2h+1 the alternative representative is at most -h.
So its degree is at least 2^h+1, contradicting the assumed threshold.
For t=1 the same conclusion says a bilinear term needs degree at least
two. Hence **every** subthreshold row monomial misses offset h.

Uniqueness of the polynomial function basis means cancellations among
other terms cannot create a missing bilinear monomial. Summing arbitrarily
many row polynomials of that degree still misses the offset, whereas (2)
requires it. This proves the lower bound in (1), including positive
exponent wraparound and zero field inputs.

## 3. A matching explicit row representation

For each nonzero c_ij in (2), choose a,b with one equal to zero and the
other at most floor(t/2), such that

~~~
a-b = i-j mod t.
~~~

Choose rho=i-a=j-b modulo t. Appropriate relative-column indices,
reduced modulo C, represent X^(2^i)Y^(2^j) by a row-rho monomial with
counts 2^a and 2^b. Negative column quotients simply wrap a positive
exponent once around its Teichmuller period; exponent zero is not used.
The formal degree is at most 2^floor(t/2)+1. Applying this to all w^2
possible bilinear terms proves the upper bound in (1).

The [acquisition module](../tools/normal_basis_and_tradeoff.py) constructs
the dual basis by binary elimination on the trace-pairing matrix, computes
(2), and emits these row polynomials. Field and normal generators are
supplied. No domain enumeration occurs in acquisition. Numeric degrees,
exponent sizes and row tables are explicit; in particular the large
degree at t=w is not hidden by writing a binary exponent.

## 4. Consequence for this parameter route

Let M=v*w be the number of bits in v such words. For t>=2,
floor(t/2)>=t/3, so S>=delta_min(t) gives

~~~
(v*w/t)*log_2(S+2) >= M/3.
~~~

For t=1 the inequality also holds, since S>=2. Thus the sufficient
sublinear-parameter condition for the generic degree/dimension bound
fails at every row choice. This is an obstruction to obtaining a saving
from **that bound alone**. An upper-bound parameter being large is not
a proof that every implementation needs its full state box or that no
additional gate structure can be exploited.

The result concerns a direct full-output AND polynomial. A factoring
representation that avoids such gates, uses a different basis/algebra,
introduces charged auxiliaries, or only evaluates a carefully structured
global predicate is outside its scope. In particular it does not prove
the native [Witt-trace division family](witt-trace-native-word-circuits.md)
hard to contract. It explains why one tempting binary-linear word packing
does not immediately supply the missing arithmetic input.

## 5. Exact fixed controls and stopping condition

The [checker](../tools/check_normal_basis_and_tradeoff.py) and
[receipt](../evidence/normal-basis-and-tradeoff.json) exhaust all 13 normal
bases in F_4,F_8,F_16. They check 2,272 word pairs against literal bitwise
AND for (2), and 6,592 word-pair evaluations of the acquired row forms
over every allowed divisor t. An independent enumeration of 141
subthreshold row monomials verifies the missing-offset control. Dual
pairing matrices, all-offset witnesses, row tables and exact degrees
remain in the receipt. Run make verify-normal-basis-and.

The proof establishes all widths/bases; the finite panel checks its
implementation and wrap conventions. No larger field or degree sweep
is needed. The registered continuation condition would have required a
counterexample giving cheap actual AND at t much larger than log delta;
the exact tradeoff closes that direct interface instead. New work should
change the arithmetic representation or the global evaluator.

Normal bases, trace duality, bilinear linearized polynomials and uniqueness
of finite-field polynomial functions are standard. The offset and
synchronization argument here was developed with OpenAI assistant support,
using the existing repository field implementation. This is an ordinary
author-checked argument with finite controls; no new Lean or independent
mathematical review is claimed.
