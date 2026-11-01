# Exact trace-constrained counts at growing precision

## Question and result

Base: 80835e8. The [registered scope](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19337426)
asks whether the earlier [square-zero construction](square-zero-trace-contraction.md)
can retain integer counts at growing precision without replacing every
word by all its binary carry columns. Divided powers supply that extension.
This changes the coefficient ring and output guarantee; it is not another
kernel search or a larger numerical screen.

Use an explicitly represented R = GR(2^H,d), with w dividing d, positive
t dividing w, C = w/t, and v summation variables on the Teichmuller copy
of F_(2^w), or its nonzero elements. Impose r predicates

~~~
Tr_(F_(2^w)/F_2)(f_i(U)) = b_i,       b_i in {0,1}.        (1)
~~~

Each f_i is an explicit ordinary polynomial of degree at most delta_i.
Its coefficient residues must lie in the summation subfield; higher
coefficient digits may lie outside the lifted subring. Let W be a scalar
synchronized product with numerator/inverse degree totals S_A,S_B, using
the [total unit inverse](precision-efficient-residue-indicators.md).
Put

~~~
J = 2H-1,
B = S_A + ((2H+1)d-1)S_B + J*t*sum_i delta_i.
~~~

The implemented compiler computes the constrained R-valued sum with cost

~~~
(2H)^O(r) * (B+2)^O(v*w/t)
          * poly(input_length,H,d,w,v,r).                (2)
~~~

Coefficient dimension, row-polynomial acquisition, ring arithmetic and
precision are charged. Field models and coefficient embeddings must be
supplied or acquired separately. For affine predicates, binary row reduction
makes r the actual constraint rank. For nonlinear predicates r is the list
length unless an equivalent reduction has been proved.

With W=1, the result is the solution count modulo 2^H. If H>v*w, its
canonical residue is the **exact integer count**, since the domain has at
most 2^(v*w)<2^H points. Smaller precision can still cancel. An arbitrary
weighted R-valued sum is not automatically an integer count or existence test.
The complete integer multiplication representation is not supplied, and
the factoring score remains **1/2**.

The later [Artin--Hasse projection](artin-hasse-trace-projection.md) handles
this polynomial trace input class through 2^r scalar Fourier queries at
precision H+r. It charges r division bits and uses classical monomialwise
splitting, giving an alternative to the growing coefficient algebra here.
Neither construction evaluates native integer digit/division circuits.

## 1. The projector needs only a Boolean residue

Retain the standard integer Hermite polynomial from the
[indicator companion](precision-efficient-residue-indicators.md):

~~~
E_H(X) = sum_(j=H)^(2H-1) binom(2H-1,j) X^j(1-X)^(2H-1-j).
~~~

The exact divisibilities X^H | E_H and (1-X)^H | 1-E_H imply that, in any
commutative ring killed by 2^H, E_H(y) equals the exact Boolean residue of y
whenever that residue is 0 or 1. Indeed y or 1-y belongs to 2R, so its Hth
power vanishes. No inverse or division by an extension degree is involved.

For predicate i use

~~~
T_i(U) = sum_(s=0)^(w-1) sigma^s(f_i(U)) - b_i,
I_i(U) = 1-E_H(T_i(U)).                                  (3)
~~~

The residue of T_i is the binary field trace minus b_i. Hence I_i is exactly
the 0/1 indicator of (1). T_i itself need not lie in Z/(2^H) when higher
coefficient digits lie outside the subring; the divisibility proof still
applies. Subtract the target bit **once, outside the trace**. Moving it
inside loses it when w is even. Apply the actual lifted automorphism to
coefficients; ordinary squaring is generally different.

## 2. A finite divided-power algebra

Let A be the free R-module with basis gamma_e for e in {0,...,J}^r, with

~~~
gamma_e gamma_f = product_i binom(e_i+f_i,e_i) gamma_(e+f)
                   if all e_i+f_i<=J,
                 0 otherwise.                           (4)
~~~

Integer binomial coefficients are mapped into R. This is an associative
commutative algebra with unit gamma_0. In one coordinate, associativity is

~~~
binom(a+b,a)binom(a+b+c,a+b)
  = binom(b+c,b)binom(a+b+c,a).
~~~

Coordinates multiply independently. If a total coordinate degree exceeds J,
both associated products vanish, so truncation preserves the identity.
The dimension is M=(J+1)^r=(2H)^r. Even naive all-pairs operations cost
M^O(1) times a polynomial in H,d,r. This growing cost is not free.
The algebra has nilpotents; where permitted, gamma_1^2=2 gamma_2,
which differs from the square-zero rule above precision one.

For scalar polynomials define

~~~
Exp_i(A) = sum_(k=0)^J A^k gamma_(k e_i).
~~~

The binomial theorem gives the integral identity

~~~
Exp_i(A) Exp_i(B) = Exp_i(A+B).                           (5)
~~~

This standard divided-power identity uses no factorial denominators.

## 3. Acquire the row packets

Write f_i=sum_mu c_mu U^mu. At row rho, in the usual vC formal symbols,
form

~~~
Q_(i,rho) = sum_(j=0)^(C-1) sum_mu
   sigma^(rho+t*j)(c_mu) product_a Y_(a,j)^(mu_a)
   - [rho=0]*b_i.
~~~

Each packet has formal degree at most delta_i. Its Teichmuller evaluation,
summed over the t rows, is T_i. Acquire the r*t explicit polynomial
factors Exp_i(Q_(i,rho)) over A, with total formal degree at most
J*t*sum_i delta_i and no additional relative column.

This acquisition is charged: multiply each packet table successively to
get its powers through J. A power has at most (J*delta_i+1)^(vC) coefficient
positions. These computations fit the polynomial-power state cost in (2);
constants have polynomial acquisition cost. Numerical degrees are parameters,
not free binary exponents, and no arithmetic-circuit closure is assumed.

By (5), their product equals product_i Exp_i(T_i). Put
ell_k equal to the coefficient of X^k in 1-E_H(X), and define the **R-linear** functional

~~~
Phi(gamma_e) = product_i ell_(e_i).
~~~

Then

~~~
Phi(product_i Exp_i(T_i)) = product_i I_i.               (6)
~~~

Phi is generally not multiplicative. Apply it after the whole contraction;
applying it to separate rows would discard the correlations needed by (6).

## 4. Why the carry and inverse compilers still apply

The [ring carry proof](teichmuller-cyclic-carries.md) uses addition,
multiplication and central scalar Teichmuller moments. Its unique integer
carry cycle is independent of the coefficient ring. Embed R and its
summation points into A. The actual weights q at zero exponent and q-1 at
positive divisible exponents remain central scalars. Thus the polynomial
calculation extends to A without treating its nilpotents as field elements.

Every supplied denominator remains R-valued. Its norm, adjugate and
collected inverse polynomial are valid in R first, then remain valid after
embedding into A. Their degree contribution is ((2H+1)d-1)S_B.
**No general A-valued denominator is inverted.**

The total degree is bounded by B. The scalar carry engine needs
(B+2)^O(vC) represented coefficient operations, with at most H+1 collected
inverse queries. Multiplying by A's operation cost, explicit packet
acquisition and the final linear functional proves (2). Linearity moves
Phi through the outer sum, so (3) and (6) prove the output guarantee.

## 5. Implementation and exact controls

The [implementation](../tools/divided_power_trace_carries.py) provides
trace_constrained_sum for polynomial predicates and scalar weights.
The coefficient representation has integer zero as its zero element and
retains modular ring addition, rather than XOR. Its coefficient automorphism
fixes the formal gamma basis; it is not ordinary squaring on A.

The affine_trace_sum frontend reuses
trace_constraint_contraction.independent_constraints. Nondegeneracy of the
residue trace pairing makes XOR elimination on the coefficient vectors
equivalent to the original affine system. Reduced vectors are lifted again;
only the residue predicates matter. Inconsistency returns zero immediately,
and punctured domains cause no change to this pointwise equivalence.
No corresponding rank reduction is inferred for nonlinear lists.

The [checker](../tools/check_divided_power_trace_carries.py) retains the
predeclared panel: 4,160 basis associativity triples, 256 exponential
identities, 13 contractions against 76 direct assignments, and six affine
controls. It independently computes field traces using binary-field
arithmetic and checks each lifted indicator pointwise before comparing sums.
Precision-one affine cases agree with the earlier square-zero evaluator.

The panel includes coupled nonlinear predicates, relative columns, mixed
punctured domains, higher coefficient digits outside the subring, an
even-width constant target, a cubic trace, exact integer recovery, a selected
scalar pole and a constant nonunit denominator. The
[receipt](../evidence/divided-power-trace-controls.json) preserves inputs,
direct points, outputs, algebra dimensions, degree invoices and source hashes.
Run make verify-divided-power-traces.

## 6. Relevance and the remaining arithmetic problem

At t=w, with H and B polynomial in an input bit bound b, the exponential
part is 2^O((v+r)log b). Thus a uniform arithmetic representation with
v+r=O(b/(log b)^2), suitable precision, and fully charged construction and
recovery would meet a score-one time budget. This is a parameter criterion,
**not** such a factoring representation. A naive bit-constraint description
can have r proportional to b and does not meet it.

The earlier operand-prefix interface has a small affine trace boundary,
but returned field-valued sums. This extension removes that particular
precision/cancellation obstruction for represented counts. It still does
not represent the complete integer phases or prove the global constructor
and recovery costs. A zero residue at insufficient precision is not a
no-solution conclusion.

The next research step must acquire an actual arithmetic representation
with the needed variable/trace budget. Larger generic precision or rank
panels are not a stronger factoring result. Existing square-zero packets,
cyclic carries, ring lifting, interpolation and scalar inverse compilation
retain their credits. These are ordinary author-checked arguments and
finite controls with substantial OpenAI assistant support; no new Lean or
independent mathematical review is claimed.
