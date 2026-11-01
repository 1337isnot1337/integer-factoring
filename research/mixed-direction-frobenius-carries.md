# Polynomial carry contraction with opposite Frobenius directions

## Question and result

The [registered mechanism](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19338422)
changes the synchronization assumption: selected word variables move
backward through Frobenius rows while the others move forward. Both carry
cycles still lie on the same underlying cycle, so their orientations can
be processed together.

Let R=GR(2^H,d), w dividing d, and let each of v variables range over the
Teichmuller copy of F_(2^w), optionally punctured. Supply a direction
epsilon_i in {1,-1} for each variable. At row rho in {0,...,w-1}, a factor
is an explicit ordinary polynomial of degree delta_f in

~~~
U_i^(2^(epsilon_i*rho mod w)),       1<=i<=v.
~~~

For S=sum_f delta_f, the [implemented polynomial extension](../tools/mixed_direction_carries.py)
computes the product sum with

~~~
(S+2)^O(v) * poly(input_length,H,d,w,v)                 (1)
~~~

bit cost, using one relative column per word. The same scalar Teichmuller
moments and actual ring arithmetic are retained. Explicit input acquisition,
numerical degree and the coefficient model remain charged. The extension
has **polynomial factors only**; it does not inherit denominator,
trace-projector or arbitrary circuit closure for free.

The registered actual-gate diagnostic gives a smaller representation of
normal-basis bitwise AND over F_16: maximum row degree four, compared with
the aligned-direction minimum five. This is an exact finite improvement
in that representation parameter, not an asymptotic family or a complete
counting/factoring algorithm. The complete factoring score remains **1/2**.

## 1. Expand only the supplied local row polynomials

For each row multiply its explicit factors into one sparse table. A global
monomial choice has nonnegative counts n_(i,rho), with total at most S.
The exponent of variable i is

~~~
E_i=sum_rho n_(i,rho)*2^(epsilon_i*rho mod w).
~~~

Let q=2^w. The scalar moment is q when E_i=0 on a full domain, q-1 when
E_i>0 and q-1 divides E_i, and zero otherwise. On a punctured domain it
is q-1 also at E_i=0. These are the actual Teichmuller sums in R,
not merely their reductions modulo two.

Let B_i be the sum, over factors, of the maximum written count of variable
i in that factor. Then B_i<=S and every monomial choice has
`sum_rho n_(i,rho)<=B_i`. The implementation charges the whole box
`product_i {0,...,B_i}`; it does not infer smaller states from cancellation.

## 2. A forward variable has the usual carry cycle

For epsilon=1 use states a_rho with

~~~
a_(rho+1)=(a_rho+n_rho)/2,       a_w=a_0.               (2)
~~~

Iteration gives `sum_rho n_rho 2^rho=(q-1)a_0`. Thus a closed integer
cycle exists exactly when the exponent is divisible by q-1. It is unique,
and its rotated formulas are

~~~
a_rho = sum_(j=0)^(w-1) n_(rho+j mod w)*2^j /(q-1).
~~~

All these numbers are nonnegative integers at a divisible exponent and
at most B. The same formula shows a_0=0 exactly when all n_rho=0.
At each scan step, parity is checked before division in (2).

## 3. Reverse the carry edge, not the scan

For epsilon=-1 scan the same rows in increasing rho, but use

~~~
a_(rho+1)=2a_rho-n_rho,          a_w=a_0.               (3)
~~~

Closing this recurrence gives

~~~
(q-1)a_0=sum_(rho=0)^(w-1) n_rho*2^(w-1-rho).
~~~

The right side is congruent to `2^(w-1) E` modulo q-1. Multiplication by
that power of two is invertible modulo the odd number q-1, so it has the
same divisibility condition as the desired reverse exponent E.
The unique rotated state formula is

~~~
a_rho=sum_(j=0)^(w-1) n_(rho+j mod w)*2^(w-1-j)/(q-1).
~~~

It is again integral, nonnegative and at most B. Its zero initial state
again occurs exactly for the all-zero exponent. The scan rejects a
negative target or one outside the charged bound, and never discards a
valid divisible monomial. Width w=1 gives a=n in either orientation,
including its zero-exponent moment.

Combining (2) or (3) independently for each variable proves that every
surviving monomial has exactly one closed joint carry cycle, with the
correct scalar moment product. Summing its row coefficients along those
cycles therefore gives exactly the original polynomial sum.

## 4. Cost and coefficient precision

Each local row table has at most `(S+1)^v` exponent positions. There are
at most `(S+1)^v` initial carry states and that many states at each scan
boundary. Even scanning every local monomial at every state gives
`(S+2)^O(v)` represented ring operations, times a polynomial in w and
the written input. The explicit ring operations have polynomial bit cost
in H,d. This proves (1) without enumerating the field domain.

For W=1 the result is its stated full-domain count modulo 2^H. For general
polynomial weights it is a weighted R-valued sum, not automatically an
existence test. Growing precision does not turn an arbitrary weight into
a Boolean predicate. Field models, embeddings and their acquisition are
still supplied or charged separately.

## 5. An actual AND diagnostic outside the aligned-row theorem

The [aligned normal-basis tradeoff](normal-basis-and-degree-tradeoff.md)
assumes both variables use the same Frobenius direction. Its coefficient
polynomial for AND has nonzero c_ij at every difference offset i-j.
That proof does not constrain opposite directions in the same way.

For forward X and reverse Y, a row monomial with counts 2^a,2^b has
exponents i=rho+a and j=-rho+b modulo w, so its **sum** offset is
i+j=a+b modulo w. An independent fixed Frobenius retiming of the Y input
shifts this offset set cyclically.

The diagnostic acquires the exact anti-diagonal support from the same
13 normal bases in F_4,F_8,F_16 as the earlier control. For a fixed
retiming, choose the least nonnegative residue s of i+j minus the retime,
then a=floor(s/2), b=ceil(s/2). This explicitly represents the monomial
with degree `2^a+2^b`. It changes the representation; it does not permute
summands or assume a new tensor identity.

All eight F_16 normal bases in this fixed panel admit maximum degree
four after retiming, below the aligned minimum five. Every word pair
agrees with literal bitwise AND. F_4 and F_8 attain degree three, equal
to their aligned bound. No asymptotic anti-diagonal support theorem is
inferred from these finite cases, and no larger basis/field scan is part
of this checkpoint.

The acquired row polynomials for AND are **added**. The contraction
algorithm multiplies its supplied factors. No automatic row-sum-to-product
conversion, composition of word gates or integer multiplier query is
claimed. Such operations require their own complete invoice.

## 6. Inverses and trace projection require a new check

Applying scalar sigma^k to a mixed factor advances the exponent indices
of **every** input variable by k. Relabeling its row would need a forward
variable to advance by k and a reverse variable to retreat by k. These
are incompatible in general. Thus the canonical norm/adjugate inverse
expansion does not preserve this one-column mixed alignment for free.

The same issue applies to the Frobenius orbit in a binary trace splitting
function. Extra columns, variables, graph states or another identity must
be charged before either extension is claimed. The new module exposes no
denominator or trace-constraint API, so it cannot accidentally invoke the
aligned compiler with an incorrect dimension estimate.

The [checker](../tools/check_mixed_direction_carries.py) and
[receipt](../evidence/mixed-direction-carry-controls.json) retain 11
contractions against 732 direct assignments, two uniform-direction
comparisons with the canonical engine, four invalid-input controls and
2,272 AND word-pair evaluations across the 13 fixed normal bases. Cases
include an output that differs from the all-forward control, a full
integer domain count at sufficient precision, zero/punctured inputs,
both width parities, a proper summation subfield and non-Teichmuller
weight coefficients. Run make verify-mixed-direction-carries.

The ring model, Teichmuller moments, cyclic-carry method and normal-basis
AND acquisition retain the preceding contributors' credit. This changed
orientation and its ordinary proof were developed with OpenAI assistant
support; there is no new Lean, independent mathematical acceptance or
stronger factoring bound. The next step needs a concrete arithmetic scalar
or a structural anti-diagonal family, with any inverse/projector costs
included, rather than a larger generic test panel.
