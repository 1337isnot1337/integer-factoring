# Linear-precision polynomials for residue roots and unit inverses

## Research question and result

Base: `ab5f29b`, with concurrent PR71 head `b9cd96f`. The Teichmuller/cyclic-carry work lifts polynomial
sums to a Galois ring. A lifted field norm is not itself an exact zero
indicator. Using the full exponent of the unit group to repair this
would have exponential degree in the precision. Standard integer Hermite
interpolation instead gives **linear degree in the precision**.

The concurrent PR71 count backend uses the same standard Hermite indicator.
This companion develops the total unit-inverse and collected-denominator
interface, and coefficient-extension packing of equations. The indicator
and polynomial carry framework retain that preceding/concurrent credit.
It supplies ordinary identities and small
explicit ring controls, with no complete factoring graph or score gain.

## 1. An integer parity polynomial

For L>=1 define the integer polynomial

```
E_L(X) = sum_(j=L)^(2L-1) binom(2L-1,j) X^j (1-X)^(2L-1-j).
```

Its degree is `2L-1`. Its two exact divisibility identities are

```
X^L divides E_L(X),
(1-X)^L divides 1-E_L(X).                               (1)
```

The first is visible termwise. The binomial theorem and the substitution
j -> 2L-1-j give `1-E_L(X)=E_L(1-X)`, proving the second. Consequently,
for every integer x,

```
E_L(x) = 0 mod 2^L if x is even,
E_L(x) = 1 mod 2^L if x is odd.                         (2)
```

The expanded coefficients are zero below L and otherwise

```
[X^k] E_L = (-1)^(k-L) binom(2L-1,k) binom(k-1,L-1),
              L<=k<=2L-1.                              (3)
```

For example `E_2=3X^2-2X^3` and `E_3=10X^3-15X^4+6X^5`.
Formula (3) follows by expanding the displayed binomial tail and using
`sum_(j=L)^k (-1)^(k-j) binom(k,j)=(-1)^(k-L)binom(k-1,L-1)`.
There are L nonzero coefficients, each with O(L) bits, computable in
polynomial time in L. This familiar interpolation construction is not
claimed as a new polynomial approximation principle.

## 2. Exact residue-root indicators over an unramified ring

Let `R=GR(2^L,w)` be explicitly represented, with Frobenius automorphism
sigma and fixed ring `Z/2^L`. For a in R put

```
n(a) = product_(h=0)^(w-1) sigma^h(a),
a#   = product_(h=1)^(w-1) sigma^h(a).
```

The empty product at w=1 is one. The norm n lies in the fixed ring,
`a*a#=n(a)`, and its reduction modulo two is the binary field norm.
It is even exactly when the residue of a is zero, equivalently when
a is a nonunit. Equations (1)-(2) therefore give the exact indicator

```
delta_L(a) = 1-E_L(n(a))
          = 1 if a mod 2 = 0, and 0 otherwise, in R.      (4)
```

This tests a **residue-field root**. An even nonzero a is deliberately
accepted. It does not test whether a is zero at full ring precision.

Let F be a polynomial over R, evaluated on the q=2^w Teichmuller lifts
in each of v variables. The sum of (4) is the number of roots of the
residue polynomial modulo `2^L`. If `L>v*w`, the root count is at most
`2^(v*w)<2^L`; its canonical residue recovers the exact nonnegative
integer count. For several polynomial equations, the product of their
indicators gives the simultaneous root count. Costs of that product
must still be charged as described below.

On Teichmuller inputs U, the exact Frobenius identity is

```
sigma^h(F(U)) = F^(sigma^h)(U_1^(2^h),...,U_v^(2^h)).    (5)
```

It uses a ring automorphism, including its action on coefficients.
Replacing it by `F(U)^(2^h)` is generally false at precision above one.

## 3. A total unit inverse with linear precision degree

Define the integer polynomial

```
J_L(X) = E_(L+1)(X)/X.                                 (6)
```

This is **formal polynomial division** before evaluation. It is allowed
because `X^(L+1)` divides the numerator. It has degree 2L, and
`X^L` divides J_L. If x is even, `J_L(x)=0 mod 2^L`. If x is odd,
`x J_L(x)=E_(L+1)(x)=1 mod 2^(L+1)`, so J_L is its inverse modulo
`2^L`. No division by a nonunit and no extra working precision is needed.
Thus

```
inv*_L(a) = a# J_L(n(a))
          = a^(-1) if a is a unit, and 0 otherwise.      (7)
```

For units, multiply (7) by a and use `a*a#=n(a)`; for nonunits the
J factor vanishes. This is a total **unit inverse**, with the nonunits
assigned zero. It is not an inverse for nonzero zero divisors. Its
reduction agrees with the field's zero-preserving inverse.

## 4. Explicit contraction invoice

Suppose F is one explicit synchronized polynomial factor of row degree
delta in the notation of the
[joint Frobenius theorem](joint-frobenius-product-contraction.md).
The w factors in (5) have the same row degree delta, with shifted rows.
For each coefficient in (3), the norm power n(F)^j is a product of these
Frobenius factors, each raised to j. Hence its total explicit degree
budget is at most

```
S_indicator <= (2L-1)*w*delta.                           (8)
```

Equation (4) is a linear combination of L such queries plus the known
constant sum. Equation (7) uses L+1 queries and the extra adjugate
product, giving

```
S_unit_inverse <= ((2L+1)*w-1)*delta.                    (9)
```

The [polynomial Galois-ring carry evaluator](teichmuller-cyclic-carries.md)
integrated from PR71 supplies the companion interface. Equations (8)-(9)
give its input budgets. For t dividing w, its
exponent parameter is `v*(w/t)`; the resulting cost has the form

```
O(L) * (S+2)^O(v*w/t) * poly(w,L,input length),          (10)
```

with polynomial bit overhead for the explicit ring operations. This uses
the credited evaluator, rather than supplying a new proof of it. It
does not assert a precision-independent S: precision appears linearly
inside the base in (10). Constants are summed with the actual cardinality
q^v, which can be nonzero modulo `2^L`.

### All denominator factors share one interpolation

There is no exponential branching cost in the number of inverse factors.
In this commutative local ring,

```
product_b inv*_L(L_b) = inv*_L(product_b L_b).            (11)
```

If any factor is a nonunit both sides vanish; otherwise both sides are
the inverse of the product. Let `B=product_b L_b`,
`A=product_a N_a`, `S_A=sum_a delta(N_a)` and
`S_B=sum_b delta(L_b)`. Keep the factors, without expanding B.
Multiplicativity of sigma and the norm gives

```
B# n(B)^j = product_b [ product_(h=1)^(w-1) sigma^h(L_b)
                        product_(h=0)^(w-1) sigma^h(L_b)^j ].
```

Every factor is still one shifted synchronized row. The entire product
`A product_b inv*_L(L_b)` is therefore a linear combination of **L+1**
polynomial carry queries, independent of the number of denominators, with

```
S <= S_A + ((2L+1)*w-1)*S_B.                            (12)
```

The individual denominators may have different rows. Equation (12)
replaces the invalid direct reuse of the field inverse exponent, while
charging precision inside S. It gives no precision-free state bound.

### Several residue equations can share one norm

Products of zero indicators need not be expanded either when their
explicit input admits a common synchronized row. For r polynomials F_i
over R, pass to an explicitly represented unramified extension
`R'=GR(2^L,r*w)` and choose beta_i whose residues form a basis over
`F_(2^w)`. The residue of `F=sum_i beta_i F_i` is zero exactly when
all residue equations vanish. Apply (4) using the norm of R'. If all F_i
have the same row and degree at most delta, F has that row and degree
at most delta; its factor budget is

```
S <= (2L-1)*r*w*delta,
```

with L+1 queries including the constant. The variables still range over
the original q-element Teichmuller set. Extension construction, basis
representation and ring bit costs must be included. Polynomials of
different rows cannot be packed into one row without justification:
the universal repacking `t=1`, with `D=w`, increases the carry exponent
to `v*w`. This is an explicit tradeoff, not a free closure assertion.

These formulas replace exponential precision degree by linear degree
for norm predicates and products of total unit inverses. Whether the remaining
degree, variable and representation costs can improve a complete factoring
algorithm is open; no such bound is claimed.

## 5. Implemented compilers and the system-count bound

The [compiler](../tools/precision_rational_carries.py) implements (11)-(12)
as `rational_sum` and the equation-packing construction as
`packed_system_count`. Its Frobenius-factor transform rotates relative
columns, wraps the binary row modulo the summation width, and applies
the actual ring automorphism to the coefficients. It does not enumerate
the Teichmuller domain. Both functions return an explicit resource invoice.

For a concrete system bound, let the explicit coefficient ring have degree
`m*w` and precision h, and let m equations in v variables have ordinary
degree at most delta. Their coefficients are supplied in the subring fixed
by `sigma^w`. The residues of

```
1, z, ..., z^(m-1)
```

form a basis over `F_(2^w)`, where z is the power-basis generator of
the supplied basic irreducible quotient. Indeed its residue generates
the whole degree-m*w field over F2, so also generates the degree-m
extension over `F_(2^w)`. A linear dependence of fewer than m powers
would contradict that minimal degree. This justifies the compiler's
particular packing basis without finding or averaging characters.

Choose `t=w`, so there is only one relative column per summation variable.
The exact modulo-`2^h` count has cost

```
(h*m*w*delta+2)^O(v) poly(input_length,h,m,w,v).          (13)
```

There are at most h+1 polynomial queries, each with
`S<=(2h-1)*m*w*delta`. The coefficient extension changes ring bit cost
and S, **without adding summation variables or precision**. The field
models and coefficient embeddings must already be supplied or acquired
at separately charged cost; (13) does not silently construct them.
For m=0 the known domain cardinality is returned. For h>v*w the count
is recovered exactly as an integer, with no nonunit precision division.

This complements PR71's additive-character compiler, whose alternative
invoice uses v+m variables and precision `h+w*(m-1)`. Both are retained.
Their relative cost depends on the available coefficient representation.
Equation (13) is a bound for explicit finite-field polynomial systems,
not for arbitrary word arithmetic. In particular trace predicates that
mix many binary rows can still enlarge the relative-column exponent;
the integer multiplication representation has not been supplied.

## 6. Exact controls

The [checker](../tools/check_precision_residue_indicators.py) expands both
versions of E_L and verifies the exact polynomial identities. It tests
all scalar residues at bounded precisions and all elements of the literal
quadratic unramified ring

```
(Z/2^L)[omega]/(omega^2+omega+1),
sigma(omega)=omega^2=-1-omega.
```

It checks norm, residue indicators, total unit inverses and direct
Teichmuller sums, including non-Teichmuller coefficients and a nonzero
even constant. Its [receipt](../evidence/precision-residue-indicator-controls.json)
retains the finite data, exact interpolation coefficients, charged
degree budgets and source hashes. Run `make verify-precision-indicators`.
No large parameter search, Lean or independent mathematical review is
included. Substantial OpenAI assistant support; the existing cyclic-carry
and active Galois-ring lift work retain their distinct source credit.

The separate [compiler checker](../tools/check_precision_rational_carries.py)
compares collected rational sums against literal ring unit inverses,
checks Frobenius-factor transforms pointwise, and compares packed system
counts with direct simultaneous root tests. Its
[receipt](../evidence/precision-rational-carry-controls.json) retains each
explicit factor, equation, output and invoice, including coefficient
extensions larger than the summation field and an exact count with
precision exceeding the domain bits. Run
`make verify-precision-rational-carries`. These are finite integration
controls for the ordinary compiler argument; no score gain is inferred.
