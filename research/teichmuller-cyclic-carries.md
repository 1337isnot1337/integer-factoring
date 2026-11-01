# Cyclic carries retain 2-adic precision on Teichmuller inputs

The polynomial part of the [joint Frobenius carry evaluator](joint-frobenius-product-contraction.md)
works over finite unramified rings with the **same carry states**. The
boundary weights must retain the actual power sums, rather than their
characteristic-two reductions. This supplies two positive interfaces:

1. Exact synchronized polynomial sums modulo `2^h`, with polynomial bit
   overhead in the precision `h` and coefficient extension degree.
2. Finite-field root counts modulo `2^h`, using an explicitly acquired
   indicator of degree `2h-1` in a lifted norm. The naive exponential
   precision degree is unnecessary.

For a residue-field system of `m` equations of degree at most `delta` in
`v` variables over `F_(2^d)`, one resulting bound is

```
(H*(delta+1)+2)^O(v+m) poly(input_length,H,d,v,m),
H=h+d*(m-1).                                             (1)
```

Field representation construction is separate and precision is a numeric
cost parameter. This is a finite-field polynomial-system interface, not
an encoding of the complete integer factoring query. It supplies no
stronger factoring exponent; the complete score remains **1/2**.

The earlier [fixed-precision network lift](fixed-precision-count-lifting.md)
applies to a broader network class with precision-dependent constants.
Here precision may grow, but the input must have the explicit polynomial
structure stated below. No resource claim about arbitrary networks follows.

## 1. Coefficient ring and exact moments

Let `O` be the unramified extension of the 2-adic integers with residue
field `F_(2^d)`, and `R=O/(2^h)`. Fix `w|d`, put `b=2^w`, and let `T_w`
be the Teichmuller representatives of the subfield `F_b`. They are zero
and the `(b-1)`st roots of unity. For a monomial exponent `K>=0`,

```
sum_(x in T_w) x^K = b                  if K=0,
                     b-1              if K>0 and b-1 divides K,
                     0                otherwise.         (2)
```

For `T_w*`, replace the first value by `b-1`. Constants in (2) mean their
actual images in `R`; neither `b=0` nor `b-1=1` is assumed.

The nonzero representatives form a cyclic group of odd order. If `K` is
not divisible by `b-1`, choose a generator `zeta`; `zeta^K-1` has nonzero
residue and is a unit. Multiplying the geometric sum by this unit proves
its vanishing. The other cases are direct counts, including the constant
monomial at zero.

Existence, multiplicative representatives and functorial Frobenius lifts
are standard facts from Rabinoff,
[*The Theory of Witt Vectors*, arXiv:1409.7445v1](https://arxiv.org/abs/1409.7445v1),
Theorem 1.2 and Example 1.3. Lemma 1.4 supplies a finite-precision power
construction of the representatives. The exact carry weights and the count
compiler below are derived here using those facts.

## 2. The state envelope is unchanged

Use exactly the explicit polynomial factors of the credited joint theorem.
Fix `t|w`, put `C=w/t`, and give each factor a row `rho` and a finite list

```
A(U)=sum_mu a_mu product_(i,j) (U_i^(2^(rho+t*j)))^(mu_(i,j)),
0<=rho<t,       0<=j<C.
```

Coefficients lie in `R`. Let `S` be the sum of the factors' maximum total
count degrees. All count vectors are explicit nonnegative integers.
No inverse factors are included in this theorem.

Group factors by row and multiply their coefficient tables as before, using
ring addition and multiplication. The same transition has

```
incoming a, outgoing z, coefficient count h=2z-a.
```

The same cyclic boundary rotates each variable's relative columns. For a
fixed monomial choice, its exponent `K` satisfies

```
K=(b-1)c_0,
```

with one unique nonnegative carry cycle if such a cycle exists. Every carry
is at most `S`. Moreover `c_0=0` exactly when `K=0`. These are integer
identities independent of the coefficient ring.

Allow zero initial carries for both domain types. Multiply each completed
path by one boundary scalar for each variable:

```
b     if its domain is T_w and its initial carry is zero;
b-1   otherwise.                                         (3)
```

Equation (2) proves that the resulting sum is exact in `R`. In precision
one these weights reduce to the credited field boundary. In higher
precision the constant term and the sign/size of `b-1` are retained.

The algorithm still uses at most `(S+1)^(vC)` carry vectors and coefficient
positions, up to the same polynomial-power matrix contraction as the
original theorem. Its cost is

```
(S+2)^O(vC) poly(input_length,h,d,w,v).                    (4)
```

Ring elements occupy `h*d` bits. Given a basic irreducible polynomial,
quotient-ring operations, Hensel acquisition of Frobenius, and Teichmuller
lifting have polynomial bit cost in `h,d`. The implementation describes
them explicitly; it does not enumerate the `2^w` summation points.
Degree and coefficient-table acquisition remain part of the input invoice.

## 3. An explicit precision-linear residue-root indicator

Take `w=d`, and let `f` be an explicit polynomial of ordinary degree at
most `delta`. Let `sigma` be the lifted Frobenius. On Teichmuller inputs,

```
N_f(U)=product_(r=0)^(d-1) sigma^r(f)(U_1^(2^r),...,U_v^(2^r))
```

is the unramified norm of `f(U)`, lying in `Z/(2^h)`: the product is fixed
by Frobenius, whose fixed ring is identified by the unique Teichmuller digit
expansion. No division by the possibly even extension degree is used.
Its parity is zero
exactly when the residue polynomial vanishes. It is generally not itself
a zero/one value in higher precision.

Define the integer polynomial

```
H_h(T)=sum_(k=h)^(2h-1) binom(2h-1,k) T^k (1-T)^(2h-1-k).
```

It has degree `2h-1`. For even `T`, every summand is divisible by `2^h`.
For odd `T`, expand `1-H_h(T)` using the complementary binomial terms;
each contains `(1-T)^h`. Therefore

```
H_h(T) = 0 mod 2^h for even T,
H_h(T) = 1 mod 2^h for odd T.                             (5)
```

Consequently `1-H_h(N_f(U))` is the exact residue-root indicator modulo
`2^h`. Its coefficient list has only `h` possibly nonzero terms, at powers
`h,...,2h-1`, with explicitly computed integer coefficients of polynomial
bit size.

Acquire each norm moment `N_f^j` by repeating the supplied row factors
`j<=2h-1` times. There is no exponential polynomial expansion before the
row-table acquisition. Each row has degree at most `j*delta`, so its carry
coordinates are at most `j*delta`: the rotated degree bound is
`j*delta*(2^d-1)`. This gives

```
(h*delta+2)^O(v) poly(input_length,h,d,v)                  (6)
```

bit cost for one residue-polynomial root count. An already supplied
synchronized polynomial weight can be included in each moment; its written
degree and row structure must also be charged.

This counts roots of `f mod 2` in the residue field. It is not a compiler
for the different congruence `f=0 mod 2^h` on arbitrary ring-valued inputs.

## 4. Several equations: pay variables and division precision

Let `f_1,...,f_m` be residue-field polynomials in `v` variables and let `Z`
be their joint zero count. For `m>=1`, introduce `m` full-field variables
`Y_i` and count roots of

```
F(U,Y)=sum_i Y_i f_i(U).
```

If the vector of `f_i(U)` is zero, there are `b^m` choices of `Y`.
Otherwise there are `b^(m-1)`. Thus its root count satisfies the integer
identity

```
T = b^(v+m-1) + b^(m-1)(b-1) Z.                         (7)
```

To obtain `Z mod 2^h`, acquire `T mod 2^H` with
`H=h+d*(m-1)`. Subtract the known first term, divide the resulting residue
by the explicit power `2^(d*(m-1))`, and invert the odd unit `b-1` modulo
`2^h`. The division is well defined at exactly that charged precision.
Applying (6) to `v+m` variables and degree `delta+1` proves (1). The case
`m=0` is the known domain size.

The extra variables and precision can be costly. No unchanged-width claim
for an arbitrary verifier or integer arithmetic program is inferred.

## 5. What this does and does not supply

This gives a concrete precision-preserving backend for **explicit
synchronized polynomial inputs**, and a residue-root count compiler with
precision growing in the cost. It is relevant to exact moments and counts
whose coefficient ring is `Z/(2^h)`, including the interface highlighted by
the [semiprime intake](sanftenberg-semiprime-intake.md).

It does not yet compile that public multiplication network into the required
polynomial form. Nor does it provide a cheap native word multiplication,
integer prefix/carry representation, or a stronger complete factoring bound.
The [elliptic word audit](elliptic-word-multiplication-cost.md) retains its
separate graph-degree and auxiliary-fibre limits.

The field inverse-value identity cannot be copied to arbitrary ring inputs.
For example, `3` has inverse `11 mod 16`; a field exponent-one rule would
return `3`. Nonunits need not have any inverse. The model's Hensel steps
invert acquired units, and (7) additionally inverts the known odd scalar
`b-1`. The later [collected-inverse companion](precision-efficient-residue-indicators.md)
supplies an inverse-factor compiler for total unit inverses, assigning zero
to nonunits, with at most `h+1` polynomial queries for the entire denominator
list and an explicit precision-dependent degree budget. Its
coefficient-extension packing also gives a second residue-system interface
without the extra summation variables or precision shift in (1). The
original additive-character interface and its exact receipt remain retained.

## 6. Reproduction and credit

Run `make verify-teichmuller-carries`. The
[ring and carry code](../tools/teichmuller_cyclic_carries.py) retains the
credited field algorithm separately. The
[checker](../tools/check_teichmuller_cyclic_carries.py) compares 384 polynomial
queries and six norm sums with 7,060 direct assignments, alongside 96
reductions to the existing field evaluator, 256 Teichmuller moments and 384
Frobenius homomorphism checks. It also checks twelve independent expansions of (5), thirteen
residue-root queries and eight system queries, including nonunit division
precision. The [receipt](../evidence/teichmuller-cyclic-carry-controls.json)
pins the seed, sources, quotient models and query digest. Regenerate with:

```
python3 tools/check_teichmuller_cyclic_carries.py --generate
```

The [clean full replay](../evidence/teichmuller-carry-full-verification.json)
pins the default and 35 focused targets at exact `b9cd96f`, with 384 Python
command blocks, 66 unit tests, five crossing tests and 1,129 Markdown files.
It predates the collected-unit-inverse implementation and the companion
precision-indicator addition; those later checks have separate scope.

The negative control over precision six is substantive: the constant
polynomial `3` has no residue roots, but the naive sum of `1-N_f` is
`32 mod 64`. The acquired polynomial (5) returns zero as required.

The cyclic-carry construction and its synchronization representation retain
credit to their existing contributors. Standard unramified/Witt facts are
credited above; the precision boundary and explicit count invoice are this
checkpoint. Substantial OpenAI assistant support was used. These are ordinary
proofs and finite controls, with no new Lean formalization or independent
mathematical review claimed.
