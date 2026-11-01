# Summing relative Frobenius factors by a vector of carries

## Status

This note extends the [rational-star primitive](frobenius-rational-trace-contraction.md)
to affine linearized polynomials with several terms. The cost is
exponential only in the relative Frobenius dimension, with a polynomial
dependence on the binary length of the base field. It gives a second
explicit family of coupled sums with a charged exact algorithm.

The result also handles vanishing denominators using a specified
zero-preserving inverse. This is useful for geometric digit boxes, whose
removable-pole value in characteristic two is zero. The result does not
establish closure of an entire multiplication network under repeated
elimination, and does not change the complete factoring exponent.

## 1. The theorem and its representation

Let `q=2^t`, `b=q^D=2^w`, where `w=tD`. Let the explicit coefficient
field `E=F_(2^d)` contain `F_b`, so `w|d`. An affine relative-linearized
factor with row `rho`, `0<=rho<t`, means

\[
 L(U)=a_*+\sum_{j=0}^{D-1}a_j U^{2^\rho q^j}.              \tag{1}
\]

Coefficients may vanish. Different factors can have different rows. Define

\[
 \iota(x)=\begin{cases}x^{-1}&x\ne0,\\0&x=0.\end{cases}
 \qquad
 F(U)=c\prod_{i=1}^{r_+}N_i(U)
             \prod_{j=1}^{r_-}\iota(L_j(U)).               \tag{2}
\]

Every `N_i` and `L_j` has form (1). Thus (2) agrees with the ordinary
rational function wherever all denominators are nonzero; elsewhere its
meaning is the explicitly stated product with `iota(0)=0`.

Put

\[
 s=r_++dr_-,\qquad
 J=D+\left\lceil\log_q(s+1)\right\rceil+1.                \tag{3}
\]

The two sums

\[
 \sum_{u\in F_b}F(u),\qquad
 \sum_{u\in F_b^\times}F(u)                              \tag{4}
\]

are computable deterministically in

\[
 (s+2)^{O(J)}\operatorname{poly}(d,t,D)                   \tag{5}
\]

field operations and polynomial additional bit overhead per operation.
The algorithm does not enumerate `q` or `b` elements. Field construction
is separate and must be charged by the caller. In particular, when
`t>=log_2(s+1)`, the relative dimension `D`, rather than `w=tD`, controls
the exponent in (5).

This theorem concerns the explicit product representation (2). It is not
a bound for an arbitrary arithmetic circuit with the same input factors.

## 2. Replace inverse values by a short Frobenius product

For every `x in E`, including zero,

\[
 \iota(x)=x^{2^{d+1}-3}.                                  \tag{6}
\]

The exponent is positive, even when `d=1`, and is congruent to `-1`
modulo `2^d-1`. Its binary expansion has exactly `d` nonzero bits: bit
zero and bits `2,...,d`, with the latter list empty for `d=1`. Consequently

\[
 \iota(L(U))=\prod_{h\in H_d}L(U)^{2^h},
 \qquad H_d=\{0\}\cup\{2,\ldots,d\}.                    \tag{7}
\]

On `F_b`, the identity `U^(2^w)=U` lets us reduce each Frobenius exponent
of `U` modulo `w`. For a factor of row `rho`, its power in (7) becomes

\[
 a_*^{2^h}+\sum_j a_j^{2^h}
        U^{2^{(\rho+h+tj)\bmod w}}.                       \tag{8}
\]

All nonconstant terms in (8) have one common row modulo `t`; their
column indices are cyclically permuted. Thus (2) has the same values on
`F_b` as a polynomial `G(U)` which is a product of `s` factors of form
(1), multiplied by `c`. This is a value identity; it need not be an
identity of formal rational functions outside that field.

Its ordinary polynomial degree is at most `s 2^(w-1)`. We can therefore
sum its values by computing only

\[
 \sum_{\ell=1}^{\lfloor s2^{w-1}/(b-1)\rfloor}
                  [U^{\ell(b-1)}]G(U).                   \tag{9}
\]

Indeed, `sum_(u in F_b)u^k` is zero for `k=0`, and is one precisely
when `k>0` is divisible by `b-1`. The number of coefficient queries in
(9) is at most `s`. The chosen `J` ensures every queried degree is below
`q^J`. If `s=0`, the function is constant and the full-field sum is zero.

The nonzero-field sum is obtained by adding `F(0)`. Formula (6) evaluates
this correction without division or a zero test.

## 3. Coefficient extraction with one carry per base-`q` digit

Consider any product of `s` factors of form (1), and a target integer
`0<=K<q^J`. For each row `rho` form the multivariate polynomial

\[
 P_\rho(Y_0,\ldots,Y_{J-1})
  =\prod_{i:\rho_i=\rho}
       \left(a_{i,*}+\sum_{j=0}^{D-1}a_{i,j}Y_j\right).   \tag{10}
\]

Its coefficient at `Y^h` records the total weight of choices contributing
`h_j` copies of the binary place `2^rho q^j`. Every coordinate of `h` is
at most `s`, and their sum is at most `s`. Direct sparse multiplication
computes (10) in `(s+2)^O(J)` field operations.

Let `k_(rho,j)` be bit `rho+tj` of `K`. A transfer at row `rho` takes an
incoming vector of integer carries `a=(a_0,...,a_(J-1))` to an outgoing
vector `z=(z_0,...,z_(J-1))`. Its weight is

\[
 T_\rho(a,z)=
 [Y^h]P_\rho(Y),\qquad
 h_j=2z_j+k_{\rho,j}-a_j,                                 \tag{11}
\]

interpreted as zero when any coordinate is negative or exceeds `s`.
All carries may be restricted to `[0,s]`, giving at most `(s+1)^J`
states. Multiply these `t` transfer matrices in row order, or propagate
each initial state through the sequence.

The boundary condition joins the last binary row of one base-`q` digit
to the first binary row of the next digit. If `a` is the initial vector
and `z` the final vector, retain exactly

\[
 a_0=0,\qquad a_{j+1}=z_j\ (0\le j<J-1),\qquad z_{J-1}=0.\tag{12}
\]

Sum the corresponding entries of `T_0...T_(t-1)`. This is `[U^K]G(U)`.

For a proof, fix all monomial choices from the original factors. They
specify the counts `h_(rho,j)`. The binary addition of their exponents
has a unique carry sequence when its initial carry is zero. Equations
(11)--(12) are exactly that sequence, with its coordinates processed in
a different order. The final zero condition forbids overflow above
`q^J`. Thus a monomial choice contributes precisely when its total degree
is `K`, and contributes its original product of coefficients exactly
once. Linearity proves the formula. There is no extra cyclic carry:
the conditions begin at digit zero and end with a zero terminal carry.

Naive matrix multiplication costs `t(s+1)^(3J)` operations. Together
with (10), the at most `s` coefficient queries in (9), and the polynomial
work forming (8), this proves (5). Faster transfer algorithms are not
needed for the stated theorem.

## 4. Geometric boxes and pole differences

For a field value `z`, let `S_b(z)=sum_(a=0)^(b-1) z^a`. Since `b` is
even, `S_b(1)=0`. Therefore the following identity holds for **every**
`z in E`, including `z=1`:

\[
 S_b(z)=(1+z^b)\,\iota(1+z).                             \tag{13}
\]

When `z` is an affine relative-linearized function of the summation
variable, both factors on the right remain of form (1) after reducing
Frobenius exponents modulo `w`. Theorem 1 therefore applies without a
nonvanishing-denominator promise. The field representation still has
to be explicit and its degree charged.

The difference of two relative Frobenius monomials,
`aU^(q^i)+cU^(q^j)`, is also a factor of form (1). Thus this theorem can
handle the new linearized denominators noted after the first rational
star contraction, provided the entire resulting expression has a
controlled product representation as in (2).

## 5. Remaining network obligation

One elimination can output a sum or circuit of many such products. The
number of terms, their field degree, and their dependence on remaining
variables must all be charged before applying the theorem again.
Nothing above proves that those quantities remain controlled through
the full factoring network.

In particular, an additive multiplication-gate phase `zeta^(sxy)` has
an exponent depending on other summation variables. It is not thereby
a fixed relative Frobenius monomial in those variables. This note does
not replace that gate by a free factor of the form (1), and does not
remove the coupled operand/frequency operation identified in the
[weighted carry contraction](weighted-dyadic-carry-contraction.md).

## 6. Exact replay

The companion checker is
[`check_relative_frobenius_trace.py`](../tools/check_relative_frobenius_trace.py).
It compares the vector-carry coefficient algorithm with ordinary exact
polynomial multiplication and compares both field sums with enumeration.
It includes nonzero cross-digit carry, final overflow rejection, zero
denominators, vanishing numerator coefficients, and `d=1`. Run:

```sh
python tools/check_relative_frobenius_trace.py --output evidence/relative-frobenius-trace.json
```

The committed [evidence](../evidence/relative-frobenius-trace.json) records
437 exact vector-carry coefficient comparisons, 240 complete relative
traces checked against 1,280 individual field evaluations, 129 evaluations
at a denominator zero, 46 zero-preserving inverse controls, 46 all-point
geometric-box controls, and three cross-digit/overflow controls. The
transfer computations checked 7,268 initial carry vectors and 9,128
row-transfer states. All comparisons passed.

The proof is self-contained; the finite replay checks the formulas and
their implementation without replacing the arbitrary-parameter proof.
