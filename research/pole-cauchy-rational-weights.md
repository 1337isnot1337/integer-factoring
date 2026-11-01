# All explicit separated rational weights

## Constructive result

Let `M` be an implicit square pole Cauchy matrix on explicit unions of
odd progressions, with its existing displacement and domain contract.
Let `A,P,B,Q` be explicitly written integer polynomials, and assume

\[
 P(a)\ne0\quad(a\in\mathcal A),\qquad
 Q(b)\ne0\quad(b\in\mathcal B).
 \tag{1}
\]

Then the complete contraction

\[
 \sum_{a\in\mathcal A,b\in\mathcal B}
      (M^{-1})_{b,a}\frac{A(a)}{P(a)}\frac{B(b)}{Q(b)}
 \tag{2}
\]

can be computed to any requested absolute 2-adic precision in polynomial
bit time. All polynomial degrees, coefficient lengths, interval endpoint
lengths, displacement, and output precision are charged. The numerical
matrix cardinality is not enumerated.

The algorithm checks (1). If a denominator vanishes on its domain, the
valuation partition rejects it with an actual integer node as witness.
Repeated roots and singular reductions modulo two are allowed. No
polynomial factorization or modular-root oracle is required.

This strictly extends the regular unit-inverse circuit contract to
arbitrary explicitly written separated rational functions. It does not
permit a compressed polynomial of exponential degree at the cost of its
short circuit, and it does not cover the joint denominator `N-ab`.
There is no new factoring exponent claim.

## 1. A polynomial-size valuation partition

Consider one nonzero denominator `P` of degree `D` on an odd domain.
Let `R` be the largest node. Set

\[
 H=\mathrm{bits}\left(\sum_{i=0}^{D}|P_i|R^i\right).
 \tag{3}
\]

Every nonzero value `P(a)` on the domain has valuation strictly below
`H`. The numerical value of `H` is polynomial in the coefficient lengths,
`D`, and the binary length of `R`.

Begin with `a=1+2X`. Remove the common power of two from the coefficients
of `P(1+2X)`. A state has the exact identity

\[
 P(\alpha+2^h X)=2^s f(X),
 \tag{4}
\]

where `0<=alpha<2^h`, `alpha` is odd, and `f` is primitive modulo two.
Its domain is the intersection with `a=alpha mod 2^h`; the size and
emptiness of that intersection are computed by linear congruence counts
on the given progression pieces.

Split `X=epsilon+2Y`, for `epsilon` zero or one. If `f(epsilon)` is odd,
the resulting polynomial has odd constant and its coefficient of degree
`i` is divisible by `2^i`. It is a unit leaf.

Otherwise remove the coefficient content `2^t` from
`f(epsilon+2Y)` and continue with that primitive polynomial. Here `t>=1`.
The number of continuing branches has an elementary degree budget.

### Root multiplicity controls the branch degree

Let `epsilon` have multiplicity `m` as a root of `f mod 2`. Write

\[
 f(\epsilon+Y)=\sum_j c_jY^j.
\]

The coefficient `c_m` is odd. In `f(epsilon+2Y)`, its valuation is
exactly `m`, so the removed content satisfies `t<=m`. Every coefficient
of degree greater than `m` has valuation greater than `m`, and therefore
vanishes modulo two after that content removal. Consequently the
normalized child's degree modulo two is at most `m`.

The multiplicities of the two possible roots sum to at most the parent's
degree modulo two. Thus the total active degree does not increase.
A normalized child of degree zero can use budget one, since it arose
from a root of multiplicity at least one; it has no further root
children. The active budget is at most `max(1,D)` at each depth.

Every continuing step increases `s` by at least one. If a nonempty state
has `s>=H`, every represented value is divisible by `2^H` while its
absolute value is smaller than `2^H`. Such a value must be zero. A
representative of the residue intersection supplies the rejected
integer witness.

Otherwise the continuing depth is at most `H`. There are
`O(max(1,D)(H+1))` processed states and leaves. Every state polynomial
has polynomial coefficient bit length: its variable substitution has
`O(H)`-bit coefficients, and its explicitly written degree remains at
most `D`.

### The final unit refinement must be retained

Content normalization alone does not ensure coefficient divisibility
by `2^i`. For example,

\[
 (a-1)^2+2=2(1+2X^2),\qquad a=1+2X.
\]

The quadratic coefficient of the unit polynomial `1+2X^2` has valuation
one. A further parity split gives the valid leaves `1+8Y^2` and
`3+8Y+8Y^2`.

The implementation stops only when the coefficient divisibilities are
explicitly satisfied. Otherwise even a unit state is split once more.
This increases the leaf count by at most a factor of two and prevents
an unjustified analytic truncation.

The [partition implementation](../tools/rational_weight_valuation_partition.py)
therefore returns disjoint nonempty residue leaves of the form

\[
 P(a)=2^s U(Z),\qquad a=\alpha+2^h Z,
 \qquad U_0\text{ odd},\quad 2^i\mid U_i.
 \tag{5}
\]

The valuation on each leaf is exactly `s`. In particular, its maximum
over all leaves is the exact maximum denominator valuation on the domain.

## 2. A local analytic quotient

On a leaf (5), write

\[
 \frac{A(a)}{P(a)}
      =2^{-s}\frac{A(\alpha+2^h Z)}{U(Z)}.
 \tag{6}
\]

Both numerator and denominator have coefficient of degree `i` divisible
by `2^i`; the numerator need not have a unit constant. Since `U_0` is
odd, formal division gives an integral series with the same coefficient
divisibility. Modulo `2^q`, only degrees below `q` are needed:

\[
 C(Z)=\sum_{i=0}^{q-1}c_i Z^i,
 \qquad
 C(Z)\equiv A(\alpha+2^h Z)/U(Z)\pmod{2^q}
 \quad(Z\in\mathbb Z_2).
 \tag{7}
\]

The coefficients are obtained by the ordinary triangular division
recurrence modulo `2^q`. Their only denominator is the odd constant
`U_0`, so this step has no uncharged division by two.

Converting `Z` directly into `(a-alpha)/2^h` would introduce many dyadic
coefficient denominators away from the leaf. The quadratic selector
gives a more useful exact finite-polynomial substitution.

## 3. Local weights become global pole polynomials

Put

\[
 c=\alpha-2^h w,\qquad
 y=\frac{2^h}{a-c},\qquad Z=y^{-1}-w.
 \tag{8}
\]

The variable `y` is integral at every integer `a`. It is a unit on the
specified residue class and has positive valuation off that class.
The cubic `y^3` has residue one on the class and zero elsewhere.

Let `P_q` denote the [Hermite selector polynomial](pole-cauchy-residue-blocks.md)
of degree `2q-1`, supported on degrees `q,...,2q-1`. Define the global
approximation contributed by this leaf as

\[
 \Phi_{\alpha,h,s}(a)
       =2^{-s}P_q(y^3)C(y^{-1}-w).
 \tag{9}
\]

This is an explicitly written proper pole polynomial at `c`. Indeed,
expanding a term of `C` gives

\[
 y^{3j}(y^{-1}-w)^i
   =\sum_{\ell=0}^i {i\choose\ell}(-w)^{i-\ell}
                        y^{3j-\ell}.
 \tag{10}
\]

Here `q<=j<=2q-1` and `i<q`, so every exponent lies between `2q+1`
and `6q-3`. All coefficient pairs are integral in the basis `1,w`.
Replacing `y^t` by `2^(ht)/(a-c)^t` leaves only the external denominator
`2^s` from (9).

On the leaf, `Z` is integral, (7) is valid, and `P_q(y^3)=1 mod 2^q`.
Thus (9) approximates the desired rational weight with absolute error
valuation at least `q-s`.

Off the leaf, every power in (10) has valuation at least `2q+1`.
Consequently (9) itself has valuation at least `2q+1-s`. This proves
the required suppression directly. It does not substitute a divergent
infinite series at an off-leaf input.

Sum (9) over all leaves. If `u=max_a v2(P(a))`, the resulting global
weight `Phi(a)` satisfies

\[
 v_2\left(\Phi(a)-\frac{A(a)}{P(a)}\right)\ge q-u
 \quad(a\in\mathcal A).
 \tag{11}
\]

At each node one leaf supplies its local approximation, and all other
leaves are suppressed by the off-leaf bound. Summing them does not lose
2-adic precision. Zero truncated numerator jets may be omitted with the
same error guarantee.

## 4. Precision for the full inverse contraction

Let `L` certify integrality of `2^L M^-1`. Put

\[
 u=\max_a v_2(P(a)),\qquad v=\max_b v_2(Q(b)),\qquad
 q=p+L+u+v.
 \tag{12}
\]

Both maxima are obtained exactly from the two valuation partitions.
The true row and column weights have valuations at least `-u` and `-v`.
The approximations in (11) preserve those lower bounds, while their
errors have valuations at least `q-u` and `q-v`. Therefore the product
weight error has valuation at least `q-u-v=p+L`. Multiplication by an
inverse entry loses at most `L` bits, and the full sum retains absolute
precision `p`.

For each pair of leaves, the
[consumer](../tools/pole_cauchy_rational_weights.py) calls the batched
pole-polynomial evaluator on the integral coefficient arrays from (10).
If the two external denominators are `2^s` and `2^t`, it requests
absolute precision `p+s+t` before dividing by `2^(s+t)`. All returned
terms are combined with a common explicit scale; no inversion of two
is treated as a unit operation.

An optional normalized-inverse mode computes (2) with `2^L M^-1` in
place of `M^-1`. It puts `2^L` directly into the row coefficients and
uses `q=p+u+v`. This saves the unnecessary local approximation precision
when the matrix has already been made integral. The weighted answer
can still be nonintegral because the rational weights may have even
denominators, so its dyadic scale remains part of the return value.

## 5. Complete cost and explicit representation boundary

The denominator partitions have polynomially many leaves by Section 1.
Their local series have degree below `q`; their pole polynomials have
degree at most `6q-3`. Formula (10) constructs their coefficients with
polynomially many integer operations and polynomial coefficient bit
lengths. Pairing the two leaf lists still gives polynomially many
batched contractions. All numerical precisions are polynomial by
(3), (12), and the existing inverse-loss certificate.

This proves polynomial bit cost for the full procedure. Polynomial
degrees are measured from the explicitly supplied coefficient arrays.
A product such as `prod_b(N-ab)` over an exponentially large domain is
not an admissible low-cost denominator description under this contract.
Allowing a joint polynomial denominator in both `a,b` would require a
different argument; the present tree is univariate and preserves the
separation of the row and column weights.

Thus the result expands the evaluable rational weight family without
claiming an implicit Smith count, global product-pole trace, or faster
factoring algorithm.

## 6. Exact controls and large implicit examples

The [checker](../tools/check_pole_cauchy_rational_weights.py) records
487 nonzero-denominator partitions and 3,833 literal leaf identities.
Ten actual domain zeros are rejected. It verifies 236 pointwise rational
weight approximations and twelve complete dense inverse contractions,
covering both normalized and unnormalized inverse modes. The final unit
refinement example from Section 1 is an explicit regression control.

A deep implicit partition uses the denominator

\[
 (a-1)^2+2^{140}
\]

on `2^200+13` odd nodes. Its maximum denominator valuation is 141.
The partition uses 71 processed primitive states and 72 leaves, with
active degree budget at most two. No nodes are listed.

The complete large contraction uses the 402-bit input
`N=2^401+45`, the same dimension, and the rational weights

\[
 \frac{a+1}{(a-1)^2+2},\qquad
 \frac{2-b}{(b-3)^2+2}.
\]

In normalized-inverse mode, the two-bit computation uses four batched
contractions and local precision four. Its scaled output is
`(1+w)/2 mod 4`, agreeing with a four-bit replay. The recorded lower
precision run takes about 0.88 seconds, with at most 39 internal bits
per batched contraction. These timings are finite implementation
controls, not factoring complexity claims.

The [evidence](../evidence/pole-cauchy-rational-weights.json) records the
partitions, exact precision budgets, results, and source hashes.

```sh
python3 tools/check_pole_cauchy_rational_weights.py \
  --output evidence/pole-cauchy-rational-weights.json
```
