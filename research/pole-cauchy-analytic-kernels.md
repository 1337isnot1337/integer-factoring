# Complete inverse contractions with joint analytic kernels

## Result and scope

The [regular weight compiler](pole-cauchy-analytic-weights.md) extends
to a circuit depending jointly on both matrix indices. In particular,
the kernel need not be given as a short sum of separated functions.
For two equal-length positive odd progressions, a complete inverse
contraction with a regular bivariate circuit has polynomial bit cost
in the encoded input, circuit and power-exponent lengths, and charged
2-adic precision. Neither progression is enumerated.

A second improvement is an exact quadratic recurrence for all mixed
polynomial moments. It replaces the separate inner sum in each moment
by one diagonal-prefix update. The required data are still only two
short interpolation jets.

Regularity excludes inverse gates whose denominators are even on the
odd input class. The singular factor-count kernel does not satisfy this
contract. The factoring exponent remains one half.

## 1. Matrix, normalization, and the contraction

Use the notation of [polynomial inverse moments](pole-cauchy-polynomial-moments.md):

\[
 t=2^kw,\qquad z=N+t,\qquad w^2+w+1=0,
 \qquad M_{a,b}=\frac{t}{z-ab}.                    \tag{1}
\]

Here `N` is positive and odd, `A,B` have equal cardinality, and `k>=1`
is large enough that every nonzero `N-ab` has absolute value below
`2^k`. Write

\[
 U=2^kM^{-1},\qquad
 \mathcal C(K)=\sum_{a\in A,b\in B}U_{b,a}K(a,b).
 \tag{2}
\]

For exact polynomial kernels, the domains may be unions of disjoint
positive odd progressions. For a regular circuit supplied through
pointwise finite-precision compilation, require that each domain is
one progression of step two, with the same length. The
[balanced conditioning theorem](pole-cauchy-conditioning.md) then makes
every entry of `U` integral. Thus pointwise congruence of kernels modulo
`2^P` implies congruence of their entire contractions modulo `2^P`.
There is no loss proportional to the number of matrix entries.

## 2. Two jets and a diagonal-prefix recurrence

Put `x_a=z/a` and define the integral formal series

\[
 R(Y)=\frac{\prod_{b\in B}(Y-b)}
            {\prod_{a\in A}(Y-x_a)}
     =\sum_{d\ge0}r_dY^d,
 \qquad
 T(Y)=\frac{\prod_{a\in A}(1-x_aY)}
            {\prod_{b\in B}(1-bY)}
     =\sum_{e\ge0}t_eY^e.                         \tag{3}
\]

Let `mu_(d,e)=C(a^d b^e)`. Reindexing the established interpolation
identity gives

\[
 \mu_{d,e}=\frac{z^{d+1}}w
    \left(\mathbf1_{d=e}-S_{d,e}\right),
 \qquad
 S_{d,e}=\sum_{\ell=0}^{\min(d,e)}r_{d-\ell}t_{e-\ell}.
 \tag{4}
\]

The sum has the recurrence

\[
 \boxed{S_{d,e}=r_dt_e+S_{d-1,e-1}},              \tag{5}
\]

where a negative index means zero. Therefore all moments with
`0<=d<=D`, `0<=e<=E` require `O((D+1)(E+1))` scalar ring operations
after the two jets are known. A row-by-row evaluation stores just
`E+1` previous values. Powers of `z` are updated once per row.
The coefficients `r_0,...,r_D` and `t_0,...,t_E` are computed by
[modular interval jets](modular-interval-jet-arithmetic.md), with
polynomial cost in their truncation degrees, precision, and the binary
lengths of progression endpoints and sizes.

Equivalently, the bivariate moment series has the compact identity

\[
 \boxed{\displaystyle
 \sum_{d,e\ge0}\mu_{d,e}X^dY^e
 =\frac zw\,\frac{1-R(zX)T(Y)}{1-zXY}.}           \tag{6}
\]

Indeed, each diagonal shift in (4) contributes one further factor
`zXY`, while the unshifted terms sum to `R(zX)T(Y)`.
Both (4) and (6) are exact formal identities. They do not assert that
an arbitrary singular kernel can be evaluated by substituting a
nonconvergent series into them.

For an explicit integer polynomial `K=sum c_(d,e) a^d b^e`, take the
corresponding finite linear combination of (4). The numerical degree
bounds are charged. The implementation accepts sparse input terms,
but its recurrence cost uses their rectangular degree bounds.

## 3. Compiling a regular joint circuit

The circuit gates are integer constants, either input variable,
addition, subtraction, multiplication, binary-encoded integer powers,
and inverses of values odd at `(a,b)=(1,1)`. Negative powers have the
same unit requirement. Substituting

\[
 a=1+2X,\qquad b=1+2Y
\]

gives, at every accepted gate,

\[
 F(1+2X,1+2Y)=\sum_{i,j\ge0}c_{i,j}X^iY^j,
 \qquad c_{i,j}\in2^{i+j}\mathbb Z_2.            \tag{7}
\]

Inputs and constants satisfy (7). Sums and products preserve it.
If a series has odd constant coefficient, its inverse does too:
solve by increasing total degree, and each nonconstant product in the
recurrence has at least the required `i+j` powers of two. Consequently
every permitted denominator is odd at every pair of odd integers.

Work modulo `2^P` and discard total degrees at least `P`. Every
discarded term vanishes modulo `2^P` on integral `X,Y`. For retained
coefficients, choose the representatives `0<=c_(i,j)<2^P` and form

\[
 F_P(a,b)=\sum_{i+j<P}\frac{c_{i,j}}{2^{i+j}}
                     (a-1)^i(b-1)^j.            \tag{8}
\]

This is an integer polynomial of total degree at most `P-1`. It agrees
with the circuit modulo `2^P` at every odd pair. Although the division
retains only `P-i-j` coefficient bits, the product of shifted variables
restores the missing `i+j` bits pointwise. No unjustified coefficient
precision is used. Expanding (8) gives the monomial terms for (4).

The implementation uses Newton inversion `B <- B(2-AB)` in the
truncated total-degree ring. Its initial constant is the exact unit
inverse modulo `2^P`; each update doubles the minimum degree of the
remaining error. Binary exponentiation handles signed encoded powers.

## 4. Precision and bit cost

For the regular kernel `F` multiplied by an external scalar `2^-r`,
the algorithm returns an integral residue `C_P` with

\[
 \sum_{a,b}(M^{-1})_{b,a}\,2^{-r}F(a,b)
 =2^{-k-r}C_P+O(2^{P-k-r}).                       \tag{9}
\]

Thus `P=p+k+r` suffices for `p` absolute output bits. The value of
the external denominator exponent is charged, even when its binary
encoding is short. Formula (9) specifies a scale and an error bound;
it does not claim that the output has exactly that valuation.

There are at most `P(P+1)/2` retained bivariate coefficients. Naive
truncated multiplication and basis conversion cost `O(P^4)` scalar
operations. If `S` is the ordinary gate count, `E` is the sum of power
exponent encoding lengths, and `I` is the number of unit inversions
including negative powers, compilation costs

\[
 O\bigl((S+E+I\log(P+1)+1)P^4\bigr)              \tag{10}
\]

operations on integers of polynomial bit length. Constants' input
lengths are charged. The two jet computations have polynomial cost;
the final moment recurrence uses at most `P^2` cells. All costs depend
on the binary encoding of the progression lengths, not on their
numerical cardinality.

## 5. Exact controls and an implicit joint example

The [checker](../tools/check_pole_cauchy_analytic_kernels.py) compiles
20 circuits at seven precisions and checks all 109,220 odd input pairs.
It compares 60 complete joint contractions against dense exact matrix
inversion and checks 30 instances of the quadratic recurrence against
the previous separate-polynomial evaluator. Three invalid or singular
circuits are rejected. The implementation is
[pole_cauchy_analytic_kernels.py](../tools/pole_cauchy_analytic_kernels.py).

The large control uses `N=2^401+45`, `m=2^200+13`,
`A={1,3,...,2m-1}`, and `B={3,5,...,2m+1}`, with the joint kernel

\[
 \begin{aligned}
 F(a,b)={}&(1+2^{64}ab)^{-(2^{512}+2^{257}+123)}\\
         &\cdot(1+2^{80}(a+b))^{-(2^{383}+2^{193}+7)}.
 \end{aligned}                                    \tag{11}
\]

At twelve absolute output bits, `k=414` and `P=426`. The compiled
total degree is twelve, with degree at most six in either variable;
the recurrence uses 49 cells. A second computation at sixteen absolute
bits agrees at the requested precision. Zero matrix nodes are
enumerated. The [evidence](../evidence/pole-cauchy-analytic-kernels.json)
records the circuit, precision, coefficients, residues, operation
counts, and source hashes.

```sh
python3 tools/check_pole_cauchy_analytic_kernels.py \
  --output evidence/pole-cauchy-analytic-kernels.json
```

The factor-count trace uses `(z-ab)^(-2)`, whose denominator is even
on odd inputs. It is not covered by the unit-inverse compilation rule.
Resolving that singular aggregation remains a separate requirement
for a factoring improvement.
