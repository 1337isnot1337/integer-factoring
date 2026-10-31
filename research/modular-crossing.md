# A normalized crossing construction in fixed characteristic

## Result and scope

The factorial multiplier in the original partition crossing is removable
in fixed characteristic. For every fixed prime `p`, this note constructs
an effective family over `F_p` whose pieces sum to the crossing tensor
with multiplier **one**. The construction extends to every commutative
`F_p`-algebra, including rings with zero divisors.

The state space has a basis of partitions with no part divisible by `p`
and with every multiplicity less than `p`. All source operator identities
survive on this quotient. Label generation and factor evaluation take
`2^(O(lambda))` bit operations for fixed `p`.

An elementary generating-function argument gives, for every fixed
`0<epsilon<1-1/sqrt(2)`, the source's width-sum condition with
`rho=1/sqrt(2)+epsilon`. This is a new normalized modular construction,
not a bound beyond exponent `1/2`. The smaller state space still has
square-root entropy. A separate actual-rank theorem would be needed for
a smaller `rho`.

There is also a homogeneous variant at prime-power energies. When
`h=p^r`, its pairing is invertible in the quotient. At `p=2`, only two
angle indices survive. Exact central ranks at `h=8,16,32` are
`2,12,261` on the full quotient spaces, whose dimensions are `2,5,23`.
These finite computations do not establish an improved asymptotic rate.

The starting operators and crossing identities are the Appendix B
construction in
[`anthropics/formal-math`, revision `16a5ce4`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The definitions are in
[`AppHalfPiecesStatement.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean)
and the integral operator recursion is in
[`Boson/Defs.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Boson/Defs.lean).
The quotient, normalization, cost analysis, and finite-field variants below
are ordinary mathematics with exact finite checks, without new Lean work.

## 1. A quotient preserved by every mode

Let `V_Z=Z[t_1,t_2,...]`, with `t_k` of weight `k`. Fix a prime `p` and set

```math
I_p=(t_{pk}:k\ge1)+(t_k^p:k\ge1),\qquad
W_p=\mathbf F_p[t_1,t_2,\ldots]/I_p.
\tag{1}
```

A basis of `W_p` consists of monomials `t^mu` for partitions `mu` having
no part divisible by `p` and all multiplicities less than `p`. Call these
partitions regular here, and write `q_p(h)` for their number at weight `h`.
This use of "regular" includes both restrictions explicitly.

The source modes are multiplication by `t_k` and `k partial/partial t_k`.
Multiplication preserves every ideal. Each derivation preserves `I_p`:

- Differentiating `t_j^p` gives zero in characteristic `p`.
- Differentiating `t_{pj}` by its matching weighted derivation gives
  the scalar `pj=0`.
- The product rule handles every multiple of these generators.

Consequently every finite normally ordered source operator preserves
`I_p`. Its action descends to `W_p`. For a fixed output weight, the source
mode sums are finite, so this also applies to every coefficient of the
operator series. The ideal is proper and homogeneous with positive
weights; projection to the quotient leaves constant coefficients intact.

The quotient here is in the **state argument** of an operator. We retain
the full polynomial `a^mu` as its operator parameter, including monomials
that belong to `I_p`. Preserving the state ideal alone would not justify
deleting such operator-parameter monomials. Neither the construction nor
its general-family checker makes that deletion.

This is a quotient of the representation on which all source identities
already hold. It does not discard a mode with a nonzero commutator.

## 2. Remove the factorial before reducing modulo `p`

For an ordinary partition `nu`, define

```math
z_\nu=\prod_{k\ge1} k^{m_k(\nu)}m_k(\nu)!.
```

The source's integral polynomial is

```math
a^\mu=\prod_{m\in\mu}
   \left((-1)^{m-1}\sum_{j=1}^m\binom mjt_j\right).
\tag{2}
```

Its stated duality, for partitions of the same weight `h`, is

```math
[1]a^\mu[-h]\left(\frac{h!}{z_\nu}t^\nu\right)
=h![\mu=\nu].
```

All operator coefficients and `a^mu` coefficients are integers.
Multiply this identity by `z_nu` and cancel the nonzero integer `h!`
**over the integers**. The resulting integral identity is

```math
[1]a^\mu[-h]t^\nu=z_\nu[\mu=\nu].
\tag{3}
```

For a regular partition, `z_nu` is a unit modulo `p`. Use the explicit
monomial polynomial below as an operator parameter, and its class in
`W_p` when it is a state:

```math
b^\nu=z_\nu^{-1}t^\nu.
```

Reducing (3) and using the preserved constant coefficient gives

```math
a^\mu[-h]b^\nu=[\mu=\nu]\,1.
\tag{4}
```

There is no cancellation of a zero factorial in `F_p`.

## 3. Exact crossing and its two factorizations

Choose any `m=2^lambda` distinct regular partitions of weight `h` and use
them as the binary labels. For `0<=i<=2h`, let the source angle operation be

```math
a^{\langle i\rangle}
=\sum_{w=0}^h\binom{h+w}{i}a_w,
```

where `a_w` is the weight-`w` part of `a`. For `-2h<=s<=2h`, form

```math
Q_s(a,b)=\sum_{i=0}^{2h}(-1)^i
 a^{\langle i\rangle}[s-h+i]\,b[h-s-i].
\tag{5}
```

The source exchange identity permits the two operators in every summand
to be reversed. These identities have integral matrix coefficients:
their validity over the rationals implies equality of their integer
coefficients, which may then be reduced modulo `p`. The quotient in
Section 1 preserves them.

Let `A_s(x,y,u,v)` be the coefficient of `t^u` in
`Q_s(a^x,b^v)t^y`. The source summation identity and (4) give

```math
\sum_{s=-2h}^{2h}A_s(x,y,u,v)
=[x=v][y=u].
\tag{6}
```

Thus the crossing multiplier is one.

In the exchanged order, the intermediate state has weight `s+i`.
In the original order, its weight is `2h-s-i`. Accordingly the two
factorizations have channel counts

```math
R_s=\sum_{i=0}^{2h}q_p(s+i),\qquad
C_s=\sum_{i=0}^{2h}q_p(2h-s-i),
\tag{7}
```

where `q_p(E)=0` for `E<0`. A zero channel count represents a zero piece;
it may be padded to a single zero channel. The factors are computed by
the same coefficient recursion as the source, with arithmetic in `F_p`
and intermediate coordinates restricted to the quotient basis.

### Extension to coefficient rings

Let `R` be any commutative `F_p`-algebra. The quotient module
`W_p tensor_(F_p) R` is free over `R` on the same basis. Equations
(4)--(7) extend by scalar extension. All factor coefficients remain
scalars from `F_p`; no ring-valued datum is inverted, and no step tests
whether such a datum is nonzero. Zero divisors in `R` cause no problem.

In particular, the construction applies to
`F_p[X]/(X^r-1)` for every positive integer `r`, including `r` divisible
by `p`. If an evaluator is compiled into a fixed arithmetic circuit,
these pieces can be substituted gate by gate over that ring. A complete
factoring algorithm must still account for the circuit's operation count,
ring arithmetic cost, and the recovery or nonzeroness test. Equation (6)
alone is not that runtime theorem.
The [scale-one transfer](scale-one-modular-transfer.md) supplies that
operation and bit-cost argument, using
[sparse divisor nonzeroness](sparse-divisor-nonzeroness.md).

## 4. Labels and factor evaluation are effective

There are enough regular partitions at quadratic energy. Here is an
explicit construction, which also avoids relying on monotonicity of
`q_p(h)`.

For `p=2`, in block `j=0,...,lambda-1`, choose one of

```math
\{8j+1,8j+7\},\qquad \{8j+3,8j+5\}.
```

The choices are disjoint odd parts, each block has sum `16j+8`, and the
`2^lambda` resulting partitions all have weight `8lambda^2`.

For `p>=3`, take the twelve integers `12j+1,...,12j+12` in block `j`.
There are six disjoint symmetric pairs, each summing to `24j+13`.
At most four integers in the block are divisible by `p`, so at least two
pairs avoid multiples of `p`. Choose two such pairs by a fixed rule.
Independent choices in all blocks give `2^lambda` distinct regular
partitions of weight `12lambda^2+lambda`.

We may therefore find the least `h` with `q_p(h)>=2^lambda` by bounded
dynamic programming. Its bound is uniform in `p`. For fixed `p`, all
ordinary partitions through energy `4h` number `2^(O(lambda))`, by the
ordinary partition upper bound. Enumerating regular partitions and their
binary prefixes takes no more time.

The source polynomial (2) can be expanded by multiplying one linear
factor at a time and collecting equal monomials. At every stage, there
are at most `sum_(E<=h) p(E)` monomials. Reduce coefficients modulo `p`
as they are produced. Intermediate integer coefficient sizes are not
part of this computation.

For an individual monomial operator, memoize the normal-order recursion
by its suffix, remaining output partition, and remaining input partition.
There are at most `(h+1)Delta^2` states, where
`Delta=sum_(E<=4h) p(E)`. Each state has polynomially many transitions.
Combining this recursion with the collected expansion costs
`poly(h,Delta)` field operations for a factor. Binomial coefficients
and the units `z_nu^-1` are computed modulo `p` in polynomial time in
`h` for fixed `p`.

Thus the number of pieces is polynomial in `lambda`, factor descriptions
use a fixed number of coefficient bits, and a uniform factor algorithm
takes `2^(O(lambda))` bit operations. All control depends on integer
indices and metadata, not on values from an extension algebra `R`.

## 5. What the state counts prove about contraction

The exact generating function is

```math
Z_p(t)=\sum_{h\ge0}q_p(h)e^{-th}
=\prod_{p\nmid k}\frac{1-e^{-pkt}}{1-e^{-kt}}.
\tag{8}
```

Write `F(t)=-sum_(k>=1) ln(1-e^(-kt))`. Cancelling the factors divisible
by `p` gives

```math
\ln Z_p(t)=F(t)-2F(pt)+F(p^2t),\qquad
A_p=\frac{\pi^2}{6}(1-1/p)^2.
```

Integral comparison with the decreasing function
`-ln(1-e^(-x))`, whose integral is `pi^2/6`, yields

```math
\ln Z_p(t)=A_p/t+O_p(\ln(1/t))\quad(t\downarrow0).
\tag{9}
```

For clarity, no individual-coefficient asymptotic is assumed here.
Give each regular partition probability `e^(-t weight)/Z_p(t)`.
Differentiating the convergent series, or summing the independent bounded
geometric occupations of each part, gives its weight mean and variance.
Direct Riemann-sum estimates give

```math
\mathbb E_t H=A_p/t^2+O_p(1/t),\qquad
\mathrm{Var}_t(H)=O_p(t^{-3}).
\tag{10}
```

Indeed, `-F'(t)=t^-1 sum_k g(kt)`, with
`g(x)=x/(e^x-1)` decreasing, integrable, and having integral `pi^2/6`.
Also `F''(t)=t^-2 sum_k f(kt)=O(t^-3)`, where
`f(x)=x^2 e^x/(e^x-1)^2` is continuous at zero and decreases
exponentially up to a polynomial factor at infinity. Applying these
estimates to the three terms in (9) proves (10).

Chebyshev's inequality puts probability at least `1/2` in the integer
window of radius `t^-7/4` about the mean, for sufficiently small `t`.
Some integer `H` in that window therefore has

```math
\ln q_p(H)\ge 2A_p/t-O_p(t^{-3/4}+\ln(1/t)),\qquad
H=A_p/t^2+o(t^{-2}).
\tag{11}
```

This follows by writing the probability of the window as
`Z_p(t)^-1 sum_H q_p(H)e^(-tH)` and bounding the number of integers
in the window. It is a lower bound for a coefficient somewhere in the
window, which is all the label-selection rule requires.

Conversely, if `Q_p(E)=sum_(0<=j<=E) q_p(j)`, then
`Q_p(E)<=e^(tE) Z_p(t)`. Taking `t=sqrt(A_p/E)` gives the uniform upper
bound

```math
\ln Q_p(E)\le 2\sqrt{A_p E}+O_p(\ln(E+2)).
\tag{12}
```

Let `h(lambda)` be the least weight with `q_p(h)>=2^lambda`.
The upper bound (12) and the coefficient-existence bound (11) imply

```math
2\sqrt{A_p h(\lambda)}
=\lambda\ln2+o_p(\lambda).
\tag{13}
```

For the upper bound on `h`, apply (11) with
`t=2A_p/((1+delta)lambda ln2)` for any fixed `delta>0`, and then let
`delta` tend to zero. For the lower bound, use
`lambda ln2<=ln q_p(h)` in (12); the already proved quadratic upper
bound makes its logarithmic error negligible. This does not assume
that `q_p(h)` increases with `h`.

Equations (7) and (12) imply, uniformly for `|s|<=2h`,

```math
\begin{aligned}
\lceil\log_2\max(1,R_s)\rceil+
\lceil\log_2\max(1,C_s)\rceil
&\le \frac{2\sqrt{A_p}}{\ln2}
       \bigl(\sqrt{2h+s}+\sqrt{2h-s}\bigr)+O_p(\ln(h+2))\\
&\le 2\sqrt2\,\lambda+o_p(\lambda).
\end{aligned}
\tag{14}
```

Consequently, for every fixed `0<epsilon<1-1/sqrt(2)`, these are effective pieces
with width sum at most
`4(1/sqrt(2)+epsilon)lambda+C log2(2lambda)`, after increasing a constant
to cover finitely many small `lambda`.

The generating function alone would not justify an `O(log lambda)`
error in (13). We claim the fixed-epsilon statement above. Neither
(14) nor the smaller entropy constant `2sqrt(A_p)` gives `rho<1/sqrt(2)`:
the label energy changes by the same constant.

## 6. A homogeneous variant when `h` is a power of `p`

Use the integral Gram operator `G_h` from
[`homogeneous-duals.md`](homogeneous-duals.md). Its rational eigenvalues
are

```math
b_j=\binom{2h-1}{h-j},\qquad 1\le j\le h,
```

with the primitive-space multiplicities proved there. Its matrix in the
monomial lattice is integral: both divided-power maps `D^k/k!` and
`U^k/k!` have integral coefficients.

The divided-power maps, not just `D` and `U`, preserve `I_p` modulo `p`.
On generators their generating series are

```math
e^{zD}t_m=\sum_{j=1}^m\binom mj z^{m-j}t_j,\qquad
e^{zU}t_m=\sum_{n\ge0}\binom{m+n-1}{n}z^n t_{m+n}.
\tag{15}
```

If `p|m`, the first coefficient is zero unless `p|j`, by the lowest
base-`p` digit in Lucas's identity. In the second formula, when `p∤n`,
that lowest digit is `binom(n mod p - 1,n mod p)=0`; the surviving terms
again have `p|(m+n)`. Frobenius handles `t_k^p` for both algebra maps.
Thus every coefficient of (15) preserves the ideal, as does `G_h`.

If `h=p^r`, all `b_j` are nonzero modulo `p`. In the lower `r` digits,
`2h-1` has only digits `p-1`, while `0<=h-j<h`; Lucas gives a product
of nonzero numbers. Therefore the determinant of the full integral
`G_h` is a unit modulo `p`. Its reduction is invertible. An invariant
quotient of a finite-dimensional invertible operator is invertible,
so `G_h` is also invertible on `W_(p,h)`.

The digit identity used here follows directly by expanding
`(1+Z)^a=product_j (1+Z^(p^j))^(a_j)` modulo `p`, where `a_j`
are the base-`p` digits. Comparing coefficients gives
`binom(a,b)=product_j binom(a_j,b_j)` modulo `p`.

The diagonal Fock form has entries `z_mu`, all units on the regular
basis. Its radical on the full characteristic-`p` monomial space is
exactly the span of the excluded monomials. The integral Gram identity
therefore descends and makes the source homogeneous pairing on regular
operator labels and regular states nonsingular. Invert its matrix over
`F_p` to obtain explicit homogeneous polynomials `a^mu` with

```math
[1]a^\mu[-h]t^\nu=[\mu=\nu].
\tag{16}
```

Ordinary elimination has polynomial cost in `q_p(h)` and uses only
`log2 p` coefficient bits. There is no exponentially long rational
denominator. For any selected subset of labels, rows of the full inverse
provide suitable duals. One may instead invert a nonsingular selected
submatrix; that is a different finite choice of dual family.

For this homogeneous family,
`a^(angle i)=binom(2h,i)a`. When `p=2` and `h` is a power of two, only
`i=0,2h` survive. For odd `p` and `h=p^r`, only `i=0,h,2h` survive, with
coefficients `1,2,1`. This substantially simplifies exact experiments.
Choosing only prime-power energies can, however, overshoot the energy
needed for a particular binary label width. Section 5's width statement
belongs to the general normalized family; it is not a claim about this
sparse choice of energies.

## 7. Why the one-bit rank saving does not tensor into a better rate

Over `F_2`, the one-bit crossing has an elementary two-piece identity.
Use the row basis `(0,0),(0,1),(1,0),(1,1)` and put

```math
v=(1,0,0,1)^T,\qquad R^{\rm row}=vv^T,\qquad C=E+R.
```

Here `E` is the crossing, whose row and column matrices are both `I_4`.
The column matrix of `C` is `ww^T`, where `w=(0,1,1,0)^T`.
Thus `R,C` have row/column ranks `(1,4)` and `(4,1)`, respectively,
and `R+C=E`. Each individual piece has width sum two at label width one.

Tensoring the individual choices creates `2^t` pieces at label width `t`.
The usual polynomial-piece alternative is to take the `t+1` coefficients

```math
A_i=[z^i](R+zC)^{\otimes t},\qquad 0\le i\le t.
\tag{17}
```

Their sum is `E^(tensor t)`, but their ranks grow more quickly than the
one-bit rank products suggest. Indeed, `v^Tv=w^Tw=0` in `F_2`.
Each cheap local matrix is rank one and square zero, hence is similar to
`J_2 direct-sum 0 direct-sum 0`, with `J_2` the nonzero nilpotent
two-by-two Jordan block.

Let `R_j` act as the local row matrix on site `j` and as the identity
elsewhere. Write `e_k` for the elementary symmetric polynomial of degree
`k` in these commuting operators. Since `R_j^2=0`,

```math
A_i^{\rm row}
=e_{t-i}(R_1,\ldots,R_t)\prod_{j=1}^t(I+R_j).
\tag{18}
```

The product on the right is invertible. This identity follows directly
from `R_j+z(I+R_j)=(I+R_j)(R_j+zI)`.
In the column flattening the cheap and expensive factors exchange roles;
its rank is the rank of `e_i` for the same local nilpotent type.

Split each local four-dimensional space into its active two-dimensional
`J_2` space and its inactive two-dimensional zero space. Tensoring these
direct sums gives invariant sectors. On a sector with exactly `k` active
sites, only the product of all `k` active operators survives in `e_k`.
It has rank one on the active factors and acts as the identity on the
inactive factors. Summing the independent sectors gives the exact lower
bound

```math
\mathrm{rank}e_k(R_1,\ldots,R_t)
\ge\binom tk\,2^{t-k}.
\tag{19}
```

For even `t` and `i=t/2`, both flattenings of `A_i` therefore have rank
at least `binom(t,t/2)2^(t/2)`. Their product is at least
`2^(3t-O(log t))`. The width-sum condition at `lambda=t` consequently
requires `rho>=3/4`, already larger than `1/sqrt(2)`.

This obstruction applies to the coefficient family (17). It does not
claim a lower bound for every decomposition over `F_2`. It does rule out
using the two one-bit ranks alone to claim `rho=1/2` for this
polynomial-piece tensor symmetrization.

## 8. Exact checks and the remaining target

Run

```sh
python3 tools/check_modular_crossing.py --output evidence/modular-crossing.json
```

The checker verifies the integral cancellation identity, preservation of
the quotient in finite operator tables, regular-partition counts, the
explicit quadratic-energy label families, normalized duality, both
orders of the crossing pieces, and their scale-one sum. It also checks
the homogeneous prime-power family and records the finite ranks above.

For `p=2`, the checker uses an independent determinant expression for a
normal-ordered coefficient. Regular monomials have distinct odd parts.
For operator parts `a`, input `sigma`, and output `pi`, choose the common
unaffected subset `eta` of size
`(length(pi)+length(sigma)-length(a))/2`. The remaining output parts
are creations; the remaining input parts are annihilations. Sum over
all such `eta`. For each choice, assignments of operator parts to these
signed modes form a permanent, which equals a determinant in `F_2`.
The matrix entries for operator part `m` and mode `k` are

```math
\binom{k-1}{m-1}\bmod2\quad\text{(creation)},\qquad
\binom{k+m-1}{m-1}\bmod2\quad\text{(annihilation)}.
```

Lucas's bit tests evaluate both entries. The checker compares this
formula with the source's integer recursion before using it for larger
finite rank calculations.

The concrete next target is a proved joint rank-width profile for the
**complete** normalized pieces that improves (14). Sparse angle support
and a smaller ambient state space do not prove that target. The scale-one
ring extension removes an arithmetic obstruction to pursuing it.

The underlying crossing construction is prior work by Anthropic. This
analysis and checker used substantial OpenAI assistant support.
