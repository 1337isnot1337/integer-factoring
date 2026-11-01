# Polynomial-circuit inverse moments and a sharp denominator barrier

## The constructive extension

Let `A,B` be equally sized sets of distinct positive odd integers, each
given as an explicitly listed union of disjoint odd progressions.
Let their common cardinality be `m`. For an odd positive integer `N`
and `k>=1`, put

```math
z=N+2^k w,\qquad t=2^k w,\qquad w^2+w+1=0,\qquad
M(a,b)=\frac{t}{z-ab}.
```

This Cauchy matrix is invertible: its nodes are distinct and its cross
denominators cannot vanish because `z` has a nonzero `w` coordinate.

Given a division-free arithmetic circuit for an integer polynomial
`h(X,Y)` and precision `P>=1`, there is a deterministic algorithm for

```math
2^k\sum_{a\in A,b\in B}(M^{-1})_{b,a}h(a,b)\pmod {2^P}
```

with bit cost polynomial in the circuit's written length, the binary
node-domain description, `P`, and the binary input parameters. The
degree of `h` need not be explicitly bounded or expanded. The algorithm
reduces the circuit in a bivariate quotient of dimension at most
`(P-1)^2`.

This extends the
[bounded-degree polynomial inverse-moment interface](pole-cauchy-polynomial-moments.md).
It is still a theorem about polynomial weights. The actual pole-square
kernel has nonunit denominators and is outside this interface. A sharp
finite-difference argument below shows that polynomial representations
of that kernel can require exponentially many denominator bits even
when polynomial degree and circuit size are unrestricted.

## 1. Cancel a common modular packet factor

Define monic polynomials

```math
P_B(Y)=\prod_{b\in B}(Y-b),\qquad
P_X(Y)=\prod_{a\in A}(Y-z/a).
```

Every root is congruent to one modulo two in the quadratic valuation
ring. Write `J=min(m,P-1)`. Expanding around `Y-1`, every term that uses
`j` root corrections has a coefficient divisible by `2^j`. Therefore

```math
P_B(Y)\equiv(Y-1)^{m-J}F(Y),\qquad
P_X(Y)\equiv(Y-1)^{m-J}G(Y)\pmod {2^P},
```

where `F,G` are monic polynomials of degree `J`. Both constant terms
are units. The products and `(Y-1)^(m-J)` have unit constant terms, so
the quotient identities give

```math
R(Y)=P_B(Y)/P_X(Y)\equiv F(Y)/G(Y)
```

as a formal power series at zero. If
`F^*(U)=U^J F(1/U)` and `G^*(U)=U^J G(1/U)`, the expansion at infinity
likewise becomes

```math
T(U)=P_X(1/U)/P_B(1/U)\equiv G^*(U)/F^*(U).
```

Both identities retain every coefficient of these infinite series
modulo `2^P`, not only their first `J` coefficients.

The packets can be constructed using only degree-`J`
[unit jets](modular-interval-jet-arithmetic.md). Compute
the normalized products `prod_b(1-Y/b)` and `prod_a(1-aY/z)` using the
existing unit-progression algorithm. Divide each by `(1-Y)^(m-J)` in
the truncated jet ring. Each result is known to be a polynomial of
degree at most `J` modulo `2^P`, by the preceding expansion. Its leading
coefficient is a unit; divide by it to obtain `F` or `G` monic. This
charges all constants, normalization, and binary powering without
enumerating either domain.

## 2. The whole moment matrix has a small Bezout numerator

Write `R(Y)=sum r_j Y^j` and `T(U)=sum t_j U^j`. The
[existing exact interpolation formula](pole-cauchy-polynomial-moments.md) says

```math
\mu_{e,d}=\sum_{a,b}b^e(M^{-1})_{b,a}a^d
 =\frac{z^{d+1}}{t}K_{e,d},\qquad
K_{e,d}=\mathbf1_{e=d}-\sum_{j=0}^{\min(d,e)}r_{d-j}t_{e-j}.
```

Hence

```math
\sum_{e,d\ge0}K_{e,d}U^eV^d
 =\frac{1-T(U)R(V)}{1-UV}
 =\frac{G(V)F^*(U)-F(V)G^*(U)}
        {(1-UV)F^*(U)G(V)}.
```

The polynomial numerator is divisible by `1-UV`. Its quotient `H(U,V)`
has degree less than `J` in each variable. This follows either by
polynomial division or by the explicit telescoping coefficient formula
for the two products. It is an identity over every commutative
coefficient ring; no field inversion is used in this division.

Let `L_(e,d)=2^k mu_(e,d)`. Incorporating the factor `z^(d+1)/w` gives

```math
\sum_{e,d\ge0}L_{e,d}U^eV^d
 \equiv\frac zw\,
 \frac{H(U,zV)}{F^*(U)G(zV)}\pmod {2^P}.
```

The numerator has degree less than `J` in both variables. Both
denominator constant terms are units. Therefore every row and column
of the normalized moment array satisfies homogeneous recurrences
starting at index `J`.

## 3. A finite polynomial quotient computes every circuit weight

The two monic recurrence polynomials are

```math
Q_B(Y)=F(Y),\qquad
Q_A(X)=\frac{X^J G(z/X)}{G(0)}.
```

Let `Lambda(h)=2^k sum_(a,b) (M^(-1))_(b,a) h(a,b)` reduced modulo
`2^P`. The recurrence identities state that

```math
\Lambda(Q_A(X)f(X,Y))=\Lambda(Q_B(Y)f(X,Y))=0
```

for every polynomial `f` over the finite quadratic coefficient ring.
Thus `Lambda` factors through

```math
(\mathbb{Z}/2^P)[w,X,Y]/(w^2+w+1,Q_A(X),Q_B(Y)).
```

Both quotient polynomials are monic. Every class has a representative
with degrees less than `J` in each variable, computed by ordinary
polynomial reduction. Evaluate the supplied arithmetic circuit directly
in this quotient. All intermediate representatives have at most `J^2`
coefficients. The result is contracted with the first `J` by `J`
normalized moment matrix, which is obtained from the existing formula
using only the first `J` origin and infinity coefficients.

For `J=0`, the generating numerator is zero and every normalized
polynomial moment is zero modulo `2^P`; the algorithm returns that
residue directly.

Naive bivariate multiplication uses `O(J^4)` coefficient operations per
node, with polynomial-time reductions. The coefficient ring has two
coordinates of `P` bits. Context construction, the initial moment
matrix, and every circuit operation therefore have polynomial bit
cost. Binary powering can be written as ordinary circuit nodes, so
monomial exponents given in binary are included. No exponentially long
coefficient list is created.

The output represents the unscaled scalar as `2^(-k)` times its residue,
with absolute error valuation at least `P-k`. For `q` absolute bits,
choose `P=q+k`. In the factoring pole setup, `k=O(n+q)`, so this is still
polynomial precision.

These quotient polynomials annihilate the polynomial functional; they
are not asserted to vanish at every original node. In particular, one
must not invert the nonunit `z-XY` in this quotient and treat that as
an evaluation of the actual pole-square kernel.

## 4. A sharp exponential denominator requirement for a factor column

Suppose a column contains a true factor pair, and its row domain is
`H=2^s` consecutive odd nodes `a_i=a_0+2i`. Choose the pole displacement
so that every nonzero `N-a_i b` has valuation below `k`. Then

```math
M(a_i,b)^2\equiv\mathbf1_{i=i_0}\pmod2
```

for the unique row `a_(i_0)=N/b`.

Let `W(X,Y)` be a polynomial of arbitrary degree over the quadratic
2-adic field with the same values modulo two on this column. Suppose
`2^L W` has integral coefficients, where `L` is a nonnegative integer.
Then necessarily

```math
L\ge(H-1)+v_2((H-1)!)=2H-2-s.
```

To prove this, take the `(H-1)`st forward difference with step two of
the univariate polynomial `W(X,b)`. Its value modulo two is
`binom(H-1,i_0)`, up to sign, hence is odd. Indeed every coefficient of
`(1+X)^(2^s-1)` is odd. For an integral polynomial, the `q`th step-two
difference is divisible by `2^q q!`: expand a monomial and use the
identity `Delta^q(j -> j^r)(0)=q! S(r,q)`. Applying this to
`2^L W(X,b)` forces the stated lower bound with `q=H-1`. Finally
`v_2(q!)=q-s` for `q=2^s-1`.

The bound is attained by the degree-`H-1` cardinal Lagrange polynomial
for row `i_0`. Its denominator is

```math
2^{H-1}(-1)^{H-1-i_0}i_0!(H-1-i_0)!.
```

Because `binom(H-1,i_0)` is odd, this denominator has exactly the
valuation in the bound. Thus allowing arbitrary degree or a succinct
polynomial circuit does not reduce the necessary denominator precision.

When `H` is exponential in the binary input length, so is this lower
bound. The result is an obstruction to polynomial-kernel approximation
with a small coefficient denominator budget. It is not a factoring
lower bound and does not exclude a direct normalized rational algorithm.
It also does not exclude a different polynomial weight chosen only to
match the final scalar contraction, rather than the kernel entrywise.

## 5. Exact checks and implementation

The [implementation](../tools/pole_cauchy_circuit_moments.py) accepts ordinary
`row`, `column`, `const`, `add`, and `mul` circuit nodes. Its context
constructor records the two quotient polynomials and the small normalized
moment matrix. All powers in the supplied examples are expanded into
addition/multiplication nodes; no exponentiation is assigned unit cost.

Run:

```sh
python3 tools/check_pole_cauchy_circuit_moments.py --output /tmp/pole-cauchy-circuit-moments.json
cmp /tmp/pole-cauchy-circuit-moments.json evidence/pole-cauchy-circuit-moments.json
```

The [checker](../tools/check_pole_cauchy_circuit_moments.py) computes 48
independent exact inverses over `Q(w)` and compares 144 circuit-weighted
contractions, including 81-bit monomial exponents and coupled high powers.
Twenty-four matrices use a quotient strictly smaller than the matrix.
Another 16-dimensional control uses precision 14 and displacement 10,
so its result includes four bits of absolute precision.

The long control has `2^100+3` rows, a 202-bit `N`, and row/column
exponents of 129 and 145 bits. Its 284-node circuit is evaluated in a
15-by-15 quotient. Complete odd residue cycles cancel from the origin
and infinity ratios because the displacement is above the working
precision: inversion followed by multiplication by `N` permutes the odd
residues. The remaining three rows and columns give an independent
literal `3`-by-`3` inverse comparison. The long domain is not enumerated.

The denominator checks verify 341 literal cardinal-polynomial values,
128 step-two finite differences, and all sharp denominator valuations
for power-of-two row counts through 512. The formula at `H=2^200` is
also evaluated from the binary count without constructing a grid or
factorial of that length. The [receipt](../evidence/pole-cauchy-circuit-moments.json)
records the parameters, residues, and hashes.

## Status

The complete factoring score remains `1/2`; no new Lean formalization
is included. The construction removes the degree-expansion cost for
polynomial weights in the inverse-moment interface. The denominator
theorem excludes direct entrywise polynomial approximation of the
pole-square kernel at polynomial coefficient precision.
