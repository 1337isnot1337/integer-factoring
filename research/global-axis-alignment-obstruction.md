# A global obstruction to aligning all low axis flags

## Question and result

Base: `8e697ab`, with PR65/67 integration separately underway. The last
three explicit high-axis kernels fail through low-prefix minors. Instead
of changing another coefficient, ask whether all axes can share a common
linear factor in their low exterior jets, which would make the higher
channels lie inside the same flag. The following homogeneous divisibility
argument rules out that alignment through a precise depth, across the
entire coefficient class, and gives a Smith-prefix bound on the quotient
by the proposed linear flag. A full-axis pressure consequence remains open.

Let k have characteristic two. For every nonempty subset T of `[d]`, let
`c_T` be a homogeneous polynomial of degree D, depending only on variables
`x_i` with `i in T`, and satisfying `c_T(1,...,1)=1`. No permutation symmetry
is needed. At axis a define, over the rational function field in the other
variables, the exterior forms

```
F_a = sum_(B subset [d]\{a}) c_({a} union B) epsilon_B,
U_a = sum_(b!=a) x_b^u epsilon_b,
```

where u is a positive integer. Take any odd set S of size h with
`h>=5` and `floor(D/(h-1))<u`. Then **some a in S** has

```
val_(x_a) [epsilon_(S\{a})](U_a F_a) <= floor(D/(h-1)).   (1)
```

Thus the low axis forms cannot all be divisible by `U_a` modulo
`x_a^r` when `floor(D/(h-1))<r<=u`. For a permutation-symmetric
coefficient family the same degree-specific failure occurs at every axis.

For even D and `u=D/2`, odd set sizes `5,7,9` force failures by
orders `floor(D/4),floor(D/6),floor(D/8)`, respectively. At D=64 these
are **16,10,8**. These conditions concern complete homogeneous polynomial
coefficients, rather than a finite set of field specializations.

## Proof by a single nonzero homogeneous polynomial

For the chosen odd S form

```
P_S = sum_(b in S) x_b^u c_(S\{b}).
```

This has degree `D+u`, and `P_S(1)=h mod 2=1`, so it is nonzero.
Every surviving monomial in it has at least one exponent exactly u:
in its originating summand the coefficient does not depend on x_b.
Choose a surviving monomial and such an index b. The other h-1 exponents
sum to D, so at least one, indexed a, is at most `floor(D/(h-1))<u`.
Consequently `val_(x_a) P_S <= floor(D/(h-1))`.

But the coefficient on the left of (1) is exactly

```
Q_(a,S) = sum_(b in S\{a}) x_b^u c_(S\{b})
        = P_S + x_a^u c_(S\{a}).
```

The added term is divisible by x_a^u and cannot cancel the chosen lower
order. This proves (1) over every characteristic-two field.

Equivalently, if all Q_(a,S) were divisible by x_a^r with r<=u, then
P_S would be divisible by the product of those x_a^r. The subset-local
monomial with one exponent exactly u would require
`D >= (h-1)r`, giving the same contradiction. Locality is essential to
this sharper bound; homogeneity alone only gives `hr<=D+u`.

Over the other-variable rational function field, a nonzero coefficient
of U_a is a unit. An invertible linear change of exterior generators
therefore sends U_a to one generator. Its annihilator equals its principal
ideal, even modulo x_a^r. Hence vanishing of U_a F_a is exactly the
claimed divisibility condition on F_a; no analytic limiting argument is
needed. Permuting axes proves the symmetric corollary.

## Relation to the linear flag and the next bound question

There is also a module consequence. Let `Abar_a` be multiplication by
F_a on the exterior quotient by `(U_a)`, of dimension `2^(d-2)`.
For the axis and set supplied by (1), put `m=2^(d-h)`. Then

```
Delta_m(Abar_a) <= m*floor(D/(h-1)).                    (2)
```

Indeed, the multiplication matrix of `U_a F_a` has a coordinate minor
of size m with determinant `Q_(a,S)^m`: use columns indexed by subsets
outside `S\{a}` and rows obtained by adjoining all of `S\{a}`. The
matrix is triangular with that coefficient on every diagonal entry.
An invertible constant change sending U_a to the first generator sends
multiplication by U_a F_a to the quotient multiplication matrix plus
zero blocks. Its nonzero determinantal-prefix orders are therefore those
of Abar_a. This proves (2) without estimating uncomputed full-axis orders.

For d=9, D=64 and u=32, the respective quotient prefixes have sizes
`16,4,1` and order sums bounded by **256,40,8**. These are uniform
consequences for the coefficient class, with an axis supplied for each
odd set size. In the symmetric class they hold at every axis.

For a permutation-symmetric two-letter coefficient of degree D, the first
axis term is a scalar times `x_a^e x_b^(D-e)`. Its linear exterior flag
is therefore of the power-linear form `U_a` with `u=D-e`. The theorem
shows why making the higher low jets all lie in this one flag is globally
incompatible with the endpoint. It applies irrespective of the higher
coefficient choices and without increasing the field or degree.

The quotient bound does **not** by itself show that these transverse coefficients
add enough independent full-axis channels. The counterexample below rules
out a lift based solely on those bounds. A successful full-prefix argument
must also use the complete crossing endpoint and its quadratic/higher-flag
compatibility. No universal high-axis pressure floor or stronger factoring
exponent is claimed. No further coefficient sampling is registered. The
score remains `1/2`.

The later [generic face-rank lemma](multiaxis-face-rank-floor.md) gives
a complementary global restriction: a symmetric kernel with any pair or
triple coefficient surviving on a coordinate face already fails the
half-score criterion in every dimension at least four. Setting all other
variables to one can annihilate these generic face minors. That lemma
does not settle the strictly positive low jets of the earlier degree-64
examples or lift the quotient bounds in (2).

## The quotient bounds alone do not lift to the required prefix

There is an exact exterior-module counterexample to that proposed
implication, so a further argument must retain more endpoint information.
On eight exterior generators over F2 consider multiplication by

```
F(z)=z^64 + z^32 epsilon_0 + z^16 epsilon_1 epsilon_2 epsilon_3
     +z^10 epsilon_1 ... epsilon_5 +z^8 epsilon_1 ... epsilon_7.
```

It has entry degree at most 64 and `F^2=z^128 I`. Quotient by
`epsilon_0`. To compute that quotient's complete profile, first set

```
B=z^16 +z^10 epsilon_4 epsilon_5
        +z^8 epsilon_4 epsilon_5 epsilon_6 epsilon_7.
```

Its sixteen-dimensional multiplication matrix has complete orders

```
b=(8,10,10,12,16^8,20,22,22,24),       sum b=256.
```

This follows by splitting the top-monomial action on each two-variable
factor into its single Jordan pair and two homology vectors: the inner
four-dimensional matrix `z^10+z^8 epsilon_6 epsilon_7` has orders
`(8,10,10,12)`, and the next split adds the complementary orders and eight
copies of 16. The quotient by epsilon_0 consequently has orders b,
96 copies of 64, and the complements `128-b`.
It attains **all three** quotient-prefix bounds `8,40,256` exactly.

The full eight-variable F has low half consisting of two copies of b
and 96 copies of 32, with the high half their complements to 128.
For a direct verification, Smith-reduce the quotient matrix A. Since
`A^2=z^128 I`, the block matrix

```
[[A,0],[z^32 I,A]]
```

reduces to two-by-two blocks with diagonal orders `a_i,128-a_i` and
off-diagonal order 32. Their low order is `min(a_i,128-a_i,32)`.
Thus the first 86 full orders sum to

```
2*256 +54*32 =2240,       3*2240 >64*86.
```

Prefix averages are nondecreasing. No prefix of more than `256/3`
channels satisfies the registered Jensen criterion. This is an algebraic
counterexample to inferring that criterion from the new quotient bounds,
the square identity, degree bound and exterior multiplication alone.

Its endpoint is **not** the required crossing coefficient table:
`F(1)=1+epsilon_0+epsilon_123+epsilon_12345+epsilon_1234567`.
It is therefore not a candidate crossing. The missing constraint is the
full endpoint and its compatibility across the remaining quadratic/higher
flags. Proving a global rate floor requires that additional information;
the divisibility and quotient theorems above remain valid.

## Exact finite controls

The [checker](../tools/check_global_axis_alignment.py) retains 96 arbitrary
homogeneous F4 coefficient tables at odd set sizes 5,7,9, in addition to
the degree-64 Frobenius fixtures. For each table it reconstructs P and all
Q, checks homogeneity/locality/endpoint, and supplies an explicit surviving
monomial, distinguished exponent and offending axis. The symmetric fixtures
attain each degree-specific valuation bound `16,10,8`; this is sharpness
of those coefficient bounds, not a claim about the whole axis spectrum.

A nonbinary fixture also performs the constant exterior coordinate change
`new epsilon_0=U`. It verifies the quotient identity coefficientwise and
a triangular sixteen-by-sixteen minor with determinant `z^256`.
For the module counterexample, an untruncated sixteen-dimensional profile
checks b, and the derived 128/256-dimensional complete profiles match the
native calculations. The algebraic block reduction above supplies the full
profile argument; the large native arrays are not its sole justification.

The [lossless compressed receipt](../evidence/global-axis-alignment-controls.json.xz) preserves
all polynomial tables and monomial witnesses, the explicit quotient minor,
complete counterexample profiles and source hashes. Run
`make verify-global-axis-alignment`. These checks do not supply a crossing
for the counterexample or a full-axis pressure theorem.

Substantial OpenAI assistant support. The preceding Cayley and exterior-flag
framework retains credit; no Lean or independent mathematical review.
