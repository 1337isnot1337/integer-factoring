# Implicit inverse moments with compressed analytic weights

## Result and scope

The [polynomial inverse-moment algorithm](pole-cauchy-polynomial-moments.md)
extends to weights represented by arithmetic circuits. These may contain
enormous binary-encoded powers and reciprocals that are odd on the odd
input class. Symbolic polynomial degree is not charged as expanded input.
For two equal-length odd progressions, the complete scalar
`g(B)^T M^-1 f(A)` is computable to specified absolute 2-adic precision
in polynomial time in the input lengths, circuit size, exponent encoding
lengths, and charged precision. Neither matrix axis is enumerated.
The singular factor-count trace remains open; the score remains one half.

## 1. Matrix and conditioning contract

Let `N` be positive and odd. Let `A,B` be positive odd progressions of
step two and equal length `m`. Choose `k>=1` so every nonzero `N-ab`
has absolute value less than `2^k`. In the unramified quadratic 2-adic
extension with `w^2+w+1=0`, set

\[
 t=2^kw,\quad z=N+t,\quad M_{a,b}=\frac{t}{z-ab}.
 \tag{1}
\]

The [balanced conditioning theorem](pole-cauchy-conditioning.md) applies:
at each odd residue level, both progressions have only the two possible
floor/ceiling counts. Permuting column residues by `b -> N/b` preserves
that property, so row and reciprocal-column counts differ by at most one.
In particular, every entry of

\[
 U=2^kM^{-1}                                        \tag{2}
\]

is integral. The stronger loss bound `k-1` is not needed. The computable
progression contract matters: this argument is not applied to arbitrary
unequal residue distributions.

If integral functions `f,F,g,G` agree pointwise modulo `2^P` on odd
integers, then

\[
 g(B)^{\mathsf T}Uf(A)
 \equiv G(B)^{\mathsf T}UF(A)\pmod{2^P}.             \tag{3}
\]

Expand the difference entry by entry and use (2). Summing any finite
number of multiples of `2^P` loses no 2-adic precision.

## 2. Regular circuit gates

The [compiler](../tools/pole_cauchy_analytic_weights.py) accepts a
topologically ordered list of gates.

| Gate | Meaning | Condition |
|---|---|---|
| `var` | Odd argument `x` | None |
| `const,c` | Encoded integer | None |
| `add,i,j` / `sub,i,j` | Sum or difference of earlier gates | None |
| `mul,i,j` | Product of earlier gates | None |
| `pow,i,e` | Binary-encoded integer power | Negative powers require a unit |
| `inverse,i` | Reciprocal | Value at `x=1` must be odd |

Write `x=1+2X`. Every accepted gate has the convergent expansion

\[
 F(1+2X)=\sum_{i\ge0}c_iX^i,\qquad c_i\in2^i\mathbb Z_2.
 \tag{4}
\]

The input and constants have this property; sums and products preserve it.
For a unit inverse the coefficients satisfy

\[
 b_0=a_0^{-1},\qquad
 b_i=-a_0^{-1}\sum_{j=1}^{i}a_jb_{i-j},              \tag{5}
\]

so induction proves the same divisibility. Its input is odd at every odd
integer because all nonconstant terms in (4) are even. Powers use repeated
squaring, applying (5) first for negative exponents. For example,
`1/(x^2+x+1)` is accepted and `1/(x-1)` is rejected.

## 3. Precision P requires degree at most P-1

Compute all gates in

\[
 (\mathbb Z/2^P\mathbb Z)[X]/(X^P).                 \tag{6}
\]

Every discarded term has valuation at least `P` on integer `X`, by (4).
Choose the representatives `0<=c_i<2^P` of the exact retained prefix.
They are divisible by `2^i`, and

\[
 F_P(x)=\sum_{i=0}^{P-1}\frac{c_i}{2^i}(x-1)^i      \tag{7}
\]

is an integer polynomial of degree at most `P-1` agreeing with the
circuit modulo `2^P` at every odd integer. Expanding (7) supplies the
monomial coefficients for the inverse-moment algorithm.

Dividing `c_i` by `2^i` retains only `P-i` coefficient bits, but
`(x-1)^i` contributes at least `i` powers of two. A different coefficient
lift changes the value by a multiple of `2^P`. The proof uses pointwise
agreement and (2), not an incorrect claim of `P` independent bits in each
divided coefficient. Trimming zero high coefficients can reduce the
compiled degree far below the worst-case bound.

## 4. Whole bilinear forms and charged denominators

Compile both weights and apply the existing two-jet evaluator to `F_P,G_P`.
Equation (3) proves the result for the entire bilinear form. External
scalars `2^-r,2^-s` may multiply the regular weights. The returned integral
residue `R` then satisfies

\[
 g(B)^{\mathsf T}M^{-1}f(A)
 =2^{-k-r-s}R+O(2^{P-k-r-s}).                       \tag{8}
\]

For absolute precision `p`, use `P=p+k+r+s`. Denominator exponent values
are charged; no enormous exact division is free. The scale is not a
claim of exact valuation after possible cancellation.

## 5. Bit complexity

Let `S` count ordinary gates and `E` be the sum of the binary encoding
lengths of all power exponents. Naive truncated arithmetic uses

\[
 O((S+E+1)P^2)                                     \tag{9}
\]

operations on polynomially many-bit integers, including inverse recurrences
and basis conversion. Modular coefficients have at most `P` bits; integer
constant bit lengths are charged as input. The polynomial inverse-moment
algorithm has polynomial bit cost in its degrees, precision, endpoint
lengths, and `log m`. Its supplied degrees are at most `P-1`, so the
combined algorithm has polynomial bit cost in these parameters and `S+E`.
It uses the stated circuit descriptions, not free access to arbitrary
exponentially long vectors. Arithmetic dependencies are those of
[modular interval jets](modular-interval-jet-arithmetic.md).

## 6. Exact controls and a large implicit example

The [checker](../tools/check_pole_cauchy_analytic_weights.py) compiles 40
circuits at ten precisions and checks all 40,920 odd residue values. It
compares 80 complete bilinear forms with dense exact inversion on 40
matrices, verifies all 512 normalized inverse entries are integral, and
checks external dyadic scales. Five invalid contracts are rejected.

The large control uses the 402-bit integer `N=2^401+45` and

\[
 m=2^{200}+13,\quad
 A=\{1,3,\ldots,2m-1\},\quad B=\{3,5,\ldots,2m+1\},
\]

with weights

\[
 f(x)=(1+2^{64}x)^{2^{512}+2^{257}+123},\qquad
 g(x)=(1+2^{80}x)^{-(2^{383}+2^{193}+7)}.             \tag{10}
\]

The power encodings have 513 and 384 bits. At the precision required for
twelve absolute output bits, their compiled degrees are six and five.
The checker computes the entire bilinear form at twelve and sixteen
absolute bits and checks agreement, enumerating zero matrix nodes.
The [evidence](../evidence/pole-cauchy-analytic-weights.json) records the
circuits, intervals, coefficients, precision budgets, residues, operation
counts, and source hashes.

```sh
python3 tools/check_pole_cauchy_analytic_weights.py \
  --output evidence/pole-cauchy-analytic-weights.json
```

The factor-count kernel has singular dependence on the product `ab`.
It is not a pair of separate regular odd-input circuits covered here.
That aggregation problem remains open.

