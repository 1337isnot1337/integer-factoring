# An exact rational pressure bound for a stabilized hot code

## Result

The stabilized twenty-four-letter physical code below supplies a uniform
balanced crossing family with parameter strictly below

```math
\rho=\frac57=0.714285714\ldots.                         \tag{1}
```

The pressure inequality is an exact integer comparison. The physical
type is also rational and explicitly recorded. This certifies an upper
bound for this particular finite-alphabet construction; it does not
claim that a numerical minimum is optimal.

The value in (1) is still larger than `1/sqrt(2)`. This checkpoint does
not improve the established factoring exponent one half. It certifies
a stronger explicit small-alphabet input for the ongoing crossing search.

## 1. The finite tensor and its complete profiles

Use interval `M=16` and retain the first twenty-four masks in increasing
order of `(natural energy, integer mask)`. This includes every mask of
energy at most 56. In the hot coordinates of
[stable-hot-reflection-codes.md](stable-hot-reflection-codes.md), each
occupied mode `a` has energy `8a+4`.

All energies are strictly below `4M=64`. The primitive tensor therefore
stabilizes under enlargement of the reflection interval, by the exact
coefficient identity proved there. Divide its complete scalar content
before computing any cut profile. The resulting polynomial degree is
112, and its endpoint is the literal twenty-four-letter crossing.

At zero, put physical height `E(x)/2` on each of the first two legs.
At infinity, put physical height `E(u)/2,E(v)/2` on the last two legs.
The infinity profile uses the raw reciprocal Laurent tensor, retaining
its negative orders. Let the complete resulting orders be

```math
\nu_1,\ldots,\nu_{576},\qquad
\mu_1,\ldots,\mu_{576}.                               \tag{2}
```

All heights and orders are integers. The checker recomputes both
complete profiles through the rectangular tensor interface and through
independently assembled square matrices. Every promised channel must
resolve; a missing pivot is an error. The strict crossing endpoint
proves the promised generic rank 576 at each cut.

The exact profiles, coefficient hash, scalar content, and all physical
labels are retained in the
[certificate record](../evidence/stable-hot-code-rational-pressure.json).

## 2. One positive rational comparison

Set

```math
z=\frac{54}{55},\qquad A=10,\quad B=7,\quad
R=\frac AB=\frac{10}{7},\quad t=A\log(1/z).             \tag{3}
```

Define positive real sums, each of which is rational here,

```math
Z_0=\sum_i z^{A\nu_i},\qquad
Z_\infty=\sum_i z^{A\mu_i},\qquad
Q=\sum_{x\in C}z^{B E(x)}.                            \tag{4}
```

The certificate verifies the strict inequality

```math
(Z_0Z_\infty)^7<Q^{20}.                               \tag{5}
```

There is no floating exponential or logarithm in this test. To evaluate
a channel sum, subtract its minimum order, clear the denominator
`55^(A(max-min))`, and sum the corresponding positive integer
monomials. The two minimum orders sum to zero, so their powers cancel
in the product. Clear the remaining powers of 55 from (5) and compare
the two positive integers exactly. The evidence stores their bit
lengths and hashes, together with a hash of their positive difference.

## 3. A fixed rational type and a strict entropy margin

Choose the same physical type on all four legs,

```math
p_x=\frac{z^{B E(x)}}Q.                               \tag{6}
```

These probabilities are positive rationals. Integer type counts are
obtained by clearing the common denominator in (6); the record contains
all of them and their sum. A tensor length divisible by that sum
realizes the type exactly. This potentially large denominator is a
fixed constant for this fixed construction.

Let `H=-sum_x p_x log p_x` and `E_p=sum_x p_x E(x)`. Then

```math
H=\log Q+B\log(1/z)E_p
  =\log Q+\frac tR E_p>0.                             \tag{7}
```

Taking logarithms of the verified positive rational inequality (5)
gives `log Z_0+log Z_infty<2R log Q`. The complete coefficient cost,
including both physical means, is therefore

```math
K=\log Z_0+\log Z_\infty+2tE_p<2RH,
\qquad \frac{K}{4H}<\frac R2=\frac57.                  \tag{8}
```

Both physical means are essential. At zero the two active legs each
contribute `E_p/2`; the infinity gauge has the same total. Dropping
either total would not certify this entropy-normalized inequality.

## 4. From the pressure margin to an effective crossing family

For a multiple `m` of the type denominator, restrict the `m`th tensor
power to the exact type class `T_p(m)` on every leg. Write its size as
`d_m`. Its polynomial coefficients `A_k`, for `0<=k<=112m`, sum to
the strict crossing on this type class.

The weighted Smith coefficient estimate gives a constant `C_t>0`
independent of `m,k` such that

```math
\begin{aligned}
\operatorname{rank}_R A_k
 &\le C_t\exp\bigl(tk+m(\log Z_0+tE_p)\bigr),\\
\operatorname{rank}_C A_k
 &\le C_t\exp\bigl(-tk+m(\log Z_\infty+tE_p)\bigr).
\end{aligned}                                         \tag{9}
```

For example, `C_t=(1-exp(-t))^(-2)` follows by bounding the nonnegative
coefficient series `P(s)^m/(1-s)^2` at `s=exp(-t)`. At infinity the
homogeneous profile adds 112 to every channel. This cancels the term
`t(112m-k)` from the complementary coefficient index, leaving exactly
the second line of (9). The raw reciprocal profile already contains
the full degree charge.

Every coefficient, including endpoints and zero coefficients, is
covered. Both ranks can also be capped at `d_m^2`. Their product is
at most `C_t^2 exp(mK)`, independently of `k`.

The type-class size satisfies `log d_m=mH+O(log m)`, effectively.
Integral bounds for `log(n!)` give, for instance,
`log d_m >= mH-24(log(m+1)+1)`. Apply the
[finite-block amplification theorem](finite-block-crossing-amplification.md)
to these at most `112m+1` pieces. Its rank-table bound is at most

```math
\frac{\log(112m+1)+\log C_t+mK/2}{2\log d_m}
  \longrightarrow\frac K{4H}<\frac57.                 \tag{10}
```

Thus one sufficiently large fixed multiple of the type denominator
already gives the strict parameter in (1). All ingredients are
effective finite-field computations. Rational bounds for the logarithms
of the explicitly positive rational margin in (5) determine a sufficient
finite block size. Amplification supplies polynomially many pieces,
explicit factorizations, and every binary label width. No continuous
capacity estimate is needed for this step.

## 5. Reproduce the certificate

```sh
OPENBLAS_NUM_THREADS=1 python tools/certify_stable_hot_code_pressure.py \
  --output evidence/stable-hot-code-rational-pressure.json
```

The [checker](../tools/certify_stable_hot_code_pressure.py) reconstructs
the full tensor from its label recipe, removes exact scalar content,
recomputes all 1,152 local channels, constructs the positive rational
type, and proves (5) by integer arithmetic. Its source and tensor hashes
are retained alongside the comparison. This certificate does not rely
on the numerical search having found a minimum.
