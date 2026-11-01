# Remove the degree restriction from the Frobenius rational-star trace

## Result

The [Frobenius rational-star theorem](frobenius-rational-trace-contraction.md)
extends to arbitrary numerator and denominator formal degrees. An
additional polynomial-size coefficient calculation accounts for the
polynomial part at infinity. The algorithm still does not enumerate
the subfield, expand its high-degree numerator, or divide by a
numerator coefficient.

Keep the original hypotheses: `b=2^w`, an explicit coefficient field
`E=F_(2^d)` with `w|d`, and

```math
R(U)=c\,
\frac{\prod_{i=1}^{r_+}(a_i+b_iU^{2^{e_i}})}
     {\prod_{j=1}^{r_-}(c_j+d_jU^{2^{f_j}})},
\qquad 0\le e_i,f_j<w.                                  \tag{1}
```

Each denominator slope `d_j` is nonzero, and every denominator root
lies outside `F_b`. Numerator coefficients may vanish. There is no
restriction comparing the two formal degrees, and `r_-=0` is allowed.

The complete and nonzero-subfield sums of (1) have deterministic
polynomial bit complexity in `r_++r_-`, `w`, and `d`, once the
explicit field is supplied. The same bound and distinct-root rational
circuit interpretation as in the original theorem continue to hold.

This extends one-variable contraction. It does not resolve the
multi-variable representation growth described in the original note.

## 1. Reverse the rational function at infinity

Define the formal degrees

```math
N=\sum_i2^{e_i},\qquad D=\sum_j2^{f_j},\qquad
\delta=N-D.
```

The series

```math
H_\infty(Z)=Z^\delta R(1/Z)
 =c\,
   \frac{\prod_i(b_i+a_iZ^{2^{e_i}})}
        {\prod_j(d_j+c_jZ^{2^{f_j}})}                    \tag{2}
```

is regular at zero. Every denominator constant in (2) is the
nonzero original slope `d_j`. Its numerator constants may vanish;
the coefficient algorithm already permits this.

The use of formal degrees is deliberate. If a leading numerator
coefficient is zero, (2) merely has additional initial zero
coefficients. No actual-degree detection or division by that leading
coefficient is required.

## 2. The missing trace contribution

Conceptually write `R=P+R_proper`, where `P` is its polynomial part.
This decomposition is used only in the proof. The algorithm does not
perform high-degree polynomial division.

For every positive integer `k`,

```math
\sum_{u\in F_b}u^k=
\begin{cases}1,&b-1\mid k,\\0,&b-1\nmid k,\end{cases}
\qquad \sum_{u\in F_b}1=0.                              \tag{3}
```

The coefficient of `U^k` in `P`, for `k>=0`, equals
`[Z^(delta-k)] H_infinity`. The proper rational part has only
positive powers after substituting `U=1/Z`, so it cannot affect
these coefficients.

Consequently the polynomial part contributes exactly

```math
I_\infty=
\sum_{\ell=1}^{\lfloor\delta/(b-1)\rfloor}
 [Z^{\delta-\ell(b-1)}]H_\infty(Z),                     \tag{4}
```

with an empty sum when `delta<b-1`. Equivalently,

```math
I_\infty=
[Z^{\delta-b+1}]
       \frac{H_\infty(Z)}{1+Z^{b-1}},                   \tag{5}
```

where a negative coefficient index means zero. Formula (5) is the
infinity-residue correction for `R(U)/(U^b+U)` in characteristic two.
Equation (3) gives an elementary proof without needing to evaluate
that residue analytically.

Add (4) to the finite-pole contributions from the original rational-star
theorem. Those contributions depend only on the principal parts at
the denominator roots and remain valid regardless of the polynomial
part. Finally, adding `R(0)` gives the nonzero-subfield sum.

## 3. Cost and uniformity

One has `delta<=r_+ b/2`, so (4) has at most `r_+` coefficient
queries. Every queried index has bit length `O(w+log(r_++1))`.
After dividing each denominator in (2) by `d_j`, each query has
exactly the Frobenius-binomial and inverse-binomial form handled by
the binary-carry algorithm.

Thus the infinity correction is polynomial work and does not change
the original polynomial bound. It also covers the denominator-free
case: the entire trace then comes from (4). For example,
`R(U)=U^(b-1)` has full-field sum one, although it has no finite pole.
Omitting the infinity term would give an incorrect zero.

On the promised domain of pairwise distinct denominator roots, the
finite-pole expression has the rational arithmetic circuit described
in the original theorem. The additional part (4) uses only inverses
of the nonzero `d_j`; its loop bounds depend on the exponent data.
Therefore the same polynomial-size rational circuit statement holds.
Over the explicit finite field, exponentiation to `2^d-2` replaces
those nonzero inversions by division-free circuits.

The coefficient field and its construction remain separately charged.
Zero denominators on the summation field are excluded here; the
[relative-linearized theorem](relative-frobenius-trace-contraction.md)
provides a different, explicitly zero-preserving interpretation when
such values are intended.

## 4. A polynomial-size circuit for zero-safe binomial quotients

There is also a separate, stronger all-coefficient interface for binomial
factors. It changes the value assigned to a vanishing denominator, rather
than treating a genuine pole as an ordinary rational value.

For `x` in the coefficient field `E=F_(2^d)`, define `iota(0)=0` and
`iota(x)=x^-1` otherwise. The positive exponent identity

```math
\iota(x)=x^{2^{d+1}-3}
        =\prod_{h\in\{0\}\cup\{2,\ldots,d\}}x^{2^h}       \tag{6}
```

holds also at `d=1`. It uses exactly `d` Frobenius factors. Apply (6)
to each denominator binomial and reduce the exponent of the scalar
summation variable modulo `w`, using `u^(2^w)=u` on `F_b`.
The resulting product has

```math
s=r_++d r_-
```

ordinary Frobenius binomial factors. Its coefficient powers are
constructed by repeated squaring. This is the binomial specialization
of the [zero-preserving inverse expansion](relative-frobenius-trace-contraction.md).

The one-dimensional binary-carry coefficient algorithm now applies
directly. The polynomial degree is at most `s2^(w-1)`, so only at most
`s` positive multiples of `b-1` contribute to its full-field sum.
The result is a fixed division-free circuit of size

```math
O\bigl((s+1)^3(w+\log(s+2))\bigr)+O(dr_-),                 \tag{7}
```

with topology determined only by the exponent data and the two field
degrees. Adding the value at zero gives the nonzero-field sum.
No root grouping, pole promise, formal degree inequality, or test of
an input coefficient is needed.

For field-valued coefficients this circuit computes exactly

```math
\sum_{u\in F_b}c\prod_i N_i(u)\prod_j\iota(D_j(u)).        \tag{8}
```

It also defines a polynomial extension to every commutative `E`-algebra.
There it computes the same sum with each `iota(x)` replaced by the
explicit positive power in (6). That power is not asserted to invert
zero divisors. The distinction matters for algebras with nilpotents.

The infinity method from Sections 1--3 remains useful: it gives the
ordinary rational trace without introducing `d` separate factors per
inverse. The zero-safe version instead gives the uniform all-coefficient
circuit in (7). The two agree whenever their field-valued denominators
are nonzero on the summation field.

## 5. Implementation and exact controls

The [module](../tools/frobenius_trace_infinity_correction.py) combines
the original finite-pole calculation with (4). The
[checker](../tools/check_frobenius_trace_infinity.py) compares the result
with literal subfield sums, including positive-degree polynomial
parts, repeated roots, canceled poles, zero leading numerator
coefficients, and denominator-free monomials.

It also compares the original and extended routines throughout their
common proper-degree domain, checks constant functions, and rejects
a genuine pole on the subfield. The
[record](../evidence/frobenius-trace-infinity-correction.json) contains
the exact source hashes and query counts.

The separate zero-safe circuit is compared with literal field-valued
evaluation, including identically zero denominator factors and `d=1`.
Over `F4[epsilon]/(epsilon^2)`, the checker evaluates its stated
polynomial extension directly. A gate-count control verifies that
coefficient values do not change the circuit schedule.

```sh
python3 tools/check_frobenius_trace_infinity.py \
  --output evidence/frobenius-trace-infinity-correction.json
```
