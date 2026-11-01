# A wider exact pressure-gain interval for the eight-letter retraction

The [explicit eight-letter nonorthogonal crossing](eight-letter-reflection-retraction.md)
improves its coordinate control for every finite scale

\[
s\ge\frac{10}{9}.
\tag{1}
\]

This strengthens the earlier interval `s>=3/2`. It uses the same crossing,
profiles, physical heights, and rational-type certificate `rho<18/25`.
It does not change the factoring exponent.

Let `c(r),z(r),i(r)` be the three exact counting polynomials in the
original note, with `r=1/s`. Define

\[
\begin{aligned}
K(r)&=r^{16}(1+r)(1-r^4),\\
A(r)&=(1+r^3)(1+r^4),\\
B(r)&=1+r+r^3+r^4+r^7.
\end{aligned}
\tag{2}
\]

The profile differences factor as

\[
z=c-KA,\qquad i=c+Kr^4B.
\tag{3}
\]

Consequently the complete pressure difference is

\[
c^2-zi
=Kc\bigl(1+r^3-r^5-r^8-r^{11}\bigr)+K^2Ar^4B.
\tag{4}
\]

For `0<r<=9/10`, the factors `K,c,A,B` are positive, and

\[
1+r^3-r^5-r^8-r^{11}
\ge1-r^8-r^{11}
\ge1-(9/10)^8-(9/10)^{11}>0.
\tag{5}
\]

The first inequality uses `r^3>=r^5`; the final comparison is the
integer inequality

\[
9^8\,10^3+9^{11}<10^{11}.
\tag{6}
\]

Equations (4)--(6) prove strict positivity and hence (1). No numerical
root estimate or pressure optimization is used. The same type-energy
polynomial occurs on both sides, so the comparison preserves the
common-type interpretation of the original result.

The short [checker](../tools/check_eight_letter_pressure_factorization.py)
reconstructs the two crossings through the existing constructor,
recomputes the four profiles, and checks every coefficient in (3)--(4)
by integer arithmetic.

```sh
python3 tools/check_eight_letter_pressure_factorization.py \
  --output evidence/eight-letter-pressure-gain-factorization.json
```

The [evidence](../evidence/eight-letter-pressure-gain-factorization.json)
records the exact positive lower bound in (5), the coefficient
identities, and source hashes. As in the original packet, the result
extends unchanged to every coefficient field of characteristic two.
