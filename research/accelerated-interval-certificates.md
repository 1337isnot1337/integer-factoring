# Exact interval pressure without enumerating the modes

## Result and scope

The interval pressure products can be enclosed using a number of exact
integer operations polynomial in the logarithm of the mode count. This
removes the linear scan through the modes from the earlier certificate
checker. It also provides a reusable way to certify long Gaussian
spectator products surrounding a fixed defect.

Eleven recorded witnesses extend the explicit interval certificates to
\(2^{1636}-1\) modes. The last rational upper bound is less than
\(2\cdot10^{-980}\) above \(1/\sqrt2\). Every bound remains above
that threshold. The factoring score is still \(1/2\); this is an
improvement in finite certificates and their computation, not a new
factoring exponent.

## 1. A finite identity with a positive remainder

For \(0\le u\le1\), the finite geometric identity gives

```math
\frac1{1+t}
=\sum_{m=0}^{N-1}\frac{(1-t)^m}{2^{m+1}}
 +\frac{(1-t)^N}{2^N(1+t)}.
```

Integrating from zero to \(u\) yields

```math
\log(1+u)
=\sum_{m=0}^{N-1}\frac1{2^{m+1}}
       \int_0^u(1-t)^m\,dt+E_N(u),
\qquad 0\le E_N(u)\le\frac{u}{2^N}.
\tag{1}
```

Now let \(u_j=xq^j\), where \(0\le x\le1\), \(0\le q<1\),
and \(0\le j<n\). Define

```math
A_k=\frac1k\sum_{j=0}^{n-1}u_j^k
 =\frac{x^k(1-q^{nk})}{k(1-q^k)},\qquad k\ge1.
\tag{2}
```

With the difference convention \(\Delta A_k=A_k-A_{k+1}\), put
\(D_m=\Delta^m A_1\). Expanding the integrand gives the exact identity

```math
D_m=\sum_{j=0}^{n-1}\int_0^{u_j}(1-t)^m\,dt,
\qquad 0\le D_m\le A_1.
\tag{3}
```

Summing (1) proves the enclosure

```math
\sum_{m=0}^{N-1}\frac{D_m}{2^{m+1}}
\le \sum_{j=0}^{n-1}\log(1+xq^j)
\le \sum_{m=0}^{N-1}\frac{D_m}{2^{m+1}}
        +\frac{A_1}{2^N}.
\tag{4}
```

The bound is uniform even when \(q\) is extremely close to one.
Computing \(q^n\) uses binary exponentiation. The moments in (2)
then require only the successive powers of \(x,q,q^n\), followed
by a triangular difference table. There is no loop of length \(n\).

## 2. Outward integer intervals

The implementation represents an interval by two signed integers
\([L,U]/2^P\). Every multiplication and division rounds the lower
endpoint downward and the upper endpoint upward. Addition, subtraction,
and multiplication by an integer are exact on this grid. All divisors
used for (2) have a certified positive lower endpoint.

Subtraction in the difference table can amplify rounding widths. The
checker retains both signed endpoints throughout, uses a precision
larger than the number of differences, and checks the final inequality
directly. It also intersects each computed \(D_m\) interval with the
proved range \([0,A_1]\). This intersection cannot discard the true
value by (3).

The last remainder in (4) is rounded upward. A negative upper endpoint
for the resulting log pressure ratio is therefore an exact proof that
the original pressure ratio is less than one. Floating logarithms,
square roots, and optimizer output are not part of this decision.

For the witnesses below, \(N=O(\log n)\) terms and
\(P=O(\log n)\) fractional bits suffice. Formula (2) and the difference
table use \(O(N^2+\log n+\log a+\log b)\) interval operations on
integers of \(O(P+\log n+\log a+\log b)\) bits. This is a
certificate-computation statement. It does not assert that the enormous
physical alphabet or a full tensor table has a compact explicit listing.

## 3. Apply the identity to the exact pressure ratio

The [original interval note](quantitative-interval-frontier.md) supplies
the polynomial crossing and its complete profiles. Write \(s=r^a\),
\(R=a/b\), with \(r=p/q>1\). For odd \(n\), its exact pressure
ratio has logarithm

```math
ab\log r
+4b\sum_{j=1}^{n}\log(1+r^{-a(2j-1)})
-a\log2
-2a\sum_{j=1}^{(n-1)/2}\log(1+r^{-8bj}).
\tag{5}
```

For even \(n\), the corresponding expression is

```math
4b\sum_{j=1}^{n}\log(1+r^{-a(2j-1)})
-2a\sum_{j=1}^{n/2}\log(1+r^{-4b(2j-1)}).
\tag{6}
```

Each sum is an instance of (4). For \(1<r\le2\), the same identity
with one term gives \(\log r=\log(1+(r-1))\), and \(x=1\)
gives \(\log2\). The checker explicitly enforces this root range
when (5) is used. Its strict negative upper bound certifies
\(\rho<a/(4b)\).

## 4. Reproducible rational witnesses

Start with \((a_1,b_1)=(3,1)\), and iterate

```math
a_{j+1}=3a_j+8b_j,\qquad b_{j+1}=a_j+3b_j.
\tag{7}
```

Thus \(a_j^2-8b_j^2=1\). At each displayed level, choose the
smallest integer \(k\ge\lceil\log_2 a_j\rceil\) satisfying
\(3(2^k-1)\ge a_j k\), and set

```math
n=2^k-1,\qquad q=\left\lfloor\frac{a_j n}{k}\right\rfloor,
\qquad r=\frac{q+1}{q}.
\tag{8}
```

These are explicit proposals; each row below is separately certified
by integer arithmetic. The last column is explanatory decimal data,
not the test used to accept a certificate.

| Pell level | Mode count | Certified upper bound | Approximate gap above \(1/\sqrt2\) |
|---:|---:|---|---:|
| 5 | \(2^{14}-1\) | \(3363/4756\) | \(3.1261\cdot10^{-8}\) |
| 6 | \(2^{17}-1\) | \(19601/27720\) | \(9.2023\cdot10^{-10}\) |
| 7 | \(2^{20}-1\) | \(114243/161564\) | \(2.7089\cdot10^{-11}\) |
| 8 | \(2^{23}-1\) | \(665857/941664\) | \(7.9743\cdot10^{-13}\) |
| 12 | \(2^{33}-1\) | \(768398401/1086679440\) | \(5.9880\cdot10^{-19}\) |
| 20 | \(2^{55}-1\) | \(1023286908188737/1447146223759344\) | \(3.3764\cdot10^{-31}\) |
| 40 | \(2^{106}-1\) | \(a_{40}/(4b_{40})\) from (7) | \(8.0613\cdot10^{-62}\) |
| 80 | \(2^{209}-1\) | \(a_{80}/(4b_{80})\) from (7) | \(4.5951\cdot10^{-123}\) |
| 160 | \(2^{414}-1\) | \(a_{160}/(4b_{160})\) from (7) | \(1.49306\cdot10^{-245}\) |
| 320 | \(2^{821}-1\) | \(a_{320}/(4b_{320})\) from (7) | \(1.57630\cdot10^{-490}\) |
| 640 | \(2^{1636}-1\) | \(a_{640}/(4b_{640})\) from (7) | \(1.75695\cdot10^{-980}\) |

Writing \(\rho_{80}=a_{80}/(4b_{80})\) for the level-80 rational upper
bound, the certificate is explicitly

```math
\rho<
\frac{8771616955900347725453582478444206703138404016040525014107137}
{12404939662976074530083669466978070217253381205037089626615744}.
```

It is above the target because its square exceeds \(1/2\) by exactly
\(1/(16b_{80}^2)\). Moreover \(\rho_{80}+1/\sqrt2>7/5\), so

```math
0<\rho_{80}-\frac1{\sqrt2}
<\frac5{112b_{80}^2}<5\cdot10^{-123}.
```

Both strict rational comparisons are reproducible without decimal
approximations. The interval floor still prevents this displayed family
from crossing the target.

The same Pell identity certifies the later gap bounds without
subtracting nearly equal numerical square roots. For level 640, exact
integer comparison gives

```math
0<\frac{a_{640}}{4b_{640}}-\frac1{\sqrt2}
<\frac5{112b_{640}^2}<2\cdot10^{-980},
```

where the final inequality is checked as
\(5\cdot10^{980}<224b_{640}^2\). The corresponding upper bounds
at levels 160 and 320 are \(2\cdot10^{-245}\) and
\(2\cdot10^{-490}\). These compare the rational certificate values
with the limiting threshold; they are not bounds below that threshold.

## 5. Evidence and reuse

Run the [checker](../tools/check_accelerated_interval_pressure.py):

```sh
python3 tools/check_accelerated_interval_pressure.py \
  --output evidence/accelerated-interval-pressure.json
```

The [evidence](../evidence/accelerated-interval-pressure.json) includes
all rational roots, interval endpoints, precision choices, and implicit
rational type probabilities. Independent `Fraction` controls compare
256 moment and difference values with direct finite sums. Four pressure
signs are also checked against the earlier direct-product interval
method, covering both parities and both signs.

The [extended record](../evidence/accelerated-interval-pressure-extended.json)
contains levels 160, 320, and 640, including strict outward interval
endpoints and exact integer checks of the displayed gap bounds. Reproduce
their pressure certificates with:

```sh
python3 -X int_max_str_digits=100000 tools/check_accelerated_interval_pressure.py \
  --levels 160 320 640 --output /tmp/accelerated-interval-pressure-extended.json
```

The Python integer-to-decimal allowance is increased only to print the
large, explicitly constructed certificate endpoints. The largest
recorded case uses 3,352 Euler terms and 20,040 fractional bits. The
argument and implementation apply to arbitrary finite
geometric products in the stated domain, so the same routine can
certify the spectator contributions in a proposed bounded-defect
pressure witness. The finite defect orders and admissible endpoint
must still be established separately.
