# Integer Smith profiles strengthen the shared-minor pressure floor

## Result and scope

Keeping the two integer Smith profiles separate gives a stronger
degree-sensitive version of the
[shared-minor bound](shared-minor-smith-floor.md).

Two consequences are:

- Every bare degree-one crossing has `rho>=3/4`, for every physical
  alphabet size `q>=2`.
- Exact certificates exclude a half-rate improvement for seven-letter
  crossings of every degree at most eleven. They also exclude several
  low-degree cases on larger alphabets, listed in Section 4.

The results hold over any coefficient field and allow arbitrary fixed
auxiliary dimension in the same rectangular-cut interface. They concern
the **bare, ungauged local Smith-pressure certificate**. Independent
endpoint gauges, adaptive cancellation or channel constructions, and
actual coefficient-rank methods remain outside the theorem. This is not
a lower bound on integer factoring.

## 1. Integer constraints on the two profiles

Remove a common scalar monomial from the polynomial crossing. Let its
actual remaining degree be `M>=1`, and put `n=q^2`. Use the full-rank
cut matrices with `n` columns; both auxiliary indices stay on the other
side of each cut. A fixed constant auxiliary contraction must give the
swap at parameter one.

Let `mu_1<=...<=mu_n` be the Smith orders of the infinity reversal,
and let `nu_1<=...<=nu_n` be those of the other cut at zero. These
are nonnegative integers. The common degree and constant entry give

```math
\mu_1=\nu_1=0,\qquad
\sum_i\mu_i\le nM,\qquad \sum_i\nu_i\le nM.       \tag{1}
```

The shared `q`-by-`q` minor gives

```math
\sum_{i=1}^q(\mu_i+\nu_i)\le qM.                 \tag{2}
```

These facts are recalled here so the exact relaxation is explicit.
The shared minor is obtained by fixing `y=u` and varying `x,v`
after the fixed auxiliary contraction. Its determinant is nonzero
at one. Its zero order is at most its degree. Minimum determinantal
valuations in the two cuts therefore imply (2), including rectangular
auxiliary matrices by Cauchy--Binet.

Put `r=s^-1` for `s>1`. The restored pressure is

```math
B(s)=r^{-M}\left(\sum_i r^{\mu_i}\right)
                 \left(\sum_i r^{\nu_i}\right).  \tag{3}
```

The degree-zero case is constant full rank and has `B=q^4`.

## 2. Balanced integer sums give a finite exact relaxation

For integers `m>=1`, `B>=0`, write `B=ma+b`, `0<=b<m`, and define

```math
J_m(B;r)=(m-b)r^a+b r^{a+1}.                     \tag{4}
```

Among `m` nonnegative integer orders with total at most `B`, this
is the smallest possible sum of their `r`-powers. One first uses the
whole budget, since `r<1`. Whenever two orders differ by at least
two, moving one unit from the larger to the smaller lowers the sum.
Repeating leaves exactly the two adjacent orders in (4).

For `0<=b<=qM`, put

```math
F_b(r)=1+J_{q-1}(b;r)+J_{n-q}(nM-b;r).            \tag{5}
```

If one profile has first-`q` sum `b`, equations (1) and (4) give
`sum_i r^(mu_i)>=F_b(r)`. This relaxes the ordering constraint
between the two groups and is therefore a valid lower bound even
when the two separately balanced groups would overlap.

Let the two actual prefix sums be `b,c`. Their sum is at most `qM`.
We may increase their budgets until the sum is exactly `qM` without
increasing the relaxed product. Indeed, while the sum is smaller,
one of them is below `qM/2`. For such an integer `b`,

```math
F_{b+1}(r)-F_b(r)
=(1-r)\left[
 r^{\lfloor(nM-b-1)/(n-q)\rfloor}
 -r^{\lfloor b/(q-1)\rfloor}\right]\le0.         \tag{6}
```

The exponent in the first term is at least the second one. This
follows already before taking floors from
`(q-1)(nM-1)>=b(n-1)` in this range. Thus the budget can be filled
one unit at a time.

We obtain the finite, degree-sensitive bound

```math
\boxed{\quad
B(s)\ge\min_{0\le b\le qM}
 r^{-M}F_b(r)F_{qM-b}(r),\qquad r=s^{-1}.
\quad}                                          \tag{7}
```

Only `b<=floor(qM/2)` needs to be checked by symmetry. This keeps
the integer allocation of each profile. Averaging their orders
first and applying real Jensen would lose precisely the rounding
information used here.

## 3. Every affine crossing has a three-quarter floor

Take `M=1`. For `1<=b<=q-1`, formula (5) simplifies to

```math
F_b(r)=(q-b)(1-r)^2+q^2r.                        \tag{8}
```

In the product `F_b F_(q-b)`, the only dependence on `b` is the
factor `b(q-b)` multiplying `(1-r)^4`. It is smallest at `b=1`
or `q-1`. The endpoint choices `b=0,q` give a larger product:
their difference from the `b=1` product is
`q^2 r(1-r)^2+(1-r)^4`, which is nonnegative.

Consequently every normalized degree-one crossing satisfies

```math
\begin{aligned}
B(s)
&\ge q^4r+q^3(1-r)^2+(q-1)\frac{(1-r)^4}{r}\\
&\ge q^3,
\end{aligned}                                    \tag{9}
```

because `q r+(1-r)^2=1+(q-2)r+r^2>=1`.
This proves

```math
\boxed{\rho_{\rm bare}\ge\frac34}
```

for every alphabet size, field, and permitted fixed auxiliary
dimension. At `q=2`, the first line of (9) is `(1+r)^4/r`, whose
minimum is the familiar `256/27`.

## 4. Exact low-degree exclusions

For the pairs `(q,M)` recorded by the checker, all real scales
obey

```math
B(s)>q^{577/204},\qquad
\rho_{\rm bare}>\frac{577}{816}>\frac1{\sqrt2}.   \tag{10}
```

The last comparison is the exact integer inequality
`2*577^2>816^2`.

The certified ranges are:

| Physical alphabet | Certified positive degrees |
| --- | --- |
| `q=7` | `1` through `11` |
| `q=8` | `1` through `5` |
| `q=9` | `1` through `3` |
| `q=10` | `1` through `3` |
| `q=11` | `1` and `3` |

Degree one for every larger alphabet is already covered by (9).
These ranges are finite exclusions, not assertions that another
degree admits a crossing. The isolated `q=11,M=3` certificate refers
to the actual normalized degree; it does not exclude degree two.

For each `q`, the checker chooses a rational `T=p/r0` satisfying
`T>q^(577/204)`. It verifies this by the integer comparison

```math
p^{204}>q^{577}r_0^{204}.                         \tag{11}
```

For every budget `b` in (7), it then proves strict positivity on
the whole interval `0<=r<=1` of the integer polynomial

```math
r_0 F_b(r)F_{qM-b}(r)-p r^M.                     \tag{12}
```

The proof uses exact rational Bernstein coefficients with dyadic
subdivision. Every leaf has strictly positive coefficients, so
there is no uncovered evaluation interval. No floating optimization
enters the replay.

For orientation, the numerical relaxation at `q=7,M=11` is about
`rho=0.70746247`. At `M=12` it can attain approximately
`0.70623045`, so integer allocation alone no longer gives the
exclusion. The formal minimizing orders there are
`(0,7 repeated six times,13 repeated forty-two times)` in both
profiles. Algebraic compatibility remains a separate requirement;
this relaxation does not construct a crossing with those orders.

## 5. Replay and search implication

Run:

```sh
python3 tools/check_integer_shared_minor_floor.py
```

The [checker](../tools/check_integer_shared_minor_floor.py) proves
the thresholds and every finite interval certificate using integers
and exact fractions. It also checks the integer balancing formula,
the budget monotonicity, and the affine simplification on finite
independent controls. The full certificate summary is saved in
[integer-shared-minor-floor.json](../evidence/integer-shared-minor-floor.json).

For example, any constant linear restriction of three elementary
binary sites has degree at most three. Its seven-letter bare
pressure is therefore excluded before enumerating physical subspaces.
Likewise an eight-letter bare restriction of four sites lies in the
certified range. Gauged and adaptive interfaces require their own
analysis and are not excluded by this argument.
