# Prime-power vacuum factors do not have finite polynomial locality

## Result and scope

A prime-power replacement of the bosonic mode set leads to the
Artin--Hasse scalar factor. Its reduction in characteristic `p` cannot
be made rational by multiplying by a nonzero Laurent polynomial. In
particular, its expansions at zero and infinity cannot agree outside
a finite range after one such multiplication.

This rules out the specific attempt to keep a fixed finite polynomial
locality relation while replacing the ordinary rational vacuum factor
by the Artin--Hasse factor. It does not exclude every Witt-vector,
prime-power, finite-cutoff, or matrix-valued construction. It does not
prove a lower bound on factoring. No crossing family or improved
factoring exponent is claimed here.

The independent all-size input is Joe Kramer-Miller's
[Theorem 1.1 in arXiv:2404.06968v1](https://arxiv.org/html/2404.06968v1):
for every prime `p`, the reduction of the Artin--Hasse exponential is
transcendental over `F_p(z)`. The paper's proof covers `p=2,3` as well
as the larger primes. The coefficient calculations below are finite
controls; they do not prove that external theorem.

## 1. Why the prime-power scalar appears

The Artin--Hasse series is

```math
E_p(z)=\exp\left(\sum_{j\ge0}\frac{z^{p^j}}{p^j}\right)
      =\sum_{n\ge0} a_n z^n.
\tag{1}
```

For comparison, including every positive integer as a mode gives

```math
\exp\left(\sum_{m\ge1}\frac{z^m}{m}\right)=\frac1{1-z}.
\tag{2}
```

The scalar relation can be seen directly with independent variables
`x_j` over the rationals. Let multiplication by `x_j` and
differentiation with respect to `x_j` be the creation and annihilation
operators. Normal ordering the two degreewise formal exponentials

```math
\exp\left(\sum_{j\ge0}\frac{w^{p^j}}{p^j}
                         \frac{\partial}{\partial x_j}\right)
\quad\hbox{and}\quad
\exp\left(\sum_{j\ge0}x_jz^{p^j}\right)
```

produces the scalar factor `E_p(wz)`. Indeed, the commutator of the
two exponents is the scalar `sum_j (wz)^(p^j)/p^j`; all further
commutators vanish. This calculation is valid coefficient by
coefficient, with only finitely many contributing modes at each
degree.

Equation (1) has coefficients in the localization `Z_(p)`, so its
reduction `e_p(z)` in `F_p[[z]]` is defined. This integrality assertion
concerns the scalar series. The displayed individual operator
exponentials have denominators and do not, by themselves, define
operators over `F_p`. Any integral operator model reproducing their
scalar factor still has to satisfy the obstruction below if it asks
for finite polynomial locality.

## 2. The elementary two-expansion obstruction

Let `K` be any field and let `f(z)=sum_(n>=0) b_n z^n` be a formal
power series. Regard

```math
f_+(z)=\sum_{n\ge0}b_nz^n,\qquad
f_-(z)=\sum_{n\ge0}b_nz^{-n}
```

as coefficientwise formal distributions on the integer powers of `z`.
Multiplication of either distribution by a Laurent polynomial is
well defined. No product of two arbitrary bilateral distributions
is used.

**Lemma.** If there are a nonzero Laurent polynomial `Q(z)` and a
Laurent polynomial `L(z)` for which

```math
Q(z)\bigl(f_+(z)-f_-(z)\bigr)=L(z),
\tag{3}
```

then `f(z)` is rational over `K(z)`.

Write `v` for the largest exponent in `Q`. The distribution `Q f_-`
has no terms of degree greater than `v`. Consequently, (3) implies
that all sufficiently high coefficients of `Q f_+` vanish. That
product is already bounded below, since `f_+` has no negative powers
and `Q` has finite support. Thus `Q f_+` is a Laurent polynomial,
and dividing by the nonzero `Q` proves rationality.

The same argument permits different nonzero multipliers at the two
ends. If

```math
Q_+(z)f_+(z)-Q_-(z)f_-(z)=L(z),\qquad Q_+\ne0,
\tag{4}
```

then the positive tail again forces `Q_+ f_+` to be a Laurent
polynomial. Finite Laurent corrections therefore cannot repair the
obstruction. Replacing `f_-` by `z^a f_-(c z)` for fixed nonzero `c`
and integer `a` also leaves the proof unchanged.

Apply the lemma with `K=F_p` and `f=e_p`. Kramer-Miller's theorem
shows that neither (3) nor (4) can hold. In particular, no finite
power of `1-z`, and no other nonzero polynomial of finite degree,
removes all tails of this scalar commutator.

This conclusion is about one stable infinite series. It does not
give a lower bound on the degree of a polynomial that cancels only
a specified finite prefix. Such cutoff-dependent constructions
need their own endpoint and rank analysis. A finite numerical
prefix cannot distinguish rationality from transcendence.

## 3. Exact finite diagnostics

Differentiating (1) gives the recurrence

```math
a_0=1,\qquad
n a_n=\sum_{p^j\le n}a_{n-p^j}.
\tag{5}
```

This computes exact rational coefficients. Their reduced denominators
are checked to be prime to `p` before modular reduction. A separate
calculation uses the Euler product

```math
E_p(z)=\prod_{\substack{m\ge1\\p\nmid m}}
                    (1-z^m)^{-\mu(m)/m}.
\tag{6}
```

For completeness, the coefficient of `z^n` in the logarithm of the
right side is

```math
\frac1n\sum_{\substack{m\mid n\\p\nmid m}}\mu(m).
```

The divisor sum is one when `n` is a power of `p` and zero otherwise.
Thus (6) and (1) agree as rational power series, without importing
the coefficients generated by (5).

For a finite coefficient list `b_0,...,b_(N-1)`, define its recurrence
length to be the smallest `d` for which some `q_0=1,q_1,...,q_d`
satisfies

```math
\sum_{j=0}^d q_j b_{n-j}=0\qquad(d\le n<N).
\tag{7}
```

Berlekamp--Massey supplies one candidate. Independent modular Gaussian
elimination verifies that (7) is feasible at that length and infeasible
at the preceding length. Feasibility is monotone in `d`, by appending
zero coefficients and dropping the earliest equation, so the two
checks certify the stated finite minimum.

The recorded run uses `p=2,3,5,7`, computes 256 coefficients in each
field, and compares the first 65 rational coefficients with (6).
It independently verifies recurrence lengths at prefix lengths
8, 16, 32, 64, and 256. A geometric series is the rational negative
control and retains recurrence length one throughout.

```sh
python3 tools/check_prime_power_vacuum.py \
  --output evidence/prime-power-vacuum-locality.json
```

The [checker](../tools/check_prime_power_vacuum.py) and
[recorded evidence](../evidence/prime-power-vacuum-locality.json)
contain the exact coefficients and recurrence denominators. No trend
in those finite lengths is used as an asymptotic rank or complexity
claim.

The argument supplements the existing
[graded-space restriction](graded-space-barrier.md): changing the
mode set must preserve the operator and endpoint identities, as well
as supply useful dimension counts. The separate
[arbitrary finite-mode crossings](fixed-mode-crossings.md) already
have valid endpoints and remain outside this vacuum-locality claim.
