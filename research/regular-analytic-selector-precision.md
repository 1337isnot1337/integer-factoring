# The sharp precision cost of regular residue selectors

## Result and scope

The [regular analytic compiler](pole-cauchy-analytic-weights.md) makes
large encoded powers and unit-denominator circuits cheap to contract
against an implicit Cauchy inverse. Can such a circuit encode a fine
residue test and recover the singular factor-count information?

For a direct pointwise parity selector, the required external dyadic
denominator has an exact lower bound. If `F` is a regular odd-input
circuit and `2^-r F(x)` is integral with parity equal to the indicator
of one odd residue class modulo `2^(h+1)`, then

\[
 \boxed{r\ge 2^{h+1}-h-2.}                         \tag{1}
\]

This bound is attained. Thus the denominator budget in the compiler
cannot generally be reduced to a polynomial in the number `h` of
selected residue bits. The result concerns this representation and
this pointwise selector contract. It is not a lower bound for
factoring, arbitrary arithmetic circuits with nonunit inverses, or
other ways to aggregate matrix information.

## 1. Finite differences of regular odd-input functions

Let `F` use integer constants, the odd variable, addition, subtraction,
multiplication, encoded integer powers, and inverses odd on the odd
input class. The compiler proof gives

\[
 H(X)=F(1+2X)=\sum_{i\ge0}c_iX^i,
 \qquad c_i\in2^i\mathbb Z_2.                     \tag{2}
\]

The series converges on every 2-adic integer. The same coefficient
property holds after translation by an integer `X_0`: the coefficient
of `X^i` in `H(X_0+X)` is
`sum_(j>=i) c_j binom(j,i) X_0^(j-i)`, a convergent sum of multiples
of `2^i`.

For `d>=1`, write the forward difference as

\[
 \Delta^dH(X_0)=
 \sum_{j=0}^d(-1)^{d-j}\binom dj H(X_0+j).         \tag{3}
\]

The elementary monomial identity
`Delta^d X^i at 0 = d! S(i,d)` uses the integer Stirling number of the
second kind, and vanishes for `i<d`. Applying it to the translated
series gives

\[
 \boxed{\displaystyle
 v_2(\Delta^dH(X_0))\ge d+v_2(d!).}               \tag{4}
\]

Indeed, each surviving term is a multiple of `2^i d!`, with `i>=d`.
The passage through the convergent series is valid because (3) is a
fixed finite linear combination of evaluations. The valuation of zero
is interpreted as infinity.

## 2. A full residue period forces the lower bound

Suppose `G(X)=2^-r H(X)` is integral and its parity is one exactly
when `X` lies in a prescribed class modulo `2^h`. Choose an integer
`X_0` in that class, and put `d=2^h-1`. Among
`X_0,X_0+1,...,X_0+d`, only the first argument has odd `G` value.
Its coefficient in (3) is `(-1)^d`, also odd. Therefore

\[
 \Delta^dG(X_0)\equiv1\pmod2,
 \qquad v_2(\Delta^dH(X_0))=r.                    \tag{5}
\]

Equations (4) and (5) imply `r>=d+v_2(d!)`. Counting the multiples
of each power of two in `d!` gives

\[
 v_2((2^h-1)!)
 =\sum_{j=1}^{h-1}(2^{h-j}-1)
 =2^h-h-1.                                      \tag{6}
\]

Combining these expressions proves (1). The argument needs parity
agreement only on one complete consecutive residue period; exact
zero-one values or a formula outside that period are not required.

## 3. An equality construction

Let `d=2^h-1` and set

\[
 G_h(X)=\binom{X+d}{d}
       =\frac{\prod_{j=1}^d(X+j)}{d!}.            \tag{7}
\]

This is integer-valued on all integers. Its parity is one precisely
for `X=0 mod 2^h`. One way to see this without a residue enumeration
is to work in `F_2[T]/(T^(2^h))`. The identity
`(1+T)^(2^h)=1` in that ring shows that the coefficient of `T^d` in
`(1+T)^n` depends only on `n mod 2^h`, including negative integers.
Writing that residue in binary factors the power into
`prod_(j with bit j set)(1+T^(2^j))`. Its coefficient of `T^(2^h-1)`
is one exactly when every one of the `h` bits is set. With `n=X+d`,
this is exactly `X=0 mod 2^h`.

Now put `X=(x-1)/2`, let `u` be the odd part of `d!`, and define

\[
 F_h(x)=u^{-1}\prod_{j=1}^d(x+2j-1),
 \qquad r=d+v_2(d!).                             \tag{8}
\]

The inverse of the odd integer `u` is a permitted regular gate, and

\[
 2^{-r}F_h(1+2X)=G_h(X).                         \tag{9}
\]

Thus equality holds in (1). Translating `x` by an even integer gives
any prescribed odd residue class. This construction establishes
sharpness; its product has `d` factors and is not claimed to have a
short circuit independent of `d`.

For example, selecting one of 256 odd residue classes modulo 512
requires at least 502 denominator bits. A parity computation through
the regular compiler consequently requires at least 503 internal
bits even before the matrix inverse normalization is charged.

## 4. Consequence for joint product selectors

Suppose a regular bivariate kernel `F(a,b)` divided by `2^r` is
integral and has parity equal to the indicator of

\[
 ab\equiv N\pmod{2^{h+1}}                        \tag{10}
\]

for every odd pair. Fix any odd `a`. It is invertible modulo
`2^(h+1)`, so (10) selects exactly one odd class of `b`. Substituting
this fixed odd integer for `a` preserves regularity of the circuit
in `b`; (1) applies unchanged. The
[joint inverse-contraction algorithm](pole-cauchy-analytic-kernels.md)
charges this external denominator and therefore does not make a
fine pointwise product selector inexpensive merely by allowing a
joint circuit or binary-encoded powers.

This does not rule out a contraction identity that recovers the
factor count without representing a pointwise parity selector.
It also does not apply to a residue test obtained by branching on
digits or by supplying a nonregular singular kernel.

## 5. Reproducible exact controls

The [checker](../tools/check_regular_analytic_selector_precision.py)
tests (4) on 40 regular circuits, twelve difference orders, and four
integer translations, including negative inputs. It then constructs
the equality examples for `h=1,...,8`, compiles them at precision
`r+1`, and checks every argument in three complete periods. It verifies
that the `d`th finite difference has valuation exactly `r` and that
the compiled polynomial agrees pointwise with the original circuit.

The [evidence](../evidence/regular-analytic-selector-precision.json)
records all eight sharp budgets, circuit sizes, finite control counts,
and source hashes. These controls supplement the general proof;
they do not substitute finite sampling for (1).

```sh
python3 tools/check_regular_analytic_selector_precision.py \
  --output evidence/regular-analytic-selector-precision.json
```
