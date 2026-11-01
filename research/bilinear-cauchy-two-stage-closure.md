# Two eliminations for a bounded-degree Cauchy interaction

## Result and scope

A restricted two-variable interaction admits a compact representation
after the first elimination and a polynomial-size exact circuit after
the second. The interaction is the zero-safe Cauchy factor
$\iota(1+\gamma UV)$, with one variable in a binary subfield and the other
on its complementary norm-one torus. The first variable may also carry
an arbitrary product of Frobenius binomials whose formal degree is
strictly below the subfield size.

The first elimination produces one factorized torus star, at most one
point correction, and a constant term for the degenerate tag. Thus the
second elimination does not require a table of first-stage values.

This theorem concerns field multiplication inside that Cauchy factor.
It does not represent the integer multiplication phase
$\zeta^{sxy}$. Section 6 identifies the additional restrictions that
prevent a whole-network conclusion.

## 1. Setup and charged circuit bound

Let $b=2^w$, $w\ge1$, and let the explicit field $E=\mathbb F_{2^d}$
contain $\mathbb F_{b^2}$, so $2w$ divides $d$. Write

$$
T=\{v\in E:v^{b+1}=1\},\qquad
\iota(x)=
\begin{cases}
x^{-1},&x\ne0,\\
0,&x=0.
\end{cases}
$$

The torus has $b+1$ elements. Fix a Frobenius-binomial product

$$
P(U)=\prod_{i=1}^{r_P}(a_i+b_iU^{2^{e_i}}),\qquad
\Delta_P=\sum_i2^{e_i}<b.                              \tag{1}
$$

Every $e_i$ is a nonnegative integer. The degree condition is formal:
it remains the same when coefficients
vanish. In particular, every occurring $e_i$ is less than $w$.
The empty product is allowed.

Let the outer weight be any zero-safe Frobenius rational star

$$
W(V)=cV^h
\prod_{i=1}^{r_+}(A_i+B_iV^{2^{f_i}})
\prod_{j=1}^{r_-}\iota(C_j+D_jV^{2^{g_j}}),              \tag{2}
$$

where $h$ is any integer, the Frobenius exponents are nonnegative
integers, and all coefficients, including the tag $\gamma$, lie in $E$.
No nonvanishing-denominator promise is imposed. Consider the joint sum

$$
I=\sum_{v\in T}W(v)
          \sum_{u\in\mathbb F_b}P(u)\iota(1+\gamma uv).   \tag{3}
$$

Put

$$
S=r_P+r_++d(r_-+1).                                    \tag{4}
$$

There is a division-free arithmetic circuit, with topology independent
of the coefficient values, computing (3) in

$$
O\!\left((S+1)^3(w+\log(S+2))
       +d(r_P+r_++r_-+1)+w\right)                       \tag{5}
$$

field additions and multiplications. Integer exponent reductions and
circuit generation have polynomial bit cost in the displayed circuit
size and the exponent input lengths. With a polynomial-bit-cost field
representation, the bit complexity is polynomial in those quantities.

The explicit field is part of the input. Constructing it is a separate
charge, for example the coarse deterministic bound $2^{O(d)}$ by
exhaustive irreducible-polynomial search. An application whose existing
field does not contain $\mathbb F_{b^2}$ must charge the required
extension. No compositum or root family is supplied at unit cost.

## 2. Exact interpolation and the exceptional value

For an indeterminate $Z$,

$$
\sum_{u\in\mathbb F_b}\frac{P(u)}{Z+u}
       =\frac{P(Z)}{Z^b+Z}.                             \tag{6}
$$

Indeed, the denominator on the right has simple roots exactly at the
subfield elements, with derivative one. Its residue at $u$ is $P(u)$.
Both sides have those principal parts and vanish at infinity because
$\deg P<b$. Their difference is therefore zero.

Suppose first that $\gamma\ne0$, and set $z=(\gamma v)^{-1}$.
Away from a subfield pole, (6) gives

$$
\sum_u P(u)\iota(1+\gamma uv)
 =\frac{zP(z)}{z^b+z}
 =\gamma^bP(\gamma^{-1}v^{-1})
             \iota(\gamma^b+\gamma v^2),                \tag{7}
$$

where the last equality uses $v^b=v^{-1}$ on $T$.

There is at most one exceptional torus point. Define the unique field
element

$$
v_0^2=\gamma^{b-1}.                                    \tag{8}
$$

Squaring is an automorphism of $E$, so $v_0$ is obtained by an inverse
Frobenius power. A pole occurs on $T$ exactly when $v_0\in T$. In that
case $z_0=(\gamma v_0)^{-1}$ lies in $\mathbb F_b^\times$.

At $Z=z_0+X$, the right side of (6) expands as

$$
\frac{P(z_0+X)}{X+X^b}
  =\frac{P(z_0)}X+P'(z_0)+\epsilon_bP(z_0)+O(X),
\qquad
\epsilon_b=
\begin{cases}1,&b=2,\\0,&b>2.\end{cases}                \tag{9}
$$

For $b>2$, the next term after $X^{-1}$ in
$1/(X+X^b)$ has positive degree. For $b=2$, it has degree zero.
The regular part on the left of (6) is the sum with the singular
summand omitted. Since $\iota(0)=0$, this is exactly the required value
at the pole.

Consequently the first-stage function is the factorized expression in
(7), plus a point atom at $v_0$ with weight

$$
z_0\bigl(P'(z_0)+\epsilon_bP(z_0)\bigr),                \tag{10}
$$

when $v_0\in T$. Expression (7) itself is zero there. The derivative in
(10) is evaluated in linear time in $r_P$, after computing the factor
values: prefix and suffix products implement the product rule.
Only factors with $e_i=0$ have nonzero derivative.

If $\gamma=0$, the first-stage function is constant:

$$
\sum_{u\in\mathbb F_b}P(u)=[U^{b-1}]P(U)=:\tau_P.        \tag{11}
$$

Under (1), this coefficient is zero if $\Delta_P<b-1$, and is
$\prod_i b_i$ if $\Delta_P=b-1$. No large coefficient table is needed.

## 3. Why the second elimination stays small

On the torus, each factor in the substituted input remains a Frobenius
binomial:

$$
a_i+b_i\gamma^{-2^{e_i}}V^{-2^{e_i}}
 =a_i+b_i\gamma^{-2^{e_i}}V^{2^{e_i+w}}.                 \tag{12}
$$

Also, for every field element, including zero,

$$
\iota(x)=x^{2^{d+1}-3}
       =\prod_{j\in\{0,2,\ldots,d\}}x^{2^j}.            \tag{13}
$$

The set in (13) has $d$ elements, with the obvious one-element
interpretation when $d=1$. Thus each inverse factor introduces exactly
$d$ Frobenius binomials. The main term in (7), multiplied by (2), has
at most $S$ factors after this expansion.

The [Fermat-torus binomial circuit](fermat-torus-binomial-contraction.md)
contracts that complete product in the cubic term of (5). The
exceptional contribution is just
$W(v_0)$ times (10), after a torus-membership test. It needs one field
point evaluation, not a second table. The $\gamma=0$ branch uses (11)
times the same torus circuit applied to $W$ alone.

This is an exact closure statement for the representation consisting
of one factorized star, one possible point atom, and the degenerate
constant term. It does not assert closure as one ordinary rational
product without point corrections.

## 4. A coefficient-independent division-free circuit

The previous algorithm may branch on $\gamma=0$ and on $v_0\in T$.
Those branches can be removed explicitly over the stated field.
Let $Q=2^d$ and define

$$
\mathcal Z(x)=1+x^{Q-1},\qquad
a=\iota(\gamma),\qquad
v_0=(\gamma^ba)^{2^{d-1}},\qquad
z_0=\iota(\gamma v_0),                                  \tag{14}
$$

$$
\chi=\mathcal Z(1+v_0^{b+1}),\qquad
\delta=\mathcal Z(\gamma).                              \tag{15}
$$

Here $\mathcal Z$ is exactly the zero indicator on $E$. Define
$\widetilde W$ at every field value using (2) with $h$ replaced by its
nonnegative residue modulo $b+1$. It agrees with $W$ on $T$ and has
a defined value at zero.

Then the single formula

$$
\begin{aligned}
I={}&
\sum_{v\in T}\gamma^bW(v)P(av^{-1})
                    \iota(\gamma^b+\gamma v^2)\\
&+\chi\,\widetilde W(v_0)\,
       z_0\bigl(P'(z_0)+\epsilon_bP(z_0)\bigr)
 +\delta\,\tau_P\sum_{v\in T}W(v)                       \tag{16}
\end{aligned}
$$

holds for all coefficient values.

For $\gamma\ne0$, it is precisely Sections 2–3. For $\gamma=0$,
the first term is zero, $v_0=z_0=0$, and $\chi=0$, leaving (11).
Every inverse in (14) and (16) is replaced by (13). All zero indicators,
square roots, and membership indicators use fixed positive powers.
The two torus circuits and point evaluations therefore give a fixed
division-free circuit with the bound (5).

This circuit statement is field-valued. It does not automatically
extend zero indicators to arbitrary algebras with nilpotents.

## 5. A minimal failure of literal Frobenius-product closure

For the constant input $P=1$, ordinary rational-function elimination
already gives

$$
\sum_{u\in\mathbb F_b}\frac1{1+uV}
       =\frac{V^{b-1}}{1+V^{b-1}}.                     \tag{17}
$$

The right side has $b-1$ distinct simple poles, at
$\mathbb F_b^\times$, with no numerator cancellation.
A single Frobenius binomial $a+cV^{2^e}$ has at most one distinct root
over an algebraic closure. Therefore any ordinary rational
single-product identity for (17) using only such denominator
binomials needs at least $b-1$ denominator factors.

The first case is $b=4$, where one bilinear denominator produces
three distinct poles after elimination. The count is exponential in
$w$. It is a statement about literal rational single-product
representations, not a lower bound for arbitrary circuits or
finite-domain function representations.

The positive result avoids this growth by using the complementary
torus equation. There, the denominator becomes the single Frobenius
binomial $\gamma^b+\gamma V^2$, with at most one exceptional point.

## 6. The remaining multiplication-network gap

The degree promise in (1) is essential. For example, $P(U)=U^b$ agrees
with $U$ on the subfield, but substituting the unreduced polynomial in
(6) introduces a nonzero polynomial part. Reducing a general product
modulo $U^b-U$ can destroy its factorized representation. This theorem
does not provide that reduction at unchanged size.

The two domains must also match: the torus has order $b+1$ for the same
$b$ as the subfield. An integer digit interval and the Fourier torus of
an independently chosen transform modulus do not automatically satisfy
this condition.

Finally, the actual multiplication-gate phase $\zeta^{sxy}$ has an
exponent controlled by other integer labels. It is not the field
product $UV$ inside (3). A replacement of that gate must preserve all
joint labels and charge its representation before this theorem can be
used. No such replacement is proved here, and the complete factoring
exponent is unchanged.

## 7. Exact replay

Run:

    python tools/check_bilinear_cauchy_closure.py --output evidence/bilinear-cauchy-closure.json

The [checker](../tools/check_bilinear_cauchy_closure.py) compares both
the branching algorithm and the fixed circuit with complete joint
enumeration. The [evidence](../evidence/bilinear-cauchy-closure.json)
records 144 exact two-stage sums and 9,504 joint reference points.
It includes 48 cases with $b=2$, 18 zero tags, 96 exceptional torus
points, 44 nonzero point corrections, and 30 cases where the possible
exceptional field point is outside the torus.

Four exact formal rational identities verify the pole-count example.
Separate countercontrols detect omission of the $b=2$ correction and
use of an input polynomial of degree $b$. Four coefficient assignments,
including all zeros, execute the same 754 multiplication gates in a
fixed test configuration.
