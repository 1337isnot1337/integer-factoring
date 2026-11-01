# A bilinear Fourier phase that defeats the naive Cauchy substitution

## Result

There is a concrete twenty-term counterexample to extending the
[bounded-degree Cauchy closure](bilinear-cauchy-two-stage-closure.md)
by inserting an encoded integer multiplication phase and leaving its
formal polynomial unreduced.

The correct joint sum is $1$. The naive substituted sum is $\omega^2$,
where $\omega^2+\omega+1=0$. Their nonzero difference is $\omega$.
The error is exactly the omitted polynomial part from an input of
degree four when the summation subfield has size four.

This is a failure of that proposed extension. It does not contradict
the proved bounded-degree theorem or rule out a different compact
representation of the coupled phase.

## 1. An exact field representation of the integer phase

Work in $E=\mathbb F_{16}$, with subfield
$\mathbb F_4=\{0,1,\omega,\omega^2\}$ and torus

$$
T=\{v\in E:v^5=1\}.
$$

Encode the integer $x=x_0+2x_1$, with $x_0,x_1\in\mathbb F_2$,
by the field element $U=x_0+x_1\omega$. The coordinate functions are

$$
L_0(U)=\omega^2U+\omega U^2,\qquad
L_1(U)=U+U^2.                                           \tag{1}
$$

They equal $x_0$ and $x_1$, respectively, on the subfield. Define

$$
P(U,V)=
\bigl(1+(1+V)L_0(U)\bigr)
\bigl(1+(1+V^2)L_1(U)\bigr).                            \tag{2}
$$

When a coordinate bit is zero, its factor is one; when it is one,
its factor is $V$ or $V^2$. Consequently

$$
P(x_0+x_1\omega,V)=V^{x_0+2x_1}.                        \tag{3}
$$

For $V=\zeta^y$ and a primitive fifth root $\zeta$, this is exactly
$\zeta^{xy}$. Thus both the integer multiplier and the frequency label
vary jointly. The phase has not been replaced by the field product $UV$.

The polynomial in (2) is a compact product of two affine linearized
factors, but its formal $U$-degree is four for every $V\in T\setminus\{1\}$.
Its leading coefficient is

$$
a_4(V)=[U^4]P(U,V)=\omega(1+V)(1+V^2).                  \tag{4}
$$

This is nonzero at all those four torus points.

## 2. The tempting substitution and its exact defect

Use the zero-preserving inverse $\iota(0)=0$ and
$\iota(z)=z^{-1}$ for $z\ne0$. The actual inner contraction is

$$
A(V)=\sum_{u\in\mathbb F_4}P(u,V)\iota(1+uV).           \tag{5}
$$

Dropping the degree restriction in the bounded-degree formula would
give

$$
B(V)=P(V^{-1},V)\iota(1+V^2).                           \tag{6}
$$

To compute the missing term, reduce the input polynomial correctly:

$$
\overline P(U,V)=P(U,V)+a_4(V)(U^4+U).                 \tag{7}
$$

It has degree at most three and agrees with $P$ at every subfield
element. For $V\ne1$, the denominator in (5) has no subfield pole:
the groups $\mathbb F_4^\times$ and $T$ have coprime orders three and
five. Applying the valid degree-bounded formula to (7) gives

$$
\begin{aligned}
A(V)+B(V)
 &=a_4(V)\frac{V^{-4}+V^{-1}}{1+V^2}\\
 &=\frac{a_4(V)}V,                                     \tag{8}
\end{aligned}
$$

where $V^5=1$ was used in the last step.

At $V=1$, equation (2) is the constant polynomial one. The true
zero-safe inner sum and (6) are both zero; $a_4(1)=0$ as well.
Thus (8) holds at every torus point. The sole denominator-zero
assignment is $(u,V)=(1,1)$, and it causes no ambiguity in this example.

Summing the defect gives

$$
\sum_{V\in T}(A(V)+B(V))
 =\omega\sum_{V\in T}(V^{-1}+1+V+V^2)
 =\omega.                                              \tag{9}
$$

The three nonconstant character sums vanish, and the odd torus size
is one in characteristic two.

## 3. The two joint values

In fact,

$$
\sum_{V\in T}A(V)=1,\qquad
\sum_{V\in T}B(V)=\omega^2.                             \tag{10}
$$

Here is a direct check of the first value. The $u=0$ term contributes
one. For $u=1$, pair each nontrivial torus element with its inverse:
each pair contributes one to $\sum V/(1+V)$, and there are two pairs,
so that contribution is zero.

For $a\in\{\omega,\omega^2\}$, ordinary torus summation gives

$$
\sum_{V\in T}\frac{V^k}{1+aV}
       =\frac{a^{5-k}}{1+a^5},\qquad 1\le k\le4.        \tag{11}
$$

For example, start with the logarithmic derivative of $Z^5+1$ to
obtain the $k=0$ sum, and use
$aV^k/(1+aV)=V^{k-1}+V^{k-1}/(1+aV)$ successively.
The $u=\omega$, $k=2$ contribution is $\omega^2$, and the
$u=\omega^2$, $k=3$ contribution is also $\omega^2$.
They cancel. This proves the first equality in (10); the second follows
from (9) and $1+\omega=\omega^2$.

## 4. Why width two is the first degree defect in this encoding

For a binary field of size $b=2^w$, choose trace-dual coordinate
functionals

$$
L_j(U)=\mathrm{Tr}_{\mathbb F_b/\mathbb F_2}(\beta_jU),
\qquad
\mathrm{enc}(U)=\sum_{j=0}^{w-1}2^jL_j(U).
$$

The trace pairing is nondegenerate: its trace polynomial has degree
$b/2<b$ and is nonzero, and multiplication by a nonzero field element
permutes the field. Hence such dual coordinates exist.

Each $L_j$ has formal degree $b/2$ with nonzero leading coefficient.
The exact phase product is

$$
V^{\mathrm{enc}(U)}
 =\prod_{j=0}^{w-1}
       \left(1+(1+V^{2^j})L_j(U)\right)                 \tag{12}
$$

as a function on the subfield. For $V\ne1$ on the torus of odd order
$b+1$, all $1+V^{2^j}$ are nonzero. The formal $U$-degree of the product
in (12) is therefore exactly

$$
\frac{wb}{2}.                                          \tag{13}
$$

This is below $b$ at width one, equals $b$ at width two, and exceeds it
thereafter. The example above is the first width at which this specific
compact encoding violates the degree hypothesis.

A reduction modulo $U^b-U$ always exists. What remains to be proved for
a fast general contraction is a controlled representation and cost for
that reduction and for the resulting jointly varying coefficients.
Equation (13) is not a lower bound on every possible representation.
The positive restricted theorem and the full multiplication-network
gap both remain as stated.

## 5. Exact replay

Run:

    python tools/check_cauchy_bilinear_phase_counterexample.py --output evidence/cauchy-bilinear-phase-counterexample.json

The [checker](../tools/check_cauchy_bilinear_phase_counterexample.py)
verifies all twenty phase encodings, the formal leading coefficient,
the exact reduction, all five inner sums and defects, and both joint
values. It uses the explicit binary polynomial $X^4+X+1$ for
$\mathbb F_{16}$. The [evidence](../evidence/cauchy-bilinear-phase-counterexample.json)
records every field value and its polynomial bit representation.
