# A sharp dyadic Bezout denominator for short arithmetic domains

## Result and scope

Let `m>=2`, let `A,B` be sets of `m` distinct positive integers, every one
congruent to one modulo four, and let `N>0` satisfy `N=3 mod 4`. Define

\[
 P(X)=\prod_{a\in A}(X-a),\qquad
 Q(X)=\prod_{b\in B}(N-bX).                       \tag{1}
\]

There is no factor pair `ab=N`, already by reduction modulo four, so the
rational gcd is one. Nevertheless, let `r` be the least nonnegative
integer for which

\[
 U(X)P(X)+V(X)Q(X)=2^r,
 \qquad U,V\in\mathbb Z_{(2)}[X].                \tag{2}
\]

Here `Z_(2)` consists of rational numbers with odd denominators. Then

\[
 m\le r\le2m-2,
 \qquad
 \boxed{r=2m-2\quad\text{if }m=2^s+1.}           \tag{3}
\]

The exact assertion holds for every choice of the nodes satisfying the
stated congruences, not only consecutive progressions. The number `2^r`
is the usual 2-adic reduced resultant; see the definition in the
[PARI library manual](https://pari.math.u-bordeaux.fr/dochtml/html-stable/usersch5.html).
The proof below computes its valuation for this family. The ordinary
resultant has valuation exactly `m^2`: since `P` is monic it is
`prod_(a,b)(N-ab)`, and every factor is two modulo four. In the sharp
family the optimal Bezout valuation is only `2m-2`, but both grow
exponentially in the binary domain length.

Taking `A=B={1,5,...,1+4(m-1)}` and `N=3` gives binary-described
arithmetic domains of encoding length `O(log m)`. Thus the necessary
dyadic denominator in a **single global polynomial Bezout certificate**
can have an exponential numerical exponent, even though the absence of
factor pairs has an elementary short certificate.

This scopes a reconstruction strategy. It is not a lower bound on
factoring, structured root separation, rational gcd decision, or algorithms
that retain enormous valuations in binary. The factoring exponent is unchanged.

## 1. The reduced inverse determines the minimum denominator

Because `P` is monic and `gcd(P,Q)=1`, there is one polynomial `V0` of
degree below `m` satisfying

\[
 V_0Q=1\pmod P.                                  \tag{4}
\]

If (2) holds, reduce `V` modulo `P`. Monic polynomial division preserves
`Z_(2)` coefficients, and uniqueness in (4) gives `V mod P=2^r V0`.
Conversely, if every coefficient of `2^r V0` belongs to `Z_(2)`, divide
`2^r-(2^r V0)Q` by the monic polynomial `P`. The quotient is in
`Z_(2)[X]` and gives (2). Thus `r` is exactly the largest dyadic
denominator exponent among the coefficients of `V0`.

This also rules out lowering `r` by using higher-degree Bezout
polynomials: their reduction has the same unique inverse.

## 2. A regular local quotient

Substitute `X=1+2Z`. The congruence assumptions give

\[
 P(1+2Z)=2^m\widetilde P(Z),\qquad
 \widetilde P(Z)=\prod_a\left(Z-\frac{a-1}{2}\right),
\]

\[
 Q(1+2Z)=2^mR(Z),\qquad
 R(Z)=\prod_b\left(\frac{N-b}{2}-bZ\right).        \tag{5}
\]

Every root of `P̃` is even, while both `(N-b)/2` and `b` are odd. Hence

\[
 \widetilde P(Z)\equiv Z^m\pmod2,\qquad
 R(Z)\equiv(1+Z)^m\pmod2.                        \tag{6}
\]

Multiplication by `R` on the free rank-`m` module
`Z_(2)[Z]/(P̃)` has odd determinant: modulo two it multiplies by a unit
in `F_2[Z]/(Z^m)`. Its inverse therefore has coefficients in `Z_(2)`.
Let

\[
 S(Z)=\sum_{i=0}^{m-1}s_iZ^i,
 \qquad SR=1\pmod{\widetilde P}.                 \tag{7}
\]

All `s_i` have nonnegative 2-adic valuation. Reducing (7) modulo two
gives the complete parity profile

\[
 s_i\equiv {m+i-1\choose i}\pmod2.               \tag{8}
\]

In particular `s_0` is odd. For `m>=2`, the top coefficient is even:

\[
 s_{m-1}\equiv{2m-2\choose m-1}=0\pmod2.
\]

Every positive central binomial coefficient is twice an integer.

## 3. Affine normalization gives the exact bound

Equations (4)--(7) give

\[
 V_0(1+Y)=2^{-m}S(Y/2)
       =\sum_{i=0}^{m-1}2^{-(m+i)}s_iY^i.         \tag{9}
\]

Translation by one is an integral invertible change of coefficient basis.
It preserves the largest dyadic denominator exponent. Therefore

\[
 r=\max_{0\le i<m,\ s_i\ne0}\bigl(m+i-v_2(s_i)\bigr). \tag{10}
\]

The constant coefficient gives `r>=m`. Every degree below `m-1` has
denominator exponent at most `2m-2`. The top coefficient is even, so it
has the same upper bound. This proves the inequalities in (3).

Now let `m=2^s+1`. At degree `m-2`, equation (8) becomes

\[
 s_{m-2}\equiv{2^{s+1}-1\choose2^s-1}=1\pmod2.   \tag{11}
\]

The parity follows from
`(1+Z)^(2^(s+1)-1)=1+Z+...+Z^(2^(s+1)-1)` over `F_2`.
This coefficient attains denominator exponent `m+(m-2)=2m-2` in (9),
proving equality. No spacing assumption beyond the residue classes is used.

## 4. Consequences for compressed elimination

The [modular interval circuits](interval-packet-circuits.md) and their
[bilinear extension](bilinear-pole-packets.md) construct long product
polynomials at polynomial precision. Affine arithmetic domains, including
this family, retain that inexpensive generation. A raw modular calculation
of the nonzero scalar in (2), however, returns zero at every precision at
most `r`.

In the sharp family, `r=2^(s+1)` while the domain endpoints have `O(s)` bits.
Raising ordinary absolute precision until the global Bezout scalar becomes
nonzero is not a polynomial-cost coprimality certificate on this family.
This concerns the scalar representation: `r` itself has a short binary
encoding, and a structured normalization may store it without allocating
`2^r` bits.

The actual rational gcd is easy to decide here: `ab=1 mod4` and `N=3 mod4`
exclude every common root. A proposed output-sensitive gcd algorithm can
use that arithmetic separation. The theorem prevents assuming that all
empty-output instances have polynomially bounded global Bezout denominators.

Compared with the earlier artificial-common-factor example, this result
determines the optimal denominator exponent, permits arbitrary Bezout
degrees, and provides a matching upper bound. It does not force every
structured extraction route to pay that denominator.

### Irregular nodes can change the finer denominator profile

For consecutive progressions, the [balanced-node interpolation theorem](balanced-bezout-denominator.md)
proves `r=m+v2((m-1)!)` for every size. That stronger formula does not hold
for general node sets. For example,

```text
m=4, N=3,
A={1,5,9,13}, B={1,5,9,17}
```

has `r=6`, whereas `m+v2((m-1)!)=5`. Its normalized inverse `S` is
`3+6Z+6Z^3 modulo 32`. This control preserves every hypothesis of (3).

## 5. Independent exact controls

The [checker](../tools/check_interval_bezout_precision.py) computes the
reduced inverse in two ways. Scaled Lagrange interpolation uses literal
values of `Q(a)` and `P'(a)`. A separate Hensel iteration in the normalized
quotient from (5) uses polynomial arithmetic and the initial parity inverse.
Every returned denominator exponent agrees. Twelve small cases additionally
use exact rational coefficient lists.

The retained [evidence](../evidence/interval-bezout-precision.json) includes
54 independent comparisons, 18 sharp power-plus-one cases, irregular node
sets, the counterexample above, and complete normalized parity profiles.
The sharp progression cases have

| `m` | Exact `r` |
| --- | --- |
| 2 | 2 |
| 3 | 4 |
| 5 | 8 |
| 9 | 16 |
| 17 | 32 |
| 33 | 64 |

The implicit `m=2^200+1` record applies the proved formula to get
`r=2^201`; it does not expand that polynomial or compute a dense inverse.

```sh
python3 tools/check_interval_bezout_precision.py \
  --output evidence/interval-bezout-precision.json
```

No new Lean formalization or stronger factoring bound is claimed.
