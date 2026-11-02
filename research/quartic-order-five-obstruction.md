# Returning order-five channels in additive-color quartic fusion

## Question, scope and disposition

The [quartic construction](quartic-charge-fusion.md) changes separate color
parity and reduces the generic order-four multiplicity. Can its higher
channels pay the degree-four restoration and cross `1/sqrt(3)`? The next
test is a generic order-five coefficient, allowing independent spatial
quartic weights rather than fitting a larger alphabet or a complete profile.

**In the restricted family below, at least `(q-1)(q-4)` order-five channels
return in every axis for `q=2^s>=8`. For every `q>=4096`, the common-scale,
unweighted three-axis profile criterion has restored pressure at least
`7/12>1/sqrt(3)`, for every common physical type.**

Thus the earlier q4096 abstract profile is incompatible with this tensor
family. Further growing-alphabet or full-profile expansion of this family
is stopped. The argument does not classify other tensors, nonvacuum
interactions, lower-valuation quartics, weighted axis scales, or all smaller
alphabets. The complete factoring exponent remains `1/2`.

## 1. The spatial freedom being allowed

Keep the pair kernel `ell^2 K`, additive colors, and distinct zero-XOR
quartic creators of the construction, with `ell=x+y+z`. For each spatial
four-set `S` allow its own weight

```math
H_S=(xyz)^4[\ell^4+xyz(a_Sx+b_Sy+c_Sz)],\qquad
a_S+b_S+c_S=1.                                      \tag{1}
```

Coefficients may lie in any characteristic-two field. This is the entire
degree-16 family divisible by `(xyz)^4` with the three prescribed leading
coefficients and zero endpoint: subtracting `ell^4` from the degree-four
quotient gives a polynomial vanishing on all coordinate planes, hence
`xyz` times a linear form. The endpoint gives the displayed sum constraint.
All construction and endpoint proofs remain valid, as does the generic
prefix `0^1,2^(q-1),4^(2(q-1))`.

Higher even creation gates of valuation at least six in each cut cannot
change the calculation through order five. We do not allow a lower-order
interaction and silently assume the same prefix.

## 2. A lowest coefficient that no allowed repair reaches

Fix the x cut and remove wholly one-side gates by the existing invertible
creation operations. Nonzero XOR charge is conserved; occupied-count
parity further separates the even and odd blocks. Fix a nonzero charge
`h`. There are

```math
r=(q-2)/2
```

unordered distinct nonzero color pairs with XOR `h`. Each has twelve
ordered placements on the four vertices of one side. Up through order
five, the even block consists of these two-occupied words: a four-occupied
word has no charge-`h` empty partner, and pairing it with a nonempty even
word costs at least order six.

Write a cross-pair entry as `x^2 A+x^3 B_uv+O(x^4)`. Directly from the
kernel,

```math
A=y^2z^2(y+z)^2,\qquad
\frac{B_{uv}}A=
 (1+u_z+v_z)\frac y{z^2}
 +(1+u_y+v_y)\frac z{y^2}.                        \tag{2}
```

Here `u_y,u_z,v_y,v_z` are spatial vertex bits; sums are in the coefficient
field. Two cross-pairs have no `x^5` term of z-degree below two. Their
coefficient at `x^5 y^9 z^2` is the sum of the four z-bits. This entry is
zero unless the row and column use the same unordered color pair. An
allowed quartic has no `x^5` coefficient below z-degree five, regardless
of its coefficients in (1).

Consequently the lowest-z coefficient of the two-occupied order-five
matrix is block diagonal with `r` identical blocks

```math
J=u\mathbf1^T+\mathbf1 v^T.                       \tag{3}
```

Both parity values occur among spatial placements, so `rank J=2` over
every characteristic-two field. The complete coefficient matrix has
rank `2r`. The order-four matrix is `A^2` times the all-ones matrix.
Eliminating one order-four pivot replaces an order-five entry `E_ij` by

```math
E_{ij}+E_{i0}+E_{0j}+E_{00}.                       \tag{4}
```

These row and column projections each lose at most one rank. Thus the
lowest-z coefficient after elimination has rank at least `2r-2=q-4`.
All other z-orders are higher, so its nonzero minor remains nonzero over
the generic field `F(y,z)`: this uses the **initial coefficient in z**,
not the generally false assertion that any coefficient bounds polynomial
matrix rank. The order-four pivot's common factor cancels in (4).

There are therefore at least `q-4` order-five Smith channels in this
charge, hence `(q-1)(q-4)` in each axis. Relabeling axes proves all three
cases. The independent q8 formal control obtains ranks six before the
projection and four after it, matching the lower bound. This does not
assert a complete generic profile or equality for every higher layer.

## 3. A dimension-and-degree tail

Let `n=q^4`, `d=q-1`, `N=d(q-4)`, and `Z_a(t)` be the complete generic
Smith generating polynomial of one axis, for `0<t<1`. Select the known
`1+3d+N=m=q^2-2q+2` channels, with orders zero, two, four and five.
Every matrix entry has degree at most 32, and the endpoint is invertible,
so the sum of all `n` Smith orders is at most `32n`.

For `q>=8`, `33m<=n` and `2m<=n`. The remaining orders have mean at most
`32n/(n-m)<=33`. Convexity of `t^v` gives the lower bound

```math
Z_a(t)\ge L_q(t):=
1+d t^2+2d t^4+N t^5+\tfrac12 q^4t^{33}.          \tag{5}
```

This is a tail estimate for actual profiles, rather than a proposed
abstract profile. We use a deliberately coarse order 33; its optimization
is unnecessary for the target exclusion.

## 4. A short all-scale pressure inequality

For `q>=4096`, put `beta=32/7`, `gamma=7/3`. We prove

```math
L_q(t)>(1+d t^\beta)^\gamma,\quad 0<t<1.           \tag{6}
```

For every `x>0`,

```math
(1+x)^{7/3}\le 1+\tfrac73x+\tfrac{14}9x^2+x^{7/3}. \tag{7}
```

Subtract the right side: its third derivative is a positive constant
times `(1+x)^(-2/3)-x^(-2/3)<0`; the function and its first two derivatives
vanish at zero. Three integrations prove (7).

For `t<=3/4`, use `x=d t^beta`. The two small-order terms in (5) supply
the linear term because `t^2+2t^4>=3t^beta`. Also

```math
\tfrac{14}9d^2t^{64/7}<\tfrac12Nt^5,
```

since `d/(q-4)<=257/256` and the exact rational check is
`(28*257/(9*256))^7*(3/4)^29<1`.
Allocate the remaining half of `Nt^5` and the entire tail to `x^(7/3)`.
Here `N>=(31/32)q^2`, so weighted Young's inequality gives

```math
\tfrac12Nt^5+\tfrac12q^4t^{33}
 \ge Cq^{7/3}t^{29/3},\qquad
C^6=\frac{6^6}{5^5}\left(\frac{31}{64}\right)^5\frac12
 >\left(\frac34\right)^6.
```

Since `t<=3/4`, this exceeds `q^(7/3)t^(32/3)`, which bounds
`d^(7/3)t^(32/3)=x^(7/3)`. This proves (6) on the lower interval.

For `t>=3/4`, the tail alone exceeds `q^(7/3)`, because
`q^(5/3)>=2^20` and `3^33>2^47`. Meanwhile
`1+d t^beta<=q`. This proves (6) on the upper interval. All constant
comparisons are exact rational or integer checks; no fitted polynomial
certificate or floating minimization is required.

## 5. Entropy restores the actual degree cost

For any common physical type `p` with entropy `H(p)>0` and occupied mass
`E`, the entropy bound with vacuum weight one and nonzero-color weight
`t^beta` is

```math
H(p)\le\log(1+d t^\beta)-\beta E\log t.
```

Apply (5)--(6) to every axis. Since `gamma beta=32/3`, the degree-four
objective satisfies

```math
\Phi_p(t)=
\frac{\frac13\sum_a\log Z_a(t)-\frac{32}3 E\log t}{4H(p)}
 >\frac\gamma4=\frac7{12}.
```

Taking the required infima gives a floor at least `7/12`.
The exact comparison `3*7^2>12^2` places it above `1/sqrt(3)`.
It excludes the listed q4096 relaxation and all growing-alphabet extensions
under (1), not just the uniform original quartic. To continue constructively
one must change an assumption that creates this low second-axis jet, rather
than enlarge the profile computation.

## Verification and credit

`make verify-quartic-order-five` checks the formal coefficient identity on
all three cuts, q8 coefficient ranks before/after projection, and the exact
rational allocations. The [receipt](../evidence/quartic-order-five.json)
pins its code and the Gaussian/quartic sources. The generic rank and entropy
arguments above are ordinary author-checked proofs, not independent review.
The construction retains the worker-4 Gaussian credit; the complementary
[global divided-square identity](charge-creation-square.md) is Sol worker
2's contribution. No new Lean or stronger complete factoring bound is claimed.
