# Homogeneous modular duals at every label width

## Result and scope

The normalized characteristic-two crossing admits a uniform homogeneous
dual construction for every binary label width. Its state labels have
one weight `h`, and all dual operators have one weight `j<=h`. Both
weights are quadratic in the label width, and

```math
h-j=O(\sqrt h\log h).
```

The construction uses a nonsingular minor of one homogeneous vacuum
pairing. Its existence follows from the already proved filtered source
duality and rank subadditivity. It does not assume that the same-weight
homogeneous Gram matrix is nonsingular, nor that finite rank patterns
continue indefinitely.

The resulting family has the existing uniform contraction coefficient
`rho=1/sqrt(2)`. It gives another complete ordinary route to exponent
`1/2`; it does **not** prove an exponent above `1/2`. The structural
gain is that every dual has a single weight and every angle is a scalar
binomial coefficient. No restriction to power-of-two state energies is
needed. A second result below gives a polynomial-time complete-profile
Pfaffian for raw homogeneous operator parameters, including `j<h`.

The source operators and crossing identities are Anthropic's Appendix B
at
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The normalized quotient, its effective transfer, and the earlier
complete-profile Pfaffian are in
[PR #22, pinned at `d94e4a99`](https://github.com/1337isnot1337/integer-factoring/tree/d94e4a99fec453768d69ae083038ba65cf31ae9e).
The [symbol note](homogeneous-central-modular-symbols.md) supplies the
individual angle formula and states the relationship to that prior work.
This selection argument, profile variant, and checker used OpenAI
assistant support. No new Lean proof is asserted.

## 1. The filtered duality contains a large homogeneous minor

Work in

```math
W=\mathbf F_2[t_1,t_3,t_5,\ldots]/(t_1^2,t_3^2,t_5^2,\ldots),
\qquad \deg t_k=k.
```

Write `q(h)=dim W_h`. Its basis consists of partitions of `h` into
distinct odd parts. For `0<=j<=h`, let `F_(j,h)` be the rectangular
matrix with rows indexed by this basis of `W_j`, columns by the basis
of `W_h`, and entries

```math
F_{j,h}(\eta,\nu)=[1]t^\eta[-h]t^\nu.                \tag{1}
```

The matrix is over `F_2`; no division by a factorial occurs.
Let `Z_(h,j)(mu,eta)` be the coefficient of `t^eta` in the weight-`j`
part of the original filtered source dual `Theta(t^mu)`. The source
duality and the valid reduction of the operator parameter give

```math
I_{q(h)}=\sum_{j=0}^h Z_{h,j}F_{j,h}.                \tag{2}
```

All Fock norms of regular characteristic-two monomials are one, so
the right side is exactly the identity. Equation (2) is a finite matrix
identity, not a statement that the individual blocks are invertible.
Taking ranks gives

```math
q(h)\le\sum_{j=0}^h\operatorname{rank}F_{j,h}.        \tag{3}
```

Fix `d=2^lambda`, with `lambda>=1`, and choose the least positive `h`
such that

```math
q(h)\ge(h+1)d.                                      \tag{4}
```

There is some `j<=h` with `rank F_(j,h)>=d`. Choose a nonsingular
`d` by `d` submatrix `B`, using actual parameter rows `eta_1,...,eta_d`
and state columns `nu_1,...,nu_d`. Set

```math
b^v=t^{\nu_v},\qquad
a^x=\sum_{k=1}^d(B^{-1})_{xk}t^{\eta_k}.             \tag{5}
```

Then every `a^x` is homogeneous of the same weight `j`, every `b^v`
has weight `h`, and

```math
a^x[-h]b^v=[x=v]1.                                 \tag{6}
```

Since `d>=2` and `q(0)=1`, a selected block always has `j>=1`.
The algorithm can select the first valid `j` in any fixed computable
order. The implementation prefers fewer set bits in `h+j`, then a
larger `j`. This favors fewer nonzero angles in finite examples; no
asymptotic advantage is assumed for that ordering.

## 2. A quadratic search bound and complete construction cost

The explicit characteristic-two block labels in
[modular-crossing.md](modular-crossing.md) give
`q(8m^2)>=2^m`. Take `m=2lambda+16`. For every `lambda>=1`,

```math
2^m\ge(8m^2+1)2^\lambda.                            \tag{7}
```

At `lambda=1` this is a direct integer inequality. When `lambda`
increases by one, the left side after division by `2^lambda` doubles,
whereas `8(2lambda+16)^2+1` increases by a factor less than two.
Thus (7) holds by induction. The least weight in (4) consequently
satisfies

```math
h\le8(2\lambda+16)^2=O(\lambda^2).                  \tag{8}
```

Compute the regular-partition count table through the right side of
(8), and find the first weight satisfying (4). Enumerate the regular
partitions through `h`. Their number is `2^O(sqrt(h))=2^O(lambda)`.
The entries (1) are normal-order coefficients; the determinant formula
for characteristic two in the modular note evaluates them within the
same bound, and in fact polynomially in the number and sizes of their
parts when the output is the vacuum.

Constructing all the rectangular blocks, finding independent rows and
columns, and inverting `B` use a fixed polynomial in the enumerated
matrix dimensions and the integer index lengths. The total cost is
`2^O(lambda)` bit operations. All retained coefficients are bits.
There is no uncharged search over all label subsets and no choice of
input-dependent advice.

## 3. The selected weights have the same leading scale

Use the distinct-odd-part asymptotic

```math
q(t)\sim
\frac{\exp(c\sqrt t)}{2\,24^{1/4}t^{3/4}},
\qquad c=\pi/\sqrt6.                               \tag{9}
```

The source and its use for the modular family are recorded in
[PR #22, Section 7](https://github.com/1337isnot1337/integer-factoring/blob/d94e4a99fec453768d69ae083038ba65cf31ae9e/research/past-half-pieces.md#7-the-full-half-rate-coefficient-over-characteristic-two).
Primary statements are Remark 3.3 of Bridges, Franke, and Garnowski,
[*Asymptotics for the twisted eta-product and applications to sign changes
in partitions*](https://doi.org/10.1007/s40687-022-00355-x), and the
self-conjugate-partition calculation in Section 5.1 of Wildon,
[*Counting partitions on the abacus*](https://arxiv.org/html/math/0609175v1).

Minimality in (4) gives `q(h-1)<h 2^lambda`. Applying (9) at `h-1`,
and using (4) at `h`, yields

```math
c\sqrt h=\lambda\ln2+O(\log(h+2))
         =\lambda\ln2+O(\log(\lambda+2)).            \tag{10}
```

This argument does not assume monotonicity of `q`. The finite search
bound (8) controls the logarithmic term.

Since `rank F_(j,h)>=d`, its number of rows satisfies `q(j)>=2^lambda`.
Equation (9), together with `j<=h`, then gives

```math
c\sqrt j\ge\lambda\ln2-O(\log(h+2)).
```

Subtracting this from (10) proves

```math
0\le\sqrt h-\sqrt j=O(\log h),\qquad
h-j=O(\sqrt h\log h).                               \tag{11}
```

The matrix selection therefore loses only polynomially many labels
relative to the full state space in the logarithmic parameter. It
does not lose a fixed fraction of the leading label entropy.

## 4. Crossing identity and uniform half-rate widths

Use (5) in the original source construction, with state labels
`t^(nu_y), t^(nu_u)`. Its hypothesis permits `a` of weight at most
`h`; hence `j<h` causes no change to the exchange or sum identity.
For `-2h<=s<=2h`, put

```math
A_s(x,y,u,v)=[t^{\nu_u}]
 \sum_{i=0}^{h+j}\binom{h+j}{i}
 a^x[s-h+i]b^v[h-s-i]t^{\nu_y}.                      \tag{12}
```

The binomial coefficient is reduced modulo two. Formula (6) and the
source identities give

```math
\sum_s A_s(x,y,u,v)=[x=v][y=u].                      \tag{13}
```

The nonzero angles are precisely the submasks of the binary integer
`h+j`. If an actual operator coefficient vanishes, that can remove
further terms; no such removal is used in the width estimate.

The two operator orders factor the piece through intermediate energies
`s+i` and `2h-s-i`. With `q(E)=0` for negative `E`, valid channel
counts are

```math
R_s=\sum_{i\subseteq h+j}q(s+i),\qquad
C_s=\sum_{i\subseteq h+j}q(2h-s-i).                  \tag{14}
```

For a nonzero piece, `s+h+j>=0` and `2h-s>=0`. Equation (9) bounds
the two sums by

```math
\begin{aligned}
\log R_s&\le c\sqrt{s+h+j}+O(\log h),\\
\log C_s&\le c\sqrt{2h-s}+O(\log h).
\end{aligned}
```

There are only `O(h)` terms in either sum, so these bounds are uniform
in the shift. Concavity and `j<=h` give

```math
\begin{aligned}
\log_2R_s+\log_2C_s
&\le\frac c{\ln2}\sqrt{2(3h+j)}+O(\log h)\\
&\le2\sqrt2\,\lambda+O(\log(\lambda+2)).             \tag{15}
\end{aligned}
```

Ceilings change the sum by at most two. A zero piece uses one zero
channel, and widths above `2lambda` may use the trivial flattening
factorization. Thus (15) is the effective piece condition with
`rho=1/sqrt(2)`.

There are `4h+1=O(lambda^2)` pieces. Their coefficients, both factors,
and the dual descriptions have the `2^O(lambda)` construction cost
from Section 2. Their scale is one in `F_2`, so they extend to every
commutative `F_2`-algebra. The established
[modular transfer](scale-one-modular-transfer.md) therefore applies and
retains the factoring time `exp(O(n/sqrt(log n)))`.

## 5. Complete profiles for raw homogeneous parameters

For exact finite calculations, intermediate spaces need not be listed.
Let `a=t^eta` be one regular monomial of weight `j<=h`, and let
`u,v,y` be regular partitions of weight `h`. Use four vertex blocks
`U_t`, `A_n`, `B_m`, `Y_t`, indexed by those parts. The alternating
matrix over `F_2[Z]` has zero within-block entries and these cross-block
entries:

| Pair | Polynomial entry |
| --- | --- |
| `U_t,A_n` | `[n-1 subseteq t-1] Z^t(1+Z)^n` |
| `A_n,Y_t` | `[n AND t=1] * (1+Z)^n` |
| `A_n,B_m` | `[n AND m=1] Z^n` |
| `U_t,B_m` | `[m-1 subseteq t-1] * (1+Z)^m` |
| `B_m,Y_t` | `[m AND t=1] Z^t(1+Z)^m` |
| `U_t,Y_r` | `[t=r] Z^r` |

Here subset means binary-digit containment. Then

```math
\operatorname{Pf}M(Z)
 =Z^{h+j}\sum_s [t^u]Q_s(a,t^v)t^y\,Z^s.            \tag{16}
```

This is a raw-parameter variant of the earlier PR #22 four-block
formula, whose first parameter is the filtered source dual. The
ordinary normal-order matching proof is the same. The raw contraction
between odd parts `n,m` is particularly simple:

```math
\sum_{k\text{ odd}>0}
 [\, (n-1)\mathbin{\&}(k-1)=0\,]
 [\,m-1\subseteq k-1\,]Z^{-k}
 =[n\mathbin{\&}m=1]
   \frac{Z^n}{(1+Z)^{n+m}}.                         \tag{17}
```

The disjointness condition in (17) uses the half-indices, equivalently
`n AND k=1`; thus the first bracket can also be written
`[(n-1) AND (k-1)=0]` for these odd integers. Write
`n=2a+1`, `m=2b+1`, `k=2c+1`. If `a AND b` is nonzero, the sum
vanishes. Otherwise the allowed `c` contain `b` and avoid `a`.
Multiplying the available binary factors gives

```math
Z^{-2b-1}\prod_{r:\,\operatorname{bit}_r(a\lor b)=0}
 (1+Z^{-2^{r+1}})
 =\frac{Z^{2a+1}}{(1+Z)^{2a+2b+2}},
```

as a Laurent series at infinity. This proves (17).

Scale each `A_n` vertex by `(1+Z)^n`, each `B_m` by `(1+Z)^m`,
and each `Y_t` by `Z^t`. Formula (17) becomes the polynomial
`A_n,B_m` entry above; the other five entries are the raw creation,
annihilation, and unchanged-state pairings with the same scalings.
The total scale is `Z^h(1+Z)^(h+j)`. The source angle sum for a
homogeneous parameter contributes `Z^(-j)(1+Z)^(h+j)`, proving the
offset `Z^(h+j)` in (16).

The matrix has `O(sqrt(h))` vertices. Every entry has degree at most
the sum of its endpoint weights, so every principal Pfaffian has
degree at most `3h+j`. Fraction-free elimination, or a polynomial
determinant followed by the characteristic-two square root, computes
the whole profile in `O(h^4)` bit operations. The final square root
means retaining the even powers of a polynomial known to be a square;
it is not a matrix square root. For a dual in (5), add the profiles
of its at most `d` constituent monomials.

## 6. Evidence and the remaining rank question

Run:

```sh
python3 tools/check_homogeneous_central_selection.py --canonical-bits 1 \
  --output evidence/homogeneous-central-selection.json
```

The checker verifies the source matrix identity (2), the rank-sum
inequality (3), actual selected inverse minors, and the raw profile
formula (16). The independent profile reference uses the original
integer normal-order recursion and explicitly sums intermediate states.
It does not use the Pfaffian blocks to obtain the expected answer.

It also builds every complete piece on selected finite label sets and
checks every entry of the summed crossing. The canonical width-one
case uses the least weight from (4), rather than a convenient fixed
weight. A larger optional run uses `--canonical-bits 2`.

The recorded file is
[homogeneous-central-selection.json](../evidence/homogeneous-central-selection.json).
These checks verify the formulas and algorithms. The all-width
existence and half-rate bound follow from (2)--(15), not from the
recorded finite ranks.

This construction leaves actual complete-piece ranks open. A separate
numerical pattern for same-weight Gram ranks is not used anywhere in
the proof. A smaller angle support, an invertible selected minor, or
a finite rank drop alone does not establish `rho<1/sqrt(2)`.
