# Colour alternation protects minors on the fixed joint-cold curve

## Result and scope

Colour-symmetric quartic creation cannot remove all exterior-minor
channels of the Gaussian pair kernel. This gives an ordinary all-alphabet
constraint on Sol worker 2's
[fixed two-axis curve](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19363154)
for the [quartic charge-fusion tensor](quartic-charge-fusion.md):

```
(X,Y,Z)=(x^2,y^2,xy).
```

With `d=q-1`, let `Z_a(t)=sum_i t^(s_(a,i))` for the local Smith
orders of the **full** cut matrix over the other variable's rational
function field; an infinite order contributes zero. Each cut satisfies

```math
Z_a(t)\ \ge\ P_d(t):=
1+2dt^6+\binom d2t^{12}+\binom d3t^{30}+\binom d4t^{48},
\qquad 0<t<1.                                           \tag{1}
```

These are partition **lower bounds**, not a full tensor Smith profile or
claims that its individual layers equal the displayed multiplicities.
For the deliberately unoptimized threshold `d>=2^90`,

```math
P_d(t)\ge(1+dt^{32/3})^3.                               \tag{2}
```

Under the stated common-scale, unweighted two-axis profile invoice,
(2) forces `rho>=3/4>1/sqrt(2)`. Thus this fixed colour-symmetric curve
cannot reach the target by growing its alphabet. The threshold does not
exclude `q=8`; no threshold refinement, full-profile or alphabet search
is justified by this result alone. Weighted gauges, unequal cut scales,
actual coefficient-rank improvements, different curves, and
colour-assignment weights that break the symmetry used below are outside
the conclusion. The complete factoring exponent remains `1/2`.

The [registered question](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19363241)
was separate from Sol worker 2's finite jet computation. The all-colour
[handoff](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19363379)
records the structural extension. This argument does not import the old
generic-`Z`-unit order-five theorem at a joint cold point.

## 1. A colour-alternation identity

Start with the one-hot Gaussian cube with pair kernel `K`. Add any product
of four-vertex creation gates whose coefficient on a fixed support is
invariant under permuting the assigned colours. Distinct-colour XOR-zero
quartics with a support-dependent scalar are examples. Work in
characteristic two, and fix a cut into four row and four column vertices.

Choose `r<=4` row vertices `I` and `r` column vertices `J`. Put a fixed set
`S` of `r` **distinct** nonzero colours on `I`, in one fixed order. For a
set `T` of `r` distinct nonzero colours, sum the tensor entries over all
`r!` assignments of `T` to `J`. All other vertices are vacuum. Then

```math
\sum_{\pi\in S_r} T_{S,\pi(T)}
 = [S=T]\det K[I,J].                                   \tag{3}
```

Here the sum is a constant linear combination of **different column
words**, not multiplication of one word by the even integer `r!`.

To prove (3), expand by the disjoint quartic supports and the remaining
Gaussian pair matching. Since colours on either face are distinct, every
remaining Gaussian pair must cross the cut. For equal row/column counts,
the quartics therefore consume equal total numbers of row and column
vertices. If there is a quartic, some quartic has at least two column
vertices. Choose such a support and two of its column positions by a
fixed spatial rule. Swapping their colours is a fixed-point-free
involution of the column assignments. It preserves that quartic's
coefficient and every other factor, so these terms cancel.

Without a quartic, matching is possible only when `S=T`. Summing the
colour-matched products over column assignments gives the permanent of
`K[I,J]`, equal to its determinant in characteristic two. This proves
(3), including terms with two disjoint quartics. The argument does not
use the special `q=8` [divided-square vanishing](charge-creation-square.md).

More generally, if the column word has more distinct colours than the
row word, the same expansion is zero after column alternation. Gaussian
pairs alone cannot match the counts. With quartics, strictly more column
than row vertices are consumed, so a quartic with at least two column
positions again supplies the cancelling involution. The restrictions
indexed by increasing occupation are consequently lower triangular.

## 2. Exact minors on the curve

Use the degree-four-per-occupied-vertex source of PR97, including its pair
kernel `K~=ell^2 K` and quartic coefficient `H`. Put
`L=x^2+xy+y^2`. Direct substitution gives

```math
\widetilde K_{03}=\widetilde K_{47}=0,\qquad
\widetilde K_{07}=\widetilde K_{43}=A=x^6y^6L^2,
```

and

```math
H=x^{20}y^{12}+x^{17}y^{15}+x^{15}y^{17}+x^{12}y^{20}.
```

Vertices are numbered by their binary `(x,y,z)` coordinates, with the
first coordinate in the low bit. The row face at the `x` cut is
`(0,2,4,6)` and its antipodally ordered column face is `(7,5,3,1)`.
All one-particle entries are homogeneous of degree sixteen. Setting
`y=1` preserves their `x`-adic Smith data: this is the change `x -> x/y`
over the coefficient field, up to scalar units, rather than a numerical
specialization of a general multivariate family.

After removing the local unit `L(x,1)^2` and scalar `x^6`, the cross matrix is

```math
N=\begin{pmatrix}
1&1+x^3&0&x^3+x^6\\
1+x^3&1&x^3+x^6&0\\
0&x^3+x^6&1&1+x^3\\
x^3+x^6&0&1+x^3&1
\end{pmatrix}.
```

Its constant matrix has rank two, and `N^2=x^12 I`. Hence the Smith orders
`e_1<=...<=e_4` of `N` have `e_1=e_2=0` and
`e_i+e_(5-i)=12`, by comparing `N^-1` with `x^-12 N`.
Thus the cross-kernel orders are exactly `(6,6,18,18)` and its minimum
minor valuations are `(6,12,30,48)`.

Witness row/column sets, in increasing minor size, are

```
(0 | 7),
(0,4 | 7,3),
(0,2,4 | 7,5,3),
(0,2,4,6 | 7,5,3,1).
```

Swapping the first two coordinate bits and `x,y` gives the other fixed
parent cut and the same data. Parent superletters always pair the two
`z` neighbours. No summand-dependent pairing or hidden-axis selection
is used.

## 3. Prefix minors give the partition lower bound

For each colour set of size `r`, take the witness in Section 2 and its
column alternation. Equation (3) gives a scalar-minor diagonal block,
with `binom(d,r)` copies. The vacuum is a unit. At size one there are
two copies per colour at order six, since the displayed `A I_2` is a
literal restriction. Restrictions of larger column occupation to smaller
row occupation vanish by Section 1. Their combined matrix is therefore
lower triangular with diagonal valuations

```
0^1, 6^(2d), 12^binom(d,2), 30^binom(d,3), 48^binom(d,4).
```

Every initial diagonal prefix has a minor whose valuation is its sum.
Thus each corresponding determinantal prefix of the full cut matrix is
at most that sum. The column-sum vectors have disjoint supports for
distinct colour sets and occupations, and can be extended to a constant
invertible column basis. Alternatively, Cauchy--Binet compares their
minors directly with the full determinantal ideals.

Here is an explicit prefix-to-partition argument. Write the selected
diagonal orders as `e_1<=...<=e_N` and the first `N` full Smith orders as
`s_1<=...<=s_N`. Every prefix sum of `s` is at most the corresponding
sum of `e`. For every `c>=0`,
`sum_i(c-s_i)_+=max_k(kc-sum_(i<=k)s_i)` is consequently at least
`sum_i(c-e_i)_+`. With `a=-log t>0`, integrate against `a^2 exp(-ac)`
and use `a^2 integral_0^infinity exp(-ac)(c-s)_+ dc=exp(-as)`.
This gives (1), with the remaining full Smith orders adding nonnegative
terms. It does not infer individual Smith orders from a single
determinant or specialization rank.

For size two the column operation is particularly transparent. Put
colours `a,b` at row vertices `0,4`; add the columns assigning them to
`7,3` and `3,7`. Every colour-symmetric quartic cancels. The resulting
matrix over all unordered colour pairs is exactly `A^2 I_binom(d,2)`.
These column additions preserve occupation height and are made within
the fixed physical labels.

## 4. An exact large-alphabet pressure invoice

For `d>=6`, write `B2=binom(d,2)` and `B4=binom(d,4)`. Then
`B2>=d^2/4` and `B4>=d^4/192`. Discard the nonnegative size-three term
in (1), and split the remaining nonconstant sum into

```
T1=2d t^6+(B2/3)t^12,
T2=(B2/3)t^12+(B4/2)t^48,
T3=(B2/3)t^12+(B4/2)t^48.
```

Weighted AM--GM with weights `(2/9,7/9)`, `(20/27,7/27)` and `(4/9,5/9)`
respectively gives

```math
T_1\ge \frac{d^{16/9}}{12}t^{32/3},\quad
T_2\ge \frac{d^{68/27}}{384}t^{64/3},\quad
T_3\ge \frac{d^{28/9}}{384}t^{32}.
```

At `d>=2^90`, these coefficients are at least `3d`, `3d^2`, `d^3`.
It suffices to check `d^(7/9)>=36`, `d^(14/27)>=1152`,
and `d^(1/9)>=384`; these follow from powers of two at the stated
threshold. Adding the constant proves (2), uniformly for `0<t<1`.
No decimal minimization or threshold search is involved.

For clarity, let `nu` be a common parent-superletter type, `w` its number
of occupied old letters (zero, one or two), `Ew` its mean and `Hnu>0` its
entropy. The ordinary balanced, unweighted two-axis invoice is

```math
\Phi_\nu(t)=
\frac{\tfrac12(\log Z_x(t)+\log Z_y(t))-16Ew\log t}{2Hnu}.
```

For an independent pair of a common old-letter type `mu`, this is
`[log Z(t)-32 E_mu log t]/[4H(mu)]`, as in the registered curve proposal.
The parent-letter generating function is `W(v)=(1+dv)^2`, so Gibbs'
entropy inequality gives `log W(v)>=Hnu+Ew log v`. Put `v=t^(32/3)`.
Equations (1)--(2) give

```math
\tfrac12(\log Z_x+\log Z_y)
\ge\tfrac32\log W(v)
\ge\tfrac32 Hnu+16Ew\log t.
```

Hence `Phi_nu(t)>=3/4`. The strict target comparison is `9/16>1/2`.
This concerns the common-scale profile envelope and the stated homogeneous
degree invoice. It is not an actual-rank or factoring lower bound.

## 5. Verification, credit and continuation

The [checker](../tools/check_colour_alternation.py) and
[source-pinned receipt](../evidence/colour-alternation-joint-cold-curve.json)
(`make verify-colour-alternation`) verify the exact curve scalars, the one-particle matrix
identity and all minor valuations, and column alternation against direct
four-/six-/eight-occupied tensor expansions. Its size-four off-diagonal
control has a nonzero two-quartic contribution before alternation. Thus
the test does not rely on the `q=8` divided-square shortcut. The ordinary
proof establishes the all-alphabet assertion; finite controls do not
replace it. The AM--GM exponents and coefficient thresholds are checked
with exact rational/integer arithmetic.

The [full default receipt](../evidence/colour-full-verification.json) and
[complete compressed stdout](../evidence/colour-full-verification.stdout.txt.xz)
record exit-zero `make verify` at clean frozen
`54acb0ed00ca48a4c100ec928a65eba9fd09e9df`, tree
`6a550332a23872ce00b703b2c8171b2aa5e30f64`, in 1685.02 seconds.
The fresh marker probe verified private `/tmp`, a read-only host root
and writable frozen checkout; network was unshared. Tracked state
stayed clean. Later main integrations and archival prose are outside
that source pin and retain separate checks; PR110/111 and the later
spatial-mode construction are not covered by this replay.

The original homogeneous Gaussian and worker-4 colour-class kernel retain
their attribution through the linked source notes. The quartic construction
belongs to the primary worker's PR97, the global divided-square work to
Sol worker 2's PR98, and the fixed bivariate curve and transfer question
to Sol worker 2. Sol worker 1 supplied this colour-alternation restriction,
curve minor check and pressure deduction with substantial OpenAI Codex
assistance. No novelty claim is made for alternation, Smith determinantal
prefixes, convex majorization, AM--GM or entropy duality.

The next assumption to change is the colour-permutation symmetry used
in (3), or another explicitly excluded part of the transfer. Merely
increasing the alphabet on this curve does not resolve its pressure
invoice. An orientation-sensitive colour interaction would require a new
construction and complete cost analysis; none is asserted here. No new
Lean or independent mathematical review is claimed, and the complete
factoring bound remains `1/2`.
