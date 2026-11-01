# A parity obstruction for homogeneous flat cube tensors

## Result and scope

Let an eight-leg binary tensor have polynomial coefficients of one common
homogeneous degree `D`, with the antipodal equality tensor as its value at
`(1,1,1)`. Suppose every odd-occupation coefficient vanishes. The tensor
need not be Gaussian, and its vacuum coefficient is unrestricted apart
from these conditions. The coefficient field is arbitrary.

For each of the three physical cuts, let `s_(a,i)` be the sixteen Smith
orders over the rational function field in the other two variables. Put

```math
Z_a(t)=\sum_{i=1}^{16}t^{s_{a,i}},\qquad
R(t)=\frac{\frac13\sum_a\log Z_a(t)-\frac D3\log t}{4\log2}.
\tag{1}
```

Then, for every `0<t<1`,

```math
R(t)>\frac{289}{500}>\frac1{\sqrt3}.
\tag{2}
```

Thus the common-degree cube criterion cannot cross its target within
the even-occupation family. This includes the degree-six flat tensors
obtained from polynomially cleared Gaussian Pfaffians. It also includes
arbitrary even coefficients and arbitrary allowed vacuum coefficients.
Odd-occupation coefficients, different degree rules, and other pressure
criteria are outside this statement. This is a criterion obstruction,
not a lower bound for integer factoring.

The proof uses shared parity minors and two exact entropy certificates.
The [standard-library checker](../tools/check_cube_flat_parity_floor.py)
replays their integer inequalities and the physical minor indexing.
The [evidence](../evidence/cube-flat-parity-floor.json) records the finite
checks. This complements the occupation-graded
[homogeneous profile budget](cube-homogeneous-profile-budget.md): here
every physical word has the same degree, so that earlier degree rule
does not apply.

## 1. Physical cuts and the common scalar factor

Write the vertices as `v in {0,1}^3`, and let `bar(v)=v xor 7`.
The coefficient indexed by a subset `U` is `S_U`. The endpoint condition is

```math
S_U(1,1,1)=[U=\bar U].
\tag{3}
```

For axis `a`, a row is a subset of the four vertices with `v_a=0`.
A column is a subset of the four vertices with `v_a=1`. The entry is
`S_(row union column)`. Odd coefficients vanish, so the cut splits into
two eight-by-eight blocks: even row and column parity, and odd row and
column parity. Equation (3) makes both blocks invertible at the endpoint.

For each axis, let `m_a` be the smallest valuation among all nonzero
coefficients. Every coefficient is divisible by
`x^m_x y^m_y z^m_z`. Divide out this common monomial and replace

```math
D\longmapsto D'=D-m_x-m_y-m_z,\qquad
s_{a,i}\longmapsto s'_{a,i}=s_{a,i}-m_a.
\tag{4}
```

The endpoint is unchanged, all coefficients remain homogeneous
polynomials, and expression (1) is unchanged. If `D'=0`, the resulting
tensor is exactly its endpoint tensor, and `R(t)=1`.
Otherwise use the divided tensor and write its degree as `D` again.
Its smallest entry valuation is zero at every axis. Consequently, at
each axis at least one parity block has first Smith order zero.

## 2. Two shared parity minors

For parity `p in {0,1}`, write the sorted orders of its block at axis `a`
as `lambda_(a,p,1),...,lambda_(a,p,8)`, and set

```math
q_{a,p}=\lambda_{a,p,1}+\lambda_{a,p,2},\qquad
d_{a,p}=\sum_{i=1}^{8}\lambda_{a,p,i}.
```

The block determinant is a nonzero homogeneous polynomial of degree
`8D`. Therefore

```math
d_{a,p}\le8D.
\tag{5}
```

Fix distinct axes `a,b`. Their common zero face contains two vertices
`v,w`. For even parity, choose the two row words
`empty` and `{v,w}`. For odd parity, choose `{v}` and `{w}`.
Use their antipodal images as the two column words.

These are the same four physical coefficients at both cuts. At either
cut they belong to its parity-`p` block. At the endpoint this two-by-two
minor is the identity, so its determinant `F_(a,b,p)` is nonzero and
homogeneous of degree `2D`. Its coordinate valuations satisfy

```math
q_{a,p}\le v_{x_a}(F_{a,b,p}),\qquad
q_{b,p}\le v_{x_b}(F_{a,b,p}).
```

The sum of two coordinate valuations cannot exceed the total degree of
a nonzero polynomial. Hence

```math
q_{a,p}+q_{b,p}\le2D.
\tag{6}
```

Adding the three axis-pair inequalities also gives

```math
\frac13\sum_a q_{a,p}\le D.
\tag{7}
```

No specialization of either other coordinate occurs in this argument.

## 3. An entropy bound for a block profile

For any probability vector `w` on the sixteen orders of one cut,
the elementary log-sum inequality gives

```math
\log Z_a(t)\ge H(w)+\log t\sum_iw_i s_{a,i},
\qquad H(w)=-\sum_iw_i\log w_i.
\tag{8}
```

We give weights, possibly depending on the axis, whose averaged
weighted order is at most `D/3`. Since `log t<0`, (8) then bounds (1)
below by the averaged entropy divided by `4 log 2`.

Choose one zero-minimum parity block at every axis. Up to permuting
axes and interchanging parity names, there are two cases.

### 3.1. The same block has minimum zero at all three axes

Use the following weights at every axis:

| Block | First order | Second order | Each of the last six orders |
|---|---:|---:|---:|
| Zero-minimum block | `86/150` | `4/150` | `2/150` |
| Other block | `18/150` | `18/150` | `2/150` |

They sum to one. Write `a=4/150`, `b=18/150`, and `c=2/150`.
The zero first order, (5), and (7) show that the averaged weighted
order is at most

```math
16Dc+(a-c)D+(b-c)D
 =D(a+b+14c)=D/3.
\tag{9}
```

The relevant entropy is `H_same` of the sixteen listed weights.

### 3.2. The zero-minimum parities have a two-plus-one split

Call the two repeated axes `A_1,A_2` and the remaining axis `B`.
Take parity zero to have minimum zero at the two `A` axes, and parity
one to have minimum zero at `B`. Use these weights:

| Axis class and parity | First order | Second order | Each of the last six orders |
|---|---:|---:|---:|
| `A`, parity zero | `285/500` | `23/500` | `6/500` |
| `A`, parity one | `57/500` | `57/500` | `7/500` |
| `B`, parity zero | `42/500` | `42/500` | `8/500` |
| `B`, parity one | `315/500` | `11/500` | `7/500` |

Each axis has total weight one. Put

```math
\begin{aligned}
a&=(q_{A_1,0}+q_{A_2,0})/2,&
e&=(q_{A_1,1}+q_{A_2,1})/2,\\
u&=q_{B,0},& d&=q_{B,1}.
\end{aligned}
```

Using the appropriate instances of (6) gives

```math
a+u\le2D,\qquad e\le D,\qquad e+d\le2D.
\tag{10}
```

After bounding the last six orders in each block with (5), the
averaged weighted order is at most

```math
\frac13\left(
  \frac{328}{500}D
 +\frac{34}{500}(a+u)
 +\frac{96}{500}e
 +\frac{4}{500}(e+d)\right).
\tag{11}
```

For example, the coefficient `34/500` agrees on the two sides because
`2(23-6)=42-8`. The remaining coefficients agree because
`2(57-7)-(11-7)=96`. All coefficients in (11) are nonnegative.
Equation (10) bounds (11) by

```math
\frac D3\frac{328+68+96+8}{500}=D/3.
\tag{12}
```

Let `H_A,H_B` denote the entropies of the two sixteen-entry weight
vectors. The averaged entropy is `(2H_A+H_B)/3`.

## 4. Exact entropy certificates

For a list of positive integer numerators `n_i` with sum `N`,

```math
H(n_i/N)=\log N-\frac1N\sum_i n_i\log n_i.
```

The checker verifies the following two strict integer comparisons:

```math
150^{750}
 >2^{1734}\left(86^{86}4^4 18^{36}2^{24}\right)^5,
\tag{13}
```

```math
500^{1500}
 >2^{3468}
  \left(285^{285}23^{23}57^{114}6^{36}7^{42}\right)^2
  \left(315^{315}42^{84}8^{48}11^{11}7^{42}\right).
\tag{14}
```

Taking logarithms proves, respectively,

```math
\frac{H_{\rm same}}{4\log2}>\frac{289}{500},\qquad
\frac{2H_A+H_B}{12\log2}>\frac{289}{500}.
\tag{15}
```

For orientation only, these two certified values are approximately
`0.5825778910` and `0.5788907290`. The theorem does not use floating
point arithmetic. Finally `3*289^2=250563>500^2`, so the rational
bound in (2) is strictly above the cube target.

## 5. Replay and next scope

Run:

```sh
python3 tools/check_cube_flat_parity_floor.py
```

The checker verifies all weight identities, the nonnegative linear
combination, both exact entropy comparisons, and every endpoint entry
in the six shared two-by-two parity minors. It also replays those
minors as explicit multivariate polynomials on twelve independently
generated even homogeneous tensors of degrees one through six.
These finite examples check the indexing and polynomial degree
argument; the proof above establishes the uniform statement.

The obstruction makes odd-occupation coefficients a necessary next
ingredient for the present flat-degree cube criterion. Merely changing
the Gaussian kernel, the even coefficients, or the vacuum polynomial
cannot cross its target under the stated hypotheses.
