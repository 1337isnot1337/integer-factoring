# A false positive from diagonal specialization of an encoded cube tensor

## Result and scope

A two-mode Gaussian kernel, restricted to the local alphabet consisting
of the vacuum and the occupied pair, gives a concrete non-Gaussian
eight-leg tensor. Its three physical cuts pass the cube target after
setting the other two variables equal. The resulting apparent bound is
strictly below `5/9 < 1/sqrt(3)` by an exact integer certificate.

The same tensor fails when the other variables are independent. An
exact finite-field specialization gives a profile dominated by the
existing Gaussian baseline, at every scale. Thus the tensor's actual
unweighted, common-type profile criterion stays above
`3/5 > 1/sqrt(3)`. This is a rejected construction, not an improved
factoring bound.

The example is useful because **every nonzero diagonal specialization**
has the false favorable profile, over arbitrarily large extension
fields. Increasing the field size does not repair a test that keeps the
two other coordinates equal. Unequal points in `F4` already reject this
example; `F256` also detects an additional valuation hidden at every
nonzero `F4` point.

The [standard-library checker](../tools/check_cube_pair_encoding_specialization.py)
and [exact evidence](../evidence/cube-pair-encoding-specialization.json)
retain the construction, direct sixteen-mode Pfaffians, three complete
Smith profiles, a small explicit minor, and the rational false-positive
certificate. The statement concerns this tensor and this profile
criterion. It does not exclude other pair encodings or optimized
physical gauges.

## 1. The encoded tensor

Work in characteristic two. Label the cube vertices by
`v in {0,1}^3`, with antipode `bar(v)=v xor 7`. Put two binary Gaussian
modes, numbered zero and one, at each vertex. Let `I` be the two-by-two
identity and let

```math
J=\begin{pmatrix}0&1\\1&0\end{pmatrix}.
```

For each nonempty set of coordinate axes, define a symmetric mode
matrix as follows. Write the three variables as `x,y,z`.

```math
\begin{aligned}
c_a&=x_a^6 I,\\
c_{ab}&=x_a^3x_b^3 I
 +(x_a^5x_b+x_ax_b^5)\begin{pmatrix}0&0\\0&1\end{pmatrix}
 +(x_a^4x_b^2+x_a^2x_b^4)\begin{pmatrix}1&0\\0&0\end{pmatrix},\\
c_{123}&=(xyz)^2 I+sJ,\\
s&=xyz(x^3+y^3+z^3+xyz).
\end{aligned}
\tag{1}
```

For distinct vertices, set

```math
K_{(v,i),(w,j)}
 =\left(\sum_{S\supseteq\operatorname{supp}(v\mathbin{\mathrm{xor}}w)}c_S\right)_{ij}.
\tag{2}
```

Set every same-vertex block to zero. This gives an alternating
sixteen-by-sixteen matrix. Each nonzero entry is homogeneous of degree
six. Every matrix `c_S` is `I` at `(1,1,1)`, so

```math
K_{(v,i),(w,j)}(1,1,1)=[w=\bar v][i=j].
\tag{3}
```

For a physical word `U subseteq {0,...,7}`, retain both modes at the
vertices in `U`, and define

```math
T_U=\operatorname{Pf}K_{U\times\{0,1\}},\qquad T_\varnothing=1.
\tag{4}
```

This is the local restriction `0 -> vacuum`, `1 -> occupied pair`.
Every nonzero `T_U` has degree exactly `6|U|`. Equation (3) gives the
required physical endpoint:

```math
T_U(1,1,1)=[U=\bar U].
\tag{5}
```

For an axis, arrange these coefficients into a sixteen-by-sixteen
matrix. Rows are subsets of the four vertices on its zero side;
columns are subsets on its one side. The entry at `(R,C)` is `T_(R union C)`.
The endpoint is an invertible antipodal permutation matrix. The
construction is symmetric under permutations of the three axes.

## 2. A smaller exact formula for the physical coefficients

Let `A` and `B` be the eight-by-eight kernels formed by modes zero and
one. The cross-mode entry between any two distinct vertices is `s`:
the off-diagonal mode part in (1) occurs only in `c_123`.

Temporarily add the same-vertex pair entry `s` at each vertex. The full
cross-mode matrix is now the rank-one matrix `s 1 1^T`. Denote the
resulting physical coefficients by `H_U`. Expanding a Pfaffian by its
cross-mode edges gives

```math
H_U=
\begin{cases}
 \operatorname{Pf}(A_U)\operatorname{Pf}(B_U),&|U|\text{ even},\\[2mm]
 s\left(\displaystyle\sum_{v\in U}\operatorname{Pf}(A_{U\setminus\{v\}})\right)
  \left(\displaystyle\sum_{w\in U}\operatorname{Pf}(B_{U\setminus\{w\}})\right),
 &|U|\text{ odd}.
\end{cases}
\tag{6}
```

Terms with two or more cross-mode edges vanish: their total coefficient
is a minor of a rank-one matrix. The remaining number of cross-mode
edges has the same parity as `|U|`, which proves (6).

Adding the local pair entries acts separately at each physical vertex:

```math
H_U=\sum_{V\subseteq U}s^{|U\setminus V|}T_V.
\tag{7}
```

The one-vertex transformation is
`[[1,0],[s,1]]`. It is polynomial and unimodular, and its square is the
identity in characteristic two. Consequently (7) preserves every
physical cut's Smith orders. The checker computes (6), applies (7)
again to recover `T`, and independently compares all 256 coefficients
with direct sixteen-mode Pfaffians at each of three points.

## 3. What the diagonal substitution reports

Setting the two other variables to one gives the exact physical Smith
orders

```math
P_{\mathrm{diag}}=
(0,7,7,14,17,17,24,24,24,24,31,31,34,41,41,48).
\tag{8}
```

These are the subset sums of `(7,7,17,17)`, so their counting function is

```math
Z_{\mathrm{diag}}(t)=(1+t^7)^2(1+t^{17})^2.
```

For physical frequency `0<p<1`, the homogeneous common-profile
expression for this degree-six physical alphabet is

```math
R(t,p)=\frac{\log Z(t)-16p\log t}{4H(p)},
\qquad H(p)=-p\log p-(1-p)\log(1-p).
\tag{9}
```

At the rational point `t=41/60`, `p=3/50`, it gives `R<5/9`.
No numerical logarithm is needed to certify this. Put

```math
A=(60^7+41^7)^2(60^{17}+41^{17})^2.
```

The checker compares the integers

```math
A^{450}3^{60}47^{940}
 <60^{21168}41^{432}50^{1000}.
\tag{10}
```

Taking logarithms of (10) is exactly `R(41/60,3/50)<5/9`.
Also `3*5^2<9^2`, so `5/9<1/sqrt(3)`.

This false positive persists at every diagonal pair `(r,r)` with
`r != 0`, in every coefficient field. Homogeneity gives

```math
T_U(x,r,r)=r^{6|U|}T_U(x/r,1,1).
\tag{11}
```

The powers of `r` split into row and column units, and the substitution
`x -> x/r` preserves the valuation. Hence (8) is unchanged.

At the diagonal point, every odd-occupation physical coefficient
vanishes. This is another specialization artifact: 120 of the 128 odd
physical words have nonzero polynomials at the `F256` point below.

## 4. Independent coordinates reject the construction

Use `F256=F2[u]/(u^8+u^4+u^3+u+1)`. Integers from zero through 255
encode polynomial-basis field elements. The checker verifies an inverse
for every nonzero element, so the field arithmetic is exact.

At `(y,z)=(2,3)`, the physical orders are

```math
P_{256}=(0,2,3,3,3,3,4,6,11,11,11,11,18,20,22,24).
\tag{12}
```

The same complete profile occurs at `(5,7)` and `(17,91)`. The two
nontrivial roots of `a^2+a+1` in this representation are 188 and 189.
For all six unequal pairs in `{1,188,189}^2`, the profile is instead

```math
P_4=(0,2,3,3,3,3,4,6,13,13,13,13,22,24,26,28).
\tag{13}
```

All three diagonal `F4` pairs retain (8). Thus unequal `F4` coordinates
already suffice to reject this example. The larger field removes
additional cancellations; it is not necessary for the rejection.

There is a particularly small witness to the failed generic valuation.
In the `x` cut, take rows `emptyset,{0}` and columns `emptyset,{7}`.
Both one-vertex coefficients vanish, so this minor is `T_{ {0,7} }`.
Directly from (1),

```math
T_{\{0,7\}}=(xyz)^4+s^2
 =x^2y^2z^2(x^6+y^6+z^6).
\tag{14}
```

At `(2,3)`, it is `31 x^2+20 x^8`. Thus the generic two-by-two
minor budget is at most two, whereas (8) would require seven.
At every nonzero `F4` pair, the coefficient of `x^2` in this particular
minor vanishes because `y^6=z^6=1`. This is the precise finite-field
blind spot; it does not mean every unequal `F4` profile is favorable.

For an all-scale rejection, let

```math
B=(0,4,8,12,16,20,20,24,24,28,28,32,36,40,44,48).
\tag{15}
```

This is the doubled physical profile of the existing
[homogeneous Gaussian seed](cube-homogeneous-gaussian.md). The checker
verifies, for every prefix length `1<=j<=16`,

```math
\sum_{i\le j}(P_{256})_i\le\sum_{i\le j}B_i,
\qquad
\sum_{i\le j}(P_4)_i\le\sum_{i\le j}B_i.
\tag{16}
```

The generic determinantal valuations are no larger than their
specialized values. Therefore its sorted orders satisfy the same
prefix bounds. For `0<t<1`, the function `a -> t^a` is decreasing and
convex. Weak majorization applied to (16) gives `Z_generic(t)>=Z_B(t)`.
The degree charge and binary entropy in (9) are unchanged.

The baseline pressure is the original Gaussian pressure after replacing
`t` by `t^2`. Its [exact all-scale floor](cube-gaussian-profile-budget.md)
is strictly greater than `3/5`. Consequently the actual unweighted,
common-type pressure of (4) is also strictly greater than `3/5`.
Since `3*3^2>5^2`, it cannot reach `1/sqrt(3)`.

## 5. Exact replay and use as a search control

Run:

```sh
python3 tools/check_cube_pair_encoding_specialization.py
```

The replay uses only the Python standard library. It checks 768 direct
sixteen-mode Pfaffians, every physical endpoint word, all nine nonzero
`F4` coordinate pairs, five further `F256` points, the explicit minor
(14), all prefix comparisons (16), and the integer inequality (10).
It also retains attaining row and column sets for every Smith prefix.

The Smith computation works in the truncated ring `F256[x]/(x^64)`.
At each step it selects an entry of minimum valuation, divides by its
unit part, and clears the remaining block. All sixteen pivots resolve
strictly below 64. Such unit operations determine the corresponding
orders of the original polynomial matrix exactly; an unresolved pivot
would be a failure, not a favorable infinite order.

For subsequent searches, a favorable specialization is only a screening
signal. Independent nonzero coordinate choices can reject it cheaply.
A surviving candidate still needs a generic multivariate certificate,
or an exact factorization that works over the appropriate rational
function field.
