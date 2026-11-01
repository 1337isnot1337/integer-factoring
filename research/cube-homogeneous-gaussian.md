# An explicit common three-cut Gaussian cube family

## Result and scope

There is an effective eight-leg crossing family satisfying the
[three-axis cube interface](cube-crossing-transfer.md) with

\[
\rho=\frac58.
\tag{1}
\]

The same polynomial coefficient tensor is factored at all three cuts.
The construction uses a homogeneous Gaussian kernel in three variables,
followed by a common constant-weight restriction on every physical leg.
Its complete three-cut profile has a rational pressure certificate.
The actual bound from the selected rational parameters is approximately
`rho <= 0.624832131518`, which is strictly smaller than `5/8`.

This gives an actual, unconditional instance of the cube interface. Its
factoring exponent is `log_3(8/5)`, approximately `0.427816`. It does
**not** improve the existing exponent one half. The local target for that
purpose remains `rho < 1/sqrt(3)`.

The construction improves the direct bound obtained by placing two
independent four-leg families on the four antipodal pairs. It is also
pointwise optimal for the common-scale, unweighted Smith-profile method
within the homogeneous Cayley class specified in section 7. An exact
positive-coefficient certificate proves a scoped floor `rho > 3/5` for
that method. Further progress therefore needs a different kernel class,
a stronger physical restriction, or additional information beyond those
profiles.

## 1. A single homogeneous kernel

Work over `F2`. Let `g_x,g_y,g_z` generate the group `(Z/2)^3`, and put
`eps_a=1+g_a`. Thus `eps_a^2=0`. Use three polynomial variables `x,y,z`.
In the group algebra define the seven nonconstant augmentation
coefficients by

\[
\begin{array}{lll}
c_x=x^6,&c_y=y^6,&c_z=z^6,\\
c_{xy}=x^3y^3,&c_{xz}=x^3z^3,&c_{yz}=y^3z^3,\\
&&c_{xyz}=x^2y^2z^2.
\end{array}
\tag{2}
\]

Set `c_empty` equal to the sum of these seven coefficients, and write

\[
f=\sum_{S\subseteq\{x,y,z\}}c_S\prod_{a\in S}\varepsilon_a
 =\sum_{s\in\{0,1\}^3}k_sg^s.
\tag{3}
\]

The coefficient `k_0` is zero. On eight modes indexed by cube vertices,
let `K_uv=k_(u xor v)`. This is an alternating matrix, with homogeneous
entries of degree six. At `(x,y,z)=(1,1,1)`, all eight coefficients `c_S`
are one, so

\[
f(1,1,1)=\prod_a(1+\varepsilon_a)=g_xg_yg_z.
\tag{4}
\]

Consequently `K(1,1,1)` is exactly the antipodal matching matrix.

For a physical binary word `b=(b_v)_v`, let `S(b)={v:b_v=1}` and define

\[
G_b(x,y,z)=\mathrm{Pf}K[S(b),S(b)].
\tag{5}
\]

The empty Pfaffian is one and odd Pfaffians are zero. Equation (4) gives
the complete endpoint

\[
G_b(1,1,1)=\prod_{\{v,\bar v\}}[b_v=b_{\bar v}].
\tag{6}
\]

Every nonzero `G_b` is homogeneous of degree `3|S(b)|`. This degree
identity is the conservation law used below. All seven nonzero group
steps have nonzero polynomial coefficients, so the one-particle kernel
has couplings across every nontrivial physical partition. In particular,
this seed is not a product of two independent Gaussian seeds on disjoint
four-leg sets.

## 2. The three exact crossblock profiles

Across the `x` plane cut, the four-by-four one-particle crossblock is
multiplication by

\[
a+b\varepsilon_y+c\varepsilon_z+d\varepsilon_y\varepsilon_z,
\quad
(a,b,c,d)=(x^6,x^3y^3,x^3z^3,x^2y^2z^2).
\tag{7}
\]

In the augmentation basis `1,eps_y,eps_z,eps_y eps_z`, its matrix is

\[
B_x=\begin{pmatrix}
a&0&0&0\\
b&a&0&0\\
c&0&a&0\\
d&c&b&a
\end{pmatrix}.
\tag{8}
\]

The group-basis and augmentation-basis matrices differ by the Boolean
zeta matrix, which is its own inverse over `F2`.

There is an exact diagonalization of (8) over
`F2[y,y^-1,z,z^-1][x]` by invertible row and column operations:

\[
\mathrm{diag}
\left(x^2y^2z^2,\ x^4y^4z^{-2},\
      x^8y^{-4}z^2,\ x^{10}y^{-2}z^{-2}\right).
\tag{9}
\]

All operation multipliers have nonnegative powers of `x`; their inverse
operations have the same property. Thus (9) is a finite Laurent-polynomial
factorization, not merely a generic-rank calculation over a rational
function field. The exact operations are retained in the certificate.

The four `x`-adic Smith orders are therefore

\[
(2,4,8,10).
\tag{10}
\]

The `y` and `z` cuts are permutations of this same construction and have
the same orders at their respective variables. They are cuts of the
single tensor (5), with no change of coefficient pieces.

There is also a short independent verification of the orders. The
smallest one-minor valuation is two and the smallest two-minor valuation
is six. Moreover `B_x^2=x^12 I`, so opposite Smith orders sum to twelve.
These facts force (10). The checker evaluates all 69 nonempty square
minors separately at each cut.

## 3. Passing from four one-particle orders to sixteen physical orders

At a cut, write the alternating kernel as

\[
K=\begin{pmatrix}A&B\\B^{\mathsf T}&C\end{pmatrix}.
\]

The matrix of the physical Gaussian tensor has the factorization

\[
M=U\,\bigoplus_{r=0}^4\bigwedge^r B\,V.
\tag{11}
\]

Here rows and columns are indexed by subsets of four modes, and

\[
U_{R,R'}=[R'\subseteq R]\mathrm{Pf}A[R\setminus R'],
\quad
V_{S',S}=[S'\subseteq S]\mathrm{Pf}C[S\setminus S'].
\tag{12}
\]

To prove (11), expand each Pfaffian by the vertices paired across the
cut. The within-side pairings give (12), and the crossing pairings give
the corresponding minor of `B`. Signs disappear in characteristic two.
Both side factors are triangular with diagonal one.

Taking exterior powers of the exact operations in (9), and then using
(11), gives a finite Laurent-polynomial factorization of the complete
sixteen-by-sixteen physical cut. Its exponents at the cut variable are
all subset sums of (10). Hence their counting polynomial is

\[
Z(t)=(1+t^2)(1+t^4)(1+t^8)(1+t^{10}).
\tag{13}
\]

For the `x` cut all factor entries are polynomials in `x` and Laurent
polynomials in `y,z`, with exponents bounded by an absolute constant.
The analogous assertion holds at the other two cuts. The checker
reconstructs all 256 entries of (11) at every cut, using the full
three-variable polynomials.

## 4. One common coefficient decomposition at every word length

Choose `m` divisible by five and put `w=m/5`. The physical alphabet is
any common subset `Q` of the weight-`w` binary words of length `m`.
For eight letters `X_v in Q`, form

\[
F((X_v)_v;x,y,z)=\prod_{j=1}^m G_{(X_{v,j})_v}(x,y,z).
\tag{14}
\]

There are exactly `8w` occupied physical bits across these eight words.
The conservation law after (6) implies that every entry in (14) is
homogeneous of degree `24w`. Write

\[
F=\sum_{k_x+k_y+k_z=24w}x^{k_x}y^{k_y}z^{k_z}T_k.
\tag{15}
\]

The endpoint of (14) is the antipodal equality tensor on `Q`, so summing
these same coefficient tensors gives that endpoint. There are at most
`binom(24w+2,2)` pieces.

We next bound the rank of one coefficient at all three cuts. At the `x`
cut, tensor the sixteen-state factorization from section 3 over the `m`
sites. A channel sequence has order equal to the sum of its `m` orders.
The left and right factors still have nonnegative `x` exponents and
other exponents bounded in absolute value by `O(m)`.

If a sequence has order greater than `k_x`, it cannot contribute to
`T_k`. For any other sequence, extracting one three-variable coefficient
from the product of its left and right factors requires only `O(m^3)`
separable terms. Restricting the rows and columns to the physical alphabet
cannot increase rank. Therefore

\[
\mathrm{rank}_x(T_k)
 \le \mathrm{poly}(m)\,
 \#\{\text{channel sequences with order at most }k_x\}.
\tag{16}
\]

For any real `0<t<1`, the right-hand count is at most
`t^(-k_x) Z(t)^m`. Apply the same argument to the other two cuts of the
same `T_k`. Using the exact conservation relation in (15) gives

\[
\prod_{a=x,y,z}\mathrm{rank}_a(T_k)
 \le \mathrm{poly}(m)\,t^{-24w}Z(t)^{3m}.
\tag{17}
\]

The Laurent coefficients cause no infinite coefficient sum: all changes
of basis used here are finite, and their tensor powers have only `O(m)`
exponent range in each variable. This is why the coefficient extraction
in (16) has a polynomial cost.

## 5. A rational pressure certificate and full effectiveness

Set `r=t^2` and define

\[
Z_0(r)=(1+r)(1+r^2)(1+r^4)(1+r^5).
\]

For `p=1/5`, let

\[
\rho_*=
\frac{\log_2 Z_0(3/7)-(4/5)\log_2(3/7)}
     {4H_2(1/5)}.
\tag{18}
\]

Since `H_2(1/5)=log_2(5)-8/5`, the strict inequality `rho_*<5/8` is
equivalent to the following rational inequality at `r=3/7`:

\[
Z_0(r)^{10}2^{40}<5^{25}r^8.
\tag{19}
\]

The checker verifies (19) by exact integer arithmetic. The ratio of its
left and right sides is approximately `0.996645566896`. No numerical
minimization is needed for the certificate.

For a requested wire width `lambda`, choose the smallest positive
multiple of five `m` such that `binom(m,m/5)>=2^lambda`, and take the
first `2^lambda` weight-`m/5` words as the common alphabet. The elementary
type-count bounds give

\[
mH_2(1/5)=\lambda+O(\log(\lambda+2)),\qquad m=O(\lambda).
\tag{20}
\]

For example, the lower type-count bound follows because the binomial
probability at its mode is at least `1/(m+1)`. The upper bound is the
corresponding probability bound by one. Minimality of `m`, with spacing
five, then gives (20).

Combining (17)--(20), and rounding the three logarithmic ranks upward,
gives channel widths

\[
r_k^x+r_k^y+r_k^z
 \le12\rho_*\lambda+O(\log(\lambda+2))
 \le12(5/8)\lambda+O(\log(\lambda+2)).
\tag{21}
\]

Every cut matrix has order `2^(4 lambda)`, so each separate width is at
most `4 lambda`.

All pieces and factors are uniformly effective within `2^O(lambda)`
bit operations. One may construct the full `2^(8 lambda)` physical
polynomial entries, multiply their `m` fixed seed factors, extract all
`O(lambda^2)` coefficients, and use ordinary binary elimination at each
of the three cuts. This charges all entries and all coefficients. It
also supplies exact factor-entry algorithms without a rank oracle or
an implicit exponentially large field. The kernel and coefficient field
are fixed independently of the factoring input.

Thus every hypothesis of the cube transfer theorem is met.

## 6. Comparison with direct products of four-leg pieces

Pair the four antipodal wires into two four-leg crossings. For either
crossing, two cube axes induce the same physical bipartition and the
third induces the other crossing bipartition. Suppose a four-leg family
has individual widths at most `2 lambda` and sum at most
`4 rho_2 lambda+O(log lambda)`.

For one such block, the three-axis width sum is of the form `2r+s` or
`r+2s`, and is bounded by
`2 lambda+4 rho_2 lambda+O(log lambda)`. The two blocks therefore give
the direct cube-interface bound

\[
\rho_3\le\frac{1+2\rho_2}{3}.
\tag{22}
\]

At the established two-cut threshold `rho_2=1/sqrt(2)`, (22) is about
`0.804737854125`. The construction above certifies `rho_3=5/8`.
Equation (22) is a bound for this direct composition, not a lower bound
for all possible products, regroupings, or joint coefficient extractions.

## 7. A complete profile boundary for homogeneous Cayley kernels

This section concerns the following specific method. Use a homogeneous
Cayley kernel of degree `D` with augmentation coefficients
`c_x=x^D,c_y=y^D,c_z=z^D`; let the three pair coefficients and the
triple coefficient be arbitrary homogeneous degree-`D` polynomials
that evaluate to one at `(1,1,1)`. Use a common physical constant-weight
type, and bound the three coefficient ranks through their unweighted
Smith counts at a common scale as in (16)--(17).

At one axis put

\[
e=\min(D,v(c_{ab}),v(c_{ac})),\qquad
f=\min(D,v(c_{abc})).
\]

For the four-by-four multiplication matrix (8), the smallest one-minor
valuation is `a=min(e,f)`. The smallest two-minor valuation is
`min(2e,D+f)`: squares of the two linear coefficients give the `2e`
term, and the remaining candidate is `c_a c_abc+c_ab c_ac`.
If its two leading terms could cancel, the square minor already gives
the same or a smaller valuation. The matrix squares to `x^(2D) I`, so
the complete ordered profile is

\[
(a,b,2D-b,2D-a),\qquad
b=\min(2e,D+f)-a.
\tag{23}
\]

Across the three axes the homogeneous degree budgets imply

\[
\sum_a a_a\le D,\qquad
\sum_a(a_a+b_a)\le3D.
\tag{24}
\]

For the first inequality, each smallest order is at most the valuation
of the triple coefficient, whose three valuations sum to at most `D`.
For the second, each pair coefficient allocates at most `D` to its two
axis valuations. Summing the three pair inequalities gives
`2 sum e_a<=3D`, and `a_a+b_a<=2e_a`.

Fix `0<t<1`. On `0<=s<=D`, the function

\[
\Phi(s)=\log(1+t^s)+\log(1+t^{2D-s})
\]

is decreasing and convex. Jensen's inequality and (24) give

\[
\sum_a\bigl(\Phi(a_a)+\Phi(b_a)\bigr)
 \ge3\bigl(\Phi(D/3)+\Phi(2D/3)\bigr).
\tag{25}
\]

Indeed, let `A=sum a_a`. Jensen first replaces each triple by its mean.
Monotonicity permits increasing the other total to `3D-A`. The function
`Phi(u)+Phi(D-u)` decreases for `0<=u<=D/2`, and `A/3<=D/3`.
This proves (25).

The left side of (25) is the logarithm of the product of all three
one-particle profile counting functions. The symmetric degree-six
construction attains equality at every scale. Thus asymmetric monomial
exponents, or polynomial pair/triple coefficients, cannot improve that
common-scale product within the stated class.

### An exact all-scale floor for this method

There is also an exact bound separating this method from the half-score
target. For every `0<r<=1`,

\[
Z_0(r)^5>(1+r^{5/3})^{12}.
\tag{26}
\]

Here is a finite integer certificate. Put `r=s^3`, subtract the right
side from the left side, and divide by `s^3`. The result `R(s)` has
degree 177. Every coefficient of

\[
(1+u)^{177}R\!\left(\frac{u}{1+u}\right)
\tag{27}
\]

is a positive integer; the smallest is five. The checker retains all
178 coefficients. Equation (27) proves positivity for `0<s<1`; direct
evaluation gives positivity at one.

The elementary binary entropy identity

\[
\max_{0\le p\le1}\{H_2(p)+p\log_2 q\}=\log_2(1+q)
\]

and (26), with `q=r^(5/3)`, imply

\[
\frac{\log_2 Z_0(r)-4p\log_2 r}{4H_2(p)}>\frac35
\quad(0<p<1).
\tag{28}
\]

Together with (25), this proves that the specified common-scale
Cayley-profile method cannot reach `rho<=3/5`, and hence cannot reach
`rho<1/sqrt(3)`. This is a boundary for that method, not for general
eight-leg families or other ways to estimate their coefficient ranks.

## 8. Retained checks and reproducibility

The [exact constructor](../tools/cube_homogeneous_gaussian.py),
[checker](../tools/check_cube_homogeneous_gaussian.py), and
[evidence](../evidence/cube-homogeneous-gaussian.json) retain:

- All 64 one-particle and 256 Gaussian endpoint entries, plus the
  homogeneous degree of every one of the 1,178 Gaussian monomials.
- Exact finite Smith transformations and all 207 crossblock minors.
- All 768 entries of the three complete Gaussian cut factorizations.
- 2,024 independently evaluated small local Smith profiles, through
  degree eight, using 139,656 literal minors.
- 545,635 complete common monomial lists through degree twelve, checking
  all three profiles and both coupled degree budgets in (24).
- One explicit five-site, weight-one binary subalphabet. Its 22 common
  coefficient pieces include exact rank factorizations at all three
  cuts, checking 16,896 factorization entries.
- The exact rational upper certificate (19) and every positive integer
  coefficient in the scoped lower certificate (27).

Run `python3 tools/check_cube_homogeneous_gaussian.py` to reproduce the
certificate. These checks verify the finite algebra supporting the proof;
the arbitrary-word-length result follows from sections 3--5.
