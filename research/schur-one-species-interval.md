# A one-species interval crossing with explicit odd Smith spectra

## Result and scope

The sparse Schur kernel has a particularly simple boundary when its second
species is empty. On an interval of modes, explicit physical gauges give
two copies of the positive odd orders in one cut and all signed odd orders
in the other cut. Their pressure is a finite product, so the remaining
type optimization is one dimensional.

The exact spectra below hold when the number of modes is either a power
of two or one less than a power of two. They give, among other rational
certificates,

| Modes | Certified rate parameter | Status |
|---:|---:|---|
| 7 | \(\rho<23/32=0.71875\) | Exact integer inequality |
| 8 | \(\rho<43/60=0.71666\ldots\) | Exact integer inequality |
| 63 | \(\rho<17/24=0.70833\ldots\) | Optional larger exact integer inequality |

These rates remain above \(1/\sqrt2\), the boundary needed to improve the
existing factoring exponent. An analytic argument proves that this
particular pressure formula stays strictly above that boundary at every
evaluation point, for both proved interval families. The family approaches
the boundary from above as the interval grows.

This is a statement about the displayed one-species restriction and
gauges. It does not rule out other gauges, a second occupied species,
different mode sets, or coefficient decompositions beyond this pressure
certificate. The derivation and checks used substantial OpenAI assistant
support.

## 1. The one-species boundary is a complete typed crossing

Use the polynomial kernel and boundary conventions of
[fixed-mode crossings](fixed-mode-crossings.md). Let \(M\) be a power of
two, let \(S=\{0,\ldots,n-1\}\subseteq\{0,\ldots,M-1\}\), and put

~~~math
z=w^2,\qquad q=1+z,\qquad c=1+w,
\qquad
P_{ab}=\mathbf1_{b\subseteq a},\quad
Q_{aa}=q^a,\quad
H_{ab}=z^{M-1-(a\mathbin{|}b)}.                         \tag{1}
~~~

The relation \(b\subseteq a\) denotes inclusion of binary digits.
Set the second species \(B\) to the empty set at every physical label.
The remaining occupied fermions are the minus species at boundaries
\(B,X\) and the plus species at boundaries \(V,Y\). Thus the tensor
entries are minors of

~~~math
\mathcal D(w)=
\begin{pmatrix}
cQP^{\mathsf T}&z^M P\\
I&cH
\end{pmatrix},                                        \tag{2}
~~~

whose rows are \((V,Y)\) and columns are \((B,X)\).

In a product of \(m\) copies, fix an occupation count \(k_a\) separately
for each remaining mode \(a\). Use this same set of words on all four
physical boundaries. The endpoint at \(w=1\) is exactly the crossing on

~~~math
D_m=\prod_{a=0}^{n-1}\binom{m}{k_a}                     \tag{3}
~~~

labels. This includes the case in which every second-species occupation
count is zero in the general fixed-mode theorem. In particular, no
limiting argument with a small positive second-species probability is
needed to define (2).

At the endpoint, the potentially nonidentity factor is a minor of
\(P\). Every nonzero term in such a minor matches a binary subset to a
superset. Fixed mode counts equalize the total index sums across the
copies, so all these inclusions must be equalities. The endpoint is
therefore diagonal with value one.

## 2. Physical gauges

List the remaining one-particle heights in boundary order
\((B_-,X_-,V_+,Y_+)\). At zero use

~~~math
h_0(a)=(0,\,-2M+1+2a,\,0,\,1+2a),                      \tag{4}
~~~

and at infinity use

~~~math
h_\infty(a)=(2a-n+1,\,2M-n,\,2a-n+1,\,-n).             \tag{5}
~~~

The valuation variable is \(w\) at zero and \(\xi=w^{-1}\) at infinity.
The row cut at zero has cross blocks \(z^M P\) and \(I\). The internal
blocks are \(cQP^{\mathsf T}\) and \(cH\).

The first internal block is regular because all its nonzero entries have
valuation zero and its two height contributions are zero. The gauged
\((a,b)\) entry of the second has valuation

~~~math
(1+2a)+(-2M+1+2b)+2(M-1-(a\mathbin{|}b))
=2(a\mathbin{\&}b)\ge0.                              \tag{6}
~~~

Consequently the internal Gaussian transformations are regular units.
The two cross blocks have the same order list

~~~math
1,3,5,\ldots,2n-1.                                    \tag{7}
~~~

For \(z^M P\), this follows by factoring out the diagonal column weights
and using the constant unit matrix \(P\). The other block is diagonal.

At infinity the column cut separates all remaining plus fermions from
all remaining minus fermions. There are no internal blocks to regularize:
the full flattening is the exterior compound of (2).

The combined restored height of the physical mode \(a\) is

~~~math
G_a=\sum_{e\in\{0,\infty\}}
 \bigl(h_{e,B}(a)+h_{e,X}(a)+h_{e,V}(a)+h_{e,Y}(a)\bigr)
=8a-4(n-1).                                           \tag{8}
~~~

## 3. Exact infinity spectrum

**Theorem.** If \(n=M\) or \(n=M-1\), the gauged infinity matrix (2)
has order list

~~~math
-(2n-1),\,-(2n-3),\ldots,-1,1,\ldots,2n-3,2n-1.       \tag{9}
~~~

Put \(t=\xi^2\), and in this section redefine

~~~math
D=\operatorname{diag}(t^a),\qquad
Q=\operatorname{diag}((1+t)^a),\qquad
O_{ab}=t^{a\mathbin{|}b}.
~~~

After removing the common factor \(\xi^{1-2n}\), multiply the first
block of rows by \((1+\xi)^{-1}\) and the second block of columns by
\(1+\xi\). These are regular units. The resulting matrix is

~~~math
F_n(t)=
\begin{pmatrix}
QP^{\mathsf T}D&DP\\
D&(1+t)O
\end{pmatrix}.                                        \tag{10}
~~~

It suffices to prove that its Smith orders in \(t\) are
\(0,1,\ldots,2n-1\).

Subtract \(QP^{\mathsf T}\) times the second block of rows from the
first, and interchange the two row blocks. The matrix becomes

~~~math
\begin{pmatrix}D&(1+t)O\\0&S\end{pmatrix},\qquad
S=DP+(1+t)QP^{\mathsf T}O.                              \tag{11}
~~~

Every entry in row \(a\) of \(O\) is divisible by \(t^a\).
Thus a regular column operation clears the upper-right block, giving
\(\operatorname{diag}(D,S)\).

### Full interval: \(n=M\)

A binary-digit expansion gives

~~~math
QP^{\mathsf T}O=(1+t)^{M-1}DP.                          \tag{12}
~~~

For a single bit, the sum vanishes unless the column bit is included
in the row bit. Otherwise its contribution is the corresponding factor
of \((1+t)^{M-1}t^a\). Since \(M\) is a power of two,
\(1+(1+t)^M=t^M\). Hence

~~~math
S=t^M DP.                                             \tag{13}
~~~

The orders of \(D\) are \(0,\ldots,n-1\), and those of \(S\) are
\(n,\ldots,2n-1\).

### One missing top mode: \(n=M-1\)

Remove the all-ones index \(M-1=n\) from the full identity (12).
Its omitted contribution is \(t^nQ\mathbf1\mathbf1^{\mathsf T}\), so

~~~math
S=t^n\bigl(tDP+(1+t)Q\mathbf1\mathbf1^{\mathsf T}\bigr).\tag{14}
~~~

The truncated inclusion matrix satisfies \(P^2=I\) and
\(\mathbf1^{\mathsf T}P=\mathbf1^{\mathsf T}\). Multiplying the
parenthesis on the right by \(P\) gives

~~~math
L=\operatorname{diag}(t,t^2,\ldots,t^n)
  +u\mathbf1^{\mathsf T},\qquad u_a=(1+t)^{a+1}.         \tag{15}
~~~

Every \(u_a\) is a unit. A \(k\)-minor of a diagonal matrix plus a
rank-one matrix contains at most one rank-one contribution. Its
valuation is therefore at least \(1+\cdots+(k-1)\). For \(k<n\),
the rows \(\{0,\ldots,k-2,n-1\}\) and columns
\(\{0,\ldots,k-1\}\) attain this bound. Their determinant is
\(u_{n-1}t^{1+\cdots+(k-1)}\).

For the full determinant, the unique term of least valuation omits
the diagonal factor \(t^n\). Its valuation is \(1+\cdots+(n-1)\).
The Smith orders of \(L\) are therefore \(0,\ldots,n-1\).
Equation (14) proves the required orders for \(S\), and hence for (10).
Replacing \(t\) by \(\xi^2\) and restoring \(\xi^{1-2n}\) gives (9).

## 4. Finite pressure and exact type optimization

For \(s>1\), the product of the two exterior profile sums is

~~~math
\mathcal B_n(s)
=s^{n^2}\prod_{j=1}^{n}(1+s^{-(2j-1)})^4.              \tag{16}
~~~

For a candidate rate \(\rho\), put \(R=4\rho\). The exact independent
type criterion from [fixed-mode crossings](fixed-mode-crossings.md) is

~~~math
\mathcal B_n(s)<\mathcal Z_{n,R}(s)^R,\qquad
\mathcal Z_{n,R}(s)
=\prod_{a=0}^{n-1}(1+s^{-G_a/R}).                       \tag{17}
~~~

The maximizing occupation probabilities are

~~~math
p_a=\frac{s^{-G_a/R}}{1+s^{-G_a/R}}.                    \tag{18}
~~~

They satisfy \(p_a+p_{n-1-a}=1\). No equal-count constraint is imposed
between the absent species and the remaining species.

When \(R=A/B\) and \(s=r^A\) with rational \(r>1\), (17) is an
integer comparison after taking the \(B\)-th power and clearing
denominators. The checker records exact witnesses for

~~~math
(n,R,r)=(7,23/8,50/49),\qquad(8,43/15,100/99),          \tag{19}
~~~

and optionally \((63,17/6,177/176)\). The probabilities (18) are then
rational as well, so common mode counts can realize them exactly along
a subsequence of copy lengths. The polynomial type losses do not alter
the strict exponential inequality.

For orientation, floating scalar optimization of the exact product
formula gives the following values. These decimals are exploratory;
the preceding rational witnesses are the certified finite statements.

| Modes | Approximate optimized \(\rho\) |
|---:|---:|
| 3 | 0.7412212202 |
| 7 | 0.7185927209 |
| 8 | 0.7166536597 |
| 15 | 0.7109577836 |
| 16 | 0.7106044045 |
| 31 | 0.7083773391 |
| 32 | 0.7083157559 |
| 63 | 0.7075175531 |
| 64 | 0.7075072113 |

## 5. The displayed formulas stay above half

First consider the displayed gauges with odd \(n=M-1\).
Write \(n=2k+1\), \(t=\log s\), and
\(f(x)=\log(1+e^{-x})\). After cancelling the filled negative
energies, the logarithmic difference in (17) is

~~~math
\begin{aligned}
\mathcal D_n(t,R)
&=\log\mathcal B_n(e^t)-R\log\mathcal Z_{n,R}(e^t)\\
&=t+4\sum_{j=1}^{n}f((2j-1)t)
 -R\log2-2R\sum_{j=1}^{k}f(8jt/R).                    \tag{20}
\end{aligned}
~~~

Let \(F(t)=\sum_{j\ge1}f(jt)\). The product definition and modular
transformation of Dedekind's eta function imply

~~~math
F(t)+F(2\pi^2/t)
=\frac{\pi^2}{12t}+\frac{t}{24}-\frac{\log2}{2}.        \tag{21}
~~~

The needed eta identities are stated in the NIST Digital Library of
Mathematical Functions, [27.14.12](https://dlmf.nist.gov/27.14.E12) and
[23.18.5](https://dlmf.nist.gov/23.18.E5). For completeness, set
\(P(t)=\sum_{j\ge1}\log(1-e^{-jt})\). Its modular transformation is
\(P(t)=-\pi^2/(6t)+\tfrac12\log(2\pi/t)+t/24+P(4\pi^2/t)\).
Subtract this identity from its version at \(2t\), and use
\(F(t)=P(2t)-P(t)\), to obtain (21).

Replacing both sums in (20) by infinite sums and applying (21) gives

~~~math
\begin{aligned}
\mathcal D_\infty(t,R)
={}&\frac{(8-R^2)\pi^2}{48t}+\frac{t}{6}\\
&+4\bigl(F(\pi^2/t)-F(2\pi^2/t)\bigr)
 +2R F(\pi^2R/(4t)).                                  \tag{22}
\end{aligned}
~~~

At \(R=2\sqrt2\), every term remaining in (22) is positive.
The finite truncation only increases the difference when \(n\ge3\).
Indeed, \(f\) is decreasing and convex. Its missing numerator tail
satisfies the midpoint bound

~~~math
\sum_{j=n+1}^{\infty}f((2j-1)t)
\le\frac1{2t}\int_{2nt}^{\infty}f(x)\,dx.              \tag{23}
~~~

The missing denominator tail, multiplied by \(R/2=\sqrt2\), satisfies
the decreasing-rectangle bound

~~~math
\frac R2\sum_{j=k+1}^{\infty}f(8jt/R)
\ge\frac1{2t}\int_{\sqrt2(n+1)t}^{\infty}f(x)\,dx.     \tag{24}
~~~

For \(n\ge3\), \(\sqrt2(n+1)\le2n\). Therefore (24) dominates
(23), and \(\mathcal D_n(t,2\sqrt2)\ge
\mathcal D_\infty(t,2\sqrt2)>0\) for every \(t>0\).
For \(n=1\), minimizing (16) gives \(256/27\), which also misses the
half threshold.

Conversely, for any fixed \(R>2\sqrt2\), the first term in (22) is
negative and dominates as \(t\) decreases. After choosing such a \(t\),
large enough odd intervals approximate the infinite sums and satisfy
(17). This proves approach to the half-rate boundary from above within
this explicit family.

For even \(n=M\), the cancelled expression is instead

~~~math
\mathcal D_n(t,R)
=4\sum_{j=1}^{n}f((2j-1)t)
 -2R\sum_{j=1}^{n/2}f(4(2j-1)t/R).                    \tag{25}
~~~

Its infinite-sum expansion has the same leading terms
\((8-R^2)\pi^2/(48t)+t/6\). The transformed remainder has a different
sign pattern, so we give a separate argument for all \(t>0\).

Put \(A(t)=\sum_{j\ge1}f((2j-1)t)\). Expansion of the logarithm gives

~~~math
A(t)=\frac12\sum_{r\ge1}\frac{(-1)^{r+1}}{r\sinh(rt)}.
~~~

The function \(tA(t)\) is strictly decreasing. To prove this, define

~~~math
h(x)=\frac{x\cosh x-\sinh x}{\sinh^2x}>0.
~~~

The sign of the derivative of \(h(x)/x\) is the opposite of the sign of

~~~math
F(x)=(x^2-1)\sinh^2x-x\sinh x\cosh x+2x^2.
~~~

For \(r\ge2\), the coefficient of \(x^{2r}\) in \(2F(x)\) is

~~~math
\frac{2^{2r-2}(4r^2-6r-4)}{(2r)!}.
~~~

It is zero at \(r=2\) and positive for every \(r\ge3\); all lower
coefficients vanish. Thus \(h(x)/x\) strictly decreases on the positive
axis. It follows that \(h(rt)/r\) is a strictly decreasing positive
sequence tending to zero. Pairing the terms of the convergent alternating
series now gives

~~~math
\frac{d}{dt}\bigl(tA(t)\bigr)
=-\frac12\sum_{r\ge1}(-1)^{r+1}\frac{h(rt)}r<0.         \tag{26}
~~~

At \(R=2\sqrt2\), the infinite version of (25) is
\(4(A(t)-\sqrt2A(\sqrt2t))>0\) by (26).
The same midpoint and decreasing-rectangle tail bounds as in (23)--(24)
apply when \(n\ge4\), so the finite expression is positive as well.

It remains to check the even size \(n=2\). At \(R=3\), put \(s=u^3\).
After clearing positive factors, the difference is the polynomial

~~~math
(u^3+1)^4(u^9+1)^4-u^{24}(u^4+1)^6.                    \tag{27}
~~~

On substituting \(u=1+v\), all 46 coefficients are positive integers;
the constant is 192. The checker verifies this expansion exactly.
Hence the rate for this displayed size is strictly above \(3/4\), and
in particular above \(1/\sqrt2\).

These arguments exclude the half threshold for the exact profiles and
height list (7)--(9). They do not establish a floor for arbitrary
one-species physical gauges or for the full two-species kernel.

## 6. Reproduction

Run

```sh
python3 tools/check_schur_one_species_interval.py \
  --output evidence/schur-one-species-interval.json
```

The checker verifies the binary-digit identities, both exact Smith
profiles, internal regularity, a complete typed endpoint, and the two
smaller integer pressure inequalities. Small cases also have independent
block-Toeplitz rank checks. Add `--large-certificate` to replay the
63-mode rational witness. The optional comparison uses large integers
without floating-point logarithms.
Its [saved report](../evidence/schur-one-species-interval-large.json)
encodes the rational occupations by their root and exponent to keep the
certificate compact.

The analytic argument in Section 5 is an ordinary proof. Finite numerical
checks do not replace its eta transformation or tail comparison.
