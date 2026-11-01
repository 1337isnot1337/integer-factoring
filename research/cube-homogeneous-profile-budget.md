# A profile floor for every homogeneous binary cube tensor

This note removes the Gaussian assumption from a scoped obstruction to the
[three-axis cube transfer](cube-crossing-transfer.md). It concerns the
common-scale, unweighted Smith-profile estimate with the same binary type
on every physical leg. It does not give a rank lower bound for individual
coefficient tensors.

**Theorem.** Let an eight-leg binary tensor have a common polynomial
deformation

\[
 T_b(x,y,z),\qquad b\in\{0,1\}^{8},
\]

over a field, with every entry homogeneous of degree \(h|b|\), where
\(h>0\). Suppose its value at \((1,1,1)\) is the product of equality
tensors on the four antipodal pairs of cube vertices. Then the
common-scale unweighted Smith-profile estimate, combined with a common
binary physical type, has

\[
 \rho>\frac35>\frac1{\sqrt3}.
\]

This includes arbitrary non-Gaussian entries and permits odd-occupation
entries. It extends the bound from
[Gaussian endpoint minors](cube-gaussian-profile-budget.md) to every
homogeneous binary tensor in this method. Multiple physical modes, larger
alphabets, unequal physical types, and stronger coefficient-rank estimates
are outside the binary conclusion here.

## 1. Shared endpoint minors

Label the cube vertices by \(v\in\{0,1\}^3\), with antipode
\(\bar v=v\mathbin\oplus111\). Let \(M_a\) be the \(16\times16\)
flattening across coordinate \(a\). Its columns can be ordered by the
antipodes of its row vertices. The endpoint matrix is then the identity.

Regard \(M_a\) as a matrix over the power-series field for coordinate
\(a\), with rational functions in the other two coordinates as
coefficients. Write its ascending Smith orders as

\[
 0=s_{a,1}\le s_{a,2}\le\cdots\le s_{a,16},\qquad
 \gamma_{a,k}=\sum_{i=1}^{k}s_{a,i}.
\]

All orders are finite and nonnegative. The constant empty entry is
exactly one, so the first order is zero. For any nonzero \(k\)-minor
\(D\), the determinantal-ideal characterization gives
\(\gamma_{a,k}\le v_a(D)\).

The same two-by-two minor belongs to all three cuts: retain only the
empty configuration and the occupied vertex on each side of the pair
\(000,111\). Its determinant is

\[
 D_2=T_\varnothing T_{\{000,111\}}
     -T_{\{000\}}T_{\{111\}}.
\]

It is homogeneous of degree \(2h\) and equals one at the endpoint. A
nonzero homogeneous polynomial of degree \(d\) satisfies
\(v_x+v_y+v_z\le d\). Consequently,

\[
 \sum_a s_{a,2}\le2h. \tag{1}
\]

For two distinct axes \(a,b\), let \(c\) be the omitted axis. Use the
two antipodal pairs with representatives \(000,e_c\). Both
representatives lie on the zero side of both cuts. Allow their two
occupation bits to vary independently on each side, and fix the other
bits to zero. This selects the same four-by-four minor in both cuts,
up to transposition or permutations. Its endpoint is the identity.

The sum of the occupation numbers over its four row configurations is
four, and the same holds for its columns. Every determinant term
therefore has degree \(8h\). Thus

\[
 \gamma_{a,4}+\gamma_{b,4}\le8h.
\]

Adding the three axis-pair inequalities gives

\[
 \sum_a\gamma_{a,4}\le12h. \tag{2}
\]

The smaller shared minor with row and column configurations
\(00,10,01\) is equally useful. Its endpoint is the three-by-three
identity and its determinant has degree \(4h\). The same axis-pair
argument gives

\[
 \sum_a\gamma_{a,3}\le6h. \tag{3}
\]

Finally, on any one cut, choose the \(k\) row configurations of lowest
occupation and their antipodally paired columns. The selected minor
has identity endpoint. If \(e_1\le\cdots\le e_{16}\) is the list of
occupation numbers of all four-bit words, every determinant term has
degree \(2h\sum_{i\le k}e_i\). Therefore, at every cut,

\[
 \gamma_{a,k}\le2h\sum_{i=1}^k e_i,
 \qquad
 (e_i)=(0,1,1,1,1,2,2,2,2,2,2,3,3,3,3,4). \tag{4}
\]

In particular, the full determinant has degree \(64h\), and
\(\sum_a\gamma_{a,16}\le192h\). Selecting arbitrary low-occupation
configurations, including the three-by-three minor, is stronger than
using only full Cartesian alphabets.

These statements use the actual common tensor and shared minors. They
do not optimize three unrelated axis kernels.

## 2. A simultaneous partition-function bound

Set

\[
 \bar s_i=\frac13\sum_a s_{a,i}.
\]

The vector \(\bar s\) is ascending. Equations (1)--(4) imply

\[
 \bar s_1=0,\qquad
 \bar s_2\le\frac{2h}{3},\qquad
 \sum_{i=1}^3\bar s_i\le2h,\qquad
 \sum_{i=1}^4\bar s_i\le4h,\qquad
 \sum_{i=1}^{16}\bar s_i\le64h.
\]

Since the cumulative sums of an ascending vector form a convex
sequence, each intermediate sum is at most the chord between the
available bounds. Combine these inequalities with every prefix bound
in (4). Every prefix sum of \(\bar s\) is then at most the
corresponding prefix sum of

\[
 \sigma=
 \left(0,\frac{2h}{3},\frac{4h}{3},2h,
       \underbrace{4h,\ldots,4h}_{7\text{ terms}},
       \underbrace{6h,\ldots,6h}_{4\text{ terms}},8h\right).
\]

For a convex decreasing function, these prefix inequalities give
\(\sum_i f(\bar s_i)\ge\sum_i f(\sigma_i)\). One can prove this
by moving mass from larger coordinates to smaller coordinates until
the target prefixes are attained: such balancing cannot increase a
convex sum, while increasing the total cannot increase a decreasing
sum. Apply this with \(f(s)=t^s\), for \(0<t<1\).

The function
\(s\mapsto\log\sum_i\exp(s_i\log t)\) is convex. Jensen's inequality
across the three cuts consequently gives

\[
 \begin{aligned}
 \frac13\sum_a\log Z_a(t)
 &\ge \log\sum_i t^{\bar s_i}\\
 &\ge\log\left(1+t^{2h/3}+t^{4h/3}+t^{2h}
                   +7t^{4h}+4t^{6h}+t^{8h}\right),
 \qquad Z_a(t)=\sum_{i=1}^{16}t^{s_{a,i}}.
 \end{aligned} \tag{5}
\]

Put \(u=t^{h/3}\). The lower bound in (5) is

\[
 Z_*(u)=1+u^2+u^4+u^6+7u^{12}+4u^{18}+u^{24}.
\]

## 3. An exact all-scale entropy certificate

For a common binary type with occupation fraction \(p\), homogeneity
puts every common coefficient of the word tensor power on

\[
 k_x+k_y+k_z=8hmp.
\]

The usual unweighted Smith estimate therefore has pressure

\[
 \rho_{\mathrm{prof}}(p,t)=
 \frac{\sum_a\log_2 Z_a(t)-8hp\log_2t}{12H_2(p)}
 \ge
 \frac{\log_2 Z_*(u)-8p\log_2u}{4H_2(p)}. \tag{6}
\]

Here \(0<p<1\); the endpoint types have zero entropy and do not encode
a growing alphabet. Polynomial factors from coefficient extraction do
not change this pressure.

Let

\[
 P(s)=(1+s^6+s^{12}+s^{18}+7s^{36}+4s^{54}+s^{72})^5
       -(1+s^{10})^{12}.
\]

Its first nonzero term has degree \(6\). Write \(Q(s)=P(s)/s^6\),
which has degree \(354\). The exact retained checker expands

\[
 (1+v)^{354}Q\left(\frac{v}{1+v}\right)
\]

and verifies that all \(355\) integer coefficients are positive; the
smallest is \(5\). Thus \(P(s)>0\) for \(0<s<1\), and setting
\(u=s^3\) gives

\[
 Z_*(u)>(1+u^{10/3})^{12/5}.
\]

The elementary entropy identity

\[
 \max_{0\le p\le1}\{H_2(p)+p\log_2 a\}=\log_2(1+a)
\]

with \(a=u^{10/3}\) now implies

\[
 \log_2 Z_*(u)>
 \frac{12}{5}H_2(p)+8p\log_2u.
\]

Substitution in (6) proves \(\rho_{\mathrm{prof}}>3/5\).
The exact comparison \(3\cdot3^2>5^2\) places this
above the transfer's half-exponent threshold.

## 4. The general alphabet interface

The shared-minor argument applies to an alphabet with nonnegative
energies \(\varepsilon_1,\ldots,\varepsilon_q\), including a zero-energy
letter. Assume each tensor entry is homogeneous of the sum of its eight
letter energies. Let \(S_r(k)\) be the sum of the \(k\) smallest
energies of words of length \(r\), with multiplicity. For the average
Smith cumulative sums over the three cuts, the exact bounds are

\[
 \bar\gamma_k\le\frac23 S_1(k)\quad(k\le q),\qquad
 \bar\gamma_k\le S_2(k)\quad(k\le q^2),\qquad
 \bar\gamma_k\le2S_4(k)\quad(k\le q^4). \tag{7}
\]

Use one antipodal pair shared by all three cuts, two pairs shared by
each axis pair, or four pairs on a single cut. In each case choose an
arbitrary set of the lowest-energy configurations and its paired
columns. Their endpoint identity proves that the determinant is
nonzero. The remaining letters are fixed to the zero-energy letter.

The lower convex hull of all the right-hand-side cumulative bounds
is an upper bound on every feasible averaged cumulative profile.
Its successive slopes therefore give an optimistic Smith partition
function by the same majorization argument. If many configurations
have equal energy, only the endpoints of their cumulative intervals
are needed; convexity supplies every interior inequality.

Equation (7) is a necessary budget for the declared profile method.
Its optimistic slopes need not belong to any actual tensor. A
numerically favorable envelope is a search guide, not a constructive
crossing certificate.

## 5. Verification and next constructions

Run:

```sh
python3 tools/check_cube_homogeneous_profile_budget.py \
  --output /tmp/cube-homogeneous-profile-budget.json
```

The checker reconstructs the shared minor selections, their degrees,
and their endpoint identities. Seeded non-Gaussian tensors include
odd-occupation entries. It verifies every exact positive coefficient
of the all-scale certificate. These finite controls support the
algebraic proof above; they are not a search-based proof of a universal
claim.

This result redirects the constructive search toward larger physical
alphabets and interactions between different physical modes. It does
not establish an improved factoring exponent or alter the current
proven \(\alpha=1/2\) benchmark.
