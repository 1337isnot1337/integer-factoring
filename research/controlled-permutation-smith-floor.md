# The Feistel Smith obstruction extends to controlled monomial matrices

## Result and predecessor

The full-cycle Feistel obstruction in
[PR #31, pinned revision 314954eb, §4](https://github.com/1337isnot1337/integer-factoring/blob/314954eb/research/modified-transpose-root-target.md)
extends from cyclic translations to arbitrary controlled permutations.
The number of independent control matrices, including linearly dependent
families of distinct permutations, is the relevant parameter.
Arbitrary nonzero entry weights are also allowed, over every field of
characteristic two, when the full cycle has total weight one.

For every power of two \(q\ge8\), the entire family below has basic Smith
certificate
\[
B(s)\ge q^3/\sqrt2,\qquad
\rho\ge17/24>1/\sqrt2.
\]
Thus this certificate cannot improve the half-rate factoring exponent
for these paths. The construction and the first-order Smith argument are
from the pinned predecessor. The extension here is the span-based fiber
bound for arbitrary control permutations.

This is a restriction on the stated sufficient certificate. It does not
bound factoring algorithms, arbitrary coefficient recombinations, or
unrelated auxiliary constructions.

## 1. Controlled permutations and the existing nilpotent path

Let \(q=2^a\), \(n=q^2\), and choose an arbitrary permutation \(L_j\) of
the \(q\) labels for each label \(j\). Regard each \(L_j\) as its
permutation matrix. On matrix units, define
\[
U(i,j)=(j,L_j(i)).
\tag{1}
\]
Assume that \(U\) is one cycle of length \(n\).
Over any field of characteristic two, \(S=I+U\) is then regular nilpotent:
the characteristic and minimal polynomials of this cyclic permutation
are both \((X+1)^n\).

Use the predecessor's reshuffling
\[
\Phi(M)(x,y,u,v)=M(E_{vy})_{xu}.
\]
The established nilpotent construction is
\[
K(z)=\Phi\left(\sum_{j=0}^{n-1}z^{n-1-j}S^jU\right).
\tag{2}
\]
It has endpoint \(K(1)=E\) and column orders at zero
\[
(0,\underbrace{n,\ldots,n}_{n-1}).
\tag{3}
\]
See the predecessor's
[nilpotent construction](https://github.com/1337isnot1337/integer-factoring/blob/314954eb/research/nilpotent-flag-crossings.md)
for the full endpoint, Smith, and effectiveness proofs. Reversing (2)
gives the row series
\[
H(w)=\Phi\bigl(U(I-wS)^{-1}\bigr)
     =\Phi(U)+w\Phi(U+U^2)+O(w^2).
\tag{4}
\]

Put
\[
r=\dim_F\mathrm{span}\{L_j:0\le j<q\}.
\tag{5}
\]
This is at most the number of distinct controls, and can be strictly
smaller. The row matrix of \(\Phi(U)\) is supported on rows \(x=y\);
the row there is the vectorization of \(L_x\). Its rank is exactly \(r\).

## 2. The new fiber bound uses the control span

The row matrix of \(\Phi(U^2)\) has exactly one nonzero entry in each row
\((x,y)\), at column
\[
\Psi(x,y)=\bigl(L_x(y),L_y^{-1}(x)\bigr).
\tag{6}
\]
Consequently
\[
\mathrm{rank}\Phi(U^2)=|\mathrm{im}\Psi|.
\tag{7}
\]

Fix a prescribed output \((b,c)\) of (6). Every possible preimage has
\[
y=L_x^{-1}(b),\qquad x=L_y(c).
\tag{8}
\]
Once \(y\) is known, the second equation fixes \(x\). Distinct possible
values of \(y\) give distinct coordinate vectors
\[
e_y=L_x^{-1}e_b=L_x^{\mathsf T}e_b.
\]
These coordinate vectors are linearly independent. They all belong to the
image, under evaluation at \(e_b\), of the \(r\)-dimensional matrix space
\(\mathrm{span}\{L_x^{\mathsf T}\}\). Therefore each fiber of
\(\Psi\) has size at most \(r\), even when distinct control matrices are
linearly dependent. Thus
\[
\boxed{\quad
\mathrm{rank}\Phi(U^2)\ge\lceil n/r\rceil .
\quad}
\tag{9}
\]

For translations \(L_x(y)=y+f(x)\), this recovers the predecessor's
counting argument. Equation (9) also applies to noncommuting control
permutations, which cannot all become cyclic translations under one
simultaneous relabeling.

## 3. Transfer to the complete Smith obstruction

Let \(\mu_i\) be the row orders of (4). There are exactly \(r\) zero
orders. Eliminate the constant pivots in \(\Phi(U)\). The first remaining
coefficient is the compression of \(\Phi(U^2)\) between a kernel and a
quotient of codimension \(r\). Deleting these two sets of directions lowers
rank by at most \(2r\). Hence
\[
\#\{i:\mu_i=1\}
\ge\max(0,\lceil n/r\rceil-2r).
\tag{10}
\]
This is the first-order argument from the predecessor, now supplied with
the stronger scope of (9).

The basic two-spectrum certificate is
\[
B(s)=\left(\sum_i s^{n-1-\mu_i}\right)
       \left(1+(n-1)s^{-n}\right),\qquad s\ge1.
\tag{11}
\]
A zero order contributes
\[
a_0(s)=s^{n-1}+(n-1)s^{-1}\ge n.
\]
An order-one pivot contributes
\[
a_1(s)=s^{n-2}+(n-1)s^{-2}\ge n-1,\qquad n\ge4.
\]
The latter inequality follows by minimizing at
\(s^n=2(n-1)/(n-2)\), followed by Bernoulli's inequality
\[
\left(1+\frac2{n-2}\right)^{n/2}
\ge1+\frac n{n-2}
=\frac{2(n-1)}{n-2}.
\]
Thus the same minimization as in the predecessor gives
\[
B(s)\ge nr+(n-1)\max(0,n/r-2r)
       \ge n^{3/2}/\sqrt2=q^3/\sqrt2.
\tag{12}
\]
For the last step, if \(r\ge\sqrt{n/2}\), the first term suffices.
Below that point the expression is
\(n(n-1)/r-(n-2)r\), which decreases with \(r\) and reaches
\(n^{3/2}/\sqrt2\) at the boundary.

For \(q\ge8\), (12) gives
\[
\rho=\frac{\log B(s)}{4\log q}
\ge\frac34-\frac{\log2}{8\log q}
\ge\frac{17}{24}>\frac1{\sqrt2}.
\tag{13}
\]
The strict final comparison is \(2\cdot17^2>24^2\).

Exchanging the cuts does not help this certificate. The column matrix
of (2) has invertible leading coefficient \(U\), so all its infinity
degrees are \(n-1\). The sum of the row orders at zero is at most the
determinant degree \(n(n-1)\). Arithmetic--geometric mean therefore gives
the opposite certificate at least \(n^2\).

The \(q=2,4\) cases are already excluded by the general degree budget in
[PR #30's Smith interface, §6](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/smith-profile-crossing-transfer.md).
That existing small-alphabet result is separate from (13).

## 4. Arbitrary monomial weights obey the same obstruction

The extension permits more than permutations. Let each \(L_j\) be an
invertible monomial matrix: for a permutation \(\pi_j\) and nonzero
weights \(\ell_j(i)\),
\[
L_j e_i=\ell_j(i)e_{\pi_j(i)}.
\]
Replace (1) by the weighted operator
\[
U(E_{ij})=\ell_j(i)E_{j,\pi_j(i)}.
\tag{14}
\]
Assume that its underlying permutation is one cycle of length \(n=q^2\),
and that the product of the \(n\) weights around that cycle is one.
Its characteristic and minimal polynomials are again \((X+1)^n\):
successive powers of a coordinate vector give a scaled coordinate basis,
and \(U^n=I\). Thus the construction and ideal column orders remain valid.

The first reshuffling still has rank
\(r=\dim_F\mathrm{span}\{L_j\}\).
The second has one nonzero entry in each row; its support map is (6)
with the permutations \(\pi_j\). Its rank is therefore the image size
of that support map.

For a fixed output label \(b\), the vectors \(L_j^{-1}e_b\) and
\(L_j^{\mathsf T}e_b\) need not be equal. They are nonzero multiples
of the same coordinate vector. Consequently distinct possible inverse
images still give linearly independent vectors in the evaluation of
\(\mathrm{span}\{L_j^{\mathsf T}\}\) at \(e_b\). The fiber size
is at most \(r\), exactly as before. No unitary condition, common phase
group, or field automorphism is required.

Equations (9)--(13), including the opposite-cut observation, therefore
hold for every such weighted cycle. This covers arbitrary monomial
weights in \(F_4\), as well as in other characteristic-two fields.

## 5. Exact controls and scope

The [standard-library checker](../tools/check_controlled_permutation_smith.py)
checks the two literal reshufflings, every fiber, the span argument, and
the first three row jets. It verifies their Toeplitz ranks using both
elimination orders. It also checks every entry of the complete crossing
endpoint, and independently recovers the ideal column orders in the
two- and four-label controls.

Its finite cases include:

- Every full-cycle member of the four-label family with a cyclic first
  control, an arbitrary second control permutation, and an arbitrary
  nontrivial binary control subset.
- Controlled-transposition cycles at \(q=2,4,8,16\). Their complete row
  orders are \((0,1,1,2)\) at two labels and
  \((0,0,1^{q^2-3},2)\) in the other three cases.
- A full 64-cycle with two noncommuting controls. Its row first-order
  multiplicity is 44, so it is outside the translation hypothesis but
  remains far from a useful spectrum.
- A full 64-cycle with four distinct control permutations of span
  dimension three. This checks the linear-dependence issue in (9).
- Weighted four- and eight-label cycles over \(F_4\). Their complete
  crossing endpoints and first three row jets are checked directly.
  The extension-field ranks are independently represented by binary
  multiplication matrices.

Run:

~~~bash
python3 tools/check_controlled_permutation_smith.py \
  --output evidence/controlled-permutation-smith.json
~~~

The [recorded evidence](../evidence/controlled-permutation-smith.json)
states the finite scope and checker hash. The all-size conclusion is
the proof of (9)--(13), not an extrapolation from these controls.

The scope extension and its checks used substantial OpenAI assistant
support. The predecessor's construction and Feistel obstruction retain
their attribution. No stronger factoring exponent is claimed.
