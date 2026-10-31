# Reflection-pair compression of the Pascal Gaussian

## Theorem and scope

Let `M>=4` be a power of two. On every boundary of the strict Pascal
Gaussian, retain exactly one occupied mode from each reflection pair
`{i,M-1-i}`. The resulting polynomial crossing has physical alphabet

\[
q=2^{M/2},
\]

instead of `2^M`. It has a literal swap endpoint. It is itself a
Gaussian, multiplied by the common monomial `w^(M^2/2)`, although
its normalized Gaussian kernel is Laurent.

There are explicit natural gauges in which both complete profiles
have the same `M` positive factors

\[
2,6,10,\ldots,4M-2.
\]

The two scalar shifts are `+M^2/2` and `-M^2/2`, respectively. The
physical energies and profile partitions give **exactly the original
M-mode Gaussian rate at twice the logarithmic scale**. Thus this is
an exact square-root reduction of the alphabet, with no new saving
exponent. It supplies a smaller polynomial seed for further searches.

The [constructor](../tools/reflection_dual_rail_gaussian.py) implements
the displayed formula. Its complete coefficients agree with independent
selected minors of the original strict kernel at `M=4,8`: respectively
256 and 65,536 coefficient comparisons. The actual polynomial degrees
are 24 and 96. These finite comparisons support the implementation;
the uniform result follows from the proof below.
The [independent checker](../tools/check_reflection_dual_rail_compression.py)
also compares original weighted cut profiles, checks the one-particle
formula through M=128, and verifies the complete pressure identity.
The [recorded evidence](../evidence/reflection-dual-rail-compression.json)
states these controls separately from the all-size proof.

## 1. The strict kernel preserves a reflection form

Put `M=2m`, `t=w^2`, and `u=1+w`. For a power-of-two size `r`, define

\[
(G_r)_{ij}=t^{i\mathbin{\&}j},\quad
(H_r)_{ij}=t^{r-1-(i\mathbin{|}j)},\quad
(J_r)_{ij}=\mathbf1_{i+j=r-1}.
\]

Then `H_r=J_rG_rJ_r` and
`G_rH_r=(1+t)^(r-1) I`. The strict one-particle matrix is

\[
D=\begin{pmatrix}uG_M&w^{2M}I\\I&uH_M\end{pmatrix}.
\]

Conjugating by `diag(I,w^M I)` gives

\[
D'=\begin{pmatrix}uG_M&w^MI\\w^MI&uH_M\end{pmatrix}.
\tag{1}
\]

Every retained row subset contains `m` modes from its second physical
block, and so does every retained column subset. The two conjugating
monomials therefore cancel in each selected minor. The restricted
coefficient tensors of `D` and `D'` are identical.

The matrix (1) preserves the alternating form
`diag(J_M,J_M)`. Indeed,
`G_M J_M G_M=(1+t)^(M-1) J_M`, the analogous identity holds for
`H_M`, and `G_M J_M=J_M H_M`. The diagonal products are

\[
u^2(1+t)^{M-1}+w^{2M}=(1+t)^M+t^M=1,
\]

and the off-diagonal products vanish in characteristic two.

It also preserves the quadratic form whose polar form is this
reflection form, with zero value on each coordinate vector. For a
column of `G_M`, each reflected pair product is the same monomial
`t^j`; there are `M/2` such terms. Their sum is zero because
`M/2` is even. The same argument applies to `H_M`. The identity
blocks in (1) give no pair contribution. Since the polar form is
preserved and every column is isotropic, the quadratic form is
preserved on all vectors.

## 2. Selected minors form a homology representation

This step gives a direct explanation of the characteristic-two
cancellations. On a space with symplectic pairs `e_i,f_i`, put

\[
\Omega=\sum_{i=1}^{M}e_i\wedge f_i,
\qquad d(a)=\Omega\wedge a.
\]

In characteristic two, `d^2=0`. For a single pair, the differential
maps the vacuum to `e_i wedge f_i`, while the two single-occupation
vectors represent the homology. Tensoring these elementary complexes
shows that the full homology is concentrated in exterior degree `M`,
with basis the `2^M` products selecting one vector from every pair.
An acyclic factor stays contractible under tensor product, which
proves this assertion without a dimension assumption.

A symplectic matrix preserves `Omega` and hence acts on this homology.
The coefficient of an admissible basis vector in an image is its usual
exterior minor. Boundaries have no admissible coefficient. Therefore
the matrix of the homology action is exactly the matrix of selected
minors. In particular, this selection is multiplicative for symplectic
matrices, even though an arbitrary submatrix selection is not.

This homology realization is standard in characteristic-two
representation theory; see, for example, Peter Sin's
[exterior-algebra and spin-module discussion](https://people.clas.ufl.edu/sin/files/mobile.pdf),
slides 16--17. The elementary argument above is all that is needed
for the present construction; no irreducibility assertion is used.

## 3. Orthogonal elimination gives a Gaussian formula

Order the low-half coordinates first and their reflected high-half
partners second. In this isotropic splitting, the upper-left block
of (1) is

\[
A=\begin{pmatrix}uG_m&t^mI\\t^mI&u t^mH_m\end{pmatrix},
\qquad
\det A=t^{m^2},
\qquad
A^{-1}=\begin{pmatrix}uH_m&I\\I&u t^{-m}G_m\end{pmatrix}.
\tag{2}
\]

The identities follow by block multiplication and
`u^2 G_mH_m=(1+t)^m=1+t^m`.

Define `R_ij=t^(i & ~j)` and `L=R^T`, with indices less than `m`.
The off-diagonal blocks are `B=diag(uR,uL)` and `C=diag(uL,uR)`.
The one-bit product identities

\[
LH_m=H_mR=RG_m=G_mL=(1+t)^{m-1}J_m
\]

give alternating matrices

\[
P=CA^{-1}=Q=A^{-1}B
 =\begin{pmatrix}
 (1+t^m)J_m&uL\\uR&(1+t^{-m})J_m
 \end{pmatrix}.
\tag{3}
\]

The orthogonal Gauss factorization is

\[
D'=\begin{pmatrix}I&0\\P&I\end{pmatrix}
\begin{pmatrix}A&0\\0&A^{-\mathsf T}\end{pmatrix}
\begin{pmatrix}I&Q\\0&I\end{pmatrix}.
\tag{4}
\]

Its homology action can be computed directly. A lower alternating
shear acts by a creation Gaussian with pair coefficients `P_ij^2`:
its transition coefficients are principal determinants of `P`, equal
to the squares of their Pfaffians. An upper shear gives the analogous
annihilation Gaussian. Complementary minors show that the middle
factor acts as

\[
\det A\;\Lambda\bigl((A^{-\mathsf T})^{\circ2}\bigr),
\]

where the superscript denotes entrywise squaring. These statements
follow from minors, including all signs, since signs coincide in
characteristic two. Multiplicativity from Section 2 then gives the
complete Gaussian coefficient tensor.

## 4. An explicit Laurent kernel on half as many physical modes

Each physical letter now has `m` logical bits. Bit zero selects the
low member `i` of its reflection pair; bit one selects `M-1-i`.
Use physical variable groups `(V,Y,U,X)`, each indexed by `0<=i<m`.
The restricted tensor is `w^(2m^2)` times the Gaussian with the
following alternating kernel. Each listed off-diagonal block also
specifies its transpose.

| Pair of groups | Kernel block |
| --- | --- |
| `V,V` and `U,U` | `(1+w^(2M)) J_m` |
| `Y,Y` and `X,X` | `(1+w^(-2M)) J_m` |
| `V,Y` and `U,X` | `(1+w^2) L^(circ 2)` |
| `V,U` | `(1+w^2) H_m^(circ 2)` |
| `V,X` and `Y,U` | `I` |
| `Y,X` | `(1+w^2) w^(-2M) G_m^(circ 2)` |

Here `G_m,H_m,L` still use `t=w^2`. A Gaussian coefficient on an even
set of variables is the Pfaffian of its kernel restriction; odd
coefficients vanish.

Although the displayed kernel has negative powers, the complete
restricted tensor is polynomial: it is the selected-minor tensor of
the original polynomial matrix `D`. Each minor has order `M`, giving
the conservative degree bound `2M^2`. No denominator is discarded.
At `w=1`, every listed block except `VX` and `YU` vanishes. The
remaining Gaussian is the literal crossing on all logical letters.

## 5. Both natural profiles are explicit

Put

\[
a_i=2M-2-4i\qquad(0\le i<m).
\]

At zero, give `X_i,Y_i` height `a_i` and `U_i,V_i` height zero.
At infinity, give `U_i,V_i` height `a_i` and `X_i,Y_i` height zero.
These are exactly the original natural heights restricted to the
reflection alphabet, after normalizing the zero letter on each leg.

Every kernel coefficient internal to either cut is regular in its
local parameter. For the reflected same-leg terms this follows from
`a_i+a_(m-1-i)=2M`. The mixed internal terms have minimum order

\[
2M-4-4(i\mathbin{|}j)\ge0.
\]

Their Gaussian multiplication maps and inverses are therefore regular.
Removing them preserves the complete Smith profiles.

Let `z` be the local parameter, `D_0=diag(z^(a_i))`, and

\[
T_{ij}=(1+z^2)z^{4(i\mathbin{\&}\neg j)}.
\]

Then both remaining one-particle cross matrices have the form

\[
\begin{pmatrix}D_0T&D_0\\D_0&D_0T\end{pmatrix},
\tag{5}
\]

up to the stated ordering of the boundary groups. The product `D_0T`
is symmetric. Both `T` and `D_0 T D_0^(-1)` are regular, and

\[
T^2=(1+z^{2M})I.
\tag{6}
\]

For (6), each bit factor of `R_m` squares to a scalar matrix;
their product is `(1+z^4)^(m-1) I`. Multiplying by
`(1+z^2)^2` gives `(1+z^4)^m I=(1+z^(2M))I`.

Regular block shears in (5) consequently produce

\[
\operatorname{diag}(D_0,z^{2M}D_0).
\]

Its orders are exactly `2,6,...,4M-2`. The full Gaussian profiles
are all subset sums of these `M` numbers. Finally, the scalar
`w^(M^2/2)` adds `M^2/2` to every zero order and subtracts the same
number from every infinity order. These scalar costs cancel in the
product of the two profile partitions.

## 6. Exact equality of the rate functions

Write `ell=log(s)` and `R=2 rho`. The complete compressed profile
product and the two physical partitions are

\[
\log B_{\rm comp}(\ell)
 =2\sum_{j=0}^{M-1}\log(1+e^{-2(2j+1)\ell}),
\]

\[
\log A_{\rm comp}(\ell/R)
 =2\sum_{i=0}^{m-1}
   \log(1+e^{-4(M-1-2i)\ell/R}).
\]

Define the pressure by `F=log B-R log A`. Pairing the positive and
negative physical energies of the full M-mode interval cancels its
common ground energy and gives directly

\[
F_{\rm comp}(R,\ell)
 =\tfrac12 F_{\rm full}(R,2\ell).
\tag{7}
\]

Thus the positive Gibbs roots agree at the corresponding scales.
The raw partition functions contain different common ground factors;
(7) includes their exact cancellation. It is not an assertion that
unnormalized partition functions always agree by a bare square root.

For illustration, the compressed `M=8` polynomial has sixteen letters
per boundary and natural optimized rate approximately
`0.716653659729`. The full eight-mode alphabet has 256 letters and
the same rate. This remains above `1/sqrt(2)`.

## 7. Further searches

The explicit Gaussian formula makes all admissible independent logical
mode gauges accessible through matrices of order `M`. Removing internal
Gaussians is allowed only inside the linear cone where every internal
coefficient remains regular. The
[mode-gauge search](../tools/search_reflection_gaussian_gauges.py)
uses this cone, complete one-particle profiles, and both physical
partitions. Its numerical searches are separate from the theorem.

Arbitrary polynomial perturbations vanishing at `w=1` can also be
applied to the complete polynomial seed. This uses the proved
[endpoint-preserving constructor](../tools/ramified_polynomial_crossings.py).
A favorable perturbed pressure still needs a strict independent
certificate before it can support a stronger factoring exponent.
