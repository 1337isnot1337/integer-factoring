# Nilpotent superoperators give exact crossing families

## Result and scope

There is an elementary construction whose column Smith orders are known
exactly before computing any tensor ranks. It accepts an arbitrary nilpotent
linear operator on a matrix algebra. A single full Jordan block attains the
abstract column spectrum

\[
(0,n,\ldots,n),\qquad n=q^2,
\]

at polynomial degree `n-1`. The other cut remains an independent obligation.
An explicit two-label instance is unimodular in the row cut, and an explicit
four-label instance comes from a regular nilpotent with a two-channel leading
row coefficient. Neither instance improves the factoring exponent.

The construction and the identities below are ordinary algebra over any
field. The supplied checker uses exact binary arithmetic for the displayed
examples. There is no new Lean proof.

## 1. A basis-free construction

Let `V=Mat_q(F)`, `n=q^2`, and let `S:V->V` be nilpotent of index `m`.
Write `E_vy` for a matrix unit. Identify a linear operator `M` on `V` with
the four-tensor

\[
\Phi(M)(x,y,u,v)=M(E_{vy})_{xu}.
\tag{1}
\]

Then `Phi(Id)` is the crossing `E(x,y,u,v)=[x=v][y=u]`. Its column
flattening, with indices `(x,u)|(y,v)`, is the matrix of `M` followed by
the fixed permutation transposing the input matrix. Consequently it has
the same Smith invariants as `M`.

Define the polynomial family of degree at most `L=m-1` by

\[
K(z)=\Phi\left(\sum_{k=0}^{m-1}z^{m-1-k}S^k(I-S)\right).
\tag{2}
\]

The endpoint is exact:

\[
K(1)=\Phi\left(\sum_{k=0}^{m-1}(S^k-S^{k+1})\right)
     =\Phi(I)=E.
\tag{3}
\]

After reversing the polynomial parameter, its row matrix is

\[
H(w):=w^{m-1}K(w^{-1})
 =\Phi\bigl((I-S)(I-wS)^{-1}\bigr).
\tag{4}
\]

The inverse in (4) is the finite sum `sum_(k=0)^(m-1) w^k S^k`.
No analytic convergence, coefficient division, or interpolation is used.
The finite matrices and their powers give a direct effective constructor.

## 2. Exact column Smith orders

Suppose the nilpotent Jordan blocks of `S` have sizes `b_1,...,b_t`.
The column orders at zero are the union, over these blocks, of

\[
\boxed{\quad m-b_j\text{ once},\qquad
                 m\text{ repeated }b_j-1\text{ times}.\quad}
\tag{5}
\]

To prove this, restrict the operator polynomial inside (2) to a Jordan
block `J_b`. It is

\[
z^m(zI-J_b)^{-1}(I-J_b).
\tag{6}
\]

The matrix `I-J_b` is constant and invertible. The matrix `zI-J_b` has
Smith form `diag(1,...,1,z^b)`: its determinant is `z^b`, and its first
`b-1` rows and last `b-1` columns give a unit minor. Inverting and
multiplying by `z^m` proves (5). Constant changes of the matrix-algebra
basis and the fixed input-transpose permutation do not affect the result.
The sum of all orders is `n(m-1)`.

For a regular nilpotent, `m=n` and there is one block. This gives precisely
`(0,n,...,n)` at degree `n-1`. This attains one of the two ideal spectra
appearing in the [Smith-profile certificate](https://github.com/1337isnot1337/integer-factoring/blob/28105165/research/smith-profile-crossing-transfer.md); it says nothing by
itself about the row orders at infinity.

## 3. The triangular matrix-basis formula

For the regular case choose a Jordan chain `A_0,...,A_(n-1)` with
`S(A_0)=0` and `S(A_i)=A_(i-1)`. Let `D_i` be its trace-dual basis:
`Tr(A_i D_j)=[i=j]`. The trace pairing on a full matrix algebra is
nondegenerate over every field, including when the characteristic divides
`q`.

Put `D_n=0` and `B_j=D_(n-1-j)-D_(n-j)`. Formula (2) becomes

\[
K(z)=\sum_{i+j\le n-1}z^{i+j}A_i\otimes B_j.
\tag{7}
\]

In the constant column bases given by the `A_i` and `D_j`, its coefficient
matrix has diagonal entries `z^(n-1)`, entries
`(1-z)z^(n-1+i-j)` for `i<j`, and zero entries below the diagonal. Calling
this matrix `T(z)` and the upper shift `J_n` gives the direct identity

\[
T(z)(zI-J_n)=z^n(I-J_n),
\tag{8}
\]

which is also a short proof of the regular spectrum.

Changing the cyclic vector of the Jordan chain does not change `K`.
Indeed

\[
\sum_i A_i\otimes D_{i+r}=\Phi(S^r),
\tag{9}
\]

so every coefficient in (7) depends only on the fixed physical operator
`S`. In particular, replacing the cyclic vector by `f(S)v`, with `f(0)`
nonzero, cannot improve either cut profile. Such a replacement is a
change of description of exactly the same four-tensor polynomial.

## 4. A two-label unimodular control

Over `F2`, take

\[
A_0=I_2,\quad A_1=E_{01},\quad A_2=E_{10},\quad
A_3=E_{00}+E_{10}.
\]

In a polynomial bit mask, bit `j` records the coefficient of `z^j`.
The row matrix of (7) is

\[
\begin{pmatrix}
13&6&10&12\\
12&5&8&10\\
12&8&5&6\\
0&12&12&13
\end{pmatrix}.
\tag{10}
\]

Exact polynomial elimination gives determinant one. The row orders at
infinity, using degree three, are `(0,0,1,11)`. The column orders at
zero are `(0,4,4,4)`. Thus row unimodularity and the ideal column spectrum
are compatible, but the row orders here are highly uneven. They are
not the ideal row spectrum `(0,4,4,4)` that would give degree spectrum
`(3,-1,-1,-1)`.

## 5. A regular nilpotent from a two-channel transpose map

For commuting matrices `A,B` and an invertible matrix `P`, define

\[
U(X)=P\left[X^{\mathsf T}(I-AB)+AX^{\mathsf T}B\right]P^{-1}.
\tag{11}
\]

It fixes the identity and preserves trace. The first assertion is
`I-AB+AB=I`; the second follows by cycling the factors inside the trace
and using `AB=BA`. Its tensor `Phi(U)` has row rank at most two, because
(11) is a sum of two maps of the form `X -> C X^T D`.
The additional condition that `S=I-U` be regular nilpotent is a genuine
condition, checked separately in the example below.

Let `q=4`, work over `F2`, let `J` be the upper shift, and let `R` be the
antidiagonal permutation matrix. Set

\[
P=R+E_{00},\quad P^{-1}=R+E_{33},\quad
A=J+J^3,\quad B=I+J^2.
\tag{12}
\]

The four row-major binary matrix masks are respectively
`4681,37448,2122,33957`. Direct exact calculation gives
`S^15(E_00)=I` and `S(I)=0`. The sixteen vectors in this chain are
independent, so `S` is regular nilpotent on `Mat_4`.
The chain is an explicit basis for (7), and its last trace-dual vector
is `I`. In this case `K(0)` is the row identity and `K(1)=E_4`.

The complete checked spectra are

\[
\begin{aligned}
\mu_\infty&=(0,0,\underbrace{1,\ldots,1}_{12},9,11),\\
\nu_0&=(0,\underbrace{16,\ldots,16}_{15}).
\end{aligned}
\tag{13}
\]

The first row coefficient has rank two. The profile nevertheless gives
no improved factoring bound. In the Smith criterion its two Laurent
sums are

\[
U_*(s)=2s^{15}+12s^{14}+s^6+s^4,\qquad
V_*(s)=1+15s^{-16}.
\tag{14}
\]

For `1<=s<=3/2`, their product is at least `180s^(-2)>=80`.
For `s>=3/2`, it is at least `2s^15>80`. Hence this certificate stays
above `80>64>4^(2sqrt(2))` for every `s>=1`. Numerical optimization
is unnecessary for the failure claim.

## 6. Finite determinant roots do not repair this example

A fractional linear parameter change can preserve zero and one while
moving a finite point `gamma!=0,1` to infinity. Clearing its scalar
denominator and normalizing at one preserves the crossing endpoint and
the column orders at zero. Its new row orders at infinity are the local
row Smith orders of the original `K` at `gamma`. This is a valid extra
search direction; checking infinity alone need not find the best profile.

For (12), fraction-free exact elimination gives

\[
\det K_{\rm row}(z)=f_7(z)^4 f_{14}(z)^8 f_{17}(z)^4,
\tag{15}
\]

where the following polynomials are irreducible over `F2`:

\[
\begin{aligned}
f_7&=z^7+z^5+z^4+z^3+1,\\
f_{14}&=z^{14}+z^{12}+z^9+z^8+z^6+z^5+z^4+z^2+1,\\
f_{17}&=z^{17}+z^{12}+z^8+z^6+z^5+z^3+1.
\end{aligned}
\]

The checker verifies both their product and their irreducibility. Every
finite root therefore has total local Smith order at most eight. For
sixteen nonnegative integer orders of total at most eight, convexity
(or moving one unit from an order at least two to an order zero) gives

\[
\sum_{i=1}^{16}s^{-\mu_i}\ge8+8s^{-1}\quad(s\ge1).
\]

The existing infinity spectrum satisfies

\[
2+12s^{-1}+s^{-9}+s^{-11}\le8+8s^{-1},
\]

since `6-4/s-s^(-9)-s^(-11)>=0` for `s>=1`.
Thus every finite-point move has a Laurent row sum at least the one
already in (14), with the same column sum. No such move improves this
particular witness, even after a field extension.

## 7. A first-order filter for the row spectrum

Write `r=rank Phi(U)` and `R_2=rank Phi(U^2)`, where `U=I-S`.
By (4),

\[
H(w)=\Phi(U)+w\Phi(U-U^2)+O(w^2).
\]

There are exactly `r` zero orders. After constant invertible row and
column operations make this leading matrix `diag(I_r,0)`, the first
coefficient on its remaining block is the compression of
`-Phi(U^2)` between two codimension-`r` spaces. Its rank is at least
`max(0,R_2-2r)`. Eliminating the unit block changes that remaining
coefficient only at order two. Therefore

\[
\#\{i:\mu_i=1\}\ge\max(0,R_2-2r).
\tag{16}
\]

For (12), the two independent ranks are `r=2` and `R_2=14`; the
actual spectrum has twelve order-one entries. The bound in (16)
guarantees ten. A small leading rank alone therefore leaves a large
first-order term to control. This is a search diagnostic for the
given family, not a lower bound on arbitrary crossings.

## 8. The exact profile for a single transpose channel

In characteristic two suppose `U(X)=P X^T P^{-1}` and `S=I+U` is
nilpotent. Put `n=q^2` and use its actual nilpotence index in (2). Then

\[
\boxed{\quad
\mu_\infty=
\begin{cases}
(0,1^{n-1}),&q\text{ odd},\\
(0,1^{n-2},2),&q\text{ even}.
\end{cases}\quad}
\tag{17}
\]

Here exponents denote multiplicities. This statement concerns the complete
row Smith profile, not just its first coefficient.

Set `C=P P^{-T}`. Then `U^2(X)=CXC^{-1}`. Nilpotence of `I+U` makes
`Ad_C` unipotent. Over an algebraic closure, all eigenvalue ratios of `C`
are therefore one, so `C` has a single eigenvalue. Also `C` is similar to
`C^{-T}`, by conjugation with `P`; that eigenvalue equals its inverse and
hence is one in characteristic two. Thus `C` is unipotent and
`Tr(C^{-1})=q` in the field.

Write

\[
R=\Phi(U)=fg^{\mathsf T},\quad A=\Phi(U^2),\quad
B=\Phi(U^3)=f_3g_3^{\mathsf T},
\]

where, using row-major vectorization,

\[
f=\operatorname{vec}P,\quad g=\operatorname{vec}(P^{-T}),\quad
f_3=\operatorname{vec}(CP),\quad g_3=\operatorname{vec}((CP)^{-T}).
\]

The matrix `A` sends a matrix `Z` to `C Z^T C^{-T}` and is its own
inverse. Direct contraction gives

\[
g^{\mathsf T}A^{-1}f=\operatorname{Tr}(C^{-1})=q,\qquad
g^{\mathsf T}A^{-1}f_3=q,\qquad
g_3^{\mathsf T}A^{-1}f=q.
\tag{18}
\]

The leading matrix `R` has rank one. Its first-order compression from
Section 7 has rank `n-1` when the first contraction in (18) is nonzero,
and rank `n-2` when it is zero. This already proves the odd case and
leaves only one exceptional order, at least two, in the even case.

For completeness, expand (4) through degree two:

\[
H(w)=R+w(R+A)+w^2(R+B)+O(w^3).
\tag{19}
\]

In the even case factor out `wA` over Laurent series. The remaining
matrix is

\[
I+(w^{-1}+1+w)A^{-1}fg^{\mathsf T}
  +wA^{-1}f_3g_3^{\mathsf T}+O(w^2).
\]

The determinant lemma for a sum of two rank-one matrices, and (18),
show that this determinant has constant term one and no pole. The omitted
`O(w^2)` term cannot change that constant: the only possible negative
power comes from a single rank-one `w^{-1}` term. Consequently
`ord_w det H=n`. The orders already found are one zero, `n-2` ones,
and one order at least two. Their sum is `n`, which forces the last
order to be exactly two and proves (17).

The small rank of `Phi(U)` alone therefore does not produce a broad
spread of large row orders in this class. Coupled transpose channels
and more general nilpotent superoperators are separate constructions.

## Verification

Run the standard-library checker:

```sh
python3 tools/check_nilpotent_flag_crossings.py \
  --output evidence/nilpotent-flag-crossings.json
```

It verifies the two displayed bases and endpoints, both complete spectra,
independent block-Toeplitz rank checks, determinant identities and
irreducible factors, cyclic-vector invariance controls, the two first-order
ranks, the single-transpose profile, and nonregular transpose examples with the Jordan orders in (5).
Power-series elimination uses the full degree bound `nL+1`; the independent
Toeplitz checks also recover every displayed order.

The construction, proofs, and exact computations used substantial OpenAI
assistant support. No claim above changes the proved factoring score `1/2`.
