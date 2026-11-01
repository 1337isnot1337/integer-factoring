# Three uniform jet identities for coherent spin-zero crossings

## Statement

Let `n` be even. Use any selected return-to-zero paths from the
[bounded-spin construction](bounded-spin-crossing.md), and compress
`K(z)^tensor n`, where `K(z)=(1-z)I+zS` and `S` swaps the two bits.
Write the resulting tensor polynomial as

```math
G_n(z)=\sum_{j=0}^n C_jz^j.
```

Over `Q`, and over every field of odd characteristic in which the paired
encoder norms are units, the following are tensor identities:

```math
\boxed{\begin{aligned}
2C_1+nC_0&=0,\\
2\sum_{j=0}^n jC_j-n\sum_{j=0}^n C_j&=0,\\
2C_{n-1}+nC_n&=0.
\end{aligned}}                                    \tag{1}
```

They also hold for the integer companion `H_n=12^n G_n`.
For even `n>=6` the three scalar linear constraints are independent in
every odd characteristic. Therefore the coefficient span has dimension
at most `n-2`. This explains the dimensions four and six observed at
lengths six and eight. It does **not** bound the ordinary cut ranks of
the coefficients or prove an improved factoring exponent.

The argument applies more generally to coherent spin-zero frames, with
any paired duals in their spin-zero span. It does not require a complete
path set or stability of the selected subspace under site permutations.

## 1. Integral spin-zero invariance

Let `E` remove a one and `F=E^T` add a one, summed over all sites.
Every return path vector `v` has weight `n/2` and satisfies `Ev=0`.
Since `[E,F]` acts by `n-2*(n/2)=0`,

```math
\lVert Fv\rVert^2=\langle v,EFv\rangle=0.
```

Thus `Fv=0` over the integers. Over `Q`, upper and lower unipotent
matrices act as the exponentials of `E` and `F`, so both fix `v`.
They generate `SL_2(Q)`. The vector is therefore fixed by `g^tensor n`
for every `g in SL_2`. This is an integral polynomial identity, so it
continues to hold after reduction modulo a prime. In particular it holds
over the algebraic closure of `F_5`, rather than only its finite group
of rational points.

This last distinction matters: finite-group invariance alone would not
justify a statement about all characteristic-five extensions.

## 2. One-site transition matrices are scalar

For two spin-zero vectors `v,w` and a site `i`, define the two-by-two
transition matrix

```math
R_i(v,w)_{ab}=\sum_{\text{all other bits}}v(\ldots,a,\ldots)
                                      w(\ldots,b,\ldots).
```

Apply `g` to `v` and `g^{-T}` to `w`. Both vectors are unchanged, and
the transformations cancel at every contracted site. Consequently
`R_i=gR_i g^{-1}` for every `g in SL_2`. Commuting with both upper
and lower unipotents forces a two-by-two matrix to be scalar. Taking
its trace gives the integral identity

```math
\boxed{2R_i(v,w)=\langle v,w\rangle I_2.}          \tag{2}
```

The proof can first be read over `Q`; both sides are integral, so the
identity reduces in every odd characteristic. It applies to distinct
path vectors as well as diagonal pairs, and hence to arbitrary linear
combinations of paths.

## 3. Derivatives at the three pairings

At `z=0`, all local kernels are the identity pairing `x=u,y=v`.
Replacing one local identity by `S` contracts the two reduced matrices
from (2) with the swapped pairing. The contraction is exactly half of
the original one: `tr(I_2)=2` rather than `tr(I_2)^2=4`.
Replacing it by `S-I` therefore contributes `-C_0/2`. Summing over
the `n` replacement sites proves the first identity in (1).

Exchanging the last two physical legs replaces `z` by `1-z`.
Apply the same calculation at the crossing pairing to obtain
`G'_n(1)=nG_n(1)/2`, which is the second identity in (1).

For the third pairing use

```math
(S-I)_{xyuv}=-\varepsilon_{xy}\varepsilon_{uv},
\qquad \varepsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
```

Since `epsilon in SL_2`, its full tensor power fixes every spin-zero
vector. The leading tensor `C_n` is therefore the product of the two
paired Gram forms, divided by the prescribed dual norms. In particular
it has rank one in the cut `xy|uv` when the path sets coincide.

Replace one factor `S-I` by `I`. The epsilon-contracted one-site
transition matrices are `epsilon*<v,w>/2`, by (2) and spin-zero
invariance. Their local identity contraction gives half the product
of the full Gram forms. There is one fewer minus sign from `S-I`,
so each replacement contributes `-C_n/2`. Summing proves
`C_{n-1}=-nC_n/2`, the third identity in (1). Equivalently this is
the first derivative at zero of the reversed polynomial
`w^n G_n(1/w)`.

For `n>=6`, the first and third equations constrain disjoint pairs
of coefficients. In the middle equation the two interior coefficients
at indices two and three have weights `4-n` and `6-n`, differing by
two. At least one is nonzero in every odd characteristic. The middle
equation is therefore independent of the other two. At `n=4` the
middle equation is their difference, giving the observed three-dimensional
span instead.

## 4. Four-qubit invariant interpretation and next question

Each contracted entry is a homogeneous degree-`n` invariant polynomial
of a local four-qubit tensor: the four path covectors are separately
fixed by `SL_2`. Luque and Thibon describe the characteristic-zero
invariant ring using generators of degrees two, four, four and six in
[*The polynomial invariants of four qubits*](https://arxiv.org/abs/quant-ph/0212069).
This provides a useful representation for the next algebra calculation;
the proof of (1) does not assume that generator theorem in positive
characteristic.

One concrete model is the orthogonal Gram matrix of the local flattening.
Put `K=aI+bS`, `s=a+b`, `t=a-b`, and `J=epsilon tensor epsilon`.
Then `K^T J K J=K^2` has eigenvalues `s^2,s^2,s^2,t^2`. Its first
three characteristic coefficients and the determinant of `K` restrict to

```math
3s^2+t^2,\quad 3s^4+3s^2t^2,\quad
s^6+3s^4t^2,\quad s^3t.
```

The pencil is `s=1,t=1-2z`. These explicit scalar polynomials model the
small coefficient spaces. A complete modular generation argument and,
more importantly, ordinary-rank formulas for the associated tensor
coefficients remain to be supplied. Counting invariants or noting the
triple Gram eigenvalue does not determine cut-rank growth.

## Verification

`make verify-bounded-spin` checks (2) over the integers on all 13,080
one-site transition matrices between return paths of lengths two through
ten. It checks all three identities on the sixteen independent literal
four-site integer polynomials. The optional NumPy target
`make verify-bounded-spin-algebra` checks (1) on every entry of the
complete eight-site coefficient tensor, in addition to its
[primary-algebra certificate](bounded-spin-rank-algebra.md).
The source-pinned receipts are
[constructor controls](../evidence/bounded-spin-crossing.json) and
[rank-algebra controls](../evidence/bounded-spin-rank-algebra.json).
This research used substantial Codex AI assistance. The complete factoring
score remains `1/2`.
