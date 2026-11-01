# Ordered qutrit phases and rank-one cold boundaries

Clock and shift matrices give an explicit ordered crossing whose entry
weights are quadratic phases. The auxiliary boundary can be normalized
to preserve the rank-one zero cut at every word length. Over `F_4`, that
normalization leaves exactly 4,096 boundary choices up to a paired phase
equivalence. Exact four-jet certificates exclude every one of them for
the **two-position bare criterion**, with all-scale product at least 738.

The endpoint and boundary normalization are uniform ordinary proofs.
The 4,096-case exclusion is finite. A separate complete coefficient-span
test rejects arbitrary F4 recombinations at two selected boundaries,
including the stabilizer boundary, without relying on Smith profiles.
Longer ordered words, new extension-field boundary coefficients, and
physical restrictions remain outside these tests. No stronger factoring
exponent or new Lean work is claimed.

## 1. A literal ordered quadratic phase

Use physical letters `0,1,2` and coefficient field
`F_4=F_2[zeta]/(zeta^2+zeta+1)`. Let

```math
X e_j=e_{j+1},\qquad Y e_j=\zeta^j e_j,\qquad YX=\zeta XY,
```

with indices modulo three. The nine matrices `W_(a,b)=X^aY^b` form a
basis of the full three-square matrix algebra. Their multiplication is

```math
(\zeta^cX^aY^b)(\zeta^{c'}X^{a'}Y^{b'})
 =\zeta^{c+c'+ba'}X^{a+a'}Y^{b+b'}.                       \tag{1}
```

All exponent arithmetic below is modulo three, including division by two;
thus `1/2=2` in exponents. Coefficient arithmetic is in characteristic two.

Let `T(t)` be the [three-letter Steiner cubic](steiner-cubic-crossings.md).
Choose any linear functional `ell_w` with
`ell_w(W_(a,b))=w_(a,b)` and `w_(0,0)=1`. Define

```math
H_m(t)=\ell_w\left(\prod_{i=1}^m
 T(t;x_i,y_i,u_i,v_i)X^{x_i}Y^{y_i-u_i}X^{-v_i}\right),   \tag{2}
```

retaining the physical order. At the crossing endpoint each local
auxiliary matrix is the identity, so `H_m(1)=E_m` exactly.

Put `a_i=x_i-v_i`, `b_i=y_i-u_i`, `A=sum_i a_i`, `B=sum_i b_i`.
Repeated application of (1) gives the complete entry formula

```math
H_m(t)=\left(\prod_i T(t;x_i,y_i,u_i,v_i)\right)
       \zeta^{Q}\,w_{A,B},\qquad
Q=-\sum_i v_i b_i+\sum_{i<j}b_i a_j.                    \tag{3}
```

The order-dependent cross term in (3) is retained. This is a special
[ordered physical gauge](ordered-physical-gauge-crossings.md) with a
general linear auxiliary closure. The closure is effective: since
`Tr(W_(a,b)^-1 W_(c,d))` is one for equal indices and zero otherwise,

```math
\ell_w(M)=\operatorname{Tr}\left(
 \left(\sum_{a,b}w_{a,b}W_{a,b}^{-1}\right)M\right).       \tag{4}
```

Left multiplication on the nine-dimensional matrix space realizes this
as a fixed boundary `f^T product(K_i)e`, with `f^T e=1`.
The degree is at most `3m`, with at most `3m+1` coefficient positions.
Fixed field, tables and auxiliary dimension retain `2^O(m)` complete
construction/factorization cost. The existing
[fixed-field descent](fixed-field-descent.md) and
[scale-one modular transfer](scale-one-modular-transfer.md) supply the
factoring interface for a family with favorable uniform rank bounds.
Those bounds remain missing here.

## 2. Classifying rank-one zero cuts at every length

At zero, the scalar tensor is the first-cut identity: `u_i=x_i,v_i=y_i`.
Then `B=-A` and the phase reduces to

```math
Q=\frac12\sum_i(x_i^2-y_i^2)-\frac12 A^2.                \tag{5}
```

After separated nonzero row/column phase factors, the second-cut kernel
is therefore

```math
h\left(\sum_i x_i-\sum_i y_i\right),\qquad
h(a)=\zeta^{-a^2/2}w_{a,-a}.                             \tag{6}
```

Only the diagonal physical pairs support it. The sum map onto `F_3` is
surjective at every positive length, so its rank is exactly that of the
three-square circulant with entries `h(s-t)`. The Fourier transform over
the third roots diagonalizes this circulant. Since `h(0)=1`, its rank is
one precisely when

```math
h(a)=\zeta^{ca}\quad\text{for one }c\in\mathbb F_3.       \tag{7}
```

Multiplying all weights by `zeta^(-c a)` multiplies (3) by
`product_i zeta^(-c x_i) zeta^(c v_i)`. This is a paired, constant
physical-leg change: it preserves the crossing and both complete cut
profiles/ranks. Thus every rank-one cold boundary is equivalent to

```math
w_{0,0}=1,\qquad w_{1,2}=w_{2,1}=\zeta^2.                \tag{8}
```

The other six entries are free. Over `F_4` this gives `4^6=4096` canonical
boundaries, representing all rank-one cold boundaries of this particular
linear clock/shift construction up to the stated equivalence. The proof
also works over larger fields containing `zeta`; new coefficients in
those fields are not part of the finite catalogue below.

For the canonical boundary the cold tensor has the separated phase
`zeta^(sum_i(x_i^2-y_i^2)/2)` on its support. In particular its second
cut has rank one for every length, regardless of the six free weights.

## 3. Exact four-jet exclusion at two positions

At length two the physical alphabet is nine, the cut matrix order is
81, and the polynomial degree is at most six. Let `r` be the first-cut
rank of the degree-six coefficient. It counts the raw infinity orders
equal to `-6`. Let `R_j` be the second-cut block-Toeplitz rank at precision
`j`. The number of its zero-place Smith orders at most three is

```math
N_4=R_4-R_3.                                            \tag{9}
```

For every real `s>=1`, the bare partitions therefore satisfy

```math
U(s)\ge r s^6,\qquad V(s)\ge N_4s^{-3},\qquad
U(s)V(s)\ge rN_4s^3\ge rN_4.                            \tag{10}
```

The complete catalogue has

```math
\min_{4096\text{ boundaries}}r(R_4-R_3)=738>9^3=729.      \tag{11}
```

Consequently every displayed boundary fails the bare target, with
`rho>=log(738)/(4 log(9))>3/4`. The comparison uses only integers and
three ranks per boundary. It does not rely on optimizing a real scale,
on a numerical LP, or on a truncated zero being treated as a missing pivot.
For example, canonical code 1694 has `r=18,R_3=27,R_4=68`, giving 738.

This is not an all-length obstruction: (2) can couple auxiliaries across
longer words, so the two-position calculation does not automatically
tensor into a lower bound for that ordered continuation. It also does
not exclude a different criterion, physical restrictions, or boundary
entries newly available over a field extension.

The later [box-budget theorem](bare-seven-to-twelve-letters.md), contributed
in [PR #55](https://github.com/1337isnot1337/integer-factoring/pull/55),
excludes the bare half target for every crossing on at most twelve letters,
at every degree. It therefore subsumes bare failure for these nine-letter
seeds. The construction, boundary classification and explicit early-jet
ranks above remain additional information. The next section's actual
coefficient-span tests are outside that bare-profile theorem.

## 4. Independent replay

The [constructor](../tools/weyl_ordered_crossing.py) builds the exact
coefficient cuts from (1)--(3). The complete boundary variation is affine
over twelve binary coordinates. It checks that those free directions
leave the cold matrix unchanged, then enumerates all 4,096 choices.

The [checker](../tools/check_weyl_ordered_boundaries.py) restricts scalars
from `F_4` to two-dimensional `F_2` matrices and uses binary elimination;
each resulting rank is exactly twice the original field rank. It replays
all **12,288** leading/Toeplitz rank entries with the opposite pivot order.
There are also 384 direct dense three-square Weyl-product and closure
controls at lengths one through six, including the separated cold phase.
[The evidence](../evidence/weyl-ordered-boundaries.json) retains every
boundary code, three ranks, the all-scale floor and source hashes.

```sh
make verify-weyl-ordered-boundaries
python3 tools/check_weyl_ordered_boundaries.py --generate
make verify
```

The replay and generation use the Python standard library. The complete
factoring score remains `1/2`; no new Lean work was used. The next
construction needs cancellation beyond this normalized two-position
boundary family, rather than another choice already covered by (11).

## 5. Actual middle-coefficient barriers at two boundaries

Write the complete two-position polynomial as

```math
H_2(t)=\sum_{j=0}^6t^j A_j,\qquad \sum_j A_j=E_2.
```

For the two canonical boundary codes below, enumerate **every** functional
`w` in `F4^7` with `w_3=1`, and form `A(w)=sum_j w_j A_j`. There are
`4^6=4096` per boundary. Exact ranks give

| Boundary code | Minimum `r(A(w)) c(A(w))` with `w_3=1` | Attaining coefficients `(w_0,...,w_6)` |
| --- | --- | --- |
| 0 | 576 | `(0,1,0,1,0,1,0)` with ranks `(36,16)` |
| 1694 | 1950 | `(0,0,0,1,0,0,0)` with ranks `(50,39)` |

Code zero is the rank-one stabilizer closure. With `zeta=2` in the
constructor's F4 encoding, it is `ell(M)=phi^T M psi`, where
`psi=(1,zeta,zeta)` and `phi=(1,zeta^2,zeta^2)`; `phi^T psi=1`.
Code 1694 attains the best four-jet floor in (11). These are **two selected
boundaries**, not all 4,096 cold-boundary classes.

Here is why the table obstructs the complete coefficient span. Scaling
by a nonzero field element does not change either rank, so the same
positive floor holds whenever `w_3!=0`. In particular, every relation
`sum_j w_j A_j=0` has `w_3=0`. Thus sending `A(w)` to `w_3` defines a
well-defined linear functional on the actual tensor span. It sends the
crossing to one. Every decomposition of the crossing into arbitrary
members of that span must therefore contain a term with nonzero middle
functional, and that term has at least the displayed rank product.
The number of pieces is unrestricted in this argument.

For the nine-letter block alphabet the finite half target requires
`r c < 9^(2 sqrt(2))`. Already 576 exceeds this threshold: the integer
comparisons `576^7>9^20` and `2*5^2>7^2` give the strict normalized floor
`rho>5/7>1/sqrt(2)`. The exact normalized lower floors are approximately
0.7231973 and 0.8619493. These are finite-block obstructions; the ordered
continuation is not their ordinary tensor power.

The [middle-coefficient checker](../tools/check_weyl_middle_coefficients.py)
generates the full F4 family by Gray-code additions and replays all
16,384 cut ranks with the opposite binary pivot order after exact scalar
restriction. Its [evidence](../evidence/weyl-middle-coefficients.json)
records every normalized direction's pair of ranks and source hashes.

```sh
make verify-weyl-middle-coefficients
python3 tools/check_weyl_middle_coefficients.py --generate
```

No floating-point decision, rank-sum shortest-path solver, or extrapolation
from monomial coefficients is used for this test. Other cold boundaries,
larger extension fields, and longer ordered words remain open for actual
coefficient recombinations.
