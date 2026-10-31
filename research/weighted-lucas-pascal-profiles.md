# Weighted Pascal profiles without polynomial elimination

## Result and scope

The Pascal-complement construction admits independent mode heights whose
profiles can be computed using only constant matrices over a prime field.
For a Lucas antichain, every height vector is allowed. Equal digit-sum
layers are concrete antichains, so the result supplies a direct profile
interface for sparse digital mode families.

The complement search also extends the binary gap-skipping traversal to
every fixed prime characteristic. It visits polynomially many radix
blocks in the bit length of the mode indices. It need not enumerate the
possibly enormous gaps between modes.

These are constructive exact profile results. They do not supply a
strict pressure witness below `1/sqrt(2)`. The original crossing degree
and the eventual coefficient count remain charged; a fast profile
calculation does not make a growing, exponentially large mode index free.

The unweighted factorization, all-field crossing, and original binary
traversal are from predecessor
`18cec06deac47c8f142ea67a60e2278ac2234d67`, in
[the all-field construction](all-field-pascal-crossings.md),
[the complement note](schur-pascal-complement.md), and
[the convolution bounds](pascal-complement-convolution.md).
The weighted reduction and prime-radix traversal below extend those
results. The work used substantial OpenAI assistant support.

## 1. An exact unit factorization of the complement block

Let `K` be a field, `S={a_1<...<a_n}` a finite mode set, and
`C={c_1<...<c_n}` its greedy independent Pascal columns outside `S`.
Write

~~~math
U_C=\left(\binom{c_j}{a_i}\right)_{i,j=1}^n,
\qquad D_C=\mathrm{diag}(t^{c_1},\ldots,t^{c_n}).
~~~

The matrix `U_C` is invertible. The complement block in the all-field
construction has the formal factorization

~~~math
K_S=Q_+\sum_{c\notin S}u_c t^c V_c(t),
\qquad
Q_+=\mathrm{diag}((1+t)^{a+1}:a\in S),
\tag{1}
~~~

where `u_c=(binom(c,a))_(a in S)` and

~~~math
V_c(0)=(-1)^c u_c^{\mathsf T}
       \mathrm{diag}((-1)^a:a\in S).
~~~

**Unit factorization.** There is an invertible power-series matrix
`W_S(t)` such that

~~~math
\boxed{K_S=Q_+ U_C D_C W_S(t).}
\tag{2}
~~~

To prove this, express every column uniquely as
`u_c=sum_j beta_(j,c) u_(c_j)`. The greedy prefix-span property gives
`beta_(j,c)=0` whenever `c<c_j`: the earlier selected columns already
span every earlier column. Therefore

~~~math
(W_S)_{j,*}=
\sum_{c\notin S}\beta_{j,c}t^{c-c_j}V_c(t)
\tag{3}
~~~

is a power series. Each fixed coefficient receives only finitely many
terms. Its constant term comes only from `c=c_j`, since the other greedy
columns have their own basis coordinates. Hence

~~~math
W_S(0)=\mathrm{diag}((-1)^{c_j})
       U_C^{\mathsf T}\mathrm{diag}((-1)^a),
\tag{4}
~~~

which is invertible. This proves (2).

The factorization is stronger than a list of unweighted Smith orders:
the left non-diagonal factor is the explicit constant matrix `U_C`.
The diagonal unit `Q_+` commutes with arbitrary diagonal mode heights.

## 2. Weighted reduction of the full infinity matrix

Use the normalized matrix

~~~math
F_S(t)=
\begin{pmatrix}
QP^{\mathsf T}D&DR\\
D&(1+t)O
\end{pmatrix},
\tag{5}
~~~

from the all-field note. Assign one rational height `h_a` to mode `a`
on each of its two row blocks, and put `H=diag(t^(h_a))`.
Rational powers mean a fixed common ramification.

Assume

~~~math
h_a\ge h_b
\quad\text{whenever }\binom ba\ne0\text{ in }K.
\tag{6}
~~~

Then regular row and column operations give

~~~math
\boxed{
\mathrm{diag}(H,H)F_S
\ \sim\
\mathrm{diag}\bigl(HD,\ H U_C D_C\bigr).
}
\tag{7}
~~~

Indeed, after interchanging the row blocks, the operation clearing the
lower-left block has coefficient `H Q P^T H^(-1)`. Its entries are
regular under (6). The usual column
clearing operation is still regular because row `a` of `O` is divisible
by `t^a`. The result is `diag(HD,HK_S)`. Apply (2), commute the
diagonal unit `Q_+` with `H`, and remove the right unit `W_S` to obtain
(7).

Consequently the order list is the multiset union of

~~~math
\{a+h_a:a\in S\}
\quad\text{and the orders of}\quad
\mathrm{diag}(t^{h_a})U_C\mathrm{diag}(t^{c_j}).
\tag{8}
~~~

Condition (6) matters. In characteristic two, take `S={0,1}` and
`h=(-2,-1)`. Applying (8) despite the failed condition predicts
`(-2,0,0,2)`; direct polynomial elimination gives `(-2,-1,1,2)`.
The implementation rejects this height vector.

### Antichains permit every height

In characteristic `p`, write each index in base `p`. The identity

~~~math
(1+X)^c=\prod_j(1+X^{p^j})^{c_j}
\quad\text{in }\mathbb F_p[X]
\tag{9}
~~~

gives Lucas's digit formula
`binom(c,a)=prod_j binom(c_j,a_j) mod p`. Thus `binom(b,a)` is
nonzero exactly when every digit of `a` is at most the corresponding
digit of `b`.

Call `S` a Lucas antichain when no two distinct modes are comparable in
this digit order. Then `P_S=I`, and (6) imposes no condition on distinct
heights. In particular, any fixed digit-sum layer is an antichain.
This includes sparse families with large gaps, and is different from
appending a common string of trailing binary ones, which the existing
[dyadic reduction](one-species-dyadic-affine-reduction.md) removes.

## 3. Weighted constant elimination computes every order

Let `A` be an invertible constant matrix and consider

~~~math
M_{ij}=A_{ij}t^{r_i+c_j}.
~~~

Choose a nonzero entry minimizing `r_i+c_j`. Every entry in its column
has at least that valuation, and every entry in its row does too.
The row and column operations clearing that pivot are therefore regular.
The Schur complement retains the same monomial form:

~~~math
A'_{k\ell}t^{r_k+c_\ell},\qquad
A'_{k\ell}=A_{k\ell}-A_{kj}A_{ij}^{-1}A_{i\ell}.
\tag{10}
~~~

Repeat on the smaller constant matrix. The pivot valuations, in order,
are the exact Smith orders. Minimum pivoting makes them nondecreasing;
ties may be resolved arbitrarily. All field operations are on constants.
No power series need be stored, no precision cutoff is selected, and
no determinant subset enumeration is used.

The routine requires `O(n^3)` field operations, plus comparisons and
additions of the height numbers. Integer and rational heights give an
exact computation. Floating heights can guide a search but require
exact rational replay before being used in a pressure certificate.

For physical infinity heights `eta_a` in the original variable `xi`,
put `h_a=eta_a/2` in (8), since the normalized variable is `t=xi^2`.
Multiply every resulting order by two and add the original common
shift `1-2n`. The restored physical mode energy increases by `2eta_a`.
The zero profile remains unchanged. Thus (8) plugs directly into the
existing product-pressure formula with independent mode heights.

## 4. Gap-skipping traversal in every prime characteristic

Fix a prime `p`. For an aligned radix block

~~~math
\mathcal I=\{B,B+1,\ldots,B+p^k-1\},\qquad p^k\mid B,
~~~

write each row index `a=a_hi p^k+a_lo`. Formula (9) gives, for
`c=B+d` in that block,

~~~math
\binom ca=
\binom{B/p^k}{a_{\rm hi}}
\binom d{a_{\rm lo}}\pmod p.
\tag{11}
~~~

For each distinct low index `ell`, form the vector whose entry at row
`a` is the first factor in (11) if `a_lo=ell`, and zero otherwise.
The span of these group vectors is exactly the span of all Pascal
columns in the block. The low-index Pascal functions are independent,
because they are distinct rows of an invertible full Pascal matrix.

Maintain the span of already accepted columns. If every group vector
is in that span, discard the entire block. Otherwise split it into its
`p` children and visit them in increasing order. At a leaf, skip indices
in `S` and accept a column precisely when it is independent. Stop after
`n` accepted columns.

This reproduces the ordinary greedy scan exactly. Every expanded path
ends at an accepted column or an excluded leaf that prevented pruning.
There are at most `n` accepted and `n` excluded leaves. For fixed `p`,
the number of visited blocks is `O(p n log_p(max S+2))`; the elementary
linear algebra and Lucas coefficient calculations give a polynomial
bit cost in `n` and the bit length of the mode indices. This is an
algorithm for the finite profile data, not a claim that all original
tensor coefficients can be enumerated in that same cost.

## 5. Verification and search interface

The reusable module is
[lucas_weighted_pascal.py](../tools/lucas_weighted_pascal.py). Run the
independent checker with

~~~sh
python tools/check_lucas_weighted_pascal.py \
  --output evidence/lucas-weighted-pascal.json
~~~

The recorded checks include:

- 13,280 comparisons of Lucas coefficients with ordinary integer
  binomial coefficients in characteristics 2, 3, 5, and 7.
- 252 complement comparisons with a separate unskipped modular scan.
- Twelve large-gap cases, with high radix positions 32, 64, and 128.
- 96 weighted profile comparisons against full polynomial Smith
  elimination, including unrestricted antichain heights and monotone
  heights on sets with comparable modes.
- The failed-regularity counterexample and an exact rational
  ramification control.

These checks support the implementation. Equations (2), (7), (10), and
(11) are the all-size constructive proofs. Searching these profiles
still requires a strict rational pressure inequality and the effective
crossing transfer before the factoring score changes.
