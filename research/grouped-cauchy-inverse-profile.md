# The grouped inverse has the same binomial precision profile

## Status

The normalized inverse in
[the Boolean Smith note](boolean-cauchy-smith-profile.md) has a complete
Smith profile across the two source groups as well as across its input
and output sides. At an odd local evaluation point the two profiles
coincide at every retained center. At either outer retained center they
coincide at every integral point, including even points and zero.

Consequently, the small leading rank at an outer center does not give
a sparse inverse decomposition at the precision needed to undo its
normalization. For a tensor product of packets, every Boolean direction
survives at that precision. This closes a specific direct-inverse
replacement of the contraction step. It does not rule out other
deformations, scalar decoders, or decompositions adapted to particular
input tensors.

## 1. Setup and result

Use the notation of the Boolean Smith note. In particular,

\[
 L=2^a,\qquad 0<|t|<L,\qquad b=|t|,\qquad
 B=\{x:x+t\notin[0,L)\},\qquad \sigma(x)=(x+t)\bmod L.
\]

Set `i=v_2(t)`. Let `G_0,G_1` partition the packet positions by bit
`i`. The permutation `sigma` interchanges these two sets. Both the
input and output ports of the normalized inverse have the source's
column orientation, so its grouping cut places both occurrences of a
position in its respective `G_g`.

For an integral local value `zeta`, the monomial-basis normalized inverse
is

\[
 R^{\mathrm{mon}}
   =U^{-1}\mathrm{diag}_I(2\zeta)^{b-|I\cap B|},\qquad
 U=S(W'(\zeta)).
 \tag{1}
\]

Here `W'` divides the physical boundary rows by `2z`, and
`W'(zeta) mod 2` is the cyclic permutation. The value-basis inverse is
obtained by local two by two zeta changes. Those changes preserve
Smith orders across the grouping cut.

**Theorem.**

* If `zeta` is a unit, the grouped matrix of (1) has orders
  `j`, with multiplicity `binom(b,j) 2^(L-b)`, for `0<=j<=b`.
* If `b=L-1` and `zeta!=0`, put `h=v_2(zeta)`. The grouped orders are
  `j(1+h)`, with multiplicity `2 binom(L-1,j)`.
* If `b=L-1` and `zeta=0`, the grouped matrix has exactly two unit
  Smith factors and all its other factors are zero.

The first statement holds more generally when `W'` is any physical
matrix congruent to the same permutation modulo two. The outer-center
statement needs only the degree-triangular structure of its Boolean
substitution inverse and its congruence to that permutation.

## 2. A weighted Boolean substitution lemma

Let `A` be an `L` by `L` matrix with `A mod 2=I`. For a coordinate
subset `C`, put `f_C(I)=|I intersect C|`. Then

\[
 v_2\bigl(S(A)_{I,J}\bigr)
      \ge |f_C(I)-f_C(J)|.
 \tag{2}
\]

The valuation of zero is infinite. To prove (2), expand a product of
the row forms indexed by `I`. Each term selects one column from each
row. If it uses `k` off-diagonal entries, it has valuation at least
`k`. Relative to the identity choices, those `k` changed choices can
change the number of distinct selected coordinates in `C` by at most
`k`. This remains true when repeated columns undergo Boolean
reduction. Thus every term contributing to `J` has the required
valuation, and summing terms preserves the bound.

Write `Q_C=diag_I 2^(f_C(I))`. Both conjugates

\[
 Q_C S(A)Q_C^{-1},\qquad Q_C^{-1}S(A)Q_C
 \tag{3}
\]

are integral. Their determinants equal the unit `det S(A)`, so they
are invertible over the same local ring. Taking inverses in (3)
proves (2) with `S(A)^(-1)` in place of `S(A)`.

For a matrix `A` congruent to the permutation `sigma`, normalize its
columns so that `A_0 mod 2=I`. Then
`S(A)=S(A_0)P_sigma`. Consequently

\[
 v_2\bigl(S(A)^{-1}_{J,I}\bigr)
   \ge\left||\sigma^{-1}(J)\cap C|-|I\cap C|\right|.
 \tag{4}
\]

This estimate holds for every choice of `C`, independently of `L`.
It is an entrywise statement about the inverse. It does not assume
that Boolean substitution commutes with matrix inversion.

## 3. Flattening and the unit-point proof

Write `J=J_0 union J_1` and `I=I_0 union I_1` using the two groups.
Flatten the inverse with row index
`r=(J_0,I_0)` and column index `c=(J_1,I_1)`.
Let `H` be the corresponding flattening of `U^(-1)`. Its reduction
modulo two is a permutation matrix `Pi`. The permutation matches

\[
 J_1=\sigma(I_0),\qquad I_1=\sigma^{-1}(J_0).
 \tag{5}
\]

Set `E=H Pi^(-1)`, so `E mod 2=I`. Let

\[
 d_0(r)=|B\cap G_0|-|I_0\cap B|,\qquad
 d_1(r)=|B\cap G_1|-|\sigma^{-1}(J_0)\cap B|.
 \tag{6}
\]

After multiplying the grouped matrix by `Pi^(-1)` on its right,
it is exactly

\[
 D_0 E D_1,\qquad
 D_g=\mathrm{diag}_r(2\zeta)^{d_g(r)}.
 \tag{7}
\]

Take `C=B intersect G_0` in (4). For the entry `E_(r,s)`, the
restriction of `sigma^(-1)(J)` to `G_0` is `I_(0,s)`, while the
restriction of `I` to `G_0` is `I_(0,r)`. It follows that

\[
 v_2(E_{r,s})\ge |d_0(r)-d_0(s)|.
 \tag{8}
\]

If `zeta` is a unit, (8) makes `D_0 E D_0^(-1)` integral.
Its determinant equals `det E`, an odd unit. Therefore (7) is
unimodularly equivalent to `D_0D_1`. Its diagonal exponents
`d_0(r)+d_1(r)` count the absent bits in a set of `b` independent
Boolean coordinates. This proves the first part of the theorem.

The argument uses that `2zeta` has valuation one. Without more
structure, (8) does not justify replacing that valuation by an
arbitrarily larger one.

## 4. The outer-center theorem at every local point

Suppose `b=L-1`. There is exactly one nonboundary position, denoted
`x_*`. Rename the two groups so that `x_*` belongs to `G_0`.
Then `B intersect G_0=G_0 minus {x_*}` and `B intersect G_1=G_1`.
Write

\[
 \epsilon(r)=\mathbf1[x_*\in I_0],\qquad
 q(r)=|J_0|-|I_0|.
 \tag{9}
\]

Boolean substitution does not increase degree. With monomials ordered
by their degree, `U` is block triangular. Each diagonal block is
invertible modulo two, so its inverse has the same triangular pattern.
After flattening and the permutation (5), this says

\[
 E_{r,s}=0\quad\text{unless}\quad q(r)\ge q(s).
 \tag{10}
\]

Indeed the input degree in that inverse entry is
`|I_(0,s)|+|J_(0,r)|`, and its output degree is
`|I_(0,r)|+|J_(0,s)|`.

At an outer center, (6) gives

\[
 d_0(r)-d_1(r)=
      (|B\cap G_0|-|B\cap G_1|)+q(r)+\epsilon(r).
 \tag{11}
\]

Order the flattened indices by increasing `d_0`, breaking ties by
increasing `d_1`. A nonzero entry `E_(r,s)` above the diagonal must
satisfy `d_1(r)<=d_1(s)`. To check this, the only problematic case
would have `d_0(r)<d_0(s)` and `d_1(r)>d_1(s)`. Their differences
are integers, so (11), together with
`epsilon(r)-epsilon(s) in {-1,0,1}`, would force `q(r)<q(s)`,
contradicting (10).

Every leading principal minor of `E` is odd because `E mod 2=I`.
Thus Gaussian elimination without row exchanges gives an integral
unit decomposition

\[
 E=L_0 A_0 U_0,
 \tag{12}
\]

where `L_0` and `U_0` are unit lower and upper triangular, and `A_0`
is diagonal with unit entries. The zero pattern (10) is preserved:
an elimination update can be nonzero only along
`q(r)>=q(k)>=q(s)`. Hence every nonzero upper entry of `U_0`
also has `d_1(r)<=d_1(s)`.

For any integral scalar `delta`, define diagonal matrices using
`delta^(d_g)` in place of `(2zeta)^(d_g)`. The two matrices

\[
 D_0 L_0 D_0^{-1},\qquad D_1^{-1}U_0D_1
 \tag{13}
\]

have only nonnegative powers of `delta` in their entries. The first
uses the increasing order of `d_0`; the second uses the property of
the upper entries just proved. They are unit triangular, including
when `delta=0`. We can therefore write (7) as an invertible
equivalence to

\[
 D_0 A_0 D_1
       =\mathrm{diag}_r\bigl(a_r\delta^{d_0(r)+d_1(r)}\bigr).
 \tag{14}
\]

Set `delta=2zeta`. The units `a_r` do not affect Smith orders, and
the exponent distribution is the same binomial distribution as
before. If `zeta=0`, only the two exponent-zero diagonal entries
survive. This proves the two remaining claims.

This argument applies separately at each integral `zeta`; the
coefficient matrix `E` is allowed to depend on that point.
The expression (13) at zero means the displayed nonnegative-power
form, with no division by zero.

## 5. Consequence for direct inverse decompositions

Consider `P` full Boolean packets whose two symbol groups lie on opposite sides
of a fixed grouping cut. Take unit local points and let their boundary
sizes be `b_1,...,b_P`. Put

\[
 M_{\mathrm{slots}}=LP,\qquad B_{\mathrm{tot}}=\sum_{k=1}^P b_k.
\]

For simplicity all packet lengths are displayed as `L`; unequal
lengths give the same argument with the sum of their lengths.
The tensor product of their grouped normalized inverses has Smith
orders

\[
 j\quad\text{with multiplicity}\quad
 \binom{B_{\mathrm{tot}}}{j}\,
                  2^{M_{\mathrm{slots}}-B_{\mathrm{tot}}},
 \qquad 0\le j\le B_{\mathrm{tot}}.
 \tag{15}
\]

This follows by taking tensor products of the exact unimodular
factorizations above. The scalar normalization used to turn these
inverses back into identities has valuation `B_tot`. Recovering even
one output bit therefore requires the normalized contraction modulo
`2^(B_tot+1)`. At that precision every one of the
`2^(M_slots)` Smith factors in (15) is nonzero.

Over `O/(2^(B_tot+1))`, the image of the diagonal matrix needs that
many generators. Each nonzero cyclic summand contributes one
dimension to `image/(2 image)` over the residue field. A sum of `r`
outer products factors through a free module of rank `r`. Over this
principal ideal ring, a submodule of that free module needs at most
`r` generators; applying the other factor cannot increase that number.
Consequently an exact sum of products of one tensor per side of this
cut requires at least

\[
 2^{M_{\mathrm{slots}}}
 \tag{16}
\]

terms at the needed precision. Using even local points at outer
centers gives the same conclusion after increasing the normalization
loss by their valuations.

This statement is about an explicit grouped factorization of the
normalized inverse itself. It does not imply that every method of
contracting that inverse with special surrounding tensors pays (16).
Fixing unused packet slots can also reduce the relevant operator.
It also does not replace the source's original mode analysis, which
decomposes the forward deformation without introducing this inverse.

## 6. Exact finite checks

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 tools/check_grouped_cauchy_inverse.py \
  --output /tmp/grouped-cauchy-inverse-profile.json
```

The checker includes randomly generated matrices `I+2A` to test both
directions of (2) and the inverse version. Its deterministic seeds
make these controls reproducible. It separately constructs grouped
inverse matrices and applies minimum-valuation elimination.

Cauchy cases include zero, odd, and even local points at outer centers,
and unit points at inner centers. Additional physical matrices
congruent to the cyclic permutation test that the unit-point theorem
does not depend on a special Cauchy identity. Small tensor products
test (15) and confirm that all Smith factors survive after the scalar
normalization loss is paid.

The counts, profiles, and source hashes are recorded in
[grouped-cauchy-inverse-profile.json](../evidence/grouped-cauchy-inverse-profile.json).
