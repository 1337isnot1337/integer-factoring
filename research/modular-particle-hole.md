# A reflected particle-hole identity for the complete modular profile

## Result and scope

There is an exact compression rule for characteristic-two labels occupying
more than half of a complete binary block of modes. It removes a common
scalar polynomial from every four-label crossing profile and replaces the
labels by explicit reflected complements. This is an all-size algebraic
identity, not an extrapolation from small ranks.

Let `K=2^a` and `F_K={1,3,...,2K-1}`. Suppose every label `mu` is a subset
of `F_K` with the same cardinality `r` and the same weight `h`. Put

```math
C_K(mu)=F_K\setminus\{2K-t:t\in mu\},\qquad
h'=K^2-2Kr+h,\qquad D_K(z)=1+z^{2K}+z^{4K}.
\tag{1}
```

The new labels have cardinality `K-r` and weight `h'`. Write
`P_mu(z)=z^(2h) sum_s A_s(x,y,u,v)z^s` for the complete normalized
profile from the folded Pfaffian formula. Then, applying `C_K` to all
four labels,

```math
P_{C_K(mu)}(z)=D_K(z)^{K-2r}P_mu(z).
\tag{2}
```

If `r>K/2`, equation (2) means an exact polynomial quotient:
`D_K^(2r-K)` divides every profile, and the quotient is precisely the
reflected-complement profile. Since `D_K(1)=1`, removing the factor
preserves the exact scale-one crossing sum. No coefficient-algebra value
is inverted.

This supplies a new effective transformation and explains some finite
rank changes under overlapping filled-core convolution. It gives no
uniform rank-width coefficient below `1/sqrt(2)` and no factoring
exponent above `1/2`.

The crossing and normal-order operators are Claude (Anthropic)'s Appendix B
at [formal-math@16a5ce4](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The folded polynomial matrix is from
[PR #22 at 735b3dd](https://github.com/1337isnot1337/integer-factoring/blob/735b3ddb37d2e63dc469cf833687355e42097dd4/research/modular-pfaffian.md).
The complete-block notation below also appears in its later
[filled-core Schur proof at c12da89](https://github.com/1337isnot1337/integer-factoring/blob/c12da897fcc6081685e8c663217ce9083db87b4a/research/fermi-sea-schur.md).
The inverse-conjugation identity and its complementary-minor consequence
below are a separate ordinary argument, developed with substantial OpenAI
assistant support. No Lean formalization is included.

## 1. The full polynomial matrix

Use the folded matrix on four copies of the full core `F_K`, in block order
`U,A,B,Y`, and call it `M_K`. Put `w=z^2`, `q=1+w`, `s=1+z`, and define
`K` by `K` matrices, for `0<=i,j<K`, by

```math
W=diag(w^i),\quad V=diag(q^i),\quad
P_{ij}=[i\subseteq j],\quad R_{ij}=[i+j=K-1],
```

```math
J_{ij}=[i\mathbin{\&}j=0],\quad
H_{ij}=w^{i\mathbin{\&}j},\quad O_{ij}=w^{i\mathbin{|}j}.
```

The subset relation refers to binary one-bits. The six upper blocks of
`M_K` are

```math
UA=zsO,\quad UB=sP^TV,\quad UY=zW,\quad
AB=zWP,\quad AY=sH,\quad BY=zsVJW.
\tag{3}
```

Opposite blocks are transposes; diagonal blocks vanish. These are direct
substitutions of `2i+1` in the folded formula. Let
`Rbar=diag(R,R,R,R)`.

**Inverse identity.** Over `F_2(z)`,

```math
M_K^{-1}=D_K^{-1} Rbar M_K Rbar.
\tag{4}
```

In particular, `M_K` is invertible as a rational-function matrix.

## 2. Prove the inverse identity by binary tensor products

All seven mode matrices are tensor products over the `a` binary digits.
For the digit with `b=w^(2^j)`, the factors are

```math
P_b=\begin{pmatrix}1&1\\0&1\end{pmatrix},\quad
R_b=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
W_b=diag(1,b),\quad V_b=diag(1,1+b),
```

```math
J_b=\begin{pmatrix}1&1\\1&0\end{pmatrix},\quad
H_b=\begin{pmatrix}1&1\\1&b\end{pmatrix},\quad
O_b=\begin{pmatrix}1&b\\b&b\end{pmatrix}.
```

The empty tensor product covers `K=1`. Two-by-two multiplication and then
tensoring prove the following six identities:

```math
\begin{aligned}
P^TVRP^TW&=WRH, & ORWP&=WRWJ^TV,\\
ORH&=P^TVRVJW, & ORP^TV&=HRWJ^TV,\\
ORW&=WPRVJW, & VPRW&=P^TWRH.
\end{aligned}
\tag{5}
```

Each side has the same digit tensor factors, with the same matrix order.
These equations make all off-diagonal blocks of `(M_K Rbar)^2` vanish:
the two terms in each off-diagonal block are equal and cancel in `F_2`.
They correspond respectively to block positions `UA,UB,UY,AB,AY,BY`.

For the diagonal blocks, the three contributions are, in an appropriate
order, `TQ I`, `Q I`, and `T I`, where `T=w^K`, `Q=q^K`.
For example, the `U,U` block is

```math
(zs)^2 OROR+s^2 P^TVRVP R+z^2 WRWR
=(TQ+Q+T)I.
```

The digit identities `OROR=(wq)^(K-1) I`, `HRHR=q^(K-1) I`,
`PRP^T=R`, `WRWR=w^(K-1) I`, and `VRVR=q^(K-1) I`
also verify the other three diagonal blocks. Since `K` is a power of two,
`Q=(1+w)^K=1+T`. Consequently

```math
(M_K Rbar)^2=(1+T+T^2)I=D_K I,
```

which proves (4). All calculations before inversion are polynomial
identities over `F_2[z]`.

## 3. The full Pfaffian and complementary minors

Taking determinants in the last identity gives
`(det M_K)^2=D_K^(4K)`. Squaring is injective in `F_2[z]`, so
`det M_K=D_K^(2K)`. Since `Pf(M_K)^2=det M_K`, another use of that
injectivity gives

```math
Pf(M_K)=D_K^K.
\tag{6}
```

For any invertible even alternating matrix `M` in characteristic two and
any even index set `I`, the complementary Pfaffian identity is

```math
Pf(M[I^c])=Pf(M) Pf(M^{-1}[I]).
\tag{7}
```

One proof uses a block Schur complement when the complementary principal
block is invertible. The inverse block on `I` is the inverse of that Schur
complement, and its Pfaffian is the reciprocal of the Schur complement's
Pfaffian. Clearing denominators gives a polynomial identity for a generic
alternating matrix; specialization proves (7) even when that principal
block is singular. Matching signs disappear in characteristic two.

Now let `S` be the `4r` vertices selected by the four labels in `M_K`.
Reflection sends their index set to `Rbar S`. The reflected-complement
labels select `(Rbar S)^c`. Apply (7), then (4):

```math
\begin{aligned}
Pf(M_K[(Rbar S)^c])
&=D_K^K Pf(M_K^{-1}[Rbar S])\\
&=D_K^K D_K^{-2r} Pf(M_K[S]).
\end{aligned}
```

The last exponent is `2r` because a Pfaffian on `4r` vertices has `2r`
factors in each matching. The two principal Pfaffians are the two complete
profiles, proving (2). This also proves the exact divisibility when
`K-2r` is negative. The label map is an involution.

## 4. Effective compression and what it does not prove

For `K=poly(lambda)`, reflecting and complementing a label is polynomial
in `lambda`. Both weights are at most `K^2`. The known Pfaffian algorithm
computes either profile in polynomial time in that weight. Polynomial
multiplication and exact division in `F_2[z]` also have polynomial cost.
Thus the transformation fits the existing effective-pieces construction
budget, including the exact quotient for high occupancy.

The variable `z` only indexes pieces. The scalar polynomial is removed
while computing their base-field coefficients, before any input table is
evaluated. Extension to coefficient algebras with zero divisors therefore
introduces no forbidden division. Evaluating at `z=1` proves that the
quotient still sums to the crossing with multiplier one.

This normalizes high occupancy and can change individual coefficient
ranks. No rank-monotonicity is asserted: exact division recombines
coefficient tensors. It does not bound the ranks of the new coefficients as `lambda`
grows. In particular it does not convert a finite local saving into a
stronger deterministic factoring exponent.

## 5. Reproduce and finite validation

```sh
python3 tools/check_modular_particle_hole.py --output evidence/modular-particle-hole.json
python3 tools/check_modular_particle_hole.py --full --output evidence/modular-particle-hole-full.json
make verify
```

The default suite checks all 24 entries in the six symbolic two-by-two
identities (5), 5,456 entries of the full inverse-conjugation identity
through `K=16`, and 978 complete four-label profiles. These include exact
quotients at occupancy above one half and the equal-occupancy case where
the profile is unchanged. The full suite extends the inverse identity
through `K=32` (21,840 entries) and checks 1,714 complete profiles.

The checker also verifies the label involution and new weights, crossing
normalization at `z=1`, and small full-core Pfaffians. Omitting the scalar
factor, or omitting the reflection when the two weights differ, supplies
negative controls. The recorded hashes pin the checker and folded matrix
implementation. These finite calculations support the ordinary all-size
proof above; they do not replace it or certify an asymptotic rank gain.
