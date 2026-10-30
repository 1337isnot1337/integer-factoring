# A rank obstruction for coordinate matching support

## Result and scope

Restricting a tensor product of the elementary binary crossing to a code
does not improve its asymptotic rate. The restriction can use any subset
of binary strings, including nonlinear codes and constant-weight sets.
The result also permits arbitrary coefficients within the specified support.

Let `S` be a set of `m` distinct binary strings of length `t`.
Work over any field. For a four-index tensor `B` on `S^4`, use row
indices `(x,y)` and column indices `(v,u)`. Its other flattening has
indices `(x,u)` and `(v,y)`. Write their ranks as `r(B)` and `c(B)`.
Define the number of nonzero crossing-diagonal entries by

```math
d(B)=|\{(x,y)\in S^2:B(x,y,y,x)\ne0\}|.
\tag{1}
```

Fix a binary orientation mask `J` of length `t`.
At each coordinate, permit the crossing pairing `x=v, y=u`.
If `J_k=0`, also permit the pairing `x=y, v=u`.
If `J_k=1`, instead permit the pairing `x=u, v=y`.
Require every nonzero entry of `B` to satisfy one of the two permitted
pairings at every coordinate. The mask is fixed for that tensor.

**Support theorem.** Every tensor with this support satisfies

```math
r(B)c(B)\ge\frac{d(B)^2}{m}.
\tag{2}
```

**Crossing consequence.** Suppose `M` such tensors sum to a nonzero
multiple of `[x=v][y=u]`. Each tensor can have its own orientation mask.
Then

```math
\max_j r(B_j)c(B_j)\ge\frac{m^3}{M^2}.
\tag{3}
```

For `m=2^lambda` and `log2 M=o(lambda)`, their largest width sum is
at least `3lambda-o(lambda)`. Thus this construction class requires
`rho>=3/4` in a bound `r_bits+c_bits<=4rho lambda+o(lambda)`.
It cannot improve on `1/sqrt(2)`.

This is a restriction on a support pattern, not on all crossing tensors,
all modular constructions, or factoring algorithms. A sum of tensors
with different masks, retained as a single piece, can leave the class.
The negative control in Section 5 shows that this distinction matters.

## 1. The row and column matrices split into XOR blocks

Every permitted coordinate pairing implies
`x xor y = v xor u`. Thus the row flattening is block diagonal.
For a binary string `a`, define

```math
X_a=\{x\in S:x\mathbin{\mathrm{xor}}a\in S\}.
```

The corresponding row block has rows `x in X_a`, representing
`(x,x xor a)`, and columns `v in X_a`, representing `(v,v xor a)`.
The column flattening has the same block index, using `a=x xor u`.

In a row block, put `b=x xor v`.
At a coordinate with `J_k=0`, a change `b_k=1` is possible only when
`a_k=0`. At a coordinate with `J_k=1`, it is possible only when `a_k=1`.
Consequently a nonzero row entry requires

```math
x|_{T_a}=v|_{T_a},\qquad T_a=\operatorname{supp}(a\mathbin{\mathrm{xor}}J).
\tag{4}
```

In the column block, the condition reverses:

```math
x|_{T_a^c}=v|_{T_a^c}.
\tag{5}
```

Equations (4)--(5) can also be checked from the two pairings directly.
They require no linear structure on `S`.

## 2. Nonzero diagonal entries force projection classes

Within block `a`, let

```math
D_a=\{x\in X_a:B(x,x\mathbin{\mathrm{xor}}a,
                   x\mathbin{\mathrm{xor}}a,x)\ne0\}.
```

By (4), the row block splits further according to projection onto `T_a`.
Each projection class occupied by `D_a` contains a nonzero diagonal
entry, so its matrix block has rank at least one.
The same argument with (5) applies to the column block.
If their ranks are `r_a,c_a`, then

```math
r_a\ge|\pi_{T_a}(D_a)|,\qquad
c_a\ge|\pi_{T_a^c}(D_a)|.
```

The two complementary projections together determine `x`.
Therefore

```math
r_ac_a\ge |D_a|.
\tag{6}
```

Ranks add over the XOR blocks. The ordinary Cauchy--Schwarz inequality
for the nonnegative real numbers `sqrt(r_a)` and `sqrt(c_a)` gives

```math
\begin{aligned}
r(B)c(B)
&=\left(\sum_a r_a\right)\left(\sum_a c_a\right)\\
&\ge\left(\sum_a\sqrt{r_ac_a}\right)^2
\ge\left(\sum_a\sqrt{|D_a|}\right)^2.
\end{aligned}
```

Each `|D_a|<=m`, so `sqrt(|D_a|)>=|D_a|/sqrt(m)`.
Also `sum_a |D_a|=d(B)`. Substitution proves (2).
The use of a real inequality compares integer ranks; it places no
restriction on the field of the tensor entries.

## 3. Cover the crossing diagonal

If `sum_j B_j` is a nonzero multiple of the crossing, every one of
its `m^2` crossing-diagonal entries is nonzero. At least one summand
is nonzero at each such entry. Hence

```math
\sum_{j=1}^M d(B_j)\ge m^2.
```

Some `j` has `d(B_j)>=m^2/M`. Equation (2) now gives (3).
This argument permits cancellation elsewhere, arbitrary field coefficients,
and a different fixed orientation mask for every summand.

Taking base-two logarithms gives the precise tradeoff

```math
\max_j\bigl(\log_2 r(B_j)+\log_2 c(B_j)\bigr)
\ge3\log_2m-2\log_2M.
\tag{7}
```

Any exact factorization needs at least the corresponding matrix rank
many channels. The same lower bound therefore applies to its bit widths.
It includes the enlarged subexponential-piece interface in
[subexponential-piece-transfer.md](subexponential-piece-transfer.md).

## 4. Application to the elementary characteristic-two family

For binary coordinates over `F_2`, put

```math
E=[x=v][y=u],\qquad R=[x=y][v=u],\qquad C=E+R.
```

The one-bit ranks are `(1,4)` for `R` and `(4,1)` for `C`.
For any coordinate subset `S`, consider the exact coefficient family

```math
A_i=[z^i]\prod_{k=1}^t(R_k+zC_k),\qquad 0\le i\le t.
\tag{8}
```

Evaluation at `z=1` proves `sum_i A_i=E` on `S^4`.
Every coefficient has the support from Section 1 with `J=0`.
Every linear combination of these coefficients has that support as well.
Thus (3) applies to any decomposition of the crossing into `M` tensors
from their linear span, not just to the displayed coefficients.

The same conclusion holds after coordinate repetitions, long coordinate
encodings, arbitrary nonlinear label restrictions, and separate fixed
coordinate exchanges in each resulting piece. It depends on `M`, not
on the length of the coordinate encoding.

### An exact coefficient formula

For a row block, write `a=x xor y=v xor u` and `b=x xor v`.
The local two-term polynomial is zero if `a` and `b` overlap.
Otherwise its product is

```math
z^{|a|}(1+z)^{|b|}.
```

Consequently

```math
A_i(x,y,u,v)
=[x\mathbin{\mathrm{xor}}y=v\mathbin{\mathrm{xor}}u]
 [a\mathbin{\mathrm{and}}b=0]
 \binom{|b|}{i-|a|}\pmod2.
\tag{9}
```

This formula computes the entire matrix without intermediate tensor
coordinates. It independently reproduces the direct local product.
At `i=|a|`, the row block is exactly the equality matrix for projection
onto `supp(a)`. The column block is exactly the equality matrix for the
complementary projection. This gives a direct special case of the proof.

Finite searches found smaller rank products for selected label sets and
coordinate repetitions. Equation (7) explains why such improvements
cannot scale through this support class with subexponentially many pieces.
The larger exhaustive search was stopped after the support theorem closed
its asymptotic target. Finite search minima are not optimality certificates.

## 5. Why mixing supports is a distinct open route

For one bit, let

```math
F(x,y,u,v)=[x\mathbin{\mathrm{xor}}y
              \mathbin{\mathrm{xor}}u\mathbin{\mathrm{xor}}v=0].
```

Both flattenings have rank two over `F_2`, while `d(F)=4` and `m=2`.
Thus `r(F)c(F)=4<8=d(F)^2/m`.
The support hypothesis fails: all three pairings occur.

Indeed, if `S_1=[x=u][v=y]`, then `F=E+R+S_1` over `F_2`.
A mixed-support sum retained as one piece can therefore evade (2).
The theorem does not exclude a scalable construction using such sums.
It does exclude treating separate support masks as a sufficient escape.

## 6. Reproduction and attribution

Run:

```sh
python3 tools/check_coordinate_support_rank.py
```

The checker forms both full flattenings and computes exact ranks.
It checks random supported tensors over several fields, arbitrary masks,
nonlinear label sets, coefficient formulas, and the complete crossing sums.
It also records the mixed-support negative control.
These finite checks supplement the all-size argument above.

The elementary one-bit crossing and its tensor coefficients are prior
constructions discussed in the pinned Anthropic Appendix H and in
[modular-crossing.md](modular-crossing.md).
This support proof and checker used substantial OpenAI assistant support.
Two separate research agents reviewed the argument. No Lean work was run.
