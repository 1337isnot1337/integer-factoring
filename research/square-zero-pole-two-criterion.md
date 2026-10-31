# An exact pole-two criterion for arbitrary square-zero physical maps

## Constructive result

Consider either the four-letter or sixteen-letter reflection crossing,
in the Bell coordinates of
[coupled-rational-physical-frames.md](coupled-rational-physical-frames.md).
Let `H` be the diagonal matrix whose entry on a letter is its bit parity:

```math
H_{ss}=\operatorname{popcount}(s)\bmod2.
```

For **every** matrix `G` over a field of characteristic two satisfying
`G^2=0`, apply `A=I+w^-2 G` to legs `x,u` and `A^-T` to legs `y,v`.
The complete transformed tensor is polynomial if and only if

```math
GHG=0.                                                   \tag{1}
```

This criterion includes arbitrary-rank and dense matrices.  Accepted
frames preserve the exact degree and strict crossing endpoint.  They
require new complete weighted profiles before any pressure improvement
can be asserted.  The theorem is instantiated here for these two fixed
reflection seeds; no all-interval low-jet identity is assumed.

## 1. Third-cut coordinates and exact low jets

Flatten the Bell-normalized tensor with rows `(x,v)` and columns `(y,u)`.
It acts on the matrix space by

```math
(M(w)X)_{xv}=\sum_{y,u}T(w;x,y,u,v)X_{yu}.
```

For a matrix `B`, write `L_B(X)=BX`, `R_B(X)=XB`, and
`ad_B=L_B+R_B`; the last expression is a commutator in characteristic
two.  The constant Bell product gives `M_0=I`.

The exact source coefficients through degree seven satisfy

```math
\begin{aligned}
M_2&=ad_H,\\
M_4&=ad_J+L_HR_H+R_H,\\
M_6&=ad_K+L_HR_J+L_JR_H+L_HR_H+R_H,\\
M_1&=M_3=M_5=M_7=0.                         \tag{2}
\end{aligned}
```

Here `J_xy=T_4(x,y,0,0)` and `K_xy=T_6(x,y,0,0)`.  All entries
of `H,J,K` lie in `F2`.  The checker extracts `J,K`, checks the parity
description of `H`, and verifies every entry of all eight coefficients.
Thus (2) is an exact finite identity for each displayed seed.  It
remains true after any extension of the coefficient field.

For the four-letter example,

```math
H=\operatorname{diag}(0,1,1,0),\qquad
J=\begin{pmatrix}
0&0&0&0\\0&0&1&0\\0&0&1&0\\0&0&0&0
\end{pmatrix}.
```

The complete `H,J,K` matrices for both alphabets are retained in the
exact evidence.  Keeping the single right-multiplication terms in
(2) matters: replacing `M_4` by only `ad_J+L_HR_H` would be false.

## 2. Conjugation gives the obstruction

Write `f=w^-2` temporarily, so `A=I+fG`.  Since `G^2=0`,
`A^-1=A`.  In the third flattening the coupled four-leg action is

```math
M\longmapsto \operatorname{Ad}_A M\operatorname{Ad}_A^{-1},
\qquad \operatorname{Ad}_A(X)=AXA^{-1}.
```

It conjugates each matrix appearing in `L`, `R`, or `ad`.  In particular

```math
H'=AHA^{-1}=H+f(GH+HG)+f^2GHG.                 \tag{3}
```

Put `B=GHG`.  In the transformed `M_4`, the term `L_H'R_H'` has
coefficient `L_B R_B` at degree four in `f`.  The terms `ad_J'` and
`R_H'` have degree at most two.  Multiplication by `w^4`, followed by
`f=w^-2`, therefore contributes

```math
[w^{-4}]M'(w)=L_B R_B.                         \tag{4}
```

No other coefficient contributes to (4).  The constant operator is
unchanged; the transformed `M_2` has degree at most two in `f`; the
odd coefficients in (2) vanish; and every coefficient with original
degree at least six has at most four factors of `f`.

If `B` is nonzero, `L_B R_B` is nonzero.  Indeed, its matrix entries
in matrix-unit coordinates are products `B_xy B_uv`; choose any two
nonzero entries of `B`.  Thus (4) is an actual pole, proving the
necessity of (1).

## 3. The same condition is sufficient

Assume `GHG=0`.  Equation (3) then has degree at most one in `f`.
The matrices `J'=AJA^-1` and `K'=AKA^-1` have degree at most two.
Using every term in (2) gives the following bounds:

| Original coefficient | Degree in `f` after conjugation | Lowest exponent after `f=w^-2` |
|---|---:|---:|
| `M_0` | 0 | 0 |
| `w^2 M_2` | at most 1 | 0 |
| `w^4 M_4` | at most 2 | 0 |
| `w^6 M_6` | at most 3 | 0 |
| `w^a M_a`, `a>=8` | at most 4 | at least 0 |

The final row holds for every tensor coefficient because the physical
action has four factors, each linear in `f`.  All odd coefficients
below eight vanish.  Every Laurent coefficient is therefore
nonnegative, proving polynomiality.

The argument applies to every square-zero `G` at once.  It does not
require diagonalizability, rank one, or parity preservation of `G`.
The endpoint and degree conclusions follow from the general coupled
frame proof: the four maps are paired inverses at `w=1`, and only the
identity term can contribute to the original highest degree.

## 4. Useful parameterizations

Condition (1) has the equivalent subspace form

```math
H(\operatorname{im}G)\subseteq\ker G.           \tag{5}
```

Any square-zero matrix commuting with `H` satisfies it.  In particular,
if `P` is an involutive permutation that preserves even and odd
letters, then `G=I+P` is square zero in characteristic two and passes
the criterion.  Perfect matchings within the two parity classes give
rank-two generators at alphabet four and rank-eight generators at
alphabet sixteen.

The family is larger than the commuting case.  At alphabet four,
`G=2E_(0,2)+3E_(1,3)` is square zero and does not commute with `H`,
but its image is `span(e_0,e_1)`, an `H`-invariant subspace contained
in its kernel.  It therefore also gives a polynomial pole-two frame.

For `G=u v^T`, condition (1) reduces to `v^T H u=0`.  In dimension
four over `F4`, exactly 1,755 nonzero rank-one square-zero generators
satisfy it.  There are ten projective `u` lying entirely within one
of the two parity planes, each with 63 nonzero annihilating `v`.
The other 75 projective `u` have nonzero parts in both planes; the
two independent conditions `v^Tu=v^THu=0` leave 15 nonzero `v`.
Hence the count is `10*63+75*15=1755`, agreeing with the complete
[rank-one catalogue](dense-rank-one-coupled-frames.md).

## 5. Exact evidence and scope

The checker verifies 526,336 low-jet coefficient entries across the
two fixed seeds.  It also compares 197,376 entries of the explicit
negative coefficient (4) for selected arbitrary-rank generators,
including accepted noncommuting maps and rejected opposite-parity
matchings.  The selected complete transformed tensors are checked
for polynomiality, unchanged degree, endpoint, and scalar content.

```sh
python3 tools/check_square_zero_second_jet.py \
  --output evidence/square-zero-pole-two-criterion.json
```

The [checker](../tools/check_square_zero_second_jet.py) and
[evidence](../evidence/square-zero-pole-two-criterion.json) retain the
exact source digests and matrices.  The construction uses
[square_zero_coupled_reflection_frames.py](../tools/square_zero_coupled_reflection_frames.py).
The criterion produces additional admissible tensors; it supplies no
pressure certificate or improved factoring machine by itself.
