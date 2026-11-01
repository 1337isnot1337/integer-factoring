# A quantum block-braid crossing compiler

The [paired-path encoders](paired-path-crossings.md) support another
polynomial crossing, built from a block braid instead of independent
binary site swaps. Its endpoint follows from constant-dimensional
cup and projector identities. This gives a different uniform family
outside the independent-site kernel. Favorable two-cut ranks are open.

Earlier [coordination work](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19324711)
checked local Temperley--Lieb relations on the classical F5 bounded
paths. This compiler uses ambient F16 matrices and explicit cup-sliding
identities to obtain its endpoint.

Work over the same fixed `F16`, with any of the twelve valid parameters
from the paired-path note. The default is `q=a^3`, of order five.
Use the height-three primal frames `F_P` and inverse duals `D_P`.

## 1. Local matrices and identities

In physical wire order, let

```math
c_q=-|10\rangle+q|01\rangle,\qquad
e_q=q^{-1}c_qc_q^T,\qquad B_q=qI-e_q.
```

Thus `e_q` vanishes on `00,11` and has mixed block
`[[q^-1,-1],[-1,q]]` in the order `10,01`. In characteristic two
the displayed minus signs equal plus signs. Direct multiplication gives

```math
e_i^2=[2]_q e_i,\quad e_i e_{i+1}e_i=e_i,\quad
B_iB_{i+1}B_i=B_{i+1}B_iB_{i+1}.                 \tag{1}
```

Disjoint operators commute. The essential cup slide is

```math
B_0B_1(c_q\otimes I)=q(I\otimes c_q).            \tag{2}
```

Products act on columns, so the rightmost gate acts first.

For `0<=h<=4`, embed the magnetic basis of `H_h` into `h` physical
wires by `i_h`: its weight-`r` column has entry
`q^(-N_01(w))` on a word `w` of weight `r`, where `N_01` counts a
zero followed later by a one. These columns have disjoint supports.
Their nonzero squared norms are

```math
g_{h,r}=q^{-r(h-r)}\frac{[h]_q!}{[r]_q![h-r]_q!}.
```

This formula follows by the recurrence obtained from appending the
last bit. All factors through `[4]_q` are units for the stated
parameters. Define `f_h=i_h diag(g_(h,r)^-1)i_h^T`.
The constant-size matrices satisfy, for `1<=h<=3`,

```math
f_{h+1}=(f_h\otimes I)
 -\frac{[h]_q}{[h+1]_q}
  (f_h\otimes I)e_{h-1}(f_h\otimes I),           \tag{3}
```

with `f_1=I`. In particular these projectors are polynomials in the
local `e_i`. The paired splitting columns obey

```math
(i_h\otimes I)U_h^+=i_{h+1},                    \tag{4}
```

and

```math
(i_h\otimes I)U_h^-
 =[h]_q(f_h\otimes I)(i_{h-1}\otimes c_q).       \tag{5}
```

For clarity, (1)--(5) are identities of ordinary finite-field
matrices. The replay checks them for every stated parameter, together
with sliding `f_h` past one wire. Its largest projector has four
physical wires. No quantum-group classification or categorical trace
is used as an unproved endpoint assumption.

## 2. The all-length block endpoint

For a retained prefix path `P` of length `n` ending at height `h`,
equations (3)--(5) express its embedding as
`F_P=A_P i_h`, where `A_P` maps `h` wires to `n` wires and is built
from cups, the projectors `f_j`, and scalar factors. This follows
inductively: an up step uses (4) and `f_(h+1)`; a down step uses
(5) and a cup. The map has `(n-h)/2` added cups.

Let `B_(n,m)` move a block of `n` wires past `m` wires. Its application
word, in zero-based adjacent indices, is

```text
for j = 0,...,m-1:
    for i = n+j-1,...,j:
        apply B_i
```

There are `nm` gates. Sliding this word through the expression for
`A_P`, using (1) to move the local projectors and (2) for every cup,
gives

```math
B_{n,m}(A_P\otimes I_m)
 =q^{m(n-h)/2}(I_m\otimes A_P)B_{h,m}.           \tag{6}
```

This is induction on the cup/projector expression: each added cup
crosses `m` wires and contributes `q^m`; length-preserving projectors
slide without a scalar. The block-word concatenations are obtained by
the braid relation and disjoint commutations. At final height zero,
`B_(0,m)=I`, so (6) holds against every input on the other block.
In particular,

```math
B_{n,n}(F_P\otimes F_Q)
 =q^{n^2/2}(F_Q\otimes F_P).                    \tag{7}
```

The inverse-dual pairing now converts (7) to the exact label crossing.

## 3. The polynomial compiler and its cost

Replace each gate in the block word by `qI-z e_i`. Define

```math
G_n(z;P,Q,R,S)=q^{-n^2/2}
 (D_R\otimes D_S)\prod(qI-z e_i)(F_P\otimes F_Q).
```

Then `G_n(1)` is the crossing, while
`G_n(0)=q^(n^2/2)[P=R][Q=S]`.
Multiplying by the affine scalar

```math
b_n(z)=q^{-n^2/2}+(1-q^{-n^2/2})z
```

makes both endpoints strict. This changes the degree to at most
`n^2+1` and preserves the complete coefficient span: `b_n(0)` is a
unit, so successive coefficient equations recover the old coefficients
from the new ones. Each direction of span containment follows.

The [implemented constructor](../tools/quantum_block_braid.py) applies
the sparse two-wire matrices to an ambient coefficient vector and
contracts with the duals. It costs `2^O(n)` bit operations, including
all `O(n^2)` coefficients. It is not claimed to be a polynomial-time
magnetic entry algorithm. Binary-selected paths have `n=O(lambda)`,
so the constructor and dense cut factorizations fit the permitted
`2^O(lambda)` effective-family budget of the
[modular interface](bounded-spin-modular-interface.md).

There are polynomially many explicit raw coefficient pieces and they
sum to the crossing. The remaining rank target is
`sum sqrt(r_j*c_j) <= poly(lambda) 2^((sqrt(2)-epsilon)lambda)`
for some fixed positive `epsilon`. No such bound is established here.
Arbitrary coefficient selection in an `O(lambda^2)`-dimensional span
also needs an effective recipe; exhaustive selection is too costly.

Run `make verify-quantum-braid`. The standard-library
[receipt](../evidence/quantum-block-braid.json) pins the local-identity
checks, small polynomial endpoint controls, and comparisons against
separate direct gate evaluations at non-endpoint parameters. The
optional `make verify-quantum-braid-profiles` constructs every entry
at four and six sites and both cut matrices of every coefficient.
Its [complete profiles](../evidence/quantum-block-braid-profiles.json)
give coefficient-span dimensions five and 33 respectively, with all
strict endpoints checked. These are structural controls below the
[universal 95-label finite exclusion](geometric-cut-rank-budget.md).
These checks
support the all-length argument above; they do not certify a favorable
rank profile or a higher score. The complete score remains `1/2`.
The construction, argument and tools used substantial Codex AI assistance.
