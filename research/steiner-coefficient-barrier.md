# A common coefficient excludes full Steiner words at size at least seven

For the [Steiner cubic crossing](steiner-cubic-crossings.md), the bare
profile obstruction can be strengthened to actual matrix ranks. At every
word length, any crossing decomposition into arbitrary local coefficient
recombinations has one piece whose two ranks are at least `(2T)^m` and
`(4T)^m`, where `T=q(q-1)/6`. Thus every nontrivial Steiner alphabet
`q>=7` is excluded at contraction parameter below `3/4` by this route.

The statement allows any number of pieces, arbitrary coefficient mixtures
across positions, independent polynomial variables at every position, and
any characteristic-two extension field. It concerns the complete word
domain. Restricted physical types and other local polynomial paths remain
outside the theorem. This is an ordinary proof with finite controls, not
a stronger factoring bound or new Lean work.

## 1. The full local coefficient product space

Fix a Steiner quasigroup of size `q`, and let `K(t)` be its cubic crossing.
At length `m`, allow a piece of the form

```math
A_L=L\bigl(K(t_1)\otimes\cdots\otimes K(t_m)\bigr),        \tag{1}
```

where `L` is any linear functional on the multivariate polynomials of
degree at most three in each variable. Equivalently, (1) allows arbitrary
linear combinations of the `4^m` local coefficient tensors. In particular,
it includes every recombination of the univariate coefficients obtained
by setting all parameters equal, and every position-symmetric subspace.

Define the common polynomial

```math
p(t_1,\ldots,t_m)=\prod_{i=1}^m t_i^3.                   \tag{2}
```

It is enough to prove two rank bounds for every piece with `L(p)!=0`.

## 2. First-cut blocks detect the common coefficient

The first-cut constant-basis decomposition in the linked note has `2T`
copies of the two-dimensional natural block

```math
M(t)=\begin{pmatrix}1+t^2&t^3\\t&1+t^2\end{pmatrix}.      \tag{3}
```

Retain only tensor-product blocks selecting (3) at every position. There
are `(2T)^m` independent copies of `M(t_1) tensor ... tensor M(t_m)`.
In each copy the entry using the upper-right entry at every position is
exactly `p`. After applying `L`, each copy is nonzero and therefore has
rank at least one. Direct-sum additivity gives

```math
\operatorname{rank}_{xy\mid uv}A_L\ge(2T)^m.             \tag{4}
```

Only constant row/column changes and restrictions were used. No
parameter-dependent Smith basis or coefficientwise rank substitution
is involved.

## 3. Second-cut blocks detect the same coefficient

Restrict the second cut at every position to off-diagonal physical pairs.
Each Steiner triple contributes the six-dimensional polynomial block

```math
t^2\bigl((1+t)(P+P^{-1})+tS\bigr),                      \tag{5}
```

with `P^3=S^2=I` and `SPS=P^-1`. The constant idempotent `P+P^-1` has
rank four. Its two-dimensional kernel contributes `t^3S`, which becomes
two scalar blocks `t^3` after a constant invertible column change.

On the four-dimensional image, (5) is `t^2 I+t^3(I+S)`. There `I+S` has
square zero and rank two: the characteristic-two standard representation
appears twice, and the transposition has one nonzero nilpotent Jordan
block in each copy. Consequently the image contributes two blocks

```math
J(t)=\begin{pmatrix}t^2&t^3\\0&t^2\end{pmatrix}.          \tag{6}
```

Thus each triple gives **four** constant-basis blocks, two scalar `t^3`
blocks and two copies of (6). There are `4T` blocks in all. Each has a
specified entry equal to `t^3`. Their tensor products give `(4T)^m`
independent blocks, each with a specified entry equal to `p`. Hence

```math
\operatorname{rank}_{xu\mid yv}A_L\ge(4T)^m.             \tag{7}
```

The diagonal/off-diagonal coupling of the complete cut is irrelevant:
the proof restricts to off-diagonal rows and columns before decomposing
them. Equations (4) and (7) detect the same functional value of the same
piece, even when `L` couples every position.

## 4. Every crossing decomposition contains a detected piece

The polynomial `p` occurs as a literal matrix entry after the constant
first-cut changes above. Therefore its pairing with `L` is well-defined
by the actual coefficient tensor, even if a coefficient presentation has
a kernel. At the crossing endpoint `t_1=...=t_m=1`, that entry is one.

If the endpoint crossing is the sum of coefficient pieces, their detected
entries sum to one. At least one piece has `L(p)!=0`. For that piece,

```math
\sqrt{\operatorname{rank}_{xy\mid uv}A_L\,
             \operatorname{rank}_{xu\mid yv}A_L}
 \ge (\sqrt8\,T)^m.                                    \tag{8}
```

For every admissible `q>=7`,

```math
8T^2=\frac{2q^2(q-1)^2}{9}\ge q^3,                      \tag{9}
```

because `2(q-1)^2>=9q` at seven and the difference increases thereafter.
The usual full-word cost parameter must therefore satisfy `rho>=3/4`.
This conclusion does not depend on the number of pieces, on a scalar
profile optimization, or on a symmetry restriction across positions.

The same proof permits a different fixed Steiner design at each position,
and even varying sizes `q_i>=7`: replace the two rank bounds by
`product_i(2T_i)` and `product_i(4T_i)`. It does not permit changing the
local coefficient space separately for each piece.

## 5. Exact controls and the remaining three-letter case

[The checker](../tools/check_steiner_coefficient_barrier.py) exhausts every
nonzero binary combination of the univariate coefficient tensors at
`q=7`, lengths one and two: 15 and 127 combinations. Among pieces detected
by the highest coefficient, the minimum ranks are exactly `14^m` and
`28^m`, so (4)--(7) are simultaneously sharp. Their minimum rank products
are `392` and `153664=392^2`; the key exact comparison is `392>7^3`.

The independently recorded complete univariate spanning thresholds are
separate finite observations. They do not replace the multivariate proof.
[The evidence](../evidence/steiner-coefficient-barrier.json) records these
values and the checker hash.

```sh
make verify-steiner-coefficients
make verify
```

At `q=3`, the same rank theorem holds with `T=1`, but its rate alone does
not exclude an improvement. The [bare cubic profiles](steiner-cubic-crossings.md)
already fail at that size; arbitrary coefficient recombinations on the
complete three-letter domain remain a distinct question. Restricted
physical types at any size also remain outside the full-word argument.
The complete factoring exponent remains `1/2`.
