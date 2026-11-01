# Binary scalar projection and a complete fused coefficient envelope

The symmetric two-copy coefficient space of one literal three-letter fused
crossing cannot supply a finite crossing certificate below `1/sqrt(2)`.
This is an exhaustive statement about every binary linear decomposition
in a 21-dimensional space, rather than a search over selected polynomial
coefficients or gauges. It does not exclude higher powers, physical type
restrictions, extension-field decompositions or the larger ordered space.
The factoring score remains `1/2`.

## 1. Literal crossing and the enlarged search space

Use the field representation and construction from
[finite fused six-vertex crossings](finite-fused-six-vertex-crossings.md),
with `m=2`, `theta=4`, and

```math
F(t)=F_0+tF_1+t^2F_2\in\operatorname{Mat}_9(\mathbb F_{16}[t]),
\qquad F(1)=I_9.
```

The checker retains the three sparse matrices literally. A separate rerun
of the existing fusion constructor matched all 243 coefficient entries,
including its invariant-subspace and strict-endpoint checks. The present
rank obstruction only needs the retained literal matrices and their
endpoint. Its definition is therefore independent of a general fusion
identity.

Define `T_i(x,y,u,v)=(F_i)_(3y+x,3u+v)`. The sum of these tensors is the
crossing `E_3(x,y,u,v)=[x=v][y=u]`. For the opposite cut the checker uses
the matrix

```math
(B_i)_{x+3u,v+3y}=(F_i)_{3y+x,3u+v}.
```

These are exactly the two tensor flattenings up to fixed row/column
permutations. On two copies the physical alphabet has size nine.

Let `chi` be the projection onto the coefficient of `1` in the binary
basis `(1,omega,u,omega*u)` of `F16`. It has `chi(1)=1`; it is not the
ordinary field trace. Project entries only after forming the complete
two-copy tensor. Let `V` be the binary span of

```math
\chi(a T_i\otimes T_i),\qquad
\chi\bigl(a(T_i\otimes T_j+T_j\otimes T_i)\bigr),
\quad i<j,\quad a\in\{1,2,4,8\}.
```

The 24 displayed generators have binary dimension 21. They contain every
scalar-projected coefficient of the same-parameter polynomial
`T(t) tensor T(t)`, and allow arbitrary further binary recombination.
They also contain `E_9`. In the retained independent-generator order its
coordinate integer is `1118481`. The evidence records every basis label.

## 2. Exact rank-span separation

For `X` in `V`, let `r_R(X),r_C(X)` be its two binary ranks, each at most
81. Exhaustively checking all `2^21` coordinate vectors gives exactly
17 nonzero vectors with

```math
r_R(X)r_C(X)\le501.
```

Their span has dimension six. The binary linear functional

```math
f(x)=x_1+x_{19}+x_{20}
```

vanishes on all 17 vectors and satisfies `f(E_9)=1`; indices start at
zero. Thus every decomposition `E_9=sum_j X_j`, with all `X_j` in `V`,
contains at least one piece with rank product at least 502. Repeated
pieces or an unbounded number of summands cannot evade this conclusion.

Its rank-sum cost consequently satisfies

```math
S=\sum_j\sqrt{r_R(X_j)r_C(X_j)}\ge\sqrt{502}>9^{\sqrt2}.
```

The last inequality uses integers only: `99^2>2*70^2` and
`502^35>9^99`. The
[finite-block amplification criterion](finite-block-crossing-amplification.md)
requires `S<9^sqrt(2)` for this nine-letter block. No decomposition in
`V` can meet that criterion.

## 3. Independent replay and limits

Both cuts split into 25 exact occupation/charge blocks. The sorted block
sizes are one 9, four 6s, four 4s, four 3s, eight 2s and four 1s. The
checker verifies the absence of off-block entries for every generator.
The native exhaustive screen adds exact block ranks and can stop once
their partial product exceeds 501: subsequent block ranks are nonnegative.
Gray-code updates visit every coordinate once.

A second pass checks the entire affine hyperplane `f=1`, comprising
`2^20` vectors, with column-scanning Gaussian elimination instead of the
first pass's row-pivot algorithm. None has rank product at most 501.
Python independently reconstructs all 17 retained complete 81 by 81
matrices and checks their 34 ranks without using the block decomposition.

```sh
make verify-fused-scalar-projection
```

The command compiles a temporary C++17 executable, reconstructs the
space and both exhaustive passes, and compares the result with
[`evidence/fused-scalar-projection.json`](../evidence/fused-scalar-projection.json).
The sources are
[`check_fused_scalar_projection.py`](../tools/check_fused_scalar_projection.py)
and [`screen_fused_scalar_projection.cpp`](../tools/screen_fused_scalar_projection.cpp).
This is a finite exact computation with an ordinary linear-separation
argument, not independent mathematical review or a new Lean proof.

The full ordered two-copy envelope has dimension 33. It allows
`chi(a T_i tensor T_j)` without symmetrizing the two positions and is
not excluded by the completed 21-dimensional certificate. More generally,
projection need not commute with a nonbinary change of primary basis;
an extension-field rank obstruction cannot silently be transferred to
the projected binary space. Higher powers and physical type restrictions
remain separate questions.

Source attribution: the literal seed is from the repository's existing
finite fusion constructor. This coefficient-envelope argument and its
replays were developed with substantial OpenAI assistant support.
