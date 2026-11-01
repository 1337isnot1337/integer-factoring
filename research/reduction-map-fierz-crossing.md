# A reduction-map Fierz crossing and its exact word bound

## An explicit two-piece boundary identity

Let `d>=2` be even and work over `F_2`, or a field extension. The
auxiliary vector space is `V=Mat_d`, of dimension `b=d^2`. Put

```math
\eta=I_d,\qquad \phi(X)=X_{00},\qquad
W(X)=X+\operatorname{Tr}(X)I_d.                    \tag{1}
```

Because `d=0` in the field, `W^2=I`, `W eta=eta`, and `phi eta=1`.
Write `L_M(X)=MX` and `R_N(X)=XN`. Define

```math
A(x,y,u,v)=L_{e_{xy}}R_{e_{uv}},\qquad
B(x,y,u,v)=W L_{e_{xu}}R_{e_{yv}}W.                \tag{2}
```

Then

```math
A\eta=[y=u]e_{xv},\qquad
B\eta=[y=u](e_{xv}+[x=v]I_d),
\qquad (A+B)\eta=E_d\eta.                        \tag{3}
```

This is a concrete boundary-stable crossing identity. The row half
algebras are original left and right multiplications. The column
half algebras are their common conjugates by `W`. They commute
across their respective cuts for every piece and every channel.

## 1. Exact extended-cut profiles

Piece `A` has one row channel, and piece `B` has one column channel.
Set `q=d(d-1)/2`. Their other exact extended-cut ranks are

```math
(r_A,c_A)=(1,d^2q),\qquad
(r_B,c_B)=(d^2q,1).                               \tag{4}
```

For example, conjugating `A` into the column factor structure gives

```math
WAW(X)=(X_{yu}+[y=u]\operatorname{Tr}X)
             (e_{xv}+[x=v]I).
```

Its coefficient at output matrix entry `(a,c)` and input entry
`(b,e)` is

```math
([a=x][c=v]+[x=v][a=c])
([b=y][e=u]+[y=u][b=e]).                            \tag{5}
```

In the extended column cut, regroup the row indices as `(x,a),(u,b)`
and the column indices as `(v,c),(y,e)`. Equation (5) is a Kronecker
product. The first factor is `I+R`, where `R` is the rank-one
diagonal-vector outer product and `R^2=dR=0`; it has rank `d^2`.
The second is `I+Swap`. It vanishes on the diagonal pairs and has
rank one on each two-element orbit of distinct letters, so its rank
is `q`. This proves `c_A=d^2q`. Exchanging the external middle legs
proves the same row rank for `B`.

Thus the direct commuting-half rank-sum base is `2d sqrt(q)`.
It is not competitive with the desired threshold.

## 2. Rank-one auxiliary words give a sharper bound

Each matrix in (2) has auxiliary rank one. With `epsilon_(a,b)`
denoting the matrix-entry functional, write

```math
A=|e_{xv}\rangle\epsilon_{yu},\qquad
B=|We_{xv}\rangle\epsilon_{uy}W.                  \tag{6}
```

The [diagonal rank-one word theorem](diagonal-rank-one-word-transfer.md)
therefore applies directly. At an adjacent pair of positions, use
the four variables `(y,u,x',v')`. The four transition kernels are

| Transition | Kernel |
|---|---|
| `AA` | `[y=x'][u=v']` |
| `AB` | `[y=x'][u=v']+[y=u][x'=v']` |
| `BA` | `[u=x'][y=v']+[y=u][x'=v']` |
| `BB` | `[u=x'][y=v']` |

Their exact cut ranks are

| Transition | Row rank | Column rank |
|---|---:|---:|
| `AA` | 1 | `d^2` |
| `AB` | `d^2` | `q` |
| `BA` | `q` | `d^2` |
| `BB` | `d^2` | 1 |

These are the same three elementary matrices just used in (5):
a rank-one pairing, the invertible `I+R`, and `I+Swap`.

The initial boundary function has rank one for `A`, and rank
`d-1` for `B`, since it is respectively
`[x=0][v=0]` and `[x=0][v=0]+[x=v]`. Both final boundary
functions are `[y=u]`, of rank `d`. Boundary ranks are the same
in both cuts. The weighted transition matrix is

```math
G=d\begin{pmatrix}1&\sqrt q\\\sqrt q&1\end{pmatrix}.
```

Every row sums to `beta=d(1+sqrt(q))`. Grouping words by their
initial symbol and all four transition counts gives at most
`2(k+1)^4` pieces and the exact upper bound

```math
\boxed{\mathcal S_k\le d^2
                  [d(1+\sqrt q)]^{k-1}.}          \tag{7}
```

The end symbol is determined by the initial symbol and transition
counts; it can instead be included explicitly without affecting
the polynomial piece count. Every factor in the scalar word uses
disjoint physical variables, so the individual rank products used
in this bound are exact.

For `d>=4`, (7) strictly improves the raw matrix-channel base
`2d sqrt(q)`. For instance the base decreases from approximately
`19.59592` to `13.79796` at `d=4`, and from `46.47580` to `29.23790`
at `d=6`. This is an actual saving for the displayed construction,
but it does not improve the factoring score: even this sharper
base is at least `d^(3/2)`, exceeding `d^(sqrt(2))`.

This last comparison describes the available upper certificate;
it is not a lower bound on all possible further recombinations of
these tensors.

## 3. Why symbol counts alone are insufficient

At `d=4`, compare the words `AAB` and `ABA`. Their total numbers
of each symbol agree. Nevertheless their exact cut ranks are
respectively `(64,384)` and `(384,384)`. The ratios differ, so
grouping those words and summing their separate geometric means
would use Cauchy--Schwarz in the wrong direction.

Recording all transition counts makes their types distinct and
restores the valid argument. The checker verifies this negative
control with exact integers.

## 4. Scope of the example

Both covectors in (6) have the same terminal projection on `eta`,
namely `[y=u]`. Also the two ket functions sum to `[x=v]eta`.
Consequently each transition row itself sums to a scalar crossing
tensor. Its rank sum is exactly `d(1+sqrt(q))`.

The common-terminal collapse theorem therefore applies: the word
bound supplies no advantage over an ordinary same-width transition
seed. This is a useful explicit noncommutative boundary identity
and a sharp improvement in its word accounting. A construction
that escapes this collapse must change the terminal projections
or use a genuinely nontrivial observable quotient endpoint.

## Verification

`check_reduction_map_fierz.py` checks the matrix-vector endpoint,
extended ranks and their exact factorizations, transition identities
and ranks, direct individual word ranks, transition-count grouping,
the weighted path formula, and the symbol-count negative control.
The evidence is recorded in `reduction-map-fierz.json`.
