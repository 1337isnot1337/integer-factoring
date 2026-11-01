# An exact binary control for three rank-one boundary pieces

## Scope

This note tests a concrete way to escape the common-terminal limitation
of [diagonal rank-one auxiliary words](diagonal-rank-one-word-transfer.md).
It permits arbitrary four-leg tensors in its three transition rows.
Over every field of characteristic two, its optimal word-rank spectral
bound has base exactly four. The finite-label target would require a
base below `2^sqrt(2)`, so this bound supplies no score improvement.
Cancellations inside grouped type pieces remain outside the obstruction.

The coefficient field has characteristic two, the external alphabet has size two,
and the three terminal functions are the distinct nonzero functions
`D,G,D+G`, where `D(y,u)=[y=u]` and `G` is neither zero nor `D`.
This is not a classification of every three-piece boundary construction.
An exhaustive `F_2` calculation provides additional finite evidence;
the proof over extension fields uses the all-field two-piece inequality.

## 1. An exact endpoint with arbitrary transition tensors

Take an auxiliary space with independent basis consisting of `eta`
and four vectors `f_(x,v)`. Write `delta(x,v)=[x=v]`, and put

```math
a_1(x,v)=f_{x,v}+\delta(x,v)\eta,
\qquad a_2(x,v)=a_3(x,v)=f_{x,v}.
```

Choose arbitrary covector tables `b_i(y,u)` subject to

```math
b_1(y,u)\eta=D(y,u),\quad
b_2(y,u)\eta=G(y,u),\quad
b_3(y,u)\eta=(D+G)(y,u).
```

Their other entries are entirely free. Define
`A_i(x,y,u,v)=a_i(x,v)b_i(y,u)`. Because the three terminal functions
sum to zero in characteristic two,

```math
\left(\sum_i A_i(x,y,u,v)\right)\eta
       =\delta(x,v)D(y,u)\eta=E_2(x,y,u,v)\eta.
\tag{1}
```

Use `phi(eta)=1` and `phi(f_(x,v))=0`. This supplies an exact
right-boundary endpoint at every amplification length. It only has
initial letter one in the scalar word expansion.

Put `B=(D,G,D+G)` and

```math
F_i(y,u,x,v)=b_i(y,u)f_{x,v},\qquad
H_i(y,u,x,v)=B_i(y,u)\delta(x,v).
```

Each `F_i` is an arbitrary tensor with sixteen field entries,
independently of the other two. The adjacent transition tensors
from letter `i` to letters one, two, and three are respectively

```math
F_i+H_i,\qquad F_i,\qquad F_i.                      \tag{2}
```

An auxiliary word factors into its first boundary, adjacent transitions,
and final boundary. At different transitions the physical variables
are disjoint in both cuts. Its two ranks therefore multiply. Grouping
words by their initial letter and ordered transition counts gives
polynomially many pieces with a rank-sum bound governed by the
nonnegative transition matrix

```math
M=\begin{pmatrix}
x_1&y_1&y_1\\x_2&y_2&y_2\\x_3&y_3&y_3
\end{pmatrix},\qquad
x_i=\sqrt{r(F_i+H_i)c(F_i+H_i)},\quad
y_i=\sqrt{r(F_i)c(F_i)}.                            \tag{3}
```

This is the direct rank-one word construction; it does not assume
cross-half commutation of general auxiliary factors. The conclusion
below concerns this exact word-rank counting bound. Cancellations
inside a grouped type piece might lower its actual rank further.

## 2. All-field transition inequality and complete F_2 frontiers

Let `H=delta(x,v)B(y,u)` and `r=rank B`, with `r=1` or `2`.
For every tensor `F`, form the pair

```math
(p,q)=\bigl(r(F+H)c(F+H),\ r(F)c(F)\bigr).
```

Independent invertible changes of the `y,u` labels transform `B`
to a rank-`r` diagonal matrix. These changes preserve both cut ranks
and bijectively permute all possible `F`. Thus the frontier depends
only on `r`.

The [all-field binary two-piece lemma](binary-two-piece-rank-sum-and-f4-boundary.md)
proves, for `H=delta(x,v)B(y,u)`, that

```math
\sqrt{r(F+H)c(F+H)}+\sqrt{r(F)c(F)}\ge2r
\quad\hbox{over every characteristic-two field}.       \tag{4}
```

For `r=2`, invertible changes of basis reduce `H` to the binary crossing,
and the lemma applies to the two pieces `F` and `F+H`. For `r=1`,
two nonzero pieces each have cost at least one; if either is zero,
the other equals `H` and has both ranks two. Thus the inequality is
independent of the field size.

Over `F_2`, exact enumeration of all 65,536 tensors for each canonical target gives
the following complete coordinatewise minimal pairs:

| Rank of `B` | Pareto frontier of `(p,q)` |
| ---: | --- |
| 1 | `(0,4), (1,1), (4,0)` |
| 2 | `(0,16), (1,9), (4,4), (9,1), (16,0)` |

Every frontier coordinate is a perfect square, and every frontier
point satisfies `sqrt(p)+sqrt(q)=2r`. Every other pair dominates
one of the frontier pairs. This independently verifies (4) over `F_2`
and describes every optimal coordinatewise rank-product profile there.

The checker constructs both four-by-four flattenings independently,
retains a witness for each frontier point, and verifies the inequality
by integer squaring. There is no floating-point rank or eigenvalue step.
The finite certificate and the field-independent proof have separate
roles; the finite enumeration is not used to infer the extension-field result.

## 3. Exact optimal spectral base

The nonzero eigenvalues of (3) are those of

```math
K=\begin{pmatrix}
x_1&y_1\\x_2+x_3&y_2+y_3
\end{pmatrix}.
\tag{5}
```

Indeed factor `M=A B` with columns of `A` equal to `x,y` and
rows of `B` equal to `(1,0,0)` and `(0,1,1)`. Then `K=B A`.
Equation (4) makes its first row sum at least four. Its second
row sum is at least

```math
2\operatorname{rank}G+2\operatorname{rank}(D+G)
 \ge2\operatorname{rank}D=4.
```

Thus `K 1>=4 1`, and iteration gives spectral radius at least four.
This also applies to the part visible from the initial boundary:
if `y_1=0`, its first state already has a loop of weight at least
four; if `y_1>0`, both other states are reached. Every terminal
function is nonzero, so none of those states is terminally invisible.

Conversely choose `F_1=F_2=F_3=0`. Then the visible first-state loop
has weight four and all later columns of `M` vanish. Its spectral
base is exactly four. This proves the claimed optimum.

There is a concrete finite target for further work: other relations among
unequal terminal projections, more pieces, or an observable
quotient that does not require the right-boundary vector identity
may fall outside this calculation. None is asserted to improve the
score without its own ranks and endpoint proof.

## Reproduction

```sh
python3 tools/check_three_piece_binary_boundary.py \
  --output evidence/three-piece-binary-boundary.json
```

The evidence records all minimal pairs, their literal tensor masks,
and the count of verified exact inequalities. The spectral conclusion
uses the row-sum proof in Section 3.
