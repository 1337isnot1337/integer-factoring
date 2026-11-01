# A four-piece F4 boundary pencil with two free ket functions

This is a separate extension of the
[binary two-piece and F4 boundary note](binary-two-piece-rank-sum-and-f4-boundary.md).
The construction has two independent free ket functions. Its exact
transition spectral minimum is still four, over F4 and every
extension field. The proof covers arbitrary free scalar tensors.

## 1. Construction and transition table

Work in characteristic two with `omega^2+omega+1=0`. Put

```math
P(y,u)=[y=u=0],\quad Q(y,u)=[y=u=1],
\quad (\beta_1,\beta_2,\beta_3,\beta_4)
 =(P,Q,P+\omega Q,P+\omega^2Q).
```

Use an auxiliary space with basis `eta`, four vectors `f(x,v)`,
and four vectors `g(x,v)`. Define

```math
\begin{aligned}
a_1&=f+g+[x=v]\eta,&
a_2&=\omega f+\omega^2g+[x=v]\eta,\\
a_3&=f,&a_4&=g.
\end{aligned}                                    \tag{1}
```

For each `i`, choose arbitrary scalar tensors `F_i,G_i` and let

```math
b_i(y,u)^T\eta=\beta_i(y,u),\quad
b_i(y,u)^Tf(x,v)=F_i(y,u,x,v),\quad
b_i(y,u)^Tg(x,v)=G_i(y,u,x,v).                    \tag{2}
```

Then `A_i=a_i b_i^T` obeys the vector endpoint

```math
\left(\sum_i A_i(x,y,u,v)\right)\eta
 =[x=v][y=u]\eta.                               \tag{3}
```

Indeed the coefficients of `f` and `g` vanish by the two displayed
relations among the `beta_i`, and the coefficient of `eta` is
`delta_(x,v)(P+Q)`.

With `phi(eta)=1` and `phi(f)=phi(g)=0`, the initial boundary ranks
are `(2,2,0,0)`, and the terminal ranks are `(1,1,2,2)`.
Writing `D_i=[x=v]beta_i(y,u)`, transition row `i` is

```math
(A,B,F,G)
 =(F_i+G_i+D_i,\ \omega F_i+\omega^2G_i+D_i,
   F_i,\ G_i).                                   \tag{4}
```

Here the letters `A,B` denote scalar transition tensors, not the
auxiliary rank-one matrices. Let `sigma(T)` be the geometric mean
of the two transition-cut ranks, as in the companion note.

## 2. The lower vector

The nonnegative transition matrix has entries `sigma(T_ij)`.
We prove the lower-vector inequality

```math
Gh\ge4h,\qquad h=(1,1,2,2).                     \tag{5}
```

### Rank-one terminal rows

For `i=1,2`, every nonzero scalar multiple of `D_i` has cost two.
If both `F,G` are nonzero, their two weighted costs already total
at least four. If both are zero, `A=B=D_i` also gives four.

If exactly one is nonzero, call it `H`. The other two transitions
are `H+D_i` and `lambda H+D_i`, with `lambda` equal to `omega`
or `omega^2`. When both are nonzero, their costs are at least one
each and the weighted cost of `H` is at least two. If one vanishes,
`H` and the remaining transition are nonzero multiples of `D_i`;
their weighted total is six. They cannot both vanish because
`lambda` is different from one. This proves the first two rows
of (5), including every zero-transition case.

### Rank-two terminal rows

For `i=3,4`, the matrix `beta_i` is invertible. The weighted
crossing `D=D_i` has both cut ranks four. In addition to (4),

```math
\begin{aligned}
D&=A+F+G=B+\omega F+\omega^2G,\\
D&=\omega A+\omega^2B+\omega^2F
  =\omega^2A+\omega B+\omega G.                  \tag{6}
\end{aligned}
```

Thus every triple among `A,B,F,G` spans `D`, with nonzero
coefficients on all three tensors.

If `A=0` or `B=0`, the two-piece inequality applied to the other
two free tensors gives `sigma(F)+sigma(G)>=4`. Their doubled
weights give the required total eight. If `F=0`, apply that
inequality to `(A,G)` and to `(B,omega^2 G)` and add; the result
is `sigma(A)+sigma(B)+2sigma(G)>=8`. The case `G=0` is identical.

Now suppose all four tensors are nonzero. If
`sigma(F)+sigma(G)>=3`, their weighted costs total at least six
and the other two costs total at least two. Otherwise both free
costs are less than two. A nonzero binary four-leg tensor with
cost less than two has rank pair only `(1,1)`, `(1,2)`, or `(2,1)`:
the pairs `(1,3),(3,1)` are impossible by the rank-one Kronecker
identity from the companion lemma.

The triple rank budgets from (6) give the following exhaustive
possibilities, allowing cut and free-tensor exchanges:

| Rank pairs of `F,G` | Lower bound on `sigma(A)+sigma(B)` | Total weighted lower bound |
| --- | ---: | ---: |
| `(1,1),(1,1)` | `4` | `8` |
| `(1,1),(1,2)` | `2+sqrt(2)` | `4+3sqrt(2)` |
| `(1,2),(1,2)` | `2sqrt(2)` | `6sqrt(2)` |
| `(1,2),(2,1)` | `2sqrt(2)` | `6sqrt(2)` |

For the first row, `A,B` each have both ranks at least two.
For the second, their row ranks are at least two and their column
ranks sum to at least three. For the third, both row ranks are at
least two. For the fourth, each cut-rank sum of `A,B` is at least
three. That last condition gives either one cost at least two
and the other at least one, or two costs at least `sqrt(2)`.
Every table entry is at least eight. This proves the remaining
rows of (5).

All terminal ranks are positive. Restricting to the vertices
reachable from the two nonzero initial boundaries preserves (5),
because edges from a reachable vertex to an unreachable vertex
are zero. The live spectral radius is consequently at least four.

Taking every `F_i=G_i=0` leaves precisely the live matrix
`[[2,2],[2,2]]` on the first two states, so four is attained.

## 3. Exact controls

`check_four_piece_f4_boundary.py` enumerates the finite integer
rank profiles allowed by all four triple budgets. It checks their
weighted lower bound with integer square-root intervals. It also
checks the F4 linear identities (3), (4), and (6), actual cut ranks,
boundary ranks, and the lower vector for zero cases and arbitrary
sampled free tensors. The universal theorem is the proof above,
not the finite sample.
