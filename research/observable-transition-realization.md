# Realizing an observable matrix table as a crossing transfer

The diagonal rank-one transfer can use any finite matrix table with
an observable crossing quotient. It does not require a common
crossing eigenvector or commuting factor halves. This note gives
an explicit realization and a three-state example whose crossing
character has no common left or right eigenvector.

The example is not an improving factoring certificate. Its purpose
is to give an exact construction outside the common-eigenvector
subclass and a finite search interface with explicit costs.

## 1. Input: a matrix table with an observable crossing quotient

Let `T(y,u,x,v)` be an `M by M` matrix over a field, with every
external label in an alphabet of size `d`. Supply vectors `c,d_0`
such that `c^T d_0=1` and, for every sequence of inputs,

```math
c^T T(q_1)\cdots T(q_k)d_0
 =\prod_{j=1}^k [y_j=u_j][x_j=v_j].               \tag{1}
```

The empty sequence is included. The notation `d_0` distinguishes
the boundary vector from the alphabet size.

This has an exact finite certificate. Starting with `(d_0,1)`,
compute the invariant reachable span under the finitely many
block matrices

```math
\operatorname{diag}(T(q),[y=u][x=v]).             \tag{2}
```

Condition (1) holds if and only if `(c^T,-1)` annihilates this span.
At most `M+1` independent vectors can be added. Equivalently one
may supply an invariant subspace containing `d_0` on which `c^T`
transforms by the scalar crossing character.

Let `r_ij,c_ij` be the two cut ranks of the scalar tensor `T_ij`,
in cuts `(y,x)|(u,v)` and `(u,x)|(y,v)`. Define

```math
G_{ij}=\sqrt{r_{ij}c_{ij}}.                       \tag{3}
```

## 2. Explicit realization

Take an auxiliary vector space of dimension `1+M d^2`, with basis
`eta` and `e_(j,x,v)` for every state `j` and external pair `(x,v)`.
Define vectors, covectors, and a boundary functional by

```math
\begin{aligned}
a_j(x,v)&=e_{j,x,v},\\
\phi(\eta)&=1,&
\phi(e_{j,x,v})&=[x=v]c_j,\\
b_i(y,u)^T\eta&=[y=u](d_0)_i,&
b_i(y,u)^T e_{j,x,v}&=T_{ij}(y,u,x,v).
\end{aligned}                                    \tag{4}
```

These assignments define honest linear maps on the stated basis;
there are no consistency relations left to solve. Put
`A_j(x,y,u,v)=a_j(x,v)b_j(y,u)^T`.

**Realization theorem.** The matrices `A_j` satisfy the all-word
crossing endpoint with boundary `(phi,eta)`. Their transition table
is exactly the prescribed `T`. In particular, for every `k>=1`
they construct an ordinary scalar crossing decomposition on `d^k`
labels with polynomially many pieces and

```math
\mathcal S_k\le s^T G^{k-1}t,
\qquad s_i=d[c_i\ne0],\quad t_i=d[(d_0)_i\ne0]. \tag{5}
```

**Proof.** The transition identity `b_i^T a_j=T_ij` follows directly
from (4). Summing the scalarizations of all ordered `k`-letter
piece words gives

```math
[x_1=v_1][y_k=u_k]\,
c^T T(y_1,u_1,x_2,v_2)\cdots
T(y_{k-1},u_{k-1},x_k,v_k)d_0.                   \tag{6}
```

For `k=1`, the matrix product is empty. Applying (1) to (6) gives
the product of the `k` crossings. The two boundary scalar matrices
are nonzero multiples of the `d by d` identity, or zero, which
gives the ranks in (5). The exact individual-word rank products
and the grouping by initial state, final state, and all ordered
transition counts are those of
[diagonal-rank-one-word-transfer.md](diagonal-rank-one-word-transfer.md).
That construction proves (5). `QED`.

Everything in (4) is explicit. For fixed finite input data, type
entries are computable in polynomial time in `k`; complete scalar
tables and rank factorizations cost `d^O(k)`. A rational upper
matrix and positive supervector provide an exact strict spectral
certificate. If the live spectral radius is below `d^(sqrt(2))`,
sufficient amplification gives an improving finite seed.

The realization neither interchanges auxiliary operators nor
assumes that the individual pieces preserve an invariant subspace.
All cancellations required by the endpoint occur in the supplied
observable matrix table.

## 3. A strict observable-quotient example

Work over F2 or an extension field and use binary labels. Let

```math
E=[x=v][y=u],\quad R=[x=y][u=v],\quad
C=[x=u][y=v].
```

These three scalar tensors are linearly independent. Consider

```math
T_0(q)=
\begin{pmatrix}
0&R(q)&0\\
0&E(q)&C(q)\\
0&0&0
\end{pmatrix},\qquad c_0=d_0=e_1.               \tag{7}
```

The index convention is `0,1,2`. The middle diagonal entry of
every ordered matrix product is the product of the middle diagonal
entries. Thus (1) holds. The reachable space from `e_1` is
`span(e_0,e_1)`, whose `e_0` direction is unobservable by `e_1^T`.

This is strictly more general than a common crossing eigenvector.
Suppose `T_0(q)v=E(q)v` for all inputs. Comparing the independent
coefficients of `E,R,C` first forces `v` into the middle coordinate,
then the `R` coefficient forces that coordinate to vanish.
There is no nonzero such `v`. The corresponding left-eigenvector
argument uses the `C` coefficient. Both conclusions persist over
every extension field.

Conjugate by the invertible matrix

```math
W=\begin{pmatrix}1&1&1\\1&1&0\\0&1&1\end{pmatrix},
\qquad T=WT_0W^{-1}.
```

Both transformed boundary vectors are `(1,1,1)`. The entry table is

```math
T=
\begin{pmatrix}
E+R+C&E+R+C&E+R\\
E+R+C&E+R+C&E+R\\
E+C&E+C&E
\end{pmatrix}.                                   \tag{8}
```

The exact two-cut rank pairs are

| Scalar tensor | Row rank | Column rank |
| --- | ---: | ---: |
| `E+R+C` | 2 | 2 |
| `E+R` | 4 | 1 |
| `E+C` | 1 | 4 |
| `E` | 4 | 4 |

Consequently the cost matrix is

```math
G=\begin{pmatrix}2&2&2\\2&2&2\\2&2&4\end{pmatrix}.
```

Its eigenvalues are `0` and `4+-2sqrt(2)`, so its spectral base
is `4+2sqrt(2)`. A simple rational upper certificate is
`G(1,1,3/2)^T <= 7(1,1,3/2)^T`. The realization has three pieces,
auxiliary dimension thirteen, and boundary ranks `(2,2,2)` on
both sides. This is an exact strict-quotient example, not a score
improvement.

## 4. A protected diagonal in one triangular basis

The same `W` gives a useful control on a small search family.
For every upper-triangular matrix table `U` with `U_11=E`,

```math
(WUW^{-1})_{22}=E.                              \tag{9}
```

Indeed row two of `W` is `(0,1,1)`, and column two of `W^{-1}`
is `(1,1,0)`. Upper triangularity makes their contraction exactly
`U_11`. The associated cost matrix has a diagonal entry four, so
its spectral radius is at least four, irrespective of its other
scalar tensor entries. Both transformed boundary coordinates at
state two are nonzero, so the diagonal is live.

Thus none of the other five upper-triangular entries can reduce the live
spectral base below four in this fixed basis, even when those entries
range over arbitrary scalar tensors. Other bases over larger fields, larger
state spaces, and nontrivial diagonal blocks are not covered by
this protected-diagonal argument.

## 5. Verification

`check_observable_transition_realization.py` checks the exact
rank table, the common-eigenvector failure, the finite reachable
certificate, all input words through length three, the explicit
thirteen-dimensional lift, its transition entries and boundary
ranks, and the protected-diagonal identity on every triangular
matrix unit. It also checks the rational upper certificate.
The proof of the realization theorem is independent of those
finite controls.


The [linear detector theorem](linear-rank-detectors-for-observable-tables.md)
gives a stronger result for entries in the binary three-pairing span:
no auxiliary size or observable quotient can lower this spectral base
below four. Additional scalar entry tensors are a separate possibility.
