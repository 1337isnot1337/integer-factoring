# Diagonal-pair rank-one auxiliary transfers

## A direct interface without commuting factor halves

Let the external alphabet have size `d`. Fix vectors
`a_i(x,v)` and covectors `b_i(y,u)^T` in an auxiliary space of
dimension `b`, and define

```math
A_i(x,y,u,v)=a_i(x,v)b_i(y,u)^T,
\qquad C(q)=\sum_{i=1}^M A_i(q).                  \tag{1}
```

The adjective rank one refers to the auxiliary matrix at each fixed
external input. The two pairs of external variables in (1) are
essential to the theorem.

Supply a common boundary `phi,eta` and an exact finite certificate
for the all-word endpoint

```math
\phi C(q_1)\cdots C(q_k)\eta=\prod_{j=1}^k E_d(q_j).
\tag{2}
```

For example, the right-eigenvector or observable-quotient certificates
in [boundary-stable-catalytic-seeds.md](boundary-stable-catalytic-seeds.md)
suffice. No global cross-half commutation assumption is needed for
this particular rank-one shape.

Define the scalar transition tensor

```math
T_{ij}(y,u,x',v')=b_i(y,u)^T a_j(x',v').            \tag{3}
```

Let `r_ij` be its rank in the cut `(y,x')|(u,v')`, and let
`c_ij` be its rank in the cut `(u,x')|(y,v')`. Also let

```math
s_i=\operatorname{rank}_{x,v}(\phi a_i(x,v)),
\qquad t_i=\operatorname{rank}_{y,u}(b_i(y,u)^T\eta),
\qquad G_{ij}=\sqrt{r_{ij}c_{ij}}.                 \tag{4}
```

All ranks are ordinary ranks over the specified coefficient field.
The transition matrix `G` is nonnegative and real, even when the
underlying tensor coefficients lie in positive characteristic.

**Theorem.** For every `k>=1`, the finite data construct an ordinary
crossing decomposition on `d^k` labels, with at most
`M^2(k+1)^(M^2)` pieces and rank sum bounded by

```math
\boxed{\mathcal S_k\le s^T G^{k-1}t.}             \tag{5}
```

The finite data and every transition rank can be checked exactly.
If the live transition matrix has spectral radius strictly below
`d^(sqrt(2))`, this gives an improving finite-seed construction after
sufficient amplification. An approximate numerical spectral radius
is not a certificate; a strict rational supervector bound as described
below is one way to verify the required inequality.

## 1. Exact word ranks

An ordered word `i_1,...,i_k` scalarizes to

```math
[\phi a_{i_1}(x_1,v_1)]
 \prod_{j=1}^{k-1}
 T_{i_j i_{j+1}}(y_j,u_j,x_{j+1},v_{j+1})
 [b_{i_k}(y_k,u_k)^T\eta].                         \tag{6}
```

This is the ordinary product rule for rank-one matrices. There is
no interchange of noncommuting matrix factors.

The factors in (6) use disjoint physical variables. In either cut,
permuting rows and columns turns its matrix into a Kronecker product
of the corresponding factor matrices. The boundary functions have
rank `s_(i_1)` and `t_(i_k)` in both cuts; the end function is merely
transposed between the two cuts. Therefore the individual word has
exact cut ranks

```math
r_w=s_{i_1}t_{i_k}\prod_{j<k}r_{i_j i_{j+1}},
\qquad
c_w=s_{i_1}t_{i_k}\prod_{j<k}c_{i_j i_{j+1}}.       \tag{7}
```

Zero factors are included: both the Kronecker product and its
predicted rank then vanish.

## 2. Group by transition counts

Group words by their initial and final symbols and by the `M^2`
ordered transition counts `n_ij`. There are polynomially many such
types for fixed `M`, and every word in one type has the same two
ranks in (7). If a type contains `N` words, rank subadditivity gives
rank bounds `N r_w` and `N c_w` for its sum. Its rank-sum contribution
is consequently at most `N sqrt(r_w c_w)`.

Summing over types gives exactly the weighted path count in (5).
Their tensors sum to the crossing by (2).

Grouping only by total symbol counts does not justify this step:
different orders can have different row/column rank ratios.
In general `sqrt((sum r_w)(sum c_w))` is greater than the sum of
the separate geometric means. The reduction-map example below has
an exact control demonstrating the failure.

Type entries can be computed by a dynamic program recording the
last symbol, initial symbol, and transition counts while multiplying
the actual local matrices in position order. For fixed finite data
this costs polynomial time in `k`. Complete scalar tables and exact
rank factors cost `d^O(k)`. Restriction to the first `2^lambda` word
labels with `k=ceil(lambda/log_2 d)` gives the usual effective all-width
interface. No exponential list of pieces is introduced.

## 3. A finite upper certificate

Remove vertices unreachable from nonzero `s` or unable to reach
nonzero `t`. A supplied rational nonnegative matrix `H` can upper
bound `G` entrywise by checking `H_ij>=0` and
`H_ij^2>=r_ij c_ij`. Supply a positive rational vector `h` and
positive rational `beta` such that `Hh<=beta h`. Then

```math
s^TG^{k-1}t
 \le(s^Th)\max_i(t_i/h_i)\,\beta^{k-1}.           \tag{8}
```

A strict certificate with `beta<d^(sqrt(2))` is sufficient. Rational
upper matrices and supervectors can approach any strict spectral
upper bound. SCC path covers may sharpen this upper estimate, as
in the corresponding channel-walk construction, while the exact
word factorization (6) remains unchanged.

## 4. A precise common-terminal collapse

Suppose the stronger endpoint `C(q)eta=E(q)eta` holds and, for every
piece, the terminal projection is the same function:

```math
b_i(y,u)^T\eta=[y=u].                              \tag{9}
```

Choosing `y=u` in the vector endpoint gives
`sum_j a_j(x,v)=[x=v]eta`. Hence, for each fixed `i`,

```math
\sum_j T_{ij}(y,u,x',v')=[y=u][x'=v'].             \tag{10}
```

After relabeling the four legs, each transition row is already an
ordinary scalar crossing seed on `d` labels, with rank sum
`sum_j G_ij`. A nonnegative matrix has spectral radius at least its
minimum row sum, so at least one of those same-width scalar seeds
is as good as the spectral base of (5).

This remains true after restricting to reachable vertices. Every
transition to an unreachable vertex is zero, and (9) makes each
reachable state able to terminate. The surviving transition rows
still satisfy (10).

The right-eigenvector hypothesis matters here. The weaker all-word
scalar endpoint or observable-quotient certificate does not by
itself imply the vector identity used in (10). Thus common terminal
projections alone do not prove this collapse for the full interface.

The [reduction-map Fierz example](reduction-map-fierz-crossing.md)
has the stronger common-terminal property. Its nontrivial path
bound therefore illustrates the exact word construction but offers
no catalytic advantage over its scalar transition-row seeds.

## 5. Two pieces with a right-eigenvector endpoint

For two pieces, the stronger vector endpoint admits a complete
terminal-span reduction. Put `B_i(y,u)=b_i(y,u)^T eta`.
Their span cannot be zero, and it must contain the identity matrix
`D(y,u)=[y=u]`.

If the span has dimension one, every nonzero `B_i` is a nonzero
scalar multiple of `D`. Rescaling the corresponding ket and
covector inversely preserves `A_i` and reduces to Section 4.
If one terminal projection is zero, label it `B_2=0` and normalize
`B_1=D`. The vector endpoint forces `a_1(x,v)=[x=v]eta`, hence
`T_21=0`. State two cannot reach the sole nonzero terminal, so it
is removed. The remaining transition is the crossing itself, with
base `d^2`.

If `B_1,B_2` are independent, write `D=c_1 B_1+c_2 B_2`.
Comparing their coefficients in the vector endpoint forces

```math
a_i(x,v)=c_i[x=v]\eta.
```

Pieces with `c_i=0` vanish. On the remaining states,
`T_ij=c_j[x'=v']B_i(y,u)`, and each of its two cut ranks is
`d rank(B_i)`. The weighted transition matrix has identical
columns and spectral radius

```math
d\sum_{i:c_i\ne0}\operatorname{rank}(B_i)\ge d^2,
```

where the inequality follows from `rank(D)=d` and rank
subadditivity. Thus no two-piece right-eigenvector construction
has a spectral advantage beyond ordinary same-width scalar seeds.
A possible new local advantage requires at least three pieces
with a genuinely intermediate terminal-span dimension, or the
weaker observable-quotient endpoint.
