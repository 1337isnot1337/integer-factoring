# Every primitive recipient set has maximal leading rank entropy

## Result and scope

The sharp surjectivity threshold in
[factoring-primitive-recipient-codec.md](factoring-primitive-recipient-codec.md)
has a complementary consequence above the threshold. For **every** set
of odd recipients, its joint message rank differs from the largest rank
allowed by the input and output cardinalities by only `O(h log h)` bits.
This includes arbitrary, irregular recipient sets of every size.

The proof uses the same negacyclic symmetry, now to bound every nonzero
kernel vector from below. It then bounds each fiber of the digit box by
Euclidean packing. This is an unconditional local rank theorem for the
existing arithmetic network. It is not a new factoring algorithm or a
lower bound on a fully contracted, fixed-product network.

Use the pinned definitions and notation in
[factoring-network-structure.md](factoring-network-structure.md):

```math
h=2^s,\quad L=2h^2,\quad B=2^L,\quad M=B^4+1,
\quad A=2h,\quad \zeta=2^{4L/h}\pmod M.
\tag{1}
```

For `J` a set of `k` distinct odd indices in `{0,...,A-1}`, define

```math
F_J(u)=\left(\sum_{i=0}^{h-1}u_i\zeta^{ij}\right)_{j\in J},
\qquad u\in\{0,\ldots,B-1\}^h.
\tag{2}
```

Gate-dependent invertible coordinate multipliers do not affect its image
size. Let `E_J` be the equality matrix between these digit settings and
the unrestricted requested message tuple. Every row is a standard basis
vector, so its rank over **every** coefficient field is `|image F_J|`.

Theorem. For `0<=k<=h`,

```math
\boxed{\quad
\frac{\min(B^h,M^k)}{(1+\sqrt h)^h}
\ \le\ \operatorname{rank}E_J
\ \le\ \min(B^h,M^k).
\quad}                                                   \tag{3}
```

In particular, writing `D_h=h log2(1+sqrt h)`,

```math
\min(hL,k\log_2M)-D_h
\le\log_2\operatorname{rank}E_J
\le\min(hL,k\log_2M).
\tag{4}
```

For the source relation `L=2h^2`, the additive error is `O(h log h)`,
whereas the full input entropy is `2h^3`. This additive estimate is
uniform in the recipient set and in its size. It complements the exact
rank `M^k` for `k<h/4` proved in the companion note.

## 1. A determinant bound for every nonzero kernel vector

For `1<=k<=h`, let

```math
\Lambda_J=\{v\in\mathbb Z^h:F_J(v)=0\pmod M\}.
\tag{5}
```

The composite-ring Vandermonde argument in the companion notes proves
that unrestricted integer evaluation onto `(Z/MZ)^k` is surjective.
Consequently `det Lambda_J=M^k`. This does not require `M` to be prime.

The negacyclic shift

```math
S(v_0,\ldots,v_{h-1})=(-v_{h-1},v_0,\ldots,v_{h-2})
```

preserves `Lambda_J`, because every recipient is odd and
`(zeta^j)^h=-1`. For any nonzero `v` in the lattice, form the integer
matrix

```math
W=(v,Sv,\ldots,S^{h-1}v).
\tag{6}
```

The irreducibility of `X^h+1` over `Q` makes `W` nonsingular, as proved
in Section 2 of the primitive-recipient note. Its columns generate a
full-rank sublattice of `Lambda_J`. Hence

```math
M^k\mid\det W,\qquad |\det W|\ge M^k.
\tag{7}
```

All columns of `W` have Euclidean norm `||v||_2`. Hadamard's determinant
inequality therefore gives

```math
M^k\le|\det W|\le\|v\|_2^h,
\qquad\boxed{\|v\|_2\ge M^{k/h}}.
\tag{8}
```

This bound holds for **every** nonzero kernel vector. It uses no
shortest-vector algorithm and no assumption about the distribution of
the evaluation points beyond distinctness and odd indices. In
particular, the lattice cannot have an exceptionally short direction
that creates large fibers in the input box.

## 2. Packing the fibers with exact volume constants

Put `delta=M^(k/h)` and fix an output tuple. Any two distinct input
vectors in its fiber differ by a nonzero vector of `Lambda_J`, so their
Euclidean distance is at least `delta`.

The input cube lies in the Euclidean ball with center
`c=((B-1)/2,...,(B-1)/2)` and radius `(B-1)sqrt(h)/2`. Open balls of
radius `delta/2` around the fiber points are disjoint. Their union lies
in the ball with center `c` and radius
`((B-1)sqrt(h)+delta)/2`. Comparing volumes cancels the dimension's
unit-ball volume and both factors of two. Thus

```math
|F_J^{-1}(w)|
\le\left(1+\frac{(B-1)\sqrt h}{\delta}\right)^h
\le\left(1+\frac{B\sqrt h}{\delta}\right)^h.
\tag{9}
```

Open balls avoid any issue when two centers are exactly `delta` apart;
their boundaries have zero volume. Dividing the number `B^h` of inputs
by the maximum fiber size gives the stronger unsimplified estimate

```math
\operatorname{rank}E_J
\ge\frac{B^h}{(1+B\sqrt h/M^{k/h})^h}.
\tag{10}
```

If `delta<=B`, the denominator in (10) is at most
`((1+sqrt h)B/delta)^h`, giving `M^k/(1+sqrt h)^h`.
If `delta>=B`, it is at most `(1+sqrt h)^h`, giving
`B^h/(1+sqrt h)^h`. These are exactly the lower bounds in (3).
The upper bound is immediate from the two finite cardinalities. For
`k=0`, the image has one element and (3) holds directly.

## 3. Exact injection at the threshold would be false

The estimate is intentionally a rank envelope, rather than a claim of
exact rank `min(B^h,M^k)`. There are actual source collisions at the
critical value `k=h/4`.

Take `h=8`, `L=128`, `D=2^64`, so `B=D^2`, `M=D^8+1`, and `zeta=D`.
For `J={1,7}`, the polynomial

```math
v(X)=-D-(B-1)X+DX^2=(X-D)(DX+1)                \tag{11}
```

vanishes at `D` and at `D^7=-D^{-1}` modulo `M`. Its coefficient vector
is a nonzero kernel vector with every coefficient between `-(B-1)` and
`B-1`. Therefore the distinct valid digit inputs

```math
(0,0,D,0,0,0,0,0),\qquad(D,B-1,0,0,0,0,0,0)
\tag{12}
```

have the same pair of messages. Their difference also illustrates the
strength of (8): its squared Euclidean norm is `B^2+1`, and its orbit
determinant has absolute value `M^2`. There is no contradiction with
the near-maximal leading entropy in (3).

## 4. A tradeoff for independently factored joint layer interfaces

The rank estimate extends one specific construction obstruction from
[half-top-structure.md](half-top-structure.md). Unlike that note, the
original wires from the unrestricted `faQ` layer may now be bundled
jointly between each pair of groups, with arbitrary recipient sets.
The following assumptions specify exactly which construction is covered.

1. Partition the original gates into `P>=2` groups, and visit every
   original local setting in each group at least once.
2. For each unordered pair of groups, factor jointly all equality
   constraints from the `faQ`--`mid` layer between that pair. These
   declared layer interfaces are factored independently of the other
   network constraints and concatenated into the clique channels.
3. Keep the dense clique interface, the explicit pad, and the level
   schedule used in the completed half-rate composition.

This allows joint channel compression within a group pair, arbitrary
rational coefficients, and arbitrary assignment of recipients to groups.
It does not cover factoring a whole contracted group tensor, using the
target-product or carry constraints to change the boundary relation, or
evaluating local sums without visiting their Cartesian setting domain.

Let `q_u` count `faQ` gates in group `u`, and let `o_u` count the odd
`mid` recipient indices in it. Then

```math
\sum_u q_u=A,\qquad\sum_u o_u=h,
\qquad q_{\max}=\max_u q_u,
\qquad H=hLq_{\max}.
\tag{13}
```

A group attaining `q_max` has at least `2^H` original local settings.
This lower bound concerns the specified enumeration procedure, even if
many settings would later be rejected by other constraints.

For one `faQ` gate in group `u`, its messages sent to `mid` gates in
group `v` project onto the `o_v` odd coordinates. By (4), the associated
equality matrix has logarithmic rank at least

```math
L\min(4o_v,h)-D_h.                            \tag{14}
```

Different `faQ` gates have independent digit settings. The requested
coordinates of different incoming vectors at `mid` gates are also
independent in this declared interface. Its matrix is therefore, up to
row and column permutations, a tensor product of the individual
fan-out equality matrices. The reverse direction between the two groups
uses disjoint setting coordinates and contributes the transpose of
another such tensor product. Ranks multiply for these factors.

If `C` is the sum over unordered group pairs of the logarithmic ranks
of their joint layer interfaces, (14) consequently implies

```math
C\ge L\sum_u q_u\sum_{v\ne u}\min(4o_v,h)
       -A(P-1)D_h.                            \tag{15}
```

The error term is conservative: pairs with no relevant odd coordinates
need no error. We retain the simpler uniform expression.

At most one group can have `o_u>3h/4`. For every other group,

```math
\sum_{v\ne u}\min(4o_v,h)
\ge\min\left(4\sum_{v\ne u}o_v,h\right)=h.
\tag{16}
```

If there is one exceptional group, its number of `faQ` gates is at most
`q_max`. If there is none, every group contributes at least `h`.
Combining (13)--(16) gives

```math
\boxed{H+C\ge AhL-A(P-1)D_h
       =L^2\left(1-\frac{(P-1)\log_2(1+\sqrt h)}L\right).}
\tag{17}
```

For example, if `h>=2` and `P<=h`, then
`log2(1+sqrt h)<=log2(h+1)<=h`, and `L=2h^2`. Thus the factor
in parentheses is greater than one half. If `n` is the source input bit
length, its capacity condition `n<=L^2` yields

```math
H+C\ge n/2.                                  \tag{18}
```

The standard source cutoff `P^40<=n<=4h^4` implies `P<=h` for `P>=2`,
so (18) holds throughout that regime. More generally, the error in
(17) is `O(P h^2 log h)`, which is `o(L^2)` whenever
`P log h=o(h^2)`.

## 5. Consequence for the unchanged explicit pad

If a common clique channel width `lambda` accommodates the independent
pair interfaces, then

```math
C\le\binom P2\lambda\le V/2,
\qquad V=P^2\lambda.
\tag{19}
```

Use the unchanged schedule `P=K_ell`,
`r=2^floor(ell/2)`, and the explicit pad exponent `q_pad>V/r` from
[half-top-structure.md](half-top-structure.md), Section 4. The complete
cost `T` includes visiting local settings and writing that pad. In the
source cutoff regime,

```math
\log_2T\ge\max(H,V/r)
\ge\frac{H+V/2}{1+r/2}
\ge\frac{n}{r+2}.                             \tag{20}
```

The schedule obeys `log2 P>=x0 r^2`, with fixed `x0>=18`. Since the
source cutoff also gives `P<=n`, equation (20) implies

```math
\log_2T\ge\frac{n}{2+\sqrt{\log_2 n/x_0}}
=\Omega\!\left(\frac n{\sqrt{\log n}}\right).
\tag{21}
```

This is a restricted obstruction for the declared independently
factored layer interfaces and explicit pad. It says nothing about the
necessity of that pad for a different evaluator. The positive next
targets remain a representation that preserves correlations across
interfaces, a different arithmetic network, or a cheaper evaluator on
a family that remains closed under every recursive step.

## 6. Reproducible finite checks

Run:

```sh
python3 tools/check_factoring_recipient_rank_envelope.py \
  --output evidence/factoring-recipient-rank-envelope.json
```

The checker exhaustively enumerates bounded-digit images for small
generalized parameter sets, checks exact orbit determinant divisibility
and Hadamard inequalities, verifies the actual source collision (12),
and exhaustively checks the integer grouping inequality underlying
(17). It uses integer bounds for the irrational-valued rank estimate.
The evidence labels reduced-width examples separately from the actual
source relation `L=2h^2`. Finite checks are controls for the proof, not a
substitute for its uniform claims or evidence of a score increase.
