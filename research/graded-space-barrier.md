# A rank restriction on slower graded-space substitutions

## Result and scope

The Appendix B crossing construction uses partition spaces with logarithmic
dimension proportional to the square root of the energy. A direct replacement
by more slowly growing spaces cannot preserve all of its interfaces.

This note proves the obstruction even when the replacement operators differ
from the original operators. The assumptions are a polynomial number of
pieces, nearly full use of the input space, and the same intermediate-energy
range. It does not prove a factoring lower bound or rule out different
crossing constructions.

The reference is Claude's Appendix B in
[`anthropics/formal-math`, revision `16a5ce4`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
Its factorization formulas and declared ranks appear in
[AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
Its assembled algebraic theorem is
[Ring/ThmHalfPieces.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Ring/ThmHalfPieces.lean).

## 1. An exact rank budget

Let `X` be a set of `m` labels. Suppose that rational-valued tensors satisfy

```math
\sum_{i=1}^{M} A_i(x,y,u,v)
=d\,[x=v][y=u],\qquad d\in\mathbb Q\setminus\{0\}.
\tag{1}
```

Suppose each tensor has a row factorization with `r_i` channels:

```math
A_i(x,y,u,v)=\sum_{j=1}^{r_i}
F_{i,j}(x,y)G_{i,j}(u,v).
\tag{2}
```

Then necessarily

```math
m^2\le\sum_{i=1}^{M}r_i.
\tag{3}
```

**Proof.** Flatten a tensor into a matrix whose rows are `(x,y)` and whose
columns are `(u,v)`. The right side of (1) is a nonzero scalar times a
permutation matrix of size `m^2`. Its rank over the rationals is `m^2`.
Equation (2) is a sum of `r_i` matrices of rank at most one. Subadditivity
of matrix rank now gives (3).

No assumption says that a factorization has minimum rank. Upper bounds on
its channel count are enough for the obstruction. Only one flattening is
needed. The common scale must remain nonzero in the field under discussion.

## 2. A general restriction on graded spaces

Let `V_0,V_1,...` be finite-dimensional rational vector spaces. Put

```math
D(H)=\sum_{E=0}^{\lfloor H\rfloor}\dim V_E.
```

Fix constants `L>0`, `C>0`, and `gamma>0`. Suppose a family of decompositions,
indexed by arbitrarily large integers `h`, has these properties:

1. The number of labels satisfies
   `log m(h)=C*h^gamma+o(h^gamma)`.
2. The cumulative dimensions satisfy
   `log D(H)=C*H^gamma+o(H^gamma)` as `H` tends to infinity.
3. There are at most `h^A` pieces, for a fixed constant `A`.
4. Each row factorization uses at most `h^B*D(L*h)` channels, for a
   fixed constant `B`.
5. The pieces sum to a nonzero scalar multiple of the crossing tensor.

Natural logarithms are used in these asymptotic assumptions. Condition 1
means that the labels use the input dimension to first order. It is an
essential assumption, not a consequence of merely choosing some labels.

**Theorem.** These conditions imply

```math
L^\gamma\ge2.
\tag{4}
```

**Proof.** The exact rank budget gives

```math
m(h)^2\le h^{A+B}D(Lh).
```

Take logarithms, divide by `C*h^gamma`, and let `h` tend to infinity.
The polynomial term contributes `O(log h)/h^gamma`, which tends to zero.
The two sides tend to `2` and `L^gamma`, respectively. This proves (4).

The same argument gives a useful form without a power-law assumption:

```math
\log D(Lh)\ge2\log m(h)-O(\log h).
\tag{5}
```

For binary labels, `m(h)=2^lambda`. Under Condition 1,
`lambda=Theta(h^gamma)` and `h=Theta(lambda^(1/gamma))`.
Thus `h` is polynomial in `lambda`, and `log h=o(lambda)`.
The piece counts in the theorem are compatible with the source's polynomial
count requirement. An exponential number of pieces in `lambda` can avoid
this rank argument's asymptotic conclusion, but violates that requirement.

## 3. The Appendix B energy range forces gamma at least one half

In Appendix B, the labels name partitions of `h`. Its pieces have indices
`s=-2h,...,2h`. The row factorization sums over

```math
i=0,\ldots,2h,\qquad \sigma\in V_{s+i}.
```

Negative energies contribute no states. Every row channel therefore lies
in a space of energy at most `4h`, and there are polynomially many repeated
copies. The column factorization has the same energy bound.

If a replacement keeps these ranges and has the first-order dimension law
in Section 2, then `L=4`. Equation (4) gives

```math
4^\gamma\ge2,
\qquad\text{hence}\qquad \gamma\ge\tfrac12.
\tag{6}
```

Consequently, a proposed replacement with
`log D(H)~C*H^(1/3)` cannot have all these properties. The failure occurs
before any analysis of the factoring program's runtime.

For the original partition construction, ordinary partition estimates give
`log p(h)=C*sqrt(h)+O(log h)`. The cumulative dimension has the same leading
term. The smallest `h` with `p(h)>=2^lambda` uses that leading dimension,
so the original construction meets the boundary `gamma=1/2`.
The relevant source bounds are in
[ClauseBound/Part/Bounds.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Part/Bounds.lean)
and [GeneralPieces/Ring/Bits.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Ring/Bits.lean).

This does not prove that the source's listed channel counts equal the
actual ranks. Rank compression inside those lists remains a different route.

## 4. Consequence for the dimension-only rate calculation

Consider the same profile as Appendix B. For a piece with `s=t*h`, the
largest row and column energies are `(2+t)h` and `(2-t)h`, respectively.
For `0<gamma<1`, dimension counting gives the leading width sum

```math
\bigl((2+t)^\gamma+(2-t)^\gamma\bigr)\lambda
\le 2^{\gamma+1}\lambda,
\qquad -2\le t\le2.
\tag{7}
```

Concavity puts the maximum at `t=0`. Polynomial multiplicities contribute
only `O(log h)=o(lambda)` additional bits.
The leading contraction coefficient from this calculation is

```math
\rho=2^{\gamma-1},\qquad
\alpha=\log_2(1/\rho)=1-\gamma.
```

The rank restriction (6) therefore makes `alpha<=1/2` for this substitution
method. Slower dimension growth cannot improve the coefficient while all
the other interfaces remain as stated.

The asymptotic assumptions above control a leading coefficient only.
An exact endpoint under Appendix A also requires its specified remainder
bound, coefficient sizes, and uniform evaluation costs. This note does not
infer those properties from a dimension asymptotic.

Increasing the energy range alone also needs care. For a symmetric profile
whose extreme energy is `Lh`, the corresponding central energy is `Lh/2`.
The dimension-only coefficient becomes

```math
\rho=\frac{(L/2)^\gamma}{2}.
```

Equation (4) forces `rho>=2^(-gamma)`. If `gamma<1/2`, this is larger than
`1/sqrt(2)`. Thus slower growth plus a proportionally enlarged symmetric
energy range still does not yield a better leading coefficient. This last
conclusion assumes that symmetric profile; it does not cover arbitrary
intermediate spaces or asymmetric factorizations.

## 5. Exact finite checks

Run:

```sh
python3 tools/check_graded_space.py
```

The checker computes the crossing flattening's rank by rational elimination
for small label sets and several nonzero scales. A zero-scale negative
control confirms why the nonzero hypothesis matters.
It also checks finite polynomial-channel budgets and exact rational tests
for the necessary condition (4).

For integer `L` and positive rational `gamma=a/b`, the test is

```math
L^a\ge2^b.
```

For example, the proposed `gamma=1/3` with range `4h` fails because
`4<8`. The original `gamma=1/2` meets the boundary because `4=4`.
These integer tests do not replace the all-`h` proof in Section 2.

The output is recorded in
[graded-space.json](../evidence/graded-space.json).

## What remains open

A better rate can use smaller actual ranks, a different energy profile,
or a decomposition outside this graded-space form. It must still satisfy
the full effective-pieces interface. This note gives a quick rejection
criterion for a broad class of apparent improvements; it supplies no
new family with `rho<1/sqrt(2)`.

The theorem and exact checks used substantial OpenAI assistant support.
They are ordinary mathematics and computation, with no new Lean proof.
