# Catalytic crossing seeds with a common boundary vector

## The weaker endpoint condition

The matrix-valued transfer does not require the local sum to equal
`E I` on the whole auxiliary space. Fix a nonzero vector `eta` and
a covector `phi` with `phi eta=1`. It is enough that every external
quadruple satisfies

```math
\left(\sum_{i=1}^M A_i(x,y,u,w)\right)\eta
 =[x=w][y=u]\eta.                                  \tag{1}
```

The same `eta` and `phi` must serve all external inputs and all
amplification lengths. The vector is an eigenvector of the sum;
individual pieces need not preserve its line.

Assume exactly the same supplied row and column factorizations and
global cross-half commutation conditions as in the
[commuting-halves transfer](commuting-halves-matrix-seed-transfer.md).
In particular every left factor entry commutes with every right
factor entry in its cut, including entries from different pieces.
Put `S=sum_i sqrt(r_i c_i)`.

**Boundary transfer.** Equation (1) constructs ordinary scalar
crossing seeds on `d^k` labels, with polynomially many type pieces
and rank sum at most `b S^k`, where `b` is the auxiliary matrix order.
All uniform construction and score conclusions of the earlier
transfer remain valid.

To prove the endpoint, set `C_j=sum_i A_i` at position `j`.
Applying the product to `eta` from the right gives

```math
C_1\cdots C_k\eta
 =E_1\cdots E_k\eta.                               \tag{2}
```

This induction uses (1) at each position; it never assumes that
individual summands preserve `eta`. Distributivity expands the left
side into the ordered piece words, and applying `phi` gives the
scalar crossing.

For a word of factor channels, the same cross-half commutation
gathers the product into `UV`. Its scalarization is

```math
\phi UV\eta=\sum_{a=1}^b(\phi U)_a(V\eta)_a.       \tag{3}
```

Thus the matrix order is paid once. Type grouping, multinomial
counts, and the proof of `b S^k` are otherwise unchanged. More
precisely, if the bilinear forms `phi UV eta` on the generated
row and column half-algebras have ranks `h_R,h_C`, the global factor
in the rank sum is `sqrt(h_R h_C)` instead of `b`.

The scalar cross-phase variant also works: its nonzero word phase
is independent of external inputs and is absorbed into one factor.
The channel-walk, SCC-path, and product-space automaton certificates
apply unchanged, with terminal functional `ell(T)=phi T eta`.
For the exact product-space automaton its terminal rank is still
the exact rank of an individual scalar channel-word function.

## 1. Finite verification and the scope of scalarization arguments

The new endpoint is a finite matrix-vector identity at every input.
Together with the supplied factorizations and commutators, it is an
exact enumerable finite-field certificate. A nonzero `eta` already
allows a normalized covector to be chosen by selecting and inverting
one nonzero coordinate. Alternatively `phi` can be supplied and
optimized for the terminal product spaces.

No new existential seed infimum is asserted: amplification still
produces ordinary finite-field seeds. The change enlarges the
algebraic construction interface.

The earlier normalized-trace and common-product-functional collapse
arguments must not be imported without another hypothesis. They
produce functionals that fix `I`. Here a replacement functional
must instead satisfy `ell(C(x,y,u,w))=E(x,y,u,w)` for the particular
endpoint family `C=sum_i A_i`. Fixing `I` alone does not imply that
condition.

## 2. An explicit Mat_4 endpoint that has no factorizing functional

The distinction is present already over `F_2` in `Mat_4`, with the
full factor pair `Mat_2 tensor I`, `I tensor Mat_2`. It does not
require a generic large-dimensional change of basis.

Write the auxiliary vector as

```math
\eta=e_0\otimes e_0+e_1\otimes e_1.
```

The space `N={T:T eta=0}` has dimension twelve. In four-by-four
matrix coordinates a basis is, for each row `a=0,...,3`,

```math
E_{a,0}+E_{a,3},\qquad E_{a,1},\qquad E_{a,2}.     \tag{4}
```

There are exactly twelve noncrossing binary external quadruples.
Assign the twelve matrices (4), in a fixed order, to those inputs.
At each of the four crossing inputs set `C=I_4`. Then

```math
C(x,y,u,w)\eta=E(x,y,u,w)\eta,                    \tag{5}
```

while `C` is not the identity-valued crossing tensor.

**Claim.** Even after extending the field, there is no linear
functional that both factorizes across this full tensor pair and
sends the entire family `C` to the scalar crossing.

Such a functional would vanish on every matrix in `N` and take
value one on `I`. Every functional vanishing on `N` has the form
`ell(T)=psi T eta`: choose any `f` with `f eta=1` and write
`T=(T eta)f+N_T`, where `N_T eta=0`. Normalization then gives
`psi eta=1`.

Under the nondegenerate trace pairing its density matrix is
`eta psi`, a rank-one matrix with image the line spanned by `eta`.
A functional that factorizes across the full tensor pair has a
product density `A tensor B`. Rank one forces both `A` and `B`
to have rank one, so its image is a product-vector line. But `eta`
has tensor rank two over every field extension. This contradiction
proves the claim.

This example rules out that particular one-step scalarization
method. It does not certify an improved rank sum. The extended-cut
ranks of the single piece `C` are `(11,11)` over `F_2`.

For a control in which the pieces do not individually preserve the
boundary, take `Q=E_(1,0)` and the two pieces `A_1=Q`, `A_2=C+Q`.
Then `Q eta=e_1` is not a multiple of `eta`, but their sum satisfies
(5). Both cuts use the same full commuting tensor halves.
Their exact extended-cut profiles are `(1,1)` and `(12,12)`.
Thus this is a valid, noncompetitive two-piece certificate with
rank sum thirteen. Its ordered word products are noncommutative.

## 3. A concrete CAR pair of factor structures

There is an explicit source of different commuting half structures
for this boundary interface. In characteristic two let `a_i,c_i`
be the annihilation and creation operators on the exterior Fock
space `F` with `m` modes. Write `n=dim F=2^m` and use
`End(F)` as the auxiliary vector space, of dimension `b=n^2`.
Let `L_A(X)=AX` and `R_A(X)=XA`; left and right multiplications
commute. Take `eta=I_F`.

Define

```math
a'_i=L_{a_i}+R_{a_i},\qquad c'_i=L_{c_i},
\qquad
a''_i=L_{c_i}+R_{c_i},\qquad c''_i=R_{a_i}.        \tag{6}
```

Each primed family is a CAR representation, and every primed
generator commutes with every double-primed generator. For example,
`[a'_i,a''_j]=delta_(ij)(L_I+R_I)=0`, while
`[a'_i,c'_j]=[a''_i,c''_j]=delta_(ij)I`. All generators square
to zero, and different modes have the required zero commutators.
Both annihilator families kill the same vector:

```math
a'_i eta=0,\qquad a''_i eta=0.                    \tag{7}
```

The vectors `prod_i(c'_i)^(p_i) prod_i(c''_i)^(q_i) eta`
are the normally ordered Clifford monomials `c_P a_Q`, which form
a basis of `End(F)`: in one mode `I,c,a,ca` are a basis of `Mat_2`,
and tensoring these bases gives the assertion for every number of
modes. Hence (6) gives two full commuting matrix
halves, each isomorphic to `Mat_n`, in a new tensor-factor structure.
In that structure `eta` is the product vacuum. In the original
left/right structure it is the maximally entangled identity vector.

The covector `phi(X)=<0|X|0>` satisfies `phi eta=1`. Its bilinear
form on the original left/right halves is `<0|AB|0>`, of rank `n`.
On the two new halves it factorizes through the product vacuum and
has rank one: `phi(c_P a_Q)` is zero unless both subsets are empty,
so `phi` is also the product vacuum covector in this basis.
This gives exact asymmetric terminal ranks for any
candidate that supplies factorizations in these two structures.

The CAR additive commutators do **not** justify gathering an
arbitrary annihilator product past an arbitrary creator product.
A useful Fierz decomposition must place its supplied factors into
the actual commuting halves above, or separately account for its
additional contractions. Equations (6)--(7) supply a concrete
algebraic interface; no improving crossing profile is claimed here.

## 4. A finite observable-quotient endpoint certificate

There is a further finite relaxation. With `C(q)=sum_i A_i(q)`,
the amplification proof only needs

```math
\phi C(q_1)\cdots C(q_k)\eta
   =\prod_{j=1}^k E(q_j)                            \tag{8}
```

for all words of external quadruples. One sufficient finite
certificate supplies an auxiliary subspace `R` containing `eta`,
invariant under every `C(q)`, such that `phi eta=1` and

```math
\phi C(q)|_R=E(q)\phi|_R.                         \tag{9}
```

All these conditions are checked on a finite basis. Applying (9)
from the left proves (8), since each remaining suffix sends `eta`
into `R`. No individual piece or factor half is required to
preserve `R`; gathering still takes place in the original ambient
matrix algebra.

There is also a direct necessary-and-sufficient finite test for
(8). On a space of dimension `b+1`, use the block diagonal matrices
`diag(C(q),E(q))` and initial vector `(eta,1)`. Generate their
reachable linear span by repeatedly applying each matrix to a
basis and adjoining independent vectors. The dimension increases
at most `b+1` times. Equation (8), including its empty-word
normalization, holds exactly when `(phi,-1)` annihilates that span.
Every reached vector is a linear combination of actual word
vectors, and closure contains every such vector, proving both
directions of the test.

The same transfer theorem, scalar-twist variants, and terminal
automata therefore accept this finite endpoint certificate too.
A scalar identity at length one alone is insufficient: off-diagonal
auxiliary transitions can leave the boundary and later return with
a nonzero contribution. The checker includes this negative control.

## Verification

`check_boundary_stable_catalytic.py` verifies the finite Mat_4 family,
both exact extended-cut rank factorizations, global cross-half
commutation, ordered type amplification through three copies, and
the failure of all 64 normalized product density matrices over
`F_2` to preserve its endpoint family. It also checks the CAR
relations, the commuting new halves, their common vacuum, and the
two terminal bilinear ranks at one and two modes. Evidence is
recorded in `boundary-stable-catalytic.json`.
