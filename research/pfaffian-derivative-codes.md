# Self-dual Pfaffian codes and uniform equal-copy crossing blocks

## Result and scope

Pfaffian derivative spaces give explicit self-dual physical crossing codes
over F2. With identical edge-mode crossings, their full natural cut matrix
splits into O(n^2) representative blocks. A cancellation identity reduces
every representative to one of two explicit binary Hamming kernels.
The tensor itself has a uniform `2^O(n)` construction, even though listing
all matching monomials individually can cost more.

The [registered diagnostic](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19339013)
computed every block at n=2,4,6. Both complete natural profiles are
established by the compression identity below; direct literal tensors
independently check both cuts at n=2,4. Sampled pressure values are
0.76023,0.82640,0.83945, respectively, all above the target. No global
scale optimum, all-scale floor or improved factoring exponent is asserted.
The complete score remains **1/2**.

This changes both the determinant code's bipartite matching structure and
its distinct-height reflection kernel. It is not a boundary retraction.
The Gaussian floor for a full ambient alphabet does not automatically
apply to the compressed physical code. Further work needs a uniform rank
mechanism from the formulas, not a blind increase of the finite n panel.

## 1. Physical codes with an exact self-pairing

Let n>=2 be even. Use one ambient binary occupation mode x_uv for every
edge `{u,v}` of the complete graph on n vertices. For every even vertex
subset S, let p_S be the Pfaffian on those vertices: the sum of the
edge-occupation monomials of all perfect matchings of S. Set p_empty=1.

Deleting an edge uv sends p_S to p_(S minus{u,v}) if u,v belong to S,
and to zero otherwise. The span is therefore closed under every bit
deletion. Different S have disjoint monomial supports. The number of
matchings on a nonempty S is `(|S|-1)!!`, an odd integer, so the ordinary
bilinear pairing over F2 gives

~~~
<p_S,p_T>=[S=T],             P^T P=I.
~~~

Both the embedding P and retraction P^T are deletion closed. The physical
alphabet has `q=2^(n-1)` labels and each label has occupation and height
`h_S=|S|/2`. Its height inventory is `sum_r binomial(n,2r)t^r`.
These statements have direct proofs for every n; the finite checks do
not supply their quantifiers. This is the Pfaffian counterpart of the
classical derivative/minor spaces, with a nondegenerate characteristic-two
self-pairing. The [determinant construction](determinant-inverse-codes.md)
has a different self-pairing because factorial counts are even.

## 2. An actual strict polynomial crossing, with all acquisition costs

For each edge mode, use one independent four-mode Gaussian crossing on
physical groups X,Y,U,V, with homogeneous alternating step coefficients

~~~
k_01=y²+xy,       k_10=x²+xy,       k_11=xy.
~~~

At x=y=1 its kernel is the matching X--V,Y--U. Its principal Pfaffians
therefore give the binary crossing. Take the tensor product over all
edge modes, and contract every physical leg with P. The identity P^T P=I
makes its endpoint exactly `[X=V][Y=U]` on the new labels. Every entry
is homogeneous of degree `h_X+h_Y+h_U+h_V`, at most 2n. All coefficients
are binary and no division or pole is used.

There is a direct `2^O(n)*poly(n)` algorithm to acquire the entire tensor.
Use a square-zero algebra on 4n vertex variables z_(leg,u). For each edge
uv, multiply its eight-term local Gaussian state, replacing an occupied
leg by `z_(leg,u) z_(leg,v)`. A state is a subset of the 4n variables; overlapping
vertex bits make its product zero. The coefficient of the four prescribed
vertex subsets is the required tensor entry. The algebra has at most
`2^(4n)` coefficient cells, and every polynomial coefficient has degree
at most 2n. Processing O(n²) edge factors therefore has the claimed time
and space bound. This avoids factorial matching enumeration, including
when the physical bit length is `log2 q=n-1`.

The [implementation](../tools/pfaffian_derivative_codes.py) includes this
subset construction and a separate literal matching expansion for small
controls. Sparse basis enumeration is used only where its finite cost is
explicit. The construction bound alone is not a favorable rank bound or
a stronger factoring algorithm.

## 3. Exact regular cut matrices

Set y=1,x=t. Removing the two within-cut Gaussian factors gives the
one-particle crossblock

~~~
B=[[d I,j I],[j I,d I]],       d=t+t²,       j=t.
~~~

Here is why this removal survives physical compression. A within-cut
factor is the invertible creation operator
`product_e (I+(1+t) C_e tensor C_e)`. Since the physical map space is
deletion closed, its row space is invariant under the adjoint creation
operators: `P^T C_e=(D_e P)^T` is a linear combination of rows of P^T.
Each induced factor is invertible, as its nilpotent term squares to zero.
The corresponding column factors are transposes and are also invertible.
They are polynomial units at the local place and do not change Smith
orders. The remaining matrix is the restriction of the full compound of
B by `P tensor P`, with a harmless physical column-pair swap.

For the opposite cut, set x=1,y=t. The same seed interchanges the two
within/cross steps, and the regular crossblock is again B. The same
self-dual maps occur on both sides. Thus **both natural cut profiles are
the same profile of this compressed compound**. This is an identity of
the actual tensor construction, not an assignment-dependent choice of
parent labeling or a relabeling of coefficient summands.

After reordering by edge, B has independent two-by-two blocks. A selected
minor is zero unless each edge has equal row/column occupation counts.
Otherwise each singly occupied edge contributes d or j, and each doubly
occupied edge contributes `d²+j²=t⁴`. This gives an independent direct
matching formula for every compressed entry.

## 4. Vertex-degree blocks and their multiplicities

For a paired physical label `(S_X,S_Y)`, let
`I=S_X intersection S_Y` and `D=S_X symmetric-difference S_Y`.
Their vertices have combined degree two and one, respectively. Equality
of edge-occupation counts in a nonzero minor preserves this vertex degree
vector. Thus the matrix splits by disjoint I,D.

Write a=|I| and delta=|D|, necessarily even delta. Inside a block the
basis is `S_X=I union A, S_Y=I union(D minus A)`, where A ranges over
subsets of D with `|A|=a mod2`. Its size is `2^(delta-1)` for delta>0.
For delta=0 it has size one if a is even, and is absent otherwise.
All edge modes use identical coefficients, so vertex permutations
identify every block with the same a,delta. Its multiplicity is

~~~
binomial(n,a)*binomial(n-a,delta).
~~~

The multiplicities times block sizes sum to q². Every block has a
permutation endpoint, full generic rank and polynomial entries of degree
at most `2a+delta`. Complete Smith computation can therefore use its
finite determinant-degree bound without interpreting a missing channel
as a zero channel.

## 5. Two closed block formulas, valid for every n

Call the representative matrix M_(a,delta), with the physical column-pair
swap included. Transpose two vertices of I in the matching expansion.
It pairs and cancels all nonfixed terms in characteristic two. A perfect
matching containing both transposed vertices is fixed by that transposition
only if it pairs them to each other. Thus each of the four matchings in
a fixed term contains their common edge. That edge contributes the
doubly occupied minor t⁴; removing it gives exactly M_(a-2,delta). Hence

~~~
M_(a,delta)=t⁴ M_(a-2,delta),       a>=2.             (1)
~~~

For a=0, every vertex has degree one. A common matching can only pair
vertices within the four cells of the row/column side labels. Their
cardinalities must all be even. Each even cell has an odd number of
matchings, and their parity condition is `|A intersection B|` even.
For even subsets A,B of D this gives

~~~
M_(0,delta)[A,B] = 0                        if |A intersection B| odd,
                  t^(delta/2)(1+t)^(|A xor B|/2)  otherwise.          (2)
~~~

For a=1 there is one internal vertex and a two-edge path through it.
If the intersection is odd, the path endpoints lie in the two diagonal
odd cells; if it is even, they lie in the two off-diagonal odd cells.
Endpoint choices have odd count. The remaining cells have odd matching
counts as before. The path contributes, after accounting for the removed
endpoints, an extra j in the former case and an extra d in the latter.
For odd subsets A,B the result is

~~~
M_(1,delta)[A,B]
 =t^(delta/2+1)(1+t)^(|A xor B|/2 + [|A intersection B| even]).        (3)
~~~

Equations (1)--(3), plus `M_(a,0)=[t^(2a)]` for even a, are the promised
uniform block construction. No fitted pattern from small n is used.

The leading matrix in (2), after removing t^(delta/2), is
`1+A dot B` on the even-weight hyperplane. Its bilinear form has radical
spanned by the all-one vector, so rank delta-2; the independent constant
adds one. Its exact leading rank is therefore delta-1, compared with
full size `2^(delta-1)`. This is a uniform local rank fact, not by itself
a favorable complete pressure estimate.

## 6. Complete finite controls and the registered stopping point

Run `make verify-pfaffian-derivative-codes`. The
[checker](../tools/check_pfaffian_derivative_codes.py) and
[receipt](../evidence/pfaffian-derivative-code-controls.json) retain:

- 1,092 orthogonal-pairing entries and 530 edge deletions;
- 1,240 direct matching entries compared with the closed formulas;
- Every complete representative Smith profile, with 108 independent
  Toeplitz hinge checks and three leading-rank controls;
- All 4,112 small tensor entries from both the literal matching and uniform
  subset constructions, their endpoint/degree checks and four complete
  literal profiles for both cuts at n=2,4.

At n=2,4,6 the exact multiplicities reconstruct profiles of length 4,64,1024.
The shared opposite profile is supplied by the proof in Section 3; only
the stated small cases have a second literal full-matrix reconstruction.
The bounded numerical scale fits are:

| n | q | Sampled natural rho | Scale |
| --- | --- | --- | --- |
| 2 | 2 | 0.7602278113 | 0.9380072 |
| 4 | 8 | 0.8263958849 | 1.1860630 |
| 6 | 32 | 0.8394548368 | 1.0177383 |

All exceed 1/sqrt2. These values do not certify global scale optima or
exclude different weights, kernels or coefficient-rank methods. No strict
witness triggers rational type/log certification. The registered finite
panel is closed; any continuation must derive new uniform rank structure
from the two kernels, rather than enlarging n merely to fit another
decimal. The construction and block identities are retained as partial
research results, not a complete factoring improvement.

Credit: classical Pfaffian derivative spaces and Gaussian cut
factorization; the existing coherent deletion-space and typed-word
interfaces. This application and its code/proofs were developed with
substantial OpenAI assistant support. No new Lean or independent
mathematical acceptance is claimed.
