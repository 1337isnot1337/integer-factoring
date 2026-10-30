# The full CNOT tensor span: exact ranks and rank-one structure

The four coefficient matrices of the three-CNOT path admit a small exact
rank representation at every tensor power. This note gives that
representation, classifies all tensors of rank one in either cut over
every characteristic-two field, and uses the classification to exclude
the full eight-label coefficient tensor space. The exclusion allows
arbitrary nonsymmetric recombinations and extension fields.

The note also gives an explicit rank-one conic and exact ranks for all
coefficients of its tensor powers. These formulas provide structured
candidate pieces for larger searches. None of these statements improves
the integer-factoring exponent. The full coefficient tensor spaces at
larger powers remain open here.

## 1. Conventions and the local spaces

Let `C0,C1,C2,C3` be the coefficient matrices in the companion
[three-CNOT construction](cnot-crossing-modules.md). Their row matrices, stored as packed binary
rows in order `00,01,10,11`, are

```text
C0 = [1,2,4,8]
C1 = [0,0,12,12]
C2 = [0,12,10,6]
C3 = [0,10,0,10].
```

The row cut has indices `(x,y)|(u,v)`. The column cut has indices
`(x,u)|(y,v)`. Exchanging the two column indices preserves every rank,
but the literal basis matrices below use the stated order.

Use the coefficient basis

\[
a=C_0+C_2,\qquad b=C_2,\qquad c=C_1,\qquad d=C_3.
\]

The crossing is `E=a+c+d`. The row basis with columns
`[1,14,12,10]` identifies the coefficient representation with

\[
 a\longmapsto I_2\oplus0,\quad
 b\longmapsto0\oplus I_2,\quad
 c\longmapsto0\oplus E_{01},\quad
 d\longmapsto0\oplus E_{10}.                              \tag{1}
\]

The coefficient letters and their scalar coordinates are distinguished
by context. A tensor with coordinates `(a,b,c,d)` has column matrix,
in the basis with columns `[6,9,1,2]`,

\[
\begin{pmatrix}
p&p+d&0&p\\
0&d&a+c&c+p\\
0&0&c&c+p\\
0&0&0&p
\end{pmatrix},\qquad p=a+b.                              \tag{2}
\]

All changes of basis are defined over `F2`, so the identities hold after
extension to any field of characteristic two.

## 2. Exact full-span row rank

Write `V=span{a,b,c,d}` and let `T` be any element of `V^tensor(t)`.
For `S` a subset of the tensor sites, define `T_S` by extracting the
coefficient `a` at every site outside `S` and substituting
`b,c,d -> I2,E01,E10` at the sites in `S`. A coefficient word with an
`a` in `S`, or a `b,c,d` outside `S`, contributes zero.

Tensoring the direct sum (1) gives the exact formula

\[
\boxed{\quad r(T)=\sum_{S\subseteq[t]}2^{t-|S|}
                    \operatorname{rank}(T_S).\quad}       \tag{3}
\]

The largest block has side length `2^t`, rather than `4^t`. The total
number of entries in all distinct blocks is
`sum_S 4^|S|=5^t`. The multiplicities in (3) are ordinary integers and
are not reduced modulo the field characteristic.

The three column characters in (2) are `p=a+b`, `c`, and `d`, with
multiplicities `2,1,1`. Evaluating their tensor products gives a cheap
column-rank lower bound. The checker implements both this lower bound
and the exact formula (3); it retains a direct column-rank routine for
small general tensors.

## 3. A second exact column formula

Put `W=span{b,c,d}` and `u=C0=a+b`. In the original row coordinates,
set

\[
X=(0,0,1,1)^\mathsf T,\qquad
Y=(0,1,0,1)^\mathsf T.
\]

Then `c=XX^T`, `d=YY^T`, and `b=XY^T+YX^T`. Define

\[
L=\begin{pmatrix}1&1\\1&0\end{pmatrix},\qquad
R=\begin{pmatrix}1&1\\0&1\end{pmatrix}.
\]

The constant two-sided column transformation
`M -> (L tensor L) M (R tensor R)` sends

\[
c\mapsto\operatorname{diag}(1,0,0,0),\quad
b\mapsto\operatorname{diag}(0,1,1,0),\quad
d\mapsto\operatorname{diag}(0,0,0,1),                    \tag{4}
\]

and sends `u` to the rank-one matrix with packed rows `[0,7,7,7]`.
This is a two-sided equivalence, not a simultaneous conjugation; no
algebra assertion is inferred from it.

For arbitrary `T` in `W^tensor(t)`, its column rank is therefore exactly
the number of nonzero entries of its active `2^t` by `2^t` row matrix.
Equivalently, each nonzero coefficient word contributes weight
`2^(number of b letters)`. Different coefficient words occupy disjoint
diagonal positions under (4), so there is no cancellation between them.

## 4. Rank-one tensors over every characteristic-two field

### Row rank one

If `r(T)=1`, formula (3) makes every proper-`S` block zero, because its
multiplicity is at least two. Thus `T` lies in `W^tensor(t)` and its
active matrix has rank one.

Locally, `W` is the space of trace-zero two by two matrices. Consequently,
every one-site partial trace of an element of `W^tensor(t)` is zero.
Write a rank-one active matrix as `xy^T`. At site `i`, split the vectors
as `x=(x0,x1)` and `y=(y0,y1)`. The partial-trace condition is

\[
x_0y_0^\mathsf T=x_1y_1^\mathsf T.                       \tag{5}
\]

If the two sides are nonzero, equality of rank-one outer products forces
`x1` to be proportional to `x0` and `y1` to be proportional to `y0`.
If both sides are zero, the nonzero vectors `x,y` force their supported
local basis coordinates to be opposite, which again gives local
factorizations. Applying this at every site makes `x` and `y` pure
tensor products.

It follows that every row-rank-one tensor is a product of local
rank-one trace-zero matrices. Up to a nonzero scalar, these have the form

\[
\begin{pmatrix}xy&x^2\\y^2&xy\end{pmatrix}
\quad\longleftrightarrow\quad xy\,b+x^2c+y^2d.            \tag{6}
\]

This parameterization works over arbitrary characteristic-two fields:
when the `c` coordinate is nonzero, take a scalar multiple of the
parameter pair `(1,b/c)`; when it is zero, the rank-one condition forces
the remaining point `d`.

The local column rank is one at the two points `c,d`, and four when
both parameters are nonzero. Thus a row-rank-one tensor has column
rank `4^k` for some `0<=k<=t`. Over `F2`, the only local projective
choices are `c,d,v=b+c+d`, and there are exactly `3^t` nonzero
row-rank-one tensors.

### Column rank one

The only local column-rank-one lines are `u,c,d`. Indeed, (2) first
forces `p=0`, because otherwise two diagonal entries are nonzero.
It also prohibits `c` and `d` from both being nonzero. The remaining
two by two minors then force either `a=0` with only `c` or `d`
nonzero, or `c=d=0`, which gives a multiple of `u`.

To tensorize this classification, let `xy^T` be a column-rank-one
matrix in `V^tensor(t)`. At one site, let `X_i` and `Y_i` be the spans
of the local slices of `x` and `y`. Contracting the other row and
column coordinates independently shows

\[
X_i\otimes Y_i\ \subseteq\ V.                            \tag{7}
\]

The three isolated rank-one lines `u,c,d` contain no product rectangle
of dimension greater than one. To see this directly, every nonzero
local rank-one vector in (7) must lie on one of those three lines;
their left directions and their right directions are pairwise distinct.
Two independent vectors on either side would therefore produce an
additional rank-one line, contradicting the local classification.
Hence both slice spaces have dimension one at every site, and the
tensor factors locally.

Every column-rank-one tensor is therefore a product of `u,c,d`, up
to scalar. Its row rank is `4^k`, where `k` counts its `u` factors.
In particular, `u^tensor(t)` is the unique column-rank-one line
with full row rank `4^t`. These conclusions hold over every
characteristic-two field, not only finite fields.

## 5. The full eight-label span misses the strict target

Consider any finite decomposition

\[
E^{\otimes3}=\sum_i T_i,\qquad T_i\in V^{\otimes3},
\]

over any characteristic-two field. Write `r_i,c_i` for the two ranks
and `S=sum_i sqrt(r_i c_i)`. The two full cut dimensions are 64.
Remove zero pieces and combine collinear copies of
`U=u^tensor(3)`; this cannot increase the cost.

Let `mathcal U` be the single line spanned by `U`. Let `mathcal V`
be the entire class of row-rank-one tensors with column rank 64.
There can be many lines in `mathcal V` over an extension field.
Every nonzero piece in either class costs 8.

The rank-one classifications imply

\[
T\notin\mathcal U\Longrightarrow r(T)/c(T)\le32,
\qquad
T\notin\mathcal V\Longrightarrow c(T)/r(T)\le32.           \tag{8}
\]

For example, a column-rank-one piece outside `mathcal U` has row
rank at most 16; any piece of column rank at least two satisfies
`r/c<=64/2`. The other implication is identical with the cuts exchanged.

There are three cases.

1. **Neither class occurs.** For `1/32<=r/c<=32`,
   `sqrt(rc)>=sqrt(32)(r+c)/33`. Rank subadditivity in each cut gives
   `sum r_i>=64` and `sum c_i>=64`, so

   \[
   S\ge128\sqrt{32}/33>19.
   \]

2. **Only one class occurs.** Select one cost-8 piece. Removing it
   leaves opposite-cut rank at least 63. If the selected piece is in
   `mathcal U`, every remaining piece is outside `mathcal V`, so
   `sqrt(rc)>=c/sqrt(32)`. If it is in `mathcal V`, use
   `sqrt(rc)>=r/sqrt(32)`, since `mathcal U` is absent. Thus

   \[
   S\ge8+63/\sqrt{32}>19.
   \]

3. **Both classes occur.** Three endpoint-class pieces already cost
   24. Otherwise there is exactly one of each. The row matrix of `U`
   is the identity. The flip `P8` satisfies
   `rank(P8+alpha I)>=8*7/2=28` for every scalar `alpha`:
   at `alpha=1`, each unordered pair contributes one independent
   rank-one block, and for `alpha!=1` the matrix is invertible.
   Removing the selected row-rank-one piece leaves row rank at
   least 27. Since every nonzero residual piece has column rank at
   least one, its total cost is at least `sqrt(27)`. Consequently,

   \[
   S\ge16+\sqrt{27}>19.
   \]

All comparisons have integer witnesses:

\[
128^2\cdot32>19^2\cdot33^2,\quad
63^2>11^2\cdot32,\quad27>3^2,
\]

and

\[
19^{70}>2^{297},\qquad99^2>2\cdot70^2.
\]

Therefore

\[
\boxed{\quad S>19>8^{99/70}>8^{\sqrt2}.\quad}              \tag{9}
\]

This settles the full third coefficient tensor power, including arbitrary
nonsymmetric recombinations and every characteristic-two extension field.
It does not tensorize automatically to higher powers.

## 6. Explicit structured pieces from the conic

The rank-one local polynomial

\[
T(z)=c+zb+z^2d=(X+zY)(X+zY)^\mathsf T                     \tag{10}
\]

has useful exact coefficient ranks. For
`T(z)^tensor(t)=sum_j T_j z^j`,

\[
\boxed{\quad
r(T_j)=\min(j+1,2t-j+1,t+1),\qquad
c(T_j)=\binom{2t}{j}.
\quad}                                                   \tag{11}
\]

For the row formula, expand `(X+zY)^tensor(t)`. Its degree-`k`
vectors are independent, because they have disjoint supports after
changing `X,Y` to coordinate vectors in their two-dimensional span.
At degree `j`, the distinct pairs `k+l=j` therefore give exactly the
stated number of independent rank-one matrix summands.

For the column formula, (4) changes `T(z)` into
`diag(1,z,z,z^2)`. Tensoring gives exactly `binom(2t,j)` nonzero
diagonal positions in degree `j`. These are integer counts of positions,
so the conclusion survives characteristic-two binomial cancellations
in scalar polynomial arithmetic.

Formula (11) supplies pieces whose row rank grows only linearly in `t`.
It does not by itself improve the decomposition cost of `v^tensor(t)`,
which already has row rank one and column rank `4^t`.

## 7. Complete column-rank-two classification

There is an exact classification at every tensor power and over every
characteristic-two field. In the coefficient basis `u,b,c,d`, apply to
the column forms in (4) the further two-sided matrices with packed rows
`[2,1,10,6]` on the left and `[2,8,11,4]` on the right. They send
`u,c,d` to `E00,E11,E22`, respectively, and send `b` to

\[
B=\begin{pmatrix}
0&0&0&1\\
0&0&0&0\\
0&0&0&1\\
1&1&0&0
\end{pmatrix}.                                           \tag{12}
\]

Thus the three letters `u,c,d` occupy three independent diagonal
positions. The letter `b` connects the fourth position to the first
and third on one side, and to the first and second on the other.

**Protected coefficient minor.** If any nonzero coefficient word has
`s` occurrences of `b`, the column rank is at least `2^s`. Fix the
other sites to their diagonal positions, and restrict the `b` sites
to positions `{0,3}`. There `b` is the two by two flip, while only `u`
survives among the three diagonal letters, as `E00`. Multiplying on
one side by the tensor flip changes the designated `b` word to a
nonzero scalar identity; every other surviving word is strictly
triangular. The resulting minor of size `2^s` is invertible.

In particular, column rank at most two forces every coefficient word
to have at most one `b`. Order matrix indices by the number of sites
equal to the fourth position. The nonzero part then has form

\[
\begin{pmatrix}D&P\\Q&0\end{pmatrix},                    \tag{13}
\]

where `D` is diagonal on the `3^t` ternary words. A nonzero one-`b`
coefficient indexed by a site and a ternary word at the other sites
contributes a column of `P` supported at exactly the two corresponding
positions `0,2`. Distinct one-`b` words have distinct two-point
supports. The corresponding row of `Q` has positions `0,1`.

For every matrix of form (13), rank is at least
`rank(P)+rank(Q)`: choose nonsingular minors of `P` and `Q`; placing
their row and column sets in opposite orders gives a block triangular
minor with their determinants on its diagonal.

If a `b` term occurs and total rank is two, both `P,Q` must have rank
one. Distinct two-point supports cannot be proportional over any field,
so precisely one one-`b` coefficient remains. The diagonal `D` must
then be zero. Indeed, outside the three involved ternary positions,
any nonzero diagonal entry adds an independent rank. Inside those
positions, pivot on the invertible two by two flip using position `0`
and the fourth position. The residual diagonal entries at `1,2` force
their coefficients to vanish, and the remaining off-diagonal Schur
entry is the coefficient at `0`, forcing it to vanish as well.

This proves that every nonzero column-rank-two tensor has exactly one
of these forms:

1. A sum of nonzero multiples of two distinct product words from
   `{u,c,d}^tensor(t)`.
2. One `b` factor and, at every other site, one of `u,c,d`, up to scalar.

Conversely, every tensor of either form has column rank two. This is
an all-field theorem; the complete binary square check is an independent
finite control, not its proof.

### Exact row ranks of all the classified pieces

For type 2 the row rank is `2*4^k`, where `k` counts the `u` factors.
For type 1, write the distinct product words as `P,Q`. All nonzero
scalar coefficients give the same ranks below.

If any site has `c` in one word and `d` in the other, both their image
spaces and their row spaces intersect trivially. Thus
`r(P+Q)=r(P)+r(Q)`.

Otherwise let `h` count common `u` sites, let `A` count sites where
only `P` has `u`, and let `B` count sites where only `Q` has `u`.
If either `A` or `B` is zero, the rank is `4^(h+max(A,B))`: after
factoring the common sites, the remaining matrix is an invertible
identity plus a square-zero rank-one map. If both are positive, the
rank is

\[
r(P+Q)=4^h(4^A+4^B-2).                                   \tag{14}
\]

To verify (14), factor the common sites and write the remainder as
`I tensor N_B + N_A tensor I`. Each `N` is a square-zero rank-one
map, hence has one Jordan block of size two and otherwise zero
blocks. The two nonzero Jordan blocks contribute rank two; their
products with the zero blocks contribute `4^A-2` and `4^B-2`.

Over `F2`, there are exactly
`binom(3^t,2)+t*3^(t-1)` column-rank-two tensors. This gives 42 at
power two, 378 at power three, and 3348 at power four. The checker
provides a generator with exact row ranks for all of them.

A useful consequence is that the only column-rank-two pieces whose
row rank exceeds half the full dimension `n=4^t` are
`alpha*u^tensor(t)+beta*P`, where `P` is a different rank-one product
word and both scalars are nonzero. Their row rank is exactly `n`.
Every other column-rank-two piece has row rank at most `n/2`.

## Reproduction

Run the standard-library checker:

```sh
python3 tools/cnot_full_span.py --output /tmp/cnot-full-span-audit.json
```

The checker verifies both local normal forms, the independent two-sided
column diagonalization, direct binary comparisons through tensor power
three, controls over `F4`, the local rank-one varieties over `F4`, and
all rank-one product controls through power four. It also checks the
exact integer inequalities in (9) and direct coefficient ranks in (11)
through power four. It exhausts all 65,535 nonzero binary tensors in
the square to check the rank-two classification, and verifies the
explicit row ranks of all 3348 classified binary fourth-power pieces.
The all-field and all-power statements follow from
the algebraic arguments above, not from finite enumeration alone.
