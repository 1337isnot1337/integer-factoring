# Symplectic relations among the Gaussian complements

## Exact result

The 65 four-particle complements of the retained five-state reflection
code split into two parts:

- 11 other reflection-admissible states, with independent slab rank 11;
- 54 nonadmissible states, with independent slab rank 26.

There are 29 explicit symplectic contraction equations among the
nonadmissible slabs. They have rank 28 and explain the total rank
`65-28=37` in the
[complete aligned module](aligned-gaussian-complement-module.md).
The equations follow uniformly from symplectic invariance. Their
independence and equality cases are also checked on the complete
retained polynomial source.

The two sectors interact through polynomial cancellation. At degree
18, capping each endpoint-zero sector separately gives 1,250 tensor
directions. Their full polynomial module sum gives 1,655. Allowing
nonzero sector endpoints to cancel raises the dimension to 1,700.
Saturating the combined rational span and imposing the endpoint gives
2,120. These are four distinct exact spaces; none is a pressure bound.

## 1. The contraction operators

Let `E` be an eight-dimensional symplectic space with coordinate pairs

\[
 \{e_i,e_{7-i}\},\qquad 0\le i<4.
\]

Write `p_i` for the two-element support of the `i`th pair. A
reflection-admissible four-subset contains one coordinate from each
pair. There are sixteen such subsets.

For any symplectic space with form `omega`, define a natural divided
contraction on an exterior product by

\[
 C_r(a_1\wedge\cdots\wedge a_k)
 =\sum_{|I|=2r}
    \operatorname{Pf}(\omega(a_i,a_j))_{i,j\in I}
    \bigwedge_{j\notin I}a_j.                          \tag{1}
\]

All signs coincide in characteristic two. The Pfaffian formula uses
no division by a factorial. Because it depends only on `omega`, it
commutes with every symplectic change of coordinates. On the indicated
coordinate basis, `C_r` removes each union of `r` complete reflection
pairs, with coefficient one.

Every reflection-admissible basis vector is annihilated by all positive
`C_r`. It also cannot be enlarged by adjoining a disjoint complete
reflection pair: each pair already contains an occupied coordinate.

## 2. Apply the contractions to a selected minor

The one-particle strict Gaussian is symplectic after the diagonal
conjugation in [reflection-pair compression](reflection-dual-rail-compression.md).
Its two physical blocks are copies of `E`. For four-particle labels on
all boundaries, the conjugating monomial cancels between the selected
row and column subsets, so their minors are unchanged.

Fix admissible labels on the three boundaries `y,u,v`, and let
`F_alpha` be the selected minor with four-subset `alpha` on the
remaining boundary `x`. The transpose exterior action also commutes
with (1). Its input, consisting of the two fixed admissible row labels,
is killed by every positive contraction.

Take the coefficient of a fixed admissible `u` after contracting the
other block. A mixed contraction involving the `u` block would have
to arise from a larger support obtained by adjoining a complete pair
disjoint from `u`. There is no such pair. Thus only contractions on
the `x` block remain.

For every two-subset `beta`, the first contraction gives

\[
 \boxed{\sum_{i:\,p_i\cap\beta=\varnothing}
               F_{\beta\cup p_i}=0.}                  \tag{2}
\]

The second divided contraction gives

\[
 \boxed{\sum_{0\le i<j<4}F_{p_i\cup p_j}=0.}          \tag{3}
\]

Every support in (2)--(3) contains a complete pair, so none is an
admissible code state. The equations therefore hold entirely within
the complementary source slabs. This proof applies to every
symplectic one-particle matrix with the three fixed admissible
boundaries, including rational matrices whose selected minors have
been made polynomial by a common factor.

## 3. Exactly 28 independent relations

The 54 nonadmissible four-subsets have two types. There are 48 with
one complete pair, one empty pair, and single occupations on the other
two pairs. There are six consisting of two complete pairs.

When `beta` chooses one coordinate from each of two different pairs,
(2) has two terms. There are 24 such equations, and their supports
partition the 48 subsets into disjoint pairs. Their rank is 24.

When `beta` is itself a complete pair, the four equations (2) are the
vertex-star vectors of the complete graph on the four pair indices.
Their six coordinates are the unions of two complete pairs. Over a
characteristic-two field, these four star vectors have rank three.
Their sum is zero, and any three are independent.

Equation (3) is the all-ones vector on those six graph edges. It is
outside the star span: every star combination has even sum around a
triangle, whereas the all-ones vector has odd sum. Thus (3) adds one
more independent equation. The total relation rank is

\[
 24+3+1=28.
\]

The resulting nonadmissible coefficient space has dimension at most
`54-28=26`. The retained source achieves rank 26 both over the binary
field, when every polynomial coefficient is packed, and over the
rational function field. Its eleven admissible complementary slabs
are independent and independent of the nonadmissible part. This
gives the exact total rank 37.

The corresponding complete four-particle primitive space has
dimension `16+26=42`. The five selected admissible states and their
37 independent complementary directions account for these 42
dimensions in the retained example.

## 4. Separate sectors, polynomial cancellation, and saturation

Let `V_A,V_N` be the two rational slab spaces and let

\[
 M_A=V_A\cap\mathbf F_2[z]^{125},\qquad
 M_N=V_N\cap\mathbf F_2[z]^{125}.
\]

Their ranks are eleven and twenty-six. The exact degree-reduced
bases have these degrees.

| Sector | Minimal generator degrees |
| --- | --- |
| Admissible | `9,12,12,13,16,17,17,17,17,18,20` |
| Nonadmissible | `4,7,8,8,8,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,12,12,13,13,13` |

Each basis has its own complete saturation certificate, using the
bounded saturation lemma: degree 19 and dimension 52 for `M_A`, and
degree 12 and dimension 88 for `M_N`.

If each endpoint-zero component is separately capped at degree 18,
formula `(1+z)M_A` gives 160 tensor directions and `(1+z)M_N` gives
1,090. Their direct sum under those separate caps has dimension 1,250.

Put `H=M_A+M_N`. The sum of the two complete endpoint-zero modules is
`(1+z)H`; it permits terms of larger individual degree to cancel.
The checker performs only
unimodular polynomial row operations on the combined 37 generators,
until their leading positions are distinct. The resulting weak Popov
basis has degree sum 335 and maximum degree 13. Distinct leading
positions imply independent leading coefficient rows, so its exact
Hilbert function gives

\[
 5\dim H_{\le17}=5(37\cdot18-335)=1655.                \tag{4}
\]

Finally, the full saturated module
`M=(V_A+V_N) intersect F_2[z]^125` has degree sum 242. Therefore

\[
 5\dim M_{\le17}=5(37\cdot18-242)=2120.                \tag{5}
\]

There is one additional distinction inside `H`. It is not saturated,
so some elements of `H` vanish at one without being divisible by
`1+z` inside `H`. The endpoint image of its 37 generators has exact
rank 28. Since all weak Popov basis degrees are at most 13, evaluation
on the degree-at-most-18 part has that same rank. Thus its complete
endpoint-zero subspace has tensor dimension

\[
 5\bigl(37\cdot19-335-28\bigr)=1700.                   \tag{6}
\]

The polynomial quotient `M/H` has finite dimension `335-242=93` over
the binary field. The 870 directions beyond the separately capped
spaces therefore split into 405 supplied by degree cancellation,
45 supplied by cancellation of the two sector endpoints, and 420
additional directions after full saturation. This distinguishes
three different ways that combining the sectors can enlarge a
bounded endpoint-zero space.

## 5. Reproduction

```sh
python3 tools/check_gaussian_complement_symplectic_relations.py
```

The [checker](../tools/check_gaussian_complement_symplectic_relations.py)
and [retained exact record](../evidence/gaussian-complement-symplectic-relations.json)
verify all 3,625 polynomial entries of (2)--(3), their exact relation
rank, both constant and rational slab ranks, the separate module
saturation certificates, and the weak Popov reduction of their sum.
They also test the contraction identities on two independently
generated binary symplectic matrices, including matrices that do not
preserve the corresponding quadratic form. These latter controls
check the implementation of the natural contraction formula; the
uniform statement follows from (1).

This derivation and implementation used substantial OpenAI assistant support.
