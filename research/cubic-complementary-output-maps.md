# Allowing a nonzero complementary output image

The cubic output map need not annihilate the complementary-mode input
family. Allowing that image gives a five-direction coefficient formula
and improves the recorded seven-label cost from `7+12 sqrt(7)` to
**36** over `F4`. This is an improvement within the construction, not a
strict finite seed: the target remains `7^sqrt(2)`, approximately 15.673.

## A five-direction formula for arbitrary dual outputs

Use the odd-dimensional cubic setup from
[the arbitrary-output identity](cubic-arbitrary-output-identity.md):
`I_i=Omega e_i`, `J_i=Omega f_i`, with `d` physical pairs and `d` logical
labels. Require only `Pi I=I_d`, and define

\[
C=\Pi J,\qquad M=I_d+CC^{\mathsf T}.
\]

There is no transpose conjugation here: the field has characteristic two,
and the displayed product is the ordinary bilinear matrix product.
Define two logical tensors

\[
(P_C)_{xy,vu}=\delta_{yv}M_{xu},\qquad
(H_C)_{xy,vu}=C_{xy}C_{uv}.
\tag{1}
\]

The physical identity in the linked note gives

\[
[pr^2]G+[qr^2]G=P_C+H_C,\qquad [q^3]G=H_C.
\tag{2}
\]

With `A=[pr^2]G`, `B=[p^2r]G`, and `eta=q+r`, the whole polynomial is

\[
\boxed{
G=\eta\bigl[r\eta E+(p\eta+qr)A+p^2B+qrP_C+q\eta H_C\bigr].
}
\tag{3}
\]

Thus the complete coefficient space lies in
`span{E,A,B,P_C,H_C}`. Equation (3) holds over every field of
characteristic two and for arbitrary degree-three output maps satisfying
the duality condition. It follows by applying `Pi` to the physical
five-polynomial identity; no new local-pattern assumption is needed.

The two distinguished pieces have exact ranks

\[
\begin{aligned}
r(P_C)&=d\,\operatorname{rank}M,
&c(P_C)&=\begin{cases}1&M\ne0,\\0&M=0,\end{cases}\\
r(H_C)&=\begin{cases}1&C\ne0,\\0&C=0,\end{cases}
&c(H_C)&=(\operatorname{rank}C)^2.
\end{aligned}
\tag{4}
\]

The crossing endpoint is unchanged. The second noncrossing endpoint is
now `H_C`; it is allowed to be nonzero. All recorded searches retain that
coefficient instead of imposing the earlier zero-endpoint condition.

## A nonbinary scalar image gives three pairing tensors

Take `C=lambda I_d`, where `lambda` is neither zero nor one. Then
`P_C=(1+lambda^2)P` and `H_C=lambda^2 H`, where

\[
P_{xy,vu}=\delta_{xu}\delta_{yv},\qquad
H_{xy,vu}=\delta_{xy}\delta_{uv}.
\]

Both scalar factors are nonzero. Consequently `P` and `H` belong to the
coefficient space, and

\[
E=P+H+(E+P+H).
\tag{5}
\]

The two noncrossing pairings have respective cut ranks `(d^2,1)` and
`(1,d^2)`. The third piece has both cut ranks

\[
\binom d2+1.
\]

For the row cut, its matrix is identity plus the swap operator plus the
outer product of the vectorized identity with itself. The identity-plus-
swap image has one vector for each unordered pair of distinct indices,
and all of those vectors have zero diagonal. The extra outer product
adds one independent diagonal direction. The column cut interchanges
the two noncrossing pairings and gives the same rank.

For every odd `d>=3` in this cubic setup, this proves the upper bound

\[
\boxed{S\le2d+\binom d2+1=\frac{(d+1)(d+2)}2.}
\tag{6}
\]

The earlier two-piece bound is also available because `P` belongs to
the space. One may take the smaller of the two bounds. Neither bound
establishes the strict seed target in the recorded cases.

### Exact seven-label optimum over `F4`

Let `F4=F2[w]/(w^2+w+1)` and take `C=w I_7`. Use cyclic coordinate
outputs

\[
\Pi_i=[e_iD_{i+1}]+w[f_iD_{i+1}],\qquad i\pmod7.
\]

The complete coefficient space has dimension five. In the raw monomial
basis `(r^3, p r^2, p^2 r, q r^2, q^3)`, an exact decomposition is:

| Tensor | Coefficient vector | Row rank | Column rank | Cost |
| --- | --- | ---: | ---: | ---: |
| `E+P+H` | `(1,w^2,0,w^2,1)` | 22 | 22 | 22 |
| `P` | `(0,w^2,0,w^2,w^2)` | 49 | 1 | 7 |
| `H` | `(0,0,0,0,w)` | 1 | 49 | 7 |

Their coefficient vectors sum to `(1,0,0,0,0)`. Their total cost is
exactly 36. The complete `F4` rank table gives an integer shortest-path
lower bound of 36, independently verified by a shortest-path potential.
Since floor square-root costs are lower bounds for the actual costs,
this proves that **36 is the exact optimum in this particular `F4`
coefficient space**. The same result holds for the conjugate choice
`C=w^2 I_7`. Optimality over larger fields is not asserted.

## Recorded searches

The binary search contains 1,195 distinct matrices `C`. It includes all
binary diagonal matrices, all binary circulant matrices, all symmetric
rank-one matrices `vv^T`, all strictly upper triangular Toeplitz
matrices, and the recorded seeded dense, rank-one, and strictly upper
triangular cases. Deduplication is global. Their coefficient dimensions
are:

| Dimension | Cases |
| ---: | ---: |
| 3 | 1 |
| 4 | 7 |
| 5 | 1,187 |

All **36,909** nonzero binary coefficient vectors are checked. Every
recorded space fails the strict seed criterion; the smallest integer
cost lower bound is 38.

The extension-field search contains 56 explicit `F4` matrices: scalar
identities, scalar projections, scalar nilpotent shifts, scalar rank-one
maps, and seeded diagonal, rank-one, and dense maps. It exhausts
**18,520 projective coefficient points** in total. Both pivot orders
are checked for every rank calculation. In every case the points with
cost below 16 span a space missing the crossing, which excludes a strict
seven-label seed. The best recorded cost is the exact value 36 above.

The independent audit replays every binary rank pair, checks 648
extension-field vectors by a separate scalar elimination implementation,
and directly expands 120,050 physical coefficient entries across ten
control maps. The 56-map list is not exhaustive over complementary
images `C`; nor does it exhaust the remaining output-map freedom for a
fixed `C`.

## Reproduction

```sh
python tools/explore_cubic_complement_outputs.py \
  --output /tmp/cubic-complement-outputs.json
python tools/explore_cubic_complement_fields.py \
  --output /tmp/cubic-complement-fields.json
python tools/check_cubic_complement_outputs.py \
  --binary evidence/cubic-complement-outputs.json \
  --field4 evidence/cubic-complement-fields.json \
  --output /tmp/cubic-complement-audit.json
```

- [Binary complement explorer](../tools/explore_cubic_complement_outputs.py).
- [Complete `F4` coefficient-space explorer](../tools/explore_cubic_complement_fields.py).
- [Independent audit](../tools/check_cubic_complement_outputs.py).
- [Binary matrices and rank tables](../evidence/cubic-complement-outputs.json).
- [`F4` matrices and rank tables](../evidence/cubic-complement-fields.json).
- [Replay certificate](../evidence/cubic-complement-audit.json).
