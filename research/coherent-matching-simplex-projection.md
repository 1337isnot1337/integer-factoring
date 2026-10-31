# A coherent matching projection with an exact quadratic cancellation

## Result and scope

The fifteen-letter [matching crossing](perfect-matching-sector-crossings.md)
has an exact coherent projection to four or five letters. The construction
uses the four-dimensional binary simplex code inside its physical matching
space and an isotropic dual. It preserves the literal crossing endpoint.

The four-letter tensor is a pencil. The five-letter tensor is cubic,
its quadratic coefficient vanishes, and its cubic coefficient is a single
pure tensor with rank one in both cuts. The complete coefficient rank
pairs are nevertheless

\[
 q=4:\ (8,8),(12,12),\qquad
 q=5:\ (8,8),(20,20),(0,0),(1,1).                         \tag{1}
\]

The canonical examples are uncompetitive in the recorded pressure search.
These exact cancellations supply a concrete coherent construction; they
do not improve the factoring exponent or exclude other physical maps.

## 1. The physical simplex subspace

Let the fifteen coordinates be the perfect matchings of six vertices.
Let `A` be the binary adjacency matrix joining matchings with exactly one
common edge, and let `J` be the all-ones matrix. Define

\[
 N=A+I+J.                                                \tag{2}
\]

Direct finite arithmetic gives `N^2=0`, `NJ=0`, and `rank N=4`.
Choose four independent columns `p_1,...,p_4`. Reading their coordinates
as four-bit column vectors lists every nonzero member of `F_2^4`
exactly once. Thus their span is the binary simplex code on fifteen
nonzero points. Its vectors are isotropic for the ordinary dot product,
and all are orthogonal to the vector `one`.

Set `P_0=(p_1,...,p_4)` and `P=(P_0,one)`. These are injections of
four- and five-dimensional physical spaces into the matching space.
In particular

\[
 P_0^T P_0=0,\qquad P^TP=\operatorname{diag}(0,0,0,0,1).
 \tag{3}
\]

All identities here concern explicit matrices of order fifteen. The
constructor and checker retain their complete coordinates; no assertion
about an arbitrary-size modular representation is needed.

## 2. An isotropic left inverse

First solve for vectors `r_i` satisfying

\[
 r_i\mathbin\cdot p_j=\delta_{ij},\qquad
 r_i\mathbin\cdot one=0.                                 \tag{4}
\]

The Gram matrix `S_ij=r_i dot r_j` is alternating, because in
characteristic two `r_i dot r_i=r_i dot one=0`. Choose the strictly
upper-triangular matrix `C` with `C+C^T=S`, and replace

\[
 r_i\longmapsto r_i+\sum_j C_{ij}p_j.                    \tag{5}
\]

The pairings (4) are unchanged. By (3), the new Gram matrix is
`S+C+C^T=0`. Append the row `one` for the five-letter example. With
these rows forming `R`,

\[
 RP=I_5,\qquad RR^T=P^TP=\operatorname{diag}(0,0,0,0,1).
 \tag{6}
\]

Omitting the last row and column gives the analogous four-letter split
with both Gram matrices zero.

Apply the paired physical maps `(R,R,P^T,P^T)` to the complete matching
tensor. The [split-projection identity](physical-split-crossings.md)
and `RP=I` prove the strict endpoint. The implementation constructs the
whole polynomial, rather than transferring profile orders through a
noninvertible map.

## 3. Why the quadratic coefficient vanishes

Write the original cubic matching tensor as `sum_(j=0)^3 T_j w^j`.
Its quadratic coefficient is exactly

\[
 (T_2)_{xyuv}=[x=y][u=v]N_{xu}.                          \tag{7}
\]

Every column of `N` belongs to the span of the four nonconstant columns
of `P`. On the fifteen nonzero points of `F_2^4`, every product of at
most three affine-linear coordinate functions with at least one
nonconstant linear factor has sum zero. To see this, include the zero
point, where the product is zero. Every resulting squarefree monomial
has degree at most three, so its sum over the full four-dimensional
cube is an even integer and hence zero in characteristic two.

Consequently

\[
 \sum_u P_c(u)P_d(u)N_{xu}=0                             \tag{8}
\]

for every `x,c,d`. This kills (7) after the two output projections,
before any choice of the left inverse. Thus the quadratic cancellation
holds for **every** left inverse of this fixed simplex embedding over
every characteristic-two extension field.

The cubic coefficient of the original tensor is `[x=y][u=v]`. After
the split it becomes

\[
 (RR^T)_{xy}(P^TP)_{uv}.                                 \tag{9}
\]

For the isotropic split (6), this is the single coordinate with all
four logical labels equal to the trivial-line label. It has rank one
in both cuts. In the four-letter projection both Gram matrices vanish,
so both the quadratic and cubic coefficients disappear.

The remaining coefficient ranks in (1) are exact finite rank
computations; no stronger cancellation is asserted from (8)--(9).

## 4. Complete profiles and verification

The two cuts of each canonical example have the same unweighted
profiles:

| Alphabet | Zero orders | Infinity orders |
| ---: | --- | --- |
| 4 | `0:8, 1:8` | `-1:12, 1:4` |
| 5 | `0:8, 1:16, 3:1` | `-3:1, -1:20, 1:4` |

The [exact checker](../tools/check_matching_coherent_projections.py)
contracts each complete coefficient independently using Kronecker maps,
checks every simplex triple product, verifies both physical Gram
matrices, and obtains the complete profiles from binary block-Toeplitz
ranks. The [verification record](../evidence/matching-coherent-projections.json)
also stores both split maps.

```sh
python3 tools/check_matching_coherent_projections.py
```

A finite full-height search of the canonical five-letter seed used
554 profile attempts, eight of which were unresolved. Its best
recorded `rho` was about `1.00000024`. This search does not certify a
global minimum; the exact construction and coefficient cancellations
are the result of this checkpoint.

The proof, construction, and exact controls used substantial OpenAI
assistant support. No Lean source was changed.
