# Cubic crossings with arbitrary dual output maps

The cubic family has a fixed noncrossing tensor even when every output
row is changed. This reduces the remaining search to two variable tensor
directions. It also gives a uniform decomposition cost upper bound.
It is **not a lower bound** excluding arbitrary output maps.

An explicit new output map gives a tensor of cut ranks `(16,15)` in the
seven-label coefficient space. Its individual cost is below the finite
seed threshold, but the complete binary space still fails the criterion.

## Setup and universal identity

Let `d>=3` be odd. Work over any field of characteristic two, with `d`
physical pairs of modes `e_i,f_i`. Put

\[
\Omega=\sum_{i=0}^{d-1}e_if_i,\qquad
I_i=\Omega e_i,\qquad J_i=\Omega f_i.
\]

Let `Pi` be an arbitrary linear map from exterior grade three to the
logical `d`-dimensional space, subject to

\[
\Pi I=I_d,\qquad \Pi J=0.
\tag{1}
\]

The map may mix physical syndromes and use nonbinary field coefficients.
Use `I` on the two input legs and `Pi` on the two output legs of the
universal Pfaffian tensor. Denote the result by `G(p,q,r)`, with logical
indices `(x,y,v,u)` as in the other crossing notes. Define

\[
E_{xy,vu}=\delta_{xv}\delta_{yu},\qquad
N_{xy,vu}=\delta_{xu}\delta_{yv}.
\]

**Theorem.** There are tensors `A,B` such that, with `eta=q+r`,

\[
\boxed{
G(p,q,r)=\eta\bigl[r\eta E+(p\eta+qr)A+p^2B+qrN\bigr].
}
\tag{2}
\]

In particular, the complete coefficient space is contained in
`span{E,N,A,B}`. The two tensors `E,N` are independent of the output
map. All output freedom is confined to `A,B`. These may depend
quadratically on `Pi`; equation (2) does not make their cut ranks small.

### Physical identity before applying the output map

Let `X,U` range over physical degree-three occupations. Let
`K(X,U;y,v)` be the Pfaffian entry after inserting `I_y,I_v` in the
inputs, before applying either output map. Its universal scalar
polynomials are

\[
r(q+r)^2,\quad p(q+r)^2,\quad p^2(q+r),\quad
qr(q+r),\quad q^2(q+r).
\tag{3}
\]

They are dual to the monomial coefficient positions
`r^3, p r^2, p^2 r, q r^2, q^3`. Write the corresponding physical
tensors as `E_0,A_0,B_0,C_0,H_0`. The exact identities are

\[
\begin{aligned}
E_0(X,U;y,v)&=I_v[X]I_y[U],\\
H_0(X,U;y,v)&=J_y[X]J_v[U],\\
(A_0+C_0)(X,U;y,v)
&=\delta_{yv}\sum_{w\in\{e_i,f_i\}}(\Omega w)[X]\,(\Omega w)[U]
  +J_y[X]J_v[U].
\end{aligned}
\tag{4}
\]

Apply `Pi` to both physical output indices. Conditions (1) send
`E_0` to `E`, send `H_0` to zero, and send `A_0+C_0` to `N`.
Substitution into (3) gives (2). Because (3)--(4) are binary polynomial
identities, this argument applies over every field of characteristic
two, including nonbinary choices of `Pi`.

### Finite verification with an all-size reduction

The proof of (3)--(4) uses an exhaustive local-pattern identity. Its
finite scope covers all odd sizes for the following reason.

Fix `X,U,y,v`. At most eight physical sites occur among the supports of
`X,U` and the two distinguished input sites `y,v`. In the expansion of
`I_y=Omega e_y` and `I_v=Omega e_v`, an additional empty site can
contribute only when both copies of `Omega` use that site. Its local
Pfaffian factor is one. If exactly one copy uses an empty site, the
local factor is zero. Thus all unused sites contribute only through
their number modulo two. For odd `d`, replace the number of sites by
nine without changing the entry.

Permute the sites other than `y,v` into a fixed order. There are exactly
**4,684** canonical pairs of degree-three output occupations, allowing
both `y=v` and `y!=v`. The checker enumerates all of them, using the
literal local Pfaffian polynomials, and verifies (3)--(4). It finds 46
nonzero canonical patterns. This is a complete finite identity check,
followed by the empty-site parity reduction, rather than a sample of
large sizes.

## A uniform two-piece decomposition

Both `N` and `E+N` belong to the coefficient space, and

\[
E=N+(E+N).
\]

The ranks are independent of `Pi`:

| Piece | Row-cut rank | Column-cut rank |
| --- | ---: | ---: |
| `N` | `d^2` | `1` |
| `E+N` | `d(d-1)/2` | `d^2-1` |

For the row cut, `N` swaps the two logical indices. The map `I+swap`
has one independent image for each unordered pair of distinct indices.
For the column cut, `N` is the outer product of the vectorized identity
matrix with itself. Its self-pairing is `d=1` in the field, so `E+N`
has rank `d^2-1`.

Consequently every map satisfying (1) admits the cost bound

\[
\boxed{S\le d+(d-1)\sqrt{\frac{d(d+1)}2}.}
\tag{5}
\]

This does not meet the strict finite-seed target at the recorded sizes.
It explains the repeated seven-label decomposition cost `7+12 sqrt(7)`
in searches that change many output rows. A better decomposition would
have to use the variable directions `A,B` more effectively.

## A new individually affordable tensor

At `d=7`, use cyclic coordinate outputs together with a common low-mode
quadratic correction. In coordinate-functional notation, set

\[
\Pi_i=[e_iD_{i+1}]+[e_i(e_0e_1+e_2e_3)],
\qquad i\pmod7,
\tag{6}
\]

where a repeated low mode gives zero. The triple-single corrections
annihilate every `I_j,J_j`, so (1) holds. Let `B` be the `p^2 r`
coefficient in (2). The exact ranks are

\[
(r(E+B),c(E+B))=(16,15),\qquad (r(B),c(B))=(42,40).
\]

Thus `E+B` has cost `sqrt(240)`, strictly below `7^sqrt(2)`. The checker
certifies this without a floating-point comparison:

\[
1393^2<2\cdot985^2,\qquad 240^{985}<7^{2786}.
\]

An individually affordable tensor is insufficient. In this case its
complement `B` is expensive, and the complete binary coefficient space
has integer shortest-path lower bound 38. It therefore does not give a
strict seed. The known two-piece upper bound remains `7+12 sqrt(7)`.

## Multirow experiment and replay

The archived experiment contains 680 exact output maps:

- 168 coherent pair corrections: choose two physical sites and eight
  low/high mode patterns, changing the other five output rows.
- 512 seeded cases alternating between two changed rows, three changed
  rows, seven changed rows, and a shared low-mode quadratic form.

Every output correction is explicitly recorded. All 680 spaces have
coefficient dimension four. They give 495 distinct binary rank tables,
but none improves the recorded two-piece cost. Each complete binary
space has integer shortest-path lower bound 38.

The independent replay uses the opposite pivot order for matrix ranks
and Bellman relaxation for integer XOR costs. It checks all 10,200
nonzero coefficient vectors, every endpoint entry, and nine direct
bilinear expansions of the output maps. The fixed tensor `N` is checked
entry by entry in every case. The search is not exhaustive over all
output maps; equations (2)--(5) have the wider scope stated above.

```sh
python tools/explore_cubic_multirow_outputs.py \
  --output /tmp/cubic-multirow-outputs.json
python tools/check_cubic_output_identity.py \
  --multirow evidence/cubic-multirow-outputs.json \
  --output /tmp/cubic-output-identity.json
```

- [Multirow construction and search](../tools/explore_cubic_multirow_outputs.py).
- [All-size identity and independent replay](../tools/check_cubic_output_identity.py).
- [All 680 output maps and rank tables](../evidence/cubic-multirow-outputs.json).
- [Identity, affordable tensor, and replay certificates](../evidence/cubic-output-identity.json).
- [Earlier untouched-output gate](cubic-output-perturbations.md).
