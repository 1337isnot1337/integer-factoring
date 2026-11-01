# Coupled rational physical frames without an added degree cost

## Constructive result and scope

The four-letter and sixteen-letter reflection crossings admit polynomial
deformations obtained from physical frames that have poles at zero.
Their poles cancel only after the actions on both paired alphabets are
combined. The strict crossing endpoint and the full degree are preserved.
The resulting tensors can leave the even Gaussian coefficient locus.

This supplies a new family for complete-profile searches. The exact
construction does not itself establish a smaller pressure bound or an
improved factoring exponent. Every numerical optimization must still use
the complete weighted profiles and one common pair of physical types.

## 1. A single paired frame cannot have a genuine pole here

Write the four physical legs in the order `x,y,u,v`. A paired change on
`x,v` acts by `A(w)` and `A(w)^(-T)`. It preserves the crossing endpoint
`[x=v][y=u]` whenever it is defined there.

The following elementary local statement explains a restriction at zero.
Let `O=F[[w]]`, and suppose the `q^2` matrices obtained by fixing `y,u`
in a tensor form an `O`-basis of `End(O^q)` in the remaining `x,v`
indices. If a matrix `A` over `F((w))` makes every conjugated slice
integral, then

\[
 A\operatorname{End}(O^q)A^{-1}\subseteq\operatorname{End}(O^q).
 \tag{1}
\]

The lattice `A^(-1) O^q` is therefore stable under every matrix unit.
The diagonal matrix units split it into fractional ideals on the
coordinate axes. The off-diagonal units force all those ideals to be
equal. Since `O` is a discrete valuation ring, the lattice is homothetic
to `O^q`. Hence `A=w^h B` for some `B in GL_q(O)`. The scalar power
cancels from its paired action.

For the reflection seeds considered here, the coefficient at zero has
full rank `q^2` in the third flattening `(x,v)|(y,u)`. The exact ranks
are 16 and 256. Its slices therefore satisfy the basis hypothesis.
This argument also survives regular invertible changes on the other
paired alphabet. It does not apply when both pairs have singular
frames simultaneously.

## 2. Normalize the leading tensor to two Bell pairs

The same constant coefficient factors exactly as

\[
 T(0;x,y,u,v)=Q_{x,y}(Q^{-1})_{u,v}.                 \tag{2}
\]

The checker extracts the two matrices, verifies every entry of (2),
and checks that they are inverses over `F2`. Applying `Q^(-1)` on `x`
and `Q^T` on `v` gives a constant physical basis change. It preserves
the strict crossing endpoint and changes (2) to

\[
 S(0;x,y,u,v)=[x=y][u=v].                           \tag{3}
\]

For example, the four-letter matrix is

\[
 Q=\begin{pmatrix}
 0&0&0&1\\0&1&1&0\\0&1&0&0\\1&0&0&1
 \end{pmatrix}.
\]

The sixteen-letter matrix and its inverse are retained in the exact
evidence. The restoration map is the exact inverse constant physical
basis change; the search can therefore return to the original physical
coordinates before optimizing gauges.

Choose distinct physical letters `i,j`, let `E=E_(i,j)`, and put

\[
 A(w)=I+c w^{-k}E,\qquad c\in F_4^*,\quad k\ge1.
 \tag{4}
\]

Because `E^2=0` and the characteristic is two, `A^(-1)=A`.
Act with `A` on `x,u` and with `A^(-T)` on `y,v`. This action preserves
both pairings (3) and `[x=v][y=u]` exactly. Preservation of (3) cancels
the poles contributed by the leading coefficient. Higher coefficients
still need the complete test below.

## 3. Test the complete Laurent polynomial

Expand the four commuting physical-leg actions in an auxiliary scalar:

\[
 (I+\epsilon cE)_x(I+\epsilon cE^T)_y
 (I+\epsilon cE)_u(I+\epsilon cE^T)_v
   =\sum_{r=0}^{4}\epsilon^r D_r.
 \tag{5}
\]

For `S(w)=sum_(a=0)^D S_a w^a`, substitution of `epsilon=w^(-k)`
gives the transformed tensor. Its coefficient at a negative exponent
`e` is

\[
 \sum_{r=0}^{4}D_r S_{e+rk},                       \tag{6}
\]

where coefficients outside the original degree range are zero. The
tensor is polynomial if and only if (6) vanishes for every
`-4k<=e<0`.

The implementation keeps the five coefficient tensors in (5), applies
one common shift `w^(4k)`, adds every term, and only then tests
divisibility by `w^(4k)`. It does not require either paired alphabet's
intermediate tensor to be polynomial. That intermediate requirement
would discard the entire intended construction.

All introduced powers are nonpositive, so the degree cannot increase.
For a nonzero leading coefficient and `k>0`, only the `r=0` term
contributes at degree `D`; the degree is exactly `D`. At `w=1`, the
four maps are finite paired dual matrices, so the endpoint is the
original strict crossing.

Each physical matrix has determinant one. Therefore the determinant
of either complete cost cut is unchanged as a rational function.
At infinity all four frames are analytic units, so the ordinary
unweighted Smith profiles there are unchanged as well. Physical
gauges can interact with the new frames; their complete weighted
profiles must be recomputed.

For these reflection seeds the column cut at zero is invertible.
Its determinant remains a unit after every admissible coupled frame.
Consequently no positive common power of `w` can appear in the final
polynomial tensor. The recorded degree preservation does not hide a
discarded scalar factor.

## 4. Exact finite families

For scalar `c=1`, the following complete searches within the displayed
pole ranges pass (6):

| Alphabet | Reflection interval | Degree | Tested pole orders | Generators with a valid pole | Valid frames |
|---:|---:|---:|:---:|---:|---:|
| 4 | 4 | 24 | `1<=k<=8` | 12 of 12 | 34 |
| 16 | 8 | 96 | `1<=k<=4` | 240 of 240 | 688 |

Every generator in both rows admits `k=1,2`. In the four-letter case,
the generators `(i,j)=(0,3),(1,2)` also admit `k=3,4`, and `(3,0)`
admits all `k=1,...,8`. These are finite admissibility statements for
the stated fields, seeds, and pole ranges. They do not assume that the
same ranges hold for arbitrary interval lengths or after another frame
has already changed the seed.

The family includes explicit tensors outside the even Gaussian locus.
After restoring the original four-letter coordinates, the frame
`(i,j,k,c)=(0,1,1,1)` has vacuum coefficient `w^8` and coefficient
`w^7+w^9` at the physical bit index one. Thus it has a nonzero odd
coefficient as well as a nonzero even coefficient. The sixteen-letter
frame `(0,8,1,1)` similarly has vacuum `w^32` and odd coefficient
`w^31+w^33` at index one. This distinguishes the family from the
even Gaussian constructions; it does not imply a pressure improvement.

The catalogue also supports ordinary single-pair polynomial frames at
infinity. A separate exact scan finds 192 scalar-labelled monomial
frames for the four-letter seed and 15,696 for the sixteen-letter seed
whose full degree does not increase. Their maximum positive powers
are eight and 32. The single-pair zero obstruction explains why that
scan contains no negative powers before the two pairs are coupled.

## 5. Independent checks and search boundary

The coupled checker independently applies the four physical-leg maps
to the globally shifted polynomial tensor. It compares every resulting
coefficient with the auxiliary-scalar expansion. It also applies the
same involution twice, restores the original physical coordinates,
and checks all endpoint entries. The selected cases include nonbinary
`F4` scalars.

The controls cover 526,336 independent coefficient comparisons,
526,336 involution entries, and 526,336 endpoint entries. Sixteen
complete weighted profiles of four-letter examples agree with the
independent fraction-free `F16` local oracle at rational gauges.
Actual pole examples are rejected by both implementations. The
infinity comparison restores the degree shift when comparing the
homogeneous exact oracle with the raw-infinity numerical convention.

The source and evidence are:

- [`reflection_laurent_physical_frames.py`](../tools/reflection_laurent_physical_frames.py)
- [`check_reflection_coupled_frames.py`](../tools/check_reflection_coupled_frames.py)
- [`reflection-coupled-physical-frames.json`](../evidence/reflection-coupled-physical-frames.json)
- [`check_reflection_laurent_frames.py`](../tools/check_reflection_laurent_frames.py)
- [`reflection-laurent-physical-frames.json`](../evidence/reflection-laurent-physical-frames.json)

```sh
python3 tools/check_reflection_coupled_frames.py \
  --output evidence/reflection-coupled-physical-frames.json
python3 tools/check_reflection_laurent_frames.py \
  --output evidence/reflection-laurent-physical-frames.json
```

The bounded sixteen-letter search is implemented in
[`search_reflection16_coupled_frames.py`](../tools/search_reflection16_coupled_frames.py).
It alternates balanced and asymmetric supporting weights and retains
the actual polynomial degree, complete physical mean shifts, and
common paired-type entropy. Its parameterization `h=tH` bounds the
height-to-scale ratios used by the native oracle. A strict numerical
supporting value below one would require independent exact replay and
the supporting-line transfer before any improved factoring claim.
