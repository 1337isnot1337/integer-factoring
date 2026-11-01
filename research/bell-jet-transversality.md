# Opposite Bell jets have a unique simultaneous eighth-order carrier

## Result

In the degree-sixteen seventh-order Bell construction, there is exactly
one constant trace carrier for which all fifteen nonunit Smith channels
have order at least eight at both designated finite places. This holds
for every pair of invertible four-by-four Bell matrices, including
nonsymmetric relative forms.

The carrier is the known eighth-power quadratic path. Any other carrier
introduces a seventh-order channel at one or both designated cuts. This
is an exact statement about this interpolation family. It does not
exclude improvements that trade earlier channels against later ones,
change physical gauges, or use additional singular places.

## 1. A tangent-space lemma

Let `q>=3`, over any field. Let two four-leg tensors be products of
invertible `q` by `q` matrices along the pairings `(x,y)|(u,v)` and
`(x,u)|(y,v)`. The tangent spaces to their respective matrix rank-one
varieties have zero intersection.

Independent invertible changes on the four legs reduce the first
tensor's two factors to identities. Further changes preserving those
two identities reduce the second tensor's factors to `I,E`, where
`E` remains invertible. These transformations preserve membership in
both tangent spaces.

A tensor in the intersection therefore satisfies

```math
\delta_{xy}H_{uv}+K_{xy}\delta_{uv}
  =\delta_{xu}L_{yv}+M_{xu}E_{yv}.                \tag{1}
```

Fix `x!=u` and view (1) as a matrix in `(y,v)`. The left side has
rank at most two. The right side is `M_(xu) E`, which has rank `q`
if that scalar is nonzero. Since `q>=3`, necessarily `M_(xu)=0`.

Now take `y=x` and `v!=u`; this makes every off-diagonal entry of `H`
zero. Taking `y!=x` and `v=u` makes every off-diagonal entry of `K`
zero. Taking `y=x,v=u` gives

```math
H_{uu}+K_{xx}=0\qquad(x\ne u).                   \tag{2}
```

A third index shows that all diagonal entries of `H` equal one scalar
`h`, while those of `K` equal `-h`. The left side of (1) is then
identically zero. This proves the lemma.

The alphabet restriction matters. For `q=2` and the two identity
pairings, the tangent intersection has dimension three. The exact
checker includes that boundary case.

## 2. The common seventh jet

Use the notation of the
[explicit Bell interpolation construction](entrywise-reflection-jet-interpolation.md).
Write

```math
S_V=e_1 C+\omega^2 e_\omega B_R(A)
                  +\omega e_{\omega^2}B_C(B)+K V,
\qquad
K=(w+1)^2(w^2+w+1)^7.                            \tag{3}
```

The carrier `V` has both paired traces equal to the identity.
At either finite point, divide the seventh Hasse coefficient by the
nonzero seventh coefficient of `K`. The four scalar coefficients of
`C,B_R,B_C,V` are respectively `1,omega,omega^2,1`. Thus the same
normalized seventh tensor occurs at both points:

```math
G=V+C+\omega B_R(A)+\omega^2 B_C(B).              \tag{4}
```

All coefficients of orders one through six vanish. At the row place,
the constant matrix has rank one. The coefficient of order seven in
its fifteen-by-fifteen Schur complement is the projection of `G`
modulo the row and column directions of that rank-one matrix.
Its vanishing is exactly membership of `G` in the corresponding
rank-one tangent space. The column place has the other tangent
condition.

Consequently all nonunit channels have order at least eight at both
places if and only if `G` belongs to both tangent spaces. The lemma
forces `G=0`, so the unique carrier is

```math
V_0=C+\omega B_R(A)+\omega^2 B_C(B).              \tag{5}
```

Both paired traces of (5) are `I`: the crossing `C` has trace zero
in characteristic two, and `omega+omega^2=1`.

## 3. The forced eighth-order behavior

With (5), the path is a matrix polynomial in `w^8`, of degree at most
two in that variable. Indeed, its constant jets through order seven
at the two finite points, together with its value at one, determine
each polynomial entry of degree at most sixteen uniquely. A polynomial
in `1,w^8,w^16` realizes those three values.

Thus there is no further carrier freedom within the simultaneous
eighth-order floor. For identity relative Bell form, the two eighth
Schur coefficients have rank six; for an alternating invertible
four-by-four relative form they have rank four. These finite ranks
are checked exactly in the companion control. They do not determine
all later Smith orders by themselves.

This leaves a meaningful search question: a useful new carrier must
accept at least one earlier seventh-order channel and compensate
elsewhere. The lemma supplies a precise restriction on the carrier
parameters, without assigning a pressure to that tradeoff.

## 4. Exact controls

The [checker](../tools/check_bell_jet_transversality.py) constructs both
projected tangent maps over `F4` and computes their complete ranks.
It covers alphabets two through five, including nonsymmetric Bell
factors and the alternating four-letter form. It independently checks
the normalized scalar seventh jets using integer binomial coefficients,
then replays complete finite jets and the two eighth Schur ranks for
three four-letter relative forms.

The [record](../evidence/bell-jet-transversality.json) contains exact
dimensions, ranks, and source hashes.

```sh
python3 tools/check_bell_jet_transversality.py \
  --output evidence/bell-jet-transversality.json
```

The uniform tangent-space proof is Section 1; the finite controls check
the implementation conventions and the exceptional two-letter case.
