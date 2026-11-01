# Homogeneous linear retractions of two reflection sites

## Result and scope

Constant linear combinations of equal-height physical words provide
crossings outside the coordinate-subset search. This note records a
complete finite test in one such family: all 120 rank-three retractions
of the four-dimensional middle-height space of two four-letter
reflection sites, while every other physical word is retained.

The resulting physical dimension is fifteen. The exact construction
produces 24 pairs of endpoint Smith profiles. A nonorthogonal choice
has strictly smaller pressure than a specified coordinate deletion at
every finite scale `s>=2`. It preserves the same type energies, so this
is a strict improvement within that puncture comparison.

The numerical type optimization remains weaker than the unrestricted
reflection seed. This family does not provide a stronger factoring
score. The point of the checkpoint is a reusable homogeneous
restriction constructor and an exact comparison demonstrating that
nonorthogonal retractions can improve a coordinate restriction.

## 1. Equal-height words can be mixed

Use the normalized four-letter reflection seed from
[the common-height fiber calculation](triangular-auxiliary-row-fiber.md),
whose letter heights are

$$h=(0,3,1,4).$$

Index a two-site physical word by `w=a+4b` and assign
`H(w)=h(a)+h(b)`. The height-four space has basis

$$\mathcal F=\{3,6,9,12\},$$

corresponding to the words `(3,0),(2,1),(1,2),(0,3)`.
All remaining height spaces are retained in full.

Let `P:F_2^15 -> F_2^16` and `L:F_2^16 -> F_2^15` be constant
maps satisfying `LP=I`, and require each map to preserve the word
heights. For the two-site crossing `T`, define

$$S=(L\otimes L)T(P\otimes P). \tag{1}$$

The swap endpoint is preserved because the same maps occur on both
physical legs. Homogeneity means the maps commute with the inherited
diagonal endpoint gauges. Thus the exact gauged endpoint matrices
can be restricted directly; no new Laurent expansion is inferred
from an approximate profile.

At zero, remove the inherited common scalar `z^8` from the two-site
matrix. After this restoration, the total type energy of a retained
letter of height `H` is

$$G=4H. \tag{2}$$

The ambient normalized zero matrix equals the opposite infinity
matrix entry by entry. Their compressed matrices need not be equal:
the zero cut uses the maps `L tensor L` and `P tensor P`, whereas
the opposite cut uses the vectorization map `A -> LAP` on both
sides. The checker computes both resulting profiles independently.

## 2. The finite family is complete

Write the four coordinates of `mathcal F` in the order displayed
above. A rank-three projector on that space is determined uniquely
by a nonzero kernel vector `v` and nonzero image functional `f`
with `f(v)=1`. It is

$$F=I+vf. \tag{3}$$

There are `15*8=120` such pairs. Choose a coordinate `j` with
`f_j=1`; for `i!=j`, the embedding and retraction are

$$P_i=e_i+f_i e_j,\qquad L_i=e_i^{\mathsf T}+v_i f. \tag{4}$$

These satisfy `LP=I` and `PL=F`. Adding the twelve unchanged
physical words gives (1). Every pair is enumerated, with no random
sampling of this middle-height family.

An internal change of basis in the retained fifteen-dimensional
space changes both endpoint matrices by constant invertible row and
column factors. Consequently the Smith profiles depend on the
projector, not on the convenient image basis in (4).

The finite checker resolves all 225 orders in both cuts for every
pair. It also verifies the full literal swap endpoint, the height
conditions, and the retraction identities. Total physical parity
splits each endpoint into two blocks; the ambient tensor's exact
parity identity justifies these blocks. No uncomputed cross-block
entries are assumed to vanish.

## 3. A strict pressure gain over coordinate deletion

Compare the following two projectors, using integer bit masks in
the middle-height basis:

| Choice | Kernel vector `v` | Image functional `f` |
|---|---:|---:|
| Coordinate deletion | 1 | 1 |
| Nonorthogonal retraction | 1 | 15 |

For each choice let `a_i` and `b_i` be its restored zero and
infinity Smith orders, and define the real counting polynomials

$$Z(r)=\sum_i r^{a_i},\qquad I(r)=\sum_i r^{b_i}.$$

The pressure is `Z(1/s)I(1/s)`. Its coordinate-minus-nonorthogonal
difference factors exactly as

$$Z_{\mathrm{coord}}(r)I_{\mathrm{coord}}(r)
 -Z_{\mathrm{mix}}(r)I_{\mathrm{mix}}(r)
 =r^9(1-r)^2 Q(r), \tag{5}$$

where `Q` is an explicit degree-67 integer polynomial with constant
coefficient one. The evidence contains every coefficient. On
`0<=r<=1/2`, keep its constant term, discard its other nonnegative
terms, and bound every negative term by its value in absolute size
at `r=1/2`. Exact integer arithmetic gives

$$Q(r)\ge
\frac{147573878836948166391}{147573952589676412928}>0. \tag{6}$$

Equations (5)--(6) prove a strict pressure improvement for every
finite `s>=2`. This is not an assertion that both separate cuts
improve: one profile gains while the other loses, and the exact
product comparison is essential.

All 120 retractions retain the same multiset of heights and hence
the same type-energy generating function. Therefore (5) improves
the associated common-type certificate at the same scale and trial
rate. This comparison would not be valid if the height charge were
silently changed.

For orientation only, numerical optimization of the exact profiles
gives approximately:

| Family member | Type-rate diagnostic | Scale |
|---|---:|---:|
| Coordinate deletion `(1,1)` | 0.73604482 | 2.03352 |
| Nonorthogonal choice `(1,15)` | 0.73595486 | 2.02961 |
| Unrestricted two-site seed | 0.73108094 | 1.93097 |

These display values are not used for the exact comparison. The
complete finite experiment found no improvement over the unrestricted
seed. It does not prove an all-scale exclusion for arbitrary
homogeneous restrictions or different height-space dimensions.

## 4. Implementation and replay

The engine
[`reflection_homogeneous_restrictions.py`](../tools/reflection_homogeneous_restrictions.py)
constructs the ambient endpoints directly from the repository's
four-letter seed. It accepts general homogeneous binary embeddings
and retractions, checks `LP=I`, and returns both exact Smith profiles.
The particular family (3) is exposed separately.

The checker
[`check_reflection_homogeneous_restrictions.py`](../tools/check_reflection_homogeneous_restrictions.py)
enumerates the 120 projectors and records their 24 distinct profile
pairs. Independent block-Toeplitz ranks check the first four jets
of both displayed choices, for a total of 32 controls. Invalid
projectors and a valid retraction that mixes different heights are
rejected. The polynomial identity (5) and lower bound (6) use only
integer arithmetic.

```sh
python3 tools/check_reflection_homogeneous_restrictions.py \
  --output evidence/reflection-homogeneous-restrictions.json
```

The [evidence](../evidence/reflection-homogeneous-restrictions.json)
contains all embedding and retraction maps, profiles, energies,
comparison coefficients, and source hashes. Use `--small-only` for
the displayed comparison and two further controls, without claiming
the complete 120-case enumeration.
