# Binary hook residues and complete Schur central ranks

## Result and scope

The complete homogeneous Schur crossing over `F_2` has two explicit
central-channel relations at every power-of-two weight `h >= 2`.
They give the all-weight bound

```math
\operatorname{rank} A_0\le p(2h)-3.                         \tag{1}
```

Both crossing flattenings have the same central rank. Exact computations
give rank `228` at weight eight and rank `8,346` at weight sixteen.
Both values attain (1). The latter calculation uses all `231` external
Schur labels and an intermediate space of dimension `8,350`.

The new structural result is a more general hook-residue identity, proved
below by binary subset tests and bordered determinants. The finite weight
sixteen calculation extends the previous complete Schur calculation from
weight eight. It does not give an improved factoring exponent. Its central
piece alone excludes the raw weight-sixteen family as an improving finite
seed. Arbitrary linear recombinations of pieces remain outside that claim.

The predecessor construction is
[PR #30 at `28105165`, Sections 5--6](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/integral-schur-crossing.md).
The mode entries here use the proved
[charged-fermion determinant formula](schur-fermion-queries.md).
The source crossing identities are Appendix B of
[`anthropics/formal-math` at `16a5ce4`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
This proof and implementation were developed with substantial OpenAI
assistant support. They are ordinary mathematics and exact finite checks,
not Lean formalization.

## 1. Central channels and their two pairings

Let `h=2^r >= 2`, and let `V_h` have the complete Schur basis at weight `h`.
Write `M(x;u,y)` for the Schur-mode entry with parameter `x`, output `u`,
and input `y`. Its mode index is `|u|-|y|`.

All homogeneous angle coefficients vanish modulo two except at indices
zero and `2h`. Thus the central piece has intermediate space

```math
K_h=\mathbb F_2 e_0\oplus V_{2h},\qquad k=1+p(2h).           \tag{2}
```

The first coordinate represents the vacuum. Define matrices with rows
indexed by pairs of weight-`h` labels and columns indexed by `K_h`:

```math
F_{(x,y),\kappa}=M(x;\kappa,y),\qquad
G_{(x,y),\kappa}=M(x;y,\kappa).                              \tag{3}
```

In one crossing flattening the unnormalized central matrix is `F G^T`.
The exchanged factorization gives `G F^T` in the other flattening,
up to row and column permutations. These matrices have the same rank.
The homogeneous vacuum Gram matrix is invertible at power-of-two weights,
by the predecessor's integral determinant formula. Applying its inverse
to the first external parameter normalizes the crossing. That invertible
one-leg transformation preserves both ranks.

The vacuum columns are

```math
F_{(x,y),0}=B_h(x,y),\qquad G_{(x,y),0}=\delta_{xy},           \tag{4}
```

where `B_h(x,y)=M(x;empty,y)`. The second equality is the vacuum creation
identity at equal weights.

## 2. The hook-residue identity

Put `N=2h`. For `0 <= a < N`, let `eta_a` be the weight-`N` hook with
Frobenius arm `a` and leg `N-1-a`:

```math
\eta_a=(a+1,1^{N-1-a}).                                    \tag{5}
```

Choose a power of two `R` with `2 <= R <= h`. For a residue `j` modulo `R`,
define the channel vector and the external residue charge

```math
u_j^{(R)}=\sum_{a\equiv j\pmod R}s_{\eta_a},
\qquad
c_j^{(R)}(y)=
\#\{a\in\operatorname{Arms}(y):a\equiv j\}
+\#\{b\in\operatorname{Legs}(y):b\equiv -1-j\}\pmod2.
                                                                  \tag{6}
```

The vacuum coordinate of `u_j^(R)` is zero.

**Theorem.** For every pair of weight-`h` Schur labels,

```math
(F u_j^{(R)})_{(x,y)}=c_j^{(R)}(y)B_h(x,y),\qquad
(G u_j^{(R)})_{(x,y)}=c_j^{(R)}(y)\delta_{xy}.                \tag{7}
```

### Proof for the first map

The charged-fermion formula uses one row and column for the output hook,
in addition to the matrix for the vacuum output. Let the hook arm be `A`
and its leg be `B=N-1-A`. Every arm and leg of `x` or `y` is less than `h`.
In characteristic two the relevant contractions are subset tests,
disjointness tests, and equalities.

Expand the bordered determinant along the hook row and column. Its terms
have three types.

1. **Both hook endpoints meet the parameter.** For parameter arm `a`
   and leg `b`, the dependence on the hook is
   `1_(a subset A) 1_(b subset B)`. Replacing `A` by `A xor h` flips
   only the highest bit of `A` and `B`. This changes neither test, since
   `a,b < h`. It also preserves `A mod R`. The terms cancel in pairs.
2. **Both endpoints meet the input.** The equality contractions require
   an input arm and leg whose sum is `N-1`. Each is less than `h`, so
   this is impossible.
3. **One endpoint meets each block.** Fix an input arm `a`.
   Equality forces `A=a`, and the remaining contraction is precisely
   the corresponding disjointness entry of the vacuum determinant.
   Its residue indicator is `a congruent j`. Summing its cofactors
   is a Laplace expansion along that input row. It gives the vacuum
   determinant multiplied by that indicator. An input leg `b` gives
   the analogous column expansion, with `b congruent -1-j`.

Adding the third type over all input arms and legs gives the first
identity in (7). Laplace expansion is valid even when the vacuum
determinant vanishes or the two blocks have unequal sizes.

### Proof for the second map

Place the hook in the input block. Both hook endpoints meeting the
parameter give
`1_(A disjoint b) 1_(B disjoint a)`. Since `B` is the bitwise complement
of `A` within `0,...,N-1`, these are the same two subset conditions
used above. The pairing `A xor h` again cancels them.

Both endpoints meeting the output are impossible by the same index bound.
The two mixed types become Laplace expansions of the vacuum-input
determinant. They count output arms of residue `j` and output legs of
residue `-1-j`. The output is `y`, and the vacuum-input determinant is
`delta_xy`. This proves the second identity.

## 3. Two universal relations and the central rank bound

For `R=2`, both residue charges vanish for every even-weight partition.
If its Durfee rank is `t`, then

```math
c_0^{(2)}(y)
=\#\{\text{even arms}\}+\#\{\text{odd legs}\}
\equiv t+\sum a_i+\sum b_i=|y|=0\pmod2.                     \tag{8}
```

The sum of the two charges is `2t=0`, so the other charge is also zero.
Consequently `F u_0=F u_1=G u_0=G u_1=0`.

Let `U=span(u_0,u_1)`. These nonzero vectors have disjoint supports.
Each support contains `h` hooks, an even number. Under the standard
coordinate pairing on `K_h`, the space `U` is therefore totally
isotropic and has dimension two. Both row spaces lie in `U^perp`.

The restriction of the coordinate pairing to `U^perp` has radical `U`.
It has rank

```math
\dim U^\perp-\dim U=(k-2)-2=k-4=p(2h)-3.                  \tag{9}
```

Every pairing of subspaces of `U^perp` has rank at most (9). Applying this
to the two row spaces proves (1). If both row spaces equal `U^perp`,
equality follows. No finite experiment is needed for the upper bound.

For larger `R`, the charges need not vanish. For example, `(2,2)` at weight
four has nonzero arm-zero charge modulo four. The checker includes this
negative control. On any label set with a common charge vector `c`,
equation (7) instead gives the common relations `u_j^(R)+c_j e_0`.
These are exact residue constraints; no asymptotic rank saving is inferred.

## 4. Exact finite ranks

The native implementation streams rows of the two channel maps. It uses
packed binary determinants for individual modes and binary row elimination.
It never constructs the external `d^2` by `d^2` central matrix.

| Weight `h` | External dimension `p(h)` | Channel dimension | Left and right channel ranks | Both central ranks |
| --- | ---: | ---: | ---: | ---: |
| 2 | 2 | 6 | 3 | 2 |
| 4 | 5 | 23 | 16 | 12 |
| 8 | 22 | 232 | 230 | 228 |
| 16 | 231 | 8,350 | 8,348 | 8,346 |

At weights eight and sixteen, the computed orthogonal complements are
exactly the two hook-parity vectors. Thus the measured row spaces equal
`U^perp`, and the central ranks attain the theorem's upper bound.

There are two exact ways to compute the pairing rank. The direct method
forms the Gram matrix of row-space bases. The second method uses

```math
\operatorname{rank}(F G^T)
=\dim A-\dim B^\perp
  +\dim((A+B^\perp)/A),                                    \tag{10}
```

where `A` and `B` are the two row spaces. It reduces a basis of `B^perp`
modulo `A`. When the codimension is two, this needs only two reductions.
Both methods are compared through weight eight. A separate complete
weight-sixteen scan also reproduced `8,346` with the direct method.

The optional full scan examines all `53,361` external pairs at weight
sixteen. A deterministic subset of `13,000` pairs already reaches both
channel ranks `8,348` and central rank `8,346`. Together with (1), this
subset proves the full central rank exactly. It is not a statistical
estimate or a guessed rank extrapolation.

At weights four, eight, and sixteen, the exact central rank `r` satisfies
`r^2 > p(h)^3`. Hence `r > p(h)^(3/2) > p(h)^sqrt(2)`. The central term
alone exceeds the strict seed threshold for the raw decomposition.
This conclusion does not restrict other decompositions, other label
families, or general factoring algorithms.

## 5. Reproduction and remaining question

The checker needs Python 3 and a C++17 compiler. Its default run checks
weights two, four, and eight. The larger run uses the sufficient subset
at weight sixteen; add `--full-h16` to scan every pair.

```sh
python3 tools/check_schur_homogeneous_central.py
python3 tools/check_schur_homogeneous_central.py --large \
  --output evidence/schur-homogeneous-central.json
python3 tools/check_schur_homogeneous_central.py --large --full-h16
```

The checker compares `6,859` small native mode entries with the independent
Python determinant implementation. It also checks `256` larger triples,
the full vacuum Gram matrices, vacuum creation identities, every binary
hook-residue identity for each scanned pair, and the explicit parity
kernel supports. The preceding query audit independently checked the
underlying formula against the original bosonic recursion.

Implementation files:
[native mode kernel](../tools/schur_fermion_core.hpp),
[central rank program](../tools/schur_homogeneous_central.cpp), and
[verification driver](../tools/check_schur_homogeneous_central.py).
The recorded output is
[the exact central-rank evidence](../evidence/schur-homogeneous-central.json).

The unresolved structural question is whether the two channel row spaces
equal `U^perp` at every power-of-two weight `h >= 8`. The present proof
establishes the two universal relations; the computations establish
saturation only at weights eight and sixteen. A proof of saturation would
describe the central piece of this family, not improve its factoring rate.
