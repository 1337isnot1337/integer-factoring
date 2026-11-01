# The first possible rectangular cofactor jet has bounded crossed rank

The rectangular cofactor gap has a stronger prefix form. On a matched
`p`-by-`b` physical subalphabet, put `n=pb`, `h=min(p,b)^2` and
`c=n-h>0`. If the normalized crossed source profile has second order `u`,
then the crossed profile `eta` of its adjugate satisfies

```math
\eta_i\ge(c-1)u\quad(i\le c),\qquad
\eta_i\ge cu\quad(i>c).                               \tag{1}
```

This holds over every field, without a degree bound or an invertible
leading channel. It strengthens the preceding
[rectangular cofactor gap](subalphabet-adjugate-caps.md). A short exact
entropy witness consequently improves the all-scale nine-letter bare
floor from `17/24` to **`709/1000 > 1/sqrt(2)`**.

These are necessary crossing constraints and a method obstruction. The
complete factoring score stays **1/2**. A tempting further tail refinement
is false even over `F2`; its small explicit counterexample is retained.

## 1. Rank of the first possible cofactor coefficient

By symmetry assume `p<b`. For a rank-one crossed tensor, its ordinary
operator on `p`-by-`b` matrices has the form
`R(M)=A M^T B^T`, with `A` of size `p`-by-`b` and `B` of size
`b`-by-`p`. Its operator rank is at most `h=p^2`.

Expand the ordinary adjugate with an independent variable `z`:

```math
adj(R+zU)=\sum_j z^j C_j(R,U).
```

Every `(n-1)`-minor needs at least `n-1-h=c-1` factors from the
perturbation, so `C_j=0` for `j<c-1`. The new claim is

```math
rank\bigl(\operatorname{cross} C_{c-1}(R,U)\bigr)\le c. \tag{2}
```

First work with generic full-rank `A,B`. Independent changes of basis on
their four physical spaces put both into the canonical inclusion/projection
forms. Then `R` is transposition on the top `p`-by-`p` corner and zero on
the complementary `p`-by-`(b-p)` rectangle. Its nonzero operator block has
size `h`, and its kernel/cokernel have size `c=p(b-p)`.

A coefficient using exactly `h` entries of `R` can occur only in a cofactor
whose omitted row and column both belong to that complementary rectangle.
Thus `C_(c-1)` is supported on the complementary ordinary block. After
crossing, only the `c` rows and `c` columns whose corresponding physical
index lies in the complement can be nonzero. Its crossed rank is at most
`c`. Invertible changes on physical legs preserve crossed matrix rank,
proving (2) on the generic locus.

Every `(c+1)`-minor in (2) is an integer polynomial in the entries of
`A,B,U`. Its vanishing over the generic characteristic-zero field is an
identity over the integers. It therefore holds over every field and after
arbitrary specialization, including singular channel factors. This is a
polynomial-rank identity; it does not assume that individual Smith orders
are preserved under coefficient specialization.

## 2. From the rank identity to complete prefixes

Work over the local DVR. If `u>0`, a unit pivot of the crossed source
gives `T=R+t^u U` with a regular rank-one crossed lift `R`. Hence

```math
adj T=t^{(c-1)u}\{C_{c-1}(R,U)+t^u V\}.                \tag{3}
```

The crossed first term in braces has ordinary rank at most `c`, by (2).
For any regular matrix `L` of rank at most `c`, the `(c+1)`th Smith order
of `L+t^u V` is at least `u`: diagonalize `L` over the DVR and reduce
modulo `t^u`, when its image needs at most `c` cyclic generators. This
is the same elementary rank-approximation argument used for individual
Smith interlacing. Equation (3) proves (1). For `u=0` the bound is immediate.

In particular, every first-`k` prefix satisfies

```math
S_k(\eta)\ge\{k(c-1)+\max(0,k-c)\}u.                 \tag{4}
```

For all `n` orders this becomes
`S_n(eta)>=c(n-1)u`, rather than the older `n(c-1)u` bound.

## 3. Lifted graded inequalities

Use the normalization and smallest physical energy subsets in the
subalphabet note. Write

```math
W=b\sum X+p\sum Y,\quad C=(n-1)d+W,
\quad L_k=k(c-1)+\max(0,k-c).
```

Applying the adjugate cell and shared-line minors with (4) gives

```math
L_k s_2+k\sum_{i<n}r_i\le kC-TopCell_k,\quad1\le k\le n, \tag{5}
```

```math
L_l s_2+l\sum_{i\le n}r_i-Top_l(r)\le lC-TopLine_l,
\quad1\le l\le\max(p,b).                               \tag{6}
```

Use the permitted longest-line energy terms exactly as in the preceding
note. Both cut orientations hold; averaging gives valid linear conditions
on the averaged profile. The absolute coefficients on `nu_1,nu_2` are
`k(n-1)-L_k` and `L_k`, both nonnegative. All other individual absolute
coefficients are also nonnegative. Individual Smith interlacing therefore
lifts these constraints to full profiles, including changed restricted
minima and all aspect ratios.

The [implementation](../tools/rectangular_adjugate_jet_caps.py) leaves the
predecessor source unchanged, preserving its source-pinned dual receipts.

## 4. A stronger nine-letter entropy witness

In the bare normalization `z=(mu+nu)/(2M)`, the full rectangular determinant
case of (5) is especially simple. Put `N=pq`, `c=N-p^2`; then

```math
J_{p,q}(z):=c(N-1)z_2+N\sum_{i<N}z_i\le N(N-1).       \tag{7}
```

Take `q=9`. Combine the old shared-prefix, square determinant and square
shared-line rows `P_9,D_9,A_9` of the
[bare entropy proof](subalphabet-bare-nine-letter-floor.md) with (7):

| Row | Integer weight |
| --- | ---: |
| `P_9` | 44835 |
| `D_9` | 7 |
| `A_9` | 484 |
| `J_(4,9)` | 56 |
| `J_(5,9)` | 59 |
| `J_(6,9)` | 38 |
| `J_(7,9)` | 25 |
| `J_(8,9)` | 3 |

The combined right side is `W=2009439/2`. Normalize each nonfirst left
coefficient by `2W`, and complete the remaining first-coordinate probability.
The checker verifies a positive rational probability vector `pi` with
`pi dot z<=1/2`. Its entropy exceeds `(709/500)*log 9` by more than
`0.00087581`, using 100-bit outward arithmetic. The same all-scale
Cauchy--Schwarz/weighted-AMGM proof therefore gives
`rho_bare>709/1000`. The fixed auxiliary projection transfer also applies.

No numerical minimizer is part of this certificate. The
[witness](../evidence/rectangular-jet-entropy-witness.json) has eight short
integer weights; the complete probability vector is stored in the replay
packet. This strengthens the nine-letter obstruction but does not settle
the bare ten-letter interface.

## 5. A false stronger tail, already over F2

It is tempting to extend (1) to the lower profile
`((c-1)u repeated c times, cu, (c+1)u repeated h-1 times)`.
That statement fails. An explicit `2`-by-`3` degree-three strict crossing
has crossed source profile

```text
0,2,2,2,2,2
```

and crossed adjugate profile

```text
2,3,5,5,6,7.
```

Here `c=2`, `u=2`. It satisfies the proved lower profile
`(2,2,4,4,4,4)` but fails the proposed `(2,2,4,6,6,6)`.
The [small binary fixture](../evidence/rectangular-binary-tail-counterexample.json)
retains all 36 polynomial entries and cofactors. The checker recomputes
them by permutation expansion and independently reconstructs both complete
profiles by fraction-free elimination. Singular kernel jets cannot be
replaced by an argument requiring an invertible leading compressed block.

## 6. Exact controls and boundary

Run `make verify-rectangular-adjugate-jets`. The
[checker](../tools/check_rectangular_adjugate_jets.py) independently
interpolates 416 cofactors across eight rectangular sources, including
singular compressed leading blocks and jet orders one and two. It checks
the full profile bound (1), 6,168 lifted inequalities on the earlier
independently profiled literal inputs, both normalized/absolute forms,
the binary counterexample and the exact nine-letter entropy witness.
The [packet](../evidence/rectangular-adjugate-jet-controls.json.xz) retains
all inputs, complete profiles, cofactors, identities and source hashes.

The ordinary proof supplies the all-field claim; these finite controls test
it without substituting for the proof. Independently gauged energy spaces,
nonlinear compatibility, larger bare alphabets and constructive kernels
remain open. Substantial OpenAI assistant support was used. No new Lean
formalization or stronger factoring algorithm is claimed.
