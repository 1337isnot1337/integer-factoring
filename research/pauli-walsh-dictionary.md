# Pauli/Walsh atoms: an exact sixty-four-label dictionary exclusion

## Result and scope

The subsequent [flagged cut-budget theorem](flagged-cut-rank-budget.md)
excludes every 64-label seed under the rank-sum criterion. The independent
dictionary certificates below retain the exact Pauli support interface and
describe why these particular search spaces fail.

Over characteristic three, the following finite dictionary cannot pass the
strict rank-sum seed criterion at `d=64`:

- The indicator atoms of all joint four-symbol types on six binary positions
  in the Pauli coefficient basis.
- Their Fourier-dual orientations.
- Any rank-preserving changes of label bases, including arbitrarily many
  globally mixing binary linear label maps, and permutations of the four
  legs that preserve or exchange the two cuts.
- Nonzero scalar multiples over any extension field of characteristic three.

For every exact crossing decomposition into these atoms,

```math
\sum_i\sqrt{r_i c_i}>64^{\sqrt2}.                     \tag{1}
```

The proof uses a complete finite rank table and the
[weighted-cut theorem's](weighted-cut-rank-budget.md) all-field inequality
`sum_i min(r_i,c_i)>=d`. It does not enumerate
possible decompositions. Linear recombinations of atoms can have different
ranks and are outside this dictionary statement. No asymptotic obstruction
for every Pauli family, and no improved factoring exponent, is claimed.

## 1. The exact support-to-rank correspondence

Let `d=2^m`, and work over any field of characteristic different from two.
For binary vectors `x,z in F_2^m`, define

```math
P_{x,z}e_t=(-1)^{z\cdot t}e_{t+x},\qquad
T_{x,z}=P_{x,z}\otimes P_{x,z}.
```

All `T_(x,z)` commute: the signs in Pauli multiplication are squared.
For a coefficient function `g`, put `B_g=sum_(x,z) g(x,z)T_(x,z)`.
The vectors `vec(P_(x,z))` form an orthogonal basis with squared norm
`d`, which is nonzero. The column-cut matrix of `B_g` is their congruence
transform of `diag(g)`. Hence

```math
c(B_g)=|\operatorname{supp}g|.                        \tag{2}
```

In the row cut, the Bell vectors

```math
w_{u,v}=\sum_t(-1)^{v\cdot t}e_t\otimes e_{t+u}
```

are a basis, and the eigenvalue of `T_(x,z)` is
`(-1)^(v dot x+u dot z)`. Consequently

```math
r(B_g)=|\operatorname{supp}\widehat g|,\qquad
\widehat g(u,v)=\sum_{x,z}(-1)^{v\cdot x+u\cdot z}g(x,z).
                                                               \tag{3}
```

An ordinary Walsh transform differs only by exchanging its two output
coordinate groups. No transpose sign or quadratic factor is missing.
Swap has Bell eigenvalue `(-1)^(u dot v)`, giving the exact identity

```math
E_d=B_g,\qquad g(x,z)=d^{-1}(-1)^{x\cdot z}.          \tag{4}
```

Applying a Walsh transform to `g` interchanges its two support sizes,
since its square is `d^2` times the identity, up to the harmless coordinate
exchange. The target quadratic sign is Walsh-self-dual up to the scalar
`d`. These facts justify including the Fourier-dual dictionary.

A common binary linear label map `M` changes Pauli indices by
`(x,z)->(Mx,M^(-T)z)`. This is invertible and preserves `x dot z`.
Both support sizes in (2)–(3) are unchanged. More general invertible
four-leg changes preserving the crossing and the two cuts likewise
preserve the cut ranks; the exclusion below therefore does not depend
on a selected finite list of label bases.

## 2. All six-position type ranks over characteristic three

For a quadruple `nu=(nu_00,nu_10,nu_01,nu_11)` summing to six, let
`g_nu(x,z)` be the indicator that the six pairs `(x_j,z_j)` have that
type. There are `binom(9,3)=84` types. An exact Walsh transform over all
`4^6=4096` indices gives the following complete table. Ranks are shown
before the optional Fourier exchange.

| Row rank | Column rank | Number of types | Rank product |
|---:|---:|---:|---:|
| 4096 | 1 | 4 | 4096 |
| 2688 | 6 | 12 | 16128 |
| 2688 | 15 | 12 | 40320 |
| 2048 | 20 | 6 | 40960 |
| 2892 | 30 | 12 | 86760 |
| 2688 | 60 | 24 | 161280 |
| 1440 | 90 | 4 | 129600 |
| 3972 | 120 | 4 | 476640 |
| 720 | 180 | 6 | 129600 |

The table is finite and exhaustive for these atoms. Reduction modulo three
is done after exact integer Walsh summation. Extension of the coefficient
field does not change any rank.

## 3. A rank-sum certificate for every choice of bases

A purported witness with total cost below `64^sqrt(2)` could contain only
atoms whose individual cost is below that number. In particular every
atom of rank product at least `129600=360^2` is excluded. The five
remaining rank pairs in the table all satisfy

```math
\frac{\max(r,c)}{\min(r,c)}\ge\frac{482}{5},
```

with equality for `(2892,30)` and its reversal. Thus every remaining atom
has

```math
\sqrt{rc}\ge\sqrt{482/5}\,\min(r,c).
```

The weighted-cut budget at `s=1` applies to any actual crossing
decomposition, over any field. Summing it gives

```math
\sum_i\sqrt{r_i c_i}\ge64\sqrt{482/5}>64^{\sqrt2},     \tag{5}
```

contradicting the proposed cost. This proves (1), including all the allowed
basis changes, Fourier orientations, and scalar coefficients. It also
explains why a small-cardinality search within this dictionary is unnecessary.

The checker makes every strict comparison with integers. It uses
`577^2-2*408^2=1`, so `577/408>sqrt(2)`, and verifies the individual
atom filter and (5) by raising positive rational quantities to integer
powers. Floating-point thresholds are not proof inputs.

## 4. Exact replay and remaining construction question

```sh
python3 tools/check_pauli_walsh_dictionary.py \
  --output evidence/pauli-walsh-dictionary.json
```

The checker independently verifies the Pauli matrices, both cut conventions,
the Bell eigenvalues, and the crossing identity in several odd
characteristics at small sizes. It then reconstructs all 84 six-position
type functions, their exact Walsh supports, and the rank-sum certificate.

A preceding finite screen of eight deterministic mixed label bases found
smaller maximum rank products after adding Fourier-dual atoms. That
maximum-product improvement does not supply a finite seed: the necessary
rank sum is the quantity in (1). The certificate above closes that atom
route at sixty-four labels independently of the number of bases.

The full space of linear recombinations across incompatible label bases
remains a separate candidate, because its elements need not have ranks
from the table. The work used substantial OpenAI assistant support.

## 5. Two-atom recombinations in eight mixed bases

A second, stronger finite screen permits a piece to be a linear combination
of at most two direct type atoms from eight explicitly stored binary linear
label bases. It then includes the Fourier transform of every such piece.
The eight matrices, and every retained rank-profile witness, are in
[`pauli-walsh-pairs.json`](../evidence/pauli-walsh-pairs.json). Here the number and choice of relative bases are
fixed by that file; the arbitrary-basis claim in the first result does not
extend to cross-basis recombinations without a new check.

All coefficient ratios over characteristic-three extension fields are
covered. After normalizing the first nonzero coefficient to one, a ratio
in `F_3` is one or two and is checked directly. If the ratio lies outside
`F_3`, an entry `a+lambda b`, with `a,b in F_3`, vanishes exactly when
both coordinates vanish. Thus the support is the union of the two supports,
in both the coefficient and Walsh representations. One additional case
per pair covers every such extension ratio.

### Why a candidate piece must cost less than 295

Every nonzero Pauli tensor obeys the support uncertainty bound

```math
r(B_g)c(B_g)\ge d^2.                                  \tag{6}
```

For a field-independent proof in odd characteristic, form the translation
matrix of the nonzero coefficient function `g`. Its column rank is
`|supp(hat g)|=r`, by Fourier diagonalization, and every column has support
size `c`. A basis of `r` columns must collectively meet every row, since
no row of the translation matrix is zero. Their supports therefore cover
all `d^2` rows, proving (6).

At `d=64`, each nonzero piece consequently costs at least 64. The exact
Pell comparison proves `64^sqrt(2)<359`. A one-piece crossing has cost
4096, so a proposed passing decomposition has at least two pieces.
Each individual piece must then have cost strictly below

```math
359-64=295.                                           \tag{7}
```

The checker retains every pair rank profile with product below `295^2`.
This deliberately enlarged threshold avoids numerical threshold decisions.

A sum of one direct atom and one Fourier-dual atom is automatically too
costly. Each direct atom has row rank at least 720 and column rank at most
180, and the dual reverses these ranks. The rank difference inequality
therefore forces both ranks of their nontrivial combination to be at least
`720-180=540`. Such a piece cannot occur under (7). Two dual atoms are the
Fourier transform of two direct atoms, so their rank pairs are reversals
of those already checked.

### The exact positive dual

The complete screen checks 676,340 nonzero pair cases, including extension
ratios. Together with the single-atom cases, all 403 retained rank profiles
satisfy the exact rational inequality

```math
\frac7{500}(r+c)+\frac{17}{4}\min(r,c)\le\sqrt{rc}.    \tag{8}
```

The checker squares the nonnegative rational left side and compares it with
`rc`. This is an exhaustive finite certificate for this dictionary, not a
universal inequality for arbitrary rank pairs.

For a crossing decomposition, the two ordinary rank budgets are
`sum r_i>=4096` and `sum c_i>=4096`, while the minimum-cut budget is
`sum min(r_i,c_i)>=64`. Summing (8) gives

```math
\sum_i\sqrt{r_i c_i}\ge
\frac7{500}\cdot8192+\frac{17}{4}\cdot64
=\frac{48336}{125}=386.688>64^{\sqrt2}.                \tag{9}
```

Again the strict comparison is certified by integer powers using the upper
Pell approximation `577/408`. Thus this larger two-atom dictionary also
fails the finite seed criterion. Larger recombinations remain outside the
certificate.

Run the additional replay with:

```sh
python3 tools/check_pauli_walsh_pairs.py \
  --output evidence/pauli-walsh-pairs.json
```
