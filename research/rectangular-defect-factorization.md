# Bounded defect blocks for the full-interval rectangular crossing

## Result and precise scope

The [mixed rectangular partial repair](rectangular-f4-repair.md) on a full
power-of-two interval has an exact finite-dimensional reduction at the
two valuation places used by the interval Gaussian gauge. The number of
active fermion modes is independent of the interval length.

Let `M` be a power of two, use all modes `0,...,M-1`, and retain the
four-by-three extreme core with constants `alpha,beta` satisfying
`1+alpha*beta != 0`. Put

```math
A=4\,2^{M-2},\qquad B=3\,2^{M-2}.
```

For `M>=8`, the infinity column cut is equivalent over its discrete
valuation ring to a `128` by `128` matrix of generic rank `96`, tensored
with explicit diagonal fermion factors. The zero row cut is equivalent
to a `512` by `512` matrix of generic rank `384`, with its own explicit
diagonal factors. Zero rows and columns introduced by the coordinate
embedding are omitted when listing finite Smith orders.

These bounds concern the following physical mode weights, in boundary
order `(x,y,u,v)`:

```math
\begin{aligned}
h_\infty(a)&=(M,-M,2a-M+1,2a-M+1),\\
h_0(a)&=(-2M+1+2a,1+2a,0,0).
\end{aligned}
\tag{1}
```

Sections 6 and 7 extend these gauges by four real bulk parameters and
arbitrary weights on the retained extreme states of `y,u,v`. They do
not assert a reduction for arbitrary physical letter weights.
The [full-core extension](rectangular-full-core-factorization.md) also
permits arbitrary extreme-state weights on `x`, using ten active modes
at each place.
The weights of a mixed `top+alpha*both` core letter are those of its
leading top state. The mismatch of the other summand is included as a
Laurent coefficient of the corresponding core operator.

This is a profile-computation reduction, not a pressure certificate or a
new factoring exponent. The all-length endpoint and polynomiality come
from the predecessor. The proof below gives the bounded support and the
rank promise without inferring a rank from truncated series.

## 1. Localized operators and the profile consequence

For a one-particle map `D`, let `Gamma(D)` denote its full exterior
compound matrix, including degree zero. A polynomial Gaussian tensor
has precisely this matrix in a cut separating its two sets of variables.
Use `c_i` and `a_i` for creation and annihilation on mode `i`.
In characteristic two their canonical relations include

```math
a_i c_j+c_j a_i=[i=j],\qquad c_i^2=a_i^2=0.
```

Every linear map on a fixed set of core modes is a polynomial in their
creation and annihilation operators. In particular, on the two `u`
extreme modes the code embedding is

```math
(I+\alpha c_{u_0} n_{u_T})(I+n_{u_0}n_{u_T}),
\qquad n_i=c_i a_i,
```

where `T=M-1`. The final factor deletes the both-occupied state. The
output mixture is `I+beta*a_(y0)*n_(yT)`, followed by the projection to
the three retained states. The twelve-dimensional correction acts only
on the four extreme `v,y` modes. Finally, multiplication of the scalar
Gaussian by the quartic contributes

```math
c_{v_T}c_{y_T}\,\Gamma(D)\,a_{u_0}a_{x_0}
```

in the infinity column cut. At the zero row cut the corresponding
operators are `c_(x0)c_(yT)` on the left and `a_(u0)a_(vT)` on the
right. Each term is still confined to the stated core operator algebra.

Suppose regular invertible one-particle maps put the Gaussian into a
diagonal form `Sigma`. Transport a right creation vector `v` through
`Gamma(Sigma)` as `Sigma*v`, and a right annihilation covector `b` as
`b*Sigma^-1`. The entire transformed defect is then a polynomial in a
finite list of left creation and annihilation vectors.

If all these vectors are supported on coordinates `J`, the full map is

```math
D_{\rm defect}\ \otimes\
 \Gamma\bigl(\mathrm{diag}(z^{\kappa_i}:i\notin J)\bigr),
\tag{2}
```

after a constant permutation of modes. Here `D_defect` includes
`Gamma(Sigma_J)`. The exterior compounds of the regular transformations
and their inverses are regular invertible matrices, so (2) preserves
every finite Smith order. No arbitrary rational change of basis is
being treated as a valuation-preserving transformation.

If the finite defect orders are `nu_j`, its exact profile contribution is

```math
\sum_\nu s^{-\nu}
=\left(\sum_j s^{-\nu_j}\right)
 \prod_{i\notin J}(1+s^{-\kappa_i}).
\tag{3}
```

The product describes actual diagonal factors. It is not a bound
obtained by discarding interactions.

### An elementary operation that preserves the Gaussian diagonal

For diagonal `Sigma=diag(z^kappa_i)`, take

```math
L=I+fE_{ij},\qquad
R=\Sigma^{-1}L^{-1}\Sigma
  =I+f z^{\kappa_j-\kappa_i}E_{ij},\qquad i\ne j.
\tag{4}
```

If both displayed coefficients have nonnegative valuation, `L` and `R`
are regular units and `L*Sigma*R=Sigma`. A creation vector changes only
in coordinate `i`, by adding `f` times coordinate `j`. An annihilation
covector changes only in coordinate `j`, by adding `f` times coordinate
`i`. This distinction permits a creation sweep followed by an
annihilation sweep without reopening the spectator coordinates.

## 2. Infinity: seven active modes

Let `z=w^-1`, `t=z^2`, and `u=1+z`, so `q=1+t=u^2`. Number the two
one-particle blocks by `0,...,M-1`, followed by another copy of that
range. Write `P_ab=[b subseteq a]` and `Q=diag(q^a)`.

The Gaussian with weights (1) admits regular invertible maps

```math
U=\begin{pmatrix}0&I\\u^{-1}I&QP^{\mathsf T}\end{pmatrix},
\qquad
V=\begin{pmatrix}I&LP\\0&uP\end{pmatrix},
\tag{5}
```

where

```math
(LP)_{ab}=t^{b-a}q^{M-b}[a\subseteq b],
\qquad L_{ab}=q\,t^{(a\mathbin{|}b)-a}.
```

Their inverses are

```math
U^{-1}=\begin{pmatrix}uQP^{\mathsf T}&uI\\I&0\end{pmatrix},
\qquad
V^{-1}=\begin{pmatrix}I&L/u\\0&P/u\end{pmatrix}.
```

The full-interval Pascal identities used in the
[odd-spectrum proof](schur-one-species-interval.md) give

```math
UD_\infty V=\Sigma_\infty,
\qquad (\Sigma_\infty)_{jj}=z^{1-2M+2j},\quad 0\le j<2M.
\tag{6}
```

Every entry of (5) and both inverses is regular at `z=0`. The negative
orders belong to the diagonal (6), not to an unpaid change of basis.

The relevant output generators are creation and annihilation at `v0,vT`
and `y0,yT`. The relevant input generators are creation and annihilation
at `u0,uT`, together with annihilation at `x0`. No `xT` operator or `x`
creation is needed, because the `x` alphabet is unrestricted and the
quartic uses only `x0` annihilation.

The transformed vectors can be read directly from (5). Their only
noncoordinate components are:

- The creation vector of `yT` has second-block entries `q^a`.
- The annihilation covector of `v0` has first-block entries `u`.
- The transported annihilation covector of `u0` has second-block
  entries `z^-1*q^(M-a)`.

All remaining components lie at index `0` or `T` in one of the two
blocks. For example,

```math
\begin{aligned}
c_{y_T}&:\quad e_T^{(1)}+\sum_{a=0}^{T}q^a e_a^{(2)},\\
a_{v_0}&:\quad u\sum_{a=0}^{T}(e_a^{(1)})^T
                         +u(e_0^{(2)})^T,\\
\widehat a_{u_0}&:\quad z^{2M-1}(e_0^{(1)})^T
                +z^{-1}\sum_{a=0}^{T}q^{M-a}(e_a^{(2)})^T.
\end{aligned}
\tag{7}
```

The hat means transport through the Gaussian diagonal.

For `M>=4`, keep the coordinates

```math
J_\infty^{(1)}=\{0,1,T\},\qquad
J_\infty^{(2)}=\{0,1,T-1,T\}.
\tag{8}
```

Use the following instances of (4), in the stated order.

1. In the second block, for `j=2,...,M-3`, add `q^(j-(M-2))` times
   row `M-2` to row `j`. This removes the spreading creation vector
   outside (8). Only the annihilation coordinate `M-2` can change.
2. In the first block, add row `j` to row `1` for `j=2,...,M-2`.
   This removes the spreading annihilation covector. Creations can
   change only in the retained row `1`.
3. In the second block, add `q^(1-j)` times row `j` to row `1`, for
   `j=2,...,M-3`. This removes the remaining spreading annihilation
   covector. Its pivot at index `1` was untouched by the first sweep.

Every shear coefficient is a unit at zero. Its row index is smaller
than its column index within a block, so the companion coefficient in
(4) is also regular. Every remaining creation and annihilation
coordinate is in (8). This proves the seven-mode assertion.

The infinity spectator orders are therefore those of (6) with
first-block indices `2,...,M-2` and second-block indices `2,...,M-3`.
For `M=4`, there is one spectator, of order `-3`. For `M=2`, keep
all four modes.

## 3. Zero: nine active modes

Now let `z=w`, `t=z^2`, `u=1+z`, and `q=u^2`. In the row cut, use
row modes `(x,y)` and column modes `(v,u)`, in that order.

The two internal Gaussian factors are regular units. The `xy` factor
has coefficients

```math
f_{ab}=u t^{a\mathbin{\&}b},
```

and the `uv` factor has coefficients

```math
g_{ab}=u q^a[a\subseteq b].
```

Removing them leaves cross blocks `diag(z^(2a+1))*P^T` and
`diag(z^(2a+1))`. A constant Pascal column change puts these into

```math
\Sigma_0=\mathrm{diag}
 (z^1,z^3,\ldots,z^{2M-1},z^1,z^3,\ldots,z^{2M-1}).
\tag{9}
```

The internal Gaussian changes mix creation and annihilation. Their
effect is retained explicitly. Removing the left internal factor sends
`a_(yi)` to

```math
a_{y_i}+u\sum_a t^{a\mathbin{\&}i}c_{x_a}.
\tag{10}
```

Removing the right factor sends a `u_i` creation to itself plus
`u*sum_(a subseteq i) q^a a_(va)`. After the Pascal change, the latter
annihilation coefficients are

```math
u\sum_{a\subseteq i\mathbin{\&}j}q^a
  =u(1+q)^{i\mathbin{\&}j}
  =u t^{i\mathbin{\&}j}.
\tag{11}
```

The checker retains both components of each such operator. It does not
replace the mixed states with occupation-preserving maps.

Only the following component spans occur after transport through (9).
Scalar Laurent coefficients common to a vector are omitted in this table.

| Block | Creation components | Annihilation components |
|---|---|---|
| First (`x`) | coordinate `e0`, the sequence `1`, and the sequence `t^a` | coordinate `eT`, the sequence `1`, and the sequence `t^-a` |
| Second (`y`) | coordinates `e0,eT` | coordinates `e0,eT` and the sequence `t^-a` |

The source generators are the left `y0,yT` creation and annihilation
operators, left `x0` creation, and all four right `u0,uT,v0,vT`
creation and annihilation operators. The joint output correction is
expanded into a left `y` matrix unit and a transposed right `v` matrix
unit. This is why both right `v` operators are included.

For `M>=8`, retain

```math
\begin{aligned}
J_0^{(1)}&=\{0,1,2,M-3,M-2,M-1\},\\
J_0^{(2)}&=\{0,M-2,M-1\}.
\end{aligned}
\tag{12}
```

Write `H=M-2`, `H'=M-3`, and let `S` be the first-block indices
outside (12).

1. For `j` in `S` and also `j=H`, add `t^(j-1)` times row `1` to
   row `j`. This removes the creation sequence `t^a` there. The
   constant creation sequence becomes `1+t^(j-1)`.
2. For `j` in `S`, add
   `(1+t^(j-1))/(1+t^(H-1))` times row `H` to row `j`.
   This removes the remaining creation component. Annihilation
   components can change only at the retained coordinates `1,H`.
3. For `j` in `S` and also `j=2`, add `t^(H'-j)` times row `j` to
   row `H'`. This removes the annihilation sequence `t^-a` at those
   columns. The constant sequence becomes `1+t^(H'-j)`.
4. For `j` in `S`, add
   `(1+t^(H'-j))/(1+t^(H'-2))` times row `j` to row `2`.
   This removes the remaining annihilation component. Creations can
   change only at the retained coordinates `2,H'`.
5. In the second block, add `t^(M-2-j)` times row `j` to row `M-2`
   for `j=1,...,M-3`. This removes its spreading annihilation
   component and preserves the two coordinate components.

Every denominator in these operations is a unit. For a shear moving
from a larger row index to a smaller column index, its positive
monomial exponent exactly pays for the opposite diagonal-order
difference in (4). The other companion coefficients have nonnegative
order directly. Thus every operation and its companion is regular.

All transported generators are supported on (12), proving the
nine-mode assertion. For `M=2,4`, retain all `2M` modes instead.

## 4. The exact generic ranks

The local rectangular crossing has full generic rank `AB` in each cut.
This follows from its proved typed endpoint, even though the local
rational matrix may have a pole at that endpoint.

Let its generic rank over the rational function field be `r`. Choose
word lengths divisible by both alphabet sizes, and restrict each
crossing pair to its uniform full-letter type. The two type dimensions
satisfy

```math
D_X(m)D_Y(m)=(AB)^m/m^{O(1)}.
```

Every restricted entry is polynomial. At `w=1` it is a crossing, of
cut rank `D_X(m)D_Y(m)`. Its generic rank is therefore at least that
number, while restriction of the `m`th tensor power has rank at most
`r^m`. Letting `m` increase proves `r>=AB`; the matrix size gives
equality.

The padded full Fock map has this same rank. Its alphabet embedding
adds zero rows and columns, rather than changing its finite orders.
Since `AB=(3/4)*2^(2M)`, dividing by the full spectator rank gives

```math
r_\infty=\tfrac34\,2^7=96,\qquad
r_0=\tfrac34\,2^9=384
```

when the corresponding active sets have seven and nine modes. With
all eight modes retained at zero for `M=4`, the rank is `192`.
These are structural promises for a finite-precision profile routine;
absence of another visible pivot is not used to prove them.

## 5. Effective construction and verification boundary

The displayed sweeps have `O(M)` steps on a constant number of vectors.
Their coefficients are explicit rational functions of polynomial degree
bounded by a polynomial in `M`. After restriction to seven or nine
active modes, ordinary finite CAR matrices construct the defect map.
Clearing a common denominator and running exact polynomial elimination
on this constant-size matrix gives a polynomial-in-`M` procedure for
its finite Smith orders. Formula (3) then gives the full profile without
constructing the exponentially large physical matrix.

This statement concerns mathematical bit complexity with a large fixed
matrix-size constant. It is not a claim that a particular truncated
series implementation has enough precision for every `M`.

The [support checker](../tools/check_rectangular_defect_support.py)
uses exact binary rational functions. Run:

```sh
python tools/check_rectangular_defect_support.py \
  --output evidence/rectangular-defect-support.json
```

Its `--artifact-directory` option writes the active rational operator
vectors and the exact shear sequence. The
[recorded support evidence](../evidence/rectangular-defect-support.json)
checks `M=2,4,8,16,32,64,128`. It includes the complete two-sided
Gaussian diagonalization and inverse regularity through `M=16`, the
Pascal identity (11), every canonical and cross pairing before and
after the sweeps, and exact vanishing of all inactive coordinates.

An independent [full Fock control](../tools/check_rectangular_defect_fock.py)
starts from the original rational seed, embeds its entire `192` by `192`
physical matrix into the `256` by `256` Fock matrix at `M=4`, and forms
the exterior compounds of the full changes of basis. At each of four
points in `F16`, it checks all 65,536 entries of the spectator identity:
both off-diagonal blocks vanish, and the second diagonal block is
`z^-3` times the first. The specialized ranks are `192` and `96`.
The four points avoid the poles and have nontrivial spectator scalars.
An additional entrywise comparison checks the compiled active operator
program at both cuts, after restoring its recorded unit denominator.
This optional replay uses NumPy and, with `--compare-active`, C++17:

```sh
python tools/check_rectangular_defect_fock.py --compare-active \
  --output evidence/rectangular-defect-fock.json
```

The [recorded finite-field checks](../evidence/rectangular-defect-fock.json)
compare 16,384 active entries at infinity and 65,536 at zero at each
point. At zero, two specializations have rank `160` and two have rank
`192`; this does not alter the generic-rank theorem. These comparisons
are implementation controls. They do not infer a rational-function
identity merely from its values at four points.

The gauge convention matters when comparing this result with the
rectangular native search tool. Here the mode weights in (1) are
absolute. That tool fixes each physical letter-zero weight to zero.
The sum of the four removed reference-letter weights is `2-2M` at
each place. Its normalized orders are therefore the absolute orders
from (3) plus `2M-2`, uniformly. The rational factors in the defect
construction already include their scalar denominators; no second
denominator restoration should be applied.

The [active profile tool](../tools/rectangular_defect_profiles.py) builds
the finite CAR program and computes its orders with a supplied generic
rank. It must retain the stated precision conditions. The operator
algebra and support argument above are the proof of the reduction; a
finite comparison is an implementation check, not a substitute for it.

## 6. Bulk charges and arbitrary extreme-state weights

The same support reduction allows two additional real bulk parameters
at each place. These statements concern weighted Smith orders in the
sense used by the rectangular profile oracle.

At infinity, add `lambda` to the weight of every occupied `y,v` mode,
and `mu` to every occupied `x,u` mode. Every change of basis in the
infinity proof is an exterior compound of a one-particle map. It
therefore preserves total occupation. The new active matrix is the
old one with row weight `lambda*|mask|` and column weight
`mu*|mask|`. Each spectator order changes from `kappa` to
`kappa+lambda+mu`.

At zero, the corresponding conserved quantity is the difference
between occupation in the first and second blocks. Add physical
mode weights

```math
(\lambda,-\lambda,-\mu,\mu)
\quad\text{in boundary order }(x,y,u,v).
\tag{13}
```

Every `xy` creation pair and `uv` annihilation pair preserves this
difference. The Pascal operation and the support sweeps stay within
their respective blocks. They therefore commute with the two diagonal
weight maps in (13). Put
`q(mask)=|mask in first block|-|mask in second block|`.
The active row and column weights are `lambda*q(mask)` and
`mu*q(mask)`. A first-block spectator order gains `lambda+mu`, while
a second-block spectator order loses that amount. Neither argument
requires the two parameters to be integral or nonnegative.

There is a further finite family of arbitrary real weights. Keep the
natural bulk dependence, and add any function of the extreme core
state on each of `y,u,v`. There are three retained states on `y,u`
and four on `v`. Removing one additive constant per leg leaves seven
independent core parameters at each place. A diagonal function of two
occupation bits is a polynomial in their two number operators. The
listed core CAR generators already contain these operators, so the
seven- and nine-mode support bounds remain unchanged.

This observation alone would require inserting the extra diagonal
operators inside the active matrix program. The next section gives
fixed regular bases in which every such weight is an ordinary row or
column weight. These bases are independent of the real weight values.

## 7. Regular coordinate bases for the core operators

Let `R` be the local valuation ring, and let `r` canonical creation and
annihilation pairs act on the `n=2^k` dimensional active Fock module.
Assume all their coefficients are in `R`. Define the vacuum projector

```math
P=\prod_{i=1}^{r}(I+c_i a_i),\qquad d=2^{k-r}.
\tag{14}
```

The canonical relations give a direct sum decomposition into the
`2^r` modules `c_S image(P)`, all of rank `d`, and

```math
P a_S c_T P=[S=T]P.
\tag{15}
```

Choose `d` independent rows `I` and `d` independent columns `J` of
`P` modulo the maximal ideal. The minor `P[I,J]` is invertible modulo
that ideal. Form square matrices by block rows and columns:

```math
L_S=E_I^{\mathsf T}P a_S,
\qquad B_T=c_T P E_J.
\tag{16}
```

Equation (15) gives

```math
LB=I_{2^r}\otimes P[I,J].
\tag{17}
```

Both matrices are integral, and their product has unit determinant.
They are consequently regular invertible changes of basis. No
division by a possibly nonunit minor is involved. The same canonical
relations give the exact intertwining identities

```math
L c_i=c_i^{\rm coordinate}L,\quad
L a_i=a_i^{\rm coordinate}L,\quad
c_i B=Bc_i^{\rm coordinate},\quad
a_i B=Ba_i^{\rm coordinate}.
\tag{18}
```

Thus `L` is a valid row coordinate map and `B` is a valid column
coordinate map. The row and column sides of the rectangular defect
use different core CAR systems; apply the construction separately to
those two systems.

At infinity the left system consists of `v0,vT,y0,yT`, and the right
system of `u0,uT`. At zero the left system is `y0,yT`, and the right
system is `v0,vT,u0,uT`. The right vectors in the support artifact
have already been transported through `Sigma`. Pull them back before
applying (14): multiply each creation component at coordinate `j` by
`z^-kappa_j` and each annihilation component by `z^kappa_j`.

These pulled-back vectors are integral. On the full Fock space they
are precisely the original right core operators conjugated by the
regular right Gaussian change of basis and its regular inverse.
Diagonal pullback does not create an inactive coordinate. The left
vectors are integral for the analogous reason. This establishes the
integrality assumption of (14) at every interval size, rather than
assuming that rational canonical pairs are regular.

After these coordinate changes, the excluded `y=both` rows and
`u=both` columns are identically zero. Deleting them leaves square
matrices of size `96` at infinity and `384` at zero, with full generic
rank. At zero with `M=4`, the size is `192`. Row core weights are the
sum of the `v,y` weights at infinity and the `y` weight at zero.
Column core weights are the `u` weight at infinity and the sum of
the `v,u` weights at zero. The physical state ordering remains
`(low,vacuum,top,both)`; coordinate occupation masks are
`(1,0,2,3)` in that order.

The bulk charges of section 6 also stay diagonal in these coordinates.
For a basis index `(S,j)`, the row charge is the charge of vacuum row
`I[j]` plus the sum of the charges of the core creations in `S`.
The column charge uses vacuum column `J[j]` in the same formula.
At infinity every core creation has charge `+1`. At zero a `v`
creation has charge `+1` and a `u` or `y` creation has charge `-1`.
This follows because `P` commutes with the bulk charge, and each core
creation is homogeneous for it. The charge must be computed by this
formula; it need not equal the occupation count of the new coordinate
index.

The [core basis checker](../tools/check_rectangular_defect_core_basis.py)
checks exact rational integrality, canonical pairings, the vacuum minor,
the ranks of both complete coordinate maps modulo `z`, and every
intertwining identity modulo `z`. Its `core_basis_data(data,side)` API
exports the regular rational generators, their ordering, the pivot
sets `I,J`, and the exact coordinate bulk charges. Run:

```sh
python tools/check_rectangular_defect_core_basis.py \
  --output evidence/rectangular-defect-core-basis.json
```

The [recorded basis evidence](../evidence/rectangular-defect-core-basis.json)
contains both sides at both places for `M=2,4,8,16,32,64,128`.
The proof of (14)--(18) supplies the identities over the valuation
ring; the finite residue calculations select and certify unit minors
for the explicit implementation.

An independent [coordinate audit](../tools/check_rectangular_defect_coordinates.py)
starts from the original full `M=4` tensor, constructs the regular
Gaussian and core basis changes, and compares every retained entry
at four `F16` points. It also compares complete weighted profiles
with direct physical matrices under arbitrary quarter-integer core
and bulk-charge weights. The
[recorded coordinate controls](../evidence/rectangular-defect-coordinates-control.json)
contain 184,320 entry comparisons and ten complete 192-order profile
comparisons, all exact or with zero numerical discrepancy as applicable.

## 8. Scalar quartic perturbations and other core corrections

The support proof depends on the location of the quartic operators,
not on their scalar coefficient. Replace the scalar `1/(1+w)^2` by

```math
f(w)=\frac{P(w)}{(1+w)^2},\qquad P(1)=1,
\tag{19}
```

where `P` is a Laurent polynomial over the coefficient field. The
ramified core endpoint is unchanged, so the same partial correction
works. The cancellation in the predecessor's denominator argument is
unchanged as well: the potentially double-pole term becomes
`P Q10+G0`, which still vanishes at `w=1`, and `P Q11` still vanishes
there. The remaining denominator at `w=1` has order at most one.
Any pole introduced at `w=0` by a Laurent monomial can be cleared by
a global power of `w`, whose value at the crossing endpoint is one.
After the usual paired-type restriction, this gives polynomial
crossings at every word length.

No new CAR vector appears in (19). The same active dimensions,
spectators, coordinate bases, and rank promise therefore apply. For
the simple choice `P(w)=w^ell`, the quartic scalar in the two local
variables is

```math
f_\infty(z)=\frac{z^{2-\ell}}{(1+z)^2},
\qquad f_0(z)=\frac{z^\ell}{(1+z)^2}.
\tag{20}
```

The effective construction remains polynomial in `M` and the degree
and coefficient encoding of `P`. An unrestricted rapidly growing
choice of `P` is not being assigned a free complexity bound.

Likewise, replacing the twelve-dimensional output correction by any
other admissible correction on the same retained `v,y` core adds no
new operator support. In particular, a strict-height triangular
postcorrection can be included in this same fixed-size calculation.
Its endpoint and polynomiality must still satisfy the paired-type
theorem; bounded support alone does not certify admissibility or
improved pressure.
