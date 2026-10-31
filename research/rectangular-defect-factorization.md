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

They do not assert a reduction for arbitrary physical letter weights.
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
This optional replay uses NumPy:

```sh
python tools/check_rectangular_defect_fock.py \
  --output evidence/rectangular-defect-fock.json
```

The [recorded finite-field checks](../evidence/rectangular-defect-fock.json)
are an implementation control. They do not infer a rational-function
identity merely from its values at four points.

The gauge convention matters when comparing this result with the
rectangular native search tool. Here the mode weights in (1) are
absolute. That tool fixes each physical letter-zero weight to zero.
The sum of the four removed reference-letter weights is `2-2M` at
each place. Its normalized orders are therefore the absolute orders
from (3) plus `2M-2`, uniformly. The rational factors in the defect
construction already include their scalar denominators; no second
denominator restoration should be applied.

This support certificate does not itself record the defect Smith
multisets. A complete numerical profile implementation must still be
compared with direct full matrices in feasible cases and must retain
the stated generic rank and precision conditions. The operator algebra
and support argument above are the proof of the reduction; a finite
comparison is an implementation check, not a substitute for it.
