# A uniform obstruction on fixed balanced-charge rectangular rays

## Result and scope

The 28 entry inequalities used by the full-core leading linear programs
now have uniform algebraic certificates. They apply to every sufficiently
large power-of-two interval, with no interpolation assumption.

For the mixed rectangular seed with `alpha=1`, `beta=omega` in
`F4=F2[omega]/(omega^2+omega+1)`, the partial output repair, and the
original quartic numerator, consider the full 25-parameter family of
[ten-mode core weights](rectangular-full-core-factorization.md).
Write its 24 raw gauge weights as

\[
\theta(M)=A\tau+B,\qquad
\tau=c\log M,\qquad \ell=\log s=\tau/M.
\tag{1}
\]

Here `A` and `B` are fixed real vectors. Require balanced bulk charges
at both orders:

\[
A_0+A_1=A_2+A_3=0,\qquad
B_0+B_1=B_2+B_3=0.                                  \tag{2}
\]

**Theorem.** If `c>1/sqrt(2)`, the critical pressure of this seed is
strictly positive for all sufficiently large power-of-two `M` on
every ray (1)--(2). In particular, these rays cannot produce a crossing
exponent below `1/sqrt(2)`. The same conclusion holds for the smaller
19-parameter family, obtained by setting the six extra `x` weights
to zero.

The statement allows arbitrary fixed core-weight slopes and offsets.
It does not cover an isolated finite-size witness, coefficients `A,B`
that change with `M`, unbalanced aggregate bulk charges, a different
repair, or the full four-by-four physical alphabet. It does not prove
a global lower bound for all rectangular crossings. The complete
factoring bound remains unchanged.

The algebraic step is more limited, and more reusable, than an exact
formula for every entry valuation: 28 regular matrix coefficients have
specified **upper bounds on their leading valuation slopes**. Those
bounds suffice for the dual inequalities. The primal optimum claims
in the [fitted-model note](rectangular-leading-face-lp.md) remain claims
about that finite model.

## 1. Keep the two scales independent

Use the local variable `z` at either valuation place, and introduce an
independent formal variable

\[
Z=z^M.
\tag{3}
\]

For power-of-two `M` in characteristic two,

\[
(1+z)^M=1+Z.
\tag{4}
\]

Every initial extreme-core CAR coefficient is a finite combination of
coordinate vectors and geometric sequences

\[
z^{e j}(1+z)^{f j},\qquad 0\le j<M,
\tag{5}
\]

with integer `e,f`. Coordinates used by the concentration procedures
have the form `j=aM+b`, with fixed integers `a,b`. Their values become
monomials in `z,1+z,Z,1+Z`.

The star-shaped concentration sweeps can be performed on whole
sequences. If

\[
r(z)=z^e(1+z)^f\ne1,
\]

then

\[
\sum_{j=0}^{M-1}r(z)^j=\frac{1+r(Z)}{1+r(z)}.
\tag{6}
\]

For the constant sequence, the sum is zero because `M` is even.
Removing finitely many exceptional coordinates subtracts their explicit
values, which is addition in characteristic two. Products of sequences
are again of the form (5), with finitely many coordinate corrections.
Thus the full concentration is a fixed rational calculation in
`F2(z,Z)`, independent of the length of the interval.

The resulting active Gaussian orders, written as `(a,b)` for `aM+b`,
are

```text
infinity: (-2,1), (-2,3), (0,-3), (0,-1),
           (0,1),  (0,3), (0,5),  (2,-5), (2,-3), (2,-1)

zero:      (0,1),  (0,3), (0,5),  (2,-5), (2,-3), (2,-1),
           (0,1),  (0,3), (2,-3), (2,-1)
```

The correction and quartic operators have fixed length in the CAR
algebra. Applying them after concentration therefore gives fixed
rational matrix coefficients in `F4(z,Z)`.

This construction uses the uniform regularity theorem for the
[four- and six-coordinate concentration sweeps](rectangular-full-core-factorization.md).
The closed sequence arithmetic replaces the `M`-dependent execution
of those same sweeps; it introduces no new basis change.

## 2. A safe specialization of the small variable

Set `z=omega`, while leaving `Z` formal. In particular, neither `Z`
nor `1+Z` is replaced by a finite-field constant in the proof.

A nonzero specialized coefficient need not always give a valid
valuation bound: a denominator whose lowest `Z` coefficient vanishes
under specialization could create a pole. Here the denominator audit
is explicit.

The two geometric ratios in the concentration are

\[
s=z^2\quad\text{or}\quad s=\frac{z^2}{(1+z)^2},
\qquad
S=s^M=Z^2\quad\text{or}\quad\frac{Z^2}{(1+Z)^2}.
\tag{7}
\]

At `z=omega`, their values are `omega^2` and `omega`, respectively.
Both are nonzero and different from one. Therefore the geometric
summation denominators `1+s` survive. Every additional geometric ratio
actually used by the sequence calculation is recorded by the checker;
its specialization is also nonzero and different from one.

The concentration denominators include

\[
1+s^{M-3},\qquad 1+s^{M-5}.
\tag{8}
\]

Their lowest `Z` coefficients are one. The two four-coordinate pivots
have the closed forms

\[
A_1=\frac{1+s^{M-2}}{1+s},
\qquad
A_H=s^{-(M-2)}
 \frac{1+s^{M-2}}{(1+s)(1+s^{M-3})},\quad H=M-2.
\tag{9}
\]

For the second identity, expand the defining sum for `A_H` and use
that `M-4` is even. Both identities are also checked in `F4(Z)`.
The lowest `Z` coefficients of both pivots remain nonzero at
`z=omega`, so their inverses are safe. All other inverses are powers
of `z`, `1+z`, `Z`, `1+Z`, or products of these already audited
factors. Addition and multiplication introduce no further inverted
polynomials.

Consequently every coefficient has a presentation `P/Q` in which
specialization preserves the lowest `Z` order of `Q`. The elementary
specialization inequality is then

\[
\nu_Z f(z,Z)\le \nu_Z f(\omega,Z)
\quad\text{whenever the specialized coefficient is nonzero.}
\tag{10}
\]

Indeed, specialization can remove lowest terms from the numerator,
but it cannot create a lower numerator term, and the denominator's
order is unchanged.

Finally, for every fixed nonzero rational function `f(z,Z)`,

\[
\nu_z f(z,z^M)=aM+b
\tag{11}
\]

for all sufficiently large `M`, where `a=nu_Z f` and `b` is an
integer. To see this, clear a common denominator and compare the
finitely many `Z` powers in its numerator and denominator. Their
smallest powers eventually dominate after (3), because their nonzero
coefficients are fixed rational functions of `z`. This is the scalar
case of the [bivariate valuation theorem](bivariate-core-pressure-limit.md).

Equations (10)--(11) turn a nonzero coefficient at specialized order
`a_0` into the uniform upper bound

\[
\nu_z f(z,z^M)\le a_0M+O(1).
\tag{12}
\]

The constant and the eventual threshold are independent of the gauge
ray. An explicit common threshold is not needed for the asymptotic
statement.

## 3. Regular coefficient maps suffice

It is unnecessary to prove that the vacuum pivot coordinates chosen
at one finite size remain a basis for every larger size.

Let `c_i,a_i` be the four regular core CAR pairs on one side of a
cut, with number operators `N_i=c_i a_i` and vacuum projector

\[
P=\prod_i(1-N_i).
\]

For a core occupation mask `S` and any ordinary Fock coordinate `r`,
use the row map

\[
e_r^{\mathsf T}P\prod_{i\in S}a_i.
\tag{13}
\]

On the other side, use a column map

\[
\left(\prod_{i\in T}c_i\right)P e_c.
\tag{14}
\]

These maps are regular over the local valuation ring, whether or not
a collection of them is a basis. The CAR relations imply that (13)
and (14) intertwine every physical core-state weight with the
corresponding scalar weight of `S` or `T`.

The core operators are also homogeneous for the bulk charge. At
infinity, every active mode has charge `+1`. At zero, the first six
active modes have charge `+1` and the last four have charge `-1`.
The left core charges at zero are `(-1,-1,+1,+1)` for `y0,yT,x0,xT`;
the right ones are `(+1,+1,-1,-1)` for `v0,vT,u0,uT`.
Thus each selected row or column has an exactly computable bulk
charge and physical-state label.

If a weighted active matrix has minimum entry valuation `gamma`,
then every coefficient obtained by regular row and column maps has
valuation at least `gamma`. A selected coefficient with upper bound
(12) therefore proves

\[
\gamma_e\le a_0+f\cdot A
\tag{15}
\]

for the leading slope at that place, where `f` is its 24-component
weight feature. No inverse of a selected vacuum minor enters this
argument.

More explicitly, compare these maps `L_1,B_1` with any regular
canonical coordinate bases `L_0,B_0` available at the same `M`. The
matrices `L_1 L_0^-1` and `B_0^-1 B_1` preserve the physical core
labels and bulk charge. They therefore remain regular after arbitrary
real diagonal weights. A uniform choice of `L_0,B_0` is unnecessary.

The right CAR pairs in (14) are first pulled back through the Gaussian
diagonal. They are regular by the right-companion theorem in the
[full-core construction](rectangular-full-core-factorization.md).
The resulting coefficients are exactly the weighted core-coordinate
coefficients used by the finite model.

## 4. The 28 nonzero coefficients

The table records the left and right core masks `S,T`, the vacuum
coordinates `r,c` in (13)--(14), the certified upper slope `a_0`, and
the nonzero coefficient of `Z^{a_0}` after `z=omega`.
All masks and coordinates are ordinary nonnegative integers in the
active-mode order above. The code represents `omega` by `2` and
`1+omega` by `3`.

| Place | `S` | `T` | `r` | `c` | `a_0` | Leading coefficient |
|---|---:|---:|---:|---:|---:|---|
| Infinity | 0 | 0 | 230 | 182 | -2 | 1 |
| Infinity | 1 | 4 | 2 | 2 | -4 | `1+omega` |
| Infinity | 2 | 8 | 2 | 2 | -2 | `1+omega` |
| Infinity | 3 | 12 | 2 | 2 | -4 | `omega` |
| Infinity | 4 | 1 | 2 | 2 | -4 | `omega` |
| Infinity | 5 | 5 | 0 | 0 | -4 | `omega` |
| Infinity | 6 | 9 | 2 | 2 | -4 | 1 |
| Infinity | 7 | 13 | 0 | 0 | -4 | 1 |
| Infinity | 8 | 2 | 2 | 2 | -2 | `1+omega` |
| Infinity | 9 | 6 | 2 | 2 | -4 | `1+omega` |
| Infinity | 10 | 10 | 6 | 18 | -4 | 1 |
| Infinity | 11 | 14 | 0 | 2 | -1 | `1+omega` |
| Infinity | 11 | 14 | 6 | 2 | -3 | 1 |
| Zero | 0 | 0 | 0 | 0 | 0 | 1 |
| Zero | 1 | 4 | 0 | 0 | 0 | `omega` |
| Zero | 2 | 8 | 130 | 132 | 0 | `omega` |
| Zero | 4 | 1 | 0 | 0 | 0 | `omega` |
| Zero | 5 | 5 | 0 | 0 | 0 | 1 |
| Zero | 6 | 9 | 0 | 4 | 0 | `1+omega` |
| Zero | 6 | 9 | 130 | 128 | 0 | `1+omega` |
| Zero | 8 | 2 | 130 | 0 | 0 | `omega` |
| Zero | 9 | 6 | 130 | 0 | 0 | `1+omega` |
| Zero | 10 | 10 | 130 | 4 | 0 | `1+omega` |
| Zero | 10 | 10 | 130 | 128 | 0 | 1 |
| Zero | 12 | 3 | 130 | 0 | 0 | `1+omega` |
| Zero | 13 | 7 | 0 | 0 | 0 | `1+omega` |
| Zero | 14 | 11 | 130 | 6 | 0 | 1 |
| Zero | 14 | 11 | 0 | 256 | 0 | `1+omega` |

Three original zero-place coordinates happened to specialize to zero
at the desired order. The replacements are
`(145,529)->(145,530)`, `(320,577)->(320,578)`, and
`(529,641)->(529,642)` in the finite coordinate indexing. Each
replacement has the same full weight feature and the required slope.
The proof uses only the replacement's explicit regular maps, so it
does not require the finite indexing to remain a basis.

The [input file](../evidence/rectangular-uniform-entry-inputs.json)
contains all masks, coordinates, features, and expected coefficients.
The checker reconstructs the features from CAR homogeneity; it does
not trust a label copied from the finite fit.

## 5. Exact duals now give uniform leading bounds

Let `gamma_infinity,gamma_0` be the actual minimum entry slopes of the
weighted active matrices. Write the physical energy slopes as

\[
\begin{aligned}
E_X&=(0,-2+v_l+x_l,\ 2+v_t+x_t,\ v_b+x_b),\\
E_Y&=(0,-2+y_l+u_l,\ 2+y_t+u_t),
\end{aligned}
\tag{16}
\]

where every displayed core weight is summed over the two valuation
places. The Gaussian active minimum and physical minimum constants
cancel. Under the balanced-charge assumptions, the coefficient of
`tau` in the critical pressure is

\[
C(A)=\min E_X+\min E_Y-\gamma_0-\gamma_\infty.
\tag{17}
\]

All nonzero entry multipliers in the 17 rational dual certificates
are supported on the 28 inequalities (15). The new checker verifies
this support and all dual coefficient identities with exact rational
arithmetic. It does not use the fitted primal points or any of the
other fitted constraints to establish the lower bounds.

For each selected physical minimum, the resulting uniform lower
bounds on `C` are:

| `X` minimum | `Y` empty | `Y` low | `Y` top |
|---|---:|---:|---:|
| Empty | 2 | 2 | 4 |
| Low | 2 | 0 | 4 |
| Top | 4 | 4 | 8 |
| Both | 4 | 2 | 4 |

Hence `C>=0`, and `C=0` forces both physical minima to be low.
On that face the five gap bounds are

\[
\begin{aligned}
E_X(\mathrm{empty})-E_X(\mathrm{low})&\ge2,\\
E_X(\mathrm{top})-E_X(\mathrm{low})&\ge4,\\
E_X(\mathrm{both})-E_X(\mathrm{low})&\ge2,\\
E_Y(\mathrm{empty})-E_Y(\mathrm{low})&\ge2,\\
E_Y(\mathrm{top})-E_Y(\mathrm{low})&\ge4.
\end{aligned}
\tag{18}
\]

These are lower bounds for the actual asymptotic family. No assertion
that every bound is attained by that family is needed.

## 6. The preserved background closes the fixed-ray argument

If `C(A)>0`, the critical pressure is `C(A) tau+O(1)` plus a Gaussian
term tending to zero. It is therefore eventually positive.

If `C(A)=0`, (18) makes every non-low physical core probability
smaller by at least

\[
O(e^{-2\tau/\sqrt2})=O(M^{-\sqrt2 c}).
\tag{19}
\]

The fixed offsets `B` change the constants in this estimate, not its
power of `M`.

The [all-low background theorem](rectangular-low-background.md)
identifies the corresponding restricted tensor exactly with a
Gaussian interval of length `L=M-2`. All frozen physical weights and
the scalar produced by the restriction cancel in the pressure.
For a row/column submatrix, the sorted Smith orders are
componentwise no smaller than those of the original matrix.
Consequently the full profile partition dominates the restricted
profile partition. Only the physical partition acquires the small
correction in (19), giving

\[
F_{\mathrm{full}}\ge F_{G,L}-O(M^{-\sqrt2 c}).
\tag{20}
\]

At the critical exponent, the
[Gaussian interval expansion](quantitative-interval-frontier.md)
has infinite-product term `ell/6+O(exp(-k/ell))`. Its positive finite
tail decays like `exp(-sqrt(2)L ell)/ell`, while the negative tail
decays like `exp(-2L ell)/ell`. Since `L ell` tends to infinity,
the positive tail eventually dominates the negative one. Thus

\[
F_{G,L}\ge\ell/12
\tag{21}
\]

for all sufficiently large `M`. The condition `c>1/sqrt(2)` also
ensures that the tails vanish and that

\[
M^{-\sqrt2 c}=o(\ell),\qquad \ell=c\log M/M.
\tag{22}
\]

Combining (20)--(22) proves eventual strict positivity. This proves
the theorem within its fixed-ray and balanced-charge scope.

## 7. Reproduction and validation

Run

```sh
python tools/check_rectangular_uniform_entry_bounds.py \
  --output evidence/rectangular-uniform-entry-bounds.json
```

The [exact checker](../tools/check_rectangular_uniform_entry_bounds.py)
uses the [closed geometric core](../tools/rectangular_geometric_core.py)
and the existing baseline CAR operator program. It computes Laurent
series over `F4` with tracked absolute precision. All 28 displayed
coefficients lie strictly below the certified precision boundary.
No floating-point optimization, native Smith calculation, or fit in
`M` is used for these inequalities.

The checker also compares the closed formulas with the original
finite concentration at `M=16,32,64`, at both places and both nontrivial
`F4` points. All 3,840 final operator components have finite values
and agree. The `M=32` comparisons test the other Frobenius branch,
where `z^M=z^2` at these field points. Although intermediate
concentration pivots can have a pole there, those poles cancel in
all final operator components used by these comparisons. These are independent
implementation controls; the uniform argument comes from the closed
geometric identities and the denominator audit above.

The [evidence file](../evidence/rectangular-uniform-entry-bounds.json)
records the coefficients, precision boundaries, geometric ratios,
pivot identities, all 17 dual conclusions, and hashes of the input
and source files.
