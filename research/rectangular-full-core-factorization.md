# Ten active modes for all four extreme cores

## Result and scope

The [bounded-defect reduction](rectangular-defect-factorization.md) extends
to arbitrary weights on the extreme states of **all four physical
legs**, including `x`. For every full power-of-two interval of length
`M>=8`, ten active modes suffice at each valuation place. After the
identically zero rows and columns are removed, each active matrix is
`768` by `768`, with full generic rank.

The remaining modes are explicit scalar Gaussian spectators. The
construction permits ten independent extreme-core weights per place,
the four real bulk-charge gauges of the predecessor, and the profile
variable. Thus there is a `25`-parameter pressure search with matrix
dimensions independent of `M`.

This is a larger admissible search family. It does not establish a
pressure below `1/sqrt(2)` or a factoring score above `1/2`. The
smaller seven- and nine-mode construction remains available when the
additional `x` weights are not needed.

## 1. The additional core operators

Keep the notation `T=M-1`, `z`, `t=z^2`, `u=1+z`, and `q=u^2`
from the predecessor. At infinity the Gaussian diagonal has orders
`kappa_j=1-2M+2j`, in two consecutive blocks of length `M`.

After the original Gaussian changes of basis and transport through
the diagonal, the additional right `x_i` operators, for `i=0,T`, are

```math
\begin{aligned}
\widehat c_{x_i}
 &=u\sum_{a=0}^{T}z^{1-2M+2(a\mathbin{|}i)}c_a^{(1)}
   +u^{-1}\sum_{a\supseteq i}z^{1+2a}c_a^{(2)},\\
\widehat a_{x_i}
 &=u\sum_{a\subseteq i}z^{-1-2a}a_a^{(2)}.
\end{aligned}
\tag{1}
```

These formulas are directly `Sigma*V^-1*e_(M+i)` and
`e_(M+i)^T*V*Sigma^-1`, using the explicit regular `V` of the
predecessor. They include the previously used `x0` annihilation.

At zero, removing the regular internal `xy` Gaussian sends the
additional left `x_i` operators to

```math
c_{x_i},\qquad
a_{x_i}+u\sum_{a=0}^{T}t^{i\mathbin{\&}a}c_{y_a}.
\tag{2}
```

Together with the previous generators, these are all creation and
annihilation operators at the low and top modes of `x,y,u,v`.
Every diagonal function of an extreme physical state is a polynomial
in the corresponding two number operators.

## 2. Two regular concentration procedures

Within either block the Gaussian orders differ by
`kappa_j-kappa_i=2(j-i)`. Let `s` be any rational function of order
two, so `1+s^d` is a unit for every positive integer `d`.
Coordinate vectors at `0,T` are retained throughout.

### Four retained coordinates and one spreading covector

Suppose creation vectors lie in the span

```math
e_0,\ e_T,\ (1)_a,\ (s^a)_a,
\tag{3}
```

and annihilation covectors lie in the span of the two coordinate
covectors and either `(1)_a` or `(s^-a)_a`. Put `H=M-2` and
`S={2,...,H-1}`. Retain `{0,1,H,T}`.

For each `j` in `S` and also `j=H`, add `s^(j-1)` times row `1`
to row `j`. The sequence `s^a` vanishes at those coordinates, and
the constant sequence becomes `1+s^(j-1)`. Next, for `j` in `S`,
add

```math
\frac{1+s^{j-1}}{1+s^{H-1}}
```

times row `H` to row `j`. All creation components outside the
retained set now vanish. The first shear coefficient has order
`2(j-1)`, exactly paying for its reverse Gaussian companion; the
second coefficient is a unit with a forward regular companion.

If the spreading annihilation covector is constant, its value at
coordinate `1` is the unit

```math
A_1=1+\sum_{j=2}^{H}s^{j-1}.
```

Its values at all `j` in `S` are still one. Adding `A_1^-1`
times row `j` to row `1` removes those annihilation components.
The coefficients and their Gaussian companions are regular.

If the spreading covector is `(s^-a)_a`, its value at coordinate
`H` is

```math
A_H=s^{-H}+\sum_{j=2}^{H-1}
 \frac{1+s^{j-1}}{1+s^{H-1}}s^{-j}.
```

The first term has uniquely lowest order, so `ord(A_H)=-2H`.
Adding `(s^-j/A_H)` times row `j` to row `H` removes coordinate
`j` of that covector. This coefficient has order `2(H-j)`, again
exactly paying for the reverse companion. Neither annihilation
sweep reopens a creation coordinate: it adds rows whose creation
components have already been made zero.

This proves the four-coordinate procedure.

### Six retained coordinates and two spreading covectors

Now allow both `(1)_a` and `(s^-a)_a` in the annihilation span.
Retain

```math
\{0,1,2,M-3,M-2,T\}.
\tag{4}
```

Use the same first two creation sweeps, but only on coordinates
outside (4), and on `H=M-2` where required for their pivot.
The annihilation coordinates `2` and `H'=M-3` are unchanged by
these sweeps.

For every outside `j`, and also `j=2`, add `s^(H'-j)` times
row `j` to row `H'`. This kills `(s^-a)_a` there; the constant
covector becomes `1+s^(H'-j)`. Then, for outside `j`, add

```math
\frac{1+s^{H'-j}}{1+s^{H'-2}}
```

times row `j` to row `2`. This removes the remaining spreading
covector. The first coefficients have precisely the order required
by the reverse companions, while the second coefficients and their
forward companions are units or regular monomials. All creation
components outside (4) remain zero.

This is the predecessor's six-coordinate procedure with an extra
retained creation vector `e_T`; that coordinate is never touched.

## 3. Apply the procedures at the two places

At infinity, the first block has creation span (3), with `s=t`,
and only the constant spreading annihilation covector. Apply the
four-coordinate procedure.

In the second block, first multiply creation coordinate `a` by
`q^-a` and annihilation coordinate `a` by `q^a`. This is a regular
diagonal unit with its regular inverse companion; it leaves the
Gaussian diagonal unchanged. The two creation sequences `q^a,t^a`
become `1,s^a`, with

```math
s=t/q,\qquad \operatorname{ord}(s)=2.
```

The spreading annihilation sequences become constant and `s^-a`,
up to common nonzero scalar factors. Apply the six-coordinate
procedure. Thus the infinity active sets are

```math
J_\infty^{(1)}=\{0,1,M-2,T\},\qquad
J_\infty^{(2)}=\{0,1,2,M-3,M-2,T\}.
\tag{5}
```

At zero, the first block has both spreading annihilation covectors
and uses the six-coordinate procedure with `s=t`. The second
block has only the inverse spreading covector and uses the
four-coordinate procedure. Therefore

```math
J_0^{(1)}=\{0,1,2,M-3,M-2,T\},\qquad
J_0^{(2)}=\{0,1,M-2,T\}.
\tag{6}
```

Every new transformation and its Gaussian companion is regular.
The full Fock matrix consequently factors as one `1024` by `1024`
active map and the scalar spectators outside (5) or (6). For
`M=2,4`, keep all `2M` modes instead.

The all-length crossing endpoint implies full generic cut rank
`(3/4)*2^(2M)`. Dividing by the spectator rank gives active rank
`(3/4)*1024=768`. This rank is supplied structurally to a truncated
series routine; it is not inferred from missing later pivots.

## 4. Regular physical core coordinates

The vacuum-projector construction in the predecessor applies
without change. Both sides now have four canonical core pairs.
The pair ordering is:

| Place | Left core pairs | Right core pairs |
| --- | --- | --- |
| Infinity | `v0,vT,y0,yT` | `u0,uT,x0,xT` |
| Zero | `y0,yT,x0,xT` | `v0,vT,u0,uT` |

Right pairs are pulled back through `Sigma` before the regular
bases are constructed. They are integral because the complete
right Gaussian transformation and its inverse are regular. Pullback
does not create an inactive coordinate.

The resulting coordinate basis makes every extreme-state weight
an ordinary row or column weight. Only the `y=both` rows and
`u=both` columns are removed. Each leaves a `768` by `768` full
rank matrix. At `M=4` the corresponding size is `192`.

There are three independent extreme-state weights on each of
`x,v`, and two on each of `y,u`, after one constant per physical
leg is removed. This gives ten per place. All four bulk charges
remain diagonal in these bases by the same homogeneity proof as
before. Scalar quartic perturbations and admissible secondary-height
postcorrections add no new operators and remain in this family.

## 5. Exact support and basis controls

The [full-core support checker](../tools/check_rectangular_full_core_support.py)
exports the rational operator images and every regular change.
It checks all canonical pairings, all cross pairings, and exact
vanishing outside the retained coordinates. It also compares (1)
directly with the original one-particle `V,V^-1,Sigma` matrices.
Run:

```sh
python tools/check_rectangular_full_core_support.py \
  --output evidence/rectangular-full-core-support.json
python tools/check_rectangular_defect_core_basis.py --full-core \
  --output evidence/rectangular-full-core-basis.json
```

The [support evidence](../evidence/rectangular-full-core-support.json)
and [basis evidence](../evidence/rectangular-full-core-basis.json)
cover `M=2,4,8,16,32,64,128`. The latter checks exact integrality
and certifies the complete left and right coordinate bases by unit
vacuum minors modulo the valuation parameter.

For callers, `data_for(place,M,full_core=True)` selects this
construction. `core_basis_data(data,side)` returns the regular
generators, their ordering, the vacuum pivot choices, and the
coordinate bulk charges. The default `full_core=False` continues
to select the smaller seven-/nine-mode calculation.

The independent [full-core coordinate control](../tools/check_rectangular_full_core_coordinates.py)
starts from the complete original `M=4` tensor. At four points of
`F16`, it reconstructs both Gaussian changes and the full rational
core basis changes, then checks all `294,912` retained entries.
It also compares ten complete 192-order profiles with direct
physical matrices under arbitrary quarter-integer weights on all
four extreme cores and both bulk charges. The
[recorded controls](../evidence/rectangular-full-core-coordinates-control.json)
have no discrepancies and verify the excluded rows and columns.

The fixed-size reduction proves polynomial mathematical complexity
in the interval length, with a larger constant than the predecessor.
An implementation must separately establish that its coefficient
precision is sufficient for the particular requested matrix.
