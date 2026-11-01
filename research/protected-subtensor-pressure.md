# Protected scalar subtensors give channel-pressure lower bounds

## Result and scope

A large scalar part of a crossing can prevent an auxiliary deformation
from improving its regular-channel certificate. The obstruction concerns
the actual tensor on a protected rectangle, so it survives every choice
of channel basis, state precision, and block length.

For the four-letter reflection seed, consider arbitrary auxiliary
deformations supported only on the physical letter `x=3`. The other
three `x` letters must retain the scalar seed times the auxiliary
identity. Auxiliary dimension and Laurent degree are unrestricted.
At the seed's natural physical gauges, every full-word regular-channel
certificate in this family has common-type comparison parameter

```math
\boxed{\quad \rho\ge\frac8{11}>\frac1{\sqrt2}.\quad}
```

This excludes the family from beating the half-exponent threshold
through the current pressure transfer. It does not lower-bound actual
coefficient ranks after direct type restriction, and it does not
exclude another physical gauge or a method that uses more information
than full-word regular-channel pressure.

The same obstruction also excludes a deformation supported only on
either other nonzero letter: the exact lower bounds are `8/11` for
`x=1` and `5/7` for `x=2`. Thus every single nonzero-letter deformation
is excluded at these gauges, with no restriction on its auxiliary
dimension, degree, or block length. The vacuum letter `x=0` is not
covered by this conclusion.

The result also gives a general screening principle for other
matrix-valued crossings: look for a regular physical projection under
which the auxiliary terminal vector spans a scalar invariant line.

## 1. A protected scalar projection

Let `K(tau;a,b)` be a `D`-by-`D` auxiliary matrix with physical row
label `a` and column label `b`. Let `f,e` be constant auxiliary
boundaries with `f^T e=1`. The ordered word tensor is

```math
\Gamma_t(\mathbf a,\mathbf b)
=f^T K(a_1,b_1)\cdots K(a_t,b_t)e.
```

Work over a local power-series ring `O`. Suppose there are regular
physical matrices `P,Q` and a scalar matrix `S` such that

```math
\left(\sum_{a,b}P_{i,a}K(a,b)Q_{b,j}\right)e
=S_{i,j}e
\qquad\text{for every }i,j.                       \tag{1}
```

Thus the projected local action preserves the terminal line and acts
on it by `S`. No condition on its other auxiliary columns is needed.
Project each site separately and expand the finite physical sums.
Starting at the rightmost auxiliary factor, (1) gives

```math
P^{\otimes t}\Gamma_t Q^{\otimes t}=S^{\otimes t}.
                                                               \tag{2}
```

The physical projections may differ at the two endpoints. A common
special case is a protected rectangle of physical labels on which
`K(a,b)=S(a,b)I_D`. Then `P,Q` are ordinary row and column selections.

For a regular matrix `M`, let its finite Smith orders be `nu_j` and
put

```math
\Phi_s(M)=\sum_j s^{-\nu_j},\qquad s>1.
```

Regular left and right multiplication cannot increase `Phi_s`.
Indeed, the binary rank after reduction modulo `tau^p` is
`sum_j max(p-nu_j,0)` and cannot increase under those multiplications.
The positive geometric summation in
[the triangular-fiber pressure argument](triangular-auxiliary-row-fiber.md#3-consequence-for-regular-channel-pressure)
then gives the assertion. The same argument works over any coefficient
field, with its ordinary vector-space dimension in place of binary rank.

By (2),

```math
\Phi_s(\Gamma_t)\ge\Phi_s(S^{\otimes t})
=\Phi_s(S)^t.                                     \tag{3}
```

The last equality follows by tensoring regular Smith factorizations:
the resulting orders are all sums of `t` local orders.

Every decomposition of `Gamma_t` into regular rank-one terms with
charged orders has weighted count at least `Phi_s(Gamma_t)`. Therefore
every exponential rate for such word counts is at least `Phi_s(S)`.
The projection bound is independent of how the decomposition was found.
It includes adaptive policies, complete DVR basis changes, larger
state sets, and refactorization of arbitrarily long blocks.

All statements remain valid for a Laurent presentation after a uniform
scalar shift makes it regular. The shift must then be removed from
both the upper count and the protected lower bound.

## 2. Protecting three reflection letters

Let `T(z;x,y,u,v)` be the degree-twelve reflection seed used in
[the row-fiber classification](triangular-auxiliary-row-fiber.md).
Its alphabet is `{0,1,2,3}`, its physical height vector is

```math
h=(0,3,1,4),
```

and its endpoint is the ordinary crossing. Consider any finite
auxiliary dimension `D` and any matrix-valued Laurent crossing
satisfying

```math
K(z;x,y,u,v)=T(z;x,y,u,v)I_D
\qquad\text{whenever }x\in\{0,1,2\}.             \tag{4}
```

For `x=3`, the entries of `K` can be arbitrary Laurent polynomials.
They may have arbitrary auxiliary couplings and arbitrarily large
degree. The usual endpoint and boundary conditions are required if
the object is to be used as a crossing, but the protected lower bound
itself needs only (4) and `f^T e=1`.

At zero, use the cut `(x,y)|(u,v)` with height
`h(x)+h(y)`. At infinity use `(x,u)|(y,v)` with height
`h(u)+h(v)+4`. The extra common infinity power is the same
regularization as in the preceding note.

In both cuts, selecting rows with `x in {0,1,2}` gives a protected
scalar matrix of shape `12` by `16`. Its exact Smith orders are

```math
\begin{aligned}
\mu_\infty&=(0,1,3,4,5,6,7,8,8,9,10,11),\\
\mu_0&=4+\mu_\infty.
\end{aligned}                                                   \tag{5}
```

Put `r=1/s` and

```math
g(r)=(1+r)(1+r^3)(1+r^5+r^7).
```

The two protected scalar rates are `s^-4 g(1/s)` and `g(1/s)`.
Undoing the extra infinity shift multiplies their product by `s^4`.
Consequently every full-word regular-channel certificate satisfies

```math
\boxed{\quad B(s)\ge g(1/s)^2.\quad}             \tag{6}
```

This is stronger than a finite search over deformations: (6) holds for
every choice in the unrestricted auxiliary family (4).

## 3. An exact all-scale comparison with the target

The natural common-type energies of the four letters are

```math
G=(0,12,4,16).
```

For a raw product `B(s)`, the standard same-type pressure comparison
has parameter `rho` determined by

```math
\log B(s)=4\rho\log\left(\sum_{a=0}^3
s^{-G(a)/(4\rho)}\right).                         \tag{7}
```

Equivalently, the type entropy calculation minimizes
`[log B(s)+sum_a pi_a G(a) log s]/[4 H(pi)]` over nontrivial types.
The right-hand side of (7) is strictly increasing in `rho`, since
its derivative is four times the entropy of the associated Gibbs
distribution. Thus it suffices to compare at `rho_0=8/11`.

At that rational parameter, (6) is strictly greater than the
right-hand side of (7) if

```math
g(r)^{11}>
(1+r^{11/8})^{16}(1+r^{33/8})^{16},
\qquad 0<r<1.                                    \tag{8}
```

Set `r=t^8`. Define the integer polynomial

```math
p(t)=\frac{
g(t^8)^{11}-(1+t^{11})^{16}(1+t^{33})^{16}}
{t^8}.                                           \tag{9}
```

The numerator is divisible by `t^8`, and `p` has degree `960`.
Expand it in the degree-960 Bernstein basis:

```math
p(t)=\sum_{k=0}^{960}b_k\binom{960}{k}
t^k(1-t)^{960-k}.                                \tag{10}
```

The exact checker proves

```math
b_k>8\qquad(0\le k\le960).                       \tag{11}
```

This uses integer arithmetic only. If `p(t)=sum_j c_j t^j`, the
checked integer numerators are

```math
960!\,b_k=
\sum_{j\le k}c_j\binom{k}{j}j!(960-j)!.
```

The Bernstein basis is nonnegative and sums to one on `[0,1]`.
Therefore `p(t)>8` throughout that interval, proving (8) for every
positive finite scale. There is no numerical grid, unproved tail
estimate, or subdivision gap.

At each scale the comparison parameter is strictly greater than
`8/11`; taking the infimum over scales gives `rho>=8/11`. Finally,

```math
\left(\frac8{11}\right)^2=\frac{64}{121}>\frac12
```

because `128>121`. The entire family remains on the wrong side of
the half-exponent threshold under this pressure transfer.

### The other two nonzero letters

Deleting `x=1` or `x=2` gives the same difference of four between
the zero and regularized infinity profiles. Their infinity profiles are

| Deformed letter | Protected infinity orders | Profile polynomial |
| --- | --- | --- |
| `1` | `0,1,3,4,5,6,7,8,8,9,12,13` | `(1+r)(1+r^5)(1+r^3+r^7)` |
| `2` | `0,1,3,4,5,7,8,8,10,11,12,15` | `(1+r^3)(1+r^7)(1+r+r^5)` |

For letter `1`, use the same `rho_0=8/11`, substitution `r=t^8`,
and powers `11` and `16` as in (9). After division by `t^8`, the
difference polynomial has degree `1136`; all `1137` Bernstein
coefficients are strictly greater than `8`.

For letter `2`, use `rho_0=5/7`. The required exact inequality is

```math
g_2(r)^7>(1+r^{7/5})^{10}(1+r^{21/5})^{10}.
```

With `r=t^5`, divide the difference by `t^5`. The resulting degree-520
polynomial has all `521` Bernstein coefficients strictly greater than
`4`. The same geometric argument proves the all-scale comparison.
Here `(5/7)^2>1/2` follows from `50>49`.

These calculations prove the stated obstruction for every single
nonzero letter. They make no assertion for a vacuum-only deformation
or one whose support spans several physical letters.

## 4. Practical meaning and limitations

The argument identifies an obstruction that a useful next deformation
must remove. Changing only one nonzero physical `x` letter preserves
a large scalar rectangular tensor. Its weighted Smith mass alone is
already too expensive. Larger auxiliary matrices and longer blocks
cannot alter this protected tensor. A next deformation at these gauges
must reach the vacuum sector or more than one physical letter.

The general projection test (1) is more flexible than deleting a
letter. It may find a regular linear combination of physical rows
and columns on which the terminal auxiliary line remains scalar,
even when no entrywise protected rectangle is apparent.

### Computing the largest regular left projection

The exact engine
[`tools/scalar_terminal_projection.py`](../tools/scalar_terminal_projection.py)
implements the case `Q=I`. Put `S_*=f^TKe` and stack the scalar
physical residual matrices

```math
K e-eS_*.
```

A physical row vector has the desired scalar property exactly when
it annihilates all these residual matrices. Transpose the stack and
factor it by minimum-valuation Schur elimination. Its regular right
factor has a unit maximal minor. Solve on that minor, taking the
remaining coordinates freely in `O`. This gives an `O`-basis `P`
of the whole saturated kernel, with regular rational entries.

Every other regular left projection with property (1) factors as
`W P` with `W` regular. Therefore its scalar matrix is `W(P S_*)`
and its weighted Smith mass cannot exceed `Phi_s(P S_*)`. The engine
thus computes the strongest lower bound available from a regular
left projection with this fixed terminal line; it does not merely
sample projection bases. A zero-dimensional kernel is an inconclusive
screen, not evidence that a crossing is competitive.

The separate checker exercises non-polynomial regular projection
coordinates and all `144` entries of a literal two-site projection
identity. A noncommuting reflection control also changes an upper
auxiliary entry at the vacuum letter, while its lower entry changes
only at `x=3`. There is no full matrix-valued scalar rectangle of
the previous form, but the protected terminal line still returns
the same twelve-row lower profile. This demonstrates the additional
scope of the projection test.

Run this engine's controls with

```sh
python3 tools/check_scalar_terminal_projections.py
```

The evidence is
[`evidence/scalar-terminal-projections.json`](../evidence/scalar-terminal-projections.json).

The [bounded boundary extraction extension](boundary-extracted-scalar-pressure.md)
also permits reachable and observable auxiliary lines other than the
actual terminal. It gives a uniform exponential-rate lower bound even
when the extraction uses suffix words of different bounded lengths.

The use of the **full-word** pressure interface matters. The protected
submatrix consists of words with no `x=3` letters. It need not survive
a direct restriction to a common type that assigns positive mass to
letter `3`. The present theorem therefore does not rule out a new
algorithm that first restricts to types and then obtains a sharper
rank decomposition. It rules out improvement through the existing
full-word channel bound followed by the same physical-gauge/type
transfer. Other gauges are also outside the stated obstruction.

## 5. Reproduction

```sh
python3 tools/check_protected_subtensor_pressure.py
```

The standard-library checker reconstructs all six rectangular Smith
factorizations, checks every local matrix entry, constructs (9) and
its two analogues from their factored formulas, and verifies all
`2,619` strict Bernstein inequalities. It records hashes of the
generated integer coefficient lists so that the polynomials are
identified without large evidence tables.

Evidence is in
[`evidence/protected-subtensor-pressure.json`](../evidence/protected-subtensor-pressure.json).
The theorem does not improve the currently proved factoring exponent.
