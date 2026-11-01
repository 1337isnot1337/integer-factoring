# Opposite filtrations force exact energy preservation

Commuting with the integer energy operator is a sufficient condition
for preserving one complete weighted Smith profile, but a larger
constant matrix group also does so: matrices preserving its filtration.
For the five-letter reference this enlarges the matrix space from
dimension 69 to 347. An exact opposite-filtration argument shows why
the extra directions cannot retain the original infinity pole.

This is a global statement over every characteristic-two field, with
no tangent-family, finite-field or leading-rank hypothesis. It justifies
the exact energy-block restriction within this larger zero-preserving
constant gauge group. It does not exclude gauges spending a worse pole
budget or otherwise preserving the zero profile. The factoring score
remains `1/2`.

## 1. Filtration preservation keeps the whole zero profile

Use the reference `A(t)=T(t)S`, integer joint energies `H_r`, and
`D_H(t)=diag(t^H_r)` of
[joint energy conjugations](joint-energy-conjugations.md). The zero
operator is `Z(t)=D_H(t)A(t)`. Let `U` be any invertible constant
matrix supported on

```math
U_{rc}\ne0\quad\Longrightarrow\quad H_r\ge H_c.         \tag{1}
```

Its inverse preserves the same filtration. Therefore
`D_H U D_H^-1` and its inverse are polynomial matrices at zero,
with nonzero constant determinant. For `A_U=UAU^-1`,

```math
Z_U=(D_HUD_H^{-1})\,Z\,U^{-1}.                        \tag{2}
```

Both factors are local units, so (2) preserves every zero Smith order.
The strict crossing endpoint also survives since `A(1)=I`.
This argument only needs (1), not exact commutation with `D_H`.

## 2. The joint infinity operator is a local unit

Write

```math
C(z)=z^8A(z^{-1})D_H(z).
```

Its entries are regular at zero. The constant **joint operator** is
invertible: the exact source satisfies `C(0)^4=I`. Its opposite-cut
realignment has rank one. In particular `C` and its formal inverse
belong to `GL_25(K[[z]])` over every characteristic-two field.

For any constant invertible `U`, without yet assuming (1),

```math
C_U(z)=U C(z)V(z),\qquad
V(z)=D_H(z)^{-1}U^{-1}D_H(z).                         \tag{3}
```

Multiplication on the left by a matrix in `GL_n(K[[z]])` preserves
the minimum entry valuation: multiplying by its inverse proves the
reverse inequality, or reduce the lowest coefficient modulo `z`.
Apply this to the local unit `U C(z)` in (3). Since

```math
V_{rc}(z)=(U^{-1})_{rc}z^{H_c-H_r},
```

the smallest unregularized infinity order after opposite-cut
realignment is exactly

```math
\boxed{-8-\delta(U),\qquad
\delta(U)=\max\{H_r-H_c:(U^{-1})_{rc}\ne0\}.}           \tag{4}
```

Realignment only permutes entries, and the smallest Smith order is
the minimum entry valuation, so it preserves this conclusion.
Also `delta(U)>=0`: `det V=det U^-1` has valuation zero, which would
be impossible if all entries had positive valuation. Thus the
reference's smallest order cannot be increased by any constant
conjugation.

## 3. Intersecting the two filtrations

Retaining the original pole means `delta(U)=0`. Equation (4) says
that `U^-1` must be supported where `H_r<=H_c`, the opposite
filtration. If `U` also obeys (1), its inverse obeys both support
conditions. Hence

```math
(U^{-1})_{rc}=0\quad\text{unless }H_r=H_c,
```

and `U` preserves the exact energy blocks. Conversely an energy-block
`U` has `delta(U)=0` and preserves the full zero profile by (2).
Therefore the pole-preserving subgroup of the 347-dimensional
filtered group is **exactly** the 69-dimensional energy-block group.
No separate component of that filtered group evades the argument.

The condition concerns the minimum infinity order. A larger leading
realignment rank or different later orders are not independently
excluded by this lemma. Gauges preserving the zero profile for a
source-specific reason outside (1), or spending an additional pole,
remain distinct possibilities.

## 4. Exact replay

The [checker](../tools/check_joint_filtered_gauges.py) verifies the
fourth-power identity, both matrix-space dimensions, and the complete
347-variable linearized no-new-pole system. Its kernel is exactly
the 69 energy-block directions, an independent entrywise control of
the weighting convention in (3).

For every strictly filtered elementary transvection
`U=I+c E_rc`, with `H_r>H_c` and `c in F4*`, it also constructs the
full polynomial seed, checks every crossing endpoint, and verifies
its minimum weighted infinity entry order is exactly
`-8-(H_r-H_c)`. There are 834 controls, with gaps from one to ten.
These finite controls test the implementation; (2)--(4) prove the
all-field and whole-group statement.

```sh
make verify-joint-filtered-gauges
python3 tools/check_joint_filtered_gauges.py --generate
```

The [evidence](../evidence/joint-filtered-gauges.json) records the
dimensions, exact gap counts, unit identity and source hashes.
Only Python's standard library is required. This work used substantial
OpenAI assistant support and retains the Gaussian reference
attribution. No new Lean, independent review, coefficient-rank
lower bound or improved factoring exponent is claimed.
