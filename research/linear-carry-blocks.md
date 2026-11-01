# Interpolation, addition, and carry blocks close under Fourier queries

## Result

The polynomial-time carry primitive extends to every block of the
[logarithmic-wire multiplication network](logarithmic-wire-network.md)
that contains only interpolation, addition, and carry gates.
The block may be disconnected, may contain the entire carry chain,
and may contain an arbitrary subset of the interpolation gates.

For each prescribed set of boundary frequencies, its Fourier tensor
entry is computed in polynomial bit time and polynomial space.
The proof includes the polynomial-size character-algebra representation;
it does not assume a large-field operation is free.

This combines the linear arithmetic layers with the carry constraints.
Middle gates, which perform polynomial multiplication, are outside this
closure theorem. An efficient query for this block still does not
evaluate its contraction with arbitrary middle or input tensors.

## 1. A general surjective linear-fiber identity

Let `M` be the arithmetic prime, and let the coefficient algebra
contain the principal `M`-th root `zeta` from the companion
construction. Its coefficient characteristic is different from `M`.
Let

```math
L:\mathbb F_M^r\longrightarrow\mathbb F_M^t
```

be a surjective linear map. For any function `F` into the coefficient
algebra and any `beta in F_M^r`, put

```math
S_L(\beta)=
\sum_{w\in\mathbb F_M^r}\zeta^{\beta^Tw}F(Lw).
\tag{1}
```

If `beta notin im L^T`, then `S_L(beta)=0`. Otherwise
there is a unique `gamma in F_M^t` with `L^T gamma=beta`,
and

```math
S_L(\beta)=
M^{r-t}\sum_{c\in\mathbb F_M^t}\zeta^{\gamma^Tc}F(c).
\tag{2}
```

To prove this, split the sum into fibers of `L`. Each fiber is an
affine translate of its `(r-t)`-dimensional kernel. The character
sum over that kernel vanishes unless `beta` annihilates the kernel.
This is equivalent to `beta in im L^T`. It follows directly
from the principal-root identity in each coordinate of a kernel basis,
so it remains valid when the coefficient algebra has zero divisors.
In the nonzero case the fiber character is `zeta^(gamma^T c)`
and the fiber has exactly `M^(r-t)` elements. Surjectivity makes
`L^T` injective and hence makes `gamma` unique.

Gaussian elimination over the arithmetic field decides the membership
test and computes `gamma` in polynomial time. The scalar factor
`M^(r-t)` is a unit in the coefficient field and must be retained.

## 2. The linear map onto the carry coefficients is surjective

Take a block containing only `IB`, `ADD`, and `CARRY` gates.
Reparameterize every included interpolation gate by its polynomial
coefficient vector rather than its incoming evaluation vector. The
Vandermonde transform is invertible over `F_M`, so this is a
bijection of gate settings.

The independent linear source variables are now:

- The `A` coefficient coordinates of each included `IB` gate.
- Each incoming boundary value of an included `ADD` gate whose
  interpolation predecessor is outside the block.
- Each coefficient boundary value of an included `CARRY` gate whose
  addition predecessor is outside the block.

There are polynomially many such variables; call their vector `w`.
Carry states and their endpoint boundary values are kept separate.
Every included linear gate is determined by `w`. Every exposed
message from the linear gates is a linear function of `w`. Thus all
its boundary characters combine into one computable phase
`zeta^(beta^T w)`.

Let `c` be the vector of coefficient inputs to the included carry
gates. It is a linear function `c=Lw`. This map is surjective.
For a carry at `j=qA+r`, choose the following pivot variable:

- If its `ADD[j]` predecessor is outside the block, its coefficient
  boundary value is itself a fresh pivot variable.
- If `ADD[j]` is inside, choose its first input. This is either a
  fresh boundary value, or coefficient `q` of the included `IB[r]`.

The chosen pivot appears in its carry coefficient with coefficient one.
No other carry coefficient uses that pivot: a coefficient of `IB[r]`
is sent to exactly one addition gate, and boundary edges are distinct.
The optional second input from `IB[A+r]` is a different source
coordinate, or the constant zero at `q=0`. Hence the chosen pivots
form an identity submatrix of `L`.

Unused interpolation coefficients are permitted. They are nuisance
variables in `w` and contribute to the fiber factor in (2), or force
vanishing through their boundary frequencies. They must not be silently
dropped.

## 3. Finish with the carry recurrence

For fixed `c`, define `F(c)` to be the weighted sum of legal
settings of the included carry gates, including the prescribed digit
checks, global first/last carry conditions when present, and carry
endpoint boundary characters.

The Fourier query of the whole block is precisely (1), with the
appropriate fixed boundary-transform normalization. Apply the linear
test from Section 1. If it fails, return zero. Otherwise use (2):
`gamma` is now the coefficient-frequency vector for the carry block.

Each carry segment is evaluated by the interval-prefix recurrence in
Section 6 of the companion note. For a segment of length `s` it
uses `O(sM)` coefficient-algebra operations and `O(M)` live
algebra elements. Disjoint segments multiply. The final answer is their
product times `M^(r-t)` and the fixed normalization factor.

The linear map, its elimination, the carry recurrence, the root algebra,
and every control counter have polynomial size or cost in the original
input length. This proves the stated bit-time and space bound.
The result remains valid after extending the coefficient algebra, for
example by the square-zero witness-bit tags.

## 4. What this changes and what remains to be composed

There are now two overlapping polynomial-time block families:

| Included gate types | Exact contraction method |
|---|---|
| `MID, IB, ADD` | Bilinear Gauss elimination and linear orthogonality |
| `IB, ADD, CARRY` | Surjective linear-fiber elimination and carry-prefix recurrence |

The overlap consists of the linear layers. This gives flexibility in
placing a cut or forming a block. A block containing both unrestricted
middle gates and carry constraints does not follow from either proof:
its map into the carry coefficients is generally bilinear rather than
linear. Input-digit interval weights are also outside both block classes.

No global exponent improvement follows solely from these query bounds.
They are exact constructive operations that a further contraction
algorithm may use.

## 5. Exact checks

The companion checker constructs the full interpolation/addition map
for a small two-stage arithmetic network. For every nonempty subset of
carry positions, it verifies surjectivity, compares the compiled
linear-fiber/carry query with independent enumeration of all original
interpolation inputs, and tests phases both inside and outside the
required row space.

The tests include disconnected carry segments, all global endpoint
configurations, nontrivial fiber scalars in odd coefficient
characteristic, and a coefficient algebra with zero divisors.
The finite arithmetic prime in these tests validates the block identity;
it is not presented as a large-input integer-convolution parameter.
