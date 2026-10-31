# Exact operator profiles over a fixed sixteen-element field

The [occupation-transfer construction](coupled-identity-operator-filters.md)
and its [Weyl reduction](identity-operator-weyl-reduction.md) extend from
\(\mathbf F_4\) to the fixed field

\[
 \mathbf F_{16}=\mathbf F_4[u]/(u^2+u+\omega),
 \qquad \omega^2+\omega+1=0.
 \tag{1}
\]

The larger field is useful for homogeneous polynomial synthesis with
eight endpoint coordinates. If each endpoint coordinate is a nonzero
linear functional on the solution space, their eight zero hyperplanes
cannot cover a vector space over a field with sixteen elements.
Consequently an endpoint with all coordinates nonzero exists over
this fixed extension. This gives a sufficient existence criterion;
over \(\mathbf F_4\), the full endpoint image must instead be checked.

No new factoring bound follows from this field extension by itself.

## Restriction of scalars preserves the exact profile

The trace of \(\omega\) from \(\mathbf F_4\) to \(\mathbf F_2\)
is one, so the quadratic in (1) is irreducible. In the basis \((1,u)\),
multiplication by \(a+ub\) has matrix

\[
 R(a+ub)=
 \begin{pmatrix}a&\omega b\\b&a+b\end{pmatrix}.
 \tag{2}
\]

Apply \(R\) coefficientwise to a complete Laurent polynomial matrix
\(A(z)\). Every matrix entry becomes a two-by-two matrix over
\(\mathbf F_4[z,z^{-1}]\). If

\[
 U A V=\operatorname{diag}(z^{\lambda_1},\ldots,z^{\lambda_d})
\]

is a Smith decomposition over \(\mathbf F_{16}[[z]]\), then
\(R(U)\) and \(R(V)\) are invertible over
\(\mathbf F_4[[z]]\): restriction of scalars also takes their
regular inverses to regular inverses. The resulting diagonal matrix
has two copies of each \(z^{\lambda_i}\). Therefore the exact
expanded \(\mathbf F_4\) Smith orders are

\[
 \lambda_1,\lambda_1,\ldots,\lambda_d,\lambda_d.
 \tag{3}
\]

The implementation checks the equality of every adjacent pair before
returning the original \(\mathbf F_{16}\) profile. This expansion
is only a way to compute the profile. It does not double the physical
alphabet or the multiplicities used in the crossing cost.

The endpoint inverse and polynomial clearing use exact field arithmetic.
A constant physical \(Y\) repair is accepted only when its endpoint
matrix preserves the energy filtration and is invertible. Polynomial
sources are retained completely before entering the native
\(\mathbf F_4\) Smith routine. Its adaptive window and source-precision
checks remain in force, including for coefficients in the extension
part of the field.

## Implementation and independent evidence

The implementation is
[f16_identity_operator_profiles.py](../tools/f16_identity_operator_profiles.py).
It provides the same `zero_profiles`, `infinity_profiles`, `profile_data`,
and `polynomial_seed` interfaces as the four-element operator wrapper.

The checker
[check_f16_identity_operator_profiles.py](../tools/check_f16_identity_operator_profiles.py)
uses three independent comparisons:

- Thirty complete Laurent matrices of orders one through three agree
  with direct evaluation of every minor over \(\mathbf F_{16}\).
- Twenty genuine extension-field diagonal filters agree with the
  separate complete-minor diagonal engine. Ten further cut profiles
  agree with the existing \(\mathbf F_4\) engine.
- Three original physical kernels agree at both cuts, before and after
  endpoint repair. These twelve original cut comparisons check 1,701
  endpoint coefficients, including a three-direction kernel and
  ramification three.

A degree-600 coefficient in the extension field gives exact orders
\((0,600)\) in a two-by-two control. It forces the complete source
past the initial precision window.

The [recorded evidence](../evidence/f16-identity-operator-profiles.json)
is reproduced by:

```sh
python3 tools/check_f16_identity_operator_profiles.py --extended \
  --output /tmp/f16-operator-replay.json
```
