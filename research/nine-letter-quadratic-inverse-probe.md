# A complete nine-letter quadratic-inverse profile

The noncommuting [quadratic-inverse constructor](inverse-degree-crossing-compatibility.md)
has an exact nine-letter path with row constant rank one, column constant
rank one after reversal, and only one order-one infinity channel. This
passes the simple first-jet screen but does **not** improve the factoring
bound. Complete profiles reveal 78 channels at order two and give an
exact all-scale exclusion.

## Construction

Work on Mat_9(F2), with U(X)=J_9(1)X^T and N=I-U. The deterministic
kernel flag has levels 1 through 32 with dimensions
5,4 repeated fifteen times,1 repeated sixteen times. Its longest Jordan
block has multiplicity one. Thus the leading-rank entropy obstruction
for repeated minimum channels does not exclude this particular path.

Assign weight one only to input kernel level 32, and weight two to levels
2 through 31. Let S and T be the corresponding split of N. Then

```math
R(w)=(I-wS-w^2T)^{-1}(I-N),\qquad R(1)=I,
\qquad K(z)=z^{61}\Phi(R(z^{-1})).
```

The inverse has degree at most two and the actual crossing degree is 61.
All entries and changes of basis are explicit finite binary matrices.
The constructor's [first-jet record](../evidence/quadratic-inverse-nine-letter-probe.json)
checks the polynomial inverse identity and endpoint, screens all 31
single-level choices, and identifies level 32 as the unique choice with
only one order-one channel.

## Complete exact profiles

The infinity row profile and zero column profile are

```math
\mu=(0,1,\underbrace{2,\ldots,2}_{78},3),\qquad
\nu=(0,31,31,31,61,62,\underbrace{63,\ldots,63}_{75}).
```

Both have 81 entries and sum nu=4,941=81*61. The shared-prefix budget
and inverse-degree inequality also pass. These necessary conditions do
not imply a small pressure.

The [complete record](../evidence/nine-letter-complete-smith-probe.json)
contains every order. The local Smith computation resolves all orders
at precision 128. Independent block-Toeplitz ranks at precisions one
through four check both returned profiles against the original tensors.

The exploratory numerical minimum is approximately rho=0.997973261 at
s=1.016203841. It is not used to prove failure at other scales.

## Exact all-scale rejection

Write B(s)=s^61(sum s^-mu)(sum s^-nu). For 1<=s<=6/5, the 78 order-two
row channels and the bound nu_i<=63 give

```math
B(s)\ge81\cdot78\,s^{-4}\ge6318(5/6)^4.
```

For s>=6/5, the zero column channel gives

```math
B(s)\ge78s^{59}\ge78(6/5)^{59}.
```

Raising these two rational bounds to the seventh power proves that both
exceed 9^20. Consequently rho>5/7>1/sqrt(2) at every scale.
This is an exact exclusion for this particular tensor; it does not
exclude arbitrary nine-letter quadratic-inverse pencils.

The useful abstract lesson is that suppressing the first coefficient
alone can move almost the entire obstruction into the second coefficient.
A surviving pencil needs a coordinated early-order profile, with its
endpoint and inverse degree retained.

## Independent local Smith engine and reproduction

The [native engine](../tools/polynomial_local_smith.cpp) works over
F2[[t]] modulo t^P, with arbitrary-size bit polynomials. It chooses a
minimum-valuation pivot, inverts its unit part by the characteristic-two
Newton identity, and clears the remaining column using regular row
operations. Divisibility holds because the pivot has minimum valuation.
The pivot row could then be cleared by regular column operations without
changing the residual block. Thus the recorded pivot valuations are
exact Smith orders below P. An unresolved zero block causes a precision
increase, never an invented high order.

The [wrapper](../tools/check_polynomial_local_smith.py) independently
checks 168 complete small square/rectangular matrices with Toeplitz ranks,
168 lower-precision truncations, and two singular controls before the
nine-letter run. These include unrestricted binary polynomial crossing
matrices of degrees two, three and four. Every pivot division is also
checked in the truncated polynomial ring.

```sh
python3 tools/check_polynomial_local_smith.py --nine-letter
```

This optional native check requires Clang and Boost headers. It uses
clang++ with -O3 -Wall -Wextra -Werror. The installed GCC 16/Boost
combination produces a Boost-header stringop-overflow warning under these
strict flags; the Clang build passes. The repository's default make verify
continues to use Python-only checks. The focused native target is
`make verify-polynomial-native`.

The construction, exact profiles, and compatibility analysis used
substantial OpenAI Codex assistance. No new Lean source or improved
factoring score is claimed.
