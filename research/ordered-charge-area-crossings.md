# Growing clocks, conserved charge and ordered area

The ordered Heisenberg construction has a normal form in terms of an
integer charge walk and its area. This yields two useful extensions:
polynomially growing clocks retain a rank-one starting cut without modular
charge aliasing, and a charge-zero, area-zero binary filter needs no
extension field at all. Both have polynomially many coefficients and
complete `2^O(label width)` construction cost.

The all-length identities and cost arguments below are ordinary proofs.
Complete binary coefficient spans at lengths two and three fail the finite
half target for the area-zero example. A favorable uniform rank upper
bound remains missing; the complete factoring score is still `1/2`.

## 1. Charge and area of a physical word

Use the [three-letter Steiner cubic](steiner-cubic-crossings.md) `T(t)`
over characteristic two. Interpret physical letters as integers `0,1,2`,
and set

```math
a_i=x_i-v_i,\quad b_i=y_i-u_i,\quad c_i=a_i+b_i,
\quad C_i=\sum_{j\le i}c_j,\quad C_0=0.
```

Define the integer area

```math
R=\sum_i\left(C_{i-1}a_i-v_ic_i\right)
 =-\sum_i v_ic_i+\sum_{i<j}c_i a_j.                     \tag{1}
```

On the nonzero scalar support, either all four local letters coincide,
or both ordered pairs are distinct elements from `{0,1,2}`. Their sums
then belong to `{1,2,3}`. Hence `|c_i|<=2`, and

```math
|C_m|\le2m,\qquad |R|\le2m(m+1).                       \tag{2}
```

The recurrence in (1) computes both statistics in polynomial time. It
retains physical order.

## 2. The growing clock and its strict endpoints

Choose `e=bit_length(2m+1)` and `L=2^e-1>2m`. Construct `F_(2^e)` with
a primitive `L`th root `zeta`. Since `2^e<=4m+2`, one may enumerate all
degree-`e` binary moduli, test irreducibility, and enumerate a primitive
element by factoring `L` with trial division and testing its prime
divisors. These are polynomial-time procedures in `m`; no number-theoretic
oracle or changing characteristic is involved.

For every field element `zeta` of this order, define

```math
G_m(t)= [C_m=0],\zeta^R\prod_i T(t;x_i,y_i,u_i,v_i).    \tag{3}
```

This is exactly a normalized ordered Weyl boundary, rather than an
independent heuristic phase. On an `L`-dimensional auxiliary, let
`X e_j=e_(j+1)`, `Y e_j=zeta^j e_j`, and let `half` mean the inverse
of two **in exponents modulo the odd integer `L`**. Put

```math
\psi_j=\zeta^{-j^2/2},\qquad \phi_j=\zeta^{j^2/2},
\quad
K_i=\zeta^{-(x_i^2-v_i^2)/2}
 T(t;x_i,y_i,u_i,v_i)X^{x_i}Y^{b_i}X^{-v_i}.            \tag{4}
```

In characteristic two the odd integer `L` represents one, so
`phi^T psi=1`. For the Fourier-chirp basis
`psi_(r,j)=zeta^(-j^2/2+rj)`, direct substitution gives the local rule

```math
\zeta^{-(x_i^2-v_i^2)/2}X^{x_i}Y^{b_i}X^{-v_i}\psi_r
 =\zeta^{-x_ic_i-r a_i}\psi_{r+c_i}.                    \tag{5}
```

The basis is invertible because `zeta` is primitive. Contracting the
ordered product from right to left with `psi_0` and `phi` forces
`sum_i c_i=0 mod L`. By (2) this is exact integer equality on the scalar
support. The accumulated exponent is

```math
-\sum_i x_i c_i-\sum_{i<j}a_i c_j=R\quad(C_m=0),        \tag{6}
```

where the equality uses `(sum_i a_i)(sum_i c_i)=0`. Equations (5)--(6)
prove (3) equals `phi^T K_1...K_m psi`.

At the crossing endpoint, `u_i=y_i,v_i=x_i`, so `a_i=b_i=c_i=0`.
At the starting identity, `u_i=x_i,v_i=y_i`, so `c_i=0` at every site.
Both cases have `R=0`. Consequently

```math
G_m(1)=[x=v][y=u],\qquad G_m(0)=[x=u][y=v].              \tag{7}
```

The second cut at zero has rank one at every length. The degree in `t`
is at most `3m`, and the clock dimension is `O(m)`. All tables, entries,
coefficients and exact rank factors can be built in `2^O(m)` bit time.

The field varies with `m`, but its extension degree is only `O(log m)`.
The [uniform polynomial-degree descent](fixed-field-descent.md#4-polynomial-degree-extensions-with-uniform-presentations)
projects **whole coefficient tensors** to `F2`, with at most this factor
in either rank. The channel-width loss is `O(log log m)`; a hypothetical
uniform rank bound at parameter `rho` would retain the same leading
parameter and the existing modular factoring transfer. No such improved
rank bound is proved here.

## 3. A formal area variable and a binary return filter

Replace `zeta` in (3) by an indeterminate `z`:

```math
F_m(t,z)=[C_m=0],z^R\prod_i T(t;x_i,y_i,u_i,v_i).        \tag{8}
```

This is a Laurent polynomial with `t` degree at most `3m` and area
exponents between `-2m(m+1)` and `2m(m+1)`. Thus it has at most `O(m^3)`
coefficient positions. At `t=1` its only nonzero entry is the crossing,
with area zero. Therefore, if `F_m=sum_s z^s F_(m,s)(t)`,

```math
F_{m,0}(1)=E_m,\qquad F_{m,s}(1)=0\quad(s\ne0).          \tag{9}
```

In particular one can discard all nonzero-area parts and retain the
binary polynomial

```math
A_m(t)=[C_m=0, R=0]\prod_i T(t;x_i,y_i,u_i,v_i).          \tag{10}
```

It still has the two strict endpoints in (7), now over `F2`, and only
`3m+1` coefficient positions. Its entries need just integer charge/area
arithmetic and scalar polynomial multiplication. Complete tabulation and
both exact rank factorizations cost `2^O(m)`.

This selection is not simultaneous site symmetrization. For local tuples
`(0,1,0,2)` and `(0,2,1,0)`, the total charge is zero in both orders,
but the area is respectively two and zero. Their nonzero scalar
polynomials are `t+t^3` and `t^2+t^3`. The corresponding entry in (10)
changes from zero to `t^3+t^4+t^5+t^6` upon swapping the positions.

More generally, any uniformly constructible weight `f_m(R)` with
`f_m(0)=1` preserves (7) when used in place of the area-zero indicator.
This is a polynomial-size table of possible area values, not a claim
that finding a useful table is free. It supplies another precise family
for future rank synthesis.

## 4. Exact finite checks and limits

For (10), the exact minimum possible largest rank product in a binary
coefficient-span decomposition is **576 at length two** and **92,338 at
length three**. The lower detectors are the coefficient coordinates
`w_3` and `w_5`, respectively. Every detected mask has at least the stated
rank product, so the detector vanishes on every zero-tensor coefficient
relation and is well-defined on the actual span. It sends the crossing
to one. Upper mask lists XOR to the endpoint with no term exceeding the
stated bound, proving finite minmax sharpness.

Both values are at least `24^m`. As in the
[subgroup-return tests](subgroup-return-crossings.md), exact integer
comparisons give the finite normalized floor `rho>5/7>1/sqrt(2)`.
These are two finite failures, not an all-length obstruction or an
optimization of total rank sum. Extension-field recombinations of (10)
are outside this binary catalogue. The growing-clock tensors (3) are a
different family and are not excluded by those two tests.

The [constructor](../tools/ordered_charge_area.py),
[checker](../tools/check_ordered_charge_area.py), and
[evidence](../evidence/ordered-charge-area.json) provide:

- Direct auxiliary-vector versus integer-area comparisons on 512 supported
  words at lengths one through eight, including both cold and endpoint
  controls; the area is independently recomputed by its double sum.
- Uniform irreducible-modulus and primitive-clock construction, plus
  explicit field-descent factorizations and sharp rank controls in
  extension degrees two through eight.
- Opposite-pivot replay of 2,304 cut ranks covering every binary mask of
  the two selected area-zero coefficient spans, with lower detectors and
  positive endpoint mask lists.

```sh
make verify-ordered-charge-area
python3 tools/check_ordered_charge_area.py --generate
make verify
```

The all-length identity and complete cost do not substitute for a uniform
rank upper bound. The practical next question is whether an area weight,
physical restriction, or explicit charge-walk factorization can retain
subcritical ranks as length grows. This construction used substantial
OpenAI assistant support and the existing Steiner seed. No new Lean work
or stronger factoring theorem is claimed.
