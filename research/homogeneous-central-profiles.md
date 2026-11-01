# Complete homogeneous profiles and their first nonzero jets

## Scope

This note records exact finite calculations for the homogeneous
characteristic-two crossing. The calculations use complete pieces,
including every angle in the source identity. They test whether the
two-angle rectangular families from the
[selection note](homogeneous-central-selection.md) have useful actual
cut ranks. None of the recorded cases improves the finite seed
criterion, so none proves a factoring exponent above `1/2`.

The useful additional observation is a distinction between small and
larger energies. The complete residual has Taylor order `16` at energies
`8` and `16`, but only order `8` at energy `32`. Consequently the
small-energy behavior does not establish a uniform order-`16` identity.
The [eighth-jet note](homogeneous-central-eighth-jet.md) gives the
underlying order-`8` identity and its first obstruction.

## 1. Arbitrary state energy and a raw basis

Let `a` be a regular monomial of weight `j<=h`, let `b` have weight
`h`, and let the input and output states `y,u` have common energy `E`.
Write

```math
F_{a,b;u,y}(Z)=\sum_s[t^u]Q_s(a,b)t^y\,Z^s,
\qquad C_{a,b;u,y}=D(a,b)[u=y],
```

where `D(a,b)=[1]a[-h]b`. The source sum identity gives `F(1)=C`.
The four-block polynomial matrix in Section 5 of the selection note
works without the restriction `E=h`:

```math
P_{a,b;u,y}(Z)=\operatorname{Pf}M(Z)
             =Z^{E+j}F_{a,b;u,y}(Z).                 \tag{1}
```

Its degree is at most `2E+h+j`, and its shift support is
`-E-j<=s<=E+h`. The proof uses the same matching expansion; the
unchanged-state factor contributes `Z^E` in place of `Z^h`.

Choose independent raw parameter rows `a_1,...,a_d` and partner
columns `b_1,...,b_d` of the vacuum matrix, giving an invertible matrix
`G=(D(a_x,b_v))`. Choose any `d` distinct states of energy `E`.
The raw summed crossing is `G(x,v)[u=y]`. Applying `G^{-1}` on the
`x` leg normalizes it to `[x=v][u=y]`. This is an invertible linear
map on one tensor leg. It preserves the rank of each complete piece
in both cuts, the rank of every Taylor coefficient, and the common
Taylor order of the residual. Therefore the checker can compute in
the raw basis without expanding inverse-Gram duals.

The experiments below take `E=j` and use all selected parameter
partitions as state labels. The nonsingularity of each selected
vacuum minor is checked. No formula for arbitrary rectangular Gram
ranks is assumed.

## 2. Complete cut ranks

For a complete shift piece let `r_s,c_s` be its two cut ranks, and
put `S=sum_s sqrt(r_s c_s)`. The seed criterion used here requires
`S<d^sqrt(2)`. Every family below has `h+j` a power of two. Its only
possible nonzero angles are `0` and `h+j`.

| `j` | `h` | `E` | Labels `d` | Maximum `r_s c_s` | Rank sum `S` | Threshold `d^sqrt(2)` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 16 | 16 | 5 | 144 | 44.325 | 9.739 |
| 32 | 32 | 32 | 23 | 68,121 | 2,713.234 | 84.290 |
| 31 | 33 | 31 | 20 | 54,860 | 3,486.488 | 69.173 |
| 30 | 34 | 30 | 18 | 29,328 | 2,106.617 | 59.597 |
| 28 | 36 | 28 | 16 | 22,506 | 1,328.737 | 50.453 |
| 24 | 40 | 24 | 11 | 4,136 | 406.309 | 29.700 |

The recorded evidence contains every nonzero shift and both ranks,
as well as the exact parameter, partner, and state domains. The
table shows that reducing the number of angles to two does not by
itself produce low ranks for the complete pieces. It is a finite
diagnostic, not an asymptotic impossibility theorem for this family.

## 3. Taylor coefficients of the residual

The checker translates the polynomial residual

```math
R(Z)=P(Z)-Z^{E+j}C=Z^{E+j}(F(Z)-C)
```

by `Z=1+T`. Multiplication by `Z^{E+j}` is a unit at `Z=1`, so it
preserves the first nonzero order. At that first order its coefficient
also agrees with the corresponding coefficient of `F-C`.
Higher coefficients of `R` and `F-C` need not agree.

| `(j,h,E)` | First nonzero residual order | Ranks at order `8` | Ranks at order `16` |
| --- | ---: | --- | --- |
| `(16,16,16)` | 16 | `(0,0)` | `(4,4)` |
| `(32,32,32)` | 8 | `(10,10)` | `(160,160)` |
| `(31,33,31)` | 4 | `(200,200)` | `(278,278)` |
| `(30,34,30)` | 8 | `(77,77)` | `(178,178)` |
| `(28,36,28)` | 8 | `(105,105)` | `(149,149)` |
| `(24,40,24)` | 16 | `(0,0)` | `(46,46)` |

The evidence also gives every power-of-two order through the degree
bound. The odd-energy row is deliberately included: even state
energy is a real hypothesis in the order-`8` identity, rather than
a cosmetic restriction.

## 4. Independent checks and reproduction

Run the small validation with:

```sh
python3 tools/check_homogeneous_central_profiles.py
```

Reproduce the six larger profiles with:

```sh
python3 tools/check_homogeneous_central_profiles.py \
  --cases 16:16,32:32,31:33,30:34,28:36,24:40 \
  --output evidence/homogeneous-central-profiles.json
```

The driver compiles the associated C++17 checker. Its arithmetic is
exact over `F_2[Z]`; it uses fraction-free Pfaffian elimination and
rejects an inexact polynomial division or a capacity overflow.
Carryless CPU multiplication is optional, with a scalar implementation
available when the feature is absent. The finite checker currently
accepts at most `64` labels and uses a bounded polynomial buffer.
These implementation bounds are explicit; the theoretical Pfaffian
formula has no such fixed label limit.

The independent Python path computes a polynomial determinant and
takes its characteristic-two square root. The small validation
compares `274` complete polynomial entries between these algorithms,
and checks their translated Taylor orders and jet ranks. A further
`20` complete polynomial entries are compared with the original
integer normal-order recursion, explicitly summing all intermediate
states and all source angles. Every larger run verifies every entry
of the summed crossing, and computes exact ranks of every complete
piece. The evidence pins the C++ source by its SHA-256 hash.

The formulas depend on the original Appendix B identities and the
prior quotient/Pfaffian results cited in the selection note. These
calculations and this explanatory note used OpenAI assistant support;
they make no Lean formalization claim.
