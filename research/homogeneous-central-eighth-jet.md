# The first characteristic-two residual jet at even state energy

## Result and scope

The homogeneous characteristic-two Appendix B profile has a uniform
eighth-order contact identity at even state energy. Its first possible
nonconstant Taylor coefficient is an explicit parameter pairing times a
diagonal residue charge on the states.

Let `a,b` be regular monomials, with respective weights `j<=h` and
`j=h mod 2`. Let the input and output states have the same even weight
`E`. Use the complete source profile from the
[homogeneous profile note](homogeneous-central-profiles.md):

```math
F_{a,b;E}(Z)=\sum_s Q_s(a,b)|_{W_E}Z^s,
\qquad C=\Gamma_0(a,b)I_{W_E}.
```

**Theorem.** Over `F_2`,

```math
F_{a,b;E}(Z)-C\ \in\ (Z-1)^8\operatorname{End}(W_E)[Z,Z^{-1}],
```

and

```math
[T^8]\bigl(F_{a,b;E}(1+T)-C\bigr)
 =\Gamma_8(a,b)\,\operatorname{diag}_{y\in W_E}\chi_8(y),
\qquad
\chi_8(y)=\#\{k\in y:k\equiv1,7\pmod8\}\pmod2.       \tag{1}
```

Section 2 defines `Gamma_8` by a finite sum of determinant minors. In
particular, (1) is an explicit formula, not a rank guessed from small
matrices. It extends bilinearly to parameter combinations of the stated
weights.

The order eight is sharp. At weight `32`, the parameters
`a=(1,3,5,23)`, `b=(1,3,13,15)` and the equal input/output state
`y=u=(1,5,9,17)` have a nonzero eighth coefficient. The complete
weight-32 residual therefore does not have uniform order sixteen.

This theorem does not supply a better factoring exponent. Removing a
fixed number of Taylor coefficients does not remove a linear fraction
of the intermediate energy, and subtracting or dividing out contact
terms requires separate control of both complete channel cuts.

The source operators, characteristic-two quotient and raw complete
Pfaffian are the ones cited in the
[homogeneous selection note](homogeneous-central-selection.md). The
new unmatched-field calculation, its minor formula and the checker
below used OpenAI assistant support. No new Lean proof is asserted.

## 1. The finite normal-order calculation

Write an odd parameter part as `n=2a+1` and an odd state mode as
`k=2c+1`. Its creation and annihilation coefficients in the quotient are

```math
C_n(k)=[a\subseteq c],\qquad A_n(k)=[a\mathbin{\&}c=0],
```

where containment refers to binary digits. Thus its normally ordered
field is

```math
\Omega_n(Z)=\sum_{k\text{ odd}>0}
 \bigl(C_n(k)t_k Z^k+A_n(k)\partial_{t_k}Z^{-k}\bigr).
                                                               \tag{2}
```

All displayed coefficients are in `F_2`. The usual factor `k` in
annihilation is one, since the surviving modes are odd. The derivative
is well defined on the square-zero quotient because the derivative of
`t_k^2` is zero in characteristic two.

The contraction of an `A` parameter part `n` with a `B` parameter
part `m`, after the source regularization, uses

```math
D(n,m)=[n\mathbin{\&}m=1].                          \tag{3}
```

More explicitly, the unregularized contraction is
`D(n,m) Z^n/(1+Z)^(n+m)`. Consider a partial matching of the parameter
parts. Let `A'` and `B'` be the unmatched parts in the two slots, let
`w_A` be the sum of the parts in `A'`, and put
`W=w_A+sum(B')`. Its contribution to the complete Laurent profile is

```math
\Gamma\,Z^{-w_A}(1+Z)^W
 :\!\prod_{n\in A'}\Omega_n(Z)
       \prod_{m\in B'}\Omega_m(1)\!:\ \bigm|_{W_E}.    \tag{4}
```

Here `Gamma` is the sum of all full contractions among the remaining
parts. Formula (4) follows by canceling the source regularizing power
against the contracted denominators. It is also immediate from the
four-block Pfaffian's matching expansion. No infinite intermediate-state
sum is evaluated at `Z=1`: after fixing the external states, every
unmatched normal product in (4) is a finite Laurent polynomial.

If `W=0`, (4) is the scalar constant `C`. Every parameter part is odd,
and the total number of parts has even parity. Therefore the number of
unmatched parts is even. Since each regular parameter contains a given
odd part at most once, four unmatched parts have total weight at least
`1+1+3+3=8`. Consequently every contribution of order below eight
comes from exactly two unmatched fields.

### Quadratic zero modes

Two unmatched fields can preserve state energy only by creating and
annihilating the same mode. Their zero mode is therefore diagonal in
the regular monomial state basis. If one field comes from each slot,
its diagonal entry on `y` is

```math
N_{n,m;y}(Z)=\sum_{k\in y}
 \bigl(C_n(k)A_m(k)Z^k+A_n(k)C_m(k)Z^{-k}\bigr).      \tag{5}
```

If both fields come from the same slot, their arguments coincide.
The diagonal entry is the value of (5) at `Z=1`, independent of that
common argument.

Because `E` is even and every part is odd, `y` contains an even number
of parts. Put `c=(k-1)/2` and denote its first two binary digits by
`c_0,c_1`. For the charge in (1),

```math
[k\equiv1,7\pmod8]=1+c_0+c_1.
```

The constant one disappears when summed over `y`. Expanding (5) by
Lucas's rule gives the following complete table needed through order
eight. Each indicated leading coefficient is `chi_8(y)`; reversing
the two fields gives the same statement.

| Unmatched parts | Their total weight `W` | Cross-slot zero mode is divisible by | Leading coefficient |
| --- | ---: | --- | --- |
| `1,1` | 2 | `T^6` | `chi_8(y)` at degree 6 |
| `1,3` | 4 | `T^4` | `chi_8(y)` at degree 4 |
| `1,5` | 6 | `T^2` | `chi_8(y)` at degree 2 |
| `3,3` | 6 | Identically zero | Zero |
| `1,7` | 8 | `T^0` | `chi_8(y)` at degree 0 |
| `3,5` | 8 | `T^0` | `chi_8(y)` at degree 0 |

Here `T=Z-1`. For example, for parts `1,1`, equation (5) is
`sum_(k in y)(Z^k+Z^(-k))`. Its coefficients through degree five are
either zero or the parity of the number of parts. At degree six the
coefficient of a single mode is `1+c_0+c_1`. The other rows follow
from the same binary-digit calculation using
`C_3(k)=c_0`, `A_3(k)=1+c_0`,
`C_5(k)=c_1`, and `A_5(k)=1+c_1`.

For same-slot fields, the pairs `1,3` and `1,5` give the parity of
the number of state parts and hence vanish. The pairs `1,7` and
`3,5` give `chi_8(y)`. Repeated same-slot parts cannot occur because
the parameters are regular monomials.

Multiplying the table by the factor `T^W Z^(-w_A)` in (4) shows that
every quadratic contribution has order at least eight. The Laurent
unit `Z^(-w_A)` leaves its first nonzero coefficient unchanged.

The only four-field configuration of total weight eight is
`1_A,3_A,1_B,3_B`. At coincident arguments its normal product contains
the square of `Omega_1`. That normally ordered square is zero: the
creation and annihilation squares vanish in the quotient, and every
mixed term occurs twice. Its contribution at order eight is therefore
zero. Contributions with larger unmatched weight cannot affect (1).
This proves the eighth-order divisibility and reduces its coefficient
to the explicitly listed unmatched pairs.

## 2. The parameter pairing is a determinant-minor sum

For two lists of odd parts `a,b`, define

```math
\Gamma_0(a,b)=
\begin{cases}
\det(D(n,m))_{n\in a,m\in b},& |a|=|b|,\\
0,&|a|\ne|b|.
\end{cases}                                        \tag{6}
```

The empty determinant is one. In characteristic two the determinant
and permanent agree, so (6) is the sum of full contraction matchings.
It also equals the source vacuum pairing.

Form the tagged vertex list consisting of the parts of `a` in slot
`A` and the parts of `b` in slot `B`. For each pair of tagged
vertices, delete them and evaluate (6) on the two remaining lists.
Include that determinant precisely in the following cases:

- The vertices occupy different slots and their unordered part pair
  is one of `1,1`; `1,3`; `1,5`; `1,7`; `3,5`.
- They occupy the same slot and their part pair is `1,7` or `3,5`.

The sum of these determinants in `F_2` is `Gamma_8(a,b)`. This is
exactly the coefficient obtained from the table in Section 1, so it
proves formula (1). The same-slot terms matter when the two original
parameters have different numbers of parts. Omitting them is not a
valid formula for arbitrary homogeneous parameters.

For parameter and state domains of size `d`, the complete tensor
`Gamma_8(x,v)[u=y]chi_8(y)` has the same rank in the two cuts:

```math
\operatorname{rank}T_8
 =\operatorname{rank}(\Gamma_8)\,
  \#\{y:\chi_8(y)=1\}.                             \tag{7}
```

This is a rank statement about the whole eighth coefficient: after
permuting rows and columns, each cut is a tensor product of the
parameter matrix and the diagonal state-charge matrix.

## 3. Sharpness, restrictions and the remaining question

For the complete regular domains at state weight `32`, there are
`23` labels. Exactly one state has charge one:

```text
(1,5,9,17).
```

The matrix of `Gamma_8` on the raw weight-32 parameters has rank
`10`. Equation (7) therefore explains the complete eighth-jet ranks
`(10,10)` recorded in the larger profile computation. The explicit
nonzero entry in the introduction shows that order sixteen fails
without a restriction.

Restricting the state labels to `chi_8=0` removes this coefficient
while preserving the source crossing identity on the restricted
labels. It does not require every intermediate operator to preserve
the charge subspace. One still needs enough state labels and a
nonsingular parameter pairing to normalize the resulting crossing.
Alternatively, parameter restrictions can make `Gamma_8` vanish.
Neither operation automatically removes the next Taylor coefficient.

The even-energy hypothesis is essential. Taking `a=b=(1)` and the
state `y=u=(1)` produces a nonzero residual coefficient below order
eight. Likewise, the small weights `8,16,24` have no state with the
displayed charge, but that finite fact does not imply that the charge
vanishes at all weights.

To turn contact identities into a stronger factoring bound, a further
construction would need to remove a linear fraction of the relevant
energy while retaining enough labels and complete factorizations in
both cuts. This theorem identifies the first contact term exactly;
it does not supply that higher-order construction.

## 4. Reproduction

Run:

```sh
python3 tools/check_homogeneous_central_eighth_jet.py \
  --output evidence/homogeneous-central-eighth-jet.json
```

The checker verifies the quadratic zero-mode table, the determinant
minor formula, complete low-weight polynomial entries, sampled larger
entries, both ranks of the factored eighth tensor, and the two negative
controls. A subset of entries is also compared with the original
integer normal-order recursion rather than the Pfaffian implementation.
The committed evidence records the exact counts and source hash.
