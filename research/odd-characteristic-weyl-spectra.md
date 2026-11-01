# Odd-characteristic Weyl spectra and a finite component certificate

## Results and scope

The weighted Weyl construction extends to characteristic `p` with
occupation weights `1,p,p^2,...`. In odd characteristic its coefficient
algebra has nontrivial spectral projectors, so the characteristic-two
additive and radical arguments do not apply unchanged.

This note derives the coefficient algebra and tests its projectors in
characteristic three. At energy 81 there are 239 labels, beyond the
general small-dimension exclusion. The original coefficient span still
fails the finite rank-sum target: every crossing decomposition contains
one piece with row rank at least 33933 and column rank at least 1617.
This holds for arbitrary recombinations and all characteristic-three
extension fields.

The result concerns the original coefficient span at the specified
energies. It does not cover the opposite symmetry families or their
mixed row products. In particular, the top rank-one pairing is no longer
a nilpotent: `R^2=dR`, and here `d=239=2 mod 3`. The mixed algebra remains
a distinct possible construction.

## 1. Exact reduction and coefficient algebra

For `h=p^R`, use occupation labels of weighted energy `h` and define

```math
P_h(z)(\alpha,y,u,\beta)=
[\alpha+u=\beta+y]\,[t^\alpha]
\prod_j(1+t_j)^{\beta_j}
 (1+(1-z)^{p^j}t_j)^{y_j}.
```

Then `P_h(1)=E`. Every label has weight-one occupation divisible
by `p`. Dividing that occupation by `p`, retaining all other
occupations, and dividing energies by `p` gives the weights

```math
1,1,p,p^2,\ldots,p^{R-1},\qquad m=h/p.
```

Lucas' identity and Frobenius show that the original polynomial
is the reduced curve with `z` replaced by `z^p`. This preserves
the label dimension and the coefficient span exactly. For
`p=3,h=81`, the reduced energy is 27 and the weights are
`[1,1,3,9,27]`.

Normalize the reduced curve by

```math
B(z)=P(1-z)=\sum_{i=0}^m z^i B_i.
```

The pair-creation and pair-annihilation proof in
[weighted-weyl-additive-row-law.md](weighted-weyl-additive-row-law.md)
is integral before reduction. With Frobenius weights `p^j`, its
general Gauss multiplication formula works over characteristic `p`:

```math
B_iB_j=
\sum_{k=0}^{\min(i,j)}
 \binom{n}{i}\binom{n}{j}\binom{H-n-1}{k}B_n,
\quad n=i+j-k,\quad H=2m+\sum_jw_j,
\tag{1}
```

where `B_n=0` for `n>m`. Thus the coefficient span is commutative
and closed under row multiplication. Its characteristic-three
specialization generally does not reduce to divided powers.

The implementation uses the basis `B_0,...,B_m`, with `B_0=I`.
The checker also verifies every basis associativity identity in
the finite coefficient-vector multiplication table.
The original coefficients are

```math
A_j=(-1)^j\sum_{i=j}^m\binom{i}{j}B_i,
\qquad P(0)=\sum_i B_i.
\tag{2}
```

## 2. Spectral projectors and the forced local component

For each tested energy the abstract finite algebra (1) splits over
`F_3` after its nilpotent part is removed. This is checked exactly.
Choose a power `q=3^s>=m+1`; the elements `B_i^q` lie in the
reduced algebra. The checker verifies `(B_i^q)^3=B_i^q` and refines
the identity with the projectors

```math
1-(B_i^q-a)^2,\qquad a\in F_3.
```

It verifies all resulting idempotents, their pairwise orthogonality,
and their sum `I`. This gives five projectors at original energy 27
and nine at original energy 81.

Let `e_0` be the projector onto the component where every positive
`B_i` has zero residue. The checker additionally verifies

```math
e_0B_i^q=0\quad(i>0),\qquad
P(0)^q=e_0,\qquad P(0)e_0=P(0).
\tag{3}
```

These are coefficient-vector equalities in the small abstract
algebra, so they hold in the actual row representation.
Consequently `P(0)` and `e_0` have the same image and the same
row rank.

For an arbitrary coefficient combination `T=sum_i a_i B_i`,
commutativity and Frobenius give

```math
(Te_0)^q=a_0^q e_0.
```

If `a_0!=0`, `T` is invertible on the image of `e_0`. Therefore

```math
\operatorname{rank}_{\rm row}T
\ge\operatorname{rank}_{\rm row}P(0).
\tag{4}
```

Since `e_0` is nonzero in the actual representation, `a_0` is an
intrinsic functional on its tensor coefficient span. It equals
one on `I`. Thus every decomposition of the crossing from this
span contains a summand with `a_0!=0`. This argument also works
over an arbitrary extension field: its coefficients need not
lie in the prime field.

## 3. Constant column entries detect the same component

Suppose the base-three digits of occupation labels `alpha,u`
have disjoint supports in every mode: at each digit position,
at least one of the two digits is zero. The diagonal crossing
entry is then the constant polynomial one:

```math
B(z)(\alpha,u,u,\alpha)
=\sum_\kappa
 \binom{\alpha}{\kappa}\binom{u}{\kappa}z^{\operatorname{wt}(\kappa)}
=1.
\tag{5}
```

Lucas' identity kills every nonzero `kappa`. The same entry of
`T=sum_i a_i B_i` is `a_0`.

Column blocks are indexed by the coordinate sum `sigma=alpha+u`.
Every distinct sum admitting such a pair gives a nonzero block
for every `T` with `a_0!=0`. The exact enumeration finds:

| Original energy | Reduced energy | Labels | Row rank of `P(0)` | Column rank of `P(0)` | Constant-one column blocks |
|---:|---:|---:|---:|---:|---:|
| 27 | 9 | 23 | 399 | 93 | 88 |
| 81 | 27 | 239 | 33933 | 1632 | 1617 |

At energy 81, equations (4)--(5) force the same summand to have
rank product at least

```math
33933\cdot1617=54869661
>239^3>239^{2\sqrt2}.
```

Its contribution alone exceeds the finite rank-sum target.
This closes the entire original coefficient span at this energy,
including spectral-projector recombinations and arbitrary numbers
of summands. No larger raw profile scan is needed.

## 4. Projector exploration and reproducibility

At energy 27 the raw, projector-only, and projector-refined
rank sums are approximately 533.49, 462.11, and 533.49,
respectively, against target 84.29. Each raw coefficient is already
supported in one of the spectral components, so the refinement
does not improve that profile. These small numbers diagnose the
algebra; they are not the proof of the energy-81 exclusion.

```sh
python3 tools/check_odd_weyl_spectrum.py --output evidence/odd-weyl-p3-h81-certificate.json
```

The driver reconstructs the abstract algebra and its local-unit
certificate, checks 21952 basis associativity triples at energy
81, enumerates the digit-disjoint witnesses, compiles the
native rank engine, and computes the two actual ranks of `P(0)`.
The native computation checks 7122902 conservation-supported
entries at energy 81. All arithmetic is exact modulo three.
Source hashes, the projector vector, and the finite inequalities
are included in the evidence. Use `--energy 27` for the smaller
independent instance.

The general projector constructor and native engine remain
available for other explicit recombinations. The proof here
does not infer asymptotic exclusion for all odd-characteristic
families from the two finite instances.

## 5. A different central projector in the mixed algebra

There is a constructive difference from the characteristic-two
rank-one radical. Let `A` now contain all coefficients of the
four curves `P,P^T,Q,Q^T`, and let `W=A v` with `v=vec(I)`.
The algebra is closed under transpose. If the standard bilinear
pairing restricted to `W` is nondegenerate, choose a basis matrix
`U` for `W`. Then

```math
e_W=U(U^{\mathsf T}U)^{-1}U^{\mathsf T}
```

is a central idempotent of `A`. Indeed, `W` and `W^perp` are
both invariant, and the displayed matrix is the projection onto
`W` along `W^perp`. It belongs to the ideal generated by
`R=v v^T`, because all outer products of orbit vectors belong
to that ideal. The ideal is the full matrix algebra on `W`.

At characteristic three, energy 27, exact complete orbit closure
gives `dim W=66` and pairing rank 66 inside the 133-dimensional
parity-zero block. Thus this construction gives an actual
central projector in the mixed algebra, outside the preceding
original-span analysis. The exact cut ranks are:

| Tensor | Row rank | Column rank |
|---|---:|---:|
| `e_W` | 66 | 225 |
| `P(0)e_W` | 48 | 220 |
| `P(0)(I-e_W)` | 351 | 209 |

The split increases the contribution of `P(0)` from about
192.63 to about 373.61. It therefore does not justify a larger
scan of this projector family. This is an effective construction
and a finite profile, not an improved crossing certificate.

```sh
python3 tools/explore_odd_weyl_orbit.py --profile 27
```

For comparison, energy three has orbit and pairing dimensions
`(1,1)`, while energy nine has `(4,3)`. The latter pairing is
degenerate, so the displayed projector formula cannot be used
there without passing to an additional quotient.
