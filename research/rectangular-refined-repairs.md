# Secondary-height repairs and scalar quartic freedom

## Result and scope

The mixed rectangular crossing admits additional polynomial output
repairs within its existing twelve-dimensional extreme core. A second
additive height permits transformations between states of equal
original height. The all-length paired-type crossing theorem and the
bounded active blocks remain valid.

With the Gaussian background row fixed, this gives `77` single-shear
directions over any characteristic-two field. Over `F4`, each direction
has three nonzero coefficients, for `231` individual variants. This
counts admissible algebraic variations. It is not a pressure bound or
a factoring-exponent improvement.

The scalar quartic factor can independently be changed to
`P(w)/(1+w)^2`, for any Laurent polynomial with `P(1)=1`. The
[bounded-defect note](rectangular-defect-factorization.md) proves that
these extensions keep the same seven- and nine-mode support bounds
and the same regular core coordinate bases.

## 1. A lexicographic additive flag

Index the retained output core by `i=v+4y`, where `v=0,...,3` and
`y=0,...,2` have physical state ordering `(low,vacuum,top,both)` and
`(low,vacuum,top)`. The original height is

```math
h(v,y)=h_X(v)+h_Y(y),\qquad
h_X=(0,0,1,1),\quad h_Y=(0,0,1).
\tag{1}
```

Choose any secondary additive height

```math
h'(v,y)=b_X(v)+b_Y(y).
\tag{2}
```

Let `V` have diagonal one, and permit a nonzero off-diagonal entry
`V_rs` only if

```math
(h(r),h'(r))<(h(s),h'(s))
\tag{3}
```

in lexicographic order. Ties between unused positions cause no
problem. Condition (3) makes `V` unitriangular after a permutation,
so it has determinant one and a finite triangular inverse.

Write `c=1+w` and `H=diag(c^h(i))`. If `C` is the mixed partial
correction, replace it by

```math
C'=VC.
\tag{4}
```

The new rational tensor is obtained from the old one by the output
map `H^-1 V H`. Its `(r,s)` entry is
`V_rs*c^(h(s)-h(r))`, a polynomial because (3) forbids an increase
in the original height. The inverse is polynomial for the same
reason. In particular, a local denominator of order at most one at
`c=0` remains of order at most one.

The ramified endpoint of the old tensor is `U=I+N`, with every
nonzero entry of `N` strictly lowering the original height. The new
endpoint is `VU`. Products of matrices satisfying (3) still satisfy
it. Thus `VU` has diagonal one and every off-diagonal entry strictly
lowers the ordered pair `(h,h')`.

For a tensor power, restrict `x,v` to one common full-letter type
and `u,y` to another. The sum of the original heights of the source
letters equals the sum for the target letters. Every surviving
endpoint term has nonnegative original-height drop, so every drop
must be zero. The same paired-type equality holds for (2). All
remaining off-diagonal terms would have a strictly positive
secondary-height drop, which also cannot sum to zero. Only diagonal
terms remain. The endpoint is the crossing at every word length.

This argument uses a second additive function only to select the
endpoint terms. No additional rational denominator or formal
parameter is introduced.

## 2. The single-shear family

For a single shear `V=I+aE_rs`, an equal-height direction with
`r!=s` always admits a secondary height. If the `v` states differ,
give the target `v` value `-1` and the source value `+1`, with the
other secondary weights zero. If only the `y` states differ, do
the same on `y`. This gives `h'(r)<h'(s)` explicitly.

The original height blocks have sizes `4,6,2`. There are

```math
4\cdot6+4\cdot2+6\cdot2=44
```

strict-height directions, and

```math
4\cdot3+6\cdot5+2\cdot1=44
```

oriented equal-height directions. Fixing row zero removes eight
strict-height directions and three equal-height directions. The
remaining counts are `36` and `41`, respectively.

Several individually admissible equal-height shears need not share
one additive secondary order. A combined correction must satisfy
(3) for one common pair `b_X,b_Y`. The constructor checks this
condition directly. Choosing one flag first, then selecting any
coefficients in its allowed positions, gives a valid dense family.

Row zero is fixed here to preserve the convenient all-low Gaussian
background. It is an additional search convention, not a requirement
of the lexicographic endpoint proof.

## 3. Constant parameter symmetries

The parameter gauge for the unmodified correction extends to this
family, provided the postcorrection is transformed with it. Define
the nominal core occupation

```math
q(v,y)=q_X(v)+q_Y(y),\qquad
q_X=(1,0,1,2),\quad q_Y=(1,0,1).
```

For a nonzero field element `k`, put `D=diag(k^q(i))`. Simultaneously
replace

```math
(\alpha,\beta,V)
\longmapsto
(k\alpha,\beta/k,DVD^{-1}).
\tag{5}
```

The resulting local tensor is the original local tensor multiplied
entrywise by

```math
k^{q(v)+q(y)-q(x)-q(u)}.
\tag{6}
```

Indeed, the base correction obeys `C'=DCD^-1`, while the Gaussian
and quartic preserve total creation and annihilation occupation.
Conjugating `V` in (5) gives the same covariance for `VC`. The
diagonal dressing by `H` commutes with `D`.

These are invertible constant diagonal physical factors. They
preserve every weighted Smith profile; their product is one on
matching paired types. The zero pattern of `V` and its admissible
flag are unchanged. Over `F4`, the previous four parameter classes
therefore still suffice when the entire corresponding shear family
is searched. The profile equality is not asserted for a fixed
untransformed `V`.

Field Frobenius also preserves profiles. It squares `alpha`,
`beta`, every entry of `V`, and every coefficient of the quartic
numerator. Thus a search using nonbinary quartic coefficients must
transform those coefficients as well.

## 4. Exact constructors and controls

The [refined repair tool](../tools/check_rectangular_refined_repairs.py)
provides:

- `single_shear(r,s,a)`, including explicit secondary-height data;
- `refined_repair(alpha,beta,V,bX,bY)`, returning the validated
  matrix `V*C`;
- `postcorrect_seed(...)`, independently applying the dressed
  polynomial output map to the original tensor coefficients;
- `conjugate_postcorrection(V,k)`, implementing (5).

Run the finite controls with:

```sh
python tools/check_rectangular_refined_repairs.py \
  --output evidence/rectangular-refined-repairs.json
```

The [recorded evidence](../evidence/rectangular-refined-repairs.json)
checks all `231` single shears at the two-mode core and three dense
secondary flags with a bulk mode. It performs `120,780` exact
two-copy paired-type polynomial and endpoint checks. It also checks
`27,648` entries of the simultaneous parameter/postcorrection gauge
and `9,216` Frobenius entries with a nonbinary quartic polynomial.
The finite checks audit the implementation; section 1 is the proof
for every tensor power.

The base [mixed repair constructor](../tools/check_rectangular_f4_repair.py)
also accepts `quartic_terms=[(coefficient,integer_exponent),...]`.
It checks `P(1)=1`, combines duplicate powers, and clears negative
Laurent powers by a recorded `global_w_power`. This global power
has value one at the crossing endpoint. It must be removed when
comparing the constructed numerator to the unscaled rational
defect: it contributes order `+global_w_power` at zero and
`-global_w_power` at infinity.

The scalar freedom changes the quartic coefficient inside the
existing finite CAR program. It does not introduce additional
spectator modes or additional physical alphabet states. Pressure
improvements from these directions still require exact profile
certificates and the uniform complexity argument required by the
factoring reduction.
