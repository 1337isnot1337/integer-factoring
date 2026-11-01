# Polynomial energy gauges are formally constant near the identity

For the five-letter Gaussian reference, a polynomial gauge with several
independent matrix coefficient directions has no new effective tangent
directions at any degree. Two exact direct-sum identities prove this
without bounding the degree or assuming a square-zero generator.
They also give a nonlinear formal theorem: every identity-based formal
polynomial energy gauge retaining the reference pole and leading rank
changes the tensor only by a constant energy conjugation.

The result holds over every characteristic-two field. It classifies
formal families through the identity; it does not classify other
components of the constraint set or improve the factoring score `1/2`.

## 1. Reference and two finite identities

Use the exact reference, energy blocks and cut conventions of
[joint energy conjugations](joint-energy-conjugations.md). Put

```math
C(z)=z^8B(z)=\sum_{j=0}^{17}z^j C_j,
\qquad C_j=B_{-8+j}.
```

Let `E` be the 69-dimensional energy-block matrix algebra, `L` the
10-dimensional centralizer of `C_0` in `E`, and `F` the four-dimensional
commutant of every `C_j` in `E`. In particular

```math
\{X\in\mathcal L:[X,C_1]=0\}=\mathcal F.               \tag{1}
```

Write `J(X)=[X,C_0]` for `X in E`, and `K(Y)=[Y,C_1]` for `Y in L`.
The exact binary matrix spaces satisfy

```math
\dim J(\mathcal E)=59,\quad \dim K(\mathcal L)=6,
\quad J(\mathcal E)\cap K(\mathcal L)=0.                \tag{2}
```

There is also a stronger statement after imposing the leading rank
condition. The realignment `R(C_0)` has rank one and `(0,0)` entry one.
Write `R(C_0)_(r,c)=a_r b_c`, with `a_0=b_0=1`. Define the linear map
`Q` on matrix perturbations by the derivatives of the pivot minors:

```math
Q(M)_{rc}=R(M)_{rc}+b_cR(M)_{r0}+a_rR(M)_{0c}
                    +a_rb_cR(M)_{00}.                 \tag{3}
```

The second direct sum is

```math
\dim QJ(\mathcal E)=53,\quad \dim QK(\mathcal L)=6,
\quad QJ(\mathcal E)\cap QK(\mathcal L)=0.              \tag{4}
```

The [checker](../tools/check_joint_polynomial_tangent.py) reconstructs
both sums from the original tensor. The joint ranks are exactly `65`
in (2) and `59` in (4). Equation (1) and all these identities extend
to arbitrary characteristic-two fields by binary linear algebra.

The previously proved constant tangent identity is

```math
\ker QJ=\mathfrak d+\mathcal L,\qquad \dim\ker QJ=16,   \tag{5}
```

where `d` consists of `D_x tensor I+I tensor D_y` with diagonal physical
matrices. The dimensions are `9+10-3=16`.

## 2. All-degree tangent theorem

Let `X(t)=sum_(i=0)^D t^i X_i`, with every `X_i in E tensor K`.
Consider the first-order gauge `I+epsilon X(t)` modulo `epsilon^2`.
Suppose `[X(z^-1),C(z)]` has no negative powers, and its constant
coefficient is tangent to realignment rank one. Then

```math
X_0\in\mathfrak d+\mathcal L,
\qquad X_i\in\mathcal F\quad(i\ge1).                  \tag{6}
```

Conversely (6) satisfies those conditions, so the complete tangent
dimension for degree at most `D` is exactly `16+4D`.

For `D>=2`, the two highest pole coefficients give

```math
[X_D,C_0]=0,\qquad [X_{D-1},C_0]+[X_D,C_1]=0.          \tag{7}
```

The first puts `X_D` in `L`. By the direct sum (2), both terms in the
second equation vanish. Equation (1) gives `X_D in F`. Subtract
`t^D X_D`, which commutes with the entire `C(z)`, and repeat.
Thus only degrees zero and one need the leading-rank condition.

At degree one, the pole equation gives `X_1 in L`. The constant
coefficient is `[X_0,C_0]+[X_1,C_1]`, and its rank condition is

```math
QJ(X_0)+QK(X_1)=0.                                   \tag{8}
```

The direct sum (4) forces both terms to vanish. Since `QK` has the
same rank six as `K`, equations (1) and (5) give (6). The degree-zero
case is exactly (5). This argument covers all degrees by induction;
no extrapolation from a degree census is used.

## 3. A precise nonlinear formal statement

Let `epsilon` and `t` be independent variables. Work in

```math
\mathscr R=K[t][[\epsilon]].
```

An element of this ring has polynomial coefficients at each epsilon
order; their degrees may grow with that order. This is different from
requiring one bounded polynomial degree for the entire formal family.
Let `U(epsilon,t)` be an energy-block matrix over `R`, with
`U(0,t)=I`. It has an inverse in the same ring by the epsilon-adic
geometric series. Suppose

```math
C_U(\epsilon,z)=U(\epsilon,z^{-1})C(z)
                                  U(\epsilon,z^{-1})^{-1} \tag{9}
```

has no negative power of `z` at any epsilon order, and every two-square
minor of `R(C_U(epsilon,0))` vanishes. Its `(0,0)` entry remains a unit,
so this is precisely leading realignment rank one. Then

```math
\boxed{U(\epsilon,t)=V(\epsilon)W(\epsilon,t),}          \tag{10}
```

where `V` is constant in `t`, energy preserving, and equals `I` at
epsilon zero; every coefficient of `W` belongs to `F`, and `W` is a
unit with the same identity specialization. Consequently

```math
C_U=VCV^{-1}.                                        \tag{11}
```

The construction keeps `W(epsilon,0)=I`, so the constant factor is
in fact `V(epsilon)=U(epsilon,0)`.

The constant factor retains leading rank one. The earlier constant
formal theorem further factors it as `V=(D_x tensor D_y)L_0`, where
the physical diagonal units and the leading-centralizer unit `L_0`
all equal the identity at epsilon zero.

### 3.1 Constant boundary solutions lift at every order

We need the truncated version of the
[constant factorization proof](joint-energy-conjugations.md#6-formal-rigidity-when-the-leading-rank-may-move).
Any constant energy matrix `V_n=I mod epsilon` retaining leading rank
one modulo `epsilon^n` factors, modulo that power, into separate physical
diagonal units and a leading-centralizer unit. The proof removes one
epsilon order at a time using (5). Left multiplication by a physical
diagonal unit preserves leading rank; right multiplication by a
leading-centralizer unit does not change the conjugated `C_0`.

Lifting the factors arbitrarily to formal units gives a constant matrix
with **exact** leading rank one and the specified residue modulo
`epsilon^n`. Thus every such truncated boundary solution has a constant
formal lift. No smoothness assertion is being assumed here.

### 3.2 Remove positive polynomial coefficients inductively

Right multiplication by a unit with coefficients in `F` leaves (9)
unchanged. Suppose such corrections have already made the current
`U_n` constant in `t` modulo `epsilon^n`; this holds initially at `n=1`.
Its constant part retains leading rank one modulo that power, since
the actual transformed leading matrix does. Choose an exact constant
lift `V_n` as in section 3.1 with `V_n=U_n mod epsilon^n`. Write

```math
U_n=V_n(I+\epsilon^n X_n(t))+O(\epsilon^{n+1}).          \tag{12}
```

At epsilon order `n`, comparison with the exact boundary solution
`V_n C V_n^-1` makes both conditions homogeneous: no negative powers
in `[X_n(z^-1),C(z)]`, and (3) applied to its constant coefficient is
zero. In this coefficient comparison `V_n mod epsilon=I`, so the
derivative is exactly at the original `C_0`. There is no unaccounted
rank defect from the lower constant terms; the exact lift removed it.

The theorem of section 2 shows that
`X_n(t)=X_(n,0)+H_n(t)`, where every coefficient of `H_n` is in `F`
and `H_n(0)=0`. Right multiplication by

```math
(I+\epsilon^n H_n(t))^{-1}
```

therefore leaves the tensor unchanged and makes `U_n` constant modulo
`epsilon^(n+1)`. This correction and its inverse are well-defined in
`R` and commute with every `C_j`. Induction converges epsilon-adically.
The limiting matrix is the constant factor `V`, and the product of
the commuting-algebra corrections gives `W` in (10). The algebra `F`
need not be commutative for this argument; it is closed under products
and every member commutes with the reference.

## 4. Meaning and limits

The theorem treats all matrix coefficient directions and arbitrary
nilpotency. It goes beyond the
[scalar-polynomial square-zero classification](joint-energy-polynomial-dressings.md#6-arbitrary-rank-scalar-polynomial-generators),
but the conclusions have different scopes: that classification covers
every individual gauge in its stated class; this theorem covers formal
families through the identity. Scaling an arbitrary isolated polynomial
gauge to the identity need not preserve the pole condition, so such
a scaling is not a proof that every gauge is covered.

This formal result does not supply an effective family over an infinite
field, classify separate components, vary the prescribed energies, or
exclude transformations allowing a worse pole or larger leading rank.
The existing `F4` constant orbit catalogue does not become an all-field
profile theorem through this reduction. The best complete factoring
score remains `1/2`.

```sh
make verify-joint-polynomial-tangent
python3 tools/check_joint_polynomial_tangent.py --generate
```

The [evidence](../evidence/joint-polynomial-tangent.json) records both
direct-sum rank triples, the effective six-dimensional next-commutator
space, the full stabilizer basis and source hashes. The checker also
replays the predecessor's constant tangent identity. The all-degree
induction and nonlinear lifting are ordinary proofs above, not claims
of Lean verification or independent review. This work used substantial
OpenAI assistant support and retains the Gaussian reference attribution.
