# Complete filtered unit curves for nonlinear crossing searches

## Construction and scope

A crossing with additive conservation can be built from an independently
chosen filtered polynomial curve at every logical state. This note gives
all such curves with specified endpoints, explicit free coordinates, and
exact formulas for both coefficient cuts. The family strictly contains
all nonlinear products of the earlier group-direction selectors, already
over the binary field at width three.

The construction is effective in singly exponential time in the logical
width. Its endpoints, generic ranks, polynomial degree, and parameter
count are proved below. It has not supplied a crossing below the
half-score threshold. The complete factoring score remains `1/2`.

Files: [constructor](../tools/filtered_unit_crossings.py),
[exact checker](../tools/check_filtered_unit_crossings.py), and
[receipt](../evidence/filtered-unit-crossings.json).

## 1. All filtered paths with two prescribed endpoints

Let `F` be any field of characteristic two, let `lambda>=1`, and put
`D=2^lambda`. Write the Boolean group algebra as

```math
 A=F[N_1,\ldots,N_\lambda]/(N_1^2,\ldots,N_\lambda^2),
 \qquad g_S=\prod_{i\in S}(1+N_i).
```

For each logical state `S`, consider

```math
 P_S(u)=\sum_{J\subseteq[\lambda]} f_{S,J}(u)N_J,
 \qquad N_J=\prod_{j\in J}N_j.
 \tag{1}
```

The complete constraints are

```math
 f_{S,\varnothing}=1,\qquad
 \deg f_{S,J}\le |J|,\qquad
 f_{S,J}(0)=0,\qquad
 f_{S,J}(1)=[J\subseteq S]\quad(J\ne\varnothing).
 \tag{2}
```

Thus `P_S(0)=1` and `P_S(1)=g_S`. For each nonempty `J`, every
allowed polynomial has the unique form

```math
 f_{S,J}(u)=[J\subseteq S]u^{|J|}
       +\sum_{j=1}^{|J|-1}a_{S,J,j}(u^j+u^{|J|}).
 \tag{3}
```

The scalars `a_(S,J,j)` are independent. In particular, the singleton
coefficients are forced to be `u[S contains i]`. There are

```math
 q_\lambda=\sum_{k=2}^{\lambda}(k-1)\binom\lambda k
          =(\lambda-2)2^{\lambda-1}+1
 \tag{4}
```

free scalars per state, and `D q_lambda` for the complete family.
For a finite field with `q` elements this gives exactly
`q^(D q_lambda)` paths. These are counts of the specified family,
not of all possible crossing tensors.

## 2. The polynomial crossing and its complete endpoints

Define the group-algebra polynomial

```math
 Q_S(z)=z^\lambda g_S P_S((1+z)/z).
 \tag{5}
```

Every monomial of (1) has `u` degree at most `lambda`, so (5) has
no negative powers and has degree at most `lambda`. Its two exact
endpoint identities are

```math
 Q_S(1)=g_S,\qquad [z^\lambda]Q_S(z)=g_S P_S(1)=1.
 \tag{6}
```

For four labels in `F_2^lambda`, set

```math
 T(z;x,y,u,v)=[x+y=u+v]\,[g_{x+u}]Q_{x+y}(z).
 \tag{7}
```

Equation (6) gives

```math
 T(1;x,y,u,v)=[x=v][y=u],\qquad
 [z^\lambda]T(z;x,y,u,v)=[x=u][y=v].
 \tag{8}
```

The crossing endpoint in (8) has rank `D^2` in each of the two
flattenings `(x,y)|(u,v)` and `(x,u)|(y,v)`. Consequently both
polynomial cut matrices have full generic rank. A finite series
profile computation may use that proved rank; a missing pivot must
still be treated as insufficient precision, rather than as a zero
Smith factor.

There is also an exact identity throughout the family. Augmentation
of `P_S` is one, so augmentation of `Q_S` is `z^lambda`. Squaring
any group-algebra element in characteristic two leaves only the
square of its augmentation, because every `g_S` squares to one.
Therefore

```math
 Q_S(z)^2=z^{2\lambda}.
 \tag{9}
```

In the monomial basis `N_J`, multiplication by `Q_S` is triangular
with diagonal `z^lambda`. The determinant of the entire row cut is
exactly `z^(lambda D^2)`. This identity is independent of every free
coefficient in (3). It does not determine the other cut's Smith list.

## 3. Nonlinear selectors are a proper subfamily

For each nonzero group element `p`, put `N_p=1+g_p`. Select an
arbitrary subset of these directions separately at every state,
subject only to their xor being `S`. Its path is

```math
 P_S(u)=\prod_{p:\theta_p(S)=1}(1+uN_p),\qquad
 \sum_p\theta_p(S)p=S.
 \tag{10}
```

No affine dependence on `S` is required. Each `N_p` has positive
nilpotent degree. A product contributing to `u^j N_J` must therefore
have `j<=|J|`; a repeated generator makes the term zero. At zero the
product is one, while at one it is the product of the selected `g_p`,
namely `g_S`. Thus (10) satisfies all of (2). The earlier
[code-section normal forms](code-section-normal-forms.md) are contained
in this larger construction after the corresponding selector choice.

The containment is strict. At width three and state zero take

```math
 P_0(u)=1+(u^2+u^3)N_1N_2N_3.
 \tag{11}
```

Its coefficient of `u` is zero. In a product (10), that coefficient
would be `sum_p theta_p N_p`. The `D-1` elements `N_p`, for nonzero
`p`, are linearly independent: the coefficient of `g_p` singles out
its own scalar. Hence every selector would have to be zero, making
the entire product one, contrary to (11). This proves strictness
without a numerical search or a parameter-count heuristic.

For (11), the group coefficients of `Q_0` are

```math
 Q_0(z)=z^3+(1+z^2)\sum_{w\in F_2^3}g_w.
 \tag{12}
```

The exact checker also enumerates the 16 binary selector products
with xor zero at this width and verifies the same distinction.

## 4. Equivalent nilpotent loop coordinates

There is a factorized description useful for searches. Start with

```math
 P_S^{\rm base}(u)=\prod_{i\in S}(1+uN_i).
```

Every allowed path can be written uniquely, after fixing an order
on subsets by increasing size, as

```math
 P_S(u)=P_S^{\rm base}(u)
       \prod_{|J|\ge2}(1+h_{S,J}(u)N_J),
 \quad \deg h_{S,J}\le|J|,\quad h_{S,J}(0)=h_{S,J}(1)=0.
 \tag{13}
```

To prove this, multiply by the inverse of the base path; that inverse
is the base path itself. Successively remove the coefficient of each
`N_J` by its factor `1+h_J N_J`, which is also its own inverse.
This operation cannot change a lower nilpotent degree. Products on
disjoint subsets add both degrees, so the filtration bound is
preserved. The remainder equals one at both endpoints, forcing
`h_J(0)=h_J(1)=0` at every step. This proves existence and uniqueness.

The coefficient formula (3) avoids the triangular change between
these loop coordinates. Both descriptions permit genuinely nonlinear
state dependence and non-Gaussian tensor entries. No Smith profile
is inferred solely from the loop factorization.

## 5. Both coefficient cuts and effectiveness

Write `K_i(S,w)=[z^i g_w]Q_S`. Additive conservation splits the
row cut into convolution blocks indexed by `S`, and the other cut
into convolution blocks indexed by `w`. Consequently

```math
 r_i=\sum_S\operatorname{rank}\operatorname{Conv}(K_i(S,\cdot)),
 \qquad
 c_i=\sum_w\operatorname{rank}\operatorname{Conv}(K_i(\cdot,w)).
 \tag{14}
```

This is a direct block reordering of the literal matrices in (7).
The implementation computes (14) and the checker independently
constructs both dense flattenings.

The input has `O(lambda 2^(2lambda))` field elements. All group
kernels, complete cut matrices, and exact coefficient rank factors
can be constructed in `2^O(lambda)` field operations. There are at
most `lambda+1` coefficient tensors. This meets the finite
construction interface; a useful asymptotic rank estimate remains
necessary. Tensoring a fixed seed also charges every multinomial
channel through the [type-amplification transfer](multiaxis-type-amplification.md).

One retained finite control has coefficient ranks

```text
row:    (1, 0, 28, 64)
column: (64, 0, 28, 1).
```

Its complete histogram cost is `sum_i sqrt(r_i c_i)=44`. At label
size eight, that gives `rho=log2(44)/6`, which is above the half-score
crossing target. This finite profile is a reproducible control,
not an improvement of the factoring bound. Full physical-gauge
searches must compute the actual matrices and charge both matched
type partitions separately.

## 6. Exact controls

The checker covers all 256 width-two families over `F_4`, every
binary selector product through width three, random filtered curves
through width six, literal group-polynomial squares, and independent
dense coefficient ranks. It retains the proper-containment example
and the stated finite rank profile. Floating optimization does not
enter these checks.

```sh
python3 tools/check_filtered_unit_crossings.py \
  --output /tmp/filtered-unit-crossings.json
cmp /tmp/filtered-unit-crossings.json evidence/filtered-unit-crossings.json
```

The ordinary constructor and its exact controls enlarge the search
space. They do not provide an all-size contraction theorem or new
Lean formalization.
