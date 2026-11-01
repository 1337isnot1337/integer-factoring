# A polynomial model for all rational joint left retractions

## Result and scope

For one explicit nine-letter hot crossing and a fixed eight-letter
coordinate embedding, every rational joint left retraction with polynomial
strict output has the form

```math
A(t)=A_0(t)+(t+1)G(t)S(t),                             \tag{1}
```

where G is an arbitrary 64-by-17 polynomial matrix. The fixed binary
17-by-64 matrix S is supplied with a polynomial right inverse of degree
11. There is no pole-order or pole-location cutoff in (1).

This gives a complete finite affine search space after imposing degree
or raw-jet bounds. It has not produced a better crossing pressure. The
result is not a classification of varying embeddings, arbitrary tensor
families, or all profile-preserving gauges.

## 1. The fixed ambient crossing

Use the exact restricted reflection tensor at interval 16 on hot words

```text
0, 1, 2, 3, 4, 5, 8, 6, 9.
```

They are ordered by height and then by integer mask. Remove its common
scalar monomial and deflate the even variable by t=w squared. Let
A_amb(t)=T(t)S_swap be its 81-by-81 crossing operator, so A_amb(1)=I.
All entries are explicitly constructed binary polynomials.

Retain the first eight physical labels. Let E insert their 64 paired
coordinates into the ambient 81-dimensional space, and let F insert
the 17 excluded coordinates. Put

```math
A_0=E^{\mathsf T}A_{\rm amb}E,\qquad
B=F^{\mathsf T}A_{\rm amb}E.                           \tag{2}
```

Every rational left inverse of E is exactly
L=E transposed+C F transposed. Thus every compressed operator is

```math
L A_{\rm amb}E=A_0+C B.                               \tag{3}
```

The map L itself may have poles, including at the crossing endpoint.
We require the output in (3) to be polynomial and its value at one
to be I. The construction evaluates that polynomial output directly;
it does not evaluate an undefined internal map at its pole.

## 2. Exact polynomial certificates

The [module receipt](../evidence/joint-hot-retraction-module.json)
supplies polynomial matrices S, R, U satisfying

```math
S R=I_{17},\qquad B=U S,\qquad \det U\ne0.             \tag{4}
```

R has size 64-by-17 and degree at most 11. The checker verifies all
289 entries of the first identity and all 1,088 entries of the second
by literal polynomial multiplication. Evaluation at the nonbinary
F4 element omega gives rank U(omega)=17, certifying that its determinant
is a nonzero binary polynomial. This remains nonzero after every
characteristic-two field extension.

The construction first saturates the excluded-row module at zero and
one, then reduces its leading rows. The final right inverse is the
additional certificate that proves saturation at every finite place,
not just those two places. A bounded search that failed to find such
an inverse would not establish (4).

## 3. Why every rational retraction reduces to (1)

Work over any characteristic-two extension field K. Identities (4)
show that B and S have the same row span over K(t).

Suppose D=C B is polynomial. Since D is in that rational row span,
write D=H S over K(t). Multiplying by R gives H=D R, which is polynomial.
Thus no rational coefficient is needed once the output is polynomial.

The strict endpoint requires D(1)=0. From S(1)R(1)=I we obtain
H(1)=D(1)R(1)=0. Every entry of H is therefore divisible by t+1,
and D=(t+1)G S. This proves necessity in (1).

Conversely, any polynomial G gives a polynomial correction vanishing
at one. Since U is invertible over K(t), choose C=(t+1)G U inverse.
Equation (3) then gives exactly (1). This proves sufficiency and gives
an explicit rational retraction realizing every member of the model.

No claim that U inverse is polynomial is required.

## 4. Degree bounds make the search finite

The leading row coefficients of S are linearly independent. Its row
degrees are

```text
5, 10, 9, 11, 15, 9, 14, 10, 6, 10, 9, 18, 15, 15, 14, 15, 19.
```

For polynomial row coefficients g_i, the degree of a nonzero combination
of the strict rows (t+1)S_i is exactly

```math
\max_i\bigl(\deg g_i+\deg S_i+1\bigr).                \tag{5}
```

At the maximal degree, independence of the leading row vectors prevents
cancellation. Hence a correction row of degree at most D needs only
coefficients of degree at most D-deg(S_i)-1. Entry-specific lower and
upper coefficient bounds are linear constraints on these finitely many
coefficients. This is a complete bounded-output model, not a guessed
Laurent window for C.

For the natural regularized cuts of this fixed seed, let H_i denote the
height of paired logical row i. Preserving k initial raw coefficients
of both cuts imposes, on each correction entry D_ij,

```math
\operatorname{val}D_{ij}\ge16+k-H_i,\qquad
\deg D_{ij}\le16+H_j-k.                               \tag{6}
```

These conditions concern the actual matrix entries. They are stronger
than equality of the first k Smith multiplicities.

## 5. The original variable permits simple endpoint zeros

The base tensor is even in w, but a new tensor need not be. Applying
(4) after t=w squared gives the complete original-variable form

```math
A(w)=A_0(w^2)+(w+1)H(w)S(w^2).                        \tag{7}
```

The right inverse remains polynomial. Thus (7) permits odd coefficients
and only a simple zero of the correction at the crossing endpoint.
Restricting to even outputs would unnecessarily require a double zero.

For even entry bounds obtained by doubling (6), the original-variable
correction space has twice the dimension of the deflated strict space.
Indeed, write a correction as E(t)+w O(t). Polynomial row-space
membership holds separately for E and O, and endpoint strictness says
E(1)=O(1). The odd part has upper degree one below the corresponding
even bound. Therefore R_1=E+O and R_2=(t+1)O belong to the deflated
strict space, and

```math
E(w^2)+wO(w^2)=R_1(w^2)+\frac{R_2(w^2)}{w+1}.         \tag{8}
```

Both summands are polynomial and vanish at one. Conversely (8) constructs
every original-variable correction from two deflated strict rows.
The two parts are independent, proving the dimension statement.

## 6. Checked raw-jet rigidity and remaining freedom

The [checker](../tools/check_joint_hot_retraction_module.py) computes
the full bounded row spaces in (6) for all paired heights from zero
through 16. At k=4 the dimensions range from zero to 105, depending
on the row height. At k=8 the maximum is nine, and at k=9 it is two.
At k=10 every row space is zero.

Consequently, within this fixed-embedding rational-retraction class,
a polynomial strict output agreeing with the first ten raw deflated
jets of both regularized cuts must equal A_0. Equation (8) also gives
rigidity when the first twenty original-variable raw jets agree.
This does not classify all retractions preserving only Smith profiles.

The nonzero lower-jet families remain constructive search inputs.
Initial targeted first-residual cancellations were unsuccessful; that
is not an exclusion of rank-preserving nonlinear moves or all useful
members of the family. No improved factoring exponent is claimed.

Run `python3 tools/check_joint_hot_retraction_module.py` to reconstruct
the basis and inverse, verify (4), check the leading-row condition,
and reproduce the bounded-space dimensions. No new Lean proof or
independent mathematical acceptance is claimed.

Substantial OpenAI assistant support was used for the normal form,
implementation, computations, and text.
