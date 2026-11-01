# A deletion-operator weight cone for coherent physical codes

## Result

The [coherent F4 retractions](coherent-hot-f4-retractions.md) admit
non-additive physical weights in a computable linear cone. Inside this
cone, both complete weighted cut profiles equal profiles of the small
selected compound matrices. The proof uses the compressed deletion
operators, not membership in a coordinate downset.

This extends the [coordinate-code integrality criterion](weighted-hot-compound-profiles.md)
to non-coordinate retracts. It supplies a search interface, not a
pressure witness below 1/sqrt(2) or a better factoring exponent.
The argument works over any field of characteristic two; the new
implementation uses F4.

## 1. Exact representations of deletion

Let P and L be homogeneous paired physical maps with LP=I. Assume
the column span of P and the column span of L transposed are invariant
under every occupation deletion D_a. Define

```math
d_a^P=L D_a P,\qquad d_a^L=P^{\mathsf T}D_a L^{\mathsf T}.
```

Then

```math
D_aP=P d_a^P,\qquad
D_aL^{\mathsf T}=L^{\mathsf T}d_a^L.
```

To prove the first identity, write D_aP=P X using invariance, then
multiply by L. The second follows by transposing LP=I and applying
the same argument. Consequently each collection d_a is a representation
of the deletion algebra: its operators commute and square to zero.

The ambient creation operator is C_a=D_a transposed. Transposing the
second intertwining identity gives

```math
L C_a=(d_a^L)^{\mathsf T}L.
```

Thus row-side creation factors act through the transpose of the
L-deletion representation. Column-side deletion factors act through
the P representation. No orthogonality assumption L=P transposed is
made.

## 2. Gaussian factors remain polynomial units

Use the stabilized hot kernel in the deflated variable t=w squared.
Its same-row Gaussian edge joining mode a on one physical leg to
mode b on the other has coefficient

```math
f_{ab}(t)=(1+t)t^{2(a\mathbin{\&}b)}.
```

The row Gaussian factor is a product of elementary factors
I+f_ab(C_a tensor C_b). Expanding the product gives the Pfaffian
matching expansion: terms with a repeated occupied mode vanish,
and disjoint edges give the matchings. Compression replaces the
creation operators by the representations in section 1. The column
factor has the corresponding deletion representation.

Define an addition edge (j,i,a) when (d_a)_{j,i} is nonzero. Its
source label is j and its target is i; deletion goes from i to j.
For two such edges, a sufficient row-weight inequality is

```math
g_x(i)-g_x(j)+g_y(k)-g_y(l)+2s(a\mathbin{\&}b)\ge0,       \tag{1}
```

where s is the positive scale multiplying t-valuations. Use the
corresponding inequality for the two column legs. Rational weights
can be cleared in one ramification variable.

Every nonzero entry of the conjugated elementary factor is regular
by (1). Its off-diagonal operator squares to zero. In characteristic
two the factor is therefore its own inverse, including after diagonal
weight conjugation. It is a polynomial unit. The product of these
factors and its inverse are regular, so removing the two outer
Gaussian factors preserves every Smith order.

This proof does not require selecting intermediate coordinate masks.
It only uses exact invariant-subspace representations and checks the
supports of their deletion matrices. Cancellations can make the cone
conservative; they cannot invalidate its sufficiency.

## 3. The two cuts and a shared sufficient cone

At zero the paired maps are L tensor L on rows and P tensor P on
columns. At infinity they are L tensor P and L tensor P in the paired
ordering specified by the coherent constructor. Exchange the y and u
physical weights at infinity.

For this mirrored choice, a common sufficient cone uses:

- x-row addition edges from d^L;
- y-row addition edges from the union of d^L and d^P supports;
- u-column addition edges from that same union;
- v-column addition edges from d^P.

Impose (1) for every permitted row pair and every permitted column
pair. The union enforces both cuts, rather than silently substituting
the zero-cut matrix for the infinity-cut matrix. The two coherent
compound matrices generally differ.

## 4. The complete paired-type pressure

Write Z_0 and Z_infinity for the two full weighted channel sums.
Their natural scalar monomial shifts cancel. Let h(x) be the physical
height sum of 2a+1 over occupied hot modes, and use t-scale s. For
already scaled additional weights g on four legs, put

```math
c_X=s h+g_x+g_v,\qquad c_Y=s h+g_y+g_u.
```

The mirrored choice doubles the additional physical means. After
dividing numerator and denominator by two, the balanced cost is

```math
\rho(p_X,p_Y)=
\frac{\tfrac12(\log Z_0+\log Z_\infty)
      +\langle p_X,c_X\rangle+\langle p_Y,c_Y\rangle}
     {H(p_X)+H(p_Y)}.
```

Its positive interior type optimum satisfies

```math
\tfrac12(\log Z_0+\log Z_\infty)
 =\rho\left(\log\sum_x e^{-c_X(x)/\rho}
             +\log\sum_y e^{-c_Y(y)/\rho}\right).
```

At zero additional weights this reduces to the previous natural
coherent pressure. The factor of two between t and w is retained.
A missing positive root or vanishing physical entropy is not accepted
as a favorable cost.

## 5. Implementation and evidence

`tools/weighted_coherent_hot_f4.py` stores both complete polynomial
coefficient planes for every occupation block and uses the existing
F4 matrix Smith kernel. Missing pivots trigger a precision increase
or an error. The [controls](../evidence/weighted-coherent-hot-f4-controls.json)
check 64 exact intertwining vectors, eight square-zero operators,
twelve commuting pairs and 288 cone inequalities. Two weighted
coordinate cases agree with the earlier independent wrapper, including
the different t/w scale convention. Thirty-three derivative coordinates
agree with finite differences in a nonsingular chamber to less than
9e-11. A weight outside the cone is rejected by the inequality test.

The derivative comparison must use a smooth chamber; at tied pivot
weights a selected chamber derivative need not equal a central finite
difference. Floating weighted pivots and pressure minimization remain
exploratory. Any favorable candidate still needs independent rational
profile replay and an exact strict pressure certificate.

Substantial OpenAI assistant support was used for the proof,
implementation, controls, and text.
