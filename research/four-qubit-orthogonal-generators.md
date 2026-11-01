# A direct odd-characteristic proof of the four-qubit generators

This records a direct proof of the invariant-generation step left open in
[bounded-spin-three-jets.md](bounded-spin-three-jets.md). It uses the
Gram model from that note and the argument shared in
[Discussion #51](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19324744).
The characteristic-zero context is Luque and Thibon's
[*The polynomial invariants of four qubits*](https://arxiv.org/abs/quant-ph/0212069).
The proof below works in odd characteristic, including five, without
specializing a characteristic-zero generation theorem.

This is an ordinary proof, not a new Lean result or an improved factoring
bound. It supplies scalar invariant generators. The trace-quotient
identification, ordinary primary multiplicities, and two cut-rank growth
estimates remain separate obligations of the bounded-spin construction.

## 1. Statement for the algebraic group

Let `k` be a field of characteristic different from two. Flatten a
four-qubit tensor as a generic four-by-four matrix `A`, with paired row
and column indices. Put

```math
\varepsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\quad
J=\varepsilon\otimes\varepsilon,\quad
Q=A^TJAJ,\quad d=\det A.
```

Thus `J^T=J`, `J^2=I`, and `det J=1`. Define `g_2,g_4,g_6` by

```math
\det(TI-Q)=T^4-g_2T^3+g_4T^2-g_6T+d^2.
```

The invariant ring for the **algebraic** action of `SL_2^4` is the free
polynomial ring

```math
\boxed{k[A]^{SL_2^4}=k[g_2,g_4,g_6,d]},             \tag{1}
```

with degrees `2,4,6,4`. Invariance under the finite group of `k`-rational
points alone is not the hypothesis in (1).

## 2. Convert the action to two orthogonal groups

Every tensor product `g tensor h`, with `g,h in SL_2`, preserves `J`
and has determinant one. Over an algebraic closure, the map

```math
SL_2\times SL_2\longrightarrow SO(J)
```

is onto. To see the dimension directly, its differential is
`(X,Y) -> X tensor I + I tensor Y`. An element in its kernel has both
`X,Y` scalar; their traces and the invertibility of two force them to
be zero. The image has dimension six, equal to that of the connected
group `SO_4`, and hence is the whole group. This uses the usual
algebraic-group image and connectedness facts in odd characteristic.

Write `G=g_1 tensor g_2`, `H=g_3 tensor g_4`. The flattening transforms
as `A -> G A H^T`. Set `B=A J`. Since `H^T J=J H^(-1)`, its action is

```math
B\longmapsto GBH^{-1},\qquad G,H\in SO(J).          \tag{2}
```

For the adjoint `B^*=J^(-1)B^T J`, one has
`B^*B=A^T J A J=Q` and `det B=det A`. Therefore the stated Gram
coefficients and determinant are global polynomial invariants.

## 3. Restrict generic orbits to a diagonal

Work first over the algebraic closure and change basis so that the
symmetric form is the identity. On the nonempty dense open set where
`B^T B` has four distinct nonzero eigenvalues, choose an orthonormal
eigenbasis. Distinct eigenvectors are orthogonal, and none is isotropic:
otherwise it would be orthogonal to the whole eigenbasis, contradicting
nondegeneracy of the form. Square roots normalize their lengths.
Changing one eigenvector's sign makes the basis matrix `H` lie in `SO_4`.

Choose square roots `x_i` of the four eigenvalues, with
`prod x_i=det B`; one sign adjustment makes this possible. Then
`G=B H diag(x_i)^(-1)` is orthogonal with determinant one, so

```math
B=G\operatorname{diag}(x_1,x_2,x_3,x_4)H^{-1}.
```

Thus every generic orbit meets the diagonal, using special orthogonal
matrices on both sides. This argument works in characteristic five;
it uses no averaging or factorial denominators.

## 4. Determine the restricted polynomials and extend back

An invariant restricted to the diagonal is fixed by all coordinate
permutations and all even sign changes. A permutation is realized by
conjugating with a signed permutation matrix of determinant one. An
even sign change is realized on one side of (2).

Because two is invertible, invariance under the even sign changes
forces each surviving monomial's four exponent parities to be equal.
Consequently the restricted polynomial has the unique form

```math
P(x_1^2,x_2^2,x_3^2,x_4^2)
 +(x_1x_2x_3x_4)R(x_1^2,x_2^2,x_3^2,x_4^2),
```

where permutation invariance makes `P,R` symmetric. The fundamental
theorem of symmetric polynomials works over every field, so both are
polynomials in the four elementary symmetric functions of `x_i^2`.
The fourth is `(prod x_i)^2`. The restriction is therefore a polynomial
in the first three and `prod x_i`, exactly the restrictions of the four
global invariants in (1).

Subtract that polynomial in the global generators from the original
invariant. It vanishes on the diagonal, hence on every generic orbit,
and hence identically by density. This proves generation.

The generators are algebraically independent. Given generic values of
their four coordinates, take roots of the corresponding quartic with
constant term `d^2`, choose square roots of its roots, and adjust one
sign so their product is `d`. This realizes a dense set of generator
values on diagonal matrices, excluding a polynomial relation.

Finally the result descends to `k`. In each homogeneous degree the
generator monomials are linearly independent and their coefficient
vectors are defined over `k`. Membership in their span after extending
the field implies membership over `k` by ordinary linear algebra.
This also gives an effective coefficient expansion in fixed degree.

## 5. Consequence for the bounded-spin pencil

Each spin-zero contraction entry is a homogeneous degree-`n` polynomial
in the sixteen local tensor entries, invariant under the algebraic
`SL_2^4` action. The integral invariance argument in the
[three-jet note](bounded-spin-three-jets.md) supplies that hypothesis in
characteristic five. Thus the entry is a weighted-homogeneous polynomial
of weight `n` in the four generators above. This statement also applies
to a selected subset of the spin-zero path labels.

For `K=aI+bS`, put `s=a+b`, `t=a-b`. The Gram matrix is `K^2`, with
eigenvalues `s^2,s^2,s^2,t^2`. Its generators restrict to

```math
g_2=3s^2+t^2,\quad g_4=3s^4+3s^2t^2,\quad
g_6=s^6+3s^4t^2,\quad d=s^3t.                     \tag{3}
```

The source crossing uses `s=1`, `t=1-2z`. Over `F_5`, writing the four
quantities in (3) as `alpha,beta,gamma,delta`, direct substitution gives

```math
\alpha\gamma-2\beta^2+\delta^2=0,\qquad
\gamma^2+(\alpha^2+\beta)\delta^2+2\beta^3=0.         \tag{4}
```

These are exact scalar identities for the pencil. The generation theorem
does not assert that every allowed invariant occurs as an entry of a
selected crossing, or identify the ordinary ranks of its tensor
coefficients. Those require the remaining representation and rank analysis.
