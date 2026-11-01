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

## 6. Complete pencil relations over F5

The two relations in (4) generate the entire kernel of the pencil map
over `F_5`. Here is a proof that also fixes the scalar Hilbert function.
Write the four generator variables as `a,b,c,d` and put

```math
F_8=ac-2b^2+d^2,\qquad
F_{12}=c^2+(a^2+b)d^2+2b^3.
```

Let `R=F_5[a,b,c,d]/(F_8,F_12)` and `L=a^2-b`.
The first relation is irreducible: viewed as a polynomial in `a`,
its two coefficients `c` and `d^2-2b^2` are relatively prime.
Modulo it, the second relation has the nonzero remainder

```math
c^2-a^3c-abc+2a^2b^2-b^3.
```

Thus the two polynomials form a regular sequence. In particular `R`
is a two-dimensional Cohen--Macaulay complete intersection with no
embedded associated primes. Setting `L=0` gives `b=a^2` and, set
theoretically, `c=a^3`, `d^2=a^4`. This quotient has dimension one,
so no two-dimensional component of `R` lies in `L=0`. Consequently `L`
is a nonzerodivisor.

In `R[L^(-1)]`, set `u=(c-ab)/(2L)`. Substitution gives the mutually
inverse presentation

```math
b=3au-u^2,\quad c=3au^2+2u^3,\quad
d^2=u^3(a+2u),\quad L=(u+a)^2.
```

The last identity has no factor two because `L=a^2-b`; only the formula
for `u` uses `2L`. The localized ring is therefore

```math
F_5[a,u,d,(u+a)^{-1}]/(d^2-u^3(a+2u)).
```

It is a domain: `u^3(a+2u)` is not a square in `F_5(a,u)`.
Since `L` is a nonzerodivisor, `R` embeds in this domain and is a domain
itself. The actual pencil image has dimension two, as is also seen from
`u=s^2`, `a=3s^2+t^2`. Hence its prime kernel has height two and is
exactly `(F_8,F_12)`. This proves completeness of the two relations.

Record half the ordinary polynomial degree by a variable `q`, so that
`a,b,c,d` have weights `1,2,3,2`. The ordinary vector-space dimensions
over `F_5` are encoded by

```math
\frac{(1-q^4)(1-q^6)}{(1-q)(1-q^2)^2(1-q^3)}
=\frac{(1+q^2)(1+q^3)}{(1-q)(1-q^2)}.
```

Its coefficient at `q^m` is `2m-2` for `m>=3`. Thus the scalar space
obtained by restricting **all** homogeneous degree-`n` invariants to
`s=1,t=1-2z` has dimension exactly `n-2` for even `n>=6`.
Setting `s=1` is injective on a fixed homogeneous degree.

For a generator of ordinary degree `e`, differentiating (3) gives
`f'(0)=-e*f(0)/2`, `f'(1)=e*f(1)/2`, and the corresponding first
identity for its reversed degree-`e` polynomial. The Leibniz rule adds
these weights. Thus every weighted-homogeneous monomial satisfies the
three [jet constraints](bounded-spin-three-jets.md). Those constraints are
independent at even lengths at least six. The universal scalar space
is therefore precisely the three-jet constraint space. A selected
crossing's actual tensor entries can span a smaller subspace; this
calculation supplies no ordinary cut-rank bound.

The kernel also has the following optional exact SymPy replay over `F_5`.
The ordinary proof above supplies the geometric and Hilbert conclusions.

```python
from sympy import symbols, groebner
u, a, b, c, d = symbols("u a b c d")
graph = [b - 3*a*u + u**2, c - 3*a*u**2 - 2*u**3,
         d**2 - a*u**3 - 2*u**4]
G = groebner(graph, u, a, b, c, d, modulus=5, order="lex")
eliminated = [p.as_expr() for p in G.polys if not p.as_expr().has(u)]
H = groebner([a*c - 2*b**2 + d**2,
              c**2 + (a**2+b)*d**2 + 2*b**3],
             a, b, c, d, modulus=5, order="lex")
assert eliminated == [p.as_expr() for p in H.polys]
```
