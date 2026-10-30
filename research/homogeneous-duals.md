# Effective homogeneous duals for the partition crossing construction

## Result and scope

This note constructs a new dual family for the Appendix B crossing identity.
Every dual polynomial has exactly weight `h`. Its integer coefficients have
`O(h^2 log(h+1))` bits. A uniform algorithm computes the family without
inverting a matrix of dimension `p(h)`.

The construction uses a polynomial of degree at most `h-1` in one sparse
operator. It gives an effective crossing family with the same declared
channel widths and the same parameter `rho=1/sqrt(2)` as Appendix B.
It does not prove a stronger factoring exponent.

The original polynomial ring, operators, crossing identity, and partition
estimates come from Claude's Appendix B at
[`anthropics/formal-math`, revision `16a5ce4`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The new work is the explicit homogeneous dual construction and its cost
analysis. This is an ordinary proof with finite exact checks, without new
Lean formalization.

## 1. The pairing to invert

Work in `V=Q[t_1,t_2,...]`, with `t_m` of weight `m`.
For a partition `mu`, write

```math
t^\mu=\prod_{m\ge1}t_m^{m_m(\mu)},\qquad
z_\mu=\prod_{m\ge1}m^{m_m(\mu)}m_m(\mu)!.
```

Equip the ring with the positive rational inner product

```math
\langle t^\mu,t^\nu\rangle=z_\mu[\mu=\nu].
\tag{1}
```

Let `V_h` be its weight-`h` subspace. The source's pairing on this subspace is

```math
\mathcal P_h(a,b)=[1]a[-h]b,
```

where `[1]` takes the constant coefficient. Its original duals use terms
of several weights. We construct replacements entirely in `V_h`.

Define the algebra map

```math
Y(t_m)=y_m=\sum_{j=1}^m\binom mjt_j.
```

For a partition of `h`, put `epsilon_mu=(-1)^(h-length(mu))`.
Then

```math
\mathcal P_h(t^\mu,t^\nu)
=\epsilon_\mu\langle Yt^\mu,Yt^\nu\rangle.
\tag{2}
```

**Proof.** Only differentiation terms contribute to the constant coefficient.
The two partitions must have the same number of parts. For each pairing of
a part `m` with a part `k`, the source operator contributes

```math
(-1)^{m-1}k\binom{k+m-1}{m-1}
=(-1)^{m-1}\frac{(m+k-1)!}{(m-1)!(k-1)!}.
```

On the other hand, Vandermonde's identity gives

```math
\langle y_m,y_k\rangle
=\sum_j j\binom mj\binom kj
=\frac{(m+k-1)!}{(m-1)!(k-1)!}.
```

Expanding products under (1) sums these products over all pairings of parts.
The factorials in `z_mu` supply the repeated-part multiplicities.
Multiplying the signs gives `epsilon_mu`, proving (2).

The source definitions and the positive-mode constant-coefficient rule are
in [AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
The integer recursion used in the checker is the source's `Boson.nop` in
[Boson/Defs.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/ClauseBound/Boson/Defs.lean).

## 2. A sparse operator diagonalizes the Gram form

Define two derivations on the polynomial ring:

```math
\mathcal D(t_1)=0,\quad
\mathcal D(t_m)=m t_{m-1}\ (m\ge2),\qquad
\mathcal U(t_m)=m t_{m+1}.
```

Let `H` multiply a homogeneous polynomial by its weight. Direct evaluation
on the generators, followed by the derivation rule, gives

```math
[\mathcal D,\mathcal U]=2H.
\tag{3}
```

The two derivations are adjoint under (1): `D*=U`.
For one shifted part, this follows from

```math
m\,m_m(\mu)z_\nu
=(m+1)m_{m+1}(\nu)z_\mu,
```

where `nu` replaces one part `m` by `m+1` in `mu`.
This proves the adjoint identity on monomials and hence on the ring.

The exponential `exp(D)` is finite on every polynomial and equals `Y`.
Indeed, its value on `t_m` is exactly the binomial sum defining `y_m`.
Both maps preserve products.

Define the positive Gram operator `G_h:V_h->V_h` by

```math
\langle G_h f,g\rangle=\langle Yf,Yg\rangle.
```

Projection back to weight `h` gives the finite formula

```math
G_h=\sum_{k\ge0}\frac{\mathcal U^k\mathcal D^k}{(k!)^2}
\quad\hbox{on }V_h.
\tag{4}
```

Terms with `k>=h` vanish when `h>=1`. Set `X_h=UD` on `V_h`.
Both `X_h` and `G_h` preserve the number of parts of a monomial.

### The eigenvalues

Put `Prim_j=ker(D:V_j->V_(j-1))`. For `j>=1`, `v in Prim_j`, and `t>=0`,
the commutator (3) gives

```math
\mathcal D\mathcal U^t v
=t(2j+t-1)\mathcal U^{t-1}v.
\tag{5}
```

This follows by induction on `t`. The positive inner product and adjoint
identity give the orthogonal decomposition

```math
V_h=\bigoplus_{j=1}^h\mathcal U^{h-j}Prim_j.
\tag{6}
```

For completeness, `ker D` is the orthogonal complement of `im U` at each
weight. Repeatedly apply that decomposition. The map `U` is injective at
every positive weight, since

```math
\|\mathcal Uv\|^2=\|\mathcal Dv\|^2+2j\|v\|^2>0
\quad(v\in V_j\setminus\{0\},\ j>0).
```

The vacuum is the exception: `U(1)=0`. It contributes no descendant to
positive weight, which explains the starting index `j=1` in (6).
Different summands are orthogonal by moving powers of `U` to powers of
`D` and then using the primitive condition.

On the summand indexed by `j`, equations (4)--(5) give

```math
\begin{aligned}
X_h&=\chi_j\,\mathrm{Id},
&\chi_j&=(h-j)(h+j-1)=h(h-1)-j(j-1),\\
G_h&=b_j\,\mathrm{Id},
&b_j&=\sum_{k=0}^{h-j}\binom{h-j}{k}\binom{h+j-1}{k}
=\binom{2h-1}{h-j}.
\end{aligned}
\tag{7}
```

The last equality is Vandermonde's identity. In particular, the Gram
operator has positive integer eigenvalues from a list of only `h` values.
No matrix of dimension `p(h)` needs to be diagonalized.

## 3. An integer polynomial for the inverse

The numbers `chi_1,...,chi_h` are distinct and decreasing. Define

```math
D_h=\left(\prod_{j=1}^h b_j\right)
    \left(\prod_{1\le j<k\le h}(\chi_j-\chi_k)\right)>0.
\tag{8}
```

Empty products equal one. For each `j`, the integer
`b_j*prod_(k!=j)(chi_j-chi_k)` divides `D_h`, up to sign.
Consequently the following polynomial has integer coefficients:

```math
R_h(Z)=\sum_{j=1}^h
\frac{D_h}{b_j\prod_{k\ne j}(\chi_j-\chi_k)}
\prod_{k\ne j}(Z-\chi_k).
\tag{9}
```

Its degree is at most `h-1`, and `R_h(chi_j)=D_h/b_j`.
The eigenspace decomposition therefore proves

```math
R_h(X_h)G_h=G_hR_h(X_h)=D_h\,\mathrm{Id}_{V_h}.
\tag{10}
```

This is an explicit inverse with a controlled common denominator.
It avoids a determinant or inverse of a matrix with exponentially many
rows in the label width.

## 4. The new dual polynomials

For every partition `mu` of `h>=1`, define

```math
\widetilde a^\mu
=\epsilon_\mu R_h(X_h)t^\mu,
\qquad
b^\mu=\frac{h!}{z_\mu}t^\mu.
\tag{11}
```

Both polynomials have integer coefficients and exactly weight `h`.
The second polynomial is the source's original choice.
The integrality of `h!/z_mu` is the usual cycle-type count, also proved
as the source's `ZDvd` statement.

**Homogeneous duality theorem.** For all partitions `mu,nu` of `h`,

```math
\widetilde a^\mu[-h]b^\nu
=D_h h!\,[\mu=\nu].
\tag{12}
```

The output has weight zero, so it is a scalar. To calculate it, use (2).
Every monomial in `R_h(X_h)t^mu` has the same number of parts as `mu`.
Its grading sign is therefore `epsilon_mu`. The two signs cancel, and
(10) leaves

```math
D_h\left\langle t^\mu,\frac{h!}{z_\nu}t^\nu\right\rangle
=D_hh!\,[\mu=\nu].
```

This proves (12).

## 5. Coefficient sizes and a uniform algorithm

The construction has polynomial bit lengths in `h`, despite the dimension
`p(h)` of the polynomial space.
The bounds `b_j<=2^(2h-1)` and `|chi_j-chi_k|<=h^2` give

```math
\log_2D_h\le h(2h-1)+h(h-1)\log_2 h.
\tag{13}
```

For a polynomial, let its coefficient norm be the sum of the absolute
values of its coefficients. Equation (9) gives

```math
\|R_h\|_1\le hD_h(1+h^2)^{h-1}.
```

On weight `h`, the coefficient norm of `D` grows by at most `h`.
On weight `h-1`, the coefficient norm of `U` grows by at most `h-1`.
Thus `||X_h f||_1<=h^2||f||_1`, and

```math
\|\widetilde a^\mu\|_1
\le hD_h(1+h^2)^{2h-2}.
\tag{14}
```

Every coefficient consequently has `O(h^2 log(h+1))` bits.
The common scale `D_h h!` obeys the same bound.

A direct algorithm computes (8)--(9) with integer arithmetic and then
applies `R_h(X_h)` by Horner's rule. Each application of `X_h` uses the two
sparse derivations. Intermediate polynomials lie in `V_h`, so they have
at most `p(h)` monomials. The coefficient bounds also cover these operations
up to another fixed polynomial factor in their bit lengths.
The total bit cost is polynomial in `p(h)` and `h`.

For a label width `lambda`, retain the source's smallest `h>=2` with
`p(h)>=2^lambda`. Its established estimates give

```math
h=O(\lambda^2),\qquad
\Delta_h=\sum_{E=0}^{4h}p(E)=2^{O(\lambda)}.
\tag{15}
```

Hence each homogeneous dual can be computed in `2^O(lambda)` bit time,
and its coefficients have polynomial bit length in `lambda`.
The scale `D_h h!` can be computed in polynomial time in `lambda`:
its formula has only polynomially many factors of polynomial bit length.
Finding `h` and computing partition numbers use integer dynamic programming.

The source's normal-order recursion evaluates one monomial operator
coefficient by memoization. Its states consist of a suffix of the operator
partition and two remaining partitions of energy at most `4h`.
There are at most `(h+1)*Delta_h^2` such states, with polynomially many
transitions per state. Summing over the at most `p(h)` monomials of a new
dual remains polynomial in `Delta_h` and `h`.

For a size bound, every nonzero normal-order term has at most `h` factors.
Its mode numbers have magnitude at most `4h`. Binomial and differentiation
coefficients therefore contribute at most `2^O(h^2)` per term.
Even the crude count `(8h)^h` of possible mode tuples adds only
`O(h log(h+1))` bits. Combined with (14), every factor and piece has
`O(h^2 log(h+1))` bits. The algorithm uses the memoized recursion, rather
than enumerating that crude number of tuples.

These are uniform algorithms with fixed polynomial degrees. They meet
the Appendix A size and evaluation requirements after increasing its
fixed constant `c`.

## 6. The resulting crossing family

Use the same binary-label injection `mu(x)` and piece indices
`s=-2h,...,2h` as the source. Replace its first dual by (11), and define

```math
\widetilde A_s(x,y,u,v)
=[t^{\mu(u)}]Q_s(\widetilde a^{\mu(x)},b^{\mu(v)})t^{\mu(y)}.
```

The source's `LemmaQ` applies to every polynomial of weight at most `h`.
It therefore applies to the new homogeneous polynomials. Equation (12)
gives the exact crossing sum

```math
\sum_{s=-2h}^{2h}\widetilde A_s(x,y,u,v)
=D_hh!\,[x=v][y=u].
\tag{16}
```

The two operator orders still give row and column factorizations through
the same intermediate spaces. Their channel counts remain

```math
R_{\rm row}(h,s)=\sum_{i=0}^{2h}p(s+i),\qquad
R_{\rm col}(h,s)=\sum_{i=0}^{2h}p(2h-s-i).
```

One simplification is exact:

```math
(\widetilde a^\mu)^{\langle i\rangle}
=\binom{2h}{i}\widetilde a^\mu.
```

There are `4h+1` pieces and a positive integer common scale.
All coefficients are integers, and Section 5 supplies their sizes and
uniform evaluation costs. The existing partition-width proof is unchanged.
Thus this is an effective family with `rho=1/sqrt(2)`.
The usual width trimming remains available because piece values can be
evaluated within the same uniform cost allowance.

The original operator identities are proved in
[Ring/LemmaQ.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Ring/LemmaQ.lean).
The effective-pieces interface is defined in
[AppRunningTimeStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppRunningTimeStatement.lean).

The complete pieces differ from the original ones. No claim identifies
their matrix ranks with their declared channel counts. A better rank
analysis of this homogeneous family remains open.

## 7. Exact checks and reproduction

Run:

```sh
python3 tools/check_homogeneous_duals.py
```

The checker computes the integer inverse polynomial and the sparse
derivations. It verifies the resulting dual pairings against the separate
source-based normal-order recursion in `check_partition_channels.py`.
It also computes both operator orders and sums all pieces at small weights.

The recorded checks cover all dual pairings at weights `1,...,8` and all
crossing entries at weights `1,...,3`. They include the commutator, grading,
coefficient-height bound, and two negative controls. Removing the grading
sign gives the wrong positive scale. Omitting the Gram inverse gives
different diagonal pairings instead of one common scale.

The output is [homogeneous-duals.json](../evidence/homogeneous-duals.json).
Finite checks cover their listed cases. Sections 1--6 prove the all-weight
construction and its costs.

The algebraic construction and checker used substantial OpenAI assistant
support. They retain the source's attribution and the repository's
Apache-2.0 license. The result is an alternative effective family, with
no new factoring exponent and no Lean formalization.
