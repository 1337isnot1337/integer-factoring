# Polynomial-time queries for the complete Schur crossing

## Result and scope

The complete integral Schur construction admits polynomial-time coefficient
queries. No power-sum expansion, character table, factorial modulus, or list
of intermediate partitions is needed for one query.

1. A mode entry with parameter, input, and output weights at most `H` is
   an integer determinant of order at most `3 sqrt(H)`. Its entries are
   signed binomial coefficients. It is computable in polynomial bit time
   over the integers, or directly over any fixed finite field.
2. A complete crossing entry on four partitions of weight `h` has a
   rational-function determinant of order at most `4 sqrt(h)`. After the
   stated prefactor, it is a polynomial of degree at most `4h`.
3. In characteristic two every nonzero matrix entry in the second formula
   is a monomial in `z` and `z+1`. Interpolation over an explicitly constructed
   extension of `F_2` recovers all its coefficients in
   `O(h^(5/2))` field operations after polynomial preprocessing. Field elements
   have `O(log h)` bits. The result descends to `F_2`.

These are entry algorithms. They do not bound the ranks of exponentially
large matrices or improve the factoring exponent above `1/2`. They permit
rank searches on the full Schur lattice without the old factorial-precision
restriction. This differs from the earlier polynomial Pfaffian evaluator
on the smaller distinct-odd-partition quotient.

The integral Schur normalization is the predecessor result in
[PR #30 at `28105165`](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/integral-schur-crossing.md).
The underlying bosonic modes and crossing identity are Appendix B of
[`anthropics/formal-math` at `16a5ce4`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
The boson--fermion correspondence and Wick determinants are classical.
The contribution here is their explicit source-compatible query formula,
the folded characteristic-two kernels, complete cost accounting, and the
independent implementation.

## 1. Source conventions and the integral state basis

Write `p_m` for the source variable `t_m`, and let `s_lambda` be a Schur
function. The Hall pairing makes the Schur basis orthonormal. For homogeneous
`f` of weight `w`, use the source-compatible field convention

```math
Y(f,z)=\sum_r f[r]z^{r-w}.                                      \tag{1}
```

The source field of `p_m` is the divided derivative of the Heisenberg
current. Its creation coefficient is `binom(k-1,m-1)`, and its annihilation
coefficient is `(-1)^(m-1) binom(k+m-1,m-1)`. Thus its polynomial-state
recursion agrees exactly with the standard Heisenberg vertex algebra,
not merely up to an uncharged change of representation.

For the classical correspondence and its state-field extension, see
Edward Frenkel, [*Vertex algebras and algebraic curves*, Sections 2.3--2.4
and 3.3](https://math.berkeley.edu/~frenkel/BOOK/bourbaki.pdf).
The Frobenius-coordinate identification with Schur states is also stated
explicitly in Linjie Shi, Na Wang, and Minru Chen,
[*The Orthogonal and Symplectic Schur Functions, Vertex Operators and
Integrable Hierarchies*, equation (36)](https://doi.org/10.2991/jnmp.k.210330.001).
Only their classical charged-fermion correspondence is used here.

Let the Frobenius coordinates of a partition `lambda` be

```math
(a_1,\ldots,a_r\mid b_1,\ldots,b_r),
\qquad |\lambda|=\sum_i(a_i+b_i+1),\qquad r^2\le|\lambda|.
```

Use two species of fermions with
`{psi^+_u,psi^-_v}=[u+v=0]`, and zero anticommutators for equal species.
Here `u,v` are half integers. Positive modes kill the vacuum. With the
paired order below, the Schur state is

```math
|\lambda\rangle
=\epsilon_\lambda
 \prod_{i=1}^r\psi^+_{-a_i-1/2}\psi^-_{-b_i-1/2}|0\rangle,
\qquad \epsilon_\lambda=(-1)^{\sum_i b_i}.                    \tag{2}
```

The field of this state replaces each creator of index `-a-1/2` by
the divided derivative `D^(a) psi(z)=(1/a!) partial_z^a psi(z)`,
and normally orders the entire product. The rational notation for divided
derivatives causes no finite-field division: every coefficient of the
derivative series is an integer binomial coefficient. The identities
are proved over the rationals first and then reduced integrally.

## 2. One Schur-mode entry is a small determinant

To compute

```math
M(\lambda;\beta,\gamma)
=[s_\beta]\,s_\lambda[|\beta|-|\gamma|]s_\gamma,               \tag{3}
```

form a list of fermions in three consecutive blocks:

- the adjoint of the creator list of `beta`, in reverse order;
- the derivative-field list of `lambda`, in its paired order;
- the creator list of `gamma`, in its paired order.

The bra and ket indices below are their nonnegative arm or leg indices.
For fields of opposite species, the only nonzero ordered contractions
are as follows. Same-species contractions and contractions within one
block are zero.

| Ordered blocks | Indices | Contraction at `z=1` |
| --- | --- | --- |
| Bra, parameter | bra `k`, derivative `a` | `binom(k,a)` |
| Parameter, ket | derivative `a`, ket `k` | `(-1)^a binom(k+a,a)` |
| Bra, ket | indices `k,l` | `[k=l]` |

These follow by differentiating the two fermion series. For example,
the ket contraction uses `binom(-k-1,a)=(-1)^a binom(k+a,a)`.
Normal ordering suppresses contractions between two fields in the parameter
block; retaining such contractions would compute a different operator.

Let `K` be the skew matrix of ordered contractions, with the entries
above used when the first field precedes the second. Wick's formula gives
the coefficient as
`epsilon_lambda epsilon_beta epsilon_gamma Pf(K)`.
It follows directly by moving the first fermion to the right using the
canonical anticommutation relation, then inducting on the list length.
The forbidden within-block pairs are exactly the terms removed by normal
ordering.

There are equally many fields of each species. Put all positive-species
fields before all negative-species fields, retaining each internal order.
Let `nu` be the parity of that permutation, and let `B` be the resulting
positive-to-negative contraction matrix. If its order is `m`, then

```math
M(\lambda;\beta,\gamma)
=\epsilon_\lambda\epsilon_\beta\epsilon_\gamma
 (-1)^{\nu+m(m-1)/2}\det B,
\qquad m=r_\lambda+r_\beta+r_\gamma.                         \tag{4}
```

All field powers of `z` have the same total exponent in this matrix
coefficient, namely `|beta|-|gamma|-|lambda|`. Therefore evaluating at
`z=1` loses no coefficient information.

For weights at most `H`, `m<=3 sqrt(H)`. Every binomial coefficient
has `O(H)` bits. Hadamard's inequality gives `O(H^(3/2)+sqrt(H) log(H+2))`
bits for the determinant and its minors. Integer fraction-free elimination
therefore has polynomial bit cost. Over a fixed finite field, elimination
uses `O(H^(3/2))` field operations, with the binomial coefficients generated
by Pascal recursion or the digit identity in that characteristic.

In characteristic two the entries reduce to bit tests:

```math
\binom{k}{a}\equiv[a\ \text{is a bit subset of }k],\qquad
\binom{k+a}{a}\equiv[a\mathbin{\&}k=0].                      \tag{5}
```

The implementation uses a binary determinant in this case. There is no
large coefficient modulus whose size grows with `h!`.

## 3. Fold the integral dual and all angle indices

The source dual transform is

```math
T(p_m)=(-1)^{m-1}\sum_{j=1}^m\binom mjp_j.
```

Write `omega(p_m)=(-1)^(m-1)p_m`, so `omega(s_lambda)=s_(lambda')`,
where the prime denotes conjugating the Young diagram. Let

```math
\mathcal D=\sum_{m\ge2}m p_{m-1}\frac{\partial}{\partial p_m}.
```

Then `T=exp(D) omega`. In the fermion model `D` is the positive Virasoro
mode `L_1`; its commutator on a creator of index `-a-1/2` has coefficient
`a` and decreases `a` by one. Consequently

```math
e^{\mathcal D}\psi^\pm_{-a-1/2}e^{-\mathcal D}
=\sum_{j=0}^a\binom aj\psi^\pm_{-j-1/2}.                    \tag{6}
```

This may also be checked from the fermionic current and its quadratic
Virasoro mode. It matches the displayed derivation on every `p_m` and
kills the vacuum. Equation (6) is integral and finite.

Fix four partitions `x,y,u,v` of weight `h`, and let `A_s(x,y,u,v)` be
the complete Schur pieces with parameters `T(s_x),s_v`, input `s_y`,
and output coefficient `[s_u]`. Let `L_0` denote weight. Expanding the
angle factor and summing its binomial coefficients gives

```math
\sum_s A_s(x,y,u,v)z^s
=(z-1)^h
 \langle u|Y((z-1)^{L_0}T(s_x),z)Y(s_v,1)|y\rangle.          \tag{7}
```

Indeed, for a weight-`w` part `a_w` of `T(s_x)`, set
`r=s-h+i` in the defining sum. The angle sum contributes
`sum_i (-1)^i binom(h+w,i) z^(h-i) = z^h(1-z^-1)^(h+w)`.
The field in (1) contributes `z^(r-w)`; multiplying by `z^w` gives
exactly `(z-1)^(h+w)`. This proves (7), including its power of `z`.

The source support theorem gives `A_s=0` outside `[-2h,2h]`. Hence

```math
P_{x,y,u,v}(z)=z^{2h}\sum_{s=-2h}^{2h}A_s(x,y,u,v)z^s          \tag{8}
```

is an integral polynomial of degree at most `4h`. Its value at one is
`[x=v][y=u]`. The matrix formula below is evaluated away from zero and one;
the polynomial identity determines those two values too.

## 4. The complete four-block determinant

Put `q=z-1`. Take the creator lists of `u`, `x'`, `v`, and `y`.
They give four ordered blocks: bra `u`; transformed field `x'` at `z`;
field `v` at one; ket `y`. The transformed derivative of index `a` is

```math
\sum_{j=0}^a\binom aj q^j D^{(j)}\psi(z).                   \tag{9}
```

The missing half weights contribute the single common factor `q^(r_x)`.
Together with (7)--(8), the full prefactor is `z^(2h) q^(h+r_x)`.

For opposite-species fields the ordered contraction table is:

| Blocks | Indices | Contraction |
| --- | --- | --- |
| Bra, transformed `x` | `k,a` | `sum_j binom(a,j) binom(k,j) q^j z^(k-j)` |
| Bra, `v` | `k,b` | `binom(k,b)` |
| Bra, ket | `k,l` | `[k=l]` |
| Transformed `x`, `v` | `a,b` | `(-1)^a binom(b,a) q^(-b-1)` |
| Transformed `x`, ket | `a,k` | `sum_j binom(a,j)(-1)^j binom(k+j,j) q^j z^(-k-j-1)` |
| `v`, ket | `b,k` | `(-1)^b binom(k+b,b)` |

Here each sum runs over `0<=j<=a`, and a binomial coefficient is zero
when its lower nonnegative index exceeds its upper one. Within-block and
same-species contractions are again zero.

The fourth row is the useful cancellation. Contract the two derivative
fields first, using `1/(z-1)` as the fermion propagator. This gives

```math
q^{-b-1}\sum_{j=0}^a(-1)^j\binom aj\binom{j+b}{j}
=(-1)^a\binom ba q^{-b-1}.                                  \tag{10}
```

For the last equality, express `binom(j+b,j)` as the coefficient of `X^b`
in `(1+X)^(j+b)`. The sum is that coefficient in `(1+X)^b(-X)^a`,
which equals `(-1)^a binom(b,a)`.

Construct the skew contraction matrix and convert it to a determinant
as in (4). Multiply by the four state signs and by
`z^(2h) q^(h+r_x)`. This is exactly (8). The determinant order is

```math
m=r_u+r_x+r_v+r_y\le4\sqrt h.                               \tag{11}
```

This proves the formula over the rationals. All entries belong to
`Z[z,z^-1,(z-1)^-1]`, and the result equals the integral polynomial
(8). Reduction modulo a prime is therefore valid without dividing by
a zero factorial or by a coefficient of an input table.

## 5. Every characteristic-two kernel entry is a monomial

In characteristic two, `q=z+1`. For nonnegative integers use `a & k`
for bitwise intersection, and `a & ~k` for the bits of `a` absent from
`k`. The first and fifth contraction rows become

```math
\begin{aligned}
\sum_j\binom aj\binom kjq^jz^{k-j}
 &=z^{k-(a\mathbin{\&}k)}=z^{k\mathbin{\&}\sim a},\\
\sum_j\binom aj\binom{k+j}{j}q^jz^{-k-j-1}
 &=z^{-k-1-(a\mathbin{\&}\sim k)}.                         \tag{12}
\end{aligned}
```

For the first identity, the surviving `j` are exactly the bit subsets
of `a & k`. Their generating polynomial is `(z+q)^(a&k)=1` after
the displayed power of `z` is removed. For the second identity,
`binom(k+j,j)` is odd exactly when `k & j=0`; the surviving subsets
are those of `a & ~k`. The same calculation applies.

The other four rows are respectively

```math
[b\subseteq k],\quad[k=l],\quad
[a\subseteq b]q^{-b-1},\quad[b\mathbin{\&}k=0].              \tag{13}
```

There are no derivative sums in this binary formula. It uses the full
Schur state space, including states discarded by the earlier regular
power-sum quotient.

## 6. Uniform coefficient queries and their cost

For characteristic two, construct `F_(2^m)` with more than `4h+3` elements
and `m=O(log(h+2))`. Enumerating monic polynomials of degree `m` and testing
irreducibility with polynomial gcds is deterministic and polynomial in
`h`. Choose `4h+1` distinct elements other than zero and one.

At each point, compute the powers of `z` and `z+1` used in (12)--(13),
build the matrix, and take its determinant. The inversions occur only
at nonzero field elements. Multiply by the prefactor from Section 4.
Interpolate a polynomial of degree at most `4h`.

Precomputing all Lagrange polynomials takes `O(h^2)` field operations:
form their common root polynomial, then use synthetic division and
evaluation for each root. The implementation follows this quadratic
procedure rather than separately multiplying every Lagrange product.
Interpolation of one entry also takes `O(h^2)` field operations.
There are `O(h)` determinants of order `O(sqrt(h))`, so their total
cost is `O(h^(5/2))` field operations. Ordinary polynomial arithmetic
represents a field operation in `poly(log(h+2))` bit time.

The resulting coefficients belong to `F_2`, because (8) already has
coefficients in that field. The implementation checks this descent and
does not round extension-field values into binary values.
The same determinant yields a polynomial-time interpolation algorithm
in every fixed characteristic, using the unsimplified table of Section 4.
Only the binary implementation and its stated operation count are supplied.

For the usual label energy `h=Theta(lambda^2)`, both one mode entry and
one complete piece entry therefore have polynomial cost in `lambda`.
Constructing an entire factor table still costs at least its output size.
Likewise, an arbitrary linear combination of exponentially many Schur
states is not a free query to this algorithm.

## 7. Independent checks

Run:

```sh
python3 tools/check_schur_fermion.py --output evidence/schur-fermion-queries.json
```

The default audit compares 3,960 integer mode entries with a separate
Jacobi--Trudi expansion and the original bosonic recurrence. It also
checks their 11,880 reductions modulo 2, 3, and 5. The complete crossing
test compares all 11,827 binary coefficient entries at weights one
through four against the predecessor's rational Schur basis change.
It checks 196 whole-polynomial values over the rationals as a separate
test of signs and prefactors. Two larger complete entry queries at
weights 16 and 32 have the correct endpoint and agree at an additional
field point outside the interpolation set.

The checks include a signed integer pairing and a wrong-conjugation
endpoint control. They do not compute the full weight-16 rank matrices
or establish an asymptotic rank improvement.
The recorded output is [schur-fermion-queries.json](../evidence/schur-fermion-queries.json),
and the implementation is [schur_fermion.py](../tools/schur_fermion.py).

The remaining task is a rank theorem or an improving exact crossing seed.
Faster entries remove a computation obstacle; they do not supply that
theorem. This work used substantial OpenAI assistant support and retains
the source and predecessor attribution above.
