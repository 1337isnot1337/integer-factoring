# Integral Schur crossing families and prime-power block composition

## Result and scope

The full ring of symmetric functions gives a scale-one integral form of
the Appendix B crossing construction. Its weight-`h` space has dimension
`p(h)`, over every field, including characteristic two. The construction
does not discard even modes or impose a square-zero quotient.

There are two effective versions:

1. At every positive weight, translated Schur functions give unit duality
   over the integers. The resulting crossing pieces retain the source's
   partition channel bounds.
2. At every prime-power weight `h=p^j`, the homogeneous vacuum pairing is
   invertible over `F_p`. Its inverse gives homogeneous unit duals.
   Tensoring a logarithmic number of these blocks supplies every binary
   label length. The piece count is `2^O((log lambda)^2)`, and the total
   channel width is `2sqrt(2)lambda+O((log lambda)^2)`.

All factor, coefficient, and index computations have uniform
`2^O(lambda)` bit cost for each fixed prime. The
[extended transfer theorem](subexponential-piece-transfer.md) therefore
applies at the exact parameter `rho=1/sqrt(2)`. The prime-power restriction
causes no loss in the resulting exponent `1/2`.

This is a new effective representation and a way to compose discrete
capacity blocks. It establishes no exponent larger than `1/2`.
Its useful additional feature is that the homogeneous angles at
prime-power weights have only two or three terms. A better asymptotic
rank bound at those weights alone would now suffice for a better rate.

The source operators and crossing identity are prior work from
[`anthropics/formal-math`, revision `16a5ce4`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
The symmetric-function facts below are classical; see Macdonald,
[*Schur Functions: Theme and Variations*, formulas (0.2), (0.12), (0.13)](https://www.mat.univie.ac.at/~slc/opapers/s28macdonald.pdf),
and [*Symmetric Functions and Hall Polynomials*, Chapter I](https://doi.org/10.1093/oso/9780198534891.001.0001).
Integral forms for vertex operators are established background, including
[Huang and Jing, *Lattice structure of modular vertex algebras*](https://arxiv.org/abs/2104.05223).
We prove the preservation statement for the exact source operators here.

## 1. The integral lattice and the mode convention

Write `p_k=t_k` for the source variables, of weight `k`. Define

```math
H(u)=\exp\!\left(\sum_{k\ge1}\frac{p_k u^k}{k}\right)
=\sum_{n\ge0}h_nu^n,
\qquad
\Lambda_{\mathbb Z}=\mathbb Z[h_1,h_2,\ldots]
\subset\mathbb Q[p_1,p_2,\ldots].
\tag{1}
```

The Schur functions `s_mu`, defined by the Jacobi--Trudi determinant,
form an integral basis. The weight-`E` part has rank `p(E)`.
The Hall pairing and the integral involution `omega` satisfy

```math
\langle p_\mu,p_\nu\rangle=z_\mu[\mu=\nu],\qquad
\langle s_\mu,s_\nu\rangle=[\mu=\nu],\qquad
\omega(p_k)=(-1)^{k-1}p_k.
\tag{2}
```

In particular, `omega(s_mu)=s_(mu')`, where `mu'` is the conjugate
partition. The lattice is not the polynomial ring `Z[p_1,p_2,...]`:
its elements may have rational power-sum coefficients. Reduction modulo
`p` will always mean reduction of integral Schur coordinates.

The source's mode `a[r]` changes the weight of a state by `r`.
For homogeneous `a` of weight `w`, use the generating-field convention

```math
\mathcal Y(a,z)=\sum_{r\in\mathbb Z}a[r]z^{r-w}.
\tag{3}
```

This shift is needed when comparing with a usual vertex-operator field.
All modes below are the source's actual modes, defined in
[AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).

## 2. Every source mode preserves the integral lattice

For an auxiliary variable `u`, the generating field has the exact formula

```math
\sum_{n\ge0}u^n\mathcal Y(h_n,z)
=M_{H(z+u)/H(z)}\circ T_{z,u},
\qquad
T_{z,u}(p_k)=p_k+z^{-k}-(z+u)^{-k}.
\tag{4}
```

Here `M` denotes multiplication, and inverse powers of `z+u` are expanded
in `u/z`. The multiplication acts after the substitution, as required
by the source's normal ordering.

To verify (4), sum the creation and annihilation coefficients in the
source definition. For a creation of `p_l`, the exponent in `H(u)` gives

```math
\sum_{n=1}^{l}\frac{u^n}{n}
  \binom{l-1}{n-1}z^{l-n}
=\frac{(z+u)^l-z^l}{l}.
```

For a derivative with respect to `p_k`, the corresponding sum is

```math
k\sum_{n\ge1}\frac{u^n}{n}(-1)^{n-1}
  \binom{k+n-1}{n-1}z^{-k-n}
=z^{-k}-(z+u)^{-k}.
```

Exponentiating the creations gives the multiplier in (4).
Exponentiating the commuting derivatives gives its substitution.

Both parts are integral on `Lambda_Z`. The multiplier is integral
because `H(z)` has constant term one and thus has an inverse with
integral symmetric-function coefficients. For the substitution, apply
it to the generators simultaneously:

```math
T_{z,u}H(v)
=H(v)\frac{1-v/(z+u)}{1-v/z}.
\tag{5}
```

Every coefficient in this expression is integral. No division by a
positive integer remains in (5).

For products `H(u_1)...H(u_l)`, the same normal-order computation gives
the multiplier `prod_j H(z+u_j)/H(z)` and substitution

```math
p_k\longmapsto p_k+
\sum_{j=1}^{l}\bigl(z^{-k}-(z+u_j)^{-k}\bigr).
\tag{6}
```

This is a generating identity for the field of a product of complete
symmetric functions. It is not an assertion that ordinary products of
the corresponding fields are already normally ordered. Extracting
`prod_j u_j^(mu_j)` supplies the field of `h_mu`. These products span
the lattice, so every mode `a[r]`, for `a in Lambda_Z`, maps
`Lambda_Z` to itself.

The formal series cause no infinite coefficient computation. For a fixed
operator coefficient of weight `w` and an input of weight `E`, negative
field powers are bounded below by `-E-w`. At a fixed field power, the
output weight is fixed, and the source's finite normal-order sum applies.
Equations (4)--(6) prove integrality for those finite coefficients.

Consequently all source identities between these modes hold over
`Lambda_Z`, and can be reduced modulo any prime or extended to any
commutative coefficient algebra. This is a change of integral lattice,
without a projection or a truncation of the operators.

## 3. Translated Schur functions have unit duality

Use the derivations from the
[homogeneous-dual construction](homogeneous-duals.md):

```math
\mathcal D(p_1)=0,\quad
\mathcal D(p_k)=kp_{k-1}\ (k\ge2),\qquad
\mathcal U(p_k)=kp_{k+1}.
```

Their divided powers preserve the integral lattice. Explicitly,

```math
e^{z\mathcal D}H(u)=H\!\left(\frac{u}{1-zu}\right),
\qquad
e^{z\mathcal U}H(u)=\frac{H(u+z)}{H(z)}.
\tag{7}
```

These identities follow by summing the derivations on each `p_k`.
They are coefficientwise integral and multiplicative, hence apply to
every element of `Lambda_Z`. In particular `T=e^D` preserves the
lattice and does not increase weight.

For homogeneous `f,g` of the same weight `h`, the source's duality lemma,
extended bilinearly from power sums, says

```math
[1](T\omega f)[-h]g=\langle f,g\rangle.
\tag{8}
```

Indeed, `T omega(p_mu)` is exactly the source's `aFam(mu)`, and its
pairing with `p_nu` is `z_mu[mu=nu]`. The source's factor `h!` comes
from its separate choice `bFam(nu)=(h!/z_nu)p_nu`.
Thus the integral families

```math
a_\mu=T\omega(s_\mu),\qquad b_\nu=s_\nu
\quad(\mu,\nu\vdash h)
```

satisfy

```math
a_\mu[-h]b_\nu=[\mu=\nu].
\tag{9}
```

The first polynomial has weight at most `h`, and the second has weight
exactly `h`. There is no factorial multiplier to invert in a finite field.

For any selection of `2^lambda` partitions of `h`, index these states
and duals by binary labels. Let `[s_mu(u)]` denote the Schur-coordinate
functional and set

```math
A_s(x,y,u,v)=
[s_{\mu(u)}]Q_s(a_{\mu(x)},b_{\mu(v)})s_{\mu(y)},
\qquad -2h\le s\le2h.
\tag{10}
```

The source's `LemmaQ` gives, over the integers and after every base change,

```math
\sum_{s=-2h}^{2h}A_s(x,y,u,v)=[x=v][y=u].
\tag{11}
```

Both operator orders factor (10) through the full intermediate Schur
spaces. The channel counts are still

```math
R_s=\sum_{i=0}^{2h}p(s+i),\qquad
C_s=\sum_{i=0}^{2h}p(2h-s-i),
\tag{12}
```

with `p(E)=0` for negative `E`. For example, the row factorization uses
`[s_sigma]a^(angle i)[s-h+i]s_mu(y)` and
`[s_mu(u)]b[h-s-i]s_sigma`, with the sign on either factor.
Every factor is integral by Section 2. The column factorization follows
from the other order in `LemmaQ`.

Choose the least `h>=2` with `p(h)>=2^lambda`.
Writing `c=pi sqrt(2/3)/log(2)`, the partition estimate
`log2 p(E)=c sqrt(E)+O(log(E+2))` yields

```math
h=O(\lambda^2),\qquad c\sqrt h=\lambda+O(\log(\lambda+2)),
```

```math
\begin{aligned}
\lceil\log_2R_s\rceil+\lceil\log_2C_s\rceil
&\le c\bigl(\sqrt{2h+s}+\sqrt{2h-s}\bigr)+O(\log(h+2))\\
&\le2\sqrt2\,\lambda+O(\log(\lambda+2)).
\end{aligned}
\tag{13}
```

This is the source's channel estimate in a different integral basis.
Each separate width can be capped at `2lambda`: if its declared channel
count is too large, factor that flattening through its own row or column
index space using an identity matrix. This operation does not increase
the sum in (13). Section 6 supplies the entry algorithm needed for this
replacement.

## 4. Homogeneous unit duals at every prime-power energy

The homogeneous pairing is also integral:

```math
B_h(a,b)=[1]a[-h]b,
\qquad a,b\in\Lambda_{\mathbb Z,h}.
```

Let `G_h=P_h e^U e^D` on weight `h`. Equation (7) shows that its matrix
in the Schur basis is integral. Over the rationals, the
[Gram decomposition](homogeneous-duals.md#2-a-sparse-operator-diagonalizes-the-gram-form)
gives

```math
G_h=\sum_{k\ge0}\frac{\mathcal U^k\mathcal D^k}{(k!)^2},
\qquad
B_h(a,b)=\langle\omega a,G_hb\rangle.
\tag{14}
```

Its eigenvalues are

```math
b_j=\binom{2h-1}{h-j},\qquad 1\le j\le h,
\tag{15}
```

with multiplicities `dim Prim_j`: one at `j=1`, and
`p(j)-p(j-1)` for `j>=2`. This determinant identity holds for the
integer Schur matrix. It can therefore be reduced modulo a prime
without asserting a primitive decomposition or semisimplicity there.

If `h=p^r`, Lucas's binomial congruence shows that every number in (15)
is a unit modulo `p`. For `0<=q<h`, the lower `r` base-`p` digits of
`2h-1` are all `p-1`, while the top digit of `q` is zero. Hence

```math
\binom{2h-1}{q}\equiv
\prod_{k=0}^{r-1}\binom{p-1}{q_k}\not\equiv0\pmod p.
\tag{16}
```

The case `h=1` is immediate. The Hall matrix is the identity in the
Schur basis, and `omega` is an invertible permutation. Equations
(14)--(16) prove that `B_h` is invertible on the **full** space
`Lambda_(F_p,h)`, of dimension `p(h)`.

Compute the inverse of this finite matrix over `F_p`. Its rows give
homogeneous elements `a_mu`, with `b_nu=s_nu`, satisfying (9).
Rows of the full inverse remain dual to any selected subset of Schur
labels. They may contain Schur functions outside that subset; no
restriction on their allowed support is imposed.

The resulting pieces satisfy (11)--(12), now with

```math
a_\mu^{\langle i\rangle}=\binom{2h}{i}a_\mu.
\tag{17}
```

In characteristic two, with `h=2^r`, only `i=0,2h` survive.
For odd `p`, with `h=p^r`, the surviving binomial coefficients at
`i=0,h,2h` are `1,2,1`. Including the outer sign in `Q_s`, their
weights are `1,-2,1`. These sparse angles are exact identities.
They do not by themselves imply small ranks for the complete pieces.

## 5. Prime-power blocks fill every label length

Here the additional piece allowance from
[subexponential-piece-transfer.md](subexponential-piece-transfer.md)
has a concrete use.

**Discrete-capacity composition theorem.** Fix a prime `p`. For each
`h_j=p^j`, set

```math
\ell_j=\lfloor\log_2p(h_j)\rfloor
```

and ignore zero capacities. Suppose a uniformly computable scale-one
family at each capacity `ell_j` has polynomially many pieces, separate
widths at most `2ell_j`, and width sum

```math
r+c\le4\rho\ell_j+C\log_2(\ell_j+2)
\tag{18}
```

for one fixed `1/2<=rho<1`. Suppose all coefficient, factor, and
index generation costs are `2^O(ell_j)`. Then there is a uniform family
at every positive integer length `lambda` with

```math
\log_2|I_\lambda|=O((\log(\lambda+2))^2),
\qquad
r+c\le4\rho\lambda+O((\log(\lambda+2))^2),
\tag{19}
```

separate widths at most `2lambda`, and complete factor and metadata
generation in `2^O(lambda)` bit time. It has the same exact crossing sum.

**Proof.** Partition asymptotics give

```math
\ell_j=c\,p^{j/2}+O(j+1).
\tag{20}
```

After removing repeated capacities, there are `O(log lambda)` capacities
at most `lambda`, and the ratio of consecutive capacities is bounded by
a constant depending only on `p`. All finitely many small exceptions
can be included in that constant.

Include a capacity-one block. It can be obtained by restricting any
fixed prime-power family with at least two labels to exactly two labels,
then using the trivial separate width caps. It remains homogeneous and
scale one. Its cost and width error are constants depending only on `p`.
This handles primes whose smallest positive `ell_j` exceeds one.

Greedily take the largest available capacity not exceeding the remaining
length. Before using capacity `a`, the remainder is smaller than the
next larger capacity. The bounded ratio therefore limits the number of
uses of `a` by a fixed constant. The final capacity-one remainder is also
bounded. This produces

```math
\lambda=\ell^{(1)}+\cdots+\ell^{(t)},\qquad
t=O(\log(\lambda+2)).
\tag{21}
```

Split each binary label into these blocks, and retain the **entire tuple**
of local piece indices. For a tuple `i=(i_1,...,i_t)`, define

```math
A_i(x,y,u,v)=\prod_{b=1}^{t}
 A^{(b)}_{i_b}(x_b,y_b,u_b,v_b).
\tag{22}
```

Finite distributivity and the local scale-one identities show that the
sum of (22) is exactly `[x=v][y=u]`. Tensor the row factorizations,
and separately tensor the column factorizations. Their widths add.
Thus (18), (21), and `ell^(b)<=lambda` imply (19).
The separate caps add to `2lambda`.

A local polynomial piece count has `O(log(ell^(b)+2))` index bits.
The complete tuple has `O(log^2(lambda+2))` bits, and the product of the
piece counts is `2^O(log^2(lambda+2))`. No coefficient grouping, shared
parameter, or cancellation among tuple indices is needed.

For uniformity, enumerate the prime powers until the first capacity
larger than `lambda`. By (20), the largest energy considered is
`O_p(lambda^2)`. Integer partition-number dynamic programming finds
every capacity and performs the greedy decomposition in polynomial time.
All local construction times sum to at most
`O(log(lambda+2)) 2^O(lambda) = 2^O(lambda)`.
Each tensor factor query splits its label and channel strings, evaluates
the local factors, and multiplies `O(log(lambda+2))` field values.
Enumerating all tuple indices takes `2^O(log^2(lambda+2))` time.

Even writing **all** coefficient and factor tables remains within
`2^O(lambda)` time. A complete piece has `2^(4lambda)` entries.
Each of its four capped factor tables has at most `2^(4lambda)` entries.
Multiplying by the number of pieces and by the per-entry generation
bound still gives `2^O(lambda)`. Index descriptors, widths, block
boundaries, and any elimination pivots are included in these bounds.
This proves the theorem.

Apply the theorem to the homogeneous full-Schur family of Section 4.
At capacity `ell_j`, (13) holds with `lambda` replaced by `ell_j`;
`floor(log2 p(h_j))` changes its error by at most one.
There are `4h_j+1=poly(ell_j+2)` pieces, and Section 6 proves the
required algorithms. Hence (19) holds with `rho=1/sqrt(2)`.
Both errors in (19) fit the extended transfer theorem; for example,
choose its piece-entropy parameter `delta=1`. It gives the complete
fixed-characteristic factoring rate

```math
T(n)=\exp\!\left(O\!\left(\frac n{\sqrt{\log n}}\right)\right).
\tag{23}
```

The contribution of this composition theorem is broader than (23):
a future uniform improvement of (18) only at prime-power energies
would apply at every length with the same improved `rho`.
The finite ranks below do not establish that hypothesis.

## 6. Effective coefficients and complete construction costs

Set `Delta_h=sum_(E=0)^(4h) p(E)=exp(O(sqrt(h)))`.
All spaces needed for (10) and its two factorizations lie in this range.
A direct algorithm has cost polynomial in `Delta_h` and `h`.
Here is an explicit route that avoids expanding exponentially many
auxiliary-variable monomials in (6).

Use the character change of basis

```math
s_\mu=\sum_{\nu\vdash E}\frac{\chi^\mu(\nu)}{z_\nu}p_\nu,
\qquad
p_\nu=\sum_{\mu\vdash E}\chi^\mu(\nu)s_\mu.
\tag{24}
```

The Murnaghan--Nakayama recurrence computes these integer character
tables in `poly(h,Delta_h)` time. One explicit dynamic program records
a remaining Young diagram and a suffix of the cycle partition.
Enumerate candidate smaller diagrams and test that the difference is
a connected rim hook with no two-by-two square. The state and transition
counts, even when repeated for every cycle partition and energy, are
polynomial in `h` and `Delta_h`. Character magnitudes are at most `E!`,
and every `z_nu` divides `E!`, so (24) uses polynomially many bits in `h`.

For the translated duals, apply the binomial substitution
`T(p_m)=sum_(k=1)^m binom(m,k)p_k` by sparse multiplication on partitions
of weight at most `h`. This uses polynomial time in the same parameters.
It increases the coefficient norm of a single power-sum monomial by
at most `2^h`. No enumeration of all choices of lowered parts is needed.

Now evaluate monomial operator coefficients using the source normal-order
recursion, memoized by an operator-partition suffix and the input and
output partitions. There are at most
`(h+1)Delta_h^2` states for each operator monomial, and polynomially many
transitions per state. Sum the resulting rational power-sum coefficients
and convert back by the second formula in (24).

All exact arithmetic here has polynomial bit lengths in `h`.
The operator has at most `h` parts; mode magnitudes are at most `4h`;
binomial and differentiation coefficients are bounded by `2^O(h)`
per part. Even the crude bound `(8h)^h` for the number of normal-order
terms adds only `O(h log(h+1))` bits. Together with (24), this bounds
every factor by `2^poly(h)` in absolute value. The algorithm uses the
memoized recurrence rather than that crude enumeration.

For computation modulo `p`, perform the indicated exact rational
calculation, use Section 2 to know that the final Schur coordinate is
an integer, and reduce that integer modulo `p`. This method never
reduces a rational denominator that is divisible by `p`.
For homogeneous duals, first compute the integral matrix of `B_h`,
reduce it, and invert it over `F_p`. Elimination takes polynomial time
in `p(h)`. Lift the inverse's Schur coefficients to `0,...,p-1` and
use the same operator algorithm. With fixed `p`, these coefficients
again have polynomial bit lengths in `h` throughout the rational
coordinate calculations.

Partition lists, label selection, intermediate channel lists, the angle
binomials, and the capping factorizations in Section 3 are all covered by
these algorithms. Since `h=O(lambda^2)` in Section 3, or
`h=O((ell_j+1)^2)` for the full prime-power blocks, their complete bit
cost is `2^O(lambda)` or `2^O(ell_j)`, respectively. The constants are
uniform in the label length and may depend on the fixed prime.

## 7. Exact checks and the remaining target

The companion checker
[`check_integral_schur_crossing.py`](../tools/check_integral_schur_crossing.py)
computes Schur coordinates through the separate power-sum recursion.
Its evidence is
[`integral-schur-crossing.json`](../evidence/integral-schur-crossing.json).
The prime-power capacity algorithm is checked independently by
[`check_prime_power_capacities.py`](../tools/check_prime_power_capacities.py),
with evidence in
[`prime-power-capacities.json`](../evidence/prime-power-capacities.json).
For primes two, three, five, and seven and lengths `1,...,256`, it checks
1,024 decompositions, 135,656 binary label encodings, and 1,263,552 piece
index encodings. Independent partition dynamic programs agree through
weight 128. Negative controls reject a missing block and, for primes
five and seven, the omission of the restricted capacity-one family.

Run:

```sh
python3 tools/check_integral_schur_crossing.py
python3 tools/check_prime_power_capacities.py
```

The recorded Schur evidence also includes the larger central-piece check:

```sh
python3 tools/check_integral_schur_crossing.py \
  --homogeneous-eight-central --output evidence/integral-schur-crossing.json
```

It records 4,896 individual integral angle factors, 552 dualities,
256,852 complete-piece integrality entries, and 11,971 comparisons of
the two operator orders. The default quick run omits the weight-eight
central piece and its dualities.

At weights two and four over `F_2`, the full homogeneous Schur family has
dimensions two and five, respectively. Its central row/column ranks are
`(2,2)` and `(12,12)`; the maximum complete-piece rank products are
four and 144. At weight eight, the 22-label homogeneous central piece
has ranks `(228,228)`; only that shift was checked at this weight.
These are finite exact facts, not asymptotic rank bounds.
The general algebra and cost claims above are proved independently of
those ranks.

The research target is now precise: prove a uniform joint rank bound
for the complete homogeneous pieces at prime-power weights that makes
(18) hold with `rho<1/sqrt(2)`. The full Schur lattice changes the
modular representation, and Section 5 removes the separate obligation
to match every label length with a single prime power. Neither change
settles the rank bound.

The construction, proof, and exact checks used substantial OpenAI
assistant support. They retain the source attribution and repository
license. No new Lean formalization is claimed or required here.
