# Homogeneous matching duals from Jack integrality

## Result

Let `M_(2ell)` be the perfect matchings of `2ell` bit positions. Encode an
edge `e={i,j}` by the odd mode `n_e=1+2(2^i+2^j)`, as in the
[matching-code lift](matching-code-lift.md). Let `G_ell` be the homogeneous
vacuum pairing on these label monomials, over `F_2`.

For every **even** `ell`, this pairing is invertible. More precisely, if
`p(ell)` denotes the ordinary partition number, then

\[
 (G_\ell-I)^{p(\ell)}=0.
 \tag{1}
\]

For odd `ell`, the corresponding statement is

\[
 G_\ell^{p(\ell)}=0.
 \tag{2}
\]

These are all-size statements. Their proof uses an integral annihilating
polynomial before reduction modulo two. It does not assume that an
eigenbasis survives reduction. This constructs homogeneous unit duals on
every even-particle matching space and an effective multigraded crossing
family with `exp(O(lambda/log lambda))` pieces. No improved joint width
bound, or factoring exponent beyond `1/2`, follows here.

## 1. The pairing as an association-scheme matrix

For matching labels `M,N`, the vacuum pairing is

\[
 G_\ell(M,N)=\det\bigl([e\cap f=\varnothing]\bigr)_{e\in M,f\in N}
 \quad\text{in }\mathbf F_2.
 \tag{3}
\]

The formula follows from the source normal ordering and Lucas' theorem.
Every mode is odd and appears once, so the normalization factors are units
and equal to one modulo two.

Write `A_k` for the integer zero-one adjacency matrix whose entries are
one when `M union N` consists of one nontrivial alternating cycle using
`k` edges of each matching, with all other edges common. The determinant
calculation in [the lift note, Section 5](matching-code-lift.md#5-the-homogeneous-pairing-on-matching-labels)
gives

\[
 G_\ell=[\ell\text{ even}]I+N_\ell\pmod2,
 \qquad
 N_\ell=\sum_{\substack{3\le k\le\ell\\k\text{ odd}}}A_k.
 \tag{4}
\]

The matrices `A_k` are members of the perfect matching association
scheme. Over characteristic zero, their common eigenspaces are indexed by
partitions `lambda` of `ell`; in particular there are at most `p(ell)`
distinct simultaneous eigenvalue lists. This is the usual multiplicity-free
decomposition of the matching permutation module; see
[Srinivasan, Section 1 and Theorem 3.1](https://arxiv.org/pdf/1807.00481).

## 2. A parity lemma for the eigenvalues

Use the integral Jack normalization `J_lambda^(alpha)`, in which the
coefficient of the squarefree monomial of degree `ell` is `ell!`. The
eigenvalue of `A_k` on the component indexed by `lambda` is exactly

\[
 \theta_{\lambda,k}
  =[p_kp_1^{\ell-k}]J_\lambda^{(2)}.
 \tag{5}
\]

For the normalization check, the valency of the relation of cycle type
`mu` is `|H_ell|/z_(2mu)`. Its eigenvalue is that valency times the zonal
spherical function. The same coefficient appears in the power-sum
expansion of the integral zonal polynomial. See
[Renteln, equations (5)–(6), and the normalization discussion after (13)](https://cayley.academic.csusb.edu/content/papers/pmdg.pdf).
At parameter two, this zonal polynomial is `J_lambda^(2)`.

Two standard Jack facts will be used:

1. Every power-sum coefficient of `J_lambda^(alpha)` belongs to
   `Z[alpha]`.
2. If the column lengths of `lambda` are `r_1,...,r_t`, then
   `J_lambda^(0)=prod_j r_j! e_(r_j)`.

For the first fact, [Knop–Sahi, Theorem 1.1](https://arxiv.org/pdf/q-alg/9610016)
gives integral polynomial coefficients in augmented monomials. Each
augmented monomial is integral in the power sums: inclusion-exclusion on
equal indices gives

\[
 \widetilde m_{(b_1,\ldots,b_s)}
 =\sum_{\pi\in\Pi_s}
   \prod_{B\in\pi}(-1)^{|B|-1}(|B|-1)!\,
   p_{\sum_{i\in B}b_i}.
 \tag{6}
\]

This proves the required power-sum integrality. It is also stated
explicitly in [Mickler–Moll, Theorem 4.1](https://arxiv.org/pdf/2211.01586).
For the second fact, see
[Dołęga–Féray, equation (12)](https://arxiv.org/pdf/1601.01501), citing
Stanley's specialization formula.

By integrality, evaluation at `alpha=2` and `alpha=0` gives the same
power-sum coefficients modulo two. In one factor `r!e_r`, the coefficient
of `p_k p_1^(r-k)` is

\[
 (-1)^{k-1}\frac{r!}{k(r-k)!}
 =(-1)^{k-1}\binom rk(k-1)!.
\]

The coefficient of `p_1^r` is one. Thus (5) satisfies

\[
 \theta_{\lambda,k}
 \equiv(-1)^{k-1}(k-1)!
       \sum_j\binom{r_j}{k}\pmod2.
 \tag{7}
\]

Consequently **every eigenvalue of `A_k` is an even integer whenever
`k>=3`**. The exception `k=2` matters: its eigenvalues need not be even.

There is a useful stronger parity identity. Put
`b(lambda)=sum_j floor(r_j/2)=lambda_2+lambda_4+...`. Then

\[
 J_\lambda^{(2)}
 \equiv p_1^{\ell-2b(\lambda)}(p_1^2+p_2)^{b(\lambda)}\pmod2.
 \tag{7a}
\]

To see this, write `R_r=r!e_r`. Newton's identity, multiplied by `(r-1)!`,
shows that modulo two
`R_r=p_1 R_(r-1)+(r-1)p_2 R_(r-2)`: every term with `p_k`, `k>=3`, has
a product of at least two consecutive integers as its coefficient. Hence
`R_(2s)=(p_1^2+p_2)^s` and `R_(2s+1)=p_1(p_1^2+p_2)^s`. Multiplying over
the columns proves (7a). In particular, every orbital whose cycle type
has a part at least three has only even eigenvalues. The remaining
eigenvalue parities, for type `(2^j,1^(ell-2j))`, are `binom(b(lambda),j)`.

## 3. Reduction of an integral annihilator

Let `eta_lambda` be the eigenvalue of the integer matrix `N_ell` in (4)
on the `lambda` component. Equation (7) shows that all `eta_lambda` are
even integers. Since `N_ell` is diagonalizable over characteristic zero,
the monic polynomial

\[
 F_\ell(X)=\prod_{\lambda\vdash\ell}(X-\eta_\lambda)\in\mathbf Z[X]
\]

annihilates it. Repeated eigenvalues cause no problem. The matrix equation
`F_ell(N_ell)=0` is an equality of integer matrices. Reducing it modulo two
gives

\[
 \overline N_\ell^{p(\ell)}=0,
\]

because every factor of `F_ell` reduces to `X`. Equations (1)–(2) now
follow from (4).

For even `ell`, the inverse is explicitly

\[
 G_\ell^{-1}=\sum_{j=0}^{p(\ell)-1}(G_\ell-I)^j.
 \tag{8}
\]

Equivalently, take a power of two `R>=p(ell)`. Then `G_ell^R=I` and
`G_ell^-1=G_ell^(R-1)`, so repeated squaring constructs the inverse. The
partition-number bound makes `log R=O(sqrt(ell))`. This argument establishes
invertibility without asserting that `G_ell^2=I`: that stronger statement
already fails at `ell=6`.

## 4. An effective homogeneous multigraded family

Take even `ell`, put `a=2ell`, and retain all matching monomials `b_M`.
Define

\[
 a_M=\sum_N(G_\ell^{-1})_{MN}b_N.
 \tag{9}
\]

All terms have common weight
`h=ell+2(2^a-1)`. Equations (3) and (9) give unit homogeneous vacuum
duality.

For clarity, one can check the multigraded crossing directly. Put
`B(e)=prod_(i in e)(1+Z_i)` and use group order `U,B,A,Y`. Before applying
the inverse in (9), the monomial Wick blocks are

| Groups | Entry for the edge blocks `e,f` in group order |
|---|---|
| `U,B` | `[e=f](1+L) B(e)` |
| `U,A` | `[e=f] L(1+L) Z^e B(e)` |
| `U,Y` | `[e=f] L Z^e` |
| `B,A` | `[e intersection f is empty] L Z^f` |
| `B,Y` | `[e intersection f is empty] L(1+L) Z^f B(e)` |
| `A,Y` | `[e intersection f is empty] (1+L) B(e)` |

Compute their Pfaffian by the division-free matching expansion. At all
variables equal to one, only `B,A` and `U,Y` survive. Their value is
`G_ell(x,v)[y=u]`. Applying `G_ell^-1` to the `x` label therefore gives
exactly `[x=v][y=u]`. This verifies the crossing evaluation independently
of any cancellation or clearing of denominators in a one-variable formula.

Each bit still occurs once in each of the four external matching labels.
The extra `Z_i` in the `U,A` block is charged to its second occurrence in
the two equal edge labels. Thus each high degree is at most four and the
low degree is at most `4ell`. Applying (9) preserves those bounds term by
term. The family has at most `(4ell+1)5^(2ell)` pieces.

The common invertible change on the `x` label preserves **both**
flattening ranks of every coefficient tensor. This permits exact rank
experiments on the raw monomial Pfaffians; no inverse is needed for their
rank computations.

For a binary label space, choose the least `a` divisible by four such
that `(a-1)!!>=2^lambda`. The overshoot is polynomial in `a`. Construct the
full inverse (9), then restrict all four external labels to the first
`2^lambda` labels. Restriction preserves the crossing and cannot increase
either rank. The entry algorithm, all coefficient tables, the inverse,
and Gaussian factors have total construction cost `2^O(lambda)`, and

\[
 \log_2\#\text{pieces}=O(\lambda/\log\lambda).
\]

No polynomial-space claim is made for these generic inverses or factors.
The missing ingredient for a stronger factoring bound is a uniform upper
bound on the sum of the two logarithmic ranks.

## 5. A collapsed curve and its complete finite profiles

Set all high variables equal to `Z` and put `L=Z^2`. In the preceding
table the two scalars become

\[
 u=LZ^2=Z^4,\qquad
 v=(1+L)(1+Z)^2=1+Z^4.
\]

For any Wick matching, the four equal group sizes force

\[
 \#(UB)=\#(AY)=i,\quad
 \#(UA)=\#(BY)=j,\quad
 \#(UY)=\#(BA)=k,\quad i+j+k=\ell.
\]

Its scalar is `u^(2(j+k))v^(2(i+j))`. Therefore the polynomial is a
polynomial in `t=Z^8`, of degree at most `2ell`, with exact crossing value
at `t=1`. This gives only `2ell+1` coefficient pieces.

For `a=4`, there are three nonzero pieces, of degrees `0,2,4`, with rank
pairs `(1,9),(4,4),(9,1)`. Their maximum product is `16`. For `a=8`, the
complete profile is:

| Degree | Row rank | Column rank |
|---:|---:|---:|
| 0 | 1 | 11,025 |
| 1 | 14 | 2,730 |
| 2 | 453 | 6,174 |
| 3 | 3,105 | 3,463 |
| 4 | 7,328 | 7,328 |
| 5 | 3,463 | 3,105 |
| 6 | 6,174 | 453 |
| 7 | 2,730 | 14 |
| 8 | 11,025 | 1 |

The latter label space has dimension `105`. Its maximum product is
`53,699,584`, so the attractive smaller profile does not persist at this
next even particle count. These finite computations make no asymptotic
rank claim.

## 6. Reproduction and limits

The checker `tools/check_matching_gram_unipotence.py` compares the cycle
formula with independent disjointness determinants, tests the predicted
nilpotence bound and the inverse on small complete matching spaces, and
checks the homogeneous lift against the unfurled monomial Wick formula.
It also records `G_6^2 != I` as a negative control, using exact cycle-orbit
convolution. The all-size result is the proof above, not an extrapolation
from those finite tests.

```sh
python tools/check_matching_gram_unipotence.py --collapsed-eight --output evidence/matching-gram-unipotence.json
```

Omitting `--collapsed-eight` runs the smaller checks. The full run uses
only the Python standard library. It computes entries on 4,154 orbits of
triples under the stabilizer of a fixed eight-vertex matching, then
reconstructs the full matrices by label permutations. Binary block
transposes produce the exact `(x,u) | (v,y)` flattening. Direct enumeration
at four vertices checks that orbit reconstruction and flattening.
