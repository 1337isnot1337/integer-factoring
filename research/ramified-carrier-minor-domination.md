# Rank requirements for perturbing a ramified crossing

## Results and scope

A low-rank constant perturbation cannot improve the local weighted Smith
cost of a ramified matrix when the ramification separates its determinant
terms by their orders modulo the ramification index. For the degree-sixteen
Bell carrier family, a perturbation of cut rank at most seven cannot improve
any determinantal prefix at either existing finite place. For every
perturbation, its first seven prefixes cannot improve.

For the alternating four-letter Bell form, this combines with a rank and
degree argument to exclude every **unweighted** improvement at either
existing place, including full-rank carrier changes. The statement concerns
the two specified local profiles. Physical weights for full-rank changes,
additional singular places, and other crossing families remain outside it.

## 1. A residue lemma for minors

Let `F` be any field and let

```math
M_0(t)=A(t^m),\qquad M(t)=M_0(t)+\eta(t)G,
\qquad \operatorname{ord}_t\eta=r,\quad \gcd(r,m)=1,
```

where `A` is a square power-series matrix, `G` is a constant matrix, and
`m>=2`. Suppose `rank(G)<m`. Then every fixed square minor satisfies

```math
\operatorname{ord}_t\det M_{I,J}
\le \operatorname{ord}_t\det (M_0)_{I,J}.                 \tag{1}
```

The order of a zero determinant is infinity. To prove (1), expand the minor
as a polynomial in the update scalar:

```math
\det M_{I,J}=\sum_{j=0}^{\operatorname{rank}G_{I,J}}
                      \eta(t)^j p_j(t^m).                \tag{2}
```

The rank bound on the update gives the stated upper limit. Each nonzero
summand has order

```math
rj+m\operatorname{ord}_z p_j(z).
```

These orders have distinct residues modulo `m`, since `0<=j<m` and
`r` is coprime to `m`. The term of minimum order cannot cancel. Thus the
order of (2) is its least summand order, which is at most the order of its
`j=0` term. If that term is zero, the inequality is automatic.

Without a rank assumption on `G`, the same argument still applies to
every minor of size less than `m`.

## 2. Physical weights and all local Chernoff sums

Apply arbitrary rational row and column monomial weights. A fixed minor
acquires the same additive weight shift on both sides of (1). Taking the
minimum over all minors of size `k` gives

```math
\Gamma_k(M)\le\Gamma_k(M_0),                             \tag{3}
```

for the corresponding weighted determinantal orders. Rational weights can
be interpreted over one fixed ramified extension of the local parameter;
the comparison does not depend on the chosen denominator.

When `M_0` is generically invertible, (1) for the full determinant also
proves that `M` is generically invertible.

Write the complete ordered Smith profiles as `sigma` and `sigma_0`.
Their prefix sums are the `Gamma_k`. For every real `a`,

```math
\sum_i(a-\sigma_i)_+
 =\max_{0\le k\le n}\{ka-\Gamma_k\}.
```

Consequently (3) increases every lower-tail hinge sum. Integrating those
hinges against `s^2 exp(-sa)` proves, for every `s>0`,

```math
\sum_i e^{-s\sigma_i(M)}
\ge \sum_i e^{-s\sigma_i(M_0)}.                           \tag{4}
```

This conclusion includes arbitrary physical monomial heights. The common
type shift at a fixed height is the same for both tensors, so it does not
change the comparison. The lemma compares local Smith bounds; it makes no
claim about actual tensor-power jet ranks below those bounds.

The condition `rank(G)<m` is necessary in general. For example, take

```math
M_0(t)=\begin{pmatrix}1&0\\0&t^2\end{pmatrix},\quad
G=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
\eta=t+t^2
```

in characteristic two. The determinant order rises from two to four.
Here the update rank equals the ramification index two, so the two
leading determinant terms can cancel.

## 3. The seventh-order Bell carrier family

Use the [opposite Bell construction](entrywise-reflection-jet-interpolation.md)
and its [unique eighth-order carrier](bell-jet-transversality.md). Write

```math
S_V=S_0+K G,\qquad G=V-V_0,\qquad
K=(w+1)^2(w^2+w+1)^7.
```

The baseline `S_0` is a matrix polynomial in `w^8`. At either designated
point, `w=a+t` gives a polynomial in `t^8`, while `K` has order seven.
Thus the residue lemma applies with `m=8,r=7`.

At either cut, an update of rank at most seven cannot improve any weighted
Smith prefix or local Chernoff sum. At every update rank, the first seven
weighted prefixes cannot improve. This excludes low-rank updates as a way
to improve these two existing local bounds. A new singular place would
require a separate complete profile and reconstruction allocation.

## 4. Complete unweighted exclusion for the alternating form

For an alternating invertible four-by-four relative Bell form, the
unweighted baseline at either designated cut has the complete profile

```math
0,\quad 8^{\times4},\quad16^{\times6},\quad24^{\times4},\quad32.
```

Here one first uses a paired constant physical change to normalize the
column Bell factor to the identity. The remaining changes that preserve
that pairing and the crossing act on the row factor by `A -> X A X^T`.
Every nondegenerate alternating form is congruent to the canonical
symplectic form. These physical changes preserve both unweighted local
Smith profiles, so the checked canonical profile applies to every such
relative form.

Its local Chernoff sum is

```math
Z_0(s)=(1+e^{-8s})^4.                                    \tag{5}
```

Let `G` be any allowed trace-carrier update. If its cut rank is at most
seven, Section 2 proves the desired comparison. Otherwise its projection
to the fifteen-by-fifteen Schur complement of the rank-one constant
matrix has rank at least `rank(G)-2>=6`. The seventh coefficient of that
Schur complement is a nonzero scalar times this projection. Therefore at
least six Smith orders are exactly seven.

The constant matrix has rank one, so there is one order-zero channel.
The determinant has degree at most `16*16=256`, and hence the sum of all
local orders is at most 256. After selecting the unit channel and six
order-seven channels, the remaining nine orders sum to at most 214.
Convexity and monotonicity give

```math
Z(s)\ge1+6e^{-7s}+9e^{-214s/9}\ge(1+e^{-8s})^4.           \tag{6}
```

The second inequality has an exact integer certificate. Set
`t=exp(-s/9)`. The difference of its two sides is

```math
P(t)=6t^{63}+9t^{214}-4t^{72}-6t^{144}-4t^{216}-t^{288}
    =t^{63}(1-t)^2 Q(t).
```

The polynomial `Q` has degree 223. All 224 coefficients of

```math
(1+y)^{223}Q\left(\frac{y}{1+y}\right)
```

are strictly positive integers; the smallest is six. Since `0<t<1`
corresponds to `y=t/(1-t)>0`, this proves the second inequality in (6).
The endpoint limits follow by continuity.

Thus no constant trace carrier improves either designated **unweighted**
local Chernoff sum for the alternating form. Within this family the common
scalar content is one: both paired traces equal `1+w^16`, and the strict
endpoint excludes its only root. There is no omitted scalar division that
would lower the degree charge in this comparison.

## 5. Exact controls

The [checker](../tools/check_ramified_carrier_minors.py) enumerates all minors
of small ramified matrices, checks weighted complete profiles of literal
Bell carrier perturbations, verifies the rank-threshold counterexample,
and constructs every integer in the positivity certificate. The
[record](../evidence/ramified-carrier-minor-controls.json) separates these
finite controls from the uniform proofs above.

```sh
python3 tools/check_ramified_carrier_minors.py \
  --output evidence/ramified-carrier-minor-controls.json
```
