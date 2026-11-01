# Finite crossing rank profiles can be realized over finite fields

## Lemma

Let a finite list of tensors `A_1,...,A_M` over an arbitrary field `K`
sum to the scale-one crossing tensor on `d` labels. For each tensor,
fix finitely many matrix flattenings and record their ranks over `K`.
There is a finite field and a crossing decomposition with exactly the
same number of pieces and exactly the same recorded ranks.

In particular the infimum of either the rational finite-seed quality
or the rank-sum quality over all fields equals its infimum over finite
fields. This is an application of ordinary finite-type specialization;
it is not a new tensor decomposition or a numerical rank improvement.

The proof is supplied in full because its use permits the field to
depend on the finite seed. A separate diagonal-selection theorem is
needed to choose one field uniformly for each entire factoring run.

## 1. Preserve identities and rank witnesses in a finitely generated ring

Discard tensors that are zero before recording `M`. For each nonzero
flattening choose one nonzero minor of size equal to its rank. Let
`delta_1,...,delta_s` be these determinants. Define the subring of `K`

```math
R=\mathbb Z[\text{all tensor entries},\delta_1^{-1},\ldots,
\delta_s^{-1}],
```

where integers are mapped into `K`. The ring is nonzero, commutative,
unital, and finitely generated as a `Z`-algebra. Choose a maximal ideal
`m` and put `F=R/m`. Section 2 proves that `F` is a finite field.

Map every tensor entry to `F`. Each crossing identity is a polynomial
identity with integer coefficients, so it remains exact. In particular
the scalar one remains one, and the sum is still the crossing tensor.

Suppose a recorded flattening had rank `r`. Every `(r+1)`-minor was
zero in the subring `R`, so these minors remain zero after reduction;
the new rank is at most `r`. Its selected `r`-minor is a unit of `R`,
since its inverse was adjoined, and therefore maps to a nonzero element
of the residue field. The new rank is at least `r`. This proves exact
rank preservation for every chosen cut and keeps every retained piece
nonzero.

One may omit the inverse determinants if only nonincreasing ranks are
needed. The inverse construction is useful because it proves equality
of attainable finite rank profiles, rather than just an inequality
between their infima.

## 2. A field finitely generated as a `Z`-algebra is finite

Let `F=Z[alpha_1,...,alpha_t]` be a field. If it had characteristic
zero, it would contain `Q` and hence also satisfy
`F=Q[alpha_1,...,alpha_t]`. Zariski's lemma, proved in Section 3,
then makes it a finite algebraic extension of `Q`.

Choose a positive integer `D` such that every `D alpha_i` is an
algebraic integer. Such a `D` exists: if an algebraic number satisfies
an integer polynomial with leading coefficient `a`, then multiplying
that number by `a` gives an algebraic integer, by directly rescaling
the polynomial to a monic one. Take a product of the finitely many
leading coefficients. Consequently

```math
\mathbb Z[\alpha_1,\ldots,\alpha_t]
 \subseteq \mathcal O_F[1/D],
```

where `O_F` is the ring of algebraic integers in `F`.

Choose a prime `q` not dividing `D`. Since the left side is a field,
it contains `1/q`. The inclusion would give `1/q=beta/D^k` for some
algebraic integer `beta`. But then `beta=D^k/q` is rational. Every
rational algebraic integer is an integer: clearing a monic equation
for a reduced fraction shows its denominator divides a power of its
coprime numerator. This would force `q | D`, a contradiction.

Thus `F` has positive characteristic `p`. It is a finitely generated
algebra over the finite field `F_p`. Zariski's lemma makes it a finite
algebraic extension of `F_p`, so it has finitely many elements. Applying
this argument to the maximal-ideal quotient in Section 1 proves the
claim used there.

## 3. Zariski's lemma in the needed form

Suppose a field `L=k[alpha_1,...,alpha_t]` is finitely generated as
an algebra over a field `k`. Choose a maximal algebraically independent
subset `x_1,...,x_r` of the generators. The remaining generators are
algebraic over `k(x_1,...,x_r)`.

Clear the finitely many denominators in their monic algebraic equations.
For some nonzero polynomial `g`, all the remaining generators are
integral over

```math
B=k[x_1,\ldots,x_r,1/g].
```

They generate `L` over `B`: the original algebra already equals `L`,
and `1/g` also belongs to the field `L`. Hence `L` is integral over
`B`. A subring over which a field is integral must itself be a field.
Indeed, for nonzero `b in B`, the element `b^(-1)` satisfies a monic
equation over `B`; multiplying that equation by a suitable power of
`b` expresses `b^(-1)` as an element of `B`.

If `r>=1`, the displayed `B` is not a field. Choose an irreducible
polynomial `h(x_1)` which does not divide `g`. There are infinitely
many irreducibles in `k[x_1]` by the usual product-plus-one argument,
while one nonzero coefficient of `g`, viewed as a polynomial in the
other variables, has only finitely many irreducible factors. The ideal
generated by `h` remains proper after localizing at `g`: its quotient
is a polynomial ring over the field `k[x_1]/(h)`, localized at a
nonzero element. Therefore `h` is not a unit of `B`, a contradiction.

It follows that `r=0`. All the generators are algebraic over `k`, and
finitely many algebraic generators give a finite extension. This proves
the form of Zariski's lemma used above.

## 4. Consequence for uniform finite-seed searches

Every finite-field decomposition is already a decomposition over a
field. Conversely, Sections 1--3 preserve the finite ranks of every
seed over any coefficient field in some finite field. Thus both
collections attain exactly the same finite rank profiles, including
all row/column rank data used by the seed-amplification criterion.

A universal finite-field enumerator can therefore search this entire
space of profiles. It need not be given, compute with, or even decide
equality among the entries of the original arbitrary field. The
specialization argument asserts existence of a finite certificate;
the enumerator will eventually print that certificate in its ordinary
finite encoding.

This is an asymptotic existence and uniformity conclusion. The residue
field may be large, and the proof gives no useful numerical threshold
for when its certificate enters a slowly growing bounded search. The
field must remain fixed within each contraction; changing characteristic
between recursive levels does not preserve the tensor identity.

## 5. Extension fields do not improve the fixed-characteristic infimum

Let `F=F_(p^e)` and choose an `F_p`-linear map
`ell:F -> F_p` with `ell(1)=1`. Such a map exists by extending the
nonzero vector one to a basis. Apply `ell` entrywise to every piece of
a scale-one crossing decomposition. The crossing sum remains exact,
because its entries are zero and one.

If a flattening `A` has rank `r` over `F`, write `A=UV` through
`r` channels. In a basis `b_1,...,b_e`, expand

```math
U_{ik}=\sum_{j=1}^e u_{ikj}b_j,\qquad u_{ikj}\in F_p.
```

Then

```math
\ell(A_{ih})=
\sum_{k=1}^r\sum_{j=1}^e u_{ikj}\,\ell(b_jV_{kh}).
```

This is an explicit prime-field factorization through `er` channels.
The same argument applies separately to both cuts. Thus projected
piece ranks satisfy `R_i'<=eR_i` and `C_i'<=eC_i`, after discarding
any zero projected pieces.

For a width-`a` seed `S`, put

```math
\mathcal S(S)=\sum_i\sqrt{R_iC_i},\qquad
\rho(S)=\frac{\log_2\mathcal S(S)}{2a}.
```

Projection therefore increases its rank sum by at most a factor `e`.
To remove this fixed factor asymptotically, first form the width-`ka`
type amplification of `S` over `F`, and only then apply `ell`.
A type with multiplicities `nu_i` has multinomial count `D_nu` and

```math
\sqrt{R_\nu C_\nu}
\le D_\nu\prod_i(\sqrt{R_iC_i})^{\nu_i}.
```

Summing over types gives rank sum at most `mathcal S(S)^k` by the
multinomial theorem. Consequently the prime-field projected seed has

```math
\rho(S_k')\le\rho(S)+\frac{\log_2 e}{2ak}.             \tag{1}
```

The map is applied once to each amplified tensor entry, so the loss is
`e`, not `e^k`. Projecting each original block before amplification
would not give (1).

Every finite extension seed can thus be approximated arbitrarily well
in rank-sum quality by finite prime-field seeds. The reverse inclusion
is immediate. Hence at each fixed characteristic

```math
\inf_{\text{seeds over }F_p}\rho
=\inf_{\text{finite extensions of }F_p}\rho.           \tag{2}
```

For the rational selector quality, a type seed has at most `(k+1)^M`
pieces. The same estimates give a largest rounded rank-width sum at
most `2k log_2 mathcal S(S)+2 log_2 e+2`; the additional piece-count
term is `O(M log(k+1))`. Dividing by `4ak` shows that its rational
quality converges from above to at most `rho(S)` as well. No numerical
logarithm oracle is required by a selector using exact ranks and bit
lengths.

Combining Sections 1--4 with (2), the rank-sum infimum over all
coefficient fields equals its infimum over prime fields. A universal
score selector may therefore restrict its field records to primes.
This does not identify one universally sufficient characteristic:
the prime may still vary with the chosen finite seed.

## 6. Commutative rings preserve factor-channel upper bounds

Suppose instead that a finite scale-one crossing identity is given over
any nonzero commutative unital ring `A`. Assume its row and column
factorizations through respectively `r_i` and `c_i` channels are
explicit finite matrix identities. Generate a `Z`-subalgebra by every
tensor entry and every factor entry appearing in these identities.
It remains nonzero, because it contains the same zero and one as `A`.

Take a maximal ideal and its finite-field residue quotient as in
Section 2. Reduction preserves both the crossing identity and both
factorizations. The resulting field ranks are therefore at most
`r_i` and `c_i`, with zero pieces discarded. In this finite-certificate
interface, a commutative coefficient ring cannot give smaller channel
bounds than are achievable over some finite field, and prime-field
seeds approximate the same rank-sum quality after Section 5.

This statement concerns supplied factor dimensions. Matrix rank over
a ring need not have the field-rank properties used earlier. In
particular, nonzero zero-divisor minors cannot automatically be
inverted while keeping a nonzero ring, so exact rank preservation is
not asserted in this ring version.

An invertible common crossing scalar can first be normalized without
changing the factor dimensions. A nonunit crossing scalar is outside
this argument: its image can vanish in a residue field, and the reduced
sum may then be the zero tensor instead of a crossing. Any algorithm
using such a scalar needs its own precision and division accounting.
