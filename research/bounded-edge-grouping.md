# An explicit bounded-edge grouping specialization

## Scope and prior result

This is an elementary derivation of the bounded-edge grouping route already
formalized in `anthropics/formal-math`, commit
`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`. The public statement is
`AppRunningTime.SomeGrouping`, its hypotheses are `AppRunningTime.Cases`,
and the theorem is `AppRunningTime.someGrouping`. The case-II proof is in
`eth/proof/rates/GeneralPieces/Group/CaseII.lean`; the assembled grouping
algorithm and value theorem are in `.../Group/Lemma.lean`.

The purpose here is to specify simple integer parameters and the exact
floor-compatible inequality needed by this factoring repository. This note
does not construct or certify a changed factoring machine. The prior formal
result and its stronger subsequent half-pieces algorithm must be credited;
the base-two grouping route is not a new discovery of this repository.

The separate [certified k6 machine result](certified-k6.md) has exponent
`alpha = 1/30` and uses `G = 5K^6`, `K_next = 25K^12`. It does not use the
asymptotically degree-two schedule discussed here. The stronger rates for
this schedule remain conditional on a uniform evaluator and its complete
machine-time proof. The [Appendix B route](appendix-b-factoring-route.md)
to `alpha = 1/2` has additional integration obligations of its own.

The factoring definitions cited below are pinned to
`IsaacRudich/SubexponentialDeterministicFactoring` at
`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`.

## Theorem with integer parameters

Let `K,D>=2`, put `L=bit_length(K)+bit_length(D)+8`, and let

```math
 g=144D^2KL.
```

Suppose `w` is a symmetric `g`-by-`g` matrix of nonnegative integer weights,
with zero diagonal and `w_uv<=B` for some positive integer `B`. Write
`d_u=sum_v w_uv`, `U=sum_u d_u`, and choose `P=gB` (or any larger bound).
Then a coloring into `K` groups exists with

```math
 Ks_a\le2g,\qquad DK^2Y_{ab}\le DU+gP\quad(a\ne b),
```

where `Y_ab` counts each edge between distinct groups once.

Unlike the dense theorem in `dense-grouping.md`, this statement does not
assume a positive lower bound on the edge weights. The price is the
explicit bound `P>=gB`; it cannot be replaced silently by the actual maximum
degree when that degree is much smaller than `gB`.

### Additive weighted tail bound

For independent Bernoulli variables, weights `0<=a_i<=A`, positive `A`,
`Z=sum a_i X_i`, and `mu=E Z`, the same moment generating function proof
as in the dense note gives, for every positive `t`,

```math
 \Pr(Z\ge\mu+t)
 \le\exp\!\left(-\frac{t^2}{2A(\mu+t)}\right).
 \tag{1}
```

When `mu>0`, take the exponential parameter `log(1+t/mu)` and use
`(1+x)log(1+x)-x>=x^2/(2(1+x))`. When `mu=0`, nonnegative `Z` is zero
almost surely, so the claimed upper-tail event has probability zero.

### Simultaneous balances

Color each vertex independently and uniformly. For every vertex `u` and
color `b`, put `X_ub=sum_{chi(v)=b}w_uv`, of mean `d_u/K`. Apply (1) with
`A=B`, `t=gB/(4DK)`. Since `d_u<=gB` and `D>=2`,

```math
 \Pr\left(X_{ub}>\frac{d_u}{K}+\frac{gB}{4DK}\right)
 \le\exp\!\left(-\frac{g}{36D^2K}\right).
 \tag{2}
```

Indeed, `2A(mu+t)<=2B(gB/K)(1+1/(4D))<=9gB^2/(4K)`.

For every color `a`, let `Z_a=sum_{chi(u)=a}d_u`, of mean `U/K`.
Its weights are at most `gB`; applying (1) with `A=gB` and
`t=g^2B/(4DK)` gives

```math
 \Pr\left(Z_a>\frac UK+\frac{g^2B}{4DK}\right)
 \le\exp\!\left(-\frac{g}{36D^2K}\right).
 \tag{3}
```

The usual binomial size bound is
`Pr(s_a>2g/K)<=exp(-g/(3K))`, no larger than the bounds in (2)--(3).
These are marginal tail bounds; no independence between the different
events is assumed. The union has `(g+2)K` events.

The common exponent is `g/(36D^2K)=4L`. Also

```math
 (g+2)K\le2gK=288D^2K^2L<2^{3L-7}<e^{3L}.
```

The first inequality uses `g>=2`, and the next uses `288<2^9`,
`K<2^{bit_length(K)}`, `D<2^{bit_length(D)}`, and `L<=2^L`.
The union failure probability is therefore less than `e^{-L}<1`.

For one coloring satisfying all balances and for distinct `a,b`,

```math
 \begin{aligned}
 Y_{ab}
 &=\sum_{u:\chi(u)=a}X_{ub}\\
 &\le\frac{Z_a}{K}+s_a\frac{gB}{4DK}\\
 &\le\frac{U}{K^2}+\frac{g^2B}{4DK^2}
       +\frac{2g}{K}\frac{gB}{4DK}\\
 &=\frac{U}{K^2}+\frac{3g^2B}{4DK^2}
 \le\frac{U+gP/D}{K^2}.
 \end{aligned}
```

This proves the theorem, including strict failure inequalities and the
ordered-total versus unordered-cross-pair convention. If `B=0`, every edge
vanishes and an ordinary balanced coloring handles the case separately.

## Deterministic use, source mapping, and rate

Enumerate the `K^g` colorings, test the final inequalities with integers, and
stop at the first valid coloring. This is uniform and has cost
`K^g poly(g,K,log(B+1),log(D+1),log(P+1))` bit operations. The probabilistic
argument is only an existence proof.

The pinned reference's `FactoringAux/Level.lean` proves
`rowNet_c_le` and `colNet_c_le`, namely `w_uv<=G^2*lambda`, and
`rowNet_fits` and `colNet_fits`, namely the degree bound `P=G^3*lambda`.
Thus `B=G^2*lambda` and `P=GB` give precisely this theorem's hypotheses
for every positive width. No new lower-bound proof for piece widths is
required. The conclusion is the same integral `GoodAccuracy` interface
used by `formal/GroupingThreshold.lean`.

Replacing 576 by 144 in the dense note's schedule preserves all its coarse
parameter bounds. In particular
`K_(d+1)=(144 D_d^2 K_d (bit_length K_d+bit_length D_d+8))^2` and
`D_d=(d+48)^2` give `log K_d=Theta(2^d)`. If a uniform evaluator and its
complete machine-time bound are established for this schedule with
contraction `rho`, the conditional rate is `alpha=log(1/rho)/log 2`.
A contraction recurrence alone does not supply that machine theorem.

Importing the existing formal result requires attention to rounding. Its
public `SomeGrouping` statement bounds `Y` by a *ceiling*. The underlying
`stCaseII` theorem proves a strict real inequality before that ceiling.
To preserve this repository's floor-based width, use that stronger result
or explicitly increase the width and charge the additive error. A ceiling
bound alone does not imply the floor-compatible inequality proved here.

The source's case-II divisibility assumption also holds: this schedule
makes `K` divide `g`. Its logarithmic hypothesis follows with ample slack
at `zeta=1/(2D)`. These observations identify a possible formal import path;
they do not establish compilation or interoperability between the two
projects by themselves.

## Sources

- [Bounded-edge grouping statement](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppRunningTimeStatement.lean).
- [Case-II proof](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Group/CaseII.lean).
- [Grouping algorithm and value theorem](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Group/Lemma.lean).
- [Factoring grid width bounds](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Level.lean).
