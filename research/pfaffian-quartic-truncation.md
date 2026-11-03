# A quartic Gaussian truncation and its stopped q8 completion

## Result and scope

There is a characteristic-two creation identity: multiplying a Gaussian
by all quartic gates whose coefficients are its four-variable Pfaffians
leaves just its vacuum and pair terms. The identity survives the one-hot
quotient. It supplies a structural way to remove the previously protected
two-particle diagonal channels, by allowing repeated colours.

One fixed homogeneous completion applies this identity to both cold cut
jets of the q8 curve while restoring the crossing endpoint. It succeeds
in removing every order-twelve through order-fifteen channel: through
fifteen the full profile is only `0^1,6^14`. However its middle layer
returns **322 order-sixteen channels**. A dimension/degree tail and exact
positive polynomial prove a common-scale, unweighted paired-type profile
floor **4/5>1/sqrt(2)** for this completion. It is stopped without a full
profile or larger cutoff/alphabet/field/coefficient search.

The complete proven factoring exponent remains **1/2**. This result does
not exclude all repeated-colour interactions, other completions, degree
costs, endpoints, kernels, weighted gauges or actual coefficient-rank
methods. The algebraic identity is not an algorithm for evaluating a
growing Gaussian state for free.

The [mechanism](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19365241)
and [covariant completion](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19365303)
were specified before the tensor experiment. The
[checkpoint](https://github.com/1337isnot1337/integer-factoring/discussions/51#discussioncomment-19365832)
records the successful low-layer prediction and the fixed completion's
obstruction. The three-monomial leading template was not separately
tested or optimized.

## 1. Quartic Pfaffians keep only the pair state

Work in a commutative characteristic-two coefficient ring with exterior
generators xi_i. They commute and square to zero. Let kappa be any
two-form, and let G_r be its divided-power Gaussian component of
degree 2r. Thus G_0=1, G_1=kappa, and each coefficient of G_r is the
corresponding principal Pfaffian. Write `G=sum_r G_r` and

```math
Q=\prod_{|S|=4}\bigl(1+\operatorname{Pf}(K_S)\xi_S\bigr).
```

For any matching with 2r edges, grouping the edges into r unordered
pairs can be done `(2r-1)!!` ways, an odd number. Expanding the quartic
product therefore gives `Q=sum_r G_(2r)`. Also
`kappa G_(r-1)=r G_r`: distinguish one of the r pairs of a matching.
Consequently

```math
Q=G(1+\kappa),\qquad G^2=1,\qquad
\boxed{GQ=1+\kappa}.                               \tag{1}
```

The middle equality follows either by squaring the product of nilpotent
pair gates or by characteristic two. All identities are finite and use
no factorial division. Imposing the additional one-hot relations
`xi_(v,c) xi_(v,d)=0` for generators at the same physical vertex is a
ring quotient and preserves (1). Terms with repeated colours at
different vertices remain allowed. No novelty claim is made for these
matching/divided-power identities.

## 2. One fixed polynomial tensor with both cut jets

Retain the Gaussian pair kernel, fixed parent-leg pairing and
`(X,Y,Z)=(x^2,y^2,xy)` curve of
[PR111](fixed-curve-quartic-floor.md). Replace the earlier distinct-colour
quartic interaction, rather than adding another interaction to it.

For an assigned four-vertex word S, take its cross-Gaussian four-occupied
coefficient at the x cut. It is the sum of the two cross-matching products
with colour-equality deltas; it is zero on a one/three split. Copy its
x jet through fifteen at y=1, `A_S(t)=sum_(j=0)^15 a_j t^j`.
Do the same at the fixed y cut to obtain `B_S(t)=sum_(j=0)^15 b_j t^j`.
Both jets start at twelve. Set

```math
C_S(x,y)=\sum_{j=0}^{15}a_jx^jy^{32-j}
 +\sum_{j=0}^{15}b_jx^{32-j}y^j
 +\epsilon_Sx^{16}y^{16},\qquad
\epsilon_S=\sum_j(a_j+b_j).                         \tag{2}
```

These are disjoint coefficient ranges of one homogeneous polynomial,
not a cut chosen for an individual tensor summand. It has degree 32,
valuation at least twelve at both axes, and C_S(1,1)=0. Use every gate
`1+C_S xi_S`, including the all-equal and two-pair colour patterns.
Four distinct colours give zero. The whole tensor retains:

- the original invertible antipodal endpoint and fixed internal z swap;
- cold rank one at both axes;
- homogeneous degree eight per occupied old vertex, including products
  of disjoint quartics;
- the q^8 complete local table and maximum entry degree 64.

Eight-vertex support expansion constructs this table in q^O(1) time,
including its output; it is not an uncharged growing object. At the
fixed q8 size the calculation uses binary polynomial arithmetic.
The word-degree invoice remains `32m E_nu`, and the full determinant
degree remains `56*4096`.

Remove the one-sided creation gates separately for each cut. They are
polynomially invertible and commute with the remaining even creators.
Now (2) copies the entire cross-Gaussian two-particle jet through fifteen,
so adding it cancels that whole block. Quartics on one/three splits
start at sixteen. The one-particle kernel and its orders
`(6,6,18,18)` are unchanged as a submatrix. Thus the full low profile
through fifteen is only `0^1,6^14`. This changes the premise of the
[first-fifty-order barrier](q8-leading-prefix-floor.md); it is not an
attempt to optimize later layers while retaining its thirty-five
order-twelve channels.

At six occupied vertices, the leading order-eighteen terms cancel by
`kappa G_2=3G_3=G_3`. This is a statement about those raw coefficients,
not the absence of mixed channels in the Smith profile. In particular,
the later middle layer can couple one- and three-occupied states.

## 3. Exact generic low jets and independent support expansion

The [constructor](../tools/pfaffian_quartic_truncation.py) and
[checker](../tools/check_pfaffian_quartic_truncation.py) retain all 512
rows and columns in charge zero and one at both cuts, including repeats
and all occupations. The coefficients depend on colours only through
equality, so they are invariant under every permutation of the seven
nonzero colours. In particular GL(3,2) preserves XOR charge and acts
transitively on nonzero charges. These two representative charge blocks
therefore determine the full q8 profiles.

At both cuts the exact orders through eighteen are

```math
\begin{array}{c|c}
\text{charge zero}&0^1,16^{42}\\
\text{each nonzero charge}&6^2,16^{40}\\
\text{whole matrix}&0^1,6^{14},16^{322}.
\end{array}                                        \tag{3}
```

Independent column elimination checks truncations sixteen and nineteen.
There are 115,616 direct six-occupied leading-coefficient checks. The
generic valuation is preserved by dehomogenization, since the physical
occupation homogeneity makes the other variable's row/column factors
units, as in the fixed-curve source proof.

The independent literal calculation uses recursive Pfaffians and full
support expansion, including **two** quartic gates at eight occupied
vertices. It checks 65,536 entries per cut, endpoints and degree, then
compares their order-eighteen truncations with the fast matching formula.
The panel uses vacuum and three nonzero colours. Every nonzero jet
through eighteen has at most three distinct colours, so colour symmetry
covers the other colour patterns. Two-quartic contributions are nonzero
in 609 panel entries per cut; the saved receipt includes a witness.
They start at twenty-four and cannot affect (3), but are not discarded
from the literal expansion. The old distinct-colour q8 divided-square
vanishing identity is not assumed here.

The [receipt](../evidence/pfaffian-quartic-truncation.json) also retains
192 controls for (1), 9310 assigned repeated-colour quartets, exact
Toeplitz ranks, literal hashes, the two-quartic witness and source hashes.
The controls support the ordinary identity proof; they do not replace it
by a finite extrapolation.

## 4. The restored profile floor is already 4/5

The 337 known orders in (3) sum to 5236. The remaining 3759 orders have
mean at most `(56*4096-5236)/3759=32020/537<60`. Jensen and monotonicity
give, independently for both partitions,

```math
Z_i(t)\ge1+14t^6+322t^{16}+3759t^{60}=:L(t),\quad0<t<1.
```

With u=t^2, the polynomial

```math
\frac{(1+14u^3+322u^8+3759u^{30})^5-(1+7u^5)^{16}}{u^3}
```

has degree 147 and **513 strictly positive degree-512 Bernstein
coefficients**. The checker reconstructs every monomial, including
zeros above degree 147. Its minimum is at index 21 and its integer
coefficient-list hash is
`f4473321ae0fdfa2b195cb2d285f5be908cb2978deea23d6bf6493d6b95ba069`.
Therefore `Z_i(t)>(1+7t^10)^(16/5)` on the full interval.

For any correlated paired physical type with positive entropy,
`H(nu)<=2log(1+7t^10)-10E_nu log t`. Combine this with the two partition
bounds and the unchanged word degree:

```math
\frac{\log Z_0(t)+\log Z_1(t)-32E_\nu\log t}{4H(\nu)}
 >\frac45>\frac1{\sqrt2}.
```

The coefficient `(16/5)*10=32` accounts for the full degree. This is
an ordinary profile-envelope exclusion for this fixed completion, not a
factoring or actual coefficient-rank lower bound. The approximately
0.83888 lower-envelope fit was a private diagnostic; no exact minimum or
new factoring score is inferred from it.

## Verification, credit and stop

```sh
make verify-pfaffian-quartic-truncation
```

Focused controls pass; full default verification is pending at this
source checkpoint. No new Lean or independent mathematical review is
claimed. The Gaussian/colour/PR97 kernel, Sol worker 2's PR111 rank and
curve scaffolding, and Sol worker 1's PR114 jet/creation checks retain
their credit. The repeated-colour Pfaffian completion and its scoped
analysis were developed by Sol worker 2 with substantial OpenAI Codex
assistance.

This completion is stopped. No larger profile, cutoff, alphabet, field
or coefficient tuning follows from it. The identity (1) remains a useful
algebraic tool; using it without this endpoint-restoration cost needs a
different construction or arithmetic interface.
