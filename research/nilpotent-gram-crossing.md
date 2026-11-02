# Fixed nilpotent Gram crossings and a row-profile obstruction

## Question, outcome and stopping rule

Can noncommuting nilpotent flags produce a degree-12 seven-letter
crossing with the deep flags needed by the
[asymmetric profile target](seven-letter-asymmetric-profile.md)?
This changes the commuting-flag assumption; it is not a larger field,
root, degree or parameter screen.

Two exact binary constructions preserve the flip endpoint and rank one
in both cold cuts. The noncommuting construction changes the column
profile, including producing odd orders, but leaves the row profile

```math
M=12,\qquad\mu=(0,4^{15},16^{33}).
```

That row profile alone excludes the bare two-place target over every
finite coefficient extension. Even the best relaxed opposite profile has
pressure greater than 345, whereas `7^(2 sqrt(2))<7^3=343`.
**These two fixed constructions are stopped.** No flag or field sweep
is justified by their profiles. General nilpotent Gram paths, different
row profiles, physical weights and collective criteria are outside this
statement. The complete factoring exponent remains `1/2`.

The single-flag control is a fourth-power reparameterization of a
cubic path, so its bare exclusion was already covered by the
[integer shared-minor degree floor](integer-shared-minor-pressure-floor.md),
which excludes seven-letter degrees one through eleven. The cycle
construction has odd column orders and actual degree twelve. The
row-specific argument here excludes that case without extrapolating
the earlier degree range.

## 1. An explicit operator-space construction

Work over `F_2`. On a 49-dimensional space, let `N` have nilpotent Jordan
blocks of lengths `(13,12,12,12)`, with `N e_i=e_(i-1)` on each chain.
Let `B_0` be the direct sum of antiidentity matrices on those blocks.
Then `N^T B_0=B_0 N`. The main chain is numbered `0,...,12`; its center
is six. Put `v=B_0 e_0=e_12`.

There are two cases. The control uses `S=I`. The second uses the following
four-chain cycle. If chain `j` starts at `a_j` and has length `l_j`, set

```math
S=I+\sum_{j=0}^3
 e_{a_j+l_j-6}\,e_{a_{j+1}+l_{j+1}-1}^{T},
```

with indices modulo four. This gives `S^2=I`, `S N^5=N^5` and `S^T v=v`,
but `NS!=SN`. Define

```math
O(z)=(zI+N)S(zI+N)^5S^{-1},\qquad
C=O(1)^{-1},\qquad G(z)=C^T O(z)^T B_0 O(z)C.             \tag{1}
```

Every factor of `O(1)` is invertible. The degree of `G` is twelve and
`G(1)=B_0`.

Identify this bilinear space with `Mat_7(F_2)` and `Tr(AB)`. The
antiidentity blocks have 24 hyperbolic pairs and one anisotropic center.
Send the first pair to `(E_01,E_10)`, the center to `I_7`, the next 20
pairs to the remaining off-diagonal matrix-unit pairs, and the final
three to a symplectic basis of the traceless diagonal subspace. The
checker constructs and verifies this isometry `P` explicitly.

In physical matrix-unit coordinates define the row-flattened tensor by

```math
K(z)=P^{-T}G(z)P^{-1}.                                  \tag{2}
```

Thus `K(1;x,y,u,v)=[x=v][y=u]`. This coupled 49-dimensional isometry is
part of the operator constructor. The opposite physical cut is recomputed
after it; no rank invariance under a coupled change of physical letters
is assumed.

The cold tensor is also explicit. Since `S N^5=N^5`, `O(0)=N^6S^-1`.
Further, `(N^6)^T B_0 N^6=B_0N^12=vv^T`, and `N^Tv=0`.
The identities above give `C^TS^-Tv=v`, so `G(0)=vv^T`.
Equation (2) therefore has exactly one nonzero entry, at row `(1,0)` and
column `(1,0)`. Both physical cuts have rank one there.

This is an arbitrary operator-basis path. It is not the triangular
quantum-Pascal ladder vacuum family: its operator determinant is `z^294`,
rather than a unit determinant. That distinction does not supply a gain;
the actual row profile below closes this fixed test separately.

## 2. Exact profiles and an independent row certificate

The [replay](../tools/check_nilpotent_gram_crossing.py) constructs both
tensors and checks 4802 endpoint entries. In both cases the row orders
are `(0,4^15,16^33)` at zero and all zero at infinity. The complete
column profiles are in the [receipt](../evidence/nilpotent-gram-crossing.json).
For the cycle, several cold column orders are odd; its infinity profile
remains `(0^45,8^4)`. These are profiles of the two fixed tensors, not a
classification of nilpotent paths.

The row orders are independently certified without the native Smith
engine by block-Toeplitz nullities:

| Precision | 1 | 4 | 5 | 16 | 17 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Nullity | 48 | 192 | 225 | 588 | 588 |

For full generic rank, nullity at precision `h` is `sum min(h,mu_i)`.
The first two values force one unit and all 48 remaining orders at least
four. The next value makes exactly 15 orders equal four. The last two
values bound every order by 16 and give total 588, forcing the other 33
to equal 16. Generic rank follows from the endpoint, and the leading row
matrix is checked invertible.

## 3. A row-only exclusion over all coefficient extensions

The following argument applies to any degree-12 crossing with this row
profile. Its nonzero row determinant is a degree-588 homogeneous section,
and its full order 588 is spent at zero. There are no other row determinant
zeros, including after coefficient extension. The tensor is nonzero at
zero and row-invertible everywhere else. Thus the other cut has at least
one unit at every place. Its local orders `nu_i` have sum at most 588,
by its own nonzero determinant section.

For `s>=1`, convexity and integer balancing therefore give

```math
\sum_i s^{-nu_i}>=1+36s^{-12}+12s^{-13}.                  \tag{3}
```

To see this, retain one zero order, enlarge the remaining total to 588,
and balance the other 48 integers. Their best values are 36 twelves and
12 thirteens. Enlarging orders decreases the sum, and balancing minimizes
it.

If the row place is zero, the pressure consequently satisfies

```math
B(s)>=(s^{12}+15s^8+33s^{-4})(1+36s^{-12}+12s^{-13})
    >=15s^8+573s^{-4}+36.
```

Put `u=s^4`. Integer AM-GM yields

```math
B(s)>=15u^2+285u^{-1}+285u^{-1}+36
 >=3(15\cdot285^2)^{1/3}+36>345>7^3>7^{2\sqrt2},         \tag{4}
```

where `15*285^2=1218375>103^3=1092727` and `8<9` certify the strict
comparisons. If the row place is elsewhere, all its orders are zero;
Jensen and `sum nu_i<=588` instead give `B(s)>=49^2=2401`.
These alternatives cover both cut orientations and all pairs of places
over every finite coefficient extension. No numerical minimization or
extra field screen is used for the exclusion.

## Reproduction and review scope

```sh
make verify-nilpotent-gram-crossing
```

The exact constructor, endpoint, cold ranks, native profiles, independent
row Toeplitz ranks and integer pressure certificate are replayed. The
opposite profiles use the existing native Smith engine; the row-only
exclusion has the separate Toeplitz control. These fixed finite checks
and the displayed ordinary proof do not constitute independent
mathematical review, a general family obstruction or a new Lean proof.

The [clean full default replay](../evidence/nilpotent-gram-crossing-full-verification.json)
passes at exact source `8b1a041d5400353bd08ac81329fd5de75b741a8e`
(tree `2707e8cfbe84986d8a02b2fd78fc0a2fd7623785`), in 1259.84 seconds.
The archive follow-up adds review prose, earlier-degree-floor credit and
verification records; checker code and the Makefile are unchanged.
