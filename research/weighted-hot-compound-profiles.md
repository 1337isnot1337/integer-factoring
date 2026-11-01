# Weighted profiles inside the hot compound integrality cone

## Result

The charge-block reduction in
[reflection-natural-cut-duality.md](reflection-natural-cut-duality.md)
also computes complete physical-weighted profiles in an explicit cone
of additional weights. The weights may depend arbitrarily on the
physical state; they need not be sums of one-particle weights. Both
cuts can then be evaluated from selected minors of the same small
one-particle matrix, without constructing the full reflection tensor.

This is an exact profile identity. Numerical optimization inside the
cone still needs a separate rational pressure certificate before it
implies an asymptotic bound.

## 1. A general triangular integrality criterion

Work over a finite field of characteristic two. If rational weights
occur, use one common ramification variable so that all valuations
are integral. Let `C_X,C_Y` be finite deletion-downset codes. The
regular natural-zero matrix has the exact restricted factorization

```math
M=U E_B V,
\quad U(R,R')=[R'\subseteq R]\operatorname{Pf}A[R\setminus R'],
\quad V(S',S)=[S'\subseteq S]\operatorname{Pf}F[S\setminus S'].
                                                               \tag{1}
```

For state weights `g_x,g_y,g_u,g_v`, put
`g_R(x,y)=g_x(x)+g_y(y)` and
`g_S(u,v)=g_u(u)+g_v(v)`. Diagonal matrices `D_R,D_S` have entries
`w^(g_R)` and `w^(g_S)`. Then

```math
D_R M D_S
 =(D_R U D_R^{-1})(D_R E_B D_S)(D_S^{-1}V D_S).        \tag{2}
```

Suppose every nonzero entry in the two outer matrices is regular.
They are triangular with diagonal one, so their inverses are also
regular: the inverse of `I+N` is a finite sum of powers of its strictly
triangular part. Therefore (2) preserves every Smith order.

A sufficient condition can be tested one Gaussian edge at a time.
For each allowed addition of an edge `{i,j}` to a selected row mask
`R'`, require

```math
g_R(R'\cup\{i,j\})-g_R(R')+\operatorname{val}A_{ij}\ge0,       \tag{3}
```

and impose the analogous condition with `g_S,F` on the column code.
Every Pfaffian matching is obtained by adding its disjoint edges one
at a time. Every intermediate mask remains selected because the codes
are downsets. Summing (3) proves regularity of every matching term;
cancellation can only raise its valuation. This proves the criterion.

## 2. Explicit cone for stabilized hot codes

Hot mode `a` has natural height `h_a=4a+2`. Assume all retained states
have energy below `4M`, so the within-leg reflection edges do not
occur. The only within-row edges join `x` mode `a` and `y` mode `b`;
their entries are `(1+w^2)w^(4(a&b))`. Thus (3) becomes

```math
g_x(x\cup\{a\})-g_x(x)
 +g_y(y\cup\{b\})-g_y(y)+4(a\mathbin{\&}b)\ge0,       \tag{4}
```

whenever both additions are allowed. The column condition is (4)
with `(x,y)` replaced by `(u,v)`. These are finitely many rational
linear inequalities. Nonnegative inclusion increments for every leg
are sufficient, but are not necessary. For example, opposite
occupation weights `g_x=c|x|,g_y=-c|y|` satisfy all row inequalities
for every rational `c`.

Outside the stabilized energy range, use the same statement after
adding the within-leg edge constraints. Such an edge joins
`a+b=M/2-1` and has valuation zero in the regular kernel. This makes
the criterion valid for every ambient power-of-two interval.

The weighted compound matrix still splits by total occupation. Within
each block, its row weight is `g_x(x)+g_y(y)` and its column weight is
`g_u(u)+g_v(v)`. No weights or channels are dropped. The endpoint
permutation proves full generic rank separately in each block.

## 3. Both cuts and the full common-type cost

Let the primitive reflection restriction have natural scalar shift
`a=M^2/2-g` as in the companion note. Write `g^0` for additional
zero-place physical weights and `g^infinity` for additional infinity
weights. The complete row-zero profile is

```math
a+\operatorname{Smith}\bigl(D_{g_x^0+g_y^0}E_B
                                D_{g_u^0+g_v^0}\bigr).          \tag{5}
```

At infinity, exchange the `y,u` legs in both the weights and their
cone constraints. The complete column-infinity profile is

```math
-a+\operatorname{Smith}\bigl(D_{g_x^\infty+g_u^\infty}E_B
                          D_{g_y^\infty+g_v^\infty}\bigr).     \tag{6}
```

The two scalar shifts cancel exactly in the sum of the logarithms
of the complete channel sums. The scalar gcd and even the primitive
tensor's actual degree are therefore unnecessary for this two-place
cost computation. This cancellation does not authorize ignoring the
degree when any other place or a different tensor is used.

For common paired types `p_X,p_Y`, the total physical mean is the
natural mean `sum_x p_X(x) E_x+sum_y p_Y(y) E_y`, plus the mean of all eight
additional leg weights. In particular, all changes to physical means
must be retained even though the triangular factors disappear. With
coefficient scale `s>0`, multiply polynomial valuations and natural
energies by `s`; additional weights may be supplied already scaled.
If compound coefficients use `t=w^2`, evaluate their local profile at
scale `2s`. This is the only variable-conversion factor.

## 4. Exact controls and the separate compound crossing

```sh
python tools/check_weighted_hot_compounds.py \
  --output evidence/weighted-hot-compound-profiles.json
```

The [checker](../tools/check_weighted_hot_compounds.py) compares complete
profiles of the original restricted Pfaffian tensor with independently
assembled compound blocks. Its 80 profiles contain 1,240 channels.
The controls include unequal paired alphabets, non-additive state
weights, negative occupation weights allowed by (4), and an interval
where within-leg reflection edges are present. All rational weights
have denominator dividing four, so every compared order is exactly
represented by the native arithmetic interface.

Forty additional full-cost comparisons retain the original physical
means and verify the variable conversion and infinity leg exchange.
An explicit weight outside the cone changes the true Smith profile;
the checker rejects that shortcut instead of treating the cone as an
optional numerical condition. The
[evidence](../evidence/weighted-hot-compound-profiles.json) stores every
complete profile and the countercontrol.

The companion constructor
[hot_compound_crossing.py](../tools/hot_compound_crossing.py) assembles
`E_B` itself as a strict polynomial crossing in variable `t=w^2`.
Every occupation block has the crossing endpoint; its vacuum entry is
one, so its scalar content is exactly one. The checker compares all
548 coefficients of small equal and unequal cases against independent
permutation determinants, and checks every endpoint. This constructor
does not inherit the source tensor's opposite-cut cost.
