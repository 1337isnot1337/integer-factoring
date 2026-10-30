# Exact flow certificates for channel-walk bounds

## Result and use

The [channel-walk transfer](twisted-channel-walk-transfer.md) bounds an
amplified rank sum using positive piece weights and two nonnegative
transition matrices. A proposed upper certificate consists of positive
rational vectors and numbers satisfying

```math
T_R(z)h_R\le\beta_Rh_R,\qquad
T_C(z^{-1})h_C\le\beta_Ch_C.
```

Its exponential base is `L=sqrt(beta_R beta_C)`. This note gives an
independent **lower certificate for every upper bound of this form**.
The certificate consists entirely of nonnegative integer edge counts.
It can prove that further optimization of the piece weights or positive
vectors cannot improve that particular spectral bound.

This is a bound on the certificate method. It is not a lower bound on
the actual tensor ranks, on other ways to count matching piece types,
or on integer factoring. The distinction is essential: Section 4 gives
graphs for which the optimized spectral bound is strictly larger than
the actual common-type count.

## 1. Balanced integer edge counts

Allow a finite directed multigraph whose edges carry piece labels
`i=1,...,M`. An edge `e:s->t` has weight `z_i`. Thus

```math
(T(z)h)_s=\sum_{e:s\to t}z_{i(e)}h_t.
```

The original channel graph is the special case in which the edge label
is the piece of its destination channel. A scalarized product-space
automaton can use the same formulation after removing states that are
unreachable or cannot lead to a nonzero terminal factor.

Supply an integer `a_e>=0` for every edge. At each vertex require

```math
a_s:=\sum_{e:s\to t}a_e
    =\sum_{e:t\to s}a_e.
\tag{1}
```

The positive total is `D=sum_e a_e=sum_s a_s`. Define the piece
totals and the rational number

```math
n_i=\sum_{e:i(e)=i}a_e,\qquad
K(a)=\frac{\prod_{s:a_s>0}a_s^{a_s}}
           {\prod_{e:a_e>0}a_e^{a_e}}.
\tag{2}
```

Zero counts make no contribution to either product. The data describe
a stationary collection of closed walks, possibly in several strongly
connected components. No probabilistic or numerical validation is needed.

### One-graph inequality

For every positive real `z,h,beta` with `T(z)h<=beta h`,

```math
\boxed{\beta^D\ge K(a)\prod_i z_i^{n_i}.}
\tag{3}
```

To prove it, apply weighted arithmetic-geometric mean to each vertex
with `a_s>0`, using probabilities `a_e/a_s` on its positive-count
outgoing edges. The sum over these edges is no larger than the full
transition sum. Therefore

```math
\beta^{a_s}
\ge a_s^{a_s}
  \prod_{e:s\to t,\ a_e>0}
       \left(\frac{z_{i(e)}h_t}{a_eh_s}\right)^{a_e}.
```

Multiply over vertices. All powers of every `h_s` cancel by (1).
The remaining factors are exactly (3).

## 2. Couple the row and column certificates

Choose flows `a^R,a^C` for the two graphs. Require equal piece
totals:

```math
n_i^R=n_i^C\quad\text{for every }i.
\tag{4}
```

Their total `D` is then equal as well. Apply (3) with weights `z`
in the row graph and `z^{-1}` in the column graph. The piece-weight
powers cancel, proving

```math
\boxed{(\beta_R\beta_C)^D
       \ge K(a^R)K(a^C),\qquad
       L^{2D}\ge K(a^R)K(a^C).}
\tag{5}
```

This holds for all positive real choices, including irrational weights
and arbitrary positive supervectors. In particular it holds for the
infimum over every rational certificate used by the transfer. A rational
flow can be put in this form by clearing denominators in both graphs
with one common multiplier.

The proof also gives the usual entropy interpretation without requiring
it for validation. Put `f_e=a_e/D` and `pi_s=a_s/D`. Then

```math
\frac{\log K(a)}D
 =-\sum_{e:s\to t}f_e\log\frac{f_e}{\pi_s}.
```

The right side is the conditional entropy of the edge choice at a
stationary vertex. Matching piece totals cancel the reciprocal weights.
All comparisons in the supplied checker instead use integer powers and
rational arithmetic.

### Practical search formulation

For upper bounds, writing `z_i=exp(y_i)` and `h_s=exp(u_s)` turns
the supervector constraints into finite log-sum-exp inequalities

```math
\sum_{e:s\to t}\exp(y_{i(e)}+u_t-u_s-b_R)\le1
```

and the corresponding column inequalities with `-y_i`. Their objective
is to minimize `b_R+b_C`. Numerical optimization can propose positive
rational data, which must then pass the exact inequalities. Integer
flows provide an independently checked lower bound. Neither a numerical
minimum nor an unverified eigenvalue is used as a proof.

## 3. Two exact optima

### The square-zero matrix example has optimal certified base eight

The example in the channel-walk note has twelve vertices: four channels
for each of three pieces. Every edge to a different piece is allowed,
and no edge to the same piece is allowed. Give each of the 96 allowed
edges count one. Every vertex has incoming and outgoing count eight.
Hence

```math
D=96,\qquad n_1=n_2=n_3=32,\qquad K(a)=8^{96}.
```

Use the same flow in both cuts. Equation (5) forces `L>=8`.
Constant weights and positive vectors already certify `L=8`, so eight
is the exact optimum of the weighted-supervector method for this graph
pair. Adjusting piece weights cannot improve it. This example remains
worse than the direct scalar crossing, as the original note states.

### Complete graphs recover the rank-sum optimum four

Take row channel counts `(1,4)` and column channel counts `(4,1)`.
The positive weights `(2,1/2)` and constant vectors certify
`beta_R=beta_C=4`.

On the row graph assign vertex probabilities `1/2` to the first piece's
single channel and `1/8` to each of the other four channels. Reverse the
roles in the column graph. For every directed edge, including loops,
set `a_(s,t)=64 pi_s pi_t`. These are integers, are stationary, and
give total 32 to each piece in both graphs. Directly from (2), both
graphs have `K=4^64`. Thus (5) forces `L>=4`, matching the upper
certificate exactly.

## 4. Scope: the spectral certificate can lose type information

Consider a row graph with two isolated loops, labelled `A` and `B`.
Let the column graph be a two-cycle whose edge labels alternate `A,B`.
In each graph give the two edges count one. Their piece totals agree,
`D=2`, and both values of `K` are one. Constant weights give spectral
base one, and (5) proves that this is the exact optimized certificate.

Nevertheless, for a walk with at least two edges, the row graph produces
only the types consisting entirely of `A` or entirely of `B`; the
column graph produces neither. The exact common-type sum
`sum_nu sqrt(R_nu C_nu)` is zero. The stationary row flow mixes two
components that a single walk cannot visit.

These abstract graphs cannot certify a nonzero crossing for every
length. They show why (5) must not be advertised as a lower bound on
actual common-type counts or tensor ranks. More detailed component or
type analysis can improve on the optimized spectral upper bound.

## Verification

```sh
python3 tools/check_channel_flow_duals.py
```

The standard-library checker validates both exact optimum certificates,
tests the inequalities with independent positive rational weights and
vectors, and rejects nonstationary flows, mismatched piece totals, and
counts assigned to missing edges. It also verifies the type-information
counterexample by explicit walk enumeration. Evidence is recorded in
`evidence/channel-flow-duals.json`. No factoring-score improvement is
claimed by these graph optimization results.
