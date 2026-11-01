# A simultaneous time and workspace transfer for crossing families

## Results and scope

The crossing interface can preserve the same logarithmic saving in
both the running-time exponent and the working-space bound. This is a
conditional theorem about a supplied crossing family. It does not
construct a family with contraction below `1/sqrt(2)`.

Fix a prime `p`, a contraction `1/2<=rho<1`, and `delta>0`. Set

```math
\alpha=\log_2(1/\rho),\qquad
\beta=\min(\alpha,\delta),\qquad
F_\beta(n)=\frac{n}{(\log(n+2))^\beta}.
```

Suppose an exact scale-one crossing on every binary label length
`lambda` has the following properties:

1. Both displayed factorizations have individual binary widths at
   most `2lambda` and total width at most

   ```math
   4\rho\lambda+
   \frac{C\lambda}{\log_2\log_2(\lambda+4)}.
   ```

2. Its number `M_lambda` of pieces satisfies

   ```math
   f(\lambda):=\left\lceil\log_2\max(2,M_\lambda)\right\rceil
   \le\frac{C\lambda}{(\log_2(2\lambda))^\delta}.
   ```

3. Complete piece enumeration, widths, and scalar factor queries
   have uniform algorithms taking `2^O(lambda)` bit time and
   `O(lambda)` working bits. A list need not be retained. The
   enumerator emits one piece at a time and has no duplicates.

Then a deterministic finite-alphabet machine computes the least prime
factor, complete prime factorization, sorted divisors, and exact
divisor-prefix counts of an `n`-bit integer using

```math
T(n)=\exp(O(F_\beta(n))),\qquad S(n)=O(F_\beta(n)).       \tag{1}
```

The input is read-only and the output is write-only. Prefix data are
also read-only and have the natural `O(n)`-bit encoding. In particular,
when `delta>=alpha`, both bounds retain the `rho`-derived exponent
`alpha`, including `alpha=1`.

There is a second, constructive interface improvement. If the original
family has polynomially many pieces, width error `O(lambda^theta)`
for some fixed `0<=theta<1`, and `2^O(lambda)` local time, then no
separate space hypothesis is needed. Type grouping of logarithmic-size
blocks produces a family satisfying all three assumptions, with the
same `rho` and with any fixed `delta`. Section 4 proves this, including
the endpoint `rho=1/2`.

All statements are ordinary algorithmic arguments. They use the
source grid recursion and the existing crossing identities; they do
not assert new Lean formalization or an improved unconditional time
exponent.

## 1. Packed piece indices and the grid accounting

Use the schedule, larger top constructor, and `n^(3/4)` early cutoff
from [sublinear-space-half-rate.md](sublinear-space-half-rate.md).
The same proof works for every fixed `rho` in the stated interval.
Writing `P=K_ell`, each grid step has

```math
\lambda_j\ge n^{1/4}/8,\qquad
V_j=K_j^2\lambda_j\le D\rho^{\ell-j}V_\ell,
\qquad V_\ell=O(n),\qquad
\rho^\ell=\Theta((\log n)^{-\alpha}).                  \tag{2}
```

The product of the width-slack factors is bounded: there are
`O(loglog n)` steps and each additional relative error is
`O(1/loglog n)`, besides the summable schedule errors. This is exactly
the geometric-envelope argument in the
[relaxed piece transfer](subexponential-piece-transfer.md).

At each crossing site retain the **ordinal** of its piece in the
enumeration, using `f(lambda)` bits. If a supplied piece descriptor is
longer, recover it by rerunning the enumerator up to that ordinal.
This costs `2^O(lambda)` time and `O(lambda)` scratch space by
assumption 3. It avoids silently identifying a cardinality bound
with a compact supplied representation.

If the emitted descriptor is itself longer than `O(lambda)` bits,
do not retain it. Supply each requested descriptor bit by rerunning
the enumerator and counting its output positions. Its output length
is at most `2^O(lambda)`, so a position needs only `O(lambda)` bits.
Simulate the factor-query algorithm with this virtual read-only
descriptor tape. Multiplying its `2^O(lambda)` steps by the
`2^O(lambda)` cost per requested bit still gives `2^O(lambda)`
time, with only a constant number of `O(lambda)` workspaces live.

The number of pieces can itself be counted by one enumeration, using
an `O(lambda)`-bit counter. Therefore the ordinal width and all packed
offsets are available uniformly. Store at most the piece ordinal and
two `O(log(lambda+2))`-bit widths at each site. Equations (2) imply

```math
\begin{aligned}
\sum_j K_j^2 f(\lambda_j)
  &\le C'\sum_j\frac{V_j}{(\log n)^\delta}
   =O\!\left(\frac n{(\log n)^\delta}\right),\\
\sum_j K_j^2\log(\lambda_j+2)
  &=O(n^{3/4}\log(n+2)).                              \tag{3}
\end{aligned}
```

The same first sum charges the logarithm of the piece-assignment
branching count. Thus compact storage does not create an uncharged
enumeration cost. In particular, the endpoint `delta=alpha` has no
extra `loglog n` factor.

Grouping metadata fit in `O(P^2 log(n+2) loglog(n+2))` bits, or can
be regenerated. Current grouped assignments, upward query arguments,
and local factor scratch have the same bound as in the sublinear
theorem:

```math
\sum_{j=1}^{\ell}\frac{V_j}{K_{j-1}}
\le Dn\rho^\ell
       \sum_{j\ge1}\frac{\rho^{-j}}{K_{j-1}}
=O\!\left(\frac n{(\log n)^\alpha}\right).             \tag{4}
```

The series converges uniformly for `1/2<=rho<1` because
`K_(j-1)>=2^(18j)`. Local `O(lambda_j)` scratch fits (4), since
`lambda_j<=V_j/K_(j-1)`. Only one such scratch area is active.
The level-zero leaf has `O(n rho^ell)` assignment bits; an earlier
leaf has at most `n^(3/4)` bits.

The [streamed top query](streamed-top-space.md) uses
`O(n/P+Af+log n)` bits, apart from its weight-control registers.
Here `Af=O(n^(3/4))` and `P=Omega(n^(1/40))`. These terms, and the
smaller grouping records, are `o(F_beta(n))` for every fixed
`0<beta<=1`.

The resulting weighted evaluator has time

```math
\exp\!\left(O\!\left(
 \frac n{(\log n)^\alpha}+
 \frac n{(\log n)^\delta}\right)\right)
```

and the same sum as its workspace bound. Local regeneration is
already covered by the parent-query term in (4). No full piece list,
factor table, or grouped table is materialized.

## 2. Fixed coefficient rings and paired descriptor recovery

Use the paired factor-input polynomial and the descriptor recovery
from [sublinear-space-half-rate.md](sublinear-space-half-rate.md).
The elementary uniform bound
`tau(N)<=2^O(n/log(n+2))` gives control descriptors of
`O(n/log(n+2))` bits. All recovered candidates remain virtual.

For characteristic two, use the fixed sixteen-element ring
`F_4[epsilon]/(epsilon^2)` and ternary evaluation hashes. For odd
`p`, use `F_p[epsilon]/(epsilon^2)` and binary evaluation hashes;
the character phases are `+1,-1`. In the latter case the projection
is multiplied by the inverse of the odd-characteristic scalar
`Q=2^k`. This scalar is nonzero modulo `p`. The hash field is only
a control representation; it is not the coefficient ring of the
weighted evaluator.

The polynomial root bound gives a singleton bucket for every actual
divisor once `Q>n(tau(N)-1)`. Each coordinate query uses `Q` weighted
evaluations. Hash only the first factor's bits, but place the single
square-zero tag on either input family to recover both factors.
Scalar-one buckets with coordinates in `{0,1}` give candidate pairs.
Verify the pair by the streamed column multiplication comparator.
Its carry is at most `n`, so it adds only logarithmic control storage.

Strict successor scans, comparisons, prime filtering by divisor-pair
products, and prime-power scans retain only a constant number of
descriptors. Their number of evaluations is polynomial in `n` and
`2^O(n/log(n+2))`. Therefore recovery adds

```math
O(n/\log(n+2))
```

to the logarithmic time and workspace bounds. Since `beta<=1`, this
is absorbed by `F_beta(n)`, proving (1). The same fixed-ring algebra
and odd-characteristic projection are established in
[bit-tagged-witness-recovery.md](bit-tagged-witness-recovery.md).

### An independent time consequence of the space bound

After all uniform query routines are inlined, the construction is a
halting deterministic finite-alphabet transducer. Use the standard
space convention that counts visited work-tape cells, with a fixed
number of tapes and heads. The read-only input tapes have `O(n)`
cells and endmarkers; their heads stay within those endmarkers.
The output is one-way write-only: its contents and current position
cannot be inspected by the machine.

An internal configuration consists of the control state, work-tape
contents and head positions, and input-head positions. With `S`
working bits there are at most `poly(n,S)c^S` such configurations
for a fixed constant `c`. Output contents and its head position do
not affect future transitions and can be omitted. A repeated internal
configuration would therefore repeat forever, contradicting termination.
Consequently the proved `S=O(F_beta(n))` bound also implies

```math
T\le\operatorname{poly}(n,S)c^S=\exp(O(F_\beta(n))).
```

Here `F_beta(n)` dominates `log(n+2)`. This supplies an independent
time corollary for the actual uniform finite-tape implementation.
It does not apply to an oracle machine before the oracle routines
and their storage have been included. The explicit operation
accounting above remains a direct verification of the construction.

## 3. What the local-space assumption covers

The enumerator in assumption 3 must itself be a streamed algorithm.
A `2^O(lambda)` time bound alone permits exponentially many work
bits and does not justify (3)--(4). Likewise, a large table of all
factor values is not an `O(lambda)`-space query interface.

Polynomial piece count with an explicitly linear-space local oracle
is a special case of the theorem for every fixed `delta`. The regular
characteristic-two family has such an oracle and gives the proved
half-rate sublinear-space theorem. The next section provides a
general way to obtain the space hypothesis from a stronger width
error but with no original space guarantee.

## 4. Type grouping regularizes logarithmic blocks, including the endpoint

Suppose the original family at label length `t` has at most
`C(t+1)^d` pieces, separate widths at most `2t`, and total width at
most

```math
4\rho t+C t^\theta,\qquad 0\le\theta<1.               \tag{5}
```

Assume all required local generation and query algorithms have bit
time at most `2^(C_0 t)`, increasing `C_0` to cover constants. Their
workspace is consequently at most this time. Choose a fixed integer
`B>2C_0`, and put

```math
L=\lceil\log_2(\lambda+2)\rceil,\qquad
r=\max(1,\min(\lambda,\lfloor L/B\rfloor)),\qquad
\lambda=br+q,\quad 0\le q<r.                           \tag{6}
```

There are `b` equal full blocks and at most one shorter tail. Let
`k=|I_r|=O(r^d)`, and enumerate the full-block pieces as
`A_1,...,A_k`. A new piece consists of a count vector
`nu=(nu_1,...,nu_k)`, with `sum nu_i=b`, and one tail piece when
`q>0`. Define

```math
B_{\nu,a}=
 \left(\sum_{\operatorname{type}(w)=\nu}
             A_{w_1}\otimes\cdots\otimes A_{w_b}\right)
 \otimes A_a^{(q)}.                                    \tag{7}
```

For `q=0`, omit the tail and use one empty factor. Summing (7) over
all types and all tail pieces recovers the full tensor product of the
block crossing sums, hence exactly `E_(2^lambda)`.

### Explicit channels and width

Write `r_i,c_i` for the original padded binary channel widths.
Put `M_nu=b!/prod_i nu_i!`. A row channel chooses a word of type
`nu`, followed by one original row channel at each position and
one tail channel. Thus its exact cardinality is

```math
R_{\nu,a}=M_\nu\,
  2^{\sum_i\nu_i r_i+r_a^{(q)}}.                       \tag{8}
```

Use products of the corresponding original left factors and products
of their right factors. Summing the channels proves (7). The column
factorization has the same construction with `c_i` in place of `r_i`.
Zero-padding (8) to a power of two adds `ceil(log_2 M_nu)` bits in
each cut. It does not add this overhead separately at each block.

The entropy bound `M_nu<=k^b` and (5) give total width

```math
\begin{aligned}
r_{\nu,a}+c_{\nu,a}
&\le4\rho\lambda+
 O\!\left(\lambda r^{\theta-1}
                  +\frac{\lambda\log(r+2)}r+1\right)\\
&=4\rho\lambda+
 O\!\left(
 \frac{\lambda}{(\log\lambda)^{1-\theta}}
 +\frac{\lambda\log\log\lambda}{\log\lambda}\right).
                                                               \tag{9}
\end{aligned}
```

This error is `o(lambda/loglog(lambda+4))`. If an individual padded
width exceeds `2lambda`, use the identity matrix times the complete
piece matrix in that cut. The entry algorithm below makes this cap
effective. Capping cannot increase the total width.

The number of new pieces is at most

```math
{b+k-1\choose k-1}\max(1,|I_q|).
```

Its binary logarithm, and the length of the count-vector descriptor,
are both

```math
O(k\log(b+1)+\log(r+2))
=O((\log(\lambda+2))^{d+1}).                           \tag{10}
```

Thus assumption 2 holds for **every fixed `delta>0`**. Type grouping
has moved the multinomial cost into the small width error (9); it
has not discarded the permutations of the blocks.

### Uniform query algorithms and working space

At a full block the original oracle uses at most
`2^(C_0 r)=O(sqrt(lambda+2))` work bits, by (6). Only one block
oracle runs at a time. Raw piece descriptors can be regenerated from
their ordinal within this workspace.

Compute multinomial counts with integer arithmetic. Even the simple
factorial implementation stores `b!` in `O(b log(b+1))=O(lambda)`
bits. Unrank a type word by testing its possible first letters with
the counts `(b-1)!/prod_i(nu_i-[i=j])!`, then repeat. Decode the
remaining local channel index by the block widths. The local channel
cardinality in (8) is independent of the order of the word, so this
is a genuine bijective channel numbering. A factor query multiplies
the corresponding block values sequentially.

For a complete-piece query, enumerate words in `[k]^b`, retain those
with the required type, and sum the products of their block entries.
A word counter uses `b ceil(log_2 k)=O(lambda)` bits, and the count
vector uses the smaller bound (10). This enumeration takes

```math
k^b\operatorname{poly}(\lambda)2^{O(r)}
 =2^{O(\lambda\log\log(\lambda+2)/\log(\lambda+2))}
 \le2^{O(\lambda)}                                      \tag{11}
```

time. It also supplies capped-factor entries. The count vector,
word, channel, factorial counters, and one active block workspace
together use `O(lambda)` bits. As usual, increasing fixed constants
covers the finitely many small lengths before the displayed
asymptotic comparisons apply.

Enumerate count vectors and tail ordinals sequentially. Equation (10)
shows that complete new-family enumeration still has `2^O(lambda)`
time and `O(lambda)` workspace, including width generation. All three
hypotheses of the transfer are now satisfied. Taking `delta>=alpha`
proves

```math
T(n)=\exp\!\left(O\!\left(\frac n{(\log n)^\alpha}\right)\right),
\qquad
S(n)=O\!\left(\frac n{(\log n)^\alpha}\right)           \tag{12}
```

from the original time-only local interface, for every fixed
`0<alpha<=1`.

### A faster scalar-entry algorithm

The word enumeration in (11) is sufficient, but can be improved.
For a fixed external entry, write its full-block coefficient polynomial
as

```math
H(Z_1,\ldots,Z_k)=\prod_{j=1}^b
       \left(\sum_{i=1}^k a_{j,i}Z_i\right),
\qquad a_{j,i}\in\mathbb F_p.
```

Choose the least extension-field size `Q=p^s` with `Q-1>b`. Every
variable degree is at most `b<Q-1`. The finite-field power-sum identity
therefore gives the exact coefficient formula

```math
[Z^\nu]H=(-1)^k
 \sum_{z\in\mathbb F_Q^k}
 H(z)\prod_{i=1}^k z_i^{Q-1-\nu_i}.                   \tag{13}
```

All exponents in the product are positive. For a monomial exponent
`m_i` between zero and `b`, the relevant field sum is zero unless
`m_i=nu_i`; in that case it is `Q-1=-1` in the field. This proves
(13), including the sign in odd characteristic. The result belongs
to the base field because it is the original coefficient.

Enumerate the field vectors sequentially, retaining `k` field elements
and a constant number of accumulators. Evaluate each factor of `H`
by regenerating its block entries. There are `Q^k` points, and
`Q=O_p(b+2)`. Construct the field by exhaustive irreducible-polynomial
testing, which has polynomial cost in `Q` and `O(log Q)` storage.
Consequently complete-piece and capped-factor queries take

```math
2^{O(k\log(b+2)+\log(\lambda+2))}
=2^{O((\log(\lambda+2))^{d+1})}                       \tag{14}
```

time. Their additional field workspace is
`O(k log(b+2))=polylog(lambda)`, and the overall workspace remains
`O(lambda)`. Uncapped factor queries already have polynomial time.
Thus the regularized family has quasipolynomial local query time,
as well as the required linear local workspace. The slower formula
(11) remains an independent direct verification of effectivity.

The literal tuple blocking in
[blocked-half-family.md](blocked-half-family.md) gives stronger
nearly linear local query time in its stated range. The present
type grouping has a different purpose: it gives a compact global
piece set and removes the endpoint obstruction for workspace. It
uses only the slack already allowed by the relaxed grid transfer.

## Verification and limits

The companion checker tests both type-grouped cuts, channel unranking,
tail blocks, power-of-two padding, identity-factor caps, complete
crossing sums, finite-field coefficient extraction in characteristics
two and three, and exact packed-cost inequalities in finite examples.
Those examples do not supply a better contracting family. The
asymptotic conclusions follow from (2)--(11), not from finite rank
extrapolation.

```sh
python tools/check_simultaneous_time_space_transfer.py \
  --output evidence/simultaneous-time-space-transfer.json
```
