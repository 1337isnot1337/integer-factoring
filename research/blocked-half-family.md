# A half-rate family with nearly linear-time factor queries

## Theorem and scope

The characteristic-two regular-partition crossing can be regularized
into a family on every binary label length `lambda` with the following
properties:

```math
\begin{aligned}
\log_2 |I_\lambda|
  &=O\!\left(\frac{\lambda}{\sqrt{\log(\lambda+2)}}\right),\\
r_i,c_i&\le2\lambda,\\
r_i+c_i
  &\le2\sqrt2\,\lambda+
       O\!\left(\frac{\lambda}
       {(\log(\lambda+2))^{1/8}
        (\log\log(\lambda+4))^{1/4}}\right).
\end{aligned}                                                   \tag{1}
```

Each piece, factor, width, label-index, and next-piece query takes

```math
\lambda\,2^{O(\sqrt{\log(\lambda+2)}
                         \log\log(\lambda+4))}
   =\lambda^{1+o(1)}                                      \tag{2}
```

bit operations and `O(lambda+1)` bits of working storage. Complete
enumeration still takes time proportional to the number of pieces.
Equation (2) is a query bound, not a claim that the entire index set
is generated in nearly linear time.

The construction is an exact scale-one crossing, with explicit factors
in both cuts. The
[relaxed transfer](subexponential-piece-transfer.md) applies at
`rho=1/sqrt(2)` and at its piece-entropy endpoint `delta=1/2`.
Consequently it retains the established half-rate factoring time.
It also fits the packed sublinear-space evaluator in
[sublinear-space-half-rate.md](sublinear-space-half-rate.md).
Neither conclusion improves the time exponent beyond one half.

The improvement is in the local interface. The original large
partition family uses `2^O(lambda)` local query time; this family has
the stronger bound (2), at the cost of a subexponential piece set.
It does not require a finite seed that individually beats the
half-rate coefficient.

## 1. The exact finite-block input

For every positive integer `t`, take the original regular `F_2`
crossing on the first `2^t` distinct-odd partitions of the least
weight `h_t` with enough labels. Use the original translated duals.
The established bounds are

```math
h_t=O(t^2),\quad |I_t|=4h_t+1,\quad
r_a,c_a\le2t,\quad
r_a+c_a\le2\sqrt2\,t+C t^{3/4}.                         \tag{3}
```

The constant `C` can be increased to cover all positive `t`.
The separate caps in (3) use the trivial identity factorization on
the binary label prefix when a raw channel count exceeds `2^(2t)`.
The individual channels are padded with zeros to powers of two.

The exact crossing and its two factorizations come from
[modular-crossing.md](modular-crossing.md). The quantitative
`O(t^(3/4))` width error is proved in
[rough-input-transfer.md](rough-input-transfer.md), Section 6.
The [linear-space factor theorem](near-linear-factor-space.md)
supplies uniform algorithms for all factors, capped piece entries,
widths, and label/channel unranking in `2^O(t)` bit time and `O(t)`
bits. No finite table of partitions or factors is assumed.

## 2. Choose the block lengths with integer arithmetic

For `lambda>=1`, define

```math
L=\left\lceil\log_2(\lambda+2)\right\rceil,\qquad
u=\lceil\sqrt L\rceil,\qquad
v=\left\lceil\log_2(L+2)\right\rceil,\qquad
r=\min(\lambda,uv).                                    \tag{4}
```

Use blocks of `r` bits and one final shorter nonempty block, if needed.
Write their lengths as `t_1,...,t_b`. Then

```math
\sum_jt_j=\lambda,\qquad t_j\le r,\qquad
b=\lceil\lambda/r\rceil,\qquad
r=\Theta(\sqrt L\log(L+2))
```

for all sufficiently large `lambda`. All rounding and small lengths
are handled by (4) itself. The square root and logarithmic ceilings
are obtained by exact integer comparisons or bit lengths.

Split each external label into these contiguous bit blocks. For a
tuple `a=(a_1,...,a_b)` of finite-block pieces, define

```math
B_a(x,y,u,v)=\prod_{j=1}^b
   A_{a_j}^{(t_j)}(x_j,y_j,u_j,v_j).                    \tag{5}
```

Finite distributivity and each block's exact crossing sum give

```math
\sum_{a\in I_{t_1}\times\cdots\times I_{t_b}}B_a
  =\prod_j E_{2^{t_j}}=E_{2^\lambda}.                   \tag{6}
```

This uses the full tuple index set. No coefficient regrouping,
cancellation assumption, or guessed rank formula is involved.

Tensor the displayed row factors of the blocks, and tensor their
column factors separately. Their binary widths are exactly the
sums of the block widths. In particular both are at most
`2 sum_j t_j=2lambda`. A channel index is the concatenation of the
block channel bit strings, including their zero padding. There is
no global mixed-radix division or subsequent global cap.

## 3. Width and piece-entropy estimates

Equation (3) gives

```math
r_a+c_a\le2\sqrt2\,\lambda+
 C\sum_jt_j^{3/4}
\le2\sqrt2\,\lambda+2C\lambda r^{-1/4}.                 \tag{7}
```

To see the last inequality, all but the last block have length `r`,
and `b r^(3/4)<=(lambda/r+1)r^(3/4)<=2lambda r^(-1/4)`.
The asymptotic size of `r` gives the width error in (1). Also

```math
\frac{\log\log(\lambda+4)}{r^{1/4}}\longrightarrow0,
```

so this error is smaller than the allowed
`O(lambda/loglog(lambda+4))` slack.

Each block has `O((t_j+1)^2)` pieces. Therefore

```math
\begin{aligned}
\log_2|I_\lambda|
 &\le C_1 b\log_2(r+2)\\
 &=O\!\left(\frac{\lambda\log(r+2)}r\right)
  =O\!\left(\frac{\lambda}{\sqrt{\log(\lambda+2)}}\right).
\end{aligned}                                                   \tag{8}
```

Here `log(r+2)=O(log(L+2))`, and the `log(L+2)` factor is
already present in the denominator `r`. This cancellation is the
reason for the extra logarithm in (4).

Thus the exact hypotheses of the relaxed transfer hold with
`rho=1/sqrt(2)` and `delta=1/2`. Its geometric volume sum explicitly
permits `delta=alpha` without an additional depth factor.

## 4. Nearly linear-time queries

Encode a piece tuple by its block shifts. For a block of length `t`,
the shift belongs to `[-2h_t,2h_t]` and takes `O(log(t+2))` bits.
The total index length is bounded by (8). Invalid padded shift codes
can be rejected. An odometer generates every valid tuple, retaining
only its packed index and scalar counters.

For a factor query, scan the piece index and the external label and
channel strings sequentially. At block `j`:

1. Read its shift and compute its two finite-block widths.
2. Read the next `t_j` bits of each external label, and the next
   appropriate block-width bits of the channel.
3. Run the finite-block factor oracle and multiply its scalar answer
   into the running product.

Every finite-block operation costs `2^O(r)` bit time and `O(r)`
scratch bits. The label-weight search and channel counts are included
in that bound. A capped block factor uses its complete-piece oracle,
which has the same bound. Complete-piece queries multiply the
block entries in (5) and have the identical estimate.

Global positions and the two width-sum counters have `O(log lambda)`
bits. Even charging `O(log lambda)` for every block advance gives

```math
O\!\left(\lambda+
     b\bigl(2^{C_2r}+\log(\lambda+2)\bigr)\right)
\le \lambda\,2^{O(\sqrt{\log(\lambda+2)}
                               \log\log(\lambda+4))}.            \tag{9}
```

Each local label or channel is copied only once into its small
workspace. Repeatedly extracting blocks by division of one large
`lambda`-bit integer would not establish (9); the sequential bit-string
representation avoids that unnecessary work.

The packed piece tuple, a constant number of input strings, and
all local scratch fit in `O(lambda+1)` bits. Width metadata can be
recomputed on a later pass instead of storing an array with a
machine-word pointer for every block. Enumeration delay and index
validation fit (9). Enumerating the complete index set costs at most
`|I_lambda|` times that bound, which is still `2^O(lambda)` as required
by the relaxed interface.

In the global grid recursion, every active label length is at least
`n^(1/4)/8`. The larger piece descriptors therefore occupy

```math
\sum_j O\!\left(K_j^2\log|I_{\lambda_j}|\right)
 =O\!\left(\frac{\sum_j V_j}{\sqrt{\log(n+2)}}\right)
 =O\!\left(\frac n{\sqrt{\log(n+2)}}\right).
```

All local workspaces remain linear in their label length. The packed
assignment counters, streamed top queries, and paired-descriptor
recovery consequently retain the same sublinear workspace bound.
Equivalently, apply the
[simultaneous transfer](simultaneous-time-space-transfer.md) with
`alpha=delta=1/2`. Thus the nearly linear local queries and the
global `O(n/sqrt(log n))` workspace can hold in the same construction.

## 5. Conditional regularization at every fixed exponent below one

The same operation is not specific to the partition family. Suppose a
fixed-field crossing family has polynomially many pieces, separate widths
at most `2t`, total width at most `4rho t+C t^beta` for a fixed `0<=beta<1`,
and all local queries in `2^O(t)` time. Put
`alpha=log2(1/rho)` and assume `0<alpha<1`.

Choose a fixed rational `gamma` with `alpha<gamma<1`, put
`L=ceil(log_2(lambda+2))`, and use blocks of length `min(lambda,ceil(L^gamma))`.
A rational power has an exact
integer implementation by powering and root comparisons. Then

```math
\log|I_\lambda|
 =O\!\left(\frac{\lambda\log\log(\lambda+2)}
                   {(\log(\lambda+2))^\gamma}\right)
 =O\!\left(\frac{\lambda}{(\log(\lambda+2))^\alpha}\right),
```

while the width slack is `O(lambda/r^(1-beta))`, smaller than
`O(lambda/logloglambda)`. The local-query time is
`lambda 2^O((loglambda)^gamma)=lambda^(1+o(1))`. The relaxed transfer
therefore preserves `alpha` at its entropy endpoint. If the original
local oracle has no better space guarantee than its running time, its
block workspace is still `2^O(r)=lambda^o(1)`; packed global inputs and
indices give `O(lambda)` total query storage.

When `alpha` itself has an exact integer power/root implementation, the
endpoint choice `r≈(loglambda)^alpha logloglambda` gives the sharper
query overhead illustrated by (4). The rational `gamma>alpha` choice
avoids any computability assumption about a real exponent.

This is a conditional interface theorem, not a new crossing rate. In
particular it supplies no family with `alpha>1/2`. The restriction
`alpha<1` is substantive for this conclusion: the displayed block
argument would no longer have `r=o(loglambda)` at the endpoint one.

## 6. Verification

The production wrapper
[blocked_modular_factors.py](../tools/blocked_modular_factors.py)
uses the exact local oracles and contiguous bit-string channels.
The checker verifies integer block selection, piece-index coding,
tensor-product channel decoding, both factor cuts, and crossing sums
in finite cases. The asymptotic rates follow from (3), (7)--(9),
not from extrapolation of those finite tests.

```sh
python3 tools/check_blocked_half_family.py \
  --output evidence/blocked-half-family.json
```
