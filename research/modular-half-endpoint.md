# The exact half endpoint for the normalized modular family

## Result and dependency

The normalized fixed-characteristic construction reaches the exact
`rho=1/sqrt(2)` endpoint when the recursion permits a sublinear width
slack. The elementary state-count argument already gives the quantitative
bound

```math
r_i+c_i<=2sqrt(2) lambda+O_p(sqrt(lambda)).
\tag{1}
```

The earlier [modular construction](modular-crossing.md) recorded only
`2sqrt(2)lambda+o(lambda)` and hence stated a fixed-epsilon consequence
for the logarithmic-slack interface. A sharper choice of its concentration
window proves (1), without using an individual partition-coefficient
asymptotic.

Combining (1) with the
[subexponential-piece transfer, PR #23 / 6dd2901](https://github.com/1337isnot1337/integer-factoring/blob/6dd290118668914d1f58a1d23d302963955d1dcb/research/subexponential-piece-transfer.md)
gives a deterministic finite-alphabet factoring machine with

```math
T(n)=exp(O_p(n/sqrt(log n))).
\tag{2}
```

Here `p` is any one fixed prime. The modular pieces have scale one and
work over every commutative `F_p`-algebra, including the cyclic polynomial
rings in the sparse-divisor transfer. All factor-generation and ring
costs are included by that transfer theorem. This is a composition of
ordinary proofs, with its wider-slack transfer dependency explicitly
pinned above. It is not a new Lean result and is not an exponent above
one half. The repo already has a separate complete ordinary half-rate
proof over the original integer construction.

The quantitative refinement and composition used substantial OpenAI
assistant support. The source operators and original quotient construction
retain their attribution in [modular-crossing.md](modular-crossing.md).

## 1. Quantitative coefficient existence

Fix `p` and let `q_p(h)` count partitions whose parts are not divisible by
`p` and whose multiplicities are less than `p`. As proved in Section 5 of
the modular note, with

```math
A=A_p=pi^2(1-1/p)^2/6>0,
```

the exponentially tilted distribution has partition function, mean, and
variance satisfying

```math
ln Z_p(t)=A/t+O_p(log(1/t)),\qquad
mu(t)=A/t^2+O_p(1/t),\qquad
Var_t(H)<=C_p t^{-3}.
\tag{3}
```

These estimates follow from the generating product and elementary
Riemann-sum bounds given there. Their constants depend only on the fixed
prime. Enlarge them to cover a fixed initial interval if necessary.

Choose a fixed `L_p` with `L_p^2>=2C_p`. Chebyshev gives probability at
least one half in the integer window

```math
|H-mu(t)|<=L_p t^{-3/2}.
```

The number of integers in that window is `O_p(t^{-3/2})`. At least one
integer `H` in it satisfies

```math
q_p(H)e^{-tH}/Z_p(t)>=1/O_p(t^{-3/2}).
```

Taking logarithms and using the lower endpoint of the window proves

```math
H=A/t^2+O_p(t^{-3/2}),\qquad
ln q_p(H)>=2A/t-O_p(t^{-1/2}+log(1/t)).
\tag{4}
```

This improves the previously used window `t^{-7/4}`. A constant multiple
of the standard-deviation bound suffices; no tail theorem beyond
Chebyshev is required.

## 2. The least admissible energy

Let `h(lambda)` be the least weight with `q_p(h)>=2^lambda`. Its uniform
algorithm is the exact bounded dynamic program already specified in the
modular note; the following real parameters are only used in its analysis.

Set

```math
1/t=lambda ln2/(2A)+K_p sqrt(lambda),
```

where `K_p` is a sufficiently large fixed constant. The error in the
second formula of (4) is `O_p(sqrt(lambda))`. The positive contribution
`2A K_p sqrt(lambda)` dominates it for large `lambda`, so `q_p(H)>=2^lambda`.
The first formula of (4) then gives

```math
h(lambda)<=H=(ln2)^2 lambda^2/(4A)+O_p(lambda^{3/2}).
\tag{5}
```

For the reverse estimate, the existing cumulative upper bound is

```math
ln Q_p(E)<=2sqrt(AE)+O_p(log(E+2)),\qquad
Q_p(E)=sum_{j<=E}q_p(j).
\tag{6}
```

Using `2^lambda<=q_p(h)<=Q_p(h)` and the already established quadratic
upper bound on `h` yields

```math
lambda ln2-O_p(log(2lambda))<=2sqrt(Ah(lambda))
<=lambda ln2+O_p(sqrt(lambda)).
\tag{7}
```

No monotonicity of `q_p(h)` is used. In particular, zero coefficients or
irregularities between admissible energies do not affect this estimate.

## 3. Channel widths and the complete machine bound

At a piece shift `s`, the construction has row and column channel counts
at most `Q_p(2h+s)` and `Q_p(2h-s)`, where `|s|<=2h`. Equations (6) and
(7) imply, uniformly over these shifts,

```math
\begin{aligned}
r_s+c_s
&<=\frac{2sqrt(A)}{ln2}
  (sqrt(2h+s)+sqrt(2h-s))+O_p(log(h+2))\\
&<=\frac{4sqrt(2Ah)}{ln2}+O_p(log(h+2))\\
&<=2sqrt(2)lambda+O_p(sqrt(lambda)).
\end{aligned}
```

Ceilings add at most two bits, absorbed by the error. Replace a separate
width above `2lambda` by the explicit trivial factorization if needed.
This only decreases the relevant width. The piece count remains
`4h+1=O_p(lambda^2)`. The modular note already accounts for exact piece,
factor, and label generation in `2^O_p(lambda)` finite-alphabet steps;
this quantitative analysis changes none of those algorithms.

Finally, for every sufficiently large `lambda`,

```math
O_p(sqrt(lambda))<=C_p lambda/log2 log2(lambda+4),\qquad
log2(4h+1)<=C_p lambda/log2(2lambda).
```

The pinned wider-slack transfer applies with `rho=1/sqrt(2)` and `delta=1`.
Its exponent is `min(log2(1/rho),delta,1)=1/2`, giving (2). Finitely many
small widths are covered by larger fixed constants and the same total
construction algorithms. The recursion includes its exact small-volume
cutoff, so no unproved lower bound on descendant widths is being assumed.

The proof of (1) is analytic. The existing modular identity, channel,
source-recursion, and sparse-divisor checkers verify only their stated
finite scopes; they are not presented as tests of this infinite-input
argument. A better factoring exponent still requires a leading width
coefficient strictly below `2sqrt(2)` or a different evaluator/reduction.
