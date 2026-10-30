# The logarithmic correction in the existing half-rate channels

## Result and scope

The partition prefactor does not hide a further logarithmic saving in the
existing declared central channels. It contributes a positive logarithmic
width correction whose coefficient can be determined exactly.

Let `p(n)` be the ordinary partition number. For positive integer `lambda`,
let `h` be the least integer at least two such that `p(h)>=2^lambda`, and put

```math
R(h)=\sum_{j=0}^{2h}p(j),\qquad
w(\lambda)=\lceil\log_2R(h)\rceil.
```

Then

```math
\boxed{w(\lambda)=\sqrt2\,\lambda+
                  (2\sqrt2-1)\log_2\lambda+O(1).}
\tag{1}
```

Both central declared channels have this width. Consequently their total is

```math
2\sqrt2\,\lambda+(4\sqrt2-2)\log_2\lambda+O(1).
\tag{2}
```

This sharpens the `O(log lambda)` remainder in the
[half-evaluator audit](half-evaluator-audit.md). It concerns the supplied
channel dimensions, not the ranks of the complete crossing pieces or a lower
bound on factoring time. It does not improve the score `1/2`.

## 1. The classical input and a cumulative asymptotic

Use the Hardy–Ramanujan asymptotic

```math
p(n)\sim\frac{e^{c\sqrt n}}{4\sqrt3\,n},\qquad
c=\pi\sqrt{2/3},\qquad \kappa=c/\ln2.
\tag{3}
```

This is the same established theorem invoked in Section 6 of
[binary-prefix-top-channel.md](binary-prefix-top-channel.md). Its primary
source is Hardy and Ramanujan, *Asymptotic formulae in combinatory analysis*
(1918), equation (1.41), [DOI: 10.1112/plms/s2-17.1.75](https://doi.org/10.1112/plms/s2-17.1.75).
The coarse bounds in the pinned evaluator alone do not give the coefficient
in (1); the full asymptotic (3) is the additional classical input.

Equation (3) implies

```math
\sum_{j=0}^{N}p(j)\sim\frac{2\sqrt N}{c}\,p(N).
\tag{4}
```

Here is a direct justification. For every fixed `A>0`, uniformly for integers
`0<=r<=A sqrt N`, equation (3) gives

```math
\frac{p(N-r)}{p(N)}=e^{-cr/(2\sqrt N)}(1+o(1))
```

in the Riemann-sum sense. Dividing the sum over this range by `sqrt N` therefore
gives the integral of `exp(-cs/2)` over `[0,A]` in the limit.
For `r<=N/2`, the same asymptotic and
`sqrt N-sqrt(N-r)>=r/(2sqrt N)` give a uniform bound
`C exp(-cr/(2sqrt N))`. For `r>N/2`, monotonicity bounds the remaining sum by
`N p(floor(N/2))`, which is exponentially small compared with `sqrt N p(N)`.
The scaled tail beyond `A sqrt N` is thus `O(exp(-cA/2))`. Letting `A` tend to
infinity proves (4), since the full integral is `2/c`.

In particular,

```math
\log_2R(h)=\kappa\sqrt{2h}-\tfrac12\log_2h+O(1).
\tag{5}
```

## 2. Invert the minimal-energy choice accurately

The neighboring ratio `p(h)/p(h-1)` tends to one by (3). Minimality of `h`
therefore gives

```math
\log_2p(h)=\lambda+o(1).
```

Applying (3) again yields

```math
\kappa\sqrt h=\lambda+\log_2h+O(1),\qquad
\log_2h=2\log_2\lambda-2\log_2\kappa+o(1).
\tag{6}
```

Substitute the first relation in (6) into (5):

```math
\log_2R(h)=\sqrt2\,\lambda+
              (\sqrt2-\tfrac12)\log_2h+O(1).
```

The second relation in (6) proves (1). Rounding up contributes less than one
bit, so it does not affect the logarithmic coefficient.

As a comparison, the single top-energy dimension satisfies

```math
\log_2p(2h)=\sqrt2\,\lambda+
                   (2\sqrt2-2)\log_2\lambda+O(1).
\tag{7}
```

Even that coefficient is positive. Equation (7) is a dimension calculation;
it does not assert that replacing the cumulative channel by its highest
energy alone preserves the crossing factorization.

## 3. Consequences for the present evaluator

A bound on these central declared widths of the form

```math
w_{row}+w_{col}\le2\sqrt2\,\lambda+C\log_2(2\lambda)+O(1)
```

requires `C>=4sqrt(2)-2`. Thus a sharper use of the partition prefactor does
not produce a negative width correction that can accumulate into a second
saving.

The separate [recurrence barrier](half-recurrence-barrier.md) already gives
`e_ell>=c/sqrt(log n)` for the square-grid coefficient template with a fixed
positive terminal coefficient. This excludes multiplication of its saving
by any function tending to infinity, including a power of `log log n`.
It is a restriction on that analysis template, not a lower bound on the
underlying machine. The [rectangular minimax theorem](rectangular-grid-frontier.md)
likewise retains a constant multiple of `1/sqrt(log K)` even with arbitrary
finite rectangular splits in its stated profile model.

A genuine additional saving would have to use information absent from these
channel counts and scalar recurrences. Examples include joint factor
compression, an input-specific cancellation, or a different representation of
the recursively generated tables. A generic balanced-tree decomposition of an
arbitrary terminal tensor does not automatically supply that information: its
standard full-rank bond counts add one linear information-volume charge per
balanced tree level. This last observation describes that construction only;
it is not an impossibility theorem for succinct table representations.

## 4. Exact finite evidence

Run from the repository root:

```bash
python3 tools/check_half_second_order.py --output /tmp/half-second-order-verification.json
```

The checker computes partition numbers by the exact pentagonal recurrence and
compares the first 81 values with an independent coin-change computation.
It verifies minimal energies and exact integer channel widths for label lengths
8 through 512. The recorded residuals are numerical diagnostics; equations
(3)–(6) supply the asymptotic proof. No Lean formalization is used.
