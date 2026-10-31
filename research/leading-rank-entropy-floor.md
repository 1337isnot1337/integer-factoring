# Entropy duals force one leading channel through ten letters

## Result and scope

For alphabet sizes 7<=q<=10, every bare ungauged polynomial crossing with
**at least two minimum-order channels in either cut** has

```math
\rho_{\rm bare}>17/24>1/\sqrt2.
```

The theorem holds over every field, at every degree, and with fixed auxiliary
indices on one physical side in each cut. Together with the
[unconditional floor through q=6](shared-minor-smith-floor.md), it says that
any potential improvement through q=10 must have exactly one minimum-order
channel in both cuts. There is no claim that this condition is sufficient.

This is an exact entropy dual for the determinant and shared-minor
constraints. No polynomial-degree enumeration, numerical minimization,
or compact interval cover is needed. Independent gauges, actual coefficient
recombinations, and other certificate refinements remain outside the scope.
It is not a factoring lower bound. The proved factoring score remains 1/2.

## 1. Separate the two leading multiplicities

Normalize as in the shared-minor proof: the actual degree is M, the
increasing infinity and zero profiles are mu and nu, each minimum is zero,
both sums are at most nM, and their q-entry prefixes sum to at most qM,
where n=q^2. If either minimum has multiplicity at least two, exchange the
names of the two profiles if necessary so mu_1=mu_2=0.

For M>0 put t=M log s and define

```math
p=M^{-1}\sum_{i=1}^q\mu_i,\qquad
b=M^{-1}\sum_{i=1}^q\nu_i,\qquad p+b\le q.
```

Jensen's inequality on the remaining profile groups gives

```math
A=2+(q-2)e^{-tp/(q-2)}+(n-q)e^{-t(n-p)/(n-q)},
```

```math
C=1+(q-1)e^{-tb/(q-1)}+(n-q)e^{-t(n-b)/(n-q)},\qquad
B(s)\ge e^t AC.
```

Extra zero orders only make these relaxed bounds weaker. If M=0, the
constant full-rank cuts have B=n^2 and the conclusion is immediate.

## 2. An exact dual with rational probabilities

For positive terms a_i and positive probabilities alpha_i summing to one,
weighted arithmetic-geometric mean proves

```math
\log\sum_i a_i\ge\sum_i\alpha_i\log(a_i/\alpha_i).
```

Choose alpha=(1/2,1/5,3/10). Define rational numbers

```math
\lambda=\frac{\alpha_1}{q-2}-\frac{\alpha_2}{n-q},\qquad
\beta_2=\frac{n-q}{n}(1-q\lambda)-\alpha_2,
```

```math
\beta_1=(q-1)\left(\lambda+\frac{\beta_2}{n-q}\right),\qquad
\beta_0=1-\beta_1-\beta_2.
```

For q=7,8,9,10 all these probabilities and lambda are positive. Direct
substitution gives

```math
\frac{\beta_1}{q-1}-\frac{\beta_2}{n-q}=\lambda,\qquad
1-\frac{n}{n-q}(\alpha_2+\beta_2)=q\lambda.
```

Apply weighted arithmetic-geometric mean to A and C. Every t,p,b term
then combines into t lambda(q-p-b), which is nonnegative. Consequently

```math
\log B(s)\ge H_q,
```

where the constant is explicit:

```math
H_q=\sum_{i=0}^2\alpha_i\log(m_i/\alpha_i)
    +\sum_{i=0}^2\beta_i\log(k_i/\beta_i),
\quad m=(2,q-2,n-q),\quad k=(1,q-1,n-q).
```

This explains the role of the probabilities: they cancel the degree costs
and charge the remaining shared-prefix slack with a positive multiplier.

## 3. Exact threshold comparisons

Let D be a common denominator of alpha and beta. Then exp(D H_q) is
the explicit rational product

```math
\prod_i(m_i/\alpha_i)^{D\alpha_i}
\prod_i(k_i/\beta_i)^{D\beta_i}.
```

The [checker](../tools/check_leading_rank_entropy_floor.py) raises this
rational number to the sixth power and proves it exceeds q^(17D) using
integer arithmetic. Thus H_q>(17/6)log q and rho>17/24. The integer
comparison 2*17^2-24^2=2 proves the strict half-rate separation.

| q | Approximate certified rho floor |
| --- | --- |
| 7 | 0.7311886613 |
| 8 | 0.7238831620 |
| 9 | 0.7174739264 |
| 10 | 0.7118434502 |

The [saved exact record](../evidence/leading-rank-entropy-floor.json)
contains every probability, positive multiplier, and denominator. The
displayed decimals play no role in the proof.

## 4. Consequence for weighted nilpotent flags

The preceding [quadratic-inverse constructor](inverse-degree-crossing-compatibility.md)
on Mat_7 used U(X)=J_7(1)X^T and N=I-U. Its longest nilpotent block
has length 16 and multiplicity three, so rank N^15=3. Give every input
kernel level any positive integer weight, in any kernel-adapted basis.

The longest polynomial path visits all levels. Every shorter path skips
a positive weight. Thus its leading operator coefficient is N^15 times
the invertible U, up to the fixed basis change; it has rank three. This
is the zero constant rank of the reversed column cut. The present theorem
therefore excludes **all such positive level weights and all such bases**,
extending the earlier finite 5,005-path first-jet screen.

This argument suggests a structural target: use a nilpotent with a unique
longest chain, so the leading column rank can be one. The Mat_9 transpose
channel has that property. A first quadratic-inverse path there has degree
61 and only one order-one infinity channel, so it survives the first-jet
screen. Its [complete exact profiles](nine-letter-quadratic-inverse-probe.md)
subsequently reveal 78 order-two channels and exclude that particular
path at every scale. Other noncommuting pencils remain outside this result.

The proof and exact checks used substantial OpenAI Codex assistance.
Run `python3 tools/check_leading_rank_entropy_floor.py` to reproduce them.
