# A uniform proof of the routing potential's lazy-walk tail

## Research question and result

Source: PR70 head `e5133d81dcf3ca92a7f0decaf48a33ee1c430ccd`, integrated
with PR72 at `aead326`. Its variable-level routing theorem uses a lazy
walk tail verified only on a finite grid. The next question is to prove
that tail uniformly, without enlarging the grid. Conditioning on the
number of nonzero steps and a tilted Rademacher law supply the proof.
This closes that named probability obligation; the source's expanding
instance and algorithm-accounting hypotheses are retained. It is not a
lower bound on arbitrary factoring algorithms or a stronger factoring score.

Let Z_w be a sum of w independent steps taking +1 and -1 each with
probability p, and zero with probability 1-2p, where `1/4<=p<=1/2`.
For every integer `w>=64` and `k>=0` with `k+1<=w/4`,

```
Pr[Z_w>k] >= (1/2) exp(-9(k+1)^2/(2w)-70).              (1)
```

This is stronger than Lemma B' of the source, which uses -72. The proof
allows every real p in that interval, including all its half-set choices
`p_r=floor(r/2)ceil(r/2)/(r(r-1))`, r>=2.

## 1. A deliberately loose Rademacher bound

If S_n is a sum of n independent fair signs, then for n>=1 and
`0<=k<n`, put `K=(k+1)/sqrt(n)`. We prove

```
Pr[S_n>k] >= exp(-3K^2/2-70).                           (2)
```

Two elementary logarithm bounds used below are

```
log 2 < 7/10,             log 19 < 3.
```

Finite positive partial sums for `exp(7/10)>2` and `exp(3)>19`
give exact rational certificates for these comparisons.

For n<=100, the all-positive outcome has probability `2^-n>exp(-70)`
and suffices. For larger n, the same outcome suffices whenever
`k+1>=7n/10`: indeed `(3/2)*(7/10)^2=147/200>7/10>log 2`.

It remains to consider n>100 and `k+1<7n/10`. Tilt the signs to have mean

```
mu=(k+1+2sqrt(n))/n < 9/10.
```

Under the tilted law, `E S_n=n*mu`, and its variance is at most n.
Chebyshev gives probability at least 3/4 to the event

```
|S_n-n*mu|<2sqrt(n).
```

This event is contained in `S_n>k`. For every outcome in it, the
likelihood ratio of the fair law to the tilted law is at least

```
exp(-n*D(mu)-2sqrt(n)*atanh(mu)),
D(mu)=((1+mu)log(1+mu)+(1-mu)log(1-mu))/2.
```

The power series gives

```
D(mu)=sum_(j>=1) mu^(2j)/(2j(2j-1)) <= (log 2)*mu^2.
```

Here the coefficients are positive and their sum at mu=1 is log 2.
Similarly `atanh(mu)/mu` is increasing, so for `mu<=9/10`,

```
atanh(mu)/mu <= (5/9)log 19 <5/3.
```

The likelihood exponent is consequently at most

```
(7/10)(K+2)^2+(10/3)(K+2)
  = (7/10)K^2+(92/15)K+142/15
  <= (3/2)K^2+191/9
  < (3/2)K^2+22.                                       (3)
```

The middle inequality is just completion of the square; its difference
has minimum zero at `K=23/6`. The resulting bound
`(3/4)exp(-3K^2/2-22)` is stronger than (2). This proves (2) for all
n and all allowed k. None of its constants comes from a finite tail grid.

## 2. Condition on the nonzero steps

Let `N~Binomial(w,2p)` count the nonzero steps. Conditional on N=n,
Z_w is distributed as S_n. Since

```
E N >=w/2,             Var N <=w/4,
```

Chebyshev gives

```
Pr[N>=w/3] >=1-9/w >=55/64 >1/2.
```

On that event `n>=w/3`, and `k+1<=w/4<=3n/4<n`.
Applying (2) gives a conditional lower bound

```
exp(-3(k+1)^2/(2n)-70)
  >=exp(-9(k+1)^2/(2w)-70).
```

Multiplication by `Pr[N>=w/3]>1/2` proves (1). In particular the
source's weaker -72 constant is now justified for all parameters, rather
than being a finite-grid observation.

## Consequence and controls

The [variable-level routing note](generic-engine-routing-rank-graphs.md)
can use Lemma B' as an ordinary proved input. Its separate expanding-primal-
graph premise, bounded variable occurrences, balanced groups and stated
measurement accounting remain hypotheses. This does not close the
large-Q sparse-label case or point-dependent measurement strategies.

The [checker](../tools/check_lazy_potential_tail.py) verifies the two exact
exponential-series certificates, completion of the square, conditioning
constants and rational lower-bound comparisons on independently generated
integer walk counts. The [receipt](../evidence/lazy-potential-tail-controls.json)
records those finite controls and source hashes. They support the ordinary
proof above; they are not its all-parameter justification. Run
`make verify-lazy-potential-tail`. Substantial OpenAI assistant support,
retaining the PR67/70 potential and MAIS source credits; no Lean or
independent mathematical review is claimed.
