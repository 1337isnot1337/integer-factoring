# Joint parity filters strictly improve the complete scalar multiplier family

## Result and scope

A filter with independent polynomial coefficients on three parity bits can
strictly outperform every filter of the form `1+R g` on the same three-bit
quotient. The comparison permits every Boolean support `g`, every polynomial
degree, every characteristic-two coefficient extension, and every allowed
ramification of the scalar multiplier `R`. It is not limited to coset
indicators, homogeneous supports, or a finite multiplier catalogue.

The separation persists for actual physical energies at every power-of-two
length `n>=16`, with all physical modes assigned to three balanced parity
groups. It holds uniformly over the choice of balanced groups. This gives a
growing, globally correlated multilevel family outside the complete scalar
multiplier envelope. Its full critical pressure remains positive, so the
result does not improve the factoring exponent or score.

The comparison concerns Boolean functions of these three group parities.
It does not compare against arbitrary Boolean supports on all physical bits
that fail to factor through a three-parity quotient.

## 1. A literal eight-state polynomial filter

Put `X=x^m`, where `m=nP`, and let `b=(b1,b2,b3)` be three disjoint group
occupation parities. Use a coefficient determined by their Hamming weight:

\[
 (F_0,F_1,F_2,F_3)
  =(1+X^2+X^4,\ 1+X+X^2,\ 1,\ 1+X+X^3).
 \tag{1}
\]

Every coefficient has endpoint 1. The literal three-bit transform is

\[
 \beta_b=\sum_{a\subseteq b}(1+X)^{|b|-|a|}F_{|a|}.
 \tag{2}
\]

It gives the following complete profiles.

| Parity weight | Coefficient | Transformed coefficient | Zero order | Infinity degree |
|---:|---|---|---:|---:|
| 0 | `1+X²+X⁴` | `1+X²+X⁴` | 0 | 4 |
| 1 | `1+X+X²` | `X³+X⁴+X⁵` | 3 | 2 |
| 2 | `1` | `X⁶` | 6 | 0 |
| 3 | `1+X+X³` | `X⁴+X⁶+X⁷` | 4 | 3 |

For a physical state `k`, the exact parity-group reduction from
[joint syndrome synthesis](joint-syndrome-polynomial-filters.md) gives

\[
 \beta_k=X^{|k|-|b(k)|}\beta_{b(k)}.
 \tag{3}
\]

Consequently the beta excesses over the baseline physical order are
`(0,2,4,1)` and the infinity degrees are `(4,2,0,3)`, indexed by parity
weight. All costs include the full degrees in the last column.

## 2. Exact separation at the uniform parity measure

Write `q=exp(-2 tau)` and initially give all eight parity states weight `1/8`.
The zero and infinity relative factors of (1) are

\[
 U=\frac{1+q+3q^2+3q^4}{8},\qquad
 V=\frac{q^{-4}+q^{-3}+3q^{-2}+3}{8}=q^{-4}U.
 \tag{4}
\]

At `q=3/4`, their product is

\[
 \mathcal R_{\rm joint}
 =\left(\frac{q^{-2}+q^{-1}+3+3q^2}{8}\right)^2
 =\frac{1261129}{1327104}.
 \tag{5}
\]

The [complete scalar-multiplier theorem](boolean-multiplier-automata.md)
reduces every multiplier `1+R g` of Boolean degree at most three to a finite
list of unbounded cost families. At `q=3/4`, the exact minimum over all 254
nonconstant Boolean supports is

\[
 \boxed{\mathcal R_{\rm scalar}\ge
         \frac{3970975}{3981312}
         >\frac{1261129}{1327104}.}
 \tag{6}
\]

Equality in the first inequality is attained by the indicator of the
all-zero parity state and `R=X²+X⁴`. The two constant Boolean supports have
cost at least 1: for a common coefficient `f`, its scalar cost is
`q^(ord(f)-deg(f))>=1`.

The finite proof has 42 unbounded automaton pieces for each nonconstant
support, hence 10,668 exact piece certificates. For a piece with first
divergence `h` and degree `h+delta`, the cost has the form

\[
 C+u q^{-h}+v q^h,\qquad C,u,v\in\mathbb Q_{\ge0}.
 \tag{7}
\]

For a continuous interval or ray, the minimum is at an endpoint or is
`C+2 sqrt(uv)`. The checker bounds that square root from below with a
100-bit dyadic rational and proves every lower bound is at least (6).
For an arithmetic progression of allowed integer `h`, convexity in `h`
leaves only the two progression points adjacent to its continuous minimum.
All comparisons between those points use exact rational arithmetic.
The checker also verifies the attaining literal multiplier directly.

This certificate applies to unbounded degree, coefficient field, and
ramification because the finite automaton theorem already dominates all
such multipliers by the listed families. The certificate does not impose
a coefficient-degree cutoff.

## 3. Transfer to actual growing physical intervals

Let `n>=16` be a power of two. Partition the physical mode energies

\[
 \lambda_i=2i+1,\qquad 0\le i<n,
 \tag{8}
\]

into three disjoint groups with sizes differing by at most one. Every group
then has at least five modes. Use aggregate charge `J=0` and scale

\[
 \tau=-\tfrac12\log(3/4),\qquad t=\tau/n.
 \tag{9}
\]

The physical and charged measures agree. A single occupation bit has bias

\[
 \theta_i=\tanh(t\lambda_i/2)<\tanh(\tau)=1/7,
 \tag{10}
\]

because `lambda_i<2n`. A group parity bias is the product of its bit
biases. The three parities are independent, and every bias is therefore
at most `epsilon=7^-5`. Thus every actual parity-state probability lies
between `(1-epsilon)^3/8` and `(1+epsilon)^3/8`.

For any positive state weights, either of the two cost expectations changes
by a multiplicative factor in `[(1-epsilon)^3,(1+epsilon)^3]` relative to
the uniform measure. Therefore

\[
 \mathcal R_{\rm joint,actual}
 \le(1+\epsilon)^6\frac{1261129}{1327104},
 \tag{11}
\]

whereas every scalar Boolean multiplier on these three parities satisfies

\[
 \mathcal R_{\rm scalar,actual}
 \ge(1-\epsilon)^6\frac{3970975}{3981312}.
 \tag{12}
\]

The right side of (12) exceeds the right side of (11) by more than `0.046`;
the checker proves the positive difference as a rational inequality.
The same argument permits the scalar comparator to use any other balanced
three-way partition: the bounds are uniform in the actual mode assignment.

For a Boolean polynomial in the three parities, substituting disjoint
linear parity forms preserves each homogeneous degree. Consequently its
homogeneous-layer pattern depends only on those parities, as required by
the scalar automaton theorem. This explains why the uniform-state
comparison applies even though each group grows with `n`.

## 4. Complete critical pressure stays positive

The global Gaussian support floor proves

\[
 G_n(0,t)>n\log(1+e^{-2nt})
 \quad\text{when}\quad e^{2nt}(1-e^{-t})\le1.
 \tag{13}
\]

Here `e^(2nt)=4/3` and `e^(-t)>3/4`, so the condition holds strictly.
Equation (13) is the previously certified theorem in
[the global homogeneous floor](global-homogeneous-boolean-filter-floor.md),
used only as a Gaussian inequality.

The exact checker proves `(1-epsilon)^6 R_joint>9/10`. Consequently

\[
 \exp\bigl(G_n(0,t)+\log\mathcal R_{\rm joint,actual}\bigr)
 >\frac9{10}\left(\frac74\right)^n>1
 \qquad(n\ge16).
 \tag{14}
\]

Thus this explicit family separation has strictly positive complete
critical pressure. The construction shows that independent polynomial
coefficients give a stronger search family; it does not supply the strict
negative-pressure witness needed for a score change.

## Verification

```sh
python tools/check_joint_boolean_multiplier_separation.py \
  --output evidence/joint-boolean-multiplier-separation.json
```

The evidence records all 254 support lower bounds, the complete automaton
and piece-certificate hashes, the exact attaining scalar multiplier, the
literal eight-state joint polynomials and endpoint checks, and the rational
physical-bias and positive-pressure comparisons.
