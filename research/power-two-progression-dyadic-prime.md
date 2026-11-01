# A narrow prime modulus in a power-of-two progression

## Explicit lemma

Let `q` be a power of two with `q>=100000`. Let `x` be real and satisfy

```math
\log x\ge40(\log q)^2.
```

For every odd residue `a`, there is a prime

```math
x<p\le2x,\qquad p\equiv a\pmod q.                 \tag{1}
```

Here and below, unqualified logarithms are natural logarithms. All
constants in this statement are explicit. In particular, with `x=2^s`
and `a=1`, the prime has exactly `s+1` binary digits and supports every
root-of-unity order dividing `q`.

The proof uses an explicit prime-number estimate and an elementary
exclusion of real exceptional zeros for power-of-two moduli. It does
not assume the generalized Riemann hypothesis.

## 1. The analytic input

Bennett, Martin, O'Bryant, and Rechnitzer, *Explicit bounds for primes
in arithmetic progressions*, Illinois Journal of Mathematics 62
(2018), 427--532, prove the following estimate in Lemmas 6.10--6.12,
printed pages 501--504:

```math
\left|\theta(y;q,a)-\frac{y}{\varphi(q)}\right|
\le \frac{14579}{10000}\,y\,X e^{-X}
    +\frac{1001}{1000}\sqrt y,\qquad
X=\sqrt{\frac{\log y}{R_1}},
\quad R_1=9.645908801,                             \tag{2}
```

provided `q>=100000`, `y>=exp(4 R_1 (log q)^2)`, and the associated
Dirichlet functions have no exceptional real zero. Here
`theta(y;q,a)` sums `log p` over primes `p<=y` in the residue class.
Their possible exceptional-zero term is absent under that hypothesis.
The estimate is effective and unconditional under the stated
zero-exclusion condition.

Primary source:
[author-hosted paper](https://personal.math.ubc.ca/~bennett/BeMaObRe-IJM-2018.pdf),
especially Definition 6.1 and Lemmas 6.10--6.12.
The original paper is also listed as
[arXiv:1802.00085](https://arxiv.org/abs/1802.00085).

## 2. Why the exceptional zero is absent

For `q=2^k`, every nonprincipal real character is induced by one of
the three primitive real characters of conductors four or eight.
Indeed, for `k>=3`, the unit group is
`C_2 times C_(2^(k-2))`. A character taking only the values `1,-1`
factors through the quotient modulo eight.

Their signs on odd residues are

| Character | Modulus | Signs in increasing odd-residue order |
|---|---:|---|
| `chi_(-4)` | 4 | `1,-1` |
| `chi_8` | 8 | `1,-1,-1,1` |
| `chi_(-8)` | 8 | `1,1,-1,-1` |

For every real `sigma>0`, the corresponding Dirichlet series
converges by the bounded-partial-sums test. Its consecutive blocks
are strictly positive.

For `chi_(-4)`, each block is

```math
(4j+1)^{-\sigma}-(4j+3)^{-\sigma}>0.
```

For `chi_(-8)`, write a block as

```math
\bigl((8j+1)^{-\sigma}-(8j+5)^{-\sigma}\bigr)
+\bigl((8j+3)^{-\sigma}-(8j+7)^{-\sigma}\bigr)>0.
```

For `chi_8`, the block is

```math
(8j+1)^{-\sigma}+(8j+7)^{-\sigma}
-(8j+3)^{-\sigma}-(8j+5)^{-\sigma}>0.              \tag{3}
```

The last inequality follows from strict convexity of `t^(-sigma)`:
the two pairs have the same midpoint, and the first pair is more
separated. Grouping consecutive blocks does not change a convergent
series. Thus each of these three `L`-functions is positive on the
entire positive real axis.

Inducing one of these characters to a larger power of two does not
change its Dirichlet series: both characters vanish on even integers
and agree on odd integers. Consequently no nonprincipal real
character modulo `q` has a positive real zero. This proves the
zero-exclusion hypothesis required by (2).

## 3. A rational error margin

Since `R_1<10`, the hypothesis `log y>=40(log q)^2` implies both the
threshold in (2) and

```math
X>2\log q.
```

The function `z e^-z` decreases for `z>=1`. Since `q>=100000`, (2)
therefore gives

```math
\frac1y\left|\theta(y;q,a)-\frac{y}{\varphi(q)}\right|
<\frac{14579}{5000}\frac{\log q}{q^2}
 +\frac{1001}{1000\sqrt y}.                       \tag{4}
```

For `q>=1`, the elementary inequality `log q<=sqrt q` follows by
maximizing `log q/sqrt q`, whose maximum is `2/e<1`. Also the
hypothesis on `y` implies `y>=q^4`. Hence the first term in (4) is
less than `1/(32q)`, because

```math
\sqrt q>100,\qquad
32\cdot14579=466528<500000.
```

The second term is also less than `1/(32q)`, because

```math
q\ge100000,\qquad
32\cdot1001=32032<100000000.
```

We have proved the explicit estimate

```math
\left|\theta(y;q,a)-\frac{y}{\varphi(q)}\right|
<\frac{y}{16q}                                  \tag{5}
```

throughout the stated range. Apply (5) at `y=x` and `y=2x`.
Because `phi(q)=q/2`,

```math
\theta(2x;q,a)-\theta(x;q,a)
>\frac{x}{\varphi(q)}-\frac{3x}{16q}
=\frac{29x}{16q}>0.
```

This proves (1).

## 4. An integer schedule for the multiplication constructor

Here is a concrete schedule requiring no unknown analytic cutoff.
Let the factoring input length satisfy `n>=2^40`, and set

```math
\begin{gathered}
u=\lceil\log_2(n+2)\rceil,\qquad
P=\left\lfloor\frac{\lfloor\sqrt n\rfloor}{16u^2}\right\rfloor,\\
A=\text{least power of two at least }Pu,\qquad
T=A^2,\quad D=T/2,\quad t=\log_2D,\\
w=2\left\lceil\frac{n+1}{2(T-A-1)}\right\rceil,
\qquad b=2^w,\qquad s=2w+t.
\end{gathered}                                                   \tag{6}
```

Then

```math
P\ge8,\qquad T\ge100000,\qquad
A^2<\frac{n}{64u^2},\qquad w>64u^2,\qquad s>128u^2.
                                                               \tag{7}
```

We give the small-cutoff details. Since `u>=41` and
`2^(u-1)<n+2<=2^u`, one has `n>2^(u-2)`. For every integer
`u>=41`,

```math
2^{u-2}>(128u^2+1)^2.                            \tag{8}
```

The base case is an integer comparison. For the induction step, the
left side doubles while the right side grows by a factor smaller
than `(42/41)^4<2`. Thus `floor(sqrt n)>=128u^2`, proving `P>=8`.
It follows that `A>=Pu>=328` and `T>=328^2>100000`.

Since `A<2Pu` and `P<=sqrt n/(16u^2)`,

```math
A<\frac{\sqrt n}{8u}.
```

Also `w>n/T`, because its denominator before even rounding is
`T-A-1<T`. These inequalities prove the remaining bounds in (7).
They also show that the denominator in (6) is positive and that
each contiguous group contains at least two indices.

Now `T<n<2^u`. Using `log 2<1`,

```math
\log(2^s)=s\log2
>128u^2\log2
>40u^2(\log2)^2
>40(\log T)^2.
```

The lemma therefore supplies a prime

```math
2^s<Q<2^{s+1},\qquad Q\equiv1\pmod T.             \tag{9}
```

The upper endpoint is a power of two and cannot be this odd prime.
In particular,

```math
Q>Db^2,\qquad
\operatorname{bitlength}(Q)=s+1=2w+t+1.           \tag{10}
```

The modulus therefore exceeds the exact convolution bound while
using only one more bit than `2w+t`. It removes the multiple-of-`D`
rounding in the Fermat-type modulus.

The schedule has

```math
P=\Theta\!\left(\frac{\sqrt n}{(\log n)^2}\right),
\qquad
A=\Theta\!\left(\frac{\sqrt n}{\log n}\right),
\qquad
w=\Theta((\log n)^2).
```

In particular `log P=(1/2+o(1))log n`. The
[prime-transform constructor](prime-transform-boundary-constructor.md)
charges all boundary labels and the separate lattice decoding cost.
The prime-existence lemma alone is not a factoring algorithm.

## 5. Deterministic prime and root construction

Search the integers congruent to one modulo `T` in the finite interval
(9). Trial division is already enough to test primality in total time
`2^O(s)` bit operations; no probabilistic primality test is needed.
Equation (9) proves that this finite search succeeds.

For the found prime `Q`, scan `g=1,2,...` until

```math
g^{(Q-1)/2}=-1\pmod Q.
```

A quadratic nonresidue exists, so this scan also costs `2^O(s)` bit
operations. If `v_2(Q-1)=r`, such a `g` has order with exact two-part
`2^r`. Since `T` is a power of two dividing `Q-1`, the element

```math
\omega=g^{(Q-1)/T}\pmod Q
```

has exact order `T`. This constructs the principal transform root
without factoring `Q-1`. Its distinct powers have unit differences
because the modulus is prime.

Under (6), the complete prime and root search costs
`exp(O((log n)^2))`. Its cost is part of the constructor, and is
smaller than the explicitly written boundary-table cost in the
companion theorem.

## 6. Exact controls

The [checker](../tools/check_power_two_progression_dyadic_prime.py)
and [record](../evidence/power-two-progression-dyadic-prime.json)
check the five rational constant comparisons, 16,368 complete
character products, 360 positive rational character blocks, and
2,363 integer schedules around binary-length boundaries from
`2^40` through `2^512`.

```sh
python3 tools/check_power_two_progression_dyadic_prime.py \
  --output evidence/power-two-progression-dyadic-prime.json
```

The prime-existence proof is the analytic estimate (2) and the
uniform arguments above. These controls verify its explicit numerical
margins and the integer schedule; they do not replace the analytic
theorem by searches for small primes.
