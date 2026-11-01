# A narrow-prime constructor at nearly square-root clique size

## Result

The transform modulus need not have the form `2^m+1`. A suitable prime
gives the same exact convolution with only `2w+log_2 T` residue bits.
Exact closest-vector decoding then recovers every bounded input box
without a dimension-sized margin in the modulus.

**Theorem.** Let `N>=2` have `n>=2^40` bits. Define

```math
u=\lceil\log_2(n+2)\rceil,\qquad r=\lfloor\sqrt n\rfloor,
\qquad P=\left\lfloor\frac{r}{16u^2}\right\rfloor. \tag{1}
```

For every operand-length class and divisor-prefix query, a deterministic
algorithm writes an exact Appendix B instance with

```math
B_{\rm tab}=1,\qquad P^2\lambda<13n.             \tag{2}
```

Its complete construction cost, including every boundary label, is

```math
n^{O(1)}P\,2^{(P-1)\lambda+O(A)}+2^{O(u^2)}
 =\exp\!\bigl(O(\sqrt n(\log n)^2)\bigr),         \tag{3}
```

where `A=Theta(sqrt(n)/log(n))` is defined below. More precisely,

```math
P^2\lambda=(12+o(1))n,\qquad
\log_2P=\tfrac12\log_2n-2\log_2\log_2n+O(1),
\qquad
\log_2\text{construction time}
 \le(12+o(1))\frac nP.                           \tag{4}
```

Each table entry takes `2^O(A) n^O(1)` bit operations. The constructor
supports a much larger initial clique than the former `P^40<=n`
interface. This changes the quantitative input scale and constructor
bounds. It does not improve the established factoring exponent `1/2`.

## 1. Dimensions and a narrow auxiliary prime

Let `A` be the least power of two at least `P u`. Put

```math
\begin{gathered}
D=A^2/2,\qquad T=A^2=2D,\qquad t=\log_2D,\\
w=2\left\lceil\frac{n+1}{2(T-A-1)}\right\rceil,
\qquad b=2^w,\qquad s=2w+t,\qquad X=2^s.          \tag{5}
\end{gathered}
```

All parameters use exact integer comparisons and integer roots. No
real approximation decides a ceiling or boundary condition. The
elementary estimates are

```math
\begin{gathered}
u\ge41,\quad P\ge2,\quad A\ge1024,\quad A\ge2P,
\quad A<\frac{\sqrt n}{8u},\\
\frac{A^2}{n}<\frac1{64u^2},\quad
w>64u^2,\quad s>128u^2.                          \tag{6}
\end{gathered}
```

For the lower bound on `A`, note that
`r>=2^floor((u-1)/2)-1`. At `u=41,42`, this exceeds
`16384u+16u^2`. The inequality persists as `u` increases by two:
the exponential term doubles and the polynomial terms grow by less
than two. Hence `P u>=1024`. The upper bound follows from
`A<2Pu` and (1); the width bound follows from `w>n/A^2`.

Choose the first prime `Q` congruent to one modulo `T` above `X`.
The [dyadic-prime lemma](power-two-progression-dyadic-prime.md) gives

```math
X<Q<2X.                                         \tag{7}
```

Indeed, `T>=1024^2>10^5`, and, with natural logarithms,

```math
\log X=s\log2>128u^2\log2
 >40u^2(\log2)^2\ge40(\log T)^2.                \tag{8}
```

The analytic lemma is unconditional and effective; its primary source
and constants are stated separately. Since `T` divides `X`, search
`X+1+jT` in order. Trial division suffices to certify primality.
Fewer than `2^s` candidates of `O(s)` bits occur below `2X`, so all
search work costs `2^O(s)` bit operations. This auxiliary search
does not require factoring `N`.

Find a quadratic nonresidue `g mod Q` by testing successive integers,
and put

```math
\omega=g^{(Q-1)/T}\pmod Q.                       \tag{9}
```

Then `omega^(T/2)=-1`, so its order is exactly `T`. This needs no
factorization of `Q-1` and costs at most `2^O(s)` more operations.
Every residue has the exact common bit width

```math
f=s+1=2w+t+1=2w+\log_2T.                        \tag{10}
```

## 2. Exact integer multiplication and factor assignments

Use the three-stage transform, operand-length classes, prefix tests,
and carry recurrence of the
[observable-boundary constructor](observable-boundary-constructor.md),
replacing its residue modulus by `Q` and its root by (9).
For first-factor length `ell`,

```math
d_1=\lceil\ell/w\rceil,\quad
d_2=\lceil(n-\ell+1)/w\rceil,\quad
d_1+d_2\le T-A,\quad \min(d_1,d_2)\le D.         \tag{11}
```

Cyclic convolution has no aliasing. Each genuine coefficient is at
most `D(b-1)^2<D b^2=X<Q`; its canonical inverse-transform residue
is its exact nonnegative integer value.

The paired first-transform heights at each residue index sum to at
most `A`. Each height is at most `A-1`: both supports in (11) are
positive, so each is at most `A^2-A-1`.

Keep `K=D(b-1)` and carry width `c=w+t`. Genuine carries stay in
`[0,K]` by the same induction. For arbitrary local settings retain
the separate rejection of an outgoing carry above `K`. Clipping a
rejected code never creates acceptance.

Check the exact first-factor length, truncated highest digits, the
second factor's `n-ell+1` bit cap, high zeros, prefix, target digits,
and final zero carry. The existing assignment bijection therefore
counts precisely the divisors in the requested class. A modular
equality alone is never used as integer factorization.

## 3. Bounded input recovery by exact closest vectors

Use the same balanced consecutive groups and carry relays. Their
maximum size is `kappa=ceil(A/P)<=A/2`; every block has at least two
indices. Fix an input gate with `a<=A-1` coefficients and hidden
consecutive block `I`. Remove unit twiddles from its observed outputs.

The unit Vandermonde construction gives the affine lattice
`x_0+Lambda` from the
[lattice-decoded constructor](lattice-decoded-boundary-constructor.md).
If `a` does not exceed the number of observations, modular interpolation
and all domain checks suffice. Otherwise its exact kernel basis is

```math
B=\begin{pmatrix}Q I_j&-V^{-1}W\bmod Q\\0&I_{a-j}\end{pmatrix},
\qquad j=A-|I|,                                 \tag{12}
```

where `[V W]` is the evaluation matrix. This is an `a`-dimensional
integer lattice with `O(log Q)`-bit basis entries. Root differences
are units, so `V` is invertible.

The cyclotomic norm proof applies unchanged:

```math
\lambda_1(\Lambda)\ge\sqrt{Q/2}>\frac{bA}{2}.    \tag{13}
```

That proof uses only a principal `A`th root, consecutive hidden
recipients, and integer cyclotomic coefficients. At each power-of-two
degree at least half the primitive recipients are external. Their
evaluation kernel has index `Q^j`. The principal ideal of a vanishing
nonzero folded polynomial has index divisible by `Q^j`; Hadamard's
bound then forces its norm to be at least `sqrt Q`. Folding increases
norm by at most `sqrt2`, and a successful cyclotomic division decreases
it. An external recipient at `1` or `-1` kills the final constant.
No special formula for the modulus enters.

The full digit box has center `z=((b-1)/2,...,(b-1)/2)` and radius
`R=(b-1)sqrt(a)/2`. Since `a<=A-1` and `A>=4`,

```math
2R<b\sqrt a\le b\sqrt{A-1}<\frac{bA}{2}
 <\lambda_1(\Lambda).                           \tag{14}
```

If a valid bounded input `x` exists, every other point in the affine
lattice is farther from `z`: its distance is at least
`lambda_1-R>R`. Therefore `x-x_0` is the unique closest vector to
`z-x_0` in `Lambda`.

Solve this exact Euclidean closest-vector problem with the deterministic
algorithm of Micciancio and Voulgaris. Its time and space are `2^O(a)`,
with polynomial dependence on input bit length. See *A Deterministic
Single Exponential Time Algorithm for Most Lattice Problems Based on
Voronoi Cell Computations*, SIAM J. Comput. 42 (2013), 1364--1391,
[DOI](https://doi.org/10.1137/100811970),
[full primary manuscript](https://eccc.weizmann.ac.il/report/2010/014/revision/1/download).
Corollary 5.6 states the exact CVP bound; the computational model on
printed page 7 explicitly allows polynomial dependence on coefficient
bit lengths. Double the lattice and target to make the center integral.

After recovery, check the full digit box, shortened highest digit,
leading bit, prefix, and every observation. If no valid bounded input
exists, this postcheck rejects. The middle, inverse, and carry settings
are then forced as before. Every entry is zero or one and takes
`2^O(A) n^O(1)` bit operations, including construction of the CVP
data. No preprocessed lattice oracle is provided for free.

## 4. Exact widths and the constant thirteen

The relay grouping has common channel width

```math
\lambda=6\kappa^2 f+2\kappa f+4\kappa c,
\qquad\kappa=\lceil A/P\rceil.                   \tag{15}
```

Put `epsilon=P/A<=1/u<=1/41`. Equations (5)--(6) imply

```math
\begin{aligned}
\frac{A^2w}{n}
&<\frac{1+1/n}{1-1/A-1/A^2}+\frac{2A^2}{n}\\
&\le\frac{1+2^{-40}}{1-2^{-10}-2^{-20}}
       +\frac{2}{64\cdot41^2}<1.001.             \tag{16}
\end{aligned}
```

Since `log_2 T<=log_2 n<u`,

```math
\frac{A^2 f}{n}<2.002+\frac1{64u}<2.003,
\qquad
\frac{A^2 c}{n}<1.001+\frac1{64u}<1.002.          \tag{17}
```

Use `P kappa/A<=1+epsilon` in (15):

```math
\begin{aligned}
\frac{P^2\lambda}{n}
&<6(1+\epsilon)^2(2.003)
 +2\epsilon(1+\epsilon)(2.003)
 +4\epsilon(1+\epsilon)(1.002)\\
&\le6(42/41)^2(2.003)
 +2(42/41^2)(2.003)+4(42/41^2)(1.002)<13.         \tag{18}
\end{aligned}
```

As `n` grows, `A^2w/n` tends to one, `f/w` tends to two,
`P kappa/A` tends to one, and the scaled relay terms tend to zero.
This proves the leading constant twelve in (4).

## 5. Written cost and the evaluator interface

Every complete row still has `2^((P-1)lambda)` labels. Enumerate all
of them, reject unused residue/carry/padding codes, apply the decoder,
and write its Boolean result. There are `P` rows.

Here `P=Theta(sqrt(n)/u^2)`, `A=Theta(sqrt(n)/u)`, and
`w,s,f=Theta(u^2)`. Prime and root search take `2^O(u^2)`;
one entry takes `2^O(sqrt(n)/u) n^O(1)`; a complete row has
`2^O(sqrt(n)u^2)` labels. Moreover,

```math
\frac{A}{n/P}=O(u^{-3})\longrightarrow0.
```

Thus CVP is lower order even in the leading written exponent.
Polynomially many length-class and prefix queries preserve these
bounds. Inputs below the finite cutoff can be handled by trial division
when composing a complete factoring machine.

The constructor permits explicitly computable local weights in any
commutative coefficient algebra. A nonzero grouped entry is the product
of the weights at its unique accepted setting. Charge the actual
arithmetic cost of that algebra. No additive `Q`th-root field is
required here: the prime is used only in the integer-message network
and its table decoder.

For an evaluator with circuit size

```math
2^{C P^2\lambda/(\log(P+2))^\alpha}
 (P(\lambda+1)+2)^{D_0P^2},\qquad0<\alpha<1,
```

the starting scale is now `log P=(1/2+o(1))log n`. The extra circuit
exponent `O(P^2 log n)=O(n/(log n)^3)` and the constructor exponent
`O(sqrt(n)(log n)^2)` are both lower order than
`n/(log n)^alpha`. If an evaluator queries an original table entry
on demand, recomputing it costs `2^O(A) n^O(1)` per query; multiplying
the circuit cost by this factor adds only `O(A+log n)` to its
logarithm. There is no unit-cost table-access assumption.
The contribution from the displayed main term is
at most

```math
(12C\,2^\alpha+o(1))\frac n{(\log n)^\alpha},    \tag{19}
```

with the same logarithm base in both expressions. The weighted sparse
recovery overhead must still be charged in its own theorem. No
crossing-family hypothesis is removed, and (19) does not supply an
exponent above one half.

## 6. Verification

The companion checker constructs actual auxiliary primes and roots,
compares the transform with integer convolution, checks the affine
observation lattices, and performs an exact finite closest-vector
search. It compares local entries with complete bounded-input
enumeration and checks parameter ceilings and the cutoff.

The finite closest-vector control is an exhaustive sphere enumerator.
It validates the reduction to CVP, without claiming to reimplement or
benchmark the general Micciancio--Voulgaris algorithm. The asymptotic
time bound invokes that published algorithm with its full cost charged.

Run

```sh
python tools/check_prime_transform_boundary_constructor.py --output evidence/prime-transform-boundary-constructor.json
```

The [checker](../tools/check_prime_transform_boundary_constructor.py)
and [evidence](../evidence/prime-transform-boundary-constructor.json)
record 1,536 exact operand-length product assignments, 80 additional
prime-transform products, 88 complete input maps with 1,720 recoveries,
and 1,237 exact finite sphere decodings. The 36 genuine and 72 mutated
local boundary entries agree with 55,902 reference assignments.
There are 5,167 exact schedule checks, including binary-length and
integer-square-root group-count boundaries. One explicit small case
fails the dyadic-prime interval below the analytic cutoff, as expected;
the checker does not extrapolate the analytic lemma to those inputs.
