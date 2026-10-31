# Polynomial exact divisor preprocessing from Nicolas's exponent theorem

## Result and scope

Write

\[
D(n)=\max_{1\le m<2^n}\tau(m),\qquad
M(n)=\min\{m<2^n:\tau(m)=D(n)\},\qquad
q(n)=\lceil\log_2D(n)\rceil.
\]

**Theorem.** There is one deterministic algorithm which computes `D(n)`,
`M(n)`, and `q(n)` in `n^O(1)` bit time and `O(log(n+2))` working bits,
with read-only input and streamed output. The exponent and constants are
absolute. The same resources suffice to stream every highly composite
number below `2^n`, in increasing order.

The new ingredient relative to
[the next-prime exchange algorithm](exchange-divisor-preprocessing.md)
is a constant-width description of all sufficiently high exponent levels.
The source is Nicolas's 1971 uniform exponent theorem. Uniform iterated
multiplication supplies the small-space exact comparisons, as in the
exchange note.

This is an **unconditional uniform complexity corollary**, with fixed
constants whose numerical values are not instantiated here. The proof
establishes existence of one finite program; it is not a supplied practical
polynomial-time divisor oracle. The constants do not vary with `n` and
are not an advice string. The executable contribution is the integer-only
threshold-center procedure and its independent exact arithmetic checks.
The previously supplied practical oracle retains its documented
`n^O(log log n)` time bound.

No improvement to the factoring exponent is claimed. The strongest
complete factoring bound remains `alpha=1/2`. Computing these quantities
as functions of the input length does not factor a supplied integer.
No novelty claim is made for the classical highly-composite estimates or
for this type of algorithmic consequence.

## 1. Precisely which classical statement is used

Jean-Louis Nicolas, *Répartition des nombres hautement composés de
Ramanujan*, Canadian Journal of Mathematics **23** (1971), 116–130,
[DOI 10.4153/CJM-1971-012-6](https://doi.org/10.4153/CJM-1971-012-6),
[primary PDF](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/4E3B60B2BC31651C55D01FFC9ECB4497/S0008414X00053037a.pdf/repartition_des_nombres_hautement_composes_de_ramanujan.pdf),
Theorem 2 and its proof on printed pages 124–125, imply the following.

There are fixed `gamma>0`, `C>0`, and a finite exceptional range such
that, for every sufficiently large highly composite integer `A`, with
**actual largest prime factor** `p`, and every prime `lambda<p`, putting
`b=v_lambda(A)` gives

\[
\log\left(1+\frac1{b+1}\right)-Cp^{-\gamma}
\ \le\ \frac{\log2\,\log\lambda}{\log p}
\ \le\ \log\left(1+\frac1b\right)+Cp^{-\gamma}.
\tag{1}
\]

The constants are uniform in `lambda` and `b`. This uniformity is
essential. Equation (1) does not merely concern each fixed exponent
separately.

Nicolas first bounds the benefit relative to a superior highly composite
number `N_epsilon`, with `x=2^(1/epsilon)`. His proof then uses the
uniform relation `p-x=O(x^tau)`, `tau<1`, and a choice
`gamma<1-tau` to replace `epsilon log lambda` by
`log(2) log(lambda)/log(p)`. The error remains `O(p^(-gamma))` because
`log lambda<=log p`. Thus the parameter in (1) really is the largest
prime of `A`; no auxiliary real superior parameter must be found by the
algorithm.

The classical largest-prime exponent is one, apart from `A=4,36`.
Therefore (1) also holds at `lambda=p` once those finite exceptions are
absorbed: its middle term is exactly `log 2`. The proof below only needs
the case `lambda<p` for levels at least two.

We may decrease `gamma` to a positive rational and increase `C` to an
integer. Those are fixed choices. No numerical effective bound for these
choices is asserted here. In particular, we do not silently extract an
implemented constant from the Baker/Feldmann or prime-gap ingredients in
Nicolas's proof.

## 2. High levels have a constant-width integer description

A least maximizer `M(n)` is highly composite: an earlier integer with
at least as many divisors would contradict maximality or the least-witness
tie rule. Its prime exponents are nonincreasing and use consecutive primes.
For any highly composite `A`, define its conjugate heights

\[
h_a=\#\{\lambda\text{ prime}:v_\lambda(A)\ge a\},
\qquad a\ge1.
\]

For a prime `p` and a positive level `a`, set

\[
\theta_a(p)
 =p^{\log(1+1/a)/\log2},\qquad
\delta(p)=\frac{C}{\log2}\,p^{-\gamma}\log p.
\tag{2}
\]

**Threshold consequence.** For every `a>=2`, all primes strictly below
`theta_a(p) exp(-delta(p))` have exponent at least `a`, and all primes
strictly above `theta_a(p) exp(delta(p))` have exponent below `a`.
For sufficiently large `p`, both thresholds relevant here are below `p`.

To prove the first assertion, suppose `b<a`. Then `b+1<=a`, so the left
side of (1) gives

\[
\frac{\log2\,\log\lambda}{\log p}
 \ge\log(1+1/a)-Cp^{-\gamma}.
\]

This forces `lambda>=theta_a exp(-delta)`. Conversely, `b>=a` and the
right side of (1) force `lambda<=theta_a exp(delta)`. Endpoints are left
inside the ambiguous interval; no strict comparison at a transcendental
threshold is required.

Choose a fixed integer `K>=2` so large that

\[
\frac{\log(1+1/(K+1))}{\log2}<\frac\gamma2.
\tag{3}
\]

For every `a>K`, `theta_a(p)<=p^(gamma/2)`. Consequently

\[
\theta_a(p)\bigl(e^{\delta(p)}-1\bigr)
 =O\!\left(p^{-\gamma/2}\log p\right)=o(1),
\tag{4}
\]

uniformly in all these levels. Fix an integer `P0` large enough that
(1) applies whenever `p>=P0`, that the thresholds for `a>=2` stay below
`p`, and that both endpoints in (4) are within `1/4` of `theta_a(p)`.
All these requirements hold eventually, so one such fixed `P0` exists.

The following arithmetic lemma supplies an integer center `j_a(p)` in
polynomial time and logarithmic workspace:

\[
|j_a(p)-\theta_a(p)|<2.
\tag{5}
\]

Since an interval between two integers at distance three contains at most
three additional primes, (4)–(5) imply the uniform **descriptor cover**

\[
\boxed{\quad
h_a\in\{\pi(j_a(p))-3,\ldots,\pi(j_a(p))+3\}
\quad(a>K,\ p\ge P0).\quad}
\tag{6}
\]

The three-unit margin is deliberately generous. It absorbs all endpoint
rounding, including the possibility that a threshold is extremely close
to an integer or prime. No lower bound for a linear form in logarithms
is needed by the arithmetic subroutine.

## 3. Integer-only centers in logarithmic workspace

Here is a concrete procedure for (5). It works for every integer `p>=2`
and `a>=1`; primality is unnecessary. Let

\[
L=\operatorname{bitlen}(p)+\operatorname{bitlen}(a)+10,
\quad F=32L,\quad \eta=2^{-F}.
\]

Represent real approximations as integers scaled by `2^F`. Round every
multiplication and division down immediately.

For `1<=r<=2`, use

\[
\log r
 =2\sum_{i\ge0}\frac{z^{2i+1}}{2i+1},\qquad
z=\frac{r-1}{r+1}\in[0,1/3].
\tag{7}
\]

Compute `z` to `F` fractional bits and retain the first `F` terms,
updating each power by one truncated multiplication by `z^2`. All terms
are nonnegative. If the output, divided by `2^F`, is `ell_F(r)`, then

\[
0\le\log r-\ell_F(r)\le(8F+4)\eta.
\tag{8}
\]

For completeness, rounding `z` costs at most `(9/4)eta` in its logarithm.
The error in each power is at most `(3/2)eta`: its recurrence contracts
by at most `1/9`, and each update introduces at most `(4/3)eta`.
Term division and the factor two cost at most `5F eta` in total. The
positive omitted tail is below `9^(-F)<=eta`. These bounds imply (8).

For an integer `v`, write `v=2^s r`, `1<=r<2`, and compute
`s*ell_F(2)+ell_F(r)`. This remains a lower approximation, with error
at most `(s+1)(8F+4)eta`. Compute `log(1+1/a)` directly from (7),
whose argument is `z=1/(2a+1)`.

Binary-search integers in `[1,p]`, testing the sign of

\[
f(j)=\log j\,\log2-\log p\,\log(1+1/a)
\tag{9}
\]

with these truncated approximations. The two products are truncated as
well. Their combined absolute error is less than

\[
2(L+1)^2(8F+4)\eta+2\eta
 <2^{-16L}<\frac1{16p}.
\tag{10}
\]

The middle inequality holds for `L>=10`, which is built into the choice
of `L`. On `[1,p]`, the derivative is
`f'(j)=log(2)/j>=1/(2p)`. Therefore an incorrect sign can occur only
within distance `1/8` of the real zero `theta_a(p)`.

Start with `lo=1`, `hi=p+1`. If the approximated value at `mid` is
nonpositive, set `lo=mid`; otherwise set `hi=mid`. Throughout the search,

\[
lo\le\theta_a+1/8,\qquad hi\ge\theta_a-1/8.
\]

When `hi=lo+1`, this gives
`theta_a-9/8<=lo<=theta_a+1/8`, proving (5).

All stored integers have `O(L)` bits, including temporary multiplication
numerators. There are `O(L)` series terms and `O(log p)` search steps.
Ordinary integer arithmetic therefore gives polynomial bit time and
`O(log(p+a))` working bits. The code
[`nicolas_threshold_centers.py`](../tools/nicolas_threshold_centers.py)
implements this procedure directly. The proof's descriptor generator can
recompute centers when required, so it need not retain a list of them.

## 4. Polynomially many candidates, each with a logarithmic descriptor

The exchange lemma already proved in the preceding note gives, for a
least maximizer with largest prime `p` and next prime `q`,

\[
p\le4n^2=:B,\qquad
h_1=\pi(p)\le n-1,\qquad
h_a=0\quad\text{for }a>E_p,
\tag{11}
\]

where we may take

\[
E_p=\min\{n-1,\;2\lceil\log_2(2p)\rceil-2\}.
\tag{12}
\]

Indeed, the sharper exchange bound is
`e_1<=2 ceil(log2 q)-2`, and Bertrand gives `q<2p`.
The powers and ceilings in (11)–(12) are computed by integer arithmetic.

Always include the empty profile, giving `M=D=1`. For each prime
`p<=B` with `pi(p)<=n-1`, enumerate the following descriptors.

- If `p<P0`, enumerate all nonincreasing positive exponent vectors on the
  `pi(p)` primes with exponents at most `E_p`. Their total number is a
  fixed finite constant depending on `P0`, independent of `n`.
- If `p>=P0`, fix `h_1=pi(p)`. For levels `2<=a<=min(K,E_p)`,
  allow any height from zero to `h_1`. For `K<a<=E_p`, choose one
  of the seven offsets in (6).

Reject any descriptor with negative or excessive heights, a nonmonotone
height sequence, `sum h_a>n-1`, or a reconstructed integer at least `2^n`.
Every least maximizer is retained by (6) and (11). The small-prime branch
covers the finite exceptional range without a table depending on `n`.

There are at most

\[
B(n+1)^{K-1}7^{\max_{p\le B}E_p}=n^{O(1)}
\tag{13}
\]

large-prime descriptors. For example,
`E_p<=4 ceil(log2 n)+6`, so the candidate count alone is bounded by
`C_K n^(K+13)` for `n>=2`. We do not use this displayed exponent as a
complete bit-runtime constant, since the arithmetic comparisons add a
fixed polynomial factor.

The largest-prime counter uses `O(log n)` bits. There are only a fixed
number `K-1` of unrestricted heights, each using `O(log n)` bits.
The remaining `O(log n)` levels use three bits per offset. Thus the whole
descriptor, its enumeration counter, and a retained best descriptor use
`O(log n)` bits. Coordinate offsets, centers, and prime counts can be
regenerated. Elementary trial division over integers at most `B` computes
prime counts and individual primes in polynomial time and logarithmic
workspace.

Choose once and for all integers `K,P0` satisfying the proved conditions.
This defines one fixed finite algorithm. Their values do not need to be
computed anew from the input. Their non-instantiation in this note makes
the result an existence theorem about uniform complexity, and does not
turn it into nonuniform advice. An explicit efficient implementation
would need numerical global constants or another certified coverage rule.
Neither is claimed by the finite experimental checks below.

## 5. Exact objective comparisons and output

From a descriptor, regenerate the exponent of the `i`th prime as

\[
e_i=\#\{a:h_a\ge i\}.
\]

Its candidate integer is a product of at most `n-1` prime factors,
with repetition, each at most `4n^2`. Before feasibility is checked, the
product has at most `O(n log n)` bits. Its divisor count is the product
of the small values `e_i+1`, with at most `n` bits under the exponent-sum
bound. Each requested bit of the list of factors can be regenerated in
polynomial time and `O(log n)` workspace from the compact descriptor.

Hesse, Allender, and Barrington, *Uniform Constant-Depth Threshold Circuits
for Division and Iterated Multiplication*, JCSS **65** (2002), 695–716,
[primary paper](https://people.cs.rutgers.edu/~allender/papers/division.pdf),
Corollary 6.3, place iterated multiplication in DLOGTIME-uniform `TC0`.
Evaluating the fixed-depth polynomial-size circuits with recomputed gate
values uses polynomial time and logarithmic workspace. Majority counts,
gate addresses, and the constant-depth evaluation stack all fit in
`O(log n)` bits. Replacing the circuit input by the bit generator above
preserves these bounds.

Therefore candidate products and divisor counts can be compared exactly,
bit by bit, within the same resources. Retain the descriptor with greatest
divisor count, breaking ties by the smaller integer. Coverage proves it
is `M(n)` and has value `D(n)`. Regenerate and stream their bits.
The bit length and power-of-two test for `D(n)` give `q(n)` in the same
bounds. The case `n=1` is handled separately.

For the list variant, the cover contains every highly composite integer
below `2^n`, not only its last one. A candidate `A` is a strict record
exactly when no smaller covered candidate has divisor count at least
`tau(A)`: any counterexample among all integers would itself be dominated
by a smaller highly composite integer, which is covered. Nested scans
therefore filter the records in polynomial time. Repeated scans selecting
the least remaining record stream them in increasing order. These are
only a fixed number of nested polynomial-sized searches, and their
constant number of descriptors still occupies `O(log n)` bits.

## 6. What the executable audit establishes

Run

```sh
python3 tools/check_nicolas_threshold_centers.py
```

The checker independently brackets logarithms by exact rational atanh
series with explicit positive tails. It certifies the threshold centers
by the signs two integers to either side, and checks the fixed-point
error bound on a separate grid. It also uses an ordinary divisor sieve
to inspect every strict record below `2^16`.

The illustrative finite cover uses cutoff `K=2` and offsets of magnitude
at most three. This is a finite control of the descriptor mechanics;
it is **not** a proof that `K=2` is a valid global Nicolas cutoff. The
recorded evidence states that limitation explicitly. No floating-point
number, numerical minimization, or inferred asymptotic fit certifies the
theorem.
