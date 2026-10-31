# Sparse divisor polynomials can be tested in one characteristic

## Result and scope

The [rough-input refinement](rough-sparse-transfer.md) improves the additional
transfer cost below to `exp(O(sqrt n))` by testing small factors first. It
preserves every fixed exponent supplied by the same evaluator hypothesis,
removing the `min(alpha,1)` ceiling of this original support bound.

Weighted divisor-prefix queries can be tested deterministically over one
fixed characteristic. Integer counts do not need to be reconstructed.
The required family of cyclic polynomial rings has
`exp(O(n/log n))` members, each of dimension `exp(O(n/log n))`,
for an integer of at most `n` bits.

This note proves that assertion, including an explicit divisor bound,
an elementary supply of primes, cancellation control, and complete
smallest-prime-factor recovery.
The arithmetic evaluator is a separate hypothesis.
The companion [scale-one modular transfer](scale-one-modular-transfer.md)
derives that evaluator from effective modular crossings.

**Conditional transfer.** Fix a prime `p`. Suppose the weighted
divisor-prefix polynomials below have uniformly generated arithmetic
circuits of size and construction cost at most `A(n)` over `F_p`.
The circuits must be correct after substitution in every commutative
unital `F_p`-algebra. Then a deterministic finite-alphabet machine
computes the smallest prime factor in time

```math
T(n)\le
\exp\!\left(
O\!\left(\log(A(n)+2)+\frac n{\log n}\right)
\right).
\tag{1}
```

The same bound holds for complete prime factorization.
In particular, an evaluator with

```math
A(n)=\exp\!\left(O\!\left(
\frac n{(\log n)^\alpha}\right)\right)
```

gives factoring time

```math
\exp\!\left(O\!\left(
\frac n{(\log n)^{\min(\alpha,1)}}\right)\right)
\qquad(\alpha>0).
\tag{2}
```

There is no additional `log log n` loss at `alpha=1`.
Equation (2) is a transfer bound for this method, not a lower bound
against faster methods.
This note does not construct an improved crossing family or assert
that an arbitrary integer Turing program can be run on ring elements.

## 1. An explicit elementary support bound

Let `N<2^n` with `n>=2`. Define the following integers:

```math
\begin{aligned}
q(n)&=\lfloor\log_2 n\rfloor,\\
\ell(n)&=\lceil\log_2 n\rceil,\\
b(n)&=\min\left\{
n,\
\min_{1\le j\le q(n)}
\left(2^j\ell(n)+\left\lceil\frac n j\right\rceil\right)
\right\},\\
K(n)&=2^{b(n)}.
\end{aligned}
\tag{3}
```

Then the number of positive divisors satisfies

```math
\tau(N)\le K(n),\qquad b(n)=O(n/\log n).
\tag{4}
```

To prove the bound, write `N=product_p p^(a_p)` and fix an allowed `j`.
For primes at most `2^j`, there are at most `2^j` factors.
Each exponent satisfies `a_p<n`, so `a_p+1<=n<=2^ell(n)`.
Their contribution to `tau(N)` is at most `2^(2^j ell(n))`.

For primes greater than `2^j`,

```math
j\sum_{p>2^j}a_p<\log_2N<n.
```

The inequality `a+1<=2^a` for positive integers `a` bounds their
contribution by `2^ceil(n/j)`. Also `tau(N)<=N<2^n`.
Taking the best of these integer bounds proves (3)--(4).
For `N=1` the assertion follows directly.

Choosing `j=floor(q(n)/2)` for sufficiently large `n` gives
`2^j<=sqrt(n)` and `j>=log2(n)/4`. Hence

```math
b(n)\le \sqrt n(\log_2 n+1)
       +\frac{4n}{\log_2 n}+1
=O(n/\log n).
```

The minimization in (3) also gives the sharper asymptotic upper bound

```math
b(n)\le
\frac n{\log_2 n}
+O\!\left(\frac{n\log\log n}{(\log n)^2}\right).
\tag{5}
```

For this, use
`j=q(n)-3 ceil(log2(q(n)+1))` once it is positive.
The first contribution is at most `n/(q(n)+1)^2`.
The second is `n/(q(n)-O(log q(n)))+O(1)`.
No distribution estimate for primes is needed.

## 2. A sparse polynomial has few vanishing cyclic reductions

Let

```math
F(X)=\sum_{d\in S}a_dX^d\in\mathbb F_p[X],
\qquad
S\subseteq\{0,\ldots,2^n-1\},
\quad
0<|S|=t\le K,
\quad
a_d\ne0.
```

For a positive integer `r`, its image in
`R_r=F_p[X]/(X^r-1)` has coefficient at residue `u` equal to
`sum_{d mod r=u} a_d`.

Fix one exponent `d_0 in S`. If the image of `F` is zero, some other
exponent must have the same residue as `d_0`; otherwise that residue
would have coefficient `a_(d_0) != 0`.
Consequently, for prime `r`, vanishing implies

```math
r\mid P,\qquad
P=\prod_{d\in S\setminus\{d_0\}}|d-d_0|.
\tag{6}
```

The empty product is one. If `t>=2`, then
`1<=P<2^(n(t-1))`, so it has fewer than `n(t-1)` distinct prime
divisors. In all cases, at most `n(K-1)` distinct primes can give a
zero image.

Thus any

```math
M=n(K-1)+1
\tag{7}
```

distinct odd primes form a deterministic nonzeroness family:

```math
F\ne0
\quad\Longleftrightarrow\quad
\text{some }F\bmod(X^r-1)\ne0.
\tag{8}
```

The algorithm need not know `d_0`. It is used only in the proof.
This argument needs to isolate one support term, so the larger product
over every pair of exponents is unnecessary.
If `F=0`, every reduction is zero, giving the other implication.

The argument also permits `r=p`. The quotient may then have nilpotents,
but its coefficient vectors still have the basis
`1,X,...,X^(r-1)`. A field assumption is not used.
One may exclude that single prime if desired.

## 3. Generate enough small primes without analytic number theory

For every `M>=1`, there are at least `M` odd primes at most

```math
Q=4(M+1)^2.
\tag{9}
```

Here is an elementary proof. For an even integer `v>=2`, the largest
binomial coefficient in row `v` is central, so

```math
\binom{v}{v/2}\ge\frac{2^v}{v+1}.
```

For each prime `r<=v`, Legendre's formula gives

```math
v_r\binom{v}{v/2}
=\sum_{a\ge1}
\left(
\left\lfloor\frac v{r^a}\right\rfloor
-2\left\lfloor\frac {v/2}{r^a}\right\rfloor
\right)
\le\lfloor\log_r v\rfloor.
```

Each summand is zero or one, and terms vanish when `r^a>v`.
Therefore `binom(v,v/2)<=v^pi(v)`, where `pi(v)` is the number
of primes at most `v`. Taking logarithms gives

```math
\pi(v)\ge\frac{v-\log_2(v+1)}{\log_2 v}.
\tag{10}
```

For `v=Q`, use `log2(M+1)<=M` to get
`log2 Q<=2(M+1)`. Also

```math
Q-1>2(M+1)(M+2)\ge(M+2)\log_2Q.
```

Since `log2(Q+1)<=log2Q+1`, equation (10) gives
`pi(Q)>M+1`. Removing the prime two leaves more than enough primes.

Enumerate odd candidates up to `Q` and test divisibility up to their
square roots. This deterministic procedure takes a fixed polynomial
in `Q` bit operations. It is sufficient for the bound; faster sieving
is optional.

Combining (3), (7), and (9) gives

```math
\begin{aligned}
M&=O(nK(n)),\\
Q&=O(n^2K(n)^2),\\
\log M+\log Q&=O(n/\log n).
\end{aligned}
\tag{11}
```

Thus both the prime list and its construction fit the required
`exp(O(n/log n))` budget.

## 4. Weighted prefixes and removal of the divisor one

Use the prefix convention of the reference multiplication network:
`k` lower bits remain free, and the fixed upper prefix is `x`.
For `0<=k<=n` and `0<=x<2^(n-k)`, define

```math
F_{N,k,x}(X)
=\sum_{\substack{d\mid N\\\lfloor d/2^k\rfloor=x}}X^d.
\tag{12}
```

The polynomial needed for factoring is

```math
G_{N,k,x}(X)
=F_{N,k,x}(X)
-[\lfloor1/2^k\rfloor=x]X.
\tag{13}
```

This removes only `d=1`. It preserves `d=N`, which is essential for
the endpoint when `N` is prime.
All remaining coefficients are one, all exponents are less than `2^n`,
and the support has at most `tau(N)-1<=K(n)` elements.
Equation (8) therefore decides whether the prefix contains a divisor
greater than one.
In characteristic two, subtraction in (13) is addition.

The circuit hypothesis in the theorem concerns the weighted polynomial

```math
C_{N,k,x}(z_0,\ldots,z_{n-1})
=\sum_{\substack{d\mid N\\\lfloor d/2^k\rfloor=x}}
  \prod_{i=0}^{n-1}z_i^{d_i},
\qquad d=\sum_i d_i2^i.
\tag{14}
```

Substitute

```math
z_i=X^{\,2^i\bmod r}\in R_r.
\tag{15}
```

Each product becomes `X^(d mod r)`, so (14) evaluates (12) in `R_r`.
The residues `2^i mod r` are generated by repeated doubling.
No polynomial of degree `2^n` is represented.

This weighted statement requires each divisor to have exactly one valid
complete network setting. An unweighted count identity alone would not
justify distributing the weights over factor digits.
The reference supplies the needed bijection in
[`MulCount.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulCount.lean),
and the gate/digit definitions are in
[`MulNet.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean).
The [companion transfer proof](scale-one-modular-transfer.md)
spells out the weighted input gates and the general-ring interpretation.

## 5. Recover the least prime factor, including prime inputs

Start with `k=n` and `x=0`.
The eligible set contains `N`, so it is nonempty for every `N>=2`.
Repeatedly do the following:

1. Test `G_(N,k-1,2x)` in the prime family.
2. If some reduction is nonzero, replace `x` by `2x`.
   Otherwise replace `x` by `2x+1`.
3. Decrease `k` by one.

The two children partition the eligible set at the current prefix.
The selected child is always nonempty.
Whenever the zero child is nonempty, every integer in it is smaller
than every integer in the one child.
After `n` steps, `x` is the least divisor of `N` greater than one.

That divisor is prime: any proper factor greater than one of a composite
divisor would be a smaller eligible divisor of `N`.
For prime `N`, the algorithm returns `N` itself.
It does not require a separate primality oracle or a square-root cutoff.

To obtain the complete prime factorization, divide by this least prime
factor and repeat until the quotient is one.
There are at most `n` such divisions and at most `n^2` prefix queries.
All quotients can be padded to the original bit width, and the same
prime family can be reused.

## 6. Charge ring representation and the arithmetic circuit

Represent `R_r` by its `r` coefficients in the fixed field.
Addition and equality testing take `O(r)` field operations.
Schoolbook cyclic convolution takes `O(r^2)`.
Since `p` is fixed, coefficients have constant bit length.
The monomials in (15) take `O(nr)` bits to write and polynomial
time in `n,r` to construct.

A uniformly generated circuit of `A(n)` gates can be evaluated on a
finite-alphabet tape machine with a fixed polynomial overhead in
`A(n),r,n`. For example, keep every gate result in a sequential store,
scan the store to read each gate's operands, and use schoolbook
multiplication. This deliberately permits a polynomial cost for wire
access; it does not assume constant-time random access to ring elements.

There are at most `M` circuit evaluations per prefix and polynomially
many prefixes. Circuit generation, storage, prime generation, all ring
operations, and zero tests therefore have total cost

```math
\operatorname{poly}\bigl(A(n)+2,\ n,\ M,\ Q\bigr)
=\exp\!\left(
O\!\left(\log(A(n)+2)+n/\log n\right)
\right).
```

This proves (1), and (2) follows by comparing the two exponents.
If the starting bound separately charges circuit construction, replace
`A(n)` by the maximum of the gate and construction bounds.

The hypothesis requires a correct weighted circuit over the cyclic
rings. A theorem about unit-cost operations only in `F_p`, or an
integer evaluator with a bounded-height parameter, is not that hypothesis.
A crossing scalar that vanishes modulo `p` cannot be divided out.
The companion result uses scale-one identities and proves the circuit
bound before choosing the ring representation.

## 7. Reproduction and limitations

Run the finite checks with:

```sh
python3 tools/check_sparse_divisor_nonzeroness.py
```

The checker uses exact integer and finite-field arithmetic.
It checks the divisor bound on a finite interval, the central-binomial
prime-supply inequalities, sparse cancellation cases, digit-weight
substitution, and complete prefix recovery on small inputs.
Negative controls show that one fixed prime, unweighted parity, or
failure to remove the divisor one can give the wrong decision.

Its divisor enumeration is a small exhaustive reference oracle.
It is not an implementation of the hypothesized fast ring evaluator.
The all-input result is supplied by the proofs above and the separately
stated evaluator hypothesis.
The recorded output is
[sparse-divisor-nonzeroness.json](../evidence/sparse-divisor-nonzeroness.json).

The multiplication network and its divisor bijection are prior work from
the pinned factoring reference.
This sparse-polynomial argument and checker used substantial OpenAI
assistant support. No new Lean formalization was run.
