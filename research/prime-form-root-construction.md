# Construct the prime-form roots within the cutoff's scale

The arithmetic notes often charge prime-form construction as `Z^O(1)`.
Under the same pinned zero-free input, the exponent of `Z` in that
particular task can be one. Given a negative quadratic-order
discriminant `Delta`, every split prime form of norm at most `Z`,
coprime to `Delta`, can be constructed deterministically in

$$Z\,\mathrm{poly}(\mathrm{bits}(|\Delta|),\log Z)\tag{1}$$

bit operations. One root orientation per prime suffices; the other
is obtained by negating the middle coefficient and reducing.

This is a construction-cost refinement. It does not find a smooth
represented norm, improve the complete factoring score, or claim a
new modular square-root algorithm. It combines the existing analytic
input with the standard Tonelli--Shanks method.

## 1. A concrete bound for a quadratic nonresidue

Let `p` be an odd prime and `b_p = bits(p)`. There is a quadratic
nonresidue modulo `p` among the integers

$$2\le z\le U_p,
 \qquad U_p=\min\left(p-1,
          [10000(b_p+1)]^{32}\right).\tag{2}$$

This uses the pinned Cureton--Thorner theorem for finite-order Hecke
characters over `K=Q`. The nonprincipal character is the Legendre
character modulo `p`; no auxiliary field or unknown integer
factorization is needed.

Here are explicit constants sufficient to obtain (2). For the
principal Dirichlet character of conductor one and a primitive
nonprincipal Dirichlet character of conductor `Q`, define

$$W_\chi(x)=\sum_{n\le x}\Lambda(n)\chi(n)
                  \frac nx\left(1-\frac nx\right).$$

The proof of the
[degree-two explicit bound](explicit-prime-form-cutoffs.md)
also gives

$$|W_\chi(x)-\delta_\chi x/6|
       \le100x^{31/32}\log(Q+3),\qquad x\ge2.\tag{3}$$

For completeness, the archimedean estimates needed to use that
constant in degree one are as follows. At `s=2`, the Hadamard
logarithmic-derivative identity and
`abs(L'/L(2,chi)) < 2` give

$$\sum_\rho\mathrm{Re}\frac1{2-\rho}
       \le\tfrac12\log Q+4.$$

The pole terms sum to at most `3/2`; the gamma term is
`psi((2+a)/2)/2 <= 1/2`, where `a` is the parity. Thus the zero
sum with kernel `1/((s+1)(s+2))` is at most

$$x^{31/32}(2\log Q+16).$$

On `Re(s)=-1/2`, the functional equation involves digamma arguments
`-1/4+it/2` and `3/4-it/2` in the even case, or
`1/4+it/2` and `5/4-it/2` in the odd case. The digamma integral
bound from the preceding note, followed by one recurrence at
`-1/4+it/2`, bounds their absolute values by
`log(abs(t)+2)+8` and `log(abs(t)+2)+4`, respectively.
The older common envelope

$$|L'/L(-1/2+it,\chi)|
 \le\log Q+24+2\log(|t|+2)$$

therefore remains valid. The left integral is at most
`x^(-1/2)(log Q+30)`. At most one trivial zero at zero is crossed.
Together these estimates are less than
`x^(31/32)(3 log Q+47)`, which is bounded by (3), also for
`Q=1`. This verifies the degree-one constants instead of assuming
that a degree-two formula has the same gamma factors.

Now suppose no integer in (2) is a nonresidue, and put

$$X=[10000(b_p+1)]^{32}.$$

Every rational prime at most `X`, other than `p`, then has Legendre
symbol one. If `X >= p`, reduce that prime modulo `p` to use the
same assertion. Therefore the principal and Legendre sums differ
only at powers of `p`, giving

$$|W_1(X)-W_\chi(X)|\le\frac{\log X}{4}.$$

But (3) and `log(p+3) <= b_p+1` imply

$$|W_1(X)-W_\chi(X)|
 \ge X/6-200X^{31/32}(b_p+1)
 =11X/75>X/8>\frac{\log X}{4},$$

a contradiction. This proves (2). Testing consecutive integers
with Euler's criterion takes polynomial bit time in `log p`.
The large degree and constant in (2) are a worst-case proof;
the test stops at the first nonresidue.

## 2. Deterministic modular square roots

For a known odd prime `p`, first apply Euler's criterion to the
desired residue `a`. A negative result rejects an inert prime.
Zero is handled separately. For a nonzero square, write
`p-1 = 2^s t`, with `t` odd, and obtain a nonresidue using (2).
The usual Tonelli--Shanks state maintains

$$x^2=a b\pmod p,\qquad \mathrm{ord}(b)\mid2^{m-1},
 \qquad \mathrm{ord}(c)=2^m.$$

Initially `x=a^((t+1)/2)`, `b=a^t`, `c=z^t`, and `m=s`.
If `b != 1`, find its order `2^i`, use
`d=c^(2^(m-i-1))`, and replace

$$x\gets xd,\quad b\gets bd^2,\quad c\gets d^2,
 \quad m\gets i.$$

The order of `b` strictly decreases. Each round takes polynomial
time in `log p`, and there are at most `s` rounds. On termination,
`x^2=a`. For `p=3 mod 4`, the direct exponent `(p+1)/4` avoids
the nonresidue search.

A standard reference for this square-root procedure is
[Menezes, van Oorschot, and Vanstone, Handbook of Applied Cryptography,
Algorithm 3.34](https://cacr.uwaterloo.ca/hac/about/chap3.pdf).
The only derandomization needed here is replacing its nonresidue
choice by the bounded scan in (2).

## 3. Build the complete bounded factor base

Enumerate the rational primes `p <= Z`. For each prime not dividing
`Delta`, find a root of `Delta` modulo `p` or reject the prime as
inert. For an odd split prime choose the root representative
`B` in `[0,2p)` with `B = Delta (mod 2)`. Then

$$\left(p,B,\frac{B^2-\Delta}{4p}\right)\tag{4}$$

is a primitive positive definite form. Reduce it and its inverse,
while retaining `p` and `B` for the ideal-word converter.
For the unramified prime `2`, splitting occurs exactly when
`Delta = 1 (mod 8)`; use `B=1` in (4).

Each residue test, bounded nonresidue scan, square-root calculation,
and form reduction has cost polynomial in `bits(abs(Delta))` and
`log Z`. Prime enumeration also fits (1). For example, testing
successive integers with deterministic polynomial-time primality
testing already gives the stated bound on a finite-alphabet
Turing machine; a sieve is preferable for the finite implementation.
This argument does not require unit-cost random-access memory.

Primes dividing `Delta` are recorded separately and excluded from
the invertible factor base. At `Delta=-4N`, a proper gcd with `N`
can of course already split the factoring input.

For the balanced cutoff in
[the multiband density note](smooth-band-density.md), setup therefore
costs `L_D[1/2,sqrt(8)]` times a polynomial bit factor. This statement
concerns setup alone. Subsequent relation search and integer linear
algebra have their own costs.

## 4. Implementation and replay

The [constructor](../tools/prime_form_factor_base.py) accepts a
discriminant and finite cutoff, and returns one bounded root and
both reduced orientations for every eligible prime. Its modular
square-root helper has a known-prime input contract. The constructor
ensures that contract by its sieve.

Run the [checker](../tools/check_prime_form_factor_base.py):

```sh
python3 tools/check_prime_form_factor_base.py
```

The [evidence](../evidence/prime-form-factor-base.json) records
complete small-prime square-root comparisons, large two-power
controls, factor-base comparisons with exhaustive roots, and an
integration check with the fixed-norm orientation decoder.
The finite tests check arithmetic correctness, not the external
analytic theorem behind the worst-case scan bound.

This note used substantial OpenAI assistant support. Its pinned
analytic dependency and credit are the same as in the prime-form
generator note.
