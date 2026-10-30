# A polynomial-size characteristic for quadratic factor recovery

## Result and scope

The exact root count in
[quadratic-congruence-factor-recovery.md](quadratic-congruence-factor-recovery.md)
is `2^(omega(N)-1)-1`. Characteristic two detects its nonzero cases
automatically. When a counting construction permits choosing its
characteristic, a small odd prime works as well, without replacing the
quadratic verifier by a prime-length tuple system.

For every integer `n>=1` and every explicitly supplied nonzero integer
`D`, one can deterministically find an odd prime `p` such that

```math
p\nmid D,\qquad \operatorname{ord}_p(2)>n,
\qquad p\le \max\{256,4(h+1)\},                                  \tag{1}
```

where

```math
h=1+\frac{n(n+1)}2+\operatorname{bitlength}(|D|).                   \tag{2}
```

The algorithm takes polynomial time in `n+bitlength(|D|)` and never
factors `D` or any number `2^j-1`. For `D=1`, the selected prime is
`O(n^2)` and its bit length is `O(log n)`. The order condition also
implies `p>n+1`.

For every input of at most `n` bits, this characteristic makes the same
one-scalar quadratic verifier a valid nonzero-count recovery oracle.
All smaller recursive inputs use the originally selected prime. There
is no rough-input trial-division stage, sparse-polynomial detector, or
tuple-length dependence on `p`.

This is a uniform field-selection and recovery theorem. It does not
construct an improved crossing tensor or prove a factoring exponent
above one half. It applies when the evaluator is available uniformly
over the selected primes, with its field arithmetic and coefficient
construction charged. A construction restricted to a prescribed fixed
finite field can instead use the unitary theorem; the present selection
does not move that construction to a new characteristic automatically.

## 1. An elementary primorial bound

Put

```math
L_m=\operatorname{lcm}(1,\ldots,m),\qquad
P_m=\prod_{\substack{p\le m\\p\ \mathrm{prime}}}p.
```

For every `0<=j<=m-1`,

```math
m{m-1\choose j}\mid L_m.                                          \tag{3}
```

To verify (3), fix a prime `p` and put `a=v_p(m)`. In the usual
factorial-floor formula for `v_p(binom(m-1,j))`, each summand is zero
or one. For powers `p^r` with `r<=a`, it is zero: `m-1` has residue
`p^r-1`, so the two lower factorial arguments have remainders summing
to exactly `p^r-1`. For `r>floor(log_p m)`, it is also zero. Thus the
valuation of the binomial coefficient is at most
`floor(log_p m)-a`. Multiplying by `m` leaves valuation at most
`floor(log_p m)`, which is the exponent of `p` in `L_m`.

At least one of the `m` binomial coefficients in row `m-1` is at least
`2^(m-1)/m`. Equation (3) gives

```math
L_m\ge2^{m-1}.                                                    \tag{4}
```

Only primes at most `sqrt(m)` occur to powers greater than one in
`L_m`. For each of them, the prime-power contribution to `L_m/P_m`
is at most `m`. There are at most `floor(sqrt(m))` such primes, so

```math
\log_2 P_m
\ge m-1-\sqrt m\log_2m.                                          \tag{5}
```

For real `m>=256`, `log_2 m<=sqrt(m)/2`: equality holds at 256 and
the difference thereafter is increasing, as differentiation verifies.
Consequently

```math
P_m\ge2^{m/2-1}\qquad(m\ge256).                                  \tag{6}
```

This deliberately coarse bound suffices. No prime number theorem or
unproved statement about primitive roots is required.

## 2. Find a prime whose base-two order is large

Define the integer only for the proof:

```math
B_n=\prod_{j=1}^n(2^j-1).
```

Every odd prime with `ord_p(2)<=n` divides `B_n`. Also,

```math
2|D|B_n<2^h                                                     \tag{7}
```

with `h` from (2). Set `X=max(256,4(h+1))`. By (6), `P_X>2^h`.
If every prime at most `X` either were two, divided `D`, or had
base-two order at most `n`, the squarefree primorial `P_X` would divide
`2|D|B_n`. Equation (7) gives a contradiction. Therefore a prime
satisfying (1) exists below `X`.

The deterministic selection algorithm enumerates candidate primes up to
`X`. For each odd candidate, test `D mod p` and compute the residues
`2^1,...,2^n mod p` by successive modular doubling. Accept the first
candidate for which the former is nonzero and none of the latter is one.

Primality can be checked by ordinary trial division, or candidates can
be generated with a sieve. Since `X=O(n^2+bitlength(|D|)+1)`, either
choice and every modular test take polynomial time in the stated input
lengths. The proof number `B_n` need not be constructed. In particular,
no factorization of Mersenne numbers or of the forbidden-prime
certificate is hidden in this selection step.

The integer `D` may be the product of a fixed finite list of exceptional
primes, or any available nonzero denominator/discriminant certificate
whose prime divisors contain the forbidden characteristics. It must be
supplied or generated with its full bit cost charged. An unspecified or
exponentially expensive exceptional-prime list cannot be used for free.

## 3. Why this field detects the quadratic witness count

Take an input `N>1` with `bit_length(N)<=n`. Divide out factors of two
and apply perfect-power preprocessing as in the quadratic note. For
each remaining current odd integer `M`, the canonical nontrivial square
roots of one are

```math
2\le z\le\lfloor M/2\rfloor,\qquad z^2=1\pmod M,
```

and their exact count is

```math
H(M)=2^{\omega(M)-1}-1.                                          \tag{8}
```

If `omega(M)=1`, this count is zero. Otherwise
`1<=omega(M)-1<=n`, since every distinct prime factor is at least two.
The selected prime satisfies `2^j!=1 (mod p)` for every exponent in
that range. Hence

```math
H(M)\ne0\pmod p\quad\Longleftrightarrow\quad\omega(M)\ge2.        \tag{9}
```

Once perfect powers are removed, this is exactly the required
prime/composite distinction. For example, merely requiring `p>n`
would not suffice: `N=3*5*7*11*13*17=255255` has 18 bits and
`H(N)=31`, which vanishes in characteristic 31. Its base-two order is
only five. Condition (1) excludes this accident.

All prefix counts satisfy `C(w)=C(w0)+C(w1)` in `F_p`. A nonzero parent
therefore has a nonzero child. Recover the scalar root, take
`gcd(M,z-1)`, and recurse on both factors. The same prime remains valid
throughout: every recursive integer has at most `n` bits and therefore
obeys the same order bound. There are `O(n^2)` queries and polynomial
ordinary arithmetic overhead.

Here `p` is the field characteristic for the *count*. The root equation
still uses the original integer modulus `M`. The auxiliary arithmetic
moduli inside the source convolution network are separate again. One
must not replace the integer congruence by its reduction modulo the
counting-field characteristic.

## 4. Uniform evaluator contract

Let `D(n)!=0` be generated with charged time and bit length, and suppose
an exact counting construction is available over every odd prime not
dividing `D(n)`, for all the prefix-pinned quadratic networks needed at
input lengths at most `n`. The same statement can use a padded common
size for all recursive queries.

Require the following explicitly:

1. The coefficient and local-table algorithms are uniform in `p`.
   Their stated bounds include all work that depends on `p`.
2. Field operations and zero tests use actual `O(log p)`-bit
   representations. Ordinary modular arithmetic costs
   `poly(log p)` bit operations; it is not a unit-cost oracle.
3. The evaluator's bound includes coefficient construction, the fused
   network or weighted-compiler construction, and all table generation.
4. The forbidden-prime certificate covers every recursive query, not
   only the top instance. A product of supplied certificates for all
   smaller sizes is legitimate only with its generation and bit costs
   included.

If `D(n)` has polynomial bit length and generation cost, then (1)
selects a prime of polynomial size, with `log p=O(log n)`, in polynomial
time. A uniform evaluator bound

```math
T_{\mathrm{eval}}(n,p)
\le\operatorname{poly}(n,\log p)\,\exp(O(F(n)))                    \tag{10}
```

therefore gives deterministic complete factoring in
`poly(n) exp(O(F(Cn)))`, including all construction costs in (10).
This preserves every fixed target `F(n)=n/(log n)^alpha` and every
fixed `F(n)=n^beta`, `0<beta<1`, whenever the required evaluator exists.
It does not introduce the earlier sparse-detector floor.

More generally, if `bitlength(D(n))` and its generation time are bounded
by `poly(n,A(n))`, selection costs `poly(n,A(n))`, and
`log p=O(log(n+A(n)))`. The resulting evaluator cost must still be
checked uniformly at that selected prime. This permits a larger
certificate when its entire cost already fits the proposed runtime;
it does not assert that an arbitrary varying-field seed admits a cheap
certificate or works over the chosen characteristic.

This result is especially useful for a family whose arithmetic formulas
work over all primes outside an effectively controlled exceptional set.
Exact characteristic-zero counting already had other routes to factor
recovery. The additional point here is a freely chosen, small, explicit
finite field while retaining the one-square verifier.

## Verification

`check_high_order_characteristic.py` checks the elementary LCM divisibility
and primorial bounds, literal prime selection with forbidden-prime
certificates, and all modular powers in each returned certificate. It
also runs complete scalar-root recovery over the selected odd fields
against independent trial division. Negative controls show that an
arbitrary small characteristic or the weaker condition `p>n` can lose
the count. These checks validate the finite transcription; the proof
above supplies the all-size theorem.

No Lean formalization or new factoring exponent is claimed.
