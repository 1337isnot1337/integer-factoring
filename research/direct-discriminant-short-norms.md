# Single-form smooth-norm search at discriminant -4N

## Result and scope

Two simplifications remove auxiliary tasks from the arithmetic route.

First, the prime sample can be absorbed into the represented norm.
The remaining search may operate on one target form, rather than return
a sample element and two bounded reduced forms. Second, the discriminant
can always be chosen as `Delta=-4N` for an odd integer with at least two
distinct prime factors. A full basis of the two-torsion at that one
discriminant suffices to split `N`.

An explicit two-dimensional ideal-HNF implementation also converts a
supplied prime-form word to an actual primitive represented norm and
performs form composition with polynomial bit cost. It avoids the
coprime-leading-coefficient search in the older small-control package.

These are ordinary proofs and exact constructive controls. Finding the
required short smooth represented norms within the desired uniform time
bound remains open. No stronger complete factoring exponent is claimed.

The starting revision is
`a5e05feb756b4ad35ae1877d54a5f2e4eb387c19`. This note builds on the
[positive relation kernel](positive-relation-kernel.md),
[prime-form generator and sampling theorem](polylog-prime-form-generators.md),
[direct smooth-norm words](smooth-norm-short-words.md), and
[two-part extraction algorithm](arithmetic-torsion-extraction.md).

## 1. Absorb the sampled prime into the norm

Let `S_R` consist of reduced forms with `Y`-smooth leading coefficient
at most `R`, coprime to `Delta`, where `4R^2<|Delta|`. Suppose its
density in the full class group is at least `delta>0`. Let a distribution
on oriented prime forms of prime norm at most `X`, coprime to `Delta`,
have character bias less than `delta/(1-delta)` when `delta<1`.

For every class `t`, there exists a primitive representation of an
integer `n` by the reduced form of `t` such that

```math
1\le n\le XR^2,\qquad \gcd(n,\Delta)=1,
\qquad n\text{ is }Z\text{-smooth},
\qquad \Omega(n)\le\lfloor\log_2(R^2)\rfloor+1,
\qquad Z=\max(X,Y).
```

### Proof

The positive-kernel lemma gives `u` in the sample and `a,b` in `S_R`
such that `t*u=a*b^-1`. Choose the integral primitive ideals represented
by these three forms. Multiplying the ideal for `a`, the conjugate ideal
for `b`, and the conjugate prime ideal for `u` gives an integral ideal
of class `t` and norm at most `XR^2`. Removing its rational-integer
content divides this norm by an integer square. The resulting ideal is
primitive, is still coprime to the conductor, and has the same class.

All its rational prime factors are at most `Z`. Its number of prime
factors, counted with multiplicity, is at most the sum for the two
leading coefficients and the one prime norm. This is at most
`floor(log2(R^2))+1`. The primitive ideal/form correspondence supplies
the asserted primitive representation. The HNF procedure below makes
this last step constructive when the three inputs are supplied.

### Single-target relation rows

Set

```math
L=\lfloor\log_2(R^2)\rfloor+1,\qquad K=L+1.
```

Include all the needed oriented prime forms of norm at most `Z` in the
generator list `g_1,...,g_m`. For each `i`, search the single reduced
form of `g_i^K` for a primitive `Z`-smooth norm of multiplicity at most
`L`. The displayed existence theorem guarantees one with `n<=XR^2`.
The search may also accept a larger norm if it satisfies the same
multiplicity budget; then `n<=Z^L` still bounds its bit length.

The [direct converter](smooth-norm-short-words.md) returns a word `v_i`
for that class with `||v_i||_1<=L`. Hence

```math
r_i=K e_i-v_i
```

is a strictly diagonally dominant relation row. The full matrix has
nonzero determinant of absolute value at most `(2K-1)^m`. Neither
the sample element nor a balanced divisor split is part of this search
output. The prime-sampling theorem is used to prove existence; the
algorithm need not materialize the sample distribution for this version
of the contract.

## 2. The fixed discriminant -4N covers every factoring direction

Let `N>1` be odd. Work with primitive reduced forms of discriminant
`-4N`, with the canonical convention

```math
-a<b\le a\le c,\qquad b\ge0\text{ if }a=c.
```

The inverse of `(a,b,c)` is the reduced form equivalent to `(a,-b,c)`.
It equals the original form exactly in the three cases
`b=0`, `b=a`, or `a=c`. Interior forms have distinct canonical inverses;
the indicated boundary cases are fixed by a proper swap or shear.
Thus these are precisely the two-torsion classes.

Their possible failures to give a proper factor can be classified
without knowing the factorization of `N`.

* If `b=0`, then `N=a*c` and `gcd(a,c)=1`. The form gives a proper
  factor unless `a=1`, the principal form `(1,0,N)`.
* If `b=a`, write `a=2u`. Primitiveness and odd `N` imply that `c` and
  `u` are odd and `gcd(u,c)=1`. Then
  `N=u*(2c-u)` is a coprime factorization. Since `c>=2u`, the only
  failure is `u=1`, giving `(2,2,(N+1)/2)`. This form exists exactly
  when `N` is `1` modulo `4`.
* If `a=c`, write `b=2v`. Then `a` is odd, `v` is even, and
  `N=(a-v)*(a+v)` is a coprime factorization. Since `0<=v<=a/2`, it
  is proper for every `N>1` in this case.

The classes that fail to split `N` therefore form either the trivial
subgroup or the order-two subgroup consisting of the principal class
and `(2,2,(N+1)/2)`. If `N` has at least two distinct prime factors,
some coprime pair `N=u*v` has both factors greater than one. The
primitive reduced form `(min(u,v),0,max(u,v))` is a two-torsion class
outside that nonfactoring subgroup.

A basis of the full two-torsion cannot lie entirely inside a proper
subgroup. Consequently at least one basis element gives a proper
factor by the displayed formulas and a gcd. This proves the required
coverage at `Delta=-4N`, including nonfundamental discriminants and
integers with repeated prime factors. No multiplier selection or
random ambiguous class is required.

The algorithm removes even factors and handles primes and perfect
powers separately with deterministic polynomial-bit arithmetic and
primality testing. These preprocessing steps do not require a
factorization of the discriminant or conductor.

## 3. Exact two-rank at this discriminant

Write `omega(N)` for the number of distinct prime factors. For odd
`N>1`, the two-rank is

```math
\dim_{\mathbf F_2} C_{-4N}[2]=
\begin{cases}
\omega(N)-1,&N\equiv3\pmod4,\\
\omega(N),&N\equiv1\pmod4.
\end{cases}
```

Here is a direct count, specialized to the required discriminant.
There are `2^(omega(N)-1)` unordered coprime factorizations `N=u*v`.
Each gives the distinct reduced form `(u,0,v)` after ordering `u<=v`.

The nonzero-middle ambiguous forms exist only for `N=1 mod 4`.
For each such coprime pair, exactly one further reduced form is

```math
\left(2u,2u,\frac{u+v}{2}\right)\quad\text{if }v>3u,
\qquad
\left(\frac{u+v}{2},v-u,\frac{u+v}{2}\right)
\quad\text{if }v<3u.
```

The equality case would force `u=1,v=3`, incompatible with `N=1 mod 4`.
Both forms are primitive because `u,v` are coprime and have the same
residue modulo `4`. Conversely, the factorizations in the preceding
section recover the unique pair from every nonzero-middle ambiguous
form. This doubles the count exactly when `N=1 mod 4`, proving the
rank formula.

This formula is used for analysis, not computed by factoring `N`.
After trial division through `2^ceil(sqrt(b))`, where `b=bits(N)`, a
remaining cofactor has `omega(N)<=sqrt(b)+O(1)`. The existing generic
two-part extraction then costs `2^O(sqrt(b))` times polynomial input
costs. Without this trial division, the elementary bound
`omega(N)<=sqrt(b)+2b/log2(b)` gives `2^O(b/log b)` extraction instead.

## 4. A polynomial-bit ideal implementation

The older finite form package searches for a coprime equivalent leading
coefficient. The new [HNF implementation](../tools/quadratic_ideal_hnf.py)
uses a fixed number of integer gcd operations instead.

Put `d=Delta mod 2` and `omega=(d+sqrt(Delta))/2`. Then
`omega^2=d*omega+(Delta-d)/4`. The ideal for `(a,b,c)` has basis

```math
(a,0),\qquad\left(-\frac{b+d}{2},1\right)
```

in coordinates relative to `(1,omega)`. Multiply two bases to obtain
four integer vectors. If their second coordinates have gcd `h`, the
extended Euclidean algorithm gives a lattice vector `(q,h)`. Let `D`
be the gcd of all their two-by-two minors. The HNF basis is

```math
(D/h,0),\qquad(q\bmod(D/h),h).
```

The determinant-gcd formula is the two-dimensional lattice index.
Closure under multiplication by `omega` implies that `h` divides both
first coordinates. Dividing the whole basis by `h` removes the largest
rational-integer content and leaves a primitive ideal with basis
`(n,0),(r,1)`. Its form is

```math
\left(n,-2r-d,\frac{(-2r-d)^2-\Delta}{4n}\right).
```

Ordinary Gauss reduction supplies its canonical form. The implementation
retains the product `S` of the proper integer shears and swaps, so if
the reduced form is `Q(z)=F(Sz)`, the primitive vector
`S^-1*(1,0)` explicitly satisfies `Q(x,y)=n`.

For a supplied oriented prime word of product norm `M`, ideal
multiplication and content removal give `n<=M`, with `M/n` an integer
square. All intermediate integers have `O(bits(M)+b)` bits. Fixed-size
HNF, binary powering, Gauss reduction, and the retained coordinate
transformation have polynomial bit cost in that written size. A product
of already-reduced forms is therefore computed in polynomial time in
`b`. No root search, coprime-leading search, or class-group enumeration
occurs in these routines.

This is standard ideal arithmetic made explicit for the current
constructors. It does not solve the search for a short prime word.

## 5. Reproduction and remaining obligation

Run:

```sh
python3 tools/check_direct_discriminant_short_norms.py
```

The [checker](../tools/check_direct_discriminant_short_norms.py) compares
the HNF operations with the independent older form implementation,
checks word-to-norm representations and the prime-absorption bound,
and checks the complete ambiguous-form classification at finite
discriminants. Its [record](../evidence/direct-discriminant-short-norms.json)
states the actual ranges and finite enumeration costs.

The open task is still to find a primitive smooth norm with the stated
multiplicity bound, for each required target form, in a uniformly bounded
time better than the current factoring bound. The existence proof and
the polynomial converter do not charge an exhaustive search as though
it were fast. Substantial OpenAI Codex assistance was used. No new Lean
formalization is added.
