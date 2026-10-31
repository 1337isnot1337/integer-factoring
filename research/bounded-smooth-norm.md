# Smooth-pair extraction is a bounded represented-norm problem

Let Delta<0 be a discriminant and let 4R^2<|Delta|. Let S_R consist of
reduced forms with leading coefficient at most R, coprime to Delta and
Y-smooth. For a reduced target Q=(A,B,C), the following are equivalent:

1. Its class lies in S_R/S_R.
2. It primitively represents a positive Y-smooth integer n coprime to
   Delta, with a divisor d satisfying d<=R and n/d<=R.

In particular n<=R^2. A supplied representation Q(x,y)=n with this
balanced-divisor property can be converted deterministically to the pair.
The search for the representation remains open. This is a precise
arithmetic replacement for the earlier smooth-pair interface, not a new
complete factoring algorithm.

## From a pair to a represented norm

Use the standard ideal interpretation of primitive quadratic forms. If
two forms have leading coefficients a,b<=R, their quotient is represented
by the product of the first ideal and the conjugate of the second. Its
norm is ab. Remove its integer content k to obtain a primitive ideal;
its norm is n=ab/k^2, and k divides gcd(a,b). These assertions follow
prime by prime: away from Delta a primitive ideal contains only one of
the two conjugate prime-ideal orientations, and common opposite
orientations cancel as rational prime factors.

The resulting primitive ideal has the target class and corresponds to an
equivalent form of leading coefficient n. Hence Q represents n by a
primitive integer vector. It is Y-smooth and coprime to Delta, and the
divisors a/k and b/k are both at most R. No enumeration of ideals is used
in this existence proof.

## Converting a supplied representation to a pair

For primitive (x,y), choose u,v with xv-yu=1. The proper change of variables
with first column (x,y) gives an equivalent form with coefficients

```math
n=Ax^2+Bxy+Cy^2,\qquad
b'=2Axu+B(xv+yu)+2Cyv.
```

For each prime p dividing n, the residue b' modulo 2p supplies the required
square root of Delta modulo 4p. Its prime form is therefore known directly;
no modular square-root search is needed. Lenstra--Pomerance Section 2.8
gives the signed prime-form factorization of the target class. Equivalently,
use the prime form oriented by that residue, with positive exponent v_p(n).

Split these exponents according to d and n/d. Put the first portion in
the left ideal and the **inverse** of the second portion in the right
ideal. Their class quotient is Q. Their norms are d and n/d. Because
both norms are <=R<sqrt(|Delta|)/2, centering the middle coefficient already
gives reduced forms with those same leading coefficients: their trailing
coefficients exceed their leading coefficients. Thus both output forms
belong to S_R. This proves the converse and gives an explicit conversion.

A nonprimitive supplied vector can first be divided by gcd(x,y). The
implementation tests the resulting primitive norm and balanced-divisor
condition, and can reject if no split exists.

## Exact balancing and its complete bit cost

Smoothness alone does not imply balanced divisibility. For example,
Q=(2,-1,125), Delta=-999 represents 125 at (0,1); 125 is 7-smooth but
has no divisor giving two factors <=12. Conversely, greedy assignment of
prime factors can miss a valid split: 210 has split 14*15 although the
simple greedy rule produces 10*21. The real-form control
Q=(2,-1,210), Delta=-1679 checks this case.

Given the prime factorization of an integer n of at most b bits, enumerating
its divisors costs 2^O(b/log b) times a polynomial in b. An elementary
uniform bound suffices: split its prime factors at sqrt(b). There are at
most sqrt(b) smaller primes, each with exponent <=b, contributing at most
sqrt(b)*log2(b+1) to log2(tau(n)). For larger primes, e+1<=2^e and their
total multiplicity is at most 2b/log2(b). Therefore

```math
\log_2\tau(n)\le\sqrt b\log_2(b+1)+2b/\log_2 b
=O(b/\log b).
```

The finite small-b range is absorbed into the constant. Testing all divisors
for the interval ceil(n/R)<=d<=R gives the exact split with this cost.
Trial division through Y supplies the smooth factorization at polynomial
cost in Y and b. With efficient standard quadratic-form operations, the
complete conversion cost is

```math
Y\,\operatorname{poly}(b)+2^{O(b/\log b)}.
```

For smaller supplied norms with Y*n<=R^2, greedy balancing already works:
assign each prime factor to the smaller load. The load ratio stays <=Y,
so the larger load is <=sqrt(Y*n)<=R. This version needs no divisor search.

Thus the exact converter fits logarithmic-saving exponent **1**, whereas
the restricted-norm version fits the smaller arithmetic scale when Y does.
This describes one component's cost. It establishes no factoring score:
finding a represented balanced smooth norm, and deterministic coverage of
factoring directions, remain unresolved. In a conditional transfer the
exact divisor step would preserve any fixed target exponent below 1.

## Compatibility with the density input

Lenstra--Pomerance Theorem 8.1 and its proof on pp.503-504 count smooth
integers up to x=(1/2)sqrt(|Delta|/d), with 1<=d<=12, and use Lemma 2.10
to produce reduced forms. Thus their counted subset already lies in the
small-leading-coefficient range used here, under the source's stated
prime-density hypotheses. No additional factor of sqrt(Y) is lost in the
radius when using exact balanced-divisor extraction. Integer boundary
cases are handled by the strict radius and coprimality condition.

## Finite controls and the remaining search

Run `python3 tools/check_bounded_smooth_norm.py`. Across eight discriminants
the [record](../evidence/bounded-smooth-norm.json) checks 2,500 conversions,
778 scaled nonprimitive controls and 56 rejections of unbalanced smooth
norms. It reaches 91 classes outside the original reduced-smooth set.
For each discriminant the set reached by exhaustive bounded-norm search
equals the independently computed full quotient set S_R/S_R. The 210
control specifically exercises exact divisor search after greedy fails.

The [implementation](../tools/bounded_smooth_norm.py) uses the small exact
quadratic-form package for controls. That package's simple composition
search is not the efficient class-group algorithm used in the cost theorem.
The theorem uses standard polynomial-bit-cost composition; the finite
implementation is not being advertised as that full asymptotic machine.

Exhaustive witness search is still expensive. Completing the square gives

```math
(2Ax+By)^2+|\Delta|y^2=4An\le4AR^2.
```

This provides exact finite coordinate bounds but can leave exponentially
many candidates in b=bits(|Delta|). The open arithmetic task is to locate
a balanced smooth represented norm within the guaranteed sampled target
family without scanning this whole region. No faster search is supplied
by the conversion theorem or its checker. Substantial OpenAI Codex
assistance was used; no new Lean formalization is added.
