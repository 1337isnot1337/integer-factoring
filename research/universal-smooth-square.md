# One fixed smooth square covers every square class in the factor base

Dense small smooth forms imply a fixed-norm existence theorem. A single
explicit smooth square is represented, possibly nonprimitively, by every
class in the square subgroup of the factor base. Its encoding length stays
at the arithmetic subexponential scale. The cap is chosen before the
relation targets, and every representation at that cap gives a short word.
Thus an algorithm finding **any** such representation would give a full-rank
relation matrix. No representation search algorithm with the required cost
is proved here.

The proved advance is a sharper conditional contract. There is no need to
search a list of candidate norm values or use a target-dependent biased
sample in this contract. Finding a representation of the one fixed norm is
exactly a bounded word problem in the squared prime-form generators. The
normal form makes that remaining search obligation explicit; it does not
solve it.

## 1. Assumptions and the fixed cap

Let G be the proper class group of primitive positive definite binary
quadratic forms of discriminant Delta < 0. Fix integers Y and R with

```math
2\le Y\le R,\qquad 4R^2<|\Delta|.
```

For each prime p <= Y coprime to Delta for which Delta is a square modulo
4p, choose the root b_p in [0,p] and the prime-form class g_p. Let P be the
set of these primes, m its size, and H the subgroup they generate. Define
S to be the reduced classes whose leading coefficient is at most R,
Y-smooth, and coprime to Delta.

Every class in S has an explicitly recoverable word in the g_p of length
at most

```math
L=\lfloor\log_2R\rfloor.
```

This is the elementary oriented prime-ideal factorization used in
[smooth-norm relations](smooth-norm-relations.md). The radius condition
ensures that each chosen prime form is already reduced and belongs to S.
The identity belongs to S as well. In particular, S is a subset of H and
generates H.

Assume a certified positive integer Q satisfies

```math
|G|\le Q|S|.
```

The following parameters depend on this density data and the factor base,
before any relation target is chosen:

```math
T=2(2Q-1)L,\qquad
M=\prod_{p\in P}p^T,\qquad
K=Tm+1.
```

The trivial subgroup case m=0 needs no relation queries. In the remaining
statements the fixed represented norm is M^2, and a relation target is
g_p^(2K).

## 2. A short-word covering theorem

Put A=S/S. The Cayley graph on H with connection set A is connected:
A contains S because the identity belongs to S, and S generates H. A is
symmetric and contains the identity.

Take a shortest path v_0,...,v_d. The translates

```math
v_0S,\ v_2S,\ v_4S,\ldots
```

are pairwise disjoint. If v_i s=v_j t with s,t in S and |i-j| >= 2,
then v_i^(-1)v_j=s/t lies in A. That would shorten the path, contradicting
the fact that every subpath of a shortest path is shortest. Consequently

```math
(\lfloor d/2\rfloor+1)|S|\le |H|\le |G|\le Q|S|,
\qquad d<2Q.
```

Each A-step has a factor-base word of length at most 2L. Every h in H
therefore has some integer exponent vector e with

```math
h=\prod_{p\in P}g_p^{e_p},\qquad
\sum_p|e_p|\le 2(2Q-1)L=T.
```

In particular |e_p| <= T for every p. This is an existence proof. Following
the graph by exhaustive enumeration, as the finite checker does, can cost
as much as enumerating H. Density alone does not make those paths
algorithmically accessible.

The same conclusion can be expressed directly in H^2 with the squared
generators: squaring a short word gives a word of length at most T in the
g_p^2 for every element of H^2.

## 3. The universal square and both directions of its word equivalence

**Theorem.** Every primitive positive form whose class belongs to H^2
represents M^2 over the integers. Conversely, every class representing
M^2 belongs to H^2. More precisely, representing M^2 is equivalent to
having a word

```math
[F]=\prod_{p\in P}g_p^{2e_p},\qquad |e_p|\le T.
```

For the forward construction, take such a word. For each p choose the
oriented prime ideal corresponding to g_p if e_p is positive, and its
conjugate if e_p is negative. Give that ideal exponent 2|e_p|. The resulting
primitive invertible ideal has class [F] and norm

```math
n_0=\prod_p p^{2|e_p|}.
```

Only one of the two primes above p occurs, so the ideal has no rational
content. The associated primitive norm form has leading coefficient n_0.
Proper equivalence to F gives a primitive vector (x_0,y_0) with
F(x_0,y_0)=n_0. Scaling by

```math
k=\prod_p p^{T-|e_p|}
```

gives F(kx_0,ky_0)=M^2. The short-word theorem supplies the required word
for every element of H^2.

For the reverse direction, suppose F(x,y)=M^2 and put
k=gcd(x,y)>0. Since k^2 divides M^2, k divides M. The primitive norm is

```math
n_0=F(x/k,y/k)=(M/k)^2.
```

Every prime divisor belongs to P, every valuation is even, and every
valuation is at most 2T. Primitive oriented ideal factorization produces
signed exponents 2e_p with |e_p| <= T and the stated class identity.
This applies to **every** returned representation, not only the witnesses
from the forward proof.

Both conversions use polynomially many bit operations in log|Delta|,
log M, and Y. The forward direction is explicit: lift the selected roots
modulo prime powers, combine them by the Chinese remainder theorem, and
reduce the form with leading coefficient n_0 while retaining its
unimodular change of variables. The arithmetic and its complete cost are
given in [factored norm search](factored-norm-search.md). The reverse
direction uses a coordinate gcd and the prime-word conversion from
[smooth-norm relations](smooth-norm-relations.md).

Thus the fixed-square oracle is polynomially equivalent to this bounded
word search in the squared generators. This equivalence prevents an
existence certificate from being counted as an efficient arithmetic
search algorithm.

### Why the square subgroup is necessary

For a general fixed smooth norm n=product p^(a_p), normalize any represented
vector by its coordinate gcd. The remaining valuation at p is congruent
to a_p modulo two. The sign of its prime-ideal orientation does not change
that parity. Every class representing n therefore lies in the single coset

```math
\left(\prod_p g_p^{a_p}\right)H^2.
```

No fixed smooth norm can cover all of H when H/H^2 is nontrivial. For
n=M^2 this coset is exactly H^2, and the theorem proves coverage of the
entire coset. Allowing nonprimitive vectors is essential: primitive
representations of M^2 permit only the extreme exponents +2T or -2T at
each prime, whereas nonprimitive vectors permit every even exponent in
[-2T,2T].

## 4. Any exact representation gives a dominant relation

For each factor-base generator g_i, compute its reduced target class

```math
F_i=g_i^{2K},\qquad K=Tm+1.
```

The target belongs to H^2, so it represents the **same** fixed norm M^2.
Suppose a proposed oracle returns any coordinates x_i,y_i with
F_i(x_i,y_i)=M^2. Normalize and convert them to the word v_i. By the
preceding theorem,

```math
[F_i]=\prod_p g_p^{(v_i)_p},\qquad
|(v_i)_p|\le 2T,\qquad
\|v_i\|_1\le 2Tm.
```

The integer row

```math
r_i=2K\mathbf e_i-v_i
```

is a relation on H, and its diagonal dominance slack is at least

```math
(r_i)_i-\sum_{j\ne i}|(r_i)_j|
\ge 2K-\|v_i\|_1\ge2.
```

All m rows therefore form a nonsingular integer matrix B. For completeness,
if Bz=0 with z nonzero, choose a coordinate of z with maximum absolute
value. Strict row diagonal dominance in that row gives a contradiction.
Since the row lattice is contained in the relation lattice of the
surjection Z^m -> H, |H| divides |det B|. Multiplying by lcm(1,...,Q)
gives a multiple of |G|, by the
[density index transfer](arithmetic-density-lift.md).

This construction cannot be satisfied merely by padding the given long
word. That word is 2K times the i-th basis vector. Its p_i valuation 2K
is greater than 2T, so it cannot fit in the fixed norm. Increasing T after
seeing the target would invalidate the contract. Here T is fixed from the
density data, and the target exponent is chosen afterward to exceed the
word budget forced by M^2. In particular none of the relation rows is zero.

The rows are even. Dividing a row by two produces an element of order at
most two, but need not produce a relation on H; no such division is used
in the determinant argument.

## 5. Bit cost and the remaining gates

Write b=ceil(log_2(|Delta|+1)). We have L=O(b), m <= Y, and

```math
\log_2 M
 =T\sum_{p\in P}\log_2p
 =O(QbY\log(Y+1)).
```

The prime list, splitting roots, parameters, targets, forward certificates,
reverse words, determinant, and index correction all have algorithms
polynomial in b, Y, and Q. If both Y and Q are
exp(O(sqrt(b log b))), these costs have that same asymptotic scale, with
changed constants. The bound is on the **bit length** of M. Its numerical
value is much larger.

An oracle with polynomial cost in b, Y, Q, and log M for these fixed-square
queries would consequently give a global class-number multiple at that
arithmetic scale. No such oracle is established. The norm is usually much
larger than the earlier R^2 cap; its larger encoding length is accounted
for above. Substituting this norm into the exact primitive orientation
search from the previous checkpoint does not give the desired bound.
Even deciding all nonprimitive orientations by direct search can require
searching up to (2T+1)^m exponent choices.

Two obligations remain distinct:

- Find any exact representation of the fixed square for the specified
  targets with the required uniform deterministic cost. Equivalently,
  find a vector in the certified exponent box for those squared classes.
- Supply the global ambiguous directions needed for factoring. A
  class-number multiple and the small-prime subgroup do not create
  directions absent from that subgroup.

The control Delta=-348, Y=7 still has |G|=6 and |H|=3. This theorem covers
all three square classes of H and supplies relation determinants and
global class-number multiples, while H has trivial two-torsion. The
missing ambiguous class (3,0,29) remains outside H. Thus the fixed-square
theorem does not close the final splitting gate or improve the proved
factoring score.

## 6. Exact finite certificates

Run `python3 tools/check_universal_smooth_square.py`. The
[evidence](../evidence/universal-smooth-square.json) enumerates small class
groups, their factor-base subgroups, the set S/S, and shortest paths. It
checks the diameter bound, builds the one fixed square before the target
exponents, constructs a representation for every square class, recovers
every exponent vector, and verifies all dominant relation rows and the
corrected class-number multiples.

The forward constructor in
[`universal_smooth_square.py`](../tools/universal_smooth_square.py) requires
a supplied bounded word. The finite checker supplies that word by
exhaustive group enumeration. This is a proof certificate and finite
control, not an implementation of the missing efficient oracle.

The recorded run checks 255 discriminants, 1,258 universal-square
representations, and 309 relation rows. The largest fixed norm has 1,599
bits. There are 66 trivial-subgroup controls and 68 controls in which H
misses global two-torsion. The largest checked S/S Cayley diameter is six,
at Delta=-43228; hence the covering checks include paths beyond one
quotient step.

This work used substantial OpenAI Codex assistance and adds no Lean code.
