# Exact divisor counts in polynomial space

## Theorem

There is a deterministic algorithm which, for an `n`-bit integer `N>=2`
and a binary prefix, computes the exact number of divisors with that
prefix in

```math
\exp\!\left(O\!\left(\frac{n}{\sqrt{\log n}}\right)\right)
\quad\text{time and}\quad n^{O(1)}\text{ bit space}.
\tag{1}
```

The least prime factor follows from `n` such queries. Complete prime
factorization adds a polynomial number of calls. This gives the same
half-rate factoring bound in polynomial space, with **polynomial
recovery overhead**. No preliminary trial division or witness-isolation
transform is needed for this integral family.

The proof uses translated integral Schur factors, queried without
partition-sized caches, and evaluates the contraction in the ring
`Z/2^(n+1)Z`. All crossing and contraction identities are integral and
division free after their factor entries have been generated. The
coefficient ring need not be a field.

This result strengthens the exact-output and space interface. The
time exponent remains `1/2`. The separate
[bit-tagged recovery theorem](bit-tagged-witness-recovery.md) remains
useful for future families that exist only in a fixed characteristic.

## 1. Stream the translated dual coefficients

The integral family in [integral-schur-crossing.md](integral-schur-crossing.md)
uses `b_x=s_x` and `a_x=T omega(s_x)`, where `T=exp(D)` and

```math
T(p_m)=\sum_{k=1}^m\binom mkp_k.
\tag{2}
```

Its unit vacuum pairing, crossing sum, and both channel factorizations
hold over the integers at every positive weight. Prime-power weights
are unnecessary for this translated version.

For a power-sum parameter `p_mu` and a target monomial `p_nu`, the
translation coefficient is zero unless `ell(mu)=ell(nu)` and
`|nu|<=|mu|`. Otherwise it is the coefficient of `X^m(nu)` in

```math
\prod_{a\in\mu}\left(\sum_{k\in\operatorname{supp}\nu}
                         \binom akX_k\right).
\tag{3}
```

This polynomial is homogeneous in its ordinary degree. Apply the
forward difference of order `m_k(nu)` in each variable at zero, then
divide the integer numerator by `prod_k m_k(nu)!`. Only the desired
monomial survives. The division is exact over the integers and is
completed before any ring reduction.

There are `prod_k(m_k(nu)+1)` grid points. These encode subpartitions
of `nu`, so their number is at most
`sum_(j<=|nu|)p(j)=exp(O(sqrt(|nu|)))`. At each point, evaluate (3)
as a product of linear forms. All binomials, products, and accumulators
have polynomial-in-`|mu|` bit length. The grid is streamed.

Thus a translation coefficient takes `exp(O(sqrt(h)))` time and
`poly(h)` space when the degrees are at most `h`. This fills the
remaining query step in the integral translated family.

## 2. Complete integral half-factor queries

The [streamed Schur algorithms](prime-power-gram-involution.md) compute
one character and one source power-sum mode coefficient in
`exp(O(sqrt(W)))` time and `poly(W)` space for degrees at most `W`.
Their algorithms use exact integer or rational arithmetic before
reducing an integral answer.

Combine these with (3) and the character change of basis

```math
s_x=\sum_{\mu\vdash h}\frac{\chi^x_\mu}{z_\mu}p_\mu,
\qquad p_\rho=\sum_u\chi^u_\rho s_u.
\tag{4}
```

For example, a translated angle entry on an input of degree `E` is
the sum over `mu` of degree `h`, lowered `nu` of degree `w<=h`,
input `kappa` of degree `E`, and output `rho` of the prescribed degree,
of

```math
\frac{(-1)^{h-\ell(\mu)}\chi^x_\mu\chi^y_\kappa}
     {z_\mu z_\kappa}
[p_\nu]T(p_\mu)
\binom{h+w}{i}
[p_\rho](p_\nu)[r]p_\kappa\,
\chi^u_\rho.
\tag{5}
```

The sign is `omega`, the binomial is the source's angle coefficient,
and the outside crossing sign `(-1)^i` is applied separately.
There is a fixed number of nested partition enumerations, all of
degrees at most `4h`. Their total count remains `exp(O(sqrt(h)))`.
A common denominator in (5) divides `h! E!`. Every individual term
and the running numerator have polynomially many bits. Integrality of
the source mode on the Schur lattice proves that the final division
is exact. The production algorithm stores no partition-indexed table.

The unshifted `b_x` factor uses the same formula without translation
and angle. These are the four factor-query interfaces required by
the two source operator orders. Their channel counts are

```math
R_s=\sum_{i=0}^{2h}p(s+i),\qquad
C_s=\sum_{i=0}^{2h}p(2h-s-i).
\tag{6}
```

Choose the least `h>=2` with `p(h)>=2^lambda`. The integral-Schur
analysis gives `h=O(lambda^2)`, polynomially many pieces, and

```math
\lceil\log_2R_s\rceil+\lceil\log_2C_s\rceil
\le2\sqrt2\lambda+O(\log(\lambda+2)).
\tag{7}
```

Cap either separate width at `2lambda` if necessary by factoring
through that flattening's own row indices, using an identity matrix.
The other factor is a complete tensor entry. Compute it by streaming
one of the factorizations in (6), which still takes `2^O(lambda)`
time and polynomial space. This cap uses no Gaussian elimination and
no field inversion.

Partition counts and rank/unrank metadata use a polynomial-size
dynamic-programming table indexed by integer weights and bounds,
not a table indexed by all partitions. Hence every factor, complete
entry, piece descriptor, and channel index has a uniform
`2^O(lambda)`-time, `poly(lambda)`-space algorithm over the integers.
Reduction modulo an arbitrary supplied integer `M>=2` adds only a
polynomial cost in `log M`.

## 3. Why the evaluator works over a ring

Use the explicit top queries and grouping from
[polynomial-top-queries.md](polynomial-top-queries.md), or the
ordinary constraint representation of
[logarithmic-wire-network.md](logarithmic-wire-network.md).
The latter may choose its arithmetic-prime parameters with the fixed
auxiliary characteristic parameter `p=2`; this is separate from the
coefficient ring used in this proof. Its Fourier representation is
not required here.

For a prefix query, these top tensors have entries zero or one. Their
integer contraction equals the number of divisors with that prefix.
Unique extension of each legal gate setting establishes that equality.

Replace each crossing by the integral scale-one identity, choose the
orientation from the actual integer channel widths, group vertices,
and recursively contract. These steps use only addition and
multiplication of coefficients. The grouping search uses integer
width data. It makes no field-rank test. The terminal enumeration is
also a sum of products. Since the crossing multiplier is one, no
scalar division occurs in the contraction.

Every identity therefore remains valid in `Z/MZ`, including when
`M` is composite and has nonzero zero divisors. The separate integer
coefficient generator in Section 2 finishes its exact divisions
before producing a ring element.

The volume and branching proof in
[subexponential-piece-transfer.md](subexponential-piece-transfer.md)
now applies to these exact operations at `rho=1/sqrt(2)`. Its estimates
are uniform in the coefficient ring once its operations are charged.
With `log M=poly(n)`, each ring operation has polynomial bit cost.
This changes only the polynomial overhead in the same recursion bound.
Consequently the time is (1).

Depth-first execution stores one branch assignment, channel counter,
and residue accumulator at a time. The recursion has `O(log log n)`
levels and polynomial-size descriptors. Each local query uses
polynomial space by Section 2, and each residue has `O(log M)` bits.
Parent queries are recomputed when needed; the existing sum of query
costs already charges that recomputation. This proves polynomial space.

## 4. Recover the exact prefix count and least factor

Fix

```math
M=2^{n+1},\qquad
D(N,k,x)=\#\{d\mid N:\lfloor d/2^k\rfloor=x\}.
```

The evaluator returns the canonical residue of `D(N,k,x)` modulo `M`.
But `0<=D(N,k,x)<=tau(N)<=N<2^n<M`, so that residue is the exact
nonnegative integer. No nonzero residue assumption, sparse polynomial
test, prime supply, or interpolation of the final count is needed.

Start with prefix `a=0`. For `k=n-1,...,0`, put `b=2a` and query

```math
c=D(N,k,b)-[\lfloor1/2^k\rfloor=b].
```

If `c>0`, set `a=b`; otherwise set `a=b+1`. The invariant is that
`a` is the smallest prefix attained by a divisor greater than one.
At termination it is the least such divisor, hence prime. This is
also correct when `N` itself is prime. There are exactly `n` queries,
and all additional recovery work is polynomial in `n`.

Repeatedly divide out the returned prime and recurse on the remaining
quotient. There are at most `n` prime factors with multiplicity. Using
the original input length as an upper bound on all work gives the
complete factorization within (1) and polynomial space.

## 5. Scope and verification

This theorem uses an integral, scale-one family with explicit
polynomial-space factors. It does not upgrade an arbitrary
characteristic-`p` rank factorization to a ring factorization. A future
family available only over `F_p` still needs the separate recovery
theorem and its hypotheses.

The algorithms are a proof of a uniform asymptotic bound. Their large
constants and direct partition enumerations are not presented as a
practical factorer for large integers.

Run:

```sh
python3 tools/check_integral_counting_transfer.py --output evidence/integral-counting-transfer.json
```

The checker compares 1,230 streamed translation coefficients with an
independent expansion, and 3,312 integral Schur factor entries with
the separate source recursion and Jacobi--Trudi implementation.
It verifies both operator orders and the crossing sum, including
745 piece comparisons in the composite rings of moduli 4, 8, 9,
16, and 25. Its exact prefix oracle recovers the least factor of
every integer from 2 through 4,096, using 45,069 prefix queries.

These finite tests check the coefficient algorithms, ring identities,
and recovery invariant. They do not benchmark the full asymptotic
recursive machine at large input sizes. The proof above supplies
the runtime, space, and uniformity statements.
