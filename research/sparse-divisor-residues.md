# Global divisor residues with a chosen trial budget

The [fixed-characteristic transfer](fixed-characteristic-transfer.md)
can be sharpened by charging a chosen amount of trial division first.
If trial division reaches 2^t without finding a factor, the remaining
divisor polynomial has at most D=2^ceil(b/t) terms. Its complete support
can then be recovered from global cyclic residue polynomials, with
polynomial cost and algebra degrees in b,D. Prefix queries are not
required for this version.

Consequently a uniform weighted evaluator at saving exponent beta gives
a complete factoring algorithm with the same beta, for **every fixed
positive rational beta**. At target scores through one, the required
cyclic-algebra degrees and query count are polynomial in b. This removes
the cap at one from the earlier transfer's overhead. A faster evaluator
remains unproved, and the established score remains `1/2`.

The support recovery is the standard Garg--Schost sparse interpolation
method, using the integer polynomial whose roots are the unknown
exponents. See [*Interpolation of polynomials given by straight-line
programs*, March 2009 preprint](https://cs.uwaterloo.ca/~eschost/publications/interp.pdf),
Section 4, Steps 0--5, pages 5--7. This note specializes that method to
divisors and supplies the trial-budget transfer, a global-query contract,
and exact implementation. It does not claim a new general sparse
interpolation algorithm. The existing weighted constructor and its source
attribution are retained.

## 1. The exact global query and recovery parameters

Fix one prime characteristic p. For every prime q the oracle returns
the canonical q-coefficient representative of

```math
A_q(Z)=\sum_{d\mid N}Z^{d\bmod q}
\quad\text{in }\mathbb F_p[Z]/(Z^q-1).                  \tag{1}
```

It must return all coefficients, with complete cost charged. The
[weighted constructor](weighted-prefix-constructor.md) supplies exactly
this target by taking the full prefix query k=b,x=0. It uses the same
network grouping and widths. The unknown faster evaluation algorithm
is not replaced by enumerating the divisors.

Assume D>=tau(N). Choose

```math
B=bD(D-1)/2,\qquad L=B+D(b+1)+3.                        \tag{2}
```

Use the first L primes. Their maximum is at most 16(L+1)^2 by the
elementary central-binomial argument already proved in the transfer
note. Finding them and writing any one residue vector have polynomial
cost in b,D. This allows q=p and does not require the quotient algebra
to be reduced.

## 2. Recovering integer exponents in a fixed characteristic

Put r=tau(N). A collision prime divides

```math
\prod_{d<e,\ d,e\mid N}(e-d),
```

whose bit length is at most B. Support size is at most r, with equality
exactly without collisions, including in small characteristic. Thus the
maximum observed size is r, identifying at least D(b+1)+3 good primes.

For those primes form

```math
C_q(X)=\prod_{a\in\operatorname{supp}(A_q)}(X-a)
       \pmod q.
```

They reduce the integer polynomial

```math
C_N(X)=\prod_{d\mid N}(X-d).                            \tag{3}
```

Its coefficient height is at most (N+1)^D<=2^(bD). The product of good
primes exceeds twice that bound, so signed coefficientwise CRT recovers
(3). Residue **positions** are integers modulo q; the coefficient
characteristic p stays fixed.

Remove X-1. The least remaining integer root is the smallest prime
factor; [exact Sturm isolation](divisor-moment-recovery.md) recovers it
with polynomial bit cost from the O(Db)-bit coefficients.

There are L=poly(b,D) queries and all q are poly(b,D). Residue inspection,
support polynomials, CRT and root isolation also have poly(b,D) bit cost.
The implementation retains the current largest support and resets its
CRT state when that size increases. Once enough CRT modulus has been
collected it still inspects later supports, so a larger support cannot
be missed. This changes storage and arithmetic only; the full L-query
bound remains charged.

## 3. Trial division makes the ring degree small

Choose 1<=t<=b and test integers through min(N,2^t). A returned first
divisor is already the smallest prime factor. If none occurs, every
prime divisor exceeds 2^t, and the elementary proof from the
[moment note](divisor-moment-recovery.md) gives

```math
\tau(N)\le D=2^{\lceil b/t\rceil}.
```

Suppose all queries (1) can be computed by a uniform deterministic
algorithm in F(b)poly(b,q) bit operations, with F(b)>=1 and fixed polynomial degrees.
The complete factoring cost is then

```math
2^t\operatorname{poly}(b)+F(b)\operatorname{poly}(b,2^{\lceil b/t\rceil}).
                                                               \tag{4}
```

This is a proved conditional transfer, with no class-group generation,
integer-count reconstruction, primitive-root construction or factoring
of a finite-field order required.

For a fixed positive rational beta, let ell=ceil(log2(b+2)) and choose

```math
t=\left\lceil\frac b{\ell^\beta}\right\rceil,
\qquad D\le2^{\lceil\ell^\beta\rceil}.
```

The exact rational-power selection rule is given in the moment note.
If F(b)=2^O(b/log^beta(b+2)), the polynomial factors in (4) have logarithm
O(log b+log^beta(b+2)), and fit within the same bound. For beta<1 they
are polynomial in b with D=b^o(1); at beta=1, D=O(b). In particular
only polynomial cyclic degrees and polynomially many complete global
queries are required throughout the currently relevant interval
1/2<beta<=1. For beta>1 the overhead remains sufficiently small and
there is no cap at one from the transfer.

If the stronger prefix oracle from the original transfer is available,
an even simpler refinement uses its existing sparse-nonzero test with
the same D and then its prefix driver. That already proves (4), with
only linear dependence on D in the number of tested primes. The global
interpolation here supplies a weaker oracle contract, preserving the
complete divisor polynomial without prefix restrictions. Its quadratic
collision bound is explicitly paid.

The hypothesis concerns complete weighted query evaluation. A fast
algorithm on unwritten tables must still pay their constructor, encoding
and scalar-ring costs, as specified in the original notes. The present
postprocessing does not establish such an algorithm.

The [weight-one coefficient realization](weight-one-divisor-coefficients.md)
constructs signed cyclic images from odd Dirichlet-character coefficients.
Its signed interpolation uses exponent bit bound b+1 and sparsity cap 2D,
with the same complete cost conclusion.

## 4. Executable controls and limits

Run `make verify-sparse-divisor-residues`, or

```sh
python3 tools/check_sparse_divisor_residues.py > evidence/sparse-divisor-residues.json
```

The [core](../tools/sparse_divisor_residues.py) accepts an exact dense
residue callback and reuses the original small-prime generator and
the exact integer root-isolation helper. The
[checker](../tools/check_sparse_divisor_residues.py) supplies residues
from explicitly charged divisor lists and compares the reconstructed
integer polynomial to their independent full product.

The [receipt](../evidence/sparse-divisor-residues.json) covers 84 small
recoveries over characteristics 2,3,5 with exact and oversupplied sparsity
bounds. It records 34 zero residues caused by characteristic cancellation.
Two large controls use the independently trial-proved 60/80-bit prime
factor data from the moment checker. They recover the known smallest
factors from 607/807 supplied global residue queries, with maximum cyclic
degrees 4463/6203. No fast residue provider is implemented, and these
finite checks establish no new factoring exponent.

The ordinary transfer has not received independent mathematical review.
This work used substantial OpenAI Codex assistance and adds no Lean.
