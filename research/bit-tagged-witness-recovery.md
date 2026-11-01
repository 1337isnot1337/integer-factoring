# Recover witnesses with bit tags and deterministic isolation

## Results and required interface

A weighted witness circuit can reveal a witness's bits, not only whether
its polynomial is zero. Combining this with deterministic isolation
removes adaptive divisor-prefix queries.

There are three useful implementations:

| Recovery method | Extra time after small-factor preprocessing | Coefficient storage |
|---|---:|---:|
| Cyclic residue hashes, any fixed characteristic | `exp(O(q+log n))` | `n^O(1) 2^O(q)` |
| Coordinate cylinders, any fixed characteristic | `exp(O(q log(n+2)))` | `poly(n)` |
| Polynomial hashes and additive Fourier projection, any fixed characteristic | `exp(O(q+log n))` | `poly(n)` |

Here preprocessing removes all prime factors at most `2^s` and
`q>=ceil(n/s)` bounds the logarithm of the remaining divisor count.
The evaluator's own time is a multiplicative factor, and its own
control-space requirement must also be charged.

In particular, if a streamed evaluator over any fixed characteristic has
time `A(n)` and polynomial space, coordinate cylinders give complete
factorization in polynomial space and time

```math
\exp\!\left(O\!\left(
\log(A(n)+2)+\sqrt{n\log(n+2)}\right)\right).
\tag{1}
```

Over **every fixed characteristic**, additive Fourier projection improves this to

```math
\exp\!\left(O\!\left(\log(A(n)+2)+\sqrt n\right)\right)
\quad\text{with polynomial space}.
\tag{2}
```

Thus every fixed logarithmic saving exponent of the evaluator is
preserved in polynomial space. These are recovery transfers, not new
crossing-rank bounds. With the existing modular half-rate evaluator,
they give the established half exponent with polynomial space.

### A bit-multivariate circuit is essential

For a witness set `S subset {0,1}^n`, the required circuit computes

```math
C_S(z_0,\ldots,z_{n-1})
=\sum_{d\in S}\prod_{j=0}^{n-1}z_j^{d_j},
\qquad d=\sum_j d_j2^j,
\tag{3}
```

with coefficient one for each distinct witness. It must be correct over
every commutative algebra over the fixed coefficient field.
An unspecified circuit for the univariate polynomial `sum_d X^d`
is insufficient for the substitutions below: ordinary integer exponent
bits cannot be tagged by replacing `X` alone.

Throughout, `A(n)` uniformly bounds arithmetic operations, scalar-factor
queries, and structural work for inputs of at most `n` bits, independently
of the target coefficient algebra. A streamed polynomial-space evaluator
uses polynomially many control bits and polynomially many live algebra
elements, with polynomial-space local factor routines. The representations
below charge the actual bit cost of every algebra operation. A mere
circuit-size bound would not imply the stated space conclusions.

The multiplication network supplies precisely (3). Give each input
factor bit its own formal weight. The explicit settings-to-divisors
bijection in `MulCount.good_card` makes each divisor contribute one
multilinear monomial. The weighted bridge and source pins are in
[scale-one-modular-transfer.md](scale-one-modular-transfer.md), Section 4.
The generator and all local setting weights have polynomial overhead.

For a digit `u` at bit offset `o`, the cyclic tagged weight used below
is exactly

```math
X^{u2^o}
\left(1+\sum_{\substack{0\le j<d\\u_j=1}}\epsilon_{o+j}\right).
\tag{4}
```

It is generally **not**
`(X^(2^o)(1+epsilon_o))^u`. The latter tags the integer `u`
as a scalar in the coefficient field, instead of tagging its binary
digits.

## 1. Square-zero tags carry the witness bits

Fix a prime `p` and a commutative `F_p`-algebra `R`. Define

```math
\mathcal D_n(R)
=R[\epsilon_0,\ldots,\epsilon_{n-1}]
 \big/(\epsilon_i\epsilon_j:0\le i,j<n).
\tag{5}
```

An element is represented by `(a,b_0,...,b_(n-1))`, and multiplication
is

```math
(a,b)(c,d)=(ac,\;ad+cb).
\tag{6}
```

The algebra is commutative and has dimension `n+1` over `R`.
Only the scalar coordinate and the linear tags survive. In particular,

```math
\prod_{j:d_j=1}(1+\epsilon_j)
=1+\sum_j d_j\epsilon_j.
\tag{7}
```

Consequently a sum of tagged witnesses records its cardinality modulo
`p` and the sum of each witness bit modulo `p`. For `p=2`
the latter is bitwise XOR.

Whenever a bucket contains one witness, its scalar coordinate is one
and its tags are exactly that witness's bits. A bucket with several
witnesses can imitate this format. Therefore the algorithm decodes only
coordinates in `{0,1}`, requires scalar coordinate one, and **verifies
the resulting candidate**. For divisors the verification is ordinary
integer division, together with the desired endpoint conditions.

This verification prevents false outputs. An isolation lemma ensures
that a true witness is eventually decoded.
Generic witness recovery therefore assumes an efficient membership
verifier for `S`; for divisors this is polynomial-time integer division.

## 2. Cyclic hashes recover every divisor without prefixes

Take `R_r=F_p[X]/(X^r-1)` and substitute in (3)

```math
z_j=X^{2^j\bmod r}(1+\epsilon_j)
\in\mathcal D_n(R_r).
\tag{8}
```

The coefficient of `X^u` in the result is

```math
\sum_{\substack{d\in S\\d\equiv u\pmod r}}
\left(1+\sum_j d_j\epsilon_j\right).
\tag{9}
```

Suppose `|S|<=K` and every exponent is less than `2^n`.
For each fixed witness `d0`, the product
`product_(d!=d0)|d-d0|` has fewer than `n(K-1)` distinct prime
divisors. Thus among any

```math
M=n(K-1)+1
\tag{10}
```

distinct primes, at least one places `d0` alone in its residue class.
This argument is applied separately to each witness. It does not claim
that one prime isolates them all at once.

Evaluate (8) for the first `M` odd primes. Scan the buckets, decode
their tagged bits as in Section 1, and verify candidates. Every true
witness appears at some prime, and every accepted candidate is a true
witness. The union of accepted candidates is exactly `S`.

For divisor recovery, use the unprefixed circuit `C_N` and verify
`1<=d<=N` and `d|N`. It is harmless to include the known divisors
one and `N`. The smallest accepted value greater than one is the
least prime factor, including the prime-input endpoint.

### Cost and deduplication

The elementary prime supply proved in
[rough-input-transfer.md](rough-input-transfer.md), Section 2, puts the
first `M` odd primes below

```math
8M\left\lceil\log_2(M+2)\right\rceil
=O(nK(\log(n+2)+\log(K+1))).
\tag{11}
```

One tagged ring element has `(n+1)r` field coordinates.
Equation (6), with schoolbook cyclic convolution, gives polynomial
operation cost in `n,r`. Prime generation, bucket scanning, integer
candidate verification, and evaluation counts are polynomial in `n,K`.
No prefix query is used.

To find the smallest witness, keep only the smallest accepted candidate
so far. To list all witnesses in increasing order without storing a
`K`-element dictionary, repeat the scan and return the smallest accepted
candidate larger than the previously output one. There are at most
`K+1` scans. This remains polynomial in `n,K` and uses only a
constant number of `n`-bit candidate registers in addition to the
evaluator workspace. The output can be streamed.

After small-factor preprocessing, `K=2^q`. The half-rate choice
`q=O(sqrt(log n))` gives `r=n^(1+o(1))`, so a tagged ring element
has `n^(2+o(1))` coordinates. The faster abstract choice
`q=Theta(sqrt(n))` permits exponential workspace in this method.
The following methods avoid those large coefficient arrays.

## 3. A small coordinate restriction isolates some witness

**Cylinder lemma.** If `S subset {0,1}^n` is nonempty and
`|S|<=2^q`, some partial assignment on at most `q` coordinates
has exactly one extension in `S`.

Start with `S`. If it has at least two members, some coordinate
distinguishes two of them. Both resulting fibers are nonempty; retain
the smaller fiber, whose size is at most half. Repeat. Each chosen
coordinate is new, because all members already agree on previously
fixed coordinates. After at most `q` steps only one member remains.
This proves the lemma.

The lemma guarantees some singleton. It does not isolate every
prespecified witness within `q` coordinates. For example, the zero
vector among `{0,e1,...,em}` requires all `m` coordinates to
isolate it, although every `e_i` is isolated by its one nonzero bit.

Enumerate all coordinate sets of size at most `q` and every zero-one
assignment on them. Their number is at most

```math
\sum_{j=0}^{\min(q,n)}\binom nj2^j
\le(2n+1)^q
\tag{12}
```

for `q>=1`. All descriptions fit in `O(n+q log(n+2))` bits.

### Compute a tagged cylinder with the same circuit

Let `Z` be the coordinates forced to zero and `O` those forced to
one. Work only in `D_n(F_p)`, which has `n+1` coordinates.
For `T subset O`, evaluate (3) with

```math
z_j=
\begin{cases}
0,&j\in Z\cup T,\\
1+\epsilon_j,&j\notin Z\cup T.
\end{cases}
```

Take the signed sum of these values with sign `(-1)^|T|`.
For a fixed witness, this sum contributes
`sum_(T subset O minus supp(d))(-1)^|T|`. It is zero unless
all forced-one bits are present, and is one otherwise.
The result is therefore exactly

```math
\sum_{\substack{d\in S\\d|_Z=0,\ d|_O=1}}
\left(1+\sum_jd_j\epsilon_j\right).
\tag{13}
```

In characteristic two all the signs are one. This identity requires no
new circuit for the partial assignment; it uses at most `2^q`
substitutions into the same unprefixed bit-multivariate circuit.

Decode and verify the vector from every cylinder. The cylinder lemma
ensures that at least one true witness is found whenever `S` is
nonempty. If every candidate is rejected, `S` is empty.
The number of circuit evaluations and all additional work are
`2^O(q log(n+2))`. Only one partial assignment, one subset `T`,
one tagged accumulator, and one candidate need be live at a time.

## 4. Polynomial-space factoring in any fixed characteristic

For an integer `N>=2`, apply Section 3 to its **proper** divisors:

```math
S_N=\{d:1<d<N,\ d\mid N\}.
\tag{14}
```

The circuit for this set is obtained from `C_N` by subtracting the
two known monomials for one and `N`. They are distinct and each can
be evaluated in `O(n)` ring operations.

A verified candidate is a proper divisor. If none is found after the
complete cylinder family, the set (14) is empty, so `N` is prime.
A found divisor need not be prime or the smallest divisor. Recursively
split it and its cofactor until every remaining integer is declared
prime. This recursion has at most `2n-1` nodes, since its prime leaves
have product `N<2^n` and every leaf is at least two.
The least resulting prime is the least prime factor.

First remove all prime factors at most `2^s` by ascending trial
division, dividing out each discovered factor completely. Every divisor
found in this ascending procedure is prime: a composite candidate would
have had a smaller prime factor already removed. The remaining integer
and all of its subsequent factors have no prime factor at most `2^s`.
By the rough-input bound, all their witness sets have size at most
`2^ceil(n/s)`, with the original `n` used throughout.

For a chosen `q>=1`, set `s=ceil(n/q)`. The resulting time is

```math
\exp\!\left(O\!\left(
\log(A(n)+2)+\frac nq+q\log(n+2)\right)\right).
\tag{15}
```

The coefficient algebra in every circuit query has only `n+1`
coordinates over a fixed field. If the arithmetic evaluator has the
streamed-space property proved for the grid evaluator in
[rough-input-transfer.md](rough-input-transfer.md), its bit-space cost is
polynomial in `n`.

For a fully integer choice, put

```math
L=\lceil\log_2(n+2)\rceil,\qquad
q=\max(1,\lfloor\sqrt{\lfloor n/L\rfloor}\rfloor).
\tag{16}
```

Then `n/q+qL=O(sqrt(nL))`. Equation (15) proves (1).
The recursive factor search can be depth first; its pending cofactors
and all loop descriptors use polynomial space. It never stores all
cylinders, queried outputs, or witnesses.

## 5. Odd characteristic allows a smaller hash family in polynomial space

Assume now that the coefficient characteristic `p` is odd.
The following hash uses a second field, of characteristic two,
**only for integer control computations**. No embedding of that field
into the evaluator's coefficient field is asserted.

Choose the least positive `k` satisfying

```math
2^k>(K-1)\left\lceil\frac nk\right\rceil.
\tag{17}
```

For `K=2^q` this gives `k=O(q+log(n+2))`; for example,
`k=q+ceil(log2 n)` already satisfies (17).
Let `H=F_(2^k)` and put `m=ceil(n/k)`.
Read an `n`-bit witness in consecutive `k`-bit chunks,
with the last chunk padded by zeros. Interpret the chunks as field
elements and form

```math
D_d(T)=\sum_{a=0}^{m-1}\operatorname{chunk}_a(d)T^a
\in\mathbb H[T].
```

For `t in H` use the hash `H_t(d)=D_d(t)`.
It is linear in the binary witness bits.

For distinct witnesses `d,e`, the difference `D_d-D_e` is a
nonzero polynomial of degree at most `m-1`. It has at most `m-1`
roots. For a fixed witness `d0`, at most `(K-1)(m-1)` values
of `t` collide with any other witness. By (17), some `t` isolates
`d0`. This proves isolation of every witness by this hash family.

### Walsh projection of one bucket

Represent elements of `H` as `k`-bit vectors in its polynomial
basis. For `v,b in {0,1}^k`, let `v dot b` denote their binary
dot product. The usual binary character identity, now evaluated over
the odd coefficient field, is

```math
[a=b]=
2^{-k}\sum_{v\in\{0,1\}^k}(-1)^{v\cdot(a+b)}.
\tag{18}
```

The scalar `2^k` is nonzero in `F_p` and hence invertible.

For fixed `t,b,v`, substitute into (3)

```math
z_j=(-1)^{\,v\cdot H_t(e_j)}(1+\epsilon_j)
\in\mathcal D_n(\mathbb F_p).
\tag{19}
```

Linearity of the hash makes the sign on a witness monomial
`(-1)^(v dot H_t(d))`. Multiply the returned vector by
`2^(-k)(-1)^(v dot b)` and sum over `v`. By (18), the result is

```math
\sum_{\substack{d\in S\\H_t(d)=b}}
\left(1+\sum_jd_j\epsilon_j\right).
\tag{20}
```

Thus one bucket can be computed while keeping only one tagged
accumulator. Enumerate `t,b,v` sequentially. The complete procedure
uses at most `2^(3k)` circuit evaluations, plus polynomial work per
evaluation. Decoding and verification recover every true witness at
some bucket. To find the least proper divisor, keep only the least
verified candidate greater than one; if none exists after excluding
the known endpoints, the input is prime.

As in Section 2, sorted witness enumeration can be performed by repeated
scans with a current lower threshold. Its extra factor `K+1` remains
`2^O(q)` and requires no exponential dictionary.

### Uniform field construction and space

Enumerate monic degree-`k` binary polynomials until finding an
irreducible one. Trial division by every monic polynomial of degree
at most `k/2` is a sufficient test. The construction takes
`2^O(k)` bit operations and polynomial space in `k`.
Finite fields guarantee that the search terminates.

Afterward, polynomial arithmetic modulo the chosen irreducible
polynomial uses `poly(k)` time and space per operation. Enumerate
field elements with `k`-bit counters; no complete field table is
needed. The vectors `H_t(e_j)` can be generated by field powers
and stored in `O(nk)` bits, or regenerated when queried.
The coefficient algebra in (19) always has only `n+1` coordinates.

With `q<=n`, all control data have polynomial length in `n`.
The complete recovery has extra time `exp(O(q+log(n+2)))` and
polynomial space. Small-factor preprocessing costs `exp(O(n/q))`.
Take `q=ceil(sqrt(n))` to prove the odd-characteristic instance of (2).

For the existing construction one may fix `p=3`. The moment estimates
in [modular-crossing.md](modular-crossing.md) hold for every fixed prime:
choosing `2A_p/t=lambda ln2+C lambda^(3/4)` gives label energy
`(lambda ln2)^2/(4A_p)+O(lambda^(7/4))` and hence total channel
width `2sqrt(2)lambda+O(lambda^(3/4))`. This is the same explicit
endpoint argument as Section 6 of
[rough-input-transfer.md](rough-input-transfer.md), with `A2` replaced
by `A3`. The relaxed-width transfer therefore gives the exact half
endpoint. The local polynomial-space interface is proved in
[modular-factor-space-odd.md](modular-factor-space-odd.md), and the
streamed grid execution is proved in the rough-input note.
Consequently the odd-characteristic recovery here applies to an
actual half-rate factoring algorithm with polynomial space.

Equation (18) cannot be used in characteristic two: the two signs
coincide and `2^k` is zero. The next section changes the control
field characteristic and obtains the same recovery bound in that case.

## 6. A different control characteristic removes the remaining restriction

Fix the evaluator characteristic `p`. Choose a fixed prime `ell!=p`,
and let `E` be a fixed finite extension of `F_p` containing a primitive
`ell`-th root of unity `omega`. Such an extension has constant degree
for fixed `p,ell`. The two choices needed here are particularly simple:

| Evaluator characteristic | Control characteristic | Coefficient field and root |
|---|---:|---|
| `p=2` | `ell=3` | `E=F_4`, `omega²+omega+1=0` |
| fixed odd `p` | `ell=2` | `E=F_p`, `omega=-1` |

The evaluator already works over every commutative `F_p`-algebra.
It can therefore be evaluated over `D_n(E)` without changing its
identity or structural operation bound. The coefficient extension only
multiplies storage and field-operation costs by constants.

Choose the least positive integer `k` with

```math
Q=\ell^k>n(K-1),\qquad \mathbb H=\mathbb F_{\ell^k}.
\tag{21}
```

There is no embedding of the control field `H` into the coefficient
field `E`. They are used for different tasks. For a binary witness,
form the polynomial and its evaluation hash

```math
D_d(T)=\sum_{j=0}^{n-1}d_jT^j\in\mathbb H[T],
\qquad H_t(d)=D_d(t).
\tag{22}
```

If `d!=e`, some coefficient of `D_d-D_e` is `1` or `-1`.
It is nonzero over `F_ell`, so this difference is a nonzero polynomial
of degree at most `n-1`. For a fixed witness, the union of collision
roots with all other witnesses has size at most `(K-1)(n-1)<Q`.
Hence some `t` isolates that witness. Every witness is isolated at
some point in the family.

### Additive characters require only coordinate dot products

Represent `H` in a polynomial basis over `F_ell` and denote its
coordinate vector by `coords`. For `v in F_ell^k`, put

```math
\chi_v(a)=\omega^{v\cdot\operatorname{coords}(a)}.
```

Coordinate additivity makes each `chi_v` a character of the additive
group of `H`. The elementary product of `k` geometric sums gives

```math
[a=b]=Q^{-1}\sum_{v\in\mathbb F_\ell^k}
\omega^{v\cdot(\operatorname{coords}(a)-\operatorname{coords}(b))}.
\tag{23}
```

The scalar `Q` is a unit in `E` because `ell!=p`. If one coordinate
of `a-b` is nonzero, its corresponding geometric sum is zero; if all
coordinates vanish, the sum is `Q`. No field-trace pairing is needed.

For fixed `t,b,v`, substitute

```math
z_j=\omega^{v\cdot\operatorname{coords}(t^j)}(1+\epsilon_j)
\quad\text{in }C_S,
```

multiply the result by
`Q^(-1) omega^(-v dot coords(b))`, and sum over `v`.
The phase on a witness monomial is exactly `chi_v(H_t(d))`.
Equation (23) therefore returns the tagged bucket sum (20).
In characteristic two, the two-dimensional coefficient field `F_4`
supplies the three distinct phases, and the normalizing scalar
`3^k` equals one. The binary-Walsh collapse is avoided.

Decode only a scalar-one vector with zero-one tags, and verify its
candidate with the membership verifier. Since every witness is isolated
somewhere, this procedure recovers every witness. The minimum, sorted
enumeration, proper-divisor endpoint subtraction, and recursive complete
factorization work exactly as in the preceding sections.

### Uniform time and polynomial space

Enumerate monic degree-`k` polynomials over `F_ell`, and test each
by trial division by all monic polynomials of degree at most `k/2`.
The first irreducible polynomial constructs `H` in `ell^O(k)` time
and polynomial space in `k`. Field operations and basis-coordinate
operations afterward have polynomial cost in `k`. The constant
coefficient extension `E` and its root can be fixed once and for all.

Stream the loops over `t,b,v`; the full scan uses at most `Q³`
circuit evaluations. Only one tagged accumulator, polynomially many
control counters, and the `n` powers `t^j` are needed. Since

```math
k=O(q+\log(n+2)),\qquad Q=2^{O(q+\log(n+2))}
```

for `K=2^q`, the recovery costs `exp(O(q+log(n+2)))` extra time
and polynomial space whenever `q<=n` and the evaluator is streamed.
No `Q`-element coefficient array or lookup table is stored.

Combining this with small-factor trial division, set
`q=ceil(sqrt(n))` and `s=ceil(n/q)`. This proves (2) over
**every fixed coefficient characteristic**. In particular the original
characteristic-two modular evaluator and its local-space proof in
[modular-factor-space.md](modular-factor-space.md) already suffice;
switching to odd characteristic is optional.

The recovery overhead is `exp(O(sqrt(n)))` in polynomial space,
so it preserves any fixed logarithmic-power saving of a streamed
evaluator. The existing half-rate evaluator consequently gives an
unconditional half-rate factoring algorithm with polynomial space.
This result removes a recovery restriction; it supplies no improved
crossing family.

## 7. Compatibility with the Fourier input-gate primitive

The cyclic bit tags do not invalidate the compact Fourier input-gate
queries from
[fourier-multiplication-gates.md](fourier-multiplication-gates.md).
For an allowed digit interval `u=a+v`, `0<=v<2^t`, with
`2^t|a`, put `z=Omega^alpha X^(2^o)`.
The scalar coefficient is `z^a S_(2^t)(z)`.

For a free bit `0<=j<t`, its tagged coefficient is

```math
z^{a+2^j}
S_{2^j}(z)
S_{2^{t-j-1}}(z^{2^{j+1}}).
\tag{24}
```

Indeed, the values with that bit set are uniquely
`v=v0+2^j+2^(j+1)w` with
`0<=v0<2^j` and `0<=w<2^(t-j-1)`.
For a fixed high bit, the tagged coefficient is either zero or the
whole scalar sum, according to that bit of `a`.

There are only polynomially many digit bits. Division-free geometric
sums therefore compute every tagged component in polynomially many
ring operations. This establishes compatibility of the new witness
weights with that local primitive. It still supplies no closed global
Fourier contraction algorithm.

The additive-character weights in Section 6 are compatible as well.
For free digit bits, a bit weight `beta_j` changes the digit interval
sum into the product `product_j(1+beta_j z^(2^j))`, times its fixed
high-bit weight. Each tagged coordinate is obtained by replacing one
factor by `beta_j z^(2^j)`. Prefix and suffix products compute all
coordinates without dividing by possibly zero factors.

## 8. Exact verification and scope

The companion checker tests the square-zero multiplication and tagging
identities, cyclic bucket decoding, cylinder isolation and
inclusion-exclusion, the binary polynomial hashes, odd-characteristic
Walsh projections, additive projections over a different control
characteristic (including coefficient characteristic two), and small
complete factor searches.
It includes controls for incorrect digit-power tagging, collision
candidates requiring verification, failure to isolate every specified
witness by a small cylinder, and characteristic-two Walsh collapse.

The polynomial-space factoring statements use the streamed evaluator
and local-factor proofs separately. A short circuit description alone
does not imply polynomial-space evaluation at its claimed time.
No assertion here improves the crossing coefficient beyond
`rho=1/sqrt(2)`.
