# Scalar recovery with a virtual rough cofactor

## Result and scope

Small-factor preprocessing can preserve a sublinear workspace bound during
**complete prime factorization**, without storing the current cofactor. A
constant-state phase path also implements the characteristic-two recovery
queries over the scalar field itself. These two observations extend the
time-only [rough-input transfer](rough-input-transfer.md) to a simultaneous
time and space interface for a scalar evaluator.

The output format matters. Complete prime factorization has a short output;
listing every divisor does not. Section 7 proves that a halting deterministic
machine with read-only input and one-way write-only output needs
`Omega(n/log n)` workspace on infinitely many inputs if it lists all divisors.
Thus an extension of this recovery interface beyond saving exponent one must
not promise complete divisor enumeration at the same time and space bound.
Section 6 instead gives an output-sensitive bound with an additional
`O(log tau(N))` working bits and a polynomial factor in `tau(N)` in time.

No faster scalar evaluator is supplied here. The achieved factoring score
remains `1/2`.

### A scalar phase-count interface

Fix a prime coefficient characteristic `p`. Put `ell=2` when `p` is odd and
`ell=3` when `p=2`. Define a scalar function on `Z/ell Z` by

```math
\psi_p(e)=
\begin{cases}
(-1)^e\in\mathbb F_p,&p\ne2,\\
1,0,1\in\mathbb F_2\text{ at }e=0,1,2,&p=2.
\end{cases}
\tag{1}
```

For an at-most-`n`-bit integer `M>=2`, a read-only sequence
`c=(c_0,...,c_(n-1))` in `Z/ell Z`, an offset `z`, and one optional bit pin,
consider the scalar query

```math
Z_p(M,c,z;\eta)=
\sum_{ab=M}
 \psi_p\!\left(\sum_{j=0}^{n-1}c_j a_j-z\right)\eta(a,b)
\quad\in\mathbb F_p.
\tag{2}
```

Here `a_j,b_j` are the padded binary digits and `eta` is either `1`, `a_j`,
or `b_j` for one specified `j`. A bit pin can therefore be implemented by
requiring the specified factor-input bit to equal one.

**Recovery theorem.** Suppose one uniform deterministic finite-alphabet
machine answers (2), for all inputs padded to width `n`, in time `A(n)` and
workspace `S(n)`, where `A(n)>=n`. Its integer and phase inputs are read-only.
For every integer `1<=s<=n`, there is a uniform algorithm that outputs the
complete prime factorization of `N<2^n`, or its least prime factor, with

```math
\begin{aligned}
T(n)&\le A(n)\exp\!\left(O\!\left(s+\frac ns+\log(n+2)\right)\right),\\
W(n)&\le S(n)+O\!\left(s+\frac ns+\log(n+2)\right).
\end{aligned}
\tag{3}
```

The input is read-only. The output tape is one-way write-only; its contents
and head position are never read. The constants are independent of `N,n,s`.
They may depend on the fixed coefficient field and the fixed evaluator.

With `s=ceil(sqrt(n))`, the extra logarithmic time and workspace in (3) are
`O(sqrt(n))`. In particular, for every fixed `alpha>0`, an evaluator with

```math
A(n)=\exp(O(F_\alpha(n))),\qquad
S(n)=O(F_\alpha(n)),\qquad
F_\alpha(n)=\frac n{(\log(n+2))^\alpha}
\tag{4}
```

gives the same two bounds for complete factorization and the least factor.
This assertion has no `alpha<=1` recovery restriction. It does not remove
the separate rank restrictions on the repository's crossing families.

### Consequence for a written scalar clique evaluator

Let a scalar clique evaluator over `F_p`, on read-only tables with `P`
vertices and `lambda` bits per edge, have time `E_p(P,lambda)` and workspace
`W_p(P,lambda)`. There is an absolute constant `C` such that, whenever
`P>=2` and `P^40<=n`, put

```math
\lambda_* = \max\left(1,\left\lceil Cn/P^2\right\rceil\right).
```

The construction in Section 2 and on-demand table generation give

```math
\begin{aligned}
T_{\rm fact}(n)
&\le (E_p(P,\lambda_*)+1)
\exp\!\left(O_p\!\left(n/P+s+n/s+\log(n+2)\right)\right),\\
W_{\rm fact}(n)
&\le W_p(P,\lambda_*)+
O_p\!\left(n/P+s+n/s+\log(n+2)\right).
\end{aligned}
\tag{5}
```

The time multiplication in (5) pays for regenerating a table symbol on every
read. The evaluator's own workspace bound is an explicit hypothesis; a
time-only theorem for a padded evaluator does not imply it. Polynomially
long virtual table addresses and all input simulation are charged below.

### Uniformity when the coefficient field varies

The graph bound `C` in (5) is independent of the characteristic. For an odd
characteristic only scalar signs are added; in characteristic two there is
always the same three-state path. This differs from the `p`-tuple witness
construction in [fixed-prime-factor-transfer.md](fixed-prime-factor-transfer.md),
whose network mass depends on `p` and which instead has polynomially many
recovery queries.

Suppose a valid finite field `F_(p^e)` is supplied in its polynomial basis,
with `b=O(e log(p+1))` bits per field element, and one uniform scalar clique
evaluator has time `E_F(P,lambda)` and workspace `W_F(P,lambda)`. The same
construction gives the uniform bounds

```math
\begin{aligned}
T_{\rm fact}
&\le (E_F(P,\lambda_*)+1)\mathrm{poly}(n,b)
 \exp(O(n/P+s+n/s)),\\
W_{\rm fact}
&\le W_F(P,\lambda_*)+O(n/P+s+n/s+\log(n+2)+b).
\end{aligned}
\tag{5a}
```

The wrapper uses only prime-field values. Every table entry and query output
is a constant polynomial in the supplied basis; other coordinates are zero.
The wrapper accumulates integers modulo `p`, computes `Q^{-1} mod p`, and
streams the canonical field encodings. Its arithmetic has polynomial cost
in `b` and `O(b)` scratch bits. It does not construct a coefficient extension
or assume a valid field representation without the stated input promise.

Thus input-dependent fields with `b=polylog(n)` introduce no additional
logarithmic-power loss when the scalar evaluator's bounds are uniform over
those fields and the field representation has a uniform construction within
the claimed time and workspace. If it is generated rather than supplied,
that construction cost must also be charged. The recovery-query overhead is
`exp(O(sqrt(n)))`, not
polynomial; that difference from the fixed-prime scalar transfer is explicit.

There is also a stronger explicit constructor. Let `S_0` be the least power
of two with `n<=4 S_0^4`, put `A=2 S_0`, `f=8 S_0^2+1`, and choose a power
of two `2<=P<=A`. With `L=A/P`, define

```math
\lambda_\#=(8L^2+1)f,
\qquad P^2\lambda_\#\le1296n.
```

Section 2 uses spare bits already present in the original carry channels.
The edge width is unchanged and table queries remain polynomial-time. In place of
(5a), this gives

```math
\begin{aligned}
T_{\rm fact}
&\le (E_F(P,\lambda_\#)+1)\mathrm{poly}(n,b)
 \exp(O(s+n/s)),\\
W_{\rm fact}
&\le W_F(P,\lambda_\#)+
 O(n/P+n^{3/4}+s+n/s+\log(n+2)+b).
\end{aligned}
\tag{5b}
```

No exhaustive grouping search or local-setting enumeration is used in this
version. The control states travel between the explicit groups while the
original arithmetic carries use their smaller valid range. For a balanced
threshold, the same spare channels can carry the entire hash state and give
bucket counts directly, without the sum over phase queries.

## 1. Remove small primes without storing an integer cofactor

Define the smooth part and rough part, for the proof, by

```math
B=2^s,\qquad
U_s(N)=\prod_{\substack{r\le B\\r\text{ prime}}}r^{v_r(N)},
\qquad R_s(N)=N/U_s(N).
\tag{6}
```

There is a uniform bit-query algorithm for `R_s(N)` with time
`exp(O(s+log(n+2)))` and workspace `O(s+log(n+2))`. It reads only the
original `N`, `s`, and the requested bit position. In particular, it stores
neither `U_s(N)` nor `R_s(N)` as an `n`-bit integer.

The arithmetic dependency is the uniform logspace powering, iterated
multiplication, and division theorem of Hesse, Allender, and Barrington;
see their [primary paper, Sections 2 and 6](https://people.cs.rutgers.edu/~allender/papers/division.pdf).
The [2014 corrigendum](https://doi.org/10.1016/j.jcss.2013.09.002) retracts a
proof-theoretic part of Corollary 6.6, not the uniform arithmetic theorem
used here. Only the weaker logspace consequence is needed.

First, a candidate integer `r<=B` occupies `s+1` bits. Test primality by
trial division, using `O(s)` workspace and time polynomial in `B`. If `r`
is prime, determine `v_r(N)` by testing `r^e | N` for `1<=e<=n`.
Powering has inputs of polynomial length in `n`, because `s<=n` and
`r^e` has at most `n(s+1)` bits. Uniform logspace powering and division
answer each divisibility test without materializing `r^e`; comparison,
subtraction, and two-number multiplication are also logspace operations.
This is a fixed number of compositions, so it needs
`O(s+log(n+2))` bits and polynomial time in `n`, in addition to the primality
test. An exponent counter has only `O(log(n+2))` bits.

For each `r=2,...,B`, define a virtual `n`-bit slot

```math
A_r=
\begin{cases}
r^{v_r(N)},&r\text{ prime},\\
1,&r\text{ composite}.
\end{cases}
\tag{7}
```

Every requested bit of a slot has the preceding time and space bounds.
The list of these slots has virtual length `(B-1)n=exp(O(s+log(n+2)))`.
Run the uniform logspace iterated-multiplication algorithm on that list,
regenerating a slot bit whenever the simulated input head requests it.
The outer arithmetic machine stores `O(log(Bn))` bits. Its suspended state
and the slot generator together still use `O(s+log(n+2))` bits. Its number
of input reads and steps is polynomial in `Bn`; after substitution the
time remains `exp(O(s+log(n+2)))`.

The product is exactly `U_s(N)<=N<2^n`. Although the generic product routine
permits a much longer output, only its first `n` bits are needed. Run
logspace division on `N/U_s(N)` with this virtual denominator. This adds
only one more fixed layer of input simulation and proves the bit-query
claim for `R_s(N)`.

The use of prime slots in (7) is essential. Multiplying the corresponding
powers of every integer would overcount prime factors: for `N=12` and
`B=4` that incorrect product is `4*3*4=48`.

The same scan prints every prime `r<=B` with its exact positive multiplicity,
in increasing order. It keeps just the current candidate and exponent.
Every remaining prime factor of `R=R_s(N)` is greater than `B`, and hence

```math
\Omega(R)<n/s,\qquad
\tau(R)=\prod_r(v_r(R)+1)\le2^{\Omega(R)}
\le 2^{\lceil n/s\rceil}.
\tag{8}
```

The case `R=1` is detected by its virtual bits and ends the algorithm.
No factoring call is applied recursively to a quotient represented by another
factoring computation. Every rough-part bit request returns directly to the
original read-only `N` through the fixed arithmetic construction above.

## 2. Character sums require only scalar fields

Choose the least `k>=1` with

```math
Q=\ell^k>n(K-1),\qquad K=2^{\lceil n/s\rceil}.
\tag{9}
```

The control field `F_Q` is represented by a monic irreducible polynomial
over `F_ell`. Exhaustive monic-polynomial search and exhaustive proper-factor
tests construct it in `Q^{O(1)}` time and `O(k)` space. Field arithmetic and
the standard coordinate dot product are then effective with `O(k)` bits.
Here

```math
k=O(n/s+\log(n+2)),\qquad Q=\exp(O(n/s+\log(n+2))).
\tag{10}
```

For `t in F_Q`, hash the first factor by

```math
h_t(a)=\sum_{j=0}^{n-1}a_jt^j\in\mathbb F_Q.
\tag{11}
```

For every divisor `a` of `R`, at most `(n-1)(tau(R)-1)` choices of `t`
collide with another divisor. This follows by applying the polynomial root
bound to each nonzero difference of two binary-digit polynomials, which
remains nonzero in either characteristic two or three. Thus (9) provides
an isolating `t` for every actual divisor.

For a chosen `t,b,v in F_Q`, supply (2) with

```math
c_j=\langle v,\mathrm{coords}(t^j)\rangle,
\qquad z=\langle v,\mathrm{coords}(b)\rangle.
\tag{12}
```

These phase inputs are virtual read-only data. Powers and dot products
regenerate one coordinate in polynomial time and `O(k+log(n+2))` bits.
The complete `n`-symbol sequence is never stored.

For every `x in F_Q`, the map `v -> <v,coords(x)>` is zero when `x=0`
and has every value in `F_ell` exactly `Q/ell` times when `x!=0`.
The scalar table in (1) has `psi_p(0)=1` and
`sum_(e in F_ell) psi_p(e)=0`. Since `ell!=p`, `Q` is invertible in
`F_p`. Consequently

```math
Q^{-1}\sum_{v\in\mathbb F_Q}
 \psi_p(\langle v,\mathrm{coords}(x)\rangle)=[x=0].
```

Summing (2) over all `v` and applying this identity gives

```math
Q^{-1}\sum_{v\in\mathbb F_Q}Z_p(R,c(v),z(v);\eta)
=\sum_{\substack{ac=R\\h_t(a)=b_0}}\eta(a,c),
\tag{13}
```

where `b_0` denotes the chosen bucket label. In characteristic two, `Q=3^k`
is one in `F_2`, so the normalizing factor is one. This proof uses only
balanced linear forms and the scalar table; no roots of unity or coefficient
field extension are needed even in the argument. Every actual evaluator call
remains scalar-valued in `F_p`. A square-zero tag is unnecessary because the
single-bit pins already give all coordinate sums.

The same identity holds over any fixed nonzero finite commutative ring.
Choose a prime `ell` coprime to its characteristic and take the scalar table
`psi(0)=1`, `psi(ell-1)=-1`, with all other values zero. Then `Q=ell^k` is
a unit and the displayed argument is unchanged. An `ell`-state path realizes
the table with constant overhead. This observation does not assert field
rank factorizations over a ring; it is solely a scalar recovery interface.

### Realize the scalar phase query as a multiplication network

For odd `p`, the sign in (2) factors over individual input bits and a
constant offset. It changes only scalar local weights.

For `p=2`, attach a three-state path to the first factor's exposed input
bits. Start in state zero and use the deterministic transition

```math
q_{j+1}=q_j+c_ja_j\pmod3.
```

At the end place scalar weight `psi_2(q_n-z)`. Two binary bits encode each
state, and the unused fourth encoding is forbidden. A transition vertex
has at most six local settings and incident width at most five bits, so
its setting exponent plus incident width is at most eight. There are `n`
such vertices. The final scalar weight and the initial state pin have
constant cost. The exposed input ports gain one equality connection;
the multiplier construction already allows its stated constant number of
copy/equality edges. An additional bit pin only changes a local weight.

For every factor pair there is exactly one path extension. Thus there is
no extra multiplicity, and the contraction is exactly (2). The existing
free-output multiplication gadget from
[fixed-prime-factor-transfer.md](fixed-prime-factor-transfer.md) has total
mass `O(n)` and largest vertex mass `O(n^(3/4))`; the same is true after
adding this path. All local weights are in the original scalar field.

The generic deterministic grouping proof in that note therefore gives
`P^2 lambda=O(n)` and setting exponent `H=O(n/P)` when `P^40<=n`.
Enumerating the settings for one table entry requires `O(n/P)` stored
bits. Fixed-setting multiplier messages use at most `O(n^(3/4))` scratch
bits, which fit within `O(n/P)`. Grouping codes and their verification
also fit: `P^8 log(n+2)=O(n/P)` in this parameter range, after increasing
a fixed constant to cover small inputs.

A table symbol can consequently be regenerated in `2^{O_p(n/P)} poly(n)`
time and `O_p(n/P)` space, plus the phase and rough-input generators.
The encoded input-head position of the scalar evaluator takes
`O(n/P+log(n+2))` bits. Suspending that evaluator and regenerating each
requested symbol proves (5). In particular, merely substituting an oracle
for a fully written table is not being counted as a free operation.

### Preserve the explicit interpolation constructor

The [original top constructor](polynomial-top-queries.md) has groups indexed
by `a=0,...,P-1`. Each group contains the input gates with indices in the
contiguous block `I_a={aL,...,(a+1)L-1}`. Its boundary determines at most
one valid local setting, recoverable by interpolation and directed
propagation. Its edge width is `(8L^2+1)f`, and its query workspace is
`O(n/P+Af+log(n+2))` by
[streamed-top-space.md](streamed-top-space.md).

For each valid local setting, compute the first factor's local phase sum

```math
u_a=\sum_{j:\,a_j\text{ belongs to an input gate in }I_a}c_j a_j
\quad\pmod\ell.
```

Bits at padded positions `j>=n` contribute zero. The reconstructed input
digits are already available one gate at a time. Accumulating this sum needs
only one residue modulo `ell`, and a selected bit pin is checked while its
own input gate is reconstructed. This does not change the unique-setting
argument or require retaining any extra array of input bits.

No extra edge bits are needed for the phase channel. Write `d=2 S_0^2`,
`B=2^d`, and `M=B^4+1`, so `f=4d+1`. In a complete valid arithmetic
setting, every coefficient lies in `[0,M)`, the initial carry is zero,
and the next carry is `floor((coefficient+carry)/B)`. The bound

```math
C_B=B^3+B^2+B+1
```

is invariant: if `carry<=C_B`, then
`coefficient+carry<=B^4+C_B=B C_B+1`, whose quotient by `B` is at most
`C_B`. Thus every global valid carry is at most `C_B`. Since `B>=4`,
`C_B<2 B^3=2^(3d+1)`. Its original `f`-bit message consequently has
`d>=2` bits available for control states.

On every adjacent-group edge `(a,a+1)`, replace the existing carry message
by the packed value `carry+2^(3d+1)q`. Decode this value before supplying
the original boundary carry to the interpolation/propagation routine.
Reject a decoded carry larger than `C_B` or a phase code `q>=ell`.
Nonadjacent edges keep their original zero padding. The restriction on
carry values can remove locally valid boundaries, but removes no complete
valid arithmetic setting, by the invariant just proved.

Let `q_0=0` be implicit. For each group before the last, require

```math
q_{a+1}=q_a+u_a\pmod\ell.
```

The incoming phase for a noninitial group is its prescribed boundary label.
At the last group, put scalar weight `psi_p(q_(P-1)+u_(P-1)-z)` instead of
an outgoing phase condition. A failed original table query or bit pin has
weight zero. An original global factor-pair setting has exactly one complete
phase path, so the augmented scalar contraction is exactly (2).

The entire packed carry and phase message still has `f` bits. The original
edge width is therefore unchanged. Since `P<=A` and the original constructor
proves `A^2 f<=144n`,

```math
P^2\lambda_\#=(8A^2+P^2)f
\le9A^2f\le1296n.
```

A scalar table query still has polynomial bit time in `n` and the coefficient
length, apart from its virtual-input accesses. Its workspace is
`O(n/P+Af+log(n+2)+b+k)`; `Af=O(n^(3/4))`. Phase-channel checks take constant
additional space. Simulate each read of the scalar evaluator's table input
by this routine, and then substitute the rough input and hash controls as
in Section 4. This proves (5b), including all table addresses and input
regeneration costs.

### Carry an entire hash state at a balanced threshold

The same codec can carry any binary control state of at most `d` bits,
not only a phase modulo two or three. This gives a direct scalar bucket
constructor. Use the binary control field `F_(2^k)`, with
`2^k>n(K-1)` and `K=2^ceil(n/s)`. When `k<=d`, let each group's state
update be addition of its partial hash

```math
u_a=\sum_{j:\,a_j\text{ belongs to group }a}a_j t^j
\quad\in\mathbb F_{2^k}.
```

Pack the binary field encoding into the spare `d` carry bits. The initial
state is zero; the terminal weight is the scalar indicator
`[q_(P-1)+u_(P-1)=b_0]`. Each valid arithmetic setting has one control
extension exactly when its first factor lies in bucket `(t,b_0)`. With an
optional bit pin, the clique contraction is therefore directly the right
side of (13), over any coefficient field. Hash-field addition and scalar
coefficient addition are separate operations; their characteristics need
not agree.

For example, choose `s=ceil(2 sqrt(n))` for `n>=1024`. Then

```math
k\le\lceil n/s\rceil+\lceil\log_2 n\rceil+1
\le\sqrt n\le d.
```

The middle inequality follows from
`sqrt(n)/2+log_2(n)+3<=sqrt(n)` in this range. The finitely many smaller
inputs can be handled directly. Thus the balanced wrapper can obtain all
bucket counts and coordinate sums by direct scalar clique queries, without
the `2^k`-term projection sum and without increasing the original edge width.
Its asymptotic bounds remain (5b), because both the old and the reduced
numbers of queries are polynomial in `n,2^k`. This improves the explicit
reduction and its query count, not the achieved factoring exponent.

## 3. Short handles for both factors

For each pair `(t,b_0)`, use (13) first with `eta=1`, and then with each
of the `2n` single-bit pins. Require the scalar count to equal one and
each returned coordinate to lie in `{0,1}`. If they do, the coordinates
define two virtual integers `d(t,b_0)` and `e(t,b_0)`. The handle consists
of `t,b_0` and the common control-field representation, using `O(k)` bits.
Every coordinate is regenerated when needed.

A count of one modulo `p` does **not** prove that a bucket is a singleton.
Check the virtual pair by the streamed multiplication comparator in
[sublinear-space-half-rate.md](sublinear-space-half-rate.md), Section 3.
It uses `O(log(n+2))` additional bits, `O(n^2)` coordinate queries, and
accepts exactly when `de=R`. Every accepted first coordinate is a divisor
of `R`; every actual divisor is accepted at an isolating handle.

For example, in characteristic two the three divisors `1,2,4` of `12`
have the same `t=1` ternary digit-sum hash. Their bucket has scalar count
one, but the coordinate sums give `d=1 xor 2 xor 4=7` and
`e=12 xor 6 xor 3=9`. The product is `63`, so verification rejects it.

Keep only a constant number of handles while scanning all `Q^2` choices.
Successor scans enumerate the actual divisors of `R` without duplicate
output, and the divisor-product tests and power scans from Sections 4–5
of the same note recover the large primes and their multiplicities.
The number of scans, queries, and virtual coordinate accesses is a fixed
polynomial in `n,Q`; no cofactor list or list of handles is retained.

All primes of `R` exceed `B`, so printing them after the small-prime scan
produces complete factorization in increasing prime order. A least-factor
version can stop at the first small prime, or retain the minimum large
prime handle when no small prime exists.

## 4. Full cost and uniformity

One scalar evaluator call on the virtual rough input and virtual phase
sequence takes at most

```math
A(n)\exp(O(s+\log(n+2)))\mathrm{poly}(k)
```

steps. Its suspended state uses `S(n)` bits; the rough-part generator and
phase generator use `O(s+k+log(n+2))` additional bits. They execute
sequentially, so their spaces add rather than multiply.

Field construction, bucket scans, bit checks, multiplication checks,
prime filtering, and multiplicity recovery have a fixed polynomial cost
in `n,Q`, besides the calls. The smooth-prime scan has
`exp(O(s+log(n+2)))` cost. Substitute (10) to obtain (3).

All parameter choices are integer operations. In particular,
`ceil(sqrt(n))` is computed by an integer square root, the control field is
constructed by a terminating finite search, and no unknown factorization
is used in selecting `K`. The input width remains the original `n` even
when the rough cofactor is shorter. The scalar evaluator must cover this
padded convention, or its bounds must be taken over all lengths at most `n`.

The only simulations of virtual inputs have a fixed number of arithmetic
layers plus one evaluator layer. The multiplication checks call bit
generators, not factoring on a newly generated input. This excludes an
unpaid chain of recursively factored virtual cofactors.

The choice `s=ceil(sqrt(n))` gives `s+n/s=O(sqrt(n))`; for any fixed
`alpha`, `sqrt(n)=o(F_alpha(n))`. This proves the stated all-exponent
recovery corollary. It concerns the listed scalar queries and their
actual time and workspace bounds, not an uncharged arithmetic-circuit
size assertion.

### An adaptive threshold from the virtual cofactor length

The overhead can also depend on the input's remaining rough part. Define
`m_s=0` if `R_s(N)=1`; otherwise let `m_s` be the binary length of `R_s(N)`.
Its value is obtained by scanning the virtual cofactor's bits, using only
an `O(log(n+2))`-bit position counter in addition to the bit generator.
The sequence `m_s` is nonincreasing. Put

```math
\sigma(N)=\min\{s\ge1:m_s\le s^2\}.
\tag{13a}
```

Test `s=1,2,...` in order and stop at the first successful test. This
terminates by `s=ceil(sqrt(n))`. Print small prime factors only on their
first encounter: at stage `s`, scan the prime candidates in
`(2^(s-1),2^s]`. The bit generator for `R_s` still starts from the original
input, so this scheduling needs no saved smooth product or recursive
cofactor. All earlier preprocessing stages together cost
`exp(O(sigma+log(n+2)))` time and `O(sigma+log(n+2))` workspace.

At the final stage the support bound strengthens to
`tau(R_sigma)<=2^ceil(m_sigma/sigma)`. Since `m_sigma<=sigma^2`,
the same recovery gives

```math
T(N)\le A(n)\exp(O(\sigma(N)+\log(n+2))),\qquad
W(N)\le S(n)+O(\sigma(N)+\log(n+2)).
\tag{13b}
```

If the rough part equals one, no evaluator call is needed. Otherwise the
evaluator may be called at its actual width `m_sigma` instead of padding
to `n`; then `A(n),S(n)` in (13b) can be replaced by
`A(m_sigma),S(m_sigma)`. The polynomial overhead still depends on the
original `n`, because each virtual bit reads that original input.

This stopping rule is within a constant factor of the best threshold for
the visible preprocessing/support tradeoff. More precisely, define

```math
C(N)=\min_{1\le t\le n}\left(t+\frac{m_t}{t}\right).
```

Then `sigma(N)<=C(N)<=2 sigma(N)`. The upper bound uses `t=sigma`.
For the lower bound, any `t>=sigma` already costs at least `sigma`.
If `t<sigma` and `sigma>=2`, monotonicity and the preceding failed test
give `m_t>=(m_(sigma-1))>(sigma-1)^2`, so
`t+m_t/t>=2 sqrt(m_t)>2(sigma-1)>=sigma`. The case `sigma=1` is immediate.
This comparison concerns this threshold schedule; it is not an optimality
claim about all factoring algorithms.

## 5. Short arithmetic outputs remain available

Encode complete factorization as the sorted prime list with repetitions
and separators. Its length is `O(n)`, because the sum of prime bit lengths
is at most `n+Omega(N)` and `Omega(N)<n`. A requested bit of this encoding
can be generated by replaying the factorization transducer and counting
its output symbols. The extra counter occupies `O(log(n+2))` bits.

Consequently any fixed deterministic logspace computation on that encoding
and the original `N` preserves (3): each simulated input read replays the
factorization, there are only polynomially many such reads, and the
additional workspace is logarithmic. This includes `tau(N)`, the Moebius
function, Euler's totient, and each fixed-power divisor sum. Their usual
prime-power formulas use uniform logspace powering, iterated multiplication,
and exact division; their output lengths are polynomial in `n`.

This closure does not provide arbitrary interval divisor counts from the
factorization in polynomially many bit operations. Such counts remain a
separate evaluator task. It also does not turn an exponentially long divisor
list into a short output.

## 6. Output-sensitive divisor enumeration

Suppose a complete-factorization transducer takes time `T_f(n)` and workspace
`W_f(n)` in the stated read-only input and write-only output model, and emits
sorted distinct primes with their exponents. The transducer above has this
form, and its output has `O(n)` symbols. The same replay argument makes its
output a virtual read-only input for subsequent arithmetic routines.

For a particular input `N`, write `D=tau(N)`. All divisors, including one and
`N`, can be emitted in increasing order with

```math
T_{\rm divisors}(N)\le T_f(n)(n+D)^{O(1)},\qquad
W_{\rm divisors}(N)\le W_f(n)+O(\log(D+1)+\log(n+2)).
\tag{14}
```

For the proof only, write `N=prod_(i=1)^r p_i^(e_i)` in increasing prime
order. The algorithm computes `D=prod_i(e_i+1)` and retains its binary
representation. This is permitted in (14), and costs `O(log(D+1))` bits.
For each integer `0<=u<D`, use mixed-radix digits

```math
u=\sum_{i=1}^r\alpha_i\prod_{j<i}(e_j+1),
\qquad 0\le\alpha_i\le e_i,
```

and define a virtual integer

```math
d(u)=\prod_{i=1}^r p_i^{\alpha_i}.
\tag{15}
```

This is a bijection from the indices onto the divisors, but the mixed-radix
order is not generally numerical order. A requested digit `alpha_i` is found
by rescanning the prime-exponent list and successively dividing a copy of `u`
by the earlier radices. This uses `O(log(D+1)+log(n+2))` bits. A requested bit
of `d(u)` is then a fixed composition of uniform logspace powering and iterated
multiplication on at most `n` virtual `n`-bit factors. Only the index `u`, a
constant number of its copies, the logarithmic arithmetic state, and one
factorization replay are live. Its time is `T_f(n) poly(n)` and workspace is
`W_f(n)+O(log(D+1)+log(n+2))`.

Retain indices for the preceding output and the best candidate in a full
scan of `0<=u<D`. Compare candidates by their virtual bit streams. Emit the
smallest value strictly larger than the preceding output, using a zero
sentinel on the first pass. The bijection proves that there are exactly `D`
successful passes. There are `O(D^2)` comparisons and a constant number of
live indices; output bits are regenerated from the winning index. This proves
(14). Unsorted enumeration needs only one scan.

Testing each selected divisor against read-only interval or prefix data and
incrementing an `O(log(D+1))`-bit counter gives exact interval and prefix
counts within (14) as well. The direct `tau(N)` computation in Section 5 is
faster because it does not enumerate divisors.

Combining (3) and (14) gives an all-exponent recovery interface with the
additional logarithmic-time and workspace term `O(log(tau(N)+1))` precisely
where explicit divisor enumeration needs it. It does not impose that term on
complete prime factorization or the short arithmetic outputs of Section 5.

## 7. An output-size obstruction for listing all divisors

Let a halting deterministic finite-alphabet machine have a fixed number of
work tapes, `S` visited working bits, an `n`-bit endmarked read-only input,
and one-way write-only output. There are at most

```math
\mathrm{poly}(n,S)c^S
\tag{16}
```

internal configurations for a fixed `c`. The output content and output-head
position cannot affect a future transition, so they do not enter this count.
If an internal configuration repeated, the machine would repeat the same
future computation and never halt. Therefore its running time and output
length are bounded by (16).

Let `N_k` be the product of the first `k` primes and let `n_k` be its bit
length. The elementary prime-supply estimate in
[rough-input-transfer.md](rough-input-transfer.md), Section 2, gives
`p_k=O(k log(k+2))`. Hence

```math
n_k=O(k\log(k+2)),\qquad k\le n_k,\qquad
\tau(N_k)=2^k=\exp(\Omega(n_k/\log(n_k+2))).
\tag{17}
```

Any explicit list containing all divisors has at least `tau(N_k)` symbols.
Combining (16) and (17) yields

```math
S(n_k)=\Omega(n_k/\log(n_k+2))
\tag{18}
```

on this infinite sequence. If `S(n_k)>n_k`, the conclusion is immediate;
otherwise the polynomial factor in (16) contributes only `O(log n_k)`
to its logarithm. This proves (18) in both cases.

It follows that for any fixed `alpha>1`, complete divisor listing cannot
have worst-case time `exp(O(F_alpha(n)))` or workspace `O(F_alpha(n))`
in this model. This is an output requirement, not a lower bound for finding
one factor, prime factorization, or short arithmetic outputs.

## Checks and attribution

Run:

```sh
python3 tools/check_output_sensitive_space_transfer.py \
  --output evidence/output-sensitive-space-transfer.json
```

The checker independently verifies the smooth/rough identity, prime-power
slots, the rough support bound, complete factorization after splitting,
and the adaptive threshold comparison on 1,027 inputs. It also verifies
the scalar phase path against balanced-form projections and bucket counts in
several fixed characteristics, the scalar filter over composite coefficient
rings, and the necessary rejection of a colliding bucket. Finite controls
for mixed-radix enumeration and primorial outputs check the elementary
numerical identities used in Sections 6–7. It also checks the explicit
constructor's padded width bound in 64,443 parameter cases and contracts
108 small grouped phase paths, including rejection of the spare state code.
These path controls concern the new channels; they do not reimplement the
existing interpolation constructor.

These checks do not implement the uniform logspace arithmetic theorem or
the asymptotic clique evaluator. The complexity claims are the written
composition arguments above, using the cited arithmetic theorem, the existing
grouping construction, and existing paired-descriptor verification.

The underlying multiplication network is due to Isaac Rudich and
Louis-Martin Rousseau; the crossing/grid evaluator used elsewhere in the
repository is due to Claude/Anthropic. This note used substantial OpenAI
assistant support. No Lean verification or improved numerical exponent is
claimed.
