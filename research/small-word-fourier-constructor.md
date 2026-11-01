# A multiplication constructor with logarithmic arithmetic words

## Result and scope

The multiplication network can use `O(log n)`-bit arithmetic residues,
`O(n/log n)` gates, and total information mass `O(n)`. Its largest gate has
mass `O(sqrt(n log n))`. An explicit grouping gives polynomial-time table
queries for every power of two

```math
2\le P\le A=\Theta\!\left(\sqrt{n/\log n}\right),
```

with `P^2 lambda=O(n)` and `O(n/P)` query workspace. The input is read-only.
The original constructor instead chooses `A=Theta(n^(1/4))` and arithmetic
words of `Theta(sqrt(n))` bits.

The change uses a prime arithmetic modulus with a prescribed root of unity.
An effective form of Linnik's theorem makes elementary searches for that
prime and root polynomial in the **original input length `n`**. The
polynomial degree can be large. No practical fast prime search is asserted.
Section 9 gives a CRT variant with polynomial-time metadata for any allowed
transform aspect ratio. A slightly larger digit size then gives
`P^2 lambda<=(36+o(1))n`. Grouping the carries with their inverse gates
reduces this to `(24+o(1))n` when `P=o(sqrt(A))`, as proved in Section 10.
These changes preserve the logarithmic saving exponent.

This note gives a separate ordinary constructor, with complete transform,
carry, grouping, and space arguments. It uses the same gate topology as
[the original explicit top constructor](polynomial-top-queries.md), due to
Isaac Rudich and Louis-Martin Rousseau. It does not establish a stronger
crossing evaluator or an achieved factoring exponent above `1/2`. In
particular, a larger allowable top-level `P` does not by itself extend any
later grid recurrence's parameter range.

## 1. Choose the parameters and arithmetic field

It suffices to treat `n>=64`; the finitely many smaller input lengths have
a direct finite algorithm. All logarithms in parameter choices are base two.
Put

```math
m=\lceil\log_2(n+2)\rceil.
```

Let `A>=2` be the largest power of two with `A^2 m<=n`, and define

```math
h=A/2,\quad T=A^2,\quad H=T/2=Ah,\quad
d=\lceil n/H\rceil,\quad B=2^d,\quad D=T B^2.
\tag{1}
```

These are integer computations. Maximality of `A` gives

```math
\frac{n}{4m}<A^2\le\frac nm,
\qquad 2m\le d<8m+1,
\qquad n\le Hd<n+H.
\tag{2}
```

In particular `A=Theta(sqrt(n/log n))`, `d=Theta(log n)`, and
`D=n^{O(1)}` as an integer, not merely as a bit length.

Let `q` be the least prime congruent to one modulo `D`. Linnik's theorem
gives absolute effective constants `C,L` such that

```math
q\le C D^L=n^{O(1)}.
\tag{3}
```

For a primary reference, Xylouris states the effective Linnik theorem and
proves an admissible exponent `5.18` in
[Acta Arithmetica 150 (2011), Theorem 1.1](https://www.impan.pl/shop/publication/transaction/download/product/83078),
[DOI 10.4064/aa150-1-4](https://doi.org/10.4064/aa150-1-4).
Only some absolute effective exponent is used here.

Search `1+D,1+2D,...` until a prime is found, testing primality by trial
division. By (3), the search has polynomial bit time in `n`. A constant
number of integers of size at most `q` and trial-division counters occupy
`O(log q)=O(log n)` bits. This argument charges the entire search; the
modulus is not advice and is not supplied by an uncharged number-theory
oracle.

Put

```math
f=\lceil\log_2 q\rceil=O(\log n).
```

The prime `q` is odd, `T` divides `q-1`, and `q>D`. Find a quadratic
nonresidue `g` by testing `g^((q-1)/2)=-1 mod q` in increasing order, then
set

```math
\omega=g^{(q-1)/T}\pmod q,\qquad \xi=\omega^A.
\tag{4}
```

A nonresidue exists because exactly half the nonzero field elements are
squares. The search uses at most `q` tests, so it again has polynomial time
in `n` and `O(log n)` workspace. It does not factor `q-1` or seek a
generator of the full multiplicative group. Since
`omega^(T/2)=-1`, `omega^T=1`, and `T` is a power of two, `omega` has
exact order `T`. Consequently `xi` has exact order `A`. Modular powers,
inverses, and ordinary arithmetic on these `f`-bit words have polynomial
bit time and `O(log n)` scratch space.

The arithmetic field `F_q` is separate from the coefficient field of the
eventual tensor contraction. Gate labels encode arithmetic residues; local
weights below are zero or one and can be evaluated over any coefficient
field.

## 2. The transform topology and exact multiplication

Write two candidate factors as base-`B` words with `H` digits,

```math
x=\sum_{j=0}^{H-1}x_j B^j,\qquad
y=\sum_{j=0}^{H-1}y_j B^j,
\qquad 0\le x_j,y_j<B.
```

Their digit arrays are zero beyond index `H-1`. A length-`T` cyclic
convolution equals the ordinary convolution because `T=2H`. Its integer
coefficients satisfy

```math
c_j=\sum_{u+v=j}x_u y_v,
\qquad 0\le c_j\le H(B-1)^2<q.
\tag{5}
```

The last inequality follows from `q>D=T B^2`. Thus an inverse transform
over `F_q`, read in canonical residues, gives the exact integer coefficients
without modular ambiguity.

Here is the complete two-dimensional transform. Use input gates
`FA_x[r],FA_y[r]`, middle gates `MID[k]`, and inverse gates `IB[r]`, with
`0<=r,k<A`. An input gate `FA_x[r]` stores the `h` digits
`x_(r+A i)`, `0<=i<h`, and sends to `MID[k]`

```math
X_{r,k}=\omega^{rk}\sum_{i=0}^{h-1}x_{r+Ai}\xi^{ik}.
\tag{6}
```

The second input family sends `Y_(r,k)` by the same formula. The setting
of `MID[k]` is its two incoming `A`-vectors. It computes

```math
\begin{aligned}
\widehat x_{k,v}&=\sum_{r=0}^{A-1}X_{r,k}\xi^{rv},\\
\widehat y_{k,v}&=\sum_{r=0}^{A-1}Y_{r,k}\xi^{rv},\\
Z_{k,r}&=A^{-1}\sum_{v=0}^{A-1}
 \widehat x_{k,v}\widehat y_{k,v}\xi^{-rv},
\end{aligned}
\tag{7}
```

and sends `Z_(k,r)` to `IB[r]`. The setting of `IB[r]` is its incoming
vector indexed by `k`; its output at index `i`, `0<=i<A`, is

```math
c_{r+Ai}=A^{-1}\sum_{k=0}^{A-1}Z_{k,r}\omega^{-rk}\xi^{-ik}.
\tag{8}
```

These formulas are ordinary length-`T` DFT multiplication. Indeed, the
frequency `k+Av` has transform

```math
\sum_j x_j\omega^{j(k+Av)}
=\sum_r\omega^{rk}\xi^{rv}\sum_i x_{r+Ai}\xi^{ik}
=\widehat x_{k,v}.
```

Expanding the length-`T` inverse transform at position `r+Ai`, first summing
over `v` and then over `k`, gives (7)–(8). The factors `A^{-1}` exist
because `q` is odd and `A` is a power of two. This proves (5) is exactly
the array produced by the network.

There are two complete input-to-middle bipartite layers and one complete
middle-to-inverse layer. Each inverse gate connects to the `A` carry gates
at positions `r+Ai`. Every arithmetic message occupies `f` bits; codes
outside `[0,q)` are rejected when used as arithmetic residues.

## 3. Carry gates and the factor-pair bijection

Let

```math
C_*=H(B-1).
\tag{9}
```

There is one carry gate at every position `0<=j<T`. Its setting is the
incoming coefficient in `[0,q)` and incoming carry in `[0,C_*]`. With
`N_j` the `j`th base-`B` digit of the given read-only input, its weight is
one exactly when

```math
(c_j+v_j)\bmod B=N_j,\qquad
v_{j+1}=\left\lfloor\frac{c_j+v_j}{B}\right\rfloor\le C_*.
\tag{10}
```

The first gate also requires `v_0=0`; the last requires `v_T=0`.
Consecutive carry gates have one equality edge carrying `v_(j+1)`.
The message in (10) always fits an `f`-bit word, even for a setting of
weight zero: `C_*<q` and `B>=2` give
`floor((q-1+C_*)/B)<q`.

For true convolution coefficients, the carry bound is invariant. Equations
(5) and (9) give

```math
c_j+v_j\le H(B-1)^2+H(B-1)=HB(B-1)=B C_*.
```

Starting from zero, every true carry lies in the required range. Since
`x,y<B^H`, their product is less than `B^T`, so its final carry is zero.

For each pair `xy=N`, the two input words have exactly one extension to
all transform and carry settings, and every weight is one. Conversely,
any valid global setting determines input words `x,y`; the transform
identities determine their exact coefficients. Multiplying (10) by `B^j`
and summing makes the carries telescope, giving

```math
xy=\sum_{j=0}^{T-1}c_j B^j
=\sum_{j=0}^{T-1}N_j B^j+v_T B^T-v_0=N.
```

The construction therefore has one valid setting per ordered factor pair.
For `N>=2`, zero candidate factors cannot give a valid setting. Every
actual positive divisor fits in the `H`-digit input word because
`N<2^n<=B^H`.

Local weights can impose a divisor prefix, pin one input bit of either
factor, or multiply prescribed bit weights. These changes preserve the
unique arithmetic extension. The contraction then gives the corresponding
weighted divisor count in the chosen coefficient field.

## 4. Information mass

Define a gate's mass as the logarithm of its setting count plus the sum
of widths of its incident messages, as in the existing grouping arguments.
The setting exponents have the following bounds.

| Gate family | Number | Setting exponent per gate |
| --- | ---: | ---: |
| Both input families | `2A` | `hd` |
| Middle | `A` | `2Af` |
| Inverse | `A` | `Af` |
| Carry | `T` | `2f` |

There are `2T` input-to-middle edges, `T` middle-to-inverse edges, `T`
inverse-to-carry edges, and `T-1` carry edges. Consequently the total
setting exponent is at most `T d+5T f`, and twice the total edge width
is less than `10T f`. The total mass is therefore at most

```math
T d+15T f=O(n).
\tag{11}
```

Here `T=A^2=Theta(n/log n)` and `d,f=O(log n)`. Each non-carry gate has
mass `O(A(d+f))=O(sqrt(n log n))`; each carry has mass `O(f)`.
The total number of gates is `T+4A=O(n/log n)`.

All local setting and message routines have polynomial bit time. Their
integer arrays have the sizes just charged. In particular, taking residues
modulo `q` inside a local gate is an ordinary operation on `O(log n)`-bit
integers, not a unit-cost operation on the full `n`-bit input.

## 5. Explicit groups and unique local settings

Choose a power of two `2<=P<=A`, put `L=A/P`, and partition `[0,A)` into
the contiguous blocks `I_a={aL,...,(a+1)L-1}`. Group `a` contains the four
gates `FA_x[r],FA_y[r],MID[r],IB[r]` for each `r in I_a`. Put carry
`j=iA+r` in the group containing `i`.

Each group contains `L<=h` middle gates. Given every boundary message,
an input gate inside the group has at least `A-L>=h` prescribed outputs
to external middle gates. Choose any `h` of them. Dividing (6) by its
nonzero twiddle gives evaluations of a polynomial of degree less than `h`
at distinct points `xi^k`. Those points are distinct because `xi` has order
`A`. Interpolation over the prime field therefore recovers at most one
input digit vector. Reject if a recovered coefficient is at least `B`,
or any other prescribed outgoing value or input-bit condition fails.

Once the contained input gates are determined, all middle, inverse, and
carry settings are determined in their directed order. Their incoming
messages either come from an already determined local predecessor or from
the prescribed boundary. Check all arithmetic ranges, local conditions,
and outgoing boundary messages. This proves that a prescribed group
boundary has at most one valid original setting. A table value is zero
or the weight of this unique setting.

The original cross-group edge count applies to this topology exactly:

```math
Y_{ab}=\left(8L^2+[|a-b|=1]\right)f.
\tag{12}
```

The two input layers contribute `4L^2`, the middle/inverse layer `2L^2`,
and the inverse/carry layer `2L^2` messages between distinct groups. The
carry segments are contiguous and contribute one message only between
consecutive groups. Pad each pair to the common width

```math
\lambda=(8L^2+1)f,
\qquad P^2\lambda=(8A^2+P^2)f\le9A^2f=O(n).
\tag{13}
```

Use a common deterministic ordering of original messages at both ends and
require all padding to be zero. This produces an ordinary clique whose
contraction is the weighted divisor count.

## 6. Polynomial-time table queries with `O(n/P)` workspace

A queried group boundary has `(P-1)lambda=O(n/P)` bits. Store that packed
boundary or read it through an equivalent indexed interface. All original
edge endpoints and offsets are generated by counters; no edge list is
stored. The `q,omega,xi` metadata use `O(log n)` bits, and regenerating
them has the polynomial time and logarithmic space from Section 1.

Interpolate one input gate at a time. A result coefficient array and a
current Lagrange numerator array each have `h` entries in `F_q`, using
`O(Af)` bits. A denominator, powers, indices, and running values require
`O(f+log n)` more bits. Evaluate or inspect the recovered digits as needed,
then discard the array. If another message needs that input gate, repeat
the interpolation. This has polynomial time because a table query uses
only polynomially many fixed-depth message computations.

In more detail, an input-output procedure does one interpolation and
evaluation. A middle-output procedure computes (7) with loop counters and
running residues, invoking input-output when the corresponding input gate
is local. An inverse-output procedure computes (8), invoking middle-output
when needed. This call graph has constant depth. It does not recursively
reconstruct a growing chain of preceding gates.

The carry gates in one group form one consecutive segment of `AL` positions.
Start from its prescribed incoming carry, or zero for the first group.
At each position compute the required coefficient by inverse-output or
read it from the boundary, check (10), and retain only the next carry.
Read one `d`-bit digit of `N` directly from the input tape. At the end check
the outgoing carry or the global final zero condition. No `n`-bit copy of
`N` or of a candidate factor is needed.

Thus a table query takes polynomial bit time and

```math
O(n/P+Af+\log n)=O(n/P)
\tag{14}
```

workspace. The last equality follows from `P<=A` and `A^2f=O(n)`.
For a separately supplied bit-weight routine, add its workspace and
multiply by its query time; only one such routine need be active at once.
At the largest allowable `P`, (14) is `O(sqrt(n log n))`.

## 7. Scalar phase controls fit the existing carry messages

Set `r=log_2 H+d`, an integer. By (9), every legal carry is less than
`HB=2^r`. Since `q>D=T B^2`,

```math
f\ge\log_2 T+2d+1,
\qquad f-r\ge d+2.
\tag{15}
```

On an adjacent-group edge, the existing carry slot can therefore hold
`carry+2^r state`, for any control state occupying at most `d+2` bits.
Decode the carry before using it in the original local arithmetic checks,
and reject invalid control codes. No clique edge becomes wider.

In particular, the scalar phase path from
[output-sensitive-space-transfer.md](output-sensitive-space-transfer.md)
fits. Accumulate each group's phase contribution while its first-factor
digits are reconstructed; require the outgoing state to equal the incoming
state plus that contribution. Start at zero and place the prescribed scalar
phase weight at the last group. Each original factor pair has exactly one
control path. The query time remains polynomial and the additional control
workspace is constant, besides the phase-bit generator.

There is also a direct bucket version when a certified support bound
`tau(N)<=K<=n` is available. Choose a binary control field with
`2^k>n(K-1)`; one can take `k<=2m+1<=d+2`. The entire hash accumulator
`sum_j a_j t^j` then fits in the spare carry bits. Use field addition for
group transitions and the scalar terminal indicator that the final hash
equals the selected bucket. A bit pin gives its paired coordinate sum.
These are direct scalar counts, with no roots of unity or extension of the
contraction's coefficient field. The separate binary control field has
polynomial size and a polynomial-time exhaustive irreducible-polynomial
construction. This support promise changes the recovery interface; it is
not a test for the promise.

## 8. Consequence and limitations

Let a uniform scalar clique evaluator on written read-only tables have
time `E(P,lambda)` and workspace `W(P,lambda)`. Combining this constructor's
phase queries with the virtual rough-cofactor recovery gives, for any
`1<=s<=n`,

```math
\begin{aligned}
T_{\rm fact}(n)
&\le(E(P,\lambda)+1)\mathrm{poly}(n)
 \exp(O(s+n/s)),\\
W_{\rm fact}(n)
&\le W(P,\lambda)+O(n/P+s+n/s+\log n),
\end{aligned}
\tag{16}
```

for the explicit `lambda` in (13), over any fixed coefficient field.
Each requested table symbol is regenerated, and that polynomial cost
multiplies the evaluator time. Table addresses have `O(n/P+log n)` bits
and are included. The evaluator's workspace remains an explicit hypothesis.
If coefficient encodings have `b` bits, charge their arithmetic and encodings
in the same way as the varying-field version of that recovery note.

Choosing `P=A` and `s=ceil(sqrt(n))` makes the additional workspace
`O(sqrt(n log n))`, and the additional logarithmic time `O(sqrt(n))`
up to polynomial factors. This is a constructor/recovery overhead statement.
There is no asserted evaluator with a matching small running time at these
parameters. The established worst-case factoring exponent remains `1/2`.

## 9. A CRT refinement and a variable transform aspect ratio

The same construction has a useful refinement: build an arithmetic modulus
from several small primes instead of searching directly modulo `D`. This
keeps prime searches polynomial in `n` even when `d` is much larger than
`log n`, and gives a sharper information-volume bound.

Choose any power of two `A>=2` with `A^2 m<=n`, retain
`T=A^2`, `H=T/2`, and put `d=ceil(n/H)`, `B=2^d`, `D=T B^2` as before.
The largest-power choice from Section 1 is one allowed option. In this
section, `D` can be exponentially large as an integer; its bit length is
still `O(n)`.

Enumerate distinct primes `q_i=1 mod T` in increasing order and stop when
their product `M` first exceeds `D`. This search has polynomial time in
`n`, for the following reason. Let `J` be the least positive integer with
`(T+1)^J>D`; then `J=O(n)`. Let `R` be the least power of three at least
`2J+3`, so `R=O(J)` and `phi(R)>=J`. Among the residue classes modulo
`TR`, there are `phi(R)` unit classes congruent to one modulo `T`.
Effective Linnik gives a distinct prime in every such class below

```math
V=C(TR)^L=n^{O(1)}.
\tag{17}
```

These primes are distinct because their classes modulo `TR` differ.
There are therefore at least `J` primes congruent to one modulo `T` below
`V`. Each exceeds `T`, so the product of the first `J` already exceeds
`D`. The actual search does not require the constants in (17) or a list of
those classes; the argument supplies its uniform polynomial bound. Trial
division, increasing prime candidates, and one prime-root search at a time
suffice.

At each prime, construct an exact `T`th root `omega_i` by the nonresidue
procedure in Section 1. Maintain only the current product `M`, its root
`omega`, and the current prime/root. When adding a new prime `q`, combine
the roots by

```math
\omega'=\omega+M\bigl((\omega_q-\omega)M^{-1}\bmod q\bigr),
\qquad M'=Mq.
\tag{18}
```

The primes are distinct, so this inverse exists. Formula (18) is already
in `[0,M')`. There is no stored list of earlier primes. At termination,
the previous product was at most `D`, and the last prime is at most `V`.
Thus the total arithmetic word length satisfies

```math
f=\lceil\log_2 M\rceil
\le 2d+(L+1)\log_2 T+L\log_2 R+O(1)
=O(d+\log n).
\tag{19}
```

All prime searches use `O(log n)` bits, while the product and CRT-root
updates use `O(d+log n)` bits. There are at most `J=O(n)` updates, all on
polynomially many bits. This proves polynomial total bit time and
`O(d+log n)` workspace, even though the final modulus itself need not be
polynomially bounded as an integer.

Use `Z/MZ` as the arithmetic ring and replace `q` by `M` in the transform
and range checks. Equations (6)–(8) give the correct DFT separately modulo
every `q_i`, hence modulo `M` by the Chinese remainder theorem. Because
`M>D`, the resulting convolution coefficients are still the unique
canonical integers in (5). The carry proof is unchanged.

Interpolation remains valid over this composite ring. At every prime
factor `q_i`, the root has exact order `T`, so `xi=omega^A` has exact
order `A`. For distinct `k,l`, `xi^k-xi^l` is nonzero modulo every `q_i`
and therefore is a unit modulo `M`. The Vandermonde determinant is a
product of such units. All required Lagrange denominators can consequently
be inverted by the integer Euclidean algorithm. Exact global order of a
root modulo a composite modulus alone would not establish this condition;
the componentwise root construction does.

The same topology, grouping, and regeneration proof now give

```math
P^2\lambda=O(n),\qquad
\text{largest gate mass}=O(n/A),\qquad
\text{table-query workspace}=O(n/P),
\tag{20}
```

for every power of two `2<=P<=A` and every allowed `A`. Indeed,
`T d<2n+T`, `T log n=O(n)`, and
`A(d+log n)=O(n/A)` under `A^2 m<=n`. The metadata cost fits the same
workspace bound. The carry/control packing in Section 7 also survives,
because `M>D` gives the same spare-bit inequality (15).

For the largest-power choice of Section 1, `d=Theta(log n)` and
`J,R=O(1)` asymptotically. Equation (19) then gives

```math
A^2f\le(L+5+o(1))n,
\qquad P^2\lambda\le(9L+45+o(1))n.
```

Here `2A^2 d<=4n+O(A^2)` and `A^2 log_2 A^2<=n` provide the two main
terms. These are upper bounds, not predictions of practical prime sizes.

A slightly larger digit size improves the asymptotic constant further.
Put `g=ceil(log_2(m+2))` and choose `A` as the largest power of two with
`A^2 m g<=n`. Then `d=Theta(mg)`, `J,R=O(g)`, and (19) yields

```math
A^2f\le(4+o(1))n,
\qquad P^2\lambda\le(36+o(1))n.
\tag{21}
```

The top-query workspace at `P=A` is now
`O(sqrt(n log n log log n))`. Prime and root searches remain polynomial
in `n`, and metadata use `O(log n log log n)` bits. This is an explicit
tradeoff between the information-volume constant and arithmetic workspace;
the logarithmic factoring saving exponent remains unchanged.

## 10. Group carry columns to reduce crossing volume

When `P` is smaller, a different carry placement improves the leading
constant. Keep the four `A`-indexed gate families in the same groups, but
place carry `j=iA+r` in the group containing `r`, rather than the group
containing `i`. Every inverse-to-carry edge is now internal. The group
still contains only `L<=h` middle gates, so the interpolation uniqueness
proof is unchanged.

For distinct groups `a,b`, the exact crossing-message count becomes

```math
6L^2+A[|a-b|=1]+(A-1)[\{a,b\}=\{0,P-1\}].
\tag{22}
```

The first two input layers and the middle/inverse layer contribute `6L^2`.
Within each of the `A` carry rows, there is one edge at each adjacent
group boundary. Between rows, there are `A-1` wrap edges from the last
group to the first. For `P=2`, the adjacent and wrap contributions lie
on the same unordered pair; for `P>=4`, they lie on different pairs.

Define

```math
e_P=\begin{cases}2A-1,&P=2,\\A,&P\ge4,\end{cases}
\qquad \lambda_{\rm col}=(6L^2+e_P)f.
```

After deterministic padding,

```math
P^2\lambda_{\rm col}=6A^2f+P^2e_Pf.
\tag{23}
```

Each group now contains `A` disjoint carry segments, each of length `L`.
Process those segments one at a time. Their initial carries are prescribed
boundary values, except for the global initial zero. At the end of each
segment, check the outgoing boundary carry or the global final zero.
All coefficient edges come from inverse gates inside the same group, whose
outputs are regenerated by the existing constant-depth routines. No carry
segment calls another segment recursively.

The table query still has polynomial time. Its workspace is
`O(A^2f/P+PAf+Af)`. If `P^2<=A`, this is `O(n/P)` with either arithmetic
modulus construction. For `P>=4` in that range, (23) is at most `7A^2f`;
the `P=2` case is at most `8A^2f`. The scalar control path also fits:
use the existing carry edge in row zero between each pair of adjacent
groups, and pack its control state as in Section 7.

With the slowly larger digit policy of Section 9, `A^2f<=(4+o(1))n`.
For any sequence of allowed powers of two with `P=o(sqrt(A))`,

```math
P^2\lambda_{\rm col}\le(24+o(1))n.
\tag{24}
```

Indeed, `P^2e_P/A^2=o(1)`. This parameter range includes every fixed
power `P<=n^delta` with `delta<1/4`, for sufficiently large `n`.
The formula is an improved top-constructor volume bound; it supplies no
new crossing identity or larger logarithmic saving exponent.

## Checks and attribution

Run the companion finite checker:

```sh
python3 tools/check_small_word_fourier_constructor.py \
  --output evidence/small-word-fourier-constructor.json
```

The checker compares the two-dimensional transform with a separate full DFT
and an integer convolution, checks carries and factor-pair uniqueness on small
instances, checks interpolation and grouped edge counts, and checks the
carry/control codec. It also checks CRT roots, unit Vandermonde differences,
composite-ring interpolation, and the variable parameter schedules. A modulus
15 control shows that a root's exact global order alone does not suffice
over a composite ring. The carry-column variant has a separate edge-count
and segment-layout check. Its small arithmetic parameters test the identities;
they do not benchmark the asymptotic Linnik search or establish its theorem.

The transform topology and grouping follow the Rudich–Rousseau multiplication
network and the repository's existing ordinary top-query proofs. The prime
selection bound is the cited effective Linnik theorem. This note used
substantial OpenAI assistant support. No Lean verification is claimed.
