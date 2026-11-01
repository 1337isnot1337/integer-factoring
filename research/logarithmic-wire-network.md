# A multiplication network with logarithmic wire widths

## Results and scope

The multiplication network can use logarithmic-size arithmetic messages
while retaining linear total information volume. This removes the
square-root-size residue modulus from the original constructor.
It also gives a polynomial-time Fourier primitive for an entire
arithmetic block.

For an `n`-bit integer, the construction below has

```math
A=\Theta\!\left(\sqrt{\frac n{\log(n+2)}}\right),\qquad
f=O(\log(n+2)),\qquad A^2f=O(n).
```

Here `A` is the number of middle gates and `f` is every wire's
bit width. For every power of two `2<=P<=A`, an explicit grouping
produces a clique with

```math
\lambda=\bigl(12(A/P)^2+1\bigr)f,\qquad
P^2\lambda\le13A^2f=O(n).
\tag{1}
```

Metadata and each grouped table query take polynomial bit time and
polynomial space. Every unweighted entry is zero or one; a weighted
entry is zero or one first-factor monomial.

In a separate Fourier representation, all local gate queries take
polynomial time over an explicitly constructed coefficient algebra of
polynomial dimension. More strongly, any block made only of middle,
interpolation, and coefficient-addition gates has a polynomial-time
Fourier boundary query, with no restriction on the number of middle
gates. Its contraction reduces to bilinear Gauss sums and linear
equations over the arithmetic prime field.
An arbitrary carry-only block also has a polynomial-time Fourier query:
a streamed interval-prefix recurrence contracts each carry segment.

These are constructive improvements to the initial instance and to a
restricted contraction class. They do not supply a contraction rule
for blocks containing arbitrary digit-interval or carry weights.
The known half exponent is preserved; no exponent above one half is
claimed.

The source inspiration is the three-layer multiplication transform in
[the pinned factoring reference](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean).
This note replaces its Fourier transform by bivariate evaluation and
interpolation. Its bijection and costs are proved below as ordinary
mathematics; no existing Lean theorem is asserted to verify this new
network automatically.

## 1. Parameters and an elementary prime search

Fix the coefficient characteristic `p`. The arithmetic prime `M`
used on wires is independent of that characteristic. Finitely many
inputs with `n<p` can be handled directly. For `n>=max(2,p)`,
put

```math
d=\lceil\log_2(n+2)\rceil,\qquad B=2^d.
```

Choose the least power of two `A>=2` for which

```math
h=A/2,\qquad m=Ah=A^2/2,\qquad md\ge n.
\tag{2}
```

Set

```math
H=m(B-1)^2+2A+p,
\qquad M=\text{the first prime greater than }H,
\qquad f=\lceil\log_2M\rceil.
\tag{3}
```

The elementary prime-supply proof in
[rough-input-transfer.md](rough-input-transfer.md), Section 2, gives
at least `H+1` odd primes at most `4(H+2)^2`. At most `H`
of them can be at most `H`, so

```math
H<M\le4(H+2)^2.
\tag{4}
```

Thus sequential trial division constructs `M` in polynomial time
and polynomial space in `n`. No least-prime-in-a-progression theorem,
unproved prime search, or primitive-root search is required.

Minimality in (2), including the finite `A=2` case, gives

```math
A^2d<8n,\qquad m<4n/d.
\tag{5}
```

Since `B<4n`, `d>=2`, and `p<=n`, equations (3)--(5)
give `H+2<64n^3` and, conservatively, `f<=16d`. Therefore

```math
13A^2f<1664n.
\tag{6}
```

All arithmetic messages are elements of `F_M`, represented by their
integers in `[0,M)`. Carry values also use this range.

## 2. A bivariate convolution computes the integer product

Let `a,b<B^m` have base-`B` digits `a_j,b_j`. Encode them as

```math
P(X,Y)=\sum_{r=0}^{A-1}\sum_{i=0}^{h-1}a_{r+Ai}X^rY^i,
\qquad
Q(X,Y)=\sum_{r=0}^{A-1}\sum_{i=0}^{h-1}b_{r+Ai}X^rY^i.
\tag{7}
```

Their product has coefficients `C_(u,v)` with
`0<=u<=2A-2` and `0<=v<=A-2`. Pad all other
coefficients by zero. Substituting `X=Z,Y=Z^A` shows that the
ordinary univariate product coefficient at `j=qA+r` is

```math
c_{qA+r}=C_{r,q}+[q>0]C_{r+A,q-1},
\quad 0\le q,r<A.
\tag{8}
```

Every `c_j` is at most `m(B-1)^2<M`. Each `C_(u,v)`
is also below `M`. Thus the field computations below recover the
actual integer coefficients, not only congruence classes.

Choose the `A` distinct evaluation points
`y_k=k in F_M`, `0<=k<A`. They are distinct because `M>A`.
The network has the following gates.

| Gate family | Number | Setting and output |
|---|---:|---|
| `FA_P[r], FA_Q[r]` | `2A` | Store `h` base-`B` digits and send their polynomial's value at each `y_k`. |
| `MID[k]` | `A` | Store two length-`A` incoming vectors; multiply the two degree-less-than-`A` polynomials in `X` and send all `2A` coefficients. |
| `IB[u]` | `2A` | Store `A` incoming evaluations; interpolate a degree-less-than-`A` polynomial in `Y` and send its coefficients. |
| `ADD[j]` | `A²` | Store two field elements and send their field sum, implementing (8). |
| `CARRY[j]` | `A²` | Store the integer coefficient and incoming carry; check the required product digit and send the outgoing carry. |

For the second input to `ADD[qA+r]`, gate `IB[A+r]` sends
its coefficient of degree `q-1` when `q>0` and sends zero
when `q=0`. This harmless zero-message edge keeps the geometry
uniform. The first input comes from degree `q` of `IB[r]`.
In particular the padded gate `IB[2A-1]` sends zeros in every
globally valid setting.

An input gate computes

```math
FA_r(k)=\sum_{i=0}^{h-1}a_{r+Ai}y_k^i.
\tag{9}
```

At point `y_k`, the middle gate computes exactly
`P(X,y_k)Q(X,y_k)`. Its output to `IB[u]` is
`sum_v C_(u,v)y_k^v`. Degree is at most `A-2`, so
interpolation at the `A` points recovers every coefficient.
Equation (8) then recovers each integer `c_j`.

The carry gate stores `(c,t) in [0,M)^2`, checks

```math
(c+t)\bmod B=\operatorname{digit}_B(N,j),
```

and sends `floor((c+t)/B)`. This output is less than `M`,
because `c+t<2M` and `B>=2`. The first incoming carry is
zero, and the last outgoing carry must be zero.

### Exact weighted witness bijection

The directed network is acyclic. Once the two input digit sequences
are fixed, every middle, interpolation, addition, and carry setting is
determined. All intermediate integer convolution coefficients are
below `M`, so no false integer product can arise from modular wraparound.

The carry equalities telescope to `ab=N`. Conversely, every divisor
`a|N` has `a,N/a<=N<2^n<=B^m` and gives a unique valid
network setting. The carry recurrence stays in `[0,M)` and has
the required final zero carry because `ab<N+1<=B^(2m)`.
Thus globally valid settings are in bijection with the positive
divisors of `N`.

Give each first-factor input bit its formal variable weight. The
network contraction is exactly the bit-multivariate polynomial
`C_N` used in
[bit-tagged-witness-recovery.md](bit-tagged-witness-recovery.md).
Optional prefix restrictions are checked on the first-factor digits
as in the original network. All of these are polynomial-time local
weight queries over any commutative coefficient algebra.

## 3. Unique grouped extensions and exact geometry

Any cluster with at most `h=A/2` middle gates has at most one
valid setting for each prescribed boundary.

Indeed, each contained input gate has at least `h` external middle
neighbors. Their boundary values give `h` evaluations of its
degree-less-than-`h` polynomial. Vandermonde interpolation over
`F_M` recovers its coefficients. Reject any recovered coefficient
outside `[0,B)`, check all remaining boundary values and prefix
conditions, then propagate middle, interpolation, addition, and carry
settings in directed order. Check every local weight and every outgoing
boundary message.

Any valid setting must follow these steps, which proves uniqueness.
Passing all checks gives a valid setting, which proves correctness of
the query procedure. Its arithmetic, data lengths, and number of gates
are polynomial in `n`. Weighted entries require only polynomially
many additional coefficient-algebra operations.

For a power of two `2<=P<=A`, write `L=A/P` and partition
`[0,A)` into contiguous blocks of length `L`. Group
`FA_P[r],FA_Q[r],MID[r],IB[r],IB[A+r]` by the block of `r`.
Group `ADD[qA+r]` and `CARRY[qA+r]` by the block of `q`.
Each group contains `L<=h` middle gates, so the query lemma applies.

For distinct groups, the exact edge counts are:

| Edge family | Edges between groups `a!=b` |
|---|---:|
| The two input families to middle | `4L²` |
| Middle to the two interpolation families | `4L²` |
| The two interpolation families to addition | `4L²` |
| Addition to carry | `0` |
| Carry chain | `1` for consecutive groups, otherwise `0` |

Therefore the pair width is

```math
Y_{ab}=\bigl(12L^2+[|a-b|=1]\bigr)f.
\tag{10}
```

Concatenate messages in a fixed common order and pad to the `lambda`
of (1), rejecting nonzero padding. The usual exact coarsening identity
then preserves the weighted contraction. Equation (6) gives a uniform
linear information-volume bound.

There are `2A²+5A` original gates and `8A²-1` original edges.
Their metadata are generated explicitly. Tables are accessed by
queries; writing every table still has its full exponential output cost.

## 4. Character algebra without field construction

Let the coefficient characteristic be the fixed prime `p`. Since
`M>p`, define the commutative algebra

```math
\mathcal R=
\mathbb F_p[Z]/(1+Z+\cdots+Z^{M-1}),
\qquad \zeta=Z\bmod(1+\cdots+Z^{M-1}).
\tag{11}
```

Its dimension is `M-1=poly(n)`. Polynomial arithmetic constructs it
and implements its operations in polynomial bit time and space.
No irreducible polynomial or primitive element is required.

We have `zeta^M=1`. Because `M` is prime, multiplication by a
nonzero residue permutes its exponents. Hence

```math
\sum_{x=0}^{M-1}\zeta^{ax}
=\begin{cases}M,&a=0\pmod M,\\0,&a\ne0\pmod M.\end{cases}
\tag{12}
```

The nonzero scalar `M in F_p` is a unit. Consequently the
normalized character transform on every `M`-valued wire is
invertible. These identities also hold after adjoining the square-zero
bit tags, or any other commutative coefficient algebra.

The quotient in (11) may have zero divisors. The argument uses (12)
and the invertibility of `M`, not a field assumption about `R`.
Using the all-ones polynomial with a composite arithmetic modulus
would not give (12) in general; primality in (3) is substantive.

## 5. A rectangular bilinear Gauss contraction

Let `J` be any `r by s` matrix over `F_M`, and let
`alpha in F_M^r, beta in F_M^s`. In `R` consider

```math
G_J(\alpha,\beta)
=\sum_{a\in\mathbb F_M^r}\sum_{b\in\mathbb F_M^s}
\zeta^{\alpha^Ta+\beta^Tb+a^TJb}.
\tag{13}
```

Let `k=rank J`. The value is zero unless both

```math
Jb_0=-\alpha\text{ is solvable},
\qquad \beta\in\operatorname{im}J^T.
\tag{14}
```

When (14) holds,

```math
G_J(\alpha,\beta)=M^{r+s-k}\zeta^{\beta^Tb_0}.
\tag{15}
```

To prove this, sum first over `a`. Equation (12) gives zero unless
`Jb=-alpha` and gives `M^r` otherwise. The solution set is
`b0+ker J`. Its character sum is zero unless `beta` annihilates
`ker J`, equivalently `beta in im J^T`. In that case every
kernel term contributes one and there are `M^(s-k)` of them.
The phase is independent of the chosen solution `b0`.

Gaussian elimination over the arithmetic field `F_M` decides
(14), computes a solution and the rank, and evaluates (15) in
polynomial time and space in `r+s+log M`, plus the explicitly
polynomial coefficient-algebra cost. No enumeration of `M^(r+s)`
states is used.

For a transformed middle gate, its output frequencies `gamma_u`
give the Hankel matrix

```math
J_{rs}=\gamma_{r+s}.
\tag{16}
```

Its raw Fourier entry is (13), with `r=s=A`. Rank-deficient
Hankel matrices and all zero-frequency cases are included in (14)--(15).
Transform normalization factors, where required by a chosen edge
convention, are powers of the unit `M` and must be retained.

## 6. Polynomial queries for an entire arithmetic block

Take any set of middle, interpolation, and addition gates, and expose
all its incident boundary wires. Fourier transform the boundary tensor,
using characters on incoming values and opposite characters on outgoing
values, with a fixed consistent convention.

Every incoming boundary value of a middle gate is a free arithmetic
variable. Its two vectors are denoted `a_k,b_k`. Every other
incoming boundary value is an additional free variable `w`.
The internal directed equalities determine every other setting:
middle outputs are bilinear in `a_k,b_k`, while interpolation and
addition outputs are linear combinations of those outputs and `w`.

Thus the phase of a boundary assignment has the form

```math
\sum_k\bigl(\alpha_k^Ta_k+\beta_k^Tb_k+a_k^TJ_kb_k\bigr)
+\delta^Tw.
\tag{17}
```

All coefficients are obtained by polynomial-time linear propagation
through the block. Crucially, no output of a middle gate feeds another
middle gate in this network. Therefore no higher-degree monomial is
created. Different middle gates have disjoint free incoming vectors.

Summing each `w` coordinate gives a linear orthogonality test.
Each middle pair is then evaluated by (14)--(15). Their product,
together with the linear cardinality and normalization factors, is the
entire Fourier boundary query. This proves the block primitive, including
disconnected blocks and zero outputs.

Input gates with digit intervals and carry gates impose additional
weights that are not included in (17). An arbitrary block containing
them need not have a bilinear phase, so this proof makes no closure
claim for the full weighted multiplication network.

### A whole carry block has a streamed polynomial-time query

Consider a consecutive segment of `L` carry gates. At position `j`,
let `m_j` be the prescribed base-`B` digit of `N` and
`alpha_j` the coefficient-input boundary frequency. Denote the
incoming and outgoing carry by `t,t' in [0,M)`. Their coefficient
input, when legal, is uniquely

```math
c=m_j+Bt'-t.
```

Thus the raw character transfer has matrix entries

```math
T_j(t,t')=
[0\le m_j+Bt'-t<M]\,
\zeta^{\alpha_j(m_j+Bt'-t)}.
\tag{17a}
```

Start the state vector with the incoming boundary character, or with
the vector supported at zero if the segment contains the global first
carry. At each step update
`v_(j+1)(t')=sum_t v_j(t)T_j(t,t')`. Finish by taking its inner
product with the outgoing boundary character, or its zero coordinate
if the segment contains the global last carry. This sums precisely
the legal settings of the segment. Signs and normalization factors
are assigned according to the same boundary-transform convention.

The update needs only `O(M)` ring operations, rather than a dense
matrix product. Put

```math
w(t)=v_j(t)\zeta^{-\alpha_jt},
\qquad
P(s)=\sum_{0\le t<s}w(t).
```

For a proposed `t'`, write

```math
\ell=\max(0,m_j+Bt'-M+1),\qquad
u=\min(M-1,m_j+Bt').
```

If `ell>u` the new coordinate is zero. Otherwise

```math
v_{j+1}(t')
=\zeta^{\alpha_j(m_j+Bt')}
\bigl(P(u+1)-P(\ell)\bigr).
\tag{17b}
```

All phases in the prefix and output loops can be generated by repeated
multiplication by fixed powers of the unit `zeta`. Constructing
the prefixes and all new coordinates therefore costs `O(M)` ring
operations and `O(M)` live ring elements. A segment of length `L`
costs `O(LM)` ring operations, plus polynomial control work.
The coefficient-algebra representation in Section 4 makes this a
polynomial bit-time and polynomial-space bound.

An arbitrary set of carry gates is a disjoint union of consecutive
segments. Its Fourier boundary tensor factors over those segments,
so the same primitive handles every carry-only block. In particular
the entire carry chain can be contracted in polynomial time for fixed
coefficient boundary frequencies.

This lemma and the arithmetic-block lemma concern separate boundary
queries. Contracting their common boundary still requires a further
algorithm; neither lemma evaluates that shared-label sum automatically.

### Other local gates

An input-gate character query factors over its `h` digits into
finite geometric sums. Prefix restrictions and independent bit weights
replace each digit sum by a product of factors
`1+beta_j z^(2^j)`, with its fixed high-bit multiplier.
Bit-tagged coordinates can be computed with prefix and suffix products.

Interpolation and addition gates are linear maps, so (12) gives
their queries by solving their explicit linear frequency constraints.
Carry gates have only `M²=poly(n)` settings, so even direct
enumeration is a polynomial-time query. The interval formulas in
[fourier-odd-characteristic.md](fourier-odd-characteristic.md),
Section 4, give a more compact division-free implementation and apply
with both value and carry domain equal to `M`.

Thus every local Fourier query, arithmetic-block query, and carry-block query has a
polynomial bound including coefficient representation and preprocessing.
This is stronger than an arithmetic-operation count that ignores a
large character field.

## 7. Compatibility with the half-rate streamed evaluator

The simplest use is to cap the top schedule at any fixed
`P<=n^c` with `c<1/2`. An adjusted sublinear-volume cutoff
then gives polynomially large labels exactly as in
[subexponential-piece-transfer.md](subexponential-piece-transfer.md).
All initial metadata and table queries are polynomial.

The existing modular family has only polynomially many pieces, so it
also permits the largest scheduled `P=K_ell<=A`. Here is the
additional accounting needed when the top label can be only logarithmic.

Stop recursion at volume `V<=n^(3/4)`, using `V^4<=n³`.
The top label satisfies `lambda>=f=Omega(log n)`. After the
first grid step, `K_(ell-1)^2<=P<=A`. Hence every later node
that performs a grid step has

```math
\lambda_j>\frac{n^{3/4}}{P}
=\Omega(n^{1/4}\sqrt{\log n}).
\tag{18}
```

Only the first step can have a merely logarithmic label. Its relative
width slack `O(1/loglog lambda)` is bounded and tends to zero.
The sum of all later relative width slacks is bounded by the existing
`O(loglog n)` depth argument. Rounding contributes at most
`P/n^(3/4)=O(n^(-1/4))` at each grid step. Thus the same
constant-factor volume envelope and parent-query series apply.

For polynomially many pieces, the top piece-assignment entropy is

```math
O\!\left(P^2\log\left(\frac{Cn}{P^2}+2\right)\right)
=O\!\left(\frac{n\log\log n}{\log n}\right),
\tag{19}
```

because `P²<=A²=O(n/log n)`. This is below the half-rate
scale. Later levels have the polynomial label lower bound (18), so
their total piece entropy is also lower order by the established
geometric-volume estimate. The grouping searches cost at most
`exp(O(n^(1/4)(log n)^2))`, and early leaves cost
`exp(O(n^(3/4)))`; both are lower order.

The schedule still has square growth. Since its largest level below
`A` has `log P=Theta(log n)`, it still gives
`rho^ell=Theta((log n)^(-log2(1/rho)))`. The modular
`rho=1/sqrt(2)` family therefore gives the exact half endpoint,
with the polynomial-space execution and recovery already proved in the
companion notes.

Equation (19) uses the **polynomial** piece bound. A generic
subexponential piece bound from the relaxed interface does not
automatically justify the largest top size; in that setting one should
retain a fixed-power cap `P<=n^c` with `c<1/2`.

## 8. Verification

The companion checker compares the bivariate pipeline with ordinary
integer convolution, checks prefix and weighted witness recovery,
verifies the exact grouping counts and information bounds, and compares
small grouped queries with exhaustive input-digit enumeration.
It also checks (13)--(15) against direct finite sums, including
rectangular and singular matrices, contracts a complete arithmetic core
from its unrestricted boundary values, compares the carry-prefix
recurrence with exhaustive short-chain settings, and tests the character
algebra in a case with zero divisors.

The finite checks validate the transcription. The proofs above supply
the uniform bounds and preserve the distinction between the arithmetic
block primitive and a full-network contraction algorithm.
