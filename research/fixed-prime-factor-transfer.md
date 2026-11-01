# Factoring from scalar counts in any fixed prime characteristic

## Result and scope

The scalar-parity transfer extends to every fixed prime characteristic.
For each prime `p`, explicit multiplication predicates have the following
properties: every composite input has a nonzero witness count modulo `p`,
every recovered witness gives a proper factor, and the corresponding
networks have total mass `O_p(n)` for an input of at most `n` bits.

The useful feature is this last, quantitative bound. A generic reduction
through prime-factorization certificates would not by itself preserve
the input-length scale needed for the factoring exponent. Here the
predicates require only a fixed number of factors, repeated `p`th
powering, and comparisons. Their combined width, local settings,
copying, and complete table construction are all charged below.

**Fixed-prime scalar transfer.** Fix a prime `p`. Suppose one uniform
deterministic finite-alphabet algorithm evaluates every written
`P`-vertex clique over `F_p`, with `lambda` bits on each edge and
arbitrary `F_p`-valued vertex tables, in time `E_p(P,lambda)`.
This time includes reading the full input.

There is a constant `C_p` such that, for `2<=N<2^n` and any integer
`P>=2` with `P^40<=n`, complete prime factorization and the smallest
prime factor can be computed in time

```math
\operatorname{poly}_p(n)\left(
  2^{C_pn/P}
  +E_p\!\left(P,\left\lfloor\frac{20000p^2n}{P^2}\right\rfloor\right)
\right).
\tag{1}
```

The field is a fixed part of the algorithm. Choosing `P` in polynomial
time gives one uniform factoring algorithm. Small lengths outside a
chosen schedule can be treated by direct trial division.

Only scalar evaluation over `F_p` is assumed. There is no need to
evaluate over polynomial rings or other `F_p`-algebras, reconstruct
an integer count, or run a probabilistic isolation procedure.
This is an improvement to the transfer interface. It does not give
a new evaluator or an achieved exponent above the current one half.

The quantitative construction uses the free-output multiplication
gadget and grouping argument in
[Scalar parity transfer](parity-factor-transfer.md), which identifies
the exact prior multiplication-network sources. This note extends
the witness family and proves the additional product, copy, comparison,
and mass bounds needed in characteristic `p`.

## 1. Two explicit witness families

Write, for the proof only,

```math
N=\prod_q q^{e_q},\qquad e_{\max}=\max_q e_q.
```

The integer primes `q` dividing `N` need not equal the fixed
coefficient characteristic `p`.

### Power-root tuples

For a positive level `j` with `r=p^j<=n-1`, let

```math
R_{p,j}(N)=
\left\{(a_1,\ldots,a_{p-1})\in\mathbb Z_{>0}^{p-1}:
 (a_1\cdots a_{p-1})^r\mid N,\quad
 (a_1,\ldots,a_{p-1})\ne(1,\ldots,1)
\right\}.
\tag{2}
```

Every tuple in (2) contains a proper divisor greater than one.
Indeed, its product `A` satisfies `2<=A<=N^(1/r)<N`, and every
nonunit coordinate divides `A` and hence `N`.

For each integer prime `q`, distributing exponents among the `p-1`
coordinates gives

```math
|R_{p,j}(N)|+1
=\prod_{q\mid N}
 \binom{\lfloor e_q/r\rfloor+p-1}{p-1}.
\tag{3}
```

The binomial coefficient counts nonnegative exponent tuples with
sum at most `floor(e_q/r)`; an additional slack coordinate gives
the usual weak-composition formula. The subtracted tuple is the
unique all-ones tuple.

### Cyclic representatives of factor tuples

Consider the ordered positive `p`-tuples

```math
\mathcal F_p(N)=\{(x_0,\ldots,x_{p-1}):x_0\cdots x_{p-1}=N\}.
```

The cyclic group of order `p` acts by rotating coordinates. Since
`p` is prime, a tuple fixed by any nonidentity rotation is constant.
Every nonconstant orbit therefore has exactly `p` elements and
has a unique lexicographically least tuple.

Let `C_p(N)` be the set of tuples that are **strictly** smaller than
each of their `p-1` nonidentity rotations, with the tuple

```math
(1,\ldots,1,N)
\tag{4}
```

removed. Strictness excludes constant tuples. The tuple in (4) is
the least representative of the orbit with just one nonunit factor.
Every remaining tuple has at least two nonunit coordinates. Each
such coordinate is a proper divisor of `N`.

The total number of ordered tuples is

```math
D_p(N)=|\mathcal F_p(N)|
=\prod_{q\mid N}\binom{e_q+p-1}{p-1}.
\tag{5}
```

There is one constant tuple exactly when `N` is a perfect `p`th
power. Consequently, for every `N>=2`,

```math
|C_p(N)|
=\frac{D_p(N)-[N\text{ is a perfect }p\text{th power}]}p-1.
\tag{6}
```

Division by `p` in (6) is an integer orbit-counting identity used
in the proof. The evaluator never divides a field element by `p`.
It counts the lexicographically selected representatives directly.

## 2. Every composite has a nonzero count modulo `p`

For `1<=e<=p-1`, define

```math
B_e=\binom{e+p-1}{p-1}=\binom{e+p-1}{e}.
```

Exactly one factor in its numerator is divisible by `p`, and its
denominator is not. More precisely,

```math
v_p(B_e)=1,\qquad
\frac{B_e}{p}
=\frac1e\prod_{h=1}^{e-1}\frac{p+h}{h}
\equiv e^{-1}\pmod p.
\tag{7}
```

The rational display has denominator prime to `p`; `B_e/p` itself
is an integer. In particular, `B_e` is zero modulo `p` whenever
`1<=e<p`, while `B_0=1`.

**Nonzero-witness lemma.** If `N` is composite, at least one of
the sets `R_(p,j)(N)` with `p^j<=n-1`, or the set `C_p(N)`, has
cardinality nonzero modulo `p`. For prime `N`, all these sets are
empty.

If `e_max>=p`, choose `j=floor(log_p(e_max))`. This level is
positive, and `r=p^j<=e_max<=n-1`. Every floor `floor(e_q/r)`
lies between zero and `p-1`, and at least one is positive.
Equations (3) and (7) give

```math
|R_{p,j}(N)|\equiv-1\pmod p.
\tag{8}
```

If `e_max<p`, there is no constant tuple in (6), because some
positive exponent is less than `p`. Let `t` be the number of
distinct integer primes dividing `N`.

If `t>=2`, equation (7) shows that `D_p(N)` is divisible by
`p^t`, hence by `p^2`. Therefore

```math
|C_p(N)|=D_p(N)/p-1\equiv-1\pmod p.
\tag{9}
```

If `t=1`, write `N=q^e` with `1<=e<p`. Then (6)--(7) give

```math
|C_p(N)|\equiv e^{-1}-1\pmod p.
\tag{10}
```

This is nonzero for `2<=e<=p-1`, and it is zero for `e=1`,
the prime case. A prime has no tuple in either witness family,
so the converse follows as well.

The algorithm tests the explicit levels and the cyclic predicate.
It does not know the exponents or use them as advice. Levels are
generated by multiplying a counter successively by the fixed `p`.
There are `O_p(log n)` of them.

For `p=2`, the cyclic representatives are exactly pairs `a<b`
with product `N`; removing `(1,N)` gives the lower-half divisor
set. Equations (8)--(10) specialize to the characteristic-two
argument. That specialization can use interval prefixes instead
of the comparison networks below and has smaller constants.

## 3. Recover a proper factor from a nonzero count

Encode each coordinate of a witness tuple by an `n`-bit binary
word and concatenate those words in a fixed order. A power-root
tuple has `(p-1)n` bits, and a cyclic tuple has `pn` bits.
For any fixed prefix, the predicate counts its valid completions
in `F_p`.

Once the root count is nonzero, query the lower binary child.
If its count is nonzero, select it. Otherwise select the upper
child. The two child counts add to the parent count, so the
selected child remains nonzero. After at most `pn` queries,
the remaining word is one valid witness tuple.

In the power-root family, select any coordinate greater than one.
Section 1 proves it is a proper divisor. In the cyclic family,
select any nonunit coordinate; at least one other coordinate is
nonunit, so the selected coordinate is also less than `N`.
The prefix search need not return the least tuple or least factor.

Recursively split both the returned factor and its quotient. If
all root counts vanish, Section 2 certifies that the current
integer is prime. As in the parity construction, the full split
tree has at most `2n-3` nodes, and all calls can use the original
padded bit width `n`. This gives complete factorization with
polynomially many scalar evaluations. Taking the minimum yields
the smallest prime factor.

The implementations below count all power-root tuples, including
the all-ones tuple, and all strict cyclic representatives,
including (4). Subtract one in `F_p` exactly when the current
prefix matches the excluded tuple. Those tuples and their
arithmetic extensions are unique, so this removes exactly one
setting. There is no division in either correction.

## 4. Arithmetic networks have total intended width `O_p(n)`

The free-output multiplier from
[the parity construction, Section 4](parity-factor-transfer.md)
uses the pinned reference's actual digit and carry definitions.
For intended input width `w`, it retains full products and has
the following bounds, including exposed input and output bits:

```math
\text{total mass}\le3200w,
\qquad
\text{maximum vertex mass}\le512w^{3/4}.
\tag{11}
```

Here mass is incident bit width plus a bound on the binary
logarithm of the number of local settings. An exposed bit vertex
can have at most two equality/copy edges in addition to its
arithmetic edge; the resulting mass bound is already included
in (11). Extra copying vertices will be charged separately.

Each valid arithmetic input pair has exactly one internal setting
and exactly its full product at the output. Connecting outputs
to later inputs with bit equality edges preserves uniqueness.
All unused input bits are zero. Every leading output bit of the
final product is checked against zero when it lies above `N`.

### Power-root arithmetic

Fix `r=p^j<=n-1` and `m=ceil(n/r)`. Each root factor `a_i`
fits in `m` bits, since `a_i^r<=N<2^n` for every valid tuple.
First multiply the `p-1` factors to form

```math
A=a_1\cdots a_{p-1}.
```

Allocate `w_0=(p-1)m` bits for this full product. The initial
product chain uses intended multiplier widths
`m,2m,...,(p-2)m`, with no multipliers in this part when `p=2`.

Then perform `j` stages of `p`th powering. At stage `i`, the
input has allocated width `w_i=p^i w_0`. Compute its powers
two through `p` by successive multiplication by that same input.
The intended multiplier widths are
`w_i,2w_i,...,(p-1)w_i`. Finally multiply `A^r` by an `n`-bit
cofactor `b`, using intended width `max(n,rw_0)`, and pin its
full output to `N`.

The sum `W` of the intended multiplier widths is

```math
W=
 \frac{(p-2)(p-1)}2m
 +\frac{p(p-1)}2\sum_{i=0}^{j-1}p^iw_0
 +\max(n,rw_0).
\tag{12}
```

Since `rm<n+r<=2n`, the middle term is at most
`p(p-1)n`, the first is at most `(p-2)(p-1)n/2`, and the last
is at most `2(p-1)n`. Thus

```math
W\le2p^2n,
\qquad \max\{\text{intended multiplier width}\}\le2pn.
\tag{13}
```

All allocated powers are retained. No congruence modulo a truncated
bit width replaces the integer equation `A^r b=N`.

### Cyclic-factor arithmetic

Multiply the `p` independent `n`-bit coordinates in sequence and
pin the full product to `N`. The intended widths are
`n,2n,...,(p-1)n`. Their sum is `p(p-1)n/2`, and the maximum
is `(p-1)n`, so (13) holds for this family too.

For any input tuple, all product and power values are deterministic.
For a power-root tuple with `A^r|N`, the cofactor is uniquely
`N/A^r`. Thus the networks count precisely the ordered tuples in
Section 1, with one setting per tuple before the known correction.
Although binary input ports can encode zero, the final product is
the positive integer `N`, so every accepted factor and cofactor
is positive. No additional positivity gadget is needed.

## 5. Copying and cyclic comparisons cost only `O_p(n)`

To supply one bit to several consumers, use a chain of one-bit
copy vertices. Each vertex stores that bit, sends it to one
consumer, and passes it to the next vertex. The first receives
the original bit. Each copy vertex has degree at most three
and one setting bit, so its mass is at most four. An existing
arithmetic bit port gains at most one connection to such a chain.
Its degree remains within the allowance in (11).

At a `p`th-power stage, at most `p` copies of each input bit
are needed. The number of added copy vertices across all stages
is at most

```math
p\sum_{i=0}^{j-1}w_i
=p m(r-1)<2pn.
\tag{14}
```

Intermediate product bits used once are connected directly.
The extra mass of the power-copy chains is therefore below `8pn`.

For cyclic representatives, compare the concatenated `pn`-bit
tuple with each of its `p-1` rotations. A comparison scans from
the most significant bit with three states: equal, less, greater.
Two state bits suffice. The transition is deterministic: retain
an unequal state, or resolve it at the first differing input pair.
Pin the initial state to equal and require the final state to be
less. This enforces strict lexicographic comparison.

One comparison-step vertex stores two input bits and a two-bit
incoming state, so its setting exponent is at most four. Its
incident width is at most six: two bits for each state edge and
one for each input. Its mass is at most ten. Invalid encodings
of the fourth state have weight zero. Once the tuple bits are
fixed, the whole state sequence is unique.

There are `p(p-1)n` comparison steps. Each tuple bit is used
`2(p-1)` times across these comparisons. Copy chains therefore
need at most `2p(p-1)n` added vertices. Together the comparisons
and their copy chains have mass at most

```math
10p(p-1)n+8p(p-1)n=18p(p-1)n.
\tag{15}
```

Prefix conditions merely pin the independent tuple bits. If a
power-root coordinate's padded prefix is inconsistent with its
`m` allocated bits, its count is zero immediately. Otherwise
the prefix adds no vertices or settings.

All predicates have `0/1` local weights, interpreted in `F_p`.
The arithmetic, copying, and state transitions have unique
extensions; none introduces unwanted multiplicity in the counts.

## 6. Uniform mass and written-clique bounds

Combining (11), (13), (14), and (15), both witness families have

```math
\sum_v(d_v+b_v)\le6500p^2n,
\qquad
\max_v(d_v+b_v)\le1024p\,n^{3/4}.
\tag{16}
```

For the first bound, multipliers contribute at most `6400p^2n`;
the relevant copy and comparison terms fit within the remaining
constant. For the second, their largest intended width is at
most `2pn`, and
`512(2p)^(3/4)<=1024p`. Copy and comparison vertices have
mass at most ten and satisfy the same bound.

The number of vertices, the network description, and each local
message computation have polynomial cost in `n` for fixed `p`.
In particular, the repeated-power stages use geometrically growing
widths, so their total width is the sum in (12), rather than
`n` times the number of stages.

Use the same generic grouping theorem as the characteristic-two
construction. Set

```math
G=3P^8,\quad U_0=6500p^2n,\quad
h=\lceil(n^3)^{1/4}\rceil,\quad
Q=\left\lceil U_0/G\right\rceil+1024p h.
\tag{17}
```

Greedily assign each vertex to the cluster of least current mass.
Equation (16) bounds every cluster's mass by `Q`. Coarsening
therefore gives incident width at most `Q`, ordered total width
at most `U0`, and setting exponent at most `Q` per cluster.

The pinned reference's
[`EngineGroup.good_exists`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineGroup.lean)
gives a grouping into `P` groups, each containing at most
`floor(2G/P)` clusters, with pair width at most

```math
\left\lfloor\frac{P U_0+GQ}{P^3}\right\rfloor.
```

Let `lambda` be the larger of that number and one, and let
`H=floor(2G/P)Q`. As before, `P^40<=n` implies `Gh<=6n`
and `G<=3n`. Hence

```math
GQ\le(6500p^2+6144p+3)n.
```

For all primes `p>=2`, the following convenient common bounds hold:

```math
P^2\lambda\le20000p^2n,
\qquad PH\le20000p^2n.
\tag{18}
```

Use the common width `lambda_star=floor(20000p^2n/P^2)`,
padding all unused bits with zeros. Search through the `P^G`
colorings to find the promised grouping. Its exponent is at most
`3n/P`, as in the parity note.

For each table entry, enumerate local original settings, check
all local weights and prescribed message equalities, and add
one modulo `p` for every accepted tuple. Internal messages are
determined by their endpoints; no extra internal-channel loop
is needed. Coarsening then preserves the exact global count
modulo `p`.

There are `P*2^((P-1)lambda_star)` entries, each obtained from
at most `2^H` setting tuples, with polynomial bit cost per check.
Equation (18) makes complete table generation and output cost
`2^O_p(n/P)`. Each entry is one fixed-size field element.
All integer parameters in (17)--(18) are computable by ordinary
integer arithmetic. There is no real logarithm or uninstantiated
number-theoretic constant in this construction.

## 7. Scalar transfer and its limits

The algorithm has `O_p(log n)` possible power levels and one
cyclic predicate at each factorization node. A nonzero witness
requires at most `pn` prefix steps. There are `O(n)` nodes.
Each query uses one fully written clique, followed by the known
unit correction when its prefix matches the excluded tuple.
The cliques use the same original padded `n`, `P`, and
`lambda_star` throughout the recursive factorization.

Their construction and evaluation costs give (1). The output
also has polynomial length. This proves the complete scalar
transfer with no assumption that the evaluator is an arithmetic
circuit or that its program can be reinterpreted over other rings.

The point of the transfer is that the construction's size is
`O_p(n)` before grouping. Generic polynomial-size witness
encodings alone would not justify preserving a bound such as
`exp(O(n/(log n)^alpha))` in the original input length.

Equation (1) states exactly what a proposed scalar evaluator
must provide, including its written input. If it gives a stronger
asymptotic bound on these instances, this reduction transfers
that bound without a sparse-polynomial support cost. The theorem
does not construct such an evaluator, remove a nonunit crossing
scalar, extend a crossing family's proved parameter range, or
reduce the internal precision of the current Cauchy evaluator.
The achieved factoring exponent remains one half.

## 8. Reproducible controls

Run

```sh
python3 tools/check_fixed_prime_factor_transfer.py \
  --output evidence/fixed-prime-factor-transfer.json
```

The checker uses integer divisor and ordered-tuple enumeration
as a finite reference oracle. It imports the already checked
free-output multiplier mass model from the characteristic-two
checker. This is evidence for the reduction, not an implementation
of a faster clique evaluator.

The recorded run covers characteristics 2, 3, 5, and 7 and verifies:

- 262,136 completeness cases, using independent prime-exponent
  formulas on every integer from 2 through 65,535 in each field;
- 4,092 complete factorizations, every input from 2 through 1,024
  in each field, with padded recursive calls;
- 2,044 direct cyclic-orbit counts and 2,635 power-root counts;
- 290,536 nonzero-prefix queries;
- 41,683 power-network mass controls and 7,996 cyclic-network
  mass controls;
- 248 exact constructor-cutoff checks and 5,460 deterministic
  comparison traces.

Controls also exclude constant tuples and the trivial factor orbit,
and illustrate why the full-orbit argument requires a prime
characteristic: a four-tuple can have a nonconstant orbit of
length two. The evidence records representative prime powers,
semiprimes, and both checker hashes. The general theorem follows
from the integer count identities, unique-extension construction,
and mass proof above.
