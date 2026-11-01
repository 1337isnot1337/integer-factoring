# Weighted binary-energy shells do not rescue the reciprocal pencil

## Result and scope

There is a genuine low-entropy graded state space behind the proposed
power-of-two mode construction. Its degree-`h` dimension is the binary
partition number. However, restricting the characteristic-two reciprocal
three-pairing pencil to that space does **not** give a leading width
contraction, even after arbitrary coefficient recombination.

For `h=2^R`, `R>=1`, let `D_h` be the size of the energy shell defined below.
Over every characteristic-two field, every exact decomposition of its
crossing into pieces from the pencil's full coefficient/symmetry span
has a piece satisfying

```math
r_{\rm row},r_{\rm col}\ge D_h^{\,2-o(1)}.
```

Equivalently, the optimal leading contraction in this span is exactly
`rho=1`. The upper bound follows by using the crossing itself as one
piece. This is a restriction on one polynomial family. It does not cover
other polynomial identities, arbitrary tensors on the shell, or a
different lacunary vertex-operator exchange relation. It gives no lower
bound on integer factoring, and it does not improve the repo's achieved
factoring exponent `1/2`.

The proof uses actual modular tilting direct sums and their torus weight
spaces. It never substitutes a state-space dimension for a matrix rank.
The representation formulas are the ones documented, with primary
references, in [SL2 coefficient ranks](https://github.com/1337isnot1337/integer-factoring/blob/05d9c1b8aed4f8494ce230f0869f9048ef879474/research/sl2-coefficient-ranks.md).
The new normalization, shell projection, entropy argument, and exact
computations used substantial OpenAI assistant support. No Lean
formalization is included.

A separate research agent independently derived the projected-entry
formula and checked the detecting block, three-chain construction,
extension-field basis argument, and multiplicity estimate. These checks
do not replace external review.

## 1. The shell and its size

For `s=0,...,R`, take `m_s=s+1` binary coordinates of weight `2^s`.
Put

```math
\mathcal S_h=\{x\in\{0,1\}^{M}:\sum_i w_i x_i=h\},\quad
M=(R+1)(R+2)/2,\quad D_h=|\mathcal S_h|.
```

Each coordinate may instead be indexed by `(j,k)`, with weight
`2^(j+k)`. Reading the bits `x_(j,k)` as the binary digits of a
nonnegative occupation number `a_j` gives a bijection

```math
\mathcal S_h\longleftrightarrow
\{(a_0,a_1,\ldots):\sum_j2^j a_j=h\}.
```

Thus `D_h` is the binary partition number, with generating function
`prod_(j>=0)(1-q^(2^j))^(-1)`. In particular,

```math
\log_2 D_h=\tfrac12R^2+O(R\log(R+2)).                 (1)
```

For an elementary upper bound, allow each occupation independently from
zero to `floor(h/2^j)`. The logarithm of their product is at most
`R^2/2+O(R)`. For a lower bound, vary occupations `a_j`,
`1<=j<=R-ceil(log2(R+2))-2`, independently up to
`floor(h/((R+2)2^j))`, and fill the remaining energy with `a_0`.
The sum of the independent contributions is at most `h`, and the
product of their numbers of choices has logarithm
`R^2/2-O(R log(R+2))`.

## 2. The pencil has a large removable scalar factor

On one binary coordinate let

```math
R(x,y,u,v)=[x=y][u=v],\quad
C(x,y,u,v)=[x=u][y=v],\quad
E(x,y,u,v)=[x=v][y=u],\quad H=R+C+E.
```

All additions here are in characteristic two. The weighted pencil is

```math
P_h(z)=\bigotimes_i(R_i+z^{w_i}H_i+z^{2w_i}C_i)
               \big|_{\mathcal S_h^4}.                 (2)
```

It satisfies `P_h(1)=E_(D_h)`. Put
`q(z)=1+z+z^2` and `W=sum_i w_i=R 2^(R+1)+1`.
Every nonzero local pattern has zero, two, or four ones. The zero- and
four-one patterns have local value `q(z)^(w_i)`. In a global entry,
let `b` be the total weight of four-one patterns. Equality of the four
label energies forces the two opposite orientations of each two-one
matching to have the same weight. Call these weights `r,c,e`.
Then

```math
b+r+c+e=h,
\quad
P_h(z)_{xyuv}
=q(z)^{W-2h}\left(q(z)^b z^{c+e}(1+z)^{r+c}\right)^2.   (3)
```

The equal-orientation assertion follows, for example, by subtracting
`energy(x)+energy(y)` from `energy(u)+energy(v)`; this gives twice
the difference between the two `R` orientations as an integer.
The analogous subtractions handle `C,E`.

Consequently there is a primitive pencil of degree at most `2h`,

```math
P_h(z)=q(z)^{W-2h}F_h(z^2),\qquad F_h(1)=E_{D_h}.       (4)
```

Its entries are exactly the expression inside parentheses in (3).
The common factor has value one at `z=1`; removing it is essential
before comparing coefficient ranks. Multiplication by a nonzero scalar
polynomial and substitution `z -> z^2` do not change the linear span of
the coefficient tensors: coefficient convolution is injective on entry
polynomials, hence its transpose is surjective. Thus the original and
primitive pencils have the same coefficient span.

Swapping the two noncrossing matchings gives the degree-`2h` reciprocal
`z^(2h)F_h(1/z)`. All eight crossing-preserving coordinate permutations
either do this or fix the pencil. Its coefficient span is already
closed under these symmetries.

## 3. Exact projected modular blocks

In the row cut `(x,y)|(u,v)`, the XOR
`g=x xor y=u xor v` is conserved. Its local matrices are

```math
M_1(z)=z^w(J+z^w I),\qquad
M_0(z)=(1+z^w)J+z^{2w}I,                              (5)
```

where `J` is the two-by-two all-one matrix. Both determinants equal
`z^(4w)`. For a fixed XOR vector, the row domain is the product of the
weight-`W(g)/2` layer on its support and the weight-`h-W(g)/2` layer
on its complement. This gives an elementary independent exact-rank
algorithm for small shells.

For the faster algorithm, let `T(n)` be the characteristic-two `SL2`
tilting module of highest weight `n`. If
`n+1=2^j+u`, `0<=u<2^j`, then

```math
T(n)=\bigotimes_{i<j}V^{(i)}
       \otimes\bigotimes_{i:u_i=1}V^{(i)},\qquad
\dim T(n)\le(n+1)^2.                                  (6)
```

Let `mu_k(n)` be its multiplicity in `V^(tensor k)`. These are genuine
direct sums, and they are recovered by descending character subtraction.
In one equal-weight class of `m` coordinates, take an on-XOR type
`T(n)` and an off-XOR type `T(p)`. Its multiplicity is

```math
\nu_m(n,p)=\sum_{k=0}^m {m\choose k}\mu_k(n)\mu_{m-k}(p). (7)
```

The corresponding `GL2` determinant twist has exponent
`(m-n-p)/2`. The separate on/off twists need not be retained, because
both matrices in (5) have the same determinant. On all classes, put

```math
\Delta=\tfrac12\sum_s2^s(m_s-n_s-p_s),\quad
A=\tfrac12\sum_s2^s n_s,\quad
B=h-\Delta-A.                                         (8)
```

The block survives shell projection exactly when the on-XOR tensor
of twisted natural modules has a state of energy `A` and the off-XOR
tensor has a state of energy `B`. Its multiplicity is
`prod_s nu_(m_s)(n_s,p_s)`. The projection commutes with the constant
intertwiner because the latter preserves all torus weights.

For two on-layer states let `i` be their weighted intersection. For two
off-layer states let `j` be their weighted intersection. After the
normalization (4), the block entry is

```math
z^{2\Delta+A}(1+z)^{i+B-j}q(z)^j.                     (9)
```

To check (9) directly, the on factor before normalization is
`z^(2A)(1+z)^(2i)`. The off factor is
`(1+z)^(2(B-j))q(z)^(W_off-2B+2j)`, and the determinant
twist contributes `z^(4Delta)`. Since
`W_off-2B=W-2h`, division and the substitution in (4) give (9).

There are `exp(O(R log(R+2)))` choices of block types. Each actual
block has dimension at most
`prod_s(m_s+1)^4=exp(O(R log(R+2)))`. Polynomial construction and
Gaussian elimination therefore compute every coefficient rank in
`exp(O(R log(R+2)))` bit operations, including multiplicities. This
is subexponential in `log D_h`. It is an exact rank algorithm; a Weyl
filtration alone would not justify it.

## 4. A central functional common to both cuts

There are two disjoint singleton labels of energy `h`, using two distinct
weight-`h` coordinates. At the crossing tuple built from them, (3) gives
the primitive entry `z^h`. Consequently evaluation at that tensor entry
defines a linear functional `ell` on the coefficient span with

```math
\ell([z^s]F_h)=[s=h],\qquad \ell(E_{D_h})=1.            (10)
```

Any crossing decomposition has a detected summand. If `z^h` belongs to
the entry-polynomial span of one projected block, every detected summand
acts nontrivially on every copy of that block. Its row rank is therefore
at least the block multiplicity. Reciprocity fixes `z^h`, giving the
same bound for the column rank of the **same summand**. This reasoning
works over every characteristic-two extension field.

The following elementary basis will produce such a block. Put
`u=1+z`, so `q=u^2+u+1`. For any `b>=0`, the `2b+2` polynomials

```math
u^t q^{b-t},\quad u^{t+1}q^{b-t},\qquad 0\le t\le b,  (11)
```

are a basis of polynomials of degree at most `2b+1`. Indeed, in a linear
relation, reduction modulo `q` forces both coefficients at `t=b` to
vanish: `u` is a unit modulo `q`, and a polynomial of degree less than
two cannot be divisible by `q` unless zero. Divide by `q` and repeat.
This argument is valid even when `q` splits over the coefficient field.

## 5. One detected block has almost the full row dimension

Fix `S=16`, put `a=2^S`, and, for all sufficiently large `R`, take

```math
K=R-\lceil\log_2(R+2)\rceil-4.                         (12)
```

For classes outside `S,...,K`, choose on type `T(0)` and off type
`T(m_s)`. At `s=S`, choose on `T(2)` and off `T(7)`. This is an
available type: the `k=10` term of (7) is at least
`binom(17,10) mu_10(2) mu_7(7)=16 binom(17,10)>0`, by (6).
For the remaining classes choose
an even on highest weight `n_s` and an off highest weight `p_s>=7`
with large multiplicity as established below.

All `Delta,A,B` in (8) are multiples of `a`. Moreover,

```math
\Delta+A=\tfrac12\sum_{s=S}^K2^s(m_s-p_s)
 \le\tfrac12\sum_{s=0}^K(s+1)2^s
 \le (K+1)2^K\le h/16.                               (13)
```

For the last inequality, use `2^K<=h/(16(R+2))` and `K+1<=R+2`.
Thus `B=h-Delta-A>=15h/16>=Delta`, and `B<=h-a` because
`A>=a`.

Every even `T(n_s)` has a state of half its total energy: the subset
sums of the powers in (6) contain every integer from zero to `n_s`.
Choose complementary half-energy states independently in all on
classes. They have intersection zero. In the prefix `T(2)`, change
the two choices to the same one of its two weight-`a` coordinates;
all other classes remain complementary. The intersection is now `a`.
The on intersection set therefore contains `0,a`.

For every energy `2^t`, `S<=t<R`, the off module has at least three
coordinates of that weight. The classes `t,t-1,t-2` each contribute
one: every chosen off type there has highest weight at least seven,
so (6) contains the three weights `1,2,4`. This includes the fixed
low prefix since `m_(S-2)=S-1>=7`.

Reserve three disjoint binary chains with one coordinate of each weight
`a,2a,...,h/2`. Each chain represents every multiple of `a` from
zero to `h-a`. For any multiple `t` of `a` in `[0,B]`, use one chain
for a common subset of weight `B-t` and the other two for disjoint
subsets of weight `t`. This gives two off-layer states of energy `B`
whose intersection is `B-t`.

Combining these pairs with the on intersections `0,a`, (9) contains
all polynomials

```math
z^{2\Delta+A}(1+z)^{t+\epsilon a}q(z)^{B-t},
\quad t=0,a,\ldots,B,\quad\epsilon\in\{0,1\}.          (14)
```

Use (11) with `z^a` in place of `z`. Because `Delta<=B`, the span in
(14) contains
`z^(2Delta+A) z^(B-Delta)=z^h`. This is the detecting block required
in Section 4.

It remains to make its multiplicity large. For every sufficiently large
`m`, choose an even `k` within two of `m/2`, with `m-k-7>=0`.
Since there are at most `k+1` summands and every summand has dimension
at most `(k+1)^2`, some even `n` satisfies

```math
\mu_k(n)\ge 2^k/(k+1)^3.
```

Similarly, some `r` in `V^(tensor(m-k-7))` has multiplicity at least
`2^(m-k-7)/(m+1)^3`. The tilting module
`T(r) tensor V^(tensor7)` contains `T(r+7)` with multiplicity one,
because its highest weight is `r+7` with coefficient one. Thus set
`p=r+7>=7`. Formula (7) gives

```math
\nu_m(n,p)\ge {m\choose k}\frac{2^{m-7}}{(m+1)^6}
             \ge \frac{4^m}{128(m+1)^9}.              (15)
```

The final coarse bound uses
`binom(m,k)>=2^m/(m+1)^3` for a central `k` within two positions.
This bound follows from the maximum binomial coefficient and at most
two adjacent coefficient ratios; a larger fixed polynomial denominator
would give the same result.

The finitely many prefix classes affect only a constant. Multiplying
(15) for `S<s<=K`, the detecting block has multiplicity `N_R` with

```math
\log_2 N_R
 \ge 2\sum_{s=S+1}^K(s+1)-O(R\log(R+2))
 =R^2-O(R\log(R+2))
 =2\log_2D_h-O(R\log(R+2)).                           (16)
```

The last step uses both elementary bounds in (1). Every decomposition
therefore has one summand with both ranks at least `N_R`. This proves
the stated `rho=1` asymptotic obstruction.

## 6. Exact finite computations

The checker computes all primitive coefficient ranks in genuine
direct-sum blocks. An independent construction in fixed-XOR blocks
agrees for every coefficient through `R=5`. The central-functional
bound sums the multiplicities of **all** blocks whose entry-polynomial
span contains `z^h`.

| `R` | `D_h` | central-functional rank floor | `D_h^sqrt(2)` | raw coefficient maximum `rho` |
|---:|---:|---:|---:|---:|
| 3 | 10 | 37 | 25.955 | 0.831379 |
| 4 | 36 | 269 | 158.834 | 0.828451 |
| 5 | 202 | 3,185 | 1,820.787 | 0.832541 |
| 6 | 1,828 | 75,559 | 41,032.783 | 0.841114 |
| 7 | 27,338 | 4,449,975 | 1,881,636.356 | 0.850701 |
| 8 | 692,004 | 796,377,453 | 181,618,796.994 | not computed |

The finite floors exclude `sum sqrt(r_i c_i)<D_h^sqrt(2)` throughout
these full coefficient/symmetry spans, over every characteristic-two
field. The separate [weighted cut-budget theorem](https://github.com/1337isnot1337/integer-factoring/blob/27339f61d41af6167dd5ebe947260907b97f604c/research/weighted-cut-rank-budget.md)
already excludes this finite criterion for all tensors through 46
labels; the entries with `D_h=202` and larger are beyond that cutoff.
The decimal comparisons in the table are explanatory; the replay uses
exact integer inequalities for its strict certificates.

The asymptotic theorem is proved in Section 5. It is not extrapolated
from this table. A broader [balanced-layer obstruction](binary-energy-symmetry-barrier.md)
for all equal-weight-position-invariant pieces gives the different
threshold `rho>=3/4`, uniformly over all fields; the stronger `rho=1`
here concerns only the characteristic-two reciprocal coefficient span.

## 7. Replay

```sh
python3 tools/check_weighted_energy.py --output evidence/weighted-energy-reciprocal-pencil.json
python3 tools/weighted_energy_pencil.py --R 5 --primitive --output evidence/energy-R5-direct.json
python3 tools/weighted_energy_pencil.py --R 5 --tilting --primitive --output evidence/energy-R5-tilting.json
python3 tools/weighted_energy_pencil.py --R 8 --central-bound --output evidence/energy-R8-central.json
```

The current exploratory script uses Python's arbitrary-precision
integers and no external packages. Full rank calculations are more
expensive than the common-functional certificate: the latter takes
about four seconds at `R=7` and about two minutes at `R=8` in the
current environment. These are measurements of the checker, not
complexity claims about a factoring implementation.

The compact replay passes 129 binary-partition counts, 65 polynomial
basis instances, 921 independently computed projected entries, all
coefficient-rank comparisons through `R=5`, central-functional bounds
through `R=7`, and explicit high-energy type constructions through
`R=512`. The optional `--full` replay also recomputes the `R=8`
finite central-functional bound. Its all-size claim rests on the proof,
not on these finite checks.
