# A multiplication constructor with directly decoded boundary tables

## Result

External first-transform values determine almost all the bounded input
digits inside a group. A cyclotomic argument proves that the bounded
input setting is unique. After it is recovered, every middle, inverse,
and carry setting is forced. All table entries are therefore Boolean,
and an explicit grouping gives fewer intergroup wires.

**Theorem.** Let `N>=2` have binary length `n`. Let `P>=2` be an integer
with `P^40<=n`. For every `1<=ell<=n` and every divisor-prefix query
`0<=s<=n`, `0<=x<2^(n-s)`, a deterministic algorithm writes an
Appendix B instance `I_(ell,s,x)` satisfying

```math
\begin{aligned}
\operatorname{zval}(P,I_{\ell,s,x})
 &=\#\{d\mid N:\operatorname{bitlength}(d)=\ell,\
                    \lfloor d/2^s\rfloor=x\},\\
\operatorname{WFe}(P,\lambda,B_{\rm tab},I_{\ell,s,x}),\qquad
P^2\lambda&<17n,\qquad B_{\rm tab}=1.             \tag{1}
\end{aligned}
```

The full written instance, including every binary table entry and all
padding, is constructed in

```math
n^{O(1)}P\,2^{(P-1)\lambda+F}
\le n^{O(1)}2^{18n/P},
\qquad P^2(F+1)<\frac{101}{100}n.                 \tag{2}
```

All polynomial degrees and constants are absolute. Uniformly as `P`
tends to infinity under the cutoff,

```math
P^2\lambda=(12+o(1))n,\qquad
P^2F=(1+o(1))n.                                  \tag{3}
```

The per-entry work has exponent `F=O(n/P^2)`. The complete row still has
`2^((P-1)lambda)` labels; (2) includes this written-output cost.
Thus this is an improved top constructor, not a claim of a
`2^O(n/P^2)` full-table algorithm or a new factoring exponent.

For the original adaptive constructor with two half-length supports,
the same decoder computes each Boolean entry in polynomial bit time.
The smaller length-conditioned width uses the enumeration cost in (2).

## 1. Arithmetic parameters and exact factor assignments

Use the [operand-length construction](operand-length-transform-constructor.md),
with the digit width rounded to an even integer. Put

```math
\begin{gathered}
h=\text{least power of two at least }P^9,\quad
A=2h,\quad D=2h^2,\quad T=A^2=2D,\\
L=T-A-1,\qquad
w=2\left\lceil\frac{n+1}{2L}\right\rceil,\quad b=2^w,\quad
t=\log_2D,\\
m=D\left\lceil\frac{2w+t}{D}\right\rceil,\quad
M=2^m+1,\quad f=m+1,\quad c=w+t,\quad K=D(b-1).
                                                               \tag{4}
\end{gathered}
```

For the length class `ell`, let

```math
d_p=\lceil\ell/w\rceil,\qquad
d_q=\lceil(n-\ell+1)/w\rceil .
```

The first factor has exactly `ell` bits; its leading digit has its
specified leading bit equal to one. The second factor is bounded above
by `2^(n-ell+1)`. Its shortened highest digit may have leading zeros.
All positions above these supports are fixed zero. The support bound

```math
d_p+d_q\le A^2-A,\qquad \min(d_p,d_q)\le D          \tag{5}
```

is unchanged by even rounding.

At residue index `r`, the two first-transform gates store respectively
`h_p(r)` and `h_q(r)` digits, in positions `r+A i`. Their coefficient
indices are initial segments, and

```math
h_p(r)+h_q(r)\le A.                              \tag{6}
```

The first transform sends, for example,

```math
a_r(j)=\omega^{rj}
       \sum_{i=0}^{h_p(r)-1}p_{r+A i}\eta^{ij},
\qquad \omega=2^{m/D}\pmod M,\quad \eta=\omega^A. \tag{7}
```

The other first transform has the same formula. Middle and inverse
gates use the exact three-stage formulas from the operand-length note.
Since `omega^D=-1`, those formulas are valid over the composite ring
`R=Z/MZ`. Each convolution coefficient is below

```math
D(b-1)^2<2^{2w+t}\le2^m<M.
```

Canonical inverse-transform residues are therefore the exact integer
coefficients. The carry recurrence, its bound `[0,K]`, and its clipped
outgoing-message convention are unchanged.

The existing assignment proof consequently applies: fixing the two
input integers forces every transform and carry setting, and acceptance
is equivalent to their product being `N`. Prefix checks and the exact
first-factor length give the count in (1). No new arithmetic
approximation or residue lift is introduced here.

The source transform is pinned to
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
The [adaptive constructor](adaptive-transform-constructor.md) proves
the principal-root identities, carry domains, and finite-alphabet
interface used in this note.

## 2. Boundary observations give an exact Vandermonde decoder

Consider any group whose first-transform, middle, and inverse families
contain the **same** index set `I subset {0,...,A-1}`. Write `k=|I|`.
There are `A-k` middle recipients outside the group.

Fix one complete boundary label. It gives a canonical residue for each
original equality edge crossing the boundary. A residue label outside
`[0,M-1]`, a carry label outside `[0,K]`, or nonzero padding makes the
entry zero immediately.

An internal first-transform gate with `a` stored digits is a polynomial
of degree below `a`, evaluated at the external recipient points in (7).
Let

```math
u=\min(a,A-k).
```

Choose the first `u` external recipient indices in a fixed order.
Enumerate the `a-u` highest bounded digits, subtract their contributions
from these observations, and solve for the first `u` coefficients.
Highest-digit boxes and prefix bits may restrict the enumeration further;
using the full digit box is an upper bound.

### Why the interpolation is valid over the composite modulus

The `A` evaluation points are powers of `eta`, where
`eta^(A/2)=-1`. For distinct indices, write their nonzero difference
modulo `A` as `2^v d`, with `d` odd. Set `q=A/2^(v+1)`. Then

```math
(1-\eta^{2^v d})\sum_{j=0}^{q-1}\eta^{2^v d j}=2.
                                                               \tag{8}
```

The modulus is odd, so two is a unit. Every difference of two distinct
evaluation points is consequently a unit. The square Vandermonde
determinant is a product of such differences, hence a unit. Lagrange
interpolation supplies an explicit inverse; no field assumption is
needed.

The factors `omega^(rj)` in (7) are also units and are removed before
interpolation. The recovered ring elements are accepted as input digits
only if their canonical representatives satisfy their actual digit
domains. Every unused external observation is checked afterwards.

Thus each free-digit choice gives at most one accepted input setting.
Conversely, every input setting matching the boundary occurs for its
own enumerated high digits and is recovered exactly by the invertible
Vandermonde system.

For the two input gates at a fixed `r in I`, the total number of free
digits is

```math
(h_p(r)-(A-k))_+ +(h_q(r)-(A-k))_+\le k.           \tag{9}
```

This follows from (6): if only one term is positive, it is at most `k`;
if both are positive, their sum is at most `2k-A<=k`.
Summing over the `k` internal residue indices gives at most `k^2`
free digits.

## 3. Every remaining internal setting is forced

Once the input gates have been decoded, an internal middle gate gets
every stored coordinate either from a decoded internal first-transform
gate or from its boundary label. Its setting is unique. Its outgoing
values are then computed and checked on all external edges.

The same reasoning fixes every internal inverse gate: each input
coordinate comes from an internal middle gate or a boundary label.
Every carry coefficient is then fixed by an internal inverse gate or
by its boundary edge.

The internal carry gates form disjoint runs along the original integer
carry chain. A run starts with incoming carry zero if it begins at
position zero. Otherwise its predecessor is external, so the incoming
carry is supplied by the boundary. The recurrence fixes every carry in
the run. Check every digit of `N`, the carry range, each external outgoing
carry, and the final zero carry where applicable.

Arbitrary residue settings can exceed the genuine convolution bound.
The existing rejection of an outgoing carry above `K` remains part of
this decoder. Clipping a rejected message to zero never creates an
accepted setting.

Every internal equality, message, digit domain, leading-bit condition,
prefix check, and local weight is checked after reconstruction. The
decoder therefore computes exactly the original grouped-table entry.
It is not only a bound on its support.

Writing `kappa=ceil(A/P)` for the largest group size used below, every
entry is a sum of at most

```math
b^{\kappa^2}=2^F,\qquad F=\kappa^2w                 \tag{10}
```

zero-or-one contributions. The next argument proves that at most one
contribution can be accepted, so the actual entry is zero or one.
All modular interpolation, transform arithmetic, checks, and additions
have bit cost polynomial in `n` per enumerated choice. The integers
and matrices have polynomially bounded dimensions and bit lengths.

### 3.1. Bounded input recovery is unique

We prove the following stronger statement for the explicit grouping:
if `I` is a consecutive interval of at most `A/2` recipient indices,
then the external observations uniquely determine every bounded input
polynomial with at most `A-1` coefficients. Equation (5) guarantees this
coefficient bound for every operand-length input gate.

Suppose two such input polynomials have identical external values.
Their difference `f` has integer coefficients of absolute value below
`b`, degree at most `A-2`, and satisfies

```math
f(\eta^j)=0\pmod M\qquad(j\notin I).                \tag{BI1}
```

For successive powers of two `L=A,A/2,...,4`, we show inductively that

```math
f(X)=S_{A/L}(X^L)g(X),\qquad
S_r(Z)=1+Z+\cdots+Z^{r-1},
\quad\deg g\le L-2,\quad |g_i|<b.                \tag{BI2}
```

The initial case is `g=f`. The primitive `L`th transform points have
recipient indices `(A/L)j`, with `j` odd. They form an arithmetic
progression of step `2A/L`. An interval of at most `A/2` indices contains
at most `L/4` of them. There are therefore at least `L/4` external
primitive points.

At a primitive `L`th point, the prefactor in (BI2) equals `A/L`.
It is a power of two and hence a unit modulo `M`. Thus `g` vanishes
at each of those external points.

Put `h=L/2` and reduce `g` modulo the cyclotomic polynomial
`Phi_L(X)=X^h+1`. Its integer remainder `v` has degree below `h` and

```math
\|v\|_2<2b\sqrt h.                                \tag{BI3}
```

This follows directly from `v_i=g_i-g_(i+h)` and the strict coefficient
bounds. We now justify the required lower bound if `v` were nonzero.

For any `k` selected primitive points, evaluation gives a surjective
homomorphism of additive groups

```math
\mathbb Z[X]/(X^h+1)\longrightarrow(\mathbb Z/M\mathbb Z)^k.
```

Surjectivity follows from the unit Vandermonde determinant on the first
`k` monomials, using (8). Its kernel consequently has index `M^k`.
If `v` vanishes at the selected points, the ideal generated by `v`
lies in that kernel.

The polynomial `X^h+1` is irreducible over the rationals: after replacing
`X` by `X+1`, Eisenstein's criterion at two applies. Multiplication by
a nonzero `v` therefore has a nonzero integer determinant on the free
module with basis `1,X,...,X^(h-1)`. The ideal index equals the absolute
value of that determinant and is divisible by `M^k`. Its columns are
the negacyclic shifts of the coefficient vector of `v`, all of the
same Euclidean norm. Hadamard's determinant bound gives

```math
\|v\|_2^h\ge M^k,\qquad
\|v\|_2\ge M^{k/h}.                               \tag{BI4}
```

Take the at least `k=L/4=h/2` external primitive points. Equations
(BI3)--(BI4) contradict

```math
2b\sqrt{L/2}\le2b\sqrt{A/2}
 \le\frac{bA}{\sqrt2}<\sqrt M,                    \tag{BI5}
```

where `A>=4` and `M>Db^2=b^2A^2/2` were proved from (4).
Hence the remainder `v` is zero, and `g` is divisible by `X^(L/2)+1`
over the integers.

The quotient has degree at most `L/2-2`. Its coefficients are the
lower-half coefficients of `g`, so their absolute values remain below
`b`. Multiplying the two geometric factors gives (BI2) with `L/2`
in place of `L`. This proves the induction without increasing the
coefficient bound.

After the last division,

```math
f(X)=c\,S_{A/2}(X^2),\qquad |c|<b.
```

At least one of recipients `0,A/2` lies outside a consecutive interval
of length at most `A/2`. At that point the geometric factor equals the
unit `A/2`. Thus `c=0 mod M`; since `|c|<b<M`, it is the integer zero.
The two bounded input polynomials are identical.

Each block in Section 4 has size at most `A/2`, and each input gate
has at most `A-1` digits. Its input setting is therefore unique.
Section 3 then forces every other setting. This proves
`B_tab=1` for all length-conditioned tables.

Uniqueness is used only to bound the entry value. The certified
algorithm still enumerates the free digits in (10) and checks every
candidate. No unproved fast lattice-recovery algorithm is assumed.

## 4. An explicit grouping with narrow pair boundaries

Partition `{0,...,A-1}` into `P` consecutive blocks `I_a`, each of size
`floor(A/P)` or `ceil(A/P)`. Every block has at least two indices.
Put the gates at index `r` in all four transform families into the group
of `r`. This is the common-index property required by the decoder.

Write carry positions as `j=j_1+A j_2`. Normally assign the carry to the
group containing `j_1`. For the first index of **every** block `I_a`,
including index zero, instead assign its carry to group

```math
j_2\bmod P.                                      \tag{11}
```

These first positions are the relays. The last position of every block
keeps its natural group because each block has at least two indices.

For any distinct pair of groups, the edge counts are bounded as follows:

| Edges | Count | Bits per edge |
|---|---:|---:|
| The two first-to-middle families and middle-to-inverse family | `6 kappa^2` | `f` |
| Inverse-to-carry edges changed by the relays | `2 kappa` | `f` |
| Carry-chain edges | `4 kappa` | `c` |

The first row counts both directions across the group pair in each of
three complete bipartite families.

For the second row, only a block's first index has a recolored carry.
There are at most `kappa` rows of either relay color, giving at most
`kappa` edges in each direction across the group pair.

Every noninternal carry-chain edge is incident to a relay. It is either
the relay-to-next-position edge or the preceding-core-to-relay edge.
At a row boundary, the last position has natural color `P-1`, and
connects to the next row's relay at index zero. This is included in the
same second incidence. Each relay color occurs in at most `kappa` rows.
For an unordered pair there are at most four relevant incidences,
each with at most `kappa` edges. Missing initial and final chain edges
only reduce the count.

Choose the common positive integer channel width

```math
\lambda=6\kappa^2 f+2\kappa f+4\kappa c,\qquad
F=\kappa^2w,\qquad B_{\rm tab}=1.                 \tag{12}
```

The exact off-diagonal codec concatenates the actual original edge
labels and pads to this width. Each original equality edge uses one
shared value; it is not counted as two independent boundary messages.
The local row rejects nonzero padding.

No coloring search is needed. The relay rule is deterministic and
preserves the common transform index set even though carry gates have
been redistributed.

## 5. Uniform constants

Put `C=Dw`. Even rounding increases the old width by at most one. The
elementary parameter bounds give

```math
\begin{gathered}
A\ge1024,\quad D\ge2^{19},\quad
D^2/n\le1/4\ (P=2),\quad D^2/n<19/50\ (P\ge3),\\
t/D\le19/2^{19},\qquad C<\frac{50101}{100000}n,\\
f\le2w+t+D,\qquad
\epsilon=P/A\le\frac1{2P^8}\le\frac1{512},\qquad
P \kappa/A\le1+\epsilon.                           \tag{13}
\end{gathered}
```

The bounds on `D^2/n` include all small group counts. At `P=2`,
`D=2^19`, so `D^2/n<=2^38/2^40=1/4`. At `P=3`,
`D=2^31`, and the exact integer inequality
`50*2^62<19*3^40` gives the stated `19/50` bound.
For `P>=4`, `D<8P^18` gives
`D^2/n<64/P^4<=1/4`.

For the bound on `C`, use
`w<(n+1)/L+2`, `D/L<1001/2000`, `D/n<=2^-20`,
and `n>=2^40`. In particular
`C<(1001/2000)(n+1)+2D<(50101/100000)n`.

Since `A^2=2D`,

```math
\begin{aligned}
A^2 f&\le4C+2Dt+2D^2,\\
A^2 c&=2C+2Dt.
\end{aligned}
```

Substitution into these inequalities yields

```math
\frac{A^2f}{n}<
\begin{cases}
2.50406,&P=2,\\
2.76407,&P\ge3,
\end{cases}
\qquad
\frac{A^2c}{n}<1.00205.                           \tag{14}
```

Use the larger first constant uniformly. Equations (12)--(14) imply

```math
\begin{aligned}
\frac{P^2\lambda}{n}
&<6(1+\epsilon)^2(2.76407)
  +2\epsilon(1+\epsilon)(2.76407)
  +4\epsilon(1+\epsilon)(1.00205)\\
&\le\frac{21846986847}{1310720000}<17.
\end{aligned}                                                    \tag{15}
```

Also `P^2/n<=2^-38` and `A^2w=2C`, so

```math
\frac{P^2(F+1)}n
<1.00202\left(\frac{513}{512}\right)^2+2^{-38}
<\frac{101}{100}.                                \tag{16}
```

As `P` tends to infinity under the cutoff, `C/n` tends to `1/2`,
`f/w` tends to two, and `P \kappa/A` tends to one. The lower-order
terms in (12) vanish after multiplying by `P^2/n`. This proves (3).

## 6. Full written cost and divisor recovery

For one group there are exactly `2^((P-1)lambda)` off-diagonal clique
labels. Enumerate them all. For each label use the decoder from
Sections 2--3 and write its zero-or-one result. Equation (10) bounds
the number of enumerated settings; Section 3.1 bounds the output by one.
Thus the full construction cost is

```math
n^{O(1)}P\,2^{(P-1)\lambda+F}.                    \tag{17}
```

It includes reading the label, rejecting unused codes, modular inverses,
all reconstruction arithmetic, summing accepted settings, writing every
Boolean entry, and writing the channel-width parameter. There is no
materialization of a smaller table followed by an uncharged expansion.
Equations (15)--(16) give (2).

The old grouped-table definition is unchanged: each entry counts
internal settings matching the specified external wires. Replacing
its enumeration by an exact decoder therefore preserves the complete
network count and the Appendix B finite-alphabet predicate.

Sum the `n` disjoint length-class counts to answer an ordinary divisor
prefix query. To recover the smallest divisor at least two, subtract
the known contribution of divisor one whenever it lies in the current
prefix interval. Descend the binary prefix tree, taking the left half
when this adjusted count is positive. There are at most `n` prefix
queries and at most `n` length-conditioned instances per query.
The result equals `N` when `N` is prime. These polynomially many calls
do not alter the exponential constructor bounds.

The leading exponent of (17) is `(12+o(1))n/P`. The per-entry
enumeration term is lower order in this full written bound. Improving
an evaluator's recursion still requires its own proof; entry
observability alone does not change that recursion rate.

## 7. Polynomial bit time per entry for half supports

Apply the same grouping and decoder to the
[first adaptive constructor](adaptive-transform-constructor.md), where
every first-transform gate has exactly `A/2` digits in each factor.
For every `P>=2`, each group has `k<=A/2`, so

```math
A-k\ge A/2.
```

All input digits are recovered by square Vandermonde interpolation
from their external observations. There are no enumerated free digits.
The remaining settings are forced as in Section 3. Each Boolean table
entry is therefore computed in polynomial bit time.

This corollary uses the half-support width choice, not the smaller
length-conditioned width in (4). Its pair width is still given by
(12), with its own `f,c`; its leading top constant is 24. The
length-conditioned theorem retains Boolean entries with the smaller
leading top constant 12, while using the enumeration bound (16).

## 8. Exact implementation controls

The companion checker reconstructs small transform networks, boundary
labels, Vandermonde solutions, and carry runs. It compares decoded
entries with literal internal-setting or input-setting counts, including
invalid codes, shortening of the highest digit, and row-boundary carry
edges. It also checks the relay edge counts and the explicit parameter
bounds with integer arithmetic.

The [exact checker](../tools/check_observable_boundary_constructor.py)
and [evidence record](../evidence/observable-boundary-constructor.json)
are reproduced by

```sh
python3 tools/check_observable_boundary_constructor.py \
  --output evidence/observable-boundary-constructor.json
```

The completed controls include 10,832 bounded input fibers; 2,000 complete
contiguous-boundary maps with 72,220 distinct images; 73 relay colorings
with 3,998,407 original edges counted; 48 valid and 96 mutated table-entry
comparisons; and 6,028 exact large-parameter cases. All 48 valid local
entries have exactly one accepted setting. The controls also exhibit
collisions when the input height reaches the excluded value `A`, when
the hidden recipient set is not an interval, and when the modulus is
too small. Thus they check the substantive hypotheses of Section 3.1.

These finite controls support the implementation. The unit-difference
identity (8), free-digit bound (9), cyclotomic injectivity proof, forced
propagation, and pair-edge count prove the uniform construction.
