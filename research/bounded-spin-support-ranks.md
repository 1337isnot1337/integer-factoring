# Uniform ordinary-rank factors from bounded physical word support

## A rank bound below the post-half scale for specified pieces

Let `n=2r>0` and use the height-three spin-zero encoders of the
[coherent construction](bounded-spin-crossing.md). Every encoder vanishes
outside the binary words whose prefix imbalance stays between -3 and 3
and ends at zero. Complementary words have proportional encoder columns.

For the [six-vertex family](bounded-spin-six-vertex.md), the three extreme
coefficient tensors are equal, and both ordinary cut ranks of their
shared tensor satisfy

```math
\boxed{r(T_{r,0,0}),c(T_{r,0,0})\le a_r,
\qquad
a_r=\frac{(2+\sqrt2)^r+(2-\sqrt2)^r}{4}.}             \tag{1}
```

The channels and their factors have explicit polynomial-time unranking
and evaluation algorithms. For a whole path alphabet of size
`L_n=F_(n-1)`, this is `O(L_n^(4/3))`, strictly below the desired
`L_n^sqrt(2)` rank-cost scale. The three extreme tensors give only one
direction and do not constitute a crossing decomposition.

The same construction bounds every near-corner coefficient. If
`A+B+C=r` and `k=B+C`, then

```math
r(T_{A,B,C}),c(T_{A,B,C})
 \le a_r {n\choose2k}{2k\choose2C}.                  \tag{2}
```

Analogous bounds hold at the other corners. A consequence is that bands
of width `k=o(n/log n)` have cheap individual factors but their combined
span cannot contain the crossing for all large `n`. Thus a construction
using these factors must add different, typically asymmetric-rank
directions. No stronger factoring bound is asserted here.

## 1. Support and its exact growth

For a word `w`, let `W_i` be the number of ones in its first `i` bits.
Along a path of current height `h_i`, the magnetic index in the encoder
is

```math
s_i=W_i-\frac{i-h_i}{2}.
```

Nonzero transitions require `0<=s_i<=h_i<=3`. Consequently
`|i-2W_i|<=3` at every prefix, and return to spin zero requires
`W_n=n/2`. This proves the support assertion over the integers, before
any reduction modulo five. Further modular cancellations may shrink
support and therefore do not invalidate the bound.

Call this word set `S_n`. Bitwise complement preserves `S_n`. Since each
spin-zero vector is fixed by `epsilon^tensor n`,

```math
v_P(\overline w)=(-1)^r v_P(w).                      \tag{3}
```

For positive length there are no complement-fixed words, so choose the
representative with first bit zero in each pair.

The imbalance walk lives at heights -3 through 3. Grouping opposite
heights, two steps act on the even absolute heights zero and two by

```math
M=\begin{pmatrix}2&1\\2&2\end{pmatrix}.
```

Thus `|S_(2r)|=(M^r)_(0,0)`. The eigenvalues are `2+sqrt(2)` and
`2-sqrt(2)`, and the two diagonal entries of `M^r` agree. Dividing the
support count by two proves (1)'s channel count. Equivalently

```math
a_1=1,\quad a_2=3,\quad
a_{r+2}=4a_{r+1}-2a_r\quad(r>=1).
```

The first values are `1,3,10,34,116,396`.
Suffix counts for the seven imbalance states unrank each representative
in polynomial bit time. The first bit is fixed to zero and the remaining
walk begins at height one. All counts have `O(n)` bits.

## 2. Extreme tensor factors

The coefficient `T_(r,0,0)` uses only `0000` and `1111` vertices, so all
four binary words agree. With the prescribed output dual norms,

```math
T_{r,0,0}(P,Q,R,S)
 =2\sum_{w\in S_n/\text{complement}}
       v_P(w)v_Q(w)\frac{v_R(w)v_S(w)}{\nu_R\nu_S}.   \tag{4}
```

Complementary terms are equal by (3), proving the factor two. Grouping
`P,Q` against `R,S` gives a rank factorization with `a_r` channels.
Grouping `P,R` against `Q,S` gives the other factorization with the same
channels, distributing the two dual norms between the factors.

For the other two extreme monomials, two of the four words are
complemented. The two factors from (3) cancel, so

```math
T_{r,0,0}=T_{0,r,0}=T_{0,0,r}.                       \tag{5}
```

This proves the stated rank bounds and the single-direction limitation.
The [factor implementation](../tools/bounded_spin_support.py) uses
zero padding for channels beyond the representative count.

To compare growth rates, put `phi=(1+sqrt(5))/2` and
`g=sqrt(2+sqrt(2))`. Then `a_r=O(g^n)` and `L_n` is within constant
factors of `phi^n`. The elementary inequalities

```math
g^6=20+14\sqrt2<41<89/2<\phi^8
```

give `g<phi^(4/3)`, and `4/3<sqrt(2)`. More precisely, the exponent
`log_phi(g)` is approximately 1.27589. This is a rate for the specified
pieces, not a factoring score.

## 3. Near-corner factors

For `T_(r-k,B,C)`, let `J` be the `2k` positions using a non-`a/f`
vertex, and let `D` be the `2C` crossing positions within `J`. Given
the first word `x`, the four words are necessarily

```math
x,\qquad x\mathbin{\mathrm{xor}}J,\qquad
x\mathbin{\mathrm{xor}}D,\qquad
x\mathbin{\mathrm{xor}}J\mathbin{\mathrm{xor}}D.       \tag{6}
```

Conversely these words contribute exactly the specified monomial when
the encoder entries are nonzero: the four weight equations impose equal
counts in each vertex pair. Words outside `S_n` automatically give zero.
Choose `x` up to complement, choose `J`, and choose `D` inside it.
Equation (6) gives a four-leg rank-one decomposition with precisely the
channel bound (2), including a factor two for complement pairing.
It therefore bounds all three flattenings, including the required two.
The other corners follow by exchanging the two output legs, which swaps
`V,W`, and by the epsilon/leg transformation exchanging `U,W` in the
[cut reciprocity proof](bounded-spin-projector.md). These transformations
preserve a four-leg rank-one decomposition.

Subset unranking and mixed-radix indexing give explicit channel words in
polynomial bit time; `corner_words` implements them. Each factor is the
product of two encoder entries, with the appropriate dual norms. Its
channel count is at most `2^O(n)`. Dense factorization can cap widths at
the physical matrix dimension where necessary, as in the
[modular interface](bounded-spin-modular-interface.md).

For a corner band with `k=o(n/log n)`, the binomial factors in (2) are
`exp(o(n))`. There are only `O(n^2)` monomials in all three bands.
Every tensor in their entire span consequently has row rank at most
`exp(o(n))*g^n`. The crossing has row rank `L_n^2`, or a constant-factor
equivalent power-of-two rank after binary selection. Since `g<phi^2`,
the crossing is outside that span for all sufficiently large `n`.
This remains true with arbitrary coefficient recombinations or field
extensions. Adding other coefficient directions is outside this exclusion.

## Reproduction

Run `make verify-bounded-spin-support` (NumPy required by the checker).
It exhausts the physical words through length twelve, verifies encoder
support and the integral complement relation, and checks complete channel
unranking. It verifies both full extreme rank factorizations at lengths
two through eight and six complete first-corner factorizations against
independent coordinate blocks. The recurrence is checked through length
forty. The [receipt](../evidence/bounded-spin-support.json) pins the
sources and the finite scope. These proofs and computations used
substantial Codex AI assistance. Uniform favorable ranks for a complete
crossing decomposition remain open.
