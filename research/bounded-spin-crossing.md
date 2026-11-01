# Coherent bounded-spin crossings over the integers and F5

## Construction and scope

There is an explicit coherent crossing over `F_5` whose physical labels
are bounded Clebsch--Gordan paths. It is defined uniformly at every
length, has a strict endpoint, and admits a finite-state coefficient
constructor. The same family has an explicit integral companion with
known crossing scale and linear coefficient bit growth.
Return-to-zero paths have Fibonacci growth; a prefix gives
exactly `2^lambda` labels with underlying length `O(lambda)`.

This changes the physical encoding by linear combinations of binary
words. It is therefore outside the coordinate-code hypothesis of the
[weighted block obstruction](difference-support-rank-floor.md). That
fact is an opening, not a rank upper bound. **No uniform post-half
contraction or improved factoring exponent is established.**

The main advance here is a complete effective family over one fixed
characteristic, with no growing denominator problem. No Lean proof or
independent mathematical review is claimed.

## 1. Integral path vectors with nonzero modular norms

The vectors use the standard `SL_2` Clebsch--Gordan steps; the elementary
argument below fixes their integral normalization and modular scope.
Let `V` have orthonormal basis `e_0,e_1`. On binary tensor words, `E`
changes one `1` to `0` and sums the results; `F` changes one `0` to `1`.
On words of length `n` and weight `k`, their commutator is

```math
[E,F]=(n-2k)I.
```

Take paths starting at height zero, with steps `+1,-1`, staying in
`{0,1,2,3}`. Associate an integral highest-weight vector to each path.
Start with `v_empty=1`. If the current path has height `d`, set

```math
v_{P+}=v_P\otimes e_0,\qquad
v_{P-}=d v_P\otimes e_1-Fv_P\otimes e_0.             \tag{1}
```

Only admissible next heights are used. The new vectors are annihilated
by the total raising operator. For a downward step this follows from
`EFv=dv`, since `Ev=0`.

Paths of a common length are mutually orthogonal over the integers.
Their squared norms satisfy

```math
\nu_{P+}=\nu_P,\qquad
\nu_{P-}=d(d+1)\nu_P.                               \tag{2}
```

Here is the inductive orthogonality argument. Different final heights
give different binary weights. Within an upward branch the old inner
product is unchanged. Within a downward branch it is multiplied by
`d^2+d`, using `F=E^T` and `EFv=dv`. At a common final height, an
upward branch and a downward branch have cross inner product
proportional to `<w,Fv>=<Ew,v>=0`.

At a downward step `d` is one, two, or three, so the factors in (2)
are `2,6,12`. None vanishes modulo five, at **any** path length.
Consequently the reduced vectors are linearly independent over `F_5`,
and the paired dual row is explicitly

```math
\widehat v_P=\nu_P^{-1}v_P.                         \tag{3}
```

This restriction avoids the degenerating norms that arise if arbitrarily
high spin paths are reduced modulo a fixed prime.

## 2. A product formula for every encoder entry

The vector need not be expanded to evaluate one coefficient. Write
`v_(P,r)=F^(r)v_P` using divided powers. The Clebsch--Gordan recursion is

```math
v_{P+,r}=v_{P,r}\otimes e_0+v_{P,r-1}\otimes e_1,
```

and

```math
v_{P-,r}=(d-r)v_{P,r}\otimes e_1
                 -(r+1)v_{P,r+1}\otimes e_0.        \tag{4}
```

Divided powers here are integral sums over subsets of positions; no
division in the finite field is required. The displayed formula also
follows by applying the total divided lowering operator to (1).

For a fixed binary word, track the current magnetic index `r`, starting
at zero. If its next bit is `b`, the transition rules are:

| Path step from height `d` | New magnetic index | Coefficient multiplier |
| --- | --- | --- |
| Up | `r+b` | `1` |
| Down, `b=1` | `r` | `d-r` |
| Down, `b=0` | `r-1` | `-r` |

Discard transitions outside the interval from zero to the new height.
Require final magnetic index zero. The product of the multipliers is
exactly `v_P(word)`. There are at most four magnetic states and the
entry takes `O(n)` arithmetic operations. Formula (2) computes the dual
normalization in the same time. Over the integers, coefficient and norm
bit lengths are `O(n)` because the allowed heights are bounded.

## 3. Compress the four legs coherently

At a binary position use the ordinary integer pencil

```math
K(z)=(1-z)I+zE,\qquad
I(x,y,u,v)=[x=u][y=v],\quad E(x,y,u,v)=[x=v][y=u].
```

Let `A,B` be any selected path sets of the same underlying length `n`.
Use `A` on the first/fourth legs and `B` on the second/third legs. Define

```math
G_n(z;P,Q,R,S)=
\sum_{x,y,u,v}v_P(x)v_Q(y)\widehat v_R(u)\widehat v_S(v)
             \prod_{i=1}^n K(z;x_i,y_i,u_i,v_i),     \tag{5}
```

over `F_5`. At one, the product in (5) is the full word crossing.
Orthogonality and (3) therefore give the strict identity

```math
\boxed{G_n(1;P,Q,R,S)=[P=S][Q=R].}                  \tag{6}
```

In particular its at most `n+1` coefficient tensors sum to the crossing,
with scale one. Different path sets on the two crossing pairs are allowed.

This is intrinsically coherent. Nonzero path vectors of positive binary
weight are annihilated by `E`, whereas no single coordinate word of that
weight is annihilated by `E`: its raised words are distinct. Thus these
path spaces do not have a coordinate-word basis. The literal coordinate
blocks used in the all-code elementary obstruction are not available
after the compression in (5).

## 4. Uniform finite-state construction and exact binary widths

Contract the four product formulas in section 2 alongside the local
pencil. A state is the tuple of four magnetic indices, at most `4^4=256`
possibilities. At each position sum over the 16 local binary four-tuples.
Polynomial degree increases by at most one. This computes every complete
entry polynomial in `O(n^2)` field operations, with an absolute constant.
Both dual norms are explicit nonzero field elements.

The two-step path-count matrix on even heights zero and two is

```math
M=\begin{pmatrix}1&1\\1&2\end{pmatrix}.
```

Hence the number of return-to-zero paths of length `2r`, `r>=1`, is
`F_(2r-1)`: the first values are `1,2,5,13,34`. These are integer
counts of distinct orthogonal vectors, not dimensions reduced modulo five.

For width `lambda>=1`, choose the least even `n>=4` with at least
`2^lambda` return paths and keep the first `2^lambda` in up-first order.
There are at least `2^floor(n/4)` return paths, by concatenating the two
four-step returns `UDUD` and `UUDD`, with a final `UD` if needed. Thus
the chosen `n` is `O(lambda)`.

Suffix walk counts give unranking in polynomial bit time: at each step,
count the admissible completions of the upward branch, compare the
binary label with that count, and select or skip the branch. Counts have
`O(n)` bits. The implementation provides `binary_path(width,label)`.

The complete crossing entry polynomial therefore has degree `O(lambda)`
and is computable in polynomial time in `lambda`. Generic dense rank
factorization over the fixed field is possible in `2^O(lambda)` time,
but a favorable uniform rank bound and factor algorithms compatible
with the evaluator remain unproved.

### Integral companion and coefficient bit bounds

Each norm in (2) divides `12^floor(n/2)`, since every downward factor
divides 12 and there are at most `floor(n/2)` downward steps. Thus

```math
H_n(z)=12^n G_n(z)
```

has integer coefficients: after computing the unnormalized four-frame
contraction, multiply it by the integer
`12^n/(nu_R nu_S)`. The same orthogonality proof over `Q` gives the
integer polynomial identity

```math
\boxed{H_n(1)=12^n E_{A,B}.}                        \tag{7}
```

The integer constructor is `integral_polynomial_entry`. It uses the same
finite-state contraction and exact integer arithmetic. A path coefficient
has absolute value at most `3^n`. There are at most `16^n` word tuples,
each local-product coefficient has magnitude at most `2^n`, and the
clearing multiplier is at most `12^n`. Therefore every coefficient of
`H_n` has magnitude at most `31104^n`, giving `O(n)=O(lambda)` bits.
The crossing scale has the same linear bit bound. The integer entry
constructor consequently takes polynomial bit time in `lambda`.

Reducing (7) modulo five gives a nonzero scalar crossing and the same
coefficient ranks as the normalized `F_5` construction. This supplies
an integer identity, **not** a lift of favorable modular rank bounds
to ranks over `Q`. A factoring improvement would still need the required
uniform rank theorem over the relevant ring or a separate valid
fixed-characteristic counting transfer, with complete costs.

## 5. Relation to fusion categories and the open rank problem

The bounded-height graph is the fusion graph of the simple object `L_2`
in the Verlinde category at prime five; the fusion rules
are stated in Etingof, Ostrik and Venkatesh,
[*Computations in symmetric fusion categories in characteristic p*](https://arxiv.org/abs/1512.02309),
section 2.1. Ostrik's
[*On symmetric fusion categories in positive characteristic*](https://arxiv.org/abs/1503.01492)
also discusses the characteristic-five Fibonacci example.

This is mathematical motivation for the encoding. Identifying all
operators in (5) with a categorical quotient, or converting categorical
dimensions into ordinary matrix ranks, would require additional proof.
Neither identification is assumed in (1)--(6).

The concrete next question is whether these explicit finite-state
coefficients admit effective factors satisfying

```math
\sum_j\sqrt{r(A_{n,j})c(A_{n,j})}
 \le\operatorname{poly}(n)|A|^{\sqrt2-\epsilon}
```

for equal-size selected path sets and a fixed `epsilon>0`, possibly after
specified coefficient recombinations. A recoupling or transfer argument
must bound the two ordinary cut ranks of the same pieces. Path growth,
category theory, and a low rank at one finite length do not establish it.

An [exact primary-algebra calculation](bounded-spin-rank-algebra.md) now
settles the full eight-site coefficient span over every extension of
`F_5`: its minimum maximum rank product is 2304, exceeding `13^3`.
The double-root sector has genuine nilpotent rank thirteen, but all
combinations below the threshold lie in a hyperplane missing the crossing.
This closes that finite block, not the uniform family. An all-length
primary decomposition or a changed kernel/selected subspace remains open.

## 6. Independent finite controls

Run `make verify-bounded-spin`. The standard-library checker compares
the magnetic product formula with the independent integral recursion for
131,891 word coefficients, verifies 6,766 integer Gram entries and 232
highest-weight vectors, and checks exact path unranking and binary widths.
Sixteen complete four-bit entry polynomials agree with literal word
sums both over the integers and modulo five, including the integer scale
`12^4=20736` and the reduction relation. Twelve additional six-bit
endpoint cases pass. Omitting dual norm
division is a negative control: the two-bit endpoint is four instead of
one modulo five.

The four-bit coefficient rank pairs are
`(4,1),(4,1),(3,3),(1,4),(1,4)`. They do not give an improved asymptotic
rate. The [receipt](../evidence/bounded-spin-crossing.json) pins the exact
source and finite scope. The construction, proof and implementation used
substantial Codex AI assistance; the strongest factoring exponent is `1/2`.
