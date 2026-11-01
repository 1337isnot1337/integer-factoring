# Exact joint channels for the multiplication network's Fourier fan-out

## Result and scope

The input fan-out in the pinned multiplication network admits an exact,
polynomial-time optimal codec when its recipients form a Fourier coset.
This is a joint compression of several original wires. It uses the actual
bounded integer digits and the composite transform modulus; it does not
replace that modulus by a field or count each wire separately.

The result also identifies a limitation of this compression. Split the
recipient layer into `P` equal Fourier cosets, where `P` is a power of two.
For the source parameters and `P>=16`, **all but four cosets have the full
possible message range**. Their equality matrices have full column rank
over every coefficient field. The four exceptional bundles have explicit
smaller ranks. As `P` grows, the sum of the logarithms of all optimal bundle
ranks approaches the uncompressed wire volume, even though the complete
fan-out taken as one bundle has a much smaller rank.

This gives a constructive local rank theorem and an obstruction to a
specified bundle-by-bundle compression. It does not bound the rank of a
whole contracted multiplication instance, prove a factoring lower bound,
or improve the factoring score beyond `1/2`. Joint changes spanning several
recipient groups and algorithms exploiting the later multiplication and
carry constraints remain outside the result.

## 1. The source map and its composite-ring interpolation

The source is Rudich and Rousseau's multiplication network at
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
The relevant definitions are `ph`, `faV`, `digZ`, `mS`, `mV`, and `mLink` in
[`FactoringAux/MulNet.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean),
and `faOut` and `pL` in
[`FactoringAux/Transform.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean).

Put

```math
h=2^s,\quad A=2h,\quad L=2h^2,\quad B=2^L,\quad
M=2^{4L}+1,\quad R=\mathbb Z/M\mathbb Z,\quad
\zeta=16^A=2^{4L/h}\pmod M.
\tag{1}
```

An input gate has a vector `u=(u_0,...,u_(h-1))` with `0<=u_i<B`.
Ignoring a known invertible coordinate multiplier, its outgoing messages
are

```math
F(u)_j=\sum_{i=0}^{h-1}u_i\zeta^{ij},\qquad 0\le j<A.
\tag{2}
```

The actual message from gate `faQ r`, or `faP r` before imposing its prefix
weight, is `16^(rj) F(u)_j`. Each multiplier is a unit modulo the odd
integer `M`. Removing or restoring these multipliers is a coordinatewise
bijection, so it changes neither range sizes nor equality-matrix ranks.
The `faQ` settings and weights are unrestricted by every divisor prefix.

### Unit-difference lemma

Every difference between two distinct elements of
`1,zeta,...,zeta^(A-1)` is a unit in `R`.

Indeed, `zeta^h=-1`. If `0<d<A`, write `d=2^v d_0` with `d_0` odd and
let `r=log2(A)-v`. Then `eta=zeta^d` satisfies
`eta^(2^(r-1))=-1`. Consequently

```math
(1-\eta)\sum_{j=0}^{2^{r-1}-1}\eta^j=2.
\tag{3}
```

The element two is a unit, so (3) supplies an explicit inverse of
`1-eta`. Multiplication by a power of `zeta` proves the assertion.
This argument never assumes that `M` is prime.

Thus every square Vandermonde matrix on distinct chosen evaluation
points is invertible over `R`: its determinant is a product of units.
Lagrange interpolation computes its inverse using polynomially many ring
operations on `O(L)`-bit residues.

### Arbitrary recipient sets

For a set `J` of `k` recipients, let `E_J` be the zero-one matrix whose
rows are all `B^h` input vectors, whose columns are all `R^J` message
vectors, and whose entry is one precisely when the messages agree with
(2). Over **every** coefficient field,

```math
\mathrm{rank}(E_J)=|F_J(\{0,\ldots,B-1\}^h)|.
\tag{4}
```

Each row is a standard coordinate vector. Selecting one row per image
element gives an identity minor, and duplicate rows give the matching
upper bound. The rank field is independent of the composite ring used to
define the messages.

For every `k`, interpolation gives

```math
B^{\min(k,h)}\le\mathrm{rank}(E_J)
   \le\min(B^h,M^k).
\tag{5}
```

For `k<=h`, vary only the first `k` input digits and set the others to
zero; the corresponding square Vandermonde map is injective. For `k>=h`,
any `h` selected messages recover all input digits, so the exact rank is
`B^h`. In particular, the full `A`-wire fan-out has rank `B^h`, and a
channel of exactly `hL` bits suffices. Its factors can encode the input
digits and check their evaluations directly in polynomial bit time.

## 2. Exact ranks on Fourier cosets

Let `P>=2` be a power of two dividing `A`. Put

```math
k=A/P,\qquad t=P/2,\qquad
J_c=\{c+Pj:0\le j<k\},\quad 0\le c<P.
\tag{6}
```

For `c=0`, define

```math
I_0=t(B-1)+1.
\tag{7}
```

For `c>0`, write `q=2^(v_2(c))` and put

```math
r=t/q,\qquad H=q(B-1)+1,\qquad C=2^{4L/r},\qquad
I_c=
\begin{cases}
H^r,&H\le C,\\
M,&H>C.
\end{cases}
\tag{8}
```

Here `r` divides `h`, so the exponent defining `C` is an integer.
The exact bundle rank is

```math
\boxed{\mathrm{rank}(E_{J_c})=I_c^k.}
\tag{9}
```

The same statement holds for any positive digit width `L` with `h|4L`,
using `zeta=2^(4L/h)`, provided the zero-coset interval in (7) lies below
`M`. The source parameters satisfy that condition. The checker uses small
members of this more general family as independent exhaustive controls.

### Proof of (9)

Write every input position uniquely as `ell+ka`, with `0<=ell<k` and
`0<=a<t`. Let `gamma=zeta^c` and `xi=zeta^P`. Removing the known source
coordinate multipliers, the `k` messages are

```math
F(u)_{c+Pj}
=\sum_{\ell=0}^{k-1}\xi^{\ell j}\gamma^\ell
       \sum_{a=0}^{t-1}u_{\ell+ka}\beta^a,
\qquad \beta=\gamma^k.
\tag{10}
```

The outer `k`-point Fourier transform is invertible over `R` by the
unit-difference lemma. The factors `gamma^ell` are units. Different
values of `ell` involve disjoint input digits. Therefore the bundle range
size is the `k`th power of the range size of

```math
\sum_{a=0}^{t-1}v_a\beta^a,\qquad 0\le v_a<B.
\tag{11}
```

If `c=0`, then `beta=1`. Every integer between zero and `t(B-1)` occurs
as the sum of `t` bounded digits. The interval does not wrap modulo `M`,
which proves (7).

If `c>0`, write `c=q d` with `d` odd. Then `beta=C^d` and `C^r=-1`.
The powers `beta^a` distribute among the `r` signed basis positions
`1,C,...,C^(r-1)`. Each position occurs exactly `q` times: multiplication
by `d` permutes the residues modulo `r`, and shifting `a` by `r`
reverses the sign. Thus (11), after a fixed translation, has exactly the
same range as

```math
\sum_{j=0}^{r-1}w_jC^j\pmod M,
\qquad 0\le w_j<H=q(B-1)+1.
\tag{12}
```

The digit sums in different positions are independent and fill their
entire integer intervals, including when some signs are negative. If
`H<=C`, positional representation is unique and the largest value is
at most `C^r-1=M-2`; hence the range size is `H^r`.

If `H>C`, consecutive digit intervals overlap. By induction on the
number of positions, the attainable integers fill the entire interval

```math
0\le z\le(H-1)(1+C+\cdots+C^{r-1}).
```

Its upper endpoint is at least `C^r=M-1`, so every residue modulo `M`
occurs. The range size is `M`. Equations (4) and (10) now prove (9).

## 3. An explicit optimal channel, including its decoder

The rank bound in (9) is attained by zero-one factors computable in
polynomial bit time. The channel is the Cartesian product of `k` sets
of size `I_c`, encoded by one mixed-radix integer below `I_c^k`.
Padding to `ceil(log2(I_c^k))` binary bits adds only zero channel values.

For each `ell`, apply the fixed translation used in (12). If `H<=C`,
encode its `r` ordinary base-`C` digits as a base-`H` integer. If `H>C`,
encode the translated residue itself, as an integer below `M`. For
`c=0`, encode the ordinary sum below `I_0`. This defines a unique channel
value for every source digit vector.

For a proposed output vector, first remove the factors `16^(r_gate*j)`.
Apply the inverse outer Fourier transform in (10), then multiply by
`gamma^(-ell)` and apply the same fixed translation. When `H<=C`, reject
unless its canonical residue has at most `r` base-`C` digits, each less
than `H`. Otherwise record those digits in base `H`. When `H>C`, every
residue is permitted. The zero-coset test is the interval bound in (7).

Let `L(u,theta)` be the indicator that this procedure encodes the source
input as `theta`, and let `R(theta,w)` be the indicator that it decodes
the proposed messages as the same `theta`. Then

```math
E_{J_c}(u,w)=\sum_{0\le\theta<I_c^k}
                   L(u,\theta)R(\theta,w).
\tag{13}
```

All normalization, translations, Fourier operations, integer digit
operations, and range tests take `poly(A,L)` bit operations. No image
enumeration or elimination in an exponentially large matrix is needed.
The source instance has `A,L` polynomially bounded by its input length,
so these are effective local factors in the required bit model.

For the unsplit full fan-out, use the original `h` digits as the channel.
Interpolation from any `h` messages recovers a proposed digit vector;
checking its digit bounds and all `A` evaluations supplies the column
factor. This also realizes its exact rank in polynomial bit time.

## 4. Fine coset partitions retain almost all original wire width

Now impose the source value `L=2h^2`. For `P>=16`, formulas (7)--(8)
simplify to

| Recipient coset | Number of cosets | One folded coordinate's range size |
| --- | ---: | --- |
| `c=0,P/2` | 2 | `(P/2)(B-1)+1` |
| `c=P/4,3P/4` | 2 | `((P/4)(B-1)+1)^2` |
| All other `c` | `P-4` | `M` |

For the last row, the root order is at least eight. At order eight,
`r=4`, `C=B`, and `q=P/8>=2`, so `H>C`. At larger orders, `r>=8` and
`C<=sqrt(B)<B<=H`. For order four, `r=2` and
`H< B^2=C`; for orders one and two the interval in (7) lies below `M`.
These inequalities hold because `P<=2h` and `B=2^(2h^2)`.

Write `W_P` for the sum of the base-two logarithms of the exact ranks of
all `P` bundles of one input gate. Then

```math
\begin{aligned}
W_P=k\bigl(& (P-4)\log_2M
 +2\log_2((P/2)(B-1)+1)\\
 &+4\log_2((P/4)(B-1)+1)\bigr).
\end{aligned}
\tag{14}
```

In particular, since every displayed nontrivial interval contains at
least `B` integers and `log2 M>4L`,

```math
W_P\ge (8-20/P)hL.                         \tag{15}
```

Conversely, elementary upper bounds in (14) give

```math
W_P=(8-20/P)hL+O\!\left(\frac{h\log P}{P}+h\right).
\tag{16}
```

The harmless `O(h)` term can be replaced by the much smaller explicit
contribution from `log2(1+2^(-4L))`. Thus for any source sequence with
`P` tending to infinity,

```math
W_P/(hL)\longrightarrow8.
\tag{17}
```

The original source writes `A(4L+1)=8hL+2h` bits for this fan-out.
Taking the full fan-out as one bundle needs only `hL` bits. Taking
independent Fourier-coset bundles loses almost all that saving in the
fine-partition limit. Every separate exact factorization of these
bundle equality matrices must respect (14), since its channel count is
at least the corresponding rank over the coefficient field.

This is a statement about the chosen independent bundle interface.
It does not add ranks across arbitrarily correlated factors of a whole
group tensor. In particular, it does not justify a lower bound after
the later multiplication, inverse-transform, or carry constraints have
been contracted. Those constraints are where a different arithmetic
method would have to obtain additional compression.

## 5. The later local transform gates have full output ranges

The source middle and inverse gates do not have a comparable hidden
reduction of their unrestricted output alphabets. Let `omega=16`, so
`zeta=omega^A`, and let `e_0` denote the residue vector with first
coordinate one and all other coordinates zero. Expanding the source
`midOut` and using Fourier orthogonality gives the exact identity

```math
\mathrm{midOut}(c,a,e_0)_j
   =A\omega^{-cj}a_j,\qquad 0\le j<A.
\tag{18}
```

Indeed, the forward transform of `e_0` is the all-ones vector; the
remaining forward transform followed by the unnormalized inverse
transform is `A` times the identity. Every multiplier in (18) is a unit
in `R`. Thus the joint output map of one middle gate is onto `R^A`,
and every selected `k`-output bundle has equality-matrix rank `M^k`
over every coefficient field. The corresponding optimal channel is
simply the `k` selected residues.

Likewise the source `ibOut` is `A^(-2)` times an invertible `A`-point
inverse Fourier transform. Its entire output range is `R^A`, and every
selected `k`-output bundle again has exact rank `M^k`. These statements
use the original gate setting domains from `mS`; they do not assume
that all simultaneous settings extend to a complete factoring witness.

There is a stronger forward-network statement for the middle gate.
For `h>=8`, fix a middle position

```math
c\notin\{0,A/4,A/2,3A/4\}.
```

Set one complete multiplication input to the integer one. Its messages
to this middle gate are exactly `e_0`. Section 4 with `P=A` proves that
each individual input gate of the other operand maps **onto `R`** at
recipient `c`. Different input gates use disjoint digit vectors, so
their messages can be chosen independently. The optimal codec above
constructs digits attaining any requested residue.

Given any desired middle output vector `w in R^A`, choose the incoming
vector by `a_j=A^(-1) omega^(cj) w_j`, and realize its coordinates through
those independent input gates. Equation (18) produces exactly `w`.
Therefore all `M^A` middle output vectors occur even after imposing
the full preceding forward-transform consistency and the original
bounded-digit input domains, with the other operand fixed to one.
The construction takes polynomial bit time in `A,L`.

This statement concerns the forward subnetwork before the target
product and carry tests. Once a fixed product `N` is imposed, the
operand-one branch has only the other operand `N`; it supplies no
lower bound for the final contracted instance. The distinction prevents
using a local full-range theorem as a claimed factoring lower bound.

## 6. Checks and remaining question

Run

```sh
python3 tools/check_factoring_network_structure.py \
  --output evidence/factoring-network-structure.json
```

The checker enumerates the original digit maps in small source and
generalized instances, compares all coset image cardinalities with
(7)--(9), verifies composite-ring unit differences and inverse Fourier
maps, checks the optimal mixed-radix encoder/decoder, and constructs
small equality matrices over two coefficient fields. The actual source
case `h=2,L=8` uses the composite modulus
`2^32+1=641*6700417`. Separate radix controls include the boundary
`H=C` and the first full-range case `H=C+1`.
The checker also verifies the original middle and inverse formulas,
constructs consistent bounded-digit forward inputs attaining arbitrary
middle outputs, and checks the all-but-four formula on larger source
parameters using exact integers. Finite controls support the
implementations; the all-size proof is the argument above.

The remaining constructive question is whether the **contracted**
multiplication and carry constraints create joint boundary relations
that can be preserved at each recursive level. The source's individual
fan-out correlations have now been computed exactly for the stated
recipient partitions. Their mere existence does not yield a smaller
logarithmic exponent.

The source transform and multiplication network are due to Rudich and
Rousseau. This ordinary rank analysis and checker were developed with
OpenAI assistant support. No Lean files are changed.
