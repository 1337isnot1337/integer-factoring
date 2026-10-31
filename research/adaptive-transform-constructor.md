# An adaptive multiplication network with smaller written top instances

## Result

The multiplication network can choose its digit count from the requested
clique size, then choose its digit width from the input length. This
removes the source's square-capacity padding and reduces the transform
modulus to the size needed for exact convolution. A separate carry width
reduces the remaining chain messages.

Theorem. Let `N>=2` have binary length `n`. For every integer `P>=2` with
`P^40<=n`, and every divisor-prefix query `0<=k<=n`,
`0<=x<2^(n-k)`, a deterministic algorithm writes an Appendix B instance
`I` with positive integer channel width `lambda` and entry bound `B_tab`
such that

```math
\begin{aligned}
\mathrm{WFe}(P,\lambda,B_{\rm tab},I),\qquad
\mathrm{zval}(P,I)
 &=\#\{d:d\mid N,\ \lfloor d/2^k\rfloor=x\},\\
P^2\lambda&\le99n,\\
B_{\rm tab}&\le117n/P.
\end{aligned}                                             \tag{1}
```

The complete construction, including grouping search and writing every
table entry in binary, takes `2^O(n/P)` bit operations. The constants and
the polynomial factor's degree are absolute.

Uniformly as `P` tends to infinity subject to `P^40<=n`, the stronger
estimates are

```math
P^2\lambda\le(36+o(1))n,
\qquad P B_{\rm tab}\le(40+o(1))n.                 \tag{2}
```

These improve the explicit top-instance bounds `4717n` and `6769n/P`
in [half-top-constructor.md](half-top-constructor.md). The new constants
include both the changed arithmetic network and a sharper count of its
actual setting domains. They are not obtained by replacing one number
in the old inequalities.

The ordinary factoring exponent remains `alpha=1/2`: the existing
evaluator still has the same level schedule and logarithmic contraction
rate. The result is a positive constructor improvement, not a new score
or a new Lean theorem.

## 1. Parameters without square-capacity padding

Let `h` be the least power of two with `h>=P^9`. Define

```math
\begin{gathered}
A=2h,\qquad D=2h^2,\qquad T=A^2=2D,\qquad G=3P^8,\\
w=\lceil n/D\rceil,\qquad t=\log_2D,\qquad b=2^w,\\
m=D\left\lceil\frac{2w+t}{D}\right\rceil,
\quad M=2^m+1,\quad f=m+1,\quad c=w+t,\\
K=D(b-1),\qquad C=Dw.
\end{gathered}                                             \tag{3}
```

Here `D` is the number of digits in each factor, `w` is the number of
bits in each digit, `C` is the resulting input capacity, `f` is a residue
message width, and `c` is a carry message width. The letter `t` denotes
only the integer `log2 D` in this note.

The parameters are explicit and require no prime search or integer
factorization. Their elementary bounds are

```math
\begin{gathered}
P^9\le h<2P^9,\qquad 2P^{18}\le D<8P^{18},\\
n\le C<n+D,\qquad 2w+t\le m\le2w+t+D-1,\\
f\le2w+t+D,\qquad c\le f.
\end{gathered}                                             \tag{4}
```

In particular `C` exceeds `n` by fewer than `D` bits. In the source
parameterization, rounding to the next permitted square capacity can
instead increase the capacity by a factor approaching sixteen.

Two useful uniform bounds under the theorem's cutoff are

```math
h\ge512,\qquad D\ge2^{19},\qquad D^2\le n/2.       \tag{5}
```

For `P>=4`, the last bound follows from
`D^2<64P^36<=n/4`. For `P=2`, the definitions give
`D^2=2^38<=n/4`. For `P=3`, they give `D^2=2^62`, and
`2^63<3^40<=n`. These cover every integer `P>=2`.

## 2. A smaller ring still gives exact convolution

Work in the composite ring `R=Z/MZ`. Set

```math
\omega=2^{m/D}\pmod M,
\qquad\omega'=\omega^{T-1},
\qquad T^{-1}=2^{(2m-1)(t+1)}\pmod M.
\tag{6}
```

Then `omega^D=-1`, `omega^T=1`, and the displayed value of `T^-1`
is correct because `2^(2m)=1` and `T=2^(t+1)`. All the relevant
power-of-two Fourier root sums have their usual values. For example,
the factorization of a power-of-two geometric sum gives zero whenever
the frequency is nonzero modulo `T`. Equivalently, differences between
distinct transform points are units: the same geometric identity
expresses a suitable multiple of each difference as two, which is a
unit modulo the odd number `M`.

The source's Cooley--Tukey identity is proved for a general commutative
ring with these root identities. Use its `faOut`, `midOut`, and `ibOut`
formulas with the new `A`, `omega`, `omega'`, and `T^-1`. The two input
sequences are supported in positions `0,...,D-1`, and `T=2D`, so their
ordinary convolution has no cyclic aliasing.

For digits `0<=p_i,q_i<b`, each coefficient satisfies

```math
0\le a_j=\sum_{i+i'=j}p_iq_{i'}
\le D(b-1)^2
<D b^2=2^{2w+t}\le2^m<M.                       \tag{7}
```

Therefore the canonical residue returned by the inverse transform is
the exact nonnegative integer coefficient `a_j`. This is the only
integer-recovery condition needed; primality of `M` is irrelevant.
The ceiling in (3) enforces both (7) and the explicit root condition
`D|m`, including when `2w+t` is itself divisible by `D`.

The pinned source is
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1),
specifically `sum_root`, `ib_mid_fa`, and the general transform formulas
in [Transform.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean).
The source's later definitions `pL`, `pm`, `omega16`, and its specialized
integer coefficient bound are replaced by (3), (6), and (7).

## 3. Carry bounds and the exact factor-assignment bijection

Retain the original gate topology: `A` gates in each of the `faP`,
`faQ`, `mid`, and `ib` families, and `T=A^2` carry gates. Use these
setting domains and messages:

| Gate | Setting domain | Outgoing messages |
| --- | --- | --- |
| `faP r`, `faQ r` | `h` independent digits in `[0,b-1]` | First transform stage, one residue to every `mid` gate |
| `mid k1` | Two independent vectors in `R^A` | Combined second forward and first inverse stages, one residue to every `ib` gate |
| `ib j1` | One vector in `R^A` | Final inverse stage, one coefficient residue to each of its `A` carry gates |
| `carry j` | A coefficient residue `a` and an integer `r` in `[0,K]` | The next carry described below |

All transform and coefficient edges have width `f`; the consecutive
carry-chain edges have width `c`. Indeed

```math
K+1=D(2^w-1)+1\le D2^w=2^c.                    \tag{8}
```

At a carry gate put `q=floor((a+r)/b)`. Its weight is one exactly when
all of the following hold: its digit `(a+r) mod b` matches the
corresponding base-`b` digit of `N`; `q<=K`; the first carry has `r=0`;
and the last carry has `q=0`. Conditions concerning the first or last
position apply only there. Its forward message is `q` when `q<=K`, and
zero otherwise. Its backward message is its stored `r`.

The clipping convention makes every message fit its declared width even
on settings with weight zero. It adds no accepted settings because the
weight separately rejects `q>K`. This detail matters: arbitrary residue
settings may exceed the exact coefficient bound (7).

For genuine convolution coefficients, start with `r_0=0` and define
`r_(j+1)=floor((a_j+r_j)/b)`. Inductively,

```math
0\le r_j\le K
\quad\Longrightarrow\quad
a_j+r_j\le D(b-1)^2+D(b-1)=bK,
\quad0\le r_{j+1}\le K.                         \tag{9}
```

Thus every actual multiplication has carries in the declared domain.
The initial carry and recurrence determine them uniquely.

Now fix all input digits, hence integers `p,q` in `[0,b^D-1]`.
Every `fa`--`mid` equality fixes one stored coordinate at a middle
gate. All middle settings are consequently forced. Their messages
force every inverse-stage setting, which in turn forces every carry
coefficient. Equation (7) identifies those coefficients with the exact
integer convolution. Carry-chain equalities and the initial carry then
force the unique carries above.

The digit and endpoint checks accept exactly when `pq=N`. One direct
proof multiplies

```math
a_j+r_j=n_j+b r_{j+1}
```

by `b^j` and sums over `0<=j<T`. The carry terms cancel because
`r_0=r_T=0`, leaving `sum_j a_j b^j=N`. Conversely, if `pq=N`, its
ordinary base-`b` carry recurrence satisfies every check, including the
last zero carry because `pq<b^T`.

The `faP` weights additionally require `floor(p/2^k)=x`. Each gate checks
the relevant bits of its own digits, with the prefix extended by zeros
through capacity `C`. This is exactly the desired condition because
`C>=n`, `N<2^n`, and the specified prefix lies below `2^(n-k)`.

Since `N>=2`, accepted `p,q` are positive. Every positive divisor of
`N` and its complementary factor lie below `2^n<=b^D`. We have therefore
proved a bijection between accepted global settings and the divisors
with the requested prefix. All intermediate transform settings and
carries are unique. There is no multiplicity from a residue lift,
an inverse root choice, or a carry variable.

## 4. Total wire volume and local setting volume

There are `4A^2` residue edges: two input-to-middle families, one
middle-to-inverse family, and the inverse-to-carry family. There are
`A^2-1` carry-chain edges. The exact unordered wire volume is

```math
E=4A^2f+(A^2-1)c.
\tag{10}
```

Use the convenient upper bound on ordered wire volume and the bound on
the sum of logarithmic setting counts

```math
U=A^2(8f+2c),\qquad S=A^2(w+4f+c).              \tag{11}
```

The first satisfies `2E<=U`. For the second, the two input families
contribute exactly `2Dw=A^2w`; the middle and inverse families contribute
at most `3A^2f`; and carry settings contribute at most `A^2(f+c)`.
No gate setting count is replaced by an unnecessary common worst case.

Using (4),

```math
\begin{aligned}
U&\le36C+20Dt+16D^2,\\
S&\le20C+10Dt+8D^2.
\end{aligned}                                             \tag{12}
```

For powers of two `D>=2^19`, `t/D<=19/2^19`. Together with
`D^2<=n/2`, `C<n+D`, and `D/n<=1/(2D)<=2^-20`, this gives

```math
\begin{aligned}
U/n&<44+36/2^{20}+190/2^{19}<45,\\
S/n&<24+20/2^{20}+95/2^{19}<25.
\end{aligned}                                             \tag{13}
```

In particular, the new network has fewer than `22.5n` unordered wire
bits throughout the theorem's cutoff. For comparison, the original
network with digit count `D_old` has

```math
(10D_{\rm old}-1)(4D_{\rm old}+1)
=40D_{\rm old}^2+6D_{\rm old}-1\ge40n
```

unordered wire bits. The comparison charges the changed input capacity,
the carry chain, and every transform edge.

As `P` tends to infinity, `D^2/n<64/P^4`, `D/n<8/P^22`, and
`Dt/n<=D^2/n` all tend to zero. Equations (4) and (12), with their
immediate matching leading lower bounds, imply

```math
U=(36+o(1))n,\qquad S=(20+o(1))n,
\qquad E=(18+o(1))n.                              \tag{14}
```

## 5. Grouping and the written table bounds

Use the same consecutive-family clustering map and generic grouping
theorem as [half-top-constructor.md](half-top-constructor.md), Section 1.
The clustering map sends an index `r` in a family of size `a` to
`floor(rG/a)`. Each cluster contains at most
`u_A=floor(A/G)+2` members of each `A`-element family and at most
`u_T=floor(A^2/G)+2` carry gates.

For one representative of each non-carry family, the weighted degrees
sum to `7Af`: the degrees are `Af`, `Af`, `3Af`, and `2Af`.
A carry gate has weighted degree at most `f+2c`. Consequently the
cluster weighted degree is bounded by

```math
Q=7u_AAf+u_T(f+2c).                              \tag{15}
```

The corresponding logarithmic setting bound is

```math
W=u_AA(w+3f)+u_T(f+c).                           \tag{16}
```

Here the two input families together use `Aw` bits per family position,
the middle family uses at most `2Af`, and the inverse family uses at most
`Af`. Carry settings use at most `f+c` bits.

The generic `good_exists` theorem in
[EngineGroup.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/EngineGroup.lean)
supplies a coloring of the `G=3P^8`
clusters into `P` groups, each with at most `floor(2G/P)` clusters,
and with every pair width at most

```math
L_0=\left\lfloor\frac{PU+GQ}{P^3}\right\rfloor.
\tag{17}
```

Find such a coloring by exhaustive search through its `P^G` codes.
Every condition is tested with integer arithmetic on widths and group
sizes. Put

```math
\lambda=\max(1,L_0),\qquad
H=W\lfloor2G/P\rfloor,\qquad B_{\rm tab}=H+1.
\tag{18}
```

The estimates needed for the constants follow directly from (15)--(16):

```math
\begin{aligned}
GQ&\le U+G R,&R&=14Af+2f+4c,\\
GW&\le S+G R_W,&R_W&=2A(w+3f)+2(f+c).
\end{aligned}                                             \tag{19}
```

Use `c<=f`, `w<=f/2`, `A>=1024`, and `G/A<=3/(2P)`. Since
`U>=8A^2f` and `S>=4A^2f`,

```math
\frac{GR}{U}
\le\frac{3(14+6/1024)}{16P}<\frac{2627}{1000P},
\qquad
\frac{GR_W}{S}
\le\frac{3(7+4/1024)}{8P}<\frac{2627}{1000P}.
\tag{20}
```

Together with (13), `P^2<=n`, and `P>=2`, this proves

```math
\begin{aligned}
P^2\lambda
&\le P^2+U(1+1/P+2627/(1000P^2))\\
&<n+45n(1+1/2+2627/4000)<99n,\\
PH&\le2S(1+2627/(1000P))\\
&<50n(1+2627/2000)<116n.
\end{aligned}                                             \tag{21}
```

Since `P<=n`, adding one to `H` gives `B_tab<=117n/P`. Letting
`P` tend to infinity in (18)--(20), using (14) and `P^2/n<=P^-38`,
gives (2).

## 6. Explicit tables, exact interface, and full bit cost

For each final group and each off-diagonal clique label, enumerate its
original local gate settings. Count those with all local weights one,
all internal messages equal, and all external messages matching the
corresponding label blocks. Append zero bits to each pair label up to
width `lambda`, and require those padding bits to be zero.

Each setting contributes zero or one. There are at most `2^H` settings
per group, so every entry is a nonnegative integer below `2^(H+1)`.
Every globally accepted network assignment has a unique shared label
tuple. Conversely, compatible accepted local settings determine one
global assignment. Section 3 therefore identifies the clique value
with the required divisor-prefix count.

The off-diagonal row codec, empty self blocks, treatment of empty groups,
and binary integer encoding are exactly those in
[half-top-constructor.md](half-top-constructor.md), Sections 2 and 4.
They apply to arbitrary finite equality-message networks with declared
per-edge widths, so the different residue and carry widths require no
new adapter argument. The generic cluster grouping theorem is also
unchanged. The arithmetic parameters, local domains, coefficient
recovery, carry checks, and their numerical bounds are the parts replaced
and proved in this note. No specialized source theorem is silently
instantiated outside its original parameter definitions.

The target interface is Appendix B at
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth),
with `PInst`, `zval`, and the entry-bound predicate `WFe` as cited in the
companion constructor note. Thus its existing evaluator can accept the
new written instances directly.

All dimensions, moduli, roots, settings, messages, and counters have
polynomial bit length in `n`. The number of gates is polynomial in `n`.
Each transform message is a finite modular sum of polynomial length;
modular powers, the displayed inverse, and all carry and prefix checks
take polynomial bit time. In particular, an operation on an `m`-bit
residue is not counted as unit cost.

For a fixed absolute constant `a`, a complete conservative cost bound is

```math
n^a\left(P^G+P\,2^{(P-1)\lambda+H}+1\right).
\tag{22}
```

This pays for each local setting separately at every written table
entry, and includes grouping search. Since
`G log2 P<=3P^9<=3n/P` and
`(P-1)lambda+H<=215n/P`, it is `2^O(n/P)`.

With `E_row=2^((P-1)lambda)`, the same explicit target encoding gives

```math
\operatorname{len}(I)
\le6+4\lambda+P\bigl(4+E_{\rm row}(3B_{\rm tab}+10)\bigr).
\tag{23}
```

Equation (22) charges this written output. The choice of dimensions
changes neither the evaluator's level growth nor its contraction rate.
Substituting (1) into the completed half-rate composition preserves its
`2^O(n/sqrt(log n))` runtime with improved construction constants.

## 7. Why a still smaller modulus needs an additional argument

Removing (7) and relying on the same carry checks would introduce false
factor assignments. For example, take `D=8`, `w=6`, but use the smaller
modulus `M=2^8+1=257` and root `omega=2`. The transform identities still
hold. Set `p=q=63`. Its only nonzero convolution coefficient is `3969`,
whose canonical residue is `114`. The carry checks for `N=114` accept
that residue sequence, although `63*63=3969!=114`.

Both factors and the target fit the digit capacity, and the carries fit
the proposed smaller domain. Thus bounded carries alone do not repair
coefficient aliasing. A further reduction would need another exact
recovery condition or a different coupled representation. The theorem
above does not assume such a condition.

## 8. Reproducible checks

Run:

```sh
python3 tools/check_adaptive_transform_constructor.py \
  --output evidence/adaptive-transform-constructor.json
```

The checker compares the actual three-stage network with independent
integer convolution, checks full digit and carry recovery, exhaustively
checks small divisor-prefix counts, and includes the false-assignment
control from Section 7. It also checks exact parameter ceilings and
the written-instance constants on boundary choices of `P,n`. Small
transform models are identified separately from the theorem's large
cutoff. The finite evidence does not replace the uniform bijection or
bit-cost proof.
