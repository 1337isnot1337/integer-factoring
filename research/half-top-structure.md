# A restricted obstruction for direct multiplication-network tables

## Result and scope

The special 0/1 multiplication network does not give a larger logarithmic
saving through separate wire compression, adaptive grouping, or faster
table writing alone, while retaining the construction and explicit pad
used in [half-rate-factoring.md](half-rate-factoring.md).

More precisely, the argument below concerns algorithms with these three
properties:

1. They partition the original multiplication gates into `P` groups and
   visit every original local setting of each group at least once. A
   single pass that scatters settings into table entries is permitted.
2. They factor original crossing-edge tensors separately and concatenate
   their channels between each pair of groups. The resulting common
   channel width is `lambda` bits. The separate factorizations may use
   arbitrary rational coefficients, including negative coefficients.
3. They write the dense Appendix B `PInst` and the explicit pad from the
   completed half-rate composition, with `P=K_ell` and
   `r=2^floor(ell/2)`.

For this class, the complete construction has a cost of at least

```math
2^{\Omega(n/\sqrt{\log n})}.
```

This is a restriction on that construction method, not a lower bound on
factoring, on all factorizations of its tables, or on the actual running
time of an unpadded evaluator. The pad used here is a sufficient pad; the
argument does not show it is necessary for these special inputs. Joint
factorization of a whole group boundary, avoiding Cartesian enumeration
by solving its local constraints, or proving a smaller sufficient pad
for a restricted family of inputs all fall outside the claim.

The result is a scoped obstruction rather than a new factoring upper
bound. The ordinary mathematical upper bound remains `alpha=1/2`.

## 1. An unrestricted layer survives every divisor-prefix query

Use the pinned factoring reference
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
For an input of binary length `n`, the top constructor chooses the least
`s>=0` with `n<=2^(4s+2)`. Write

```math
S=2^s,\qquad h=S,\qquad A=2S,\qquad d=2S^2,
\qquad M=2^{4d}+1.
\tag{1}
```

Here `d` is the source's `pL s`, not a clique label width. The source's
parameters are `ph=h`, `pA=A`, `pM=M`, and `fw=4d+1`. In particular,

```math
n\le4S^4=Ahd.
\tag{2}
```

The definitions `mS`, `mLink`, `mV`, and `mW` in
[MulNet.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean)
give the following concrete layer:

| Gate or edge | Source definition |
| --- | --- |
| `faQ a`, for `a<A` | A vector of `h` independently chosen `d`-bit digits |
| `mid b`, for `b<A` | Two independent vectors of `A` residues modulo `M` |
| `faQ a` to `mid b` | Present for every pair `(a,b)` |
| `faQ` and `mid` weights | Identically one, for every `N,k,x` |
| Message from `faQ a` | `faV s a u b` |
| Return message from `mid b` | Coordinate `a` of its second residue vector |

There are therefore exactly `A^2` edges in this complete bipartite layer,
and every `faQ` gate has exactly `2^(hd)` settings. The divisor-prefix
condition affects the `faP` weights. It does not restrict the typed
settings or the weights of either endpoint in the layer above.

The original edge tensor is the integer equality indicator of its two
messages. This is `Fac.numTest` and `eval_numTest` in
[Gates.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Gates.lean).
Its message-range condition is automatic here: both messages are
residues below `M`, and `M<2^(4d+1)`.

## 2. Each separate wire needs at least d channel bits

Fix an edge from `faQ a` to `mid b`. For each `0<=v<2^d`, choose a
`faQ a` setting whose digit zero is `v` and whose other digits are zero.
The formula for `faOut` in
[Transform.lean](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean)
and the zero extension in `digZ` give its outgoing message as

```math
\mu_v=16^{ab}v\pmod M.
\tag{3}
```

The modulus `M` is odd. Hence multiplication by `16^(ab)` is a bijection
of its residues. Also `2^d<M`, so the `2^d` messages in (3) are distinct.
This step does not assume `M` is prime.

For each `0<=w<2^d`, choose a `mid b` setting whose second-vector
coordinate `a` is `mu_w`, with every other coordinate zero. These are
valid settings because each coordinate is an arbitrary residue. On
these chosen rows and columns the edge tensor is

```math
E(u^{(v)},z^{(w)})=[\mu_v=\mu_w]=[v=w].
\tag{4}
```

Thus the unrestricted original edge tensor contains an identity
submatrix of order `2^d`, and its rank over the rationals is at least
`2^d`. If it factors through a `c`-bit channel,

```math
E(u,z)=\sum_{\theta\in\{0,1\}^c} L(u,\theta)R(z,\theta),
```

then its rational rank is at most `2^c`. Consequently `c>=d`.
Integer factorizations are included by embedding their coefficients
in the rationals. The same conclusion holds for a general finite
channel alphabet: it must have at least `2^d` elements.

The claim is about factoring each original edge before taking products
of edge channels. Rank is not being added across correlated wires.
In particular, this does not assert that the product of several such
edge factors has rank equal to the product of their separate ranks.

## 3. Grouping trades local enumeration against crossing volume

Partition the original gates into any `P>=2` groups. No balance condition
and no use of the original cluster map is required. For group `u`, let
`q_u` and `m_u` count its `faQ` and `mid` gates. Put

```math
q_{\max}=\max_u q_u,\qquad H_Q=hdq_{\max}.
\tag{5}
```

There are `A` gates of each kind, so `sum_u m_u=A`. The number of internal
edges from the chosen layer is exactly `sum_u q_u m_u`, and is at most
`A q_max`. At least `A^2-A q_max` layer edges therefore cross groups.

Let `lambda` be a common clique channel width that accommodates the
concatenation of separate original-edge channels between every pair of
groups. Section 2 forces at least `d` bits for each crossing edge of this
layer. Other crossing edges can only add to the required widths. Hence

```math
\binom P2\lambda\ge d(A^2-Aq_{\max}).
\tag{6}
```

Substituting `A=2h` and `d=2S^2` makes the local-setting contribution
cancel exactly:

```math
\begin{aligned}
H_Q+\frac{P(P-1)\lambda}{4}
&\ge hdq_{\max}+\frac{d(A^2-Aq_{\max})}{2}\\
&=\frac{dA^2}{2}=4S^4\ge n.
\end{aligned}
\tag{7}
```

In particular, for the volume `V=P^2 lambda` used by the evaluator,

```math
H_Q+V/4\ge n.
\tag{8}
```

A group with `q_max` input gates has at least `2^(hd q_max)=2^H_Q`
original local settings. All other gate-setting sets are nonempty.
Visiting every local setting once therefore costs at least `2^H_Q`
steps, even if many settings are subsequently rejected by internal
equalities or zero weights. This lower bound applies to the literal
enumeration procedure, not to every possible method of evaluating the
same local sum.

One consequence of (8) is the simpler dichotomy: either local enumeration
already costs at least `2^(n/2)`, or `P^2 lambda>=2n`. Reducing the number
of crossing edges by concentrating gates has a corresponding enumeration
cost. Better balancing constants cannot remove this tradeoff.

## 4. Consequence for the explicit pad and the fixed level schedule

The completed composition uses the Appendix B schedule

```math
K_0=2^{x_0},\qquad
K_{\ell+1}=(2^{3\ell+x_0}K_\ell)^2,
\qquad P=K_\ell,\quad r=2^{\lfloor\ell/2\rfloor},
\tag{9}
```

where `x0>=18` is fixed. Its declared, computable pad exponent is

```math
q=\left\lfloor\frac{V}{r}\right\rfloor+1+
C_0P^2\mathrm{size}\bigl(P(\lambda+B)\bigr)>V/r,
\tag{10}
```

with fixed nonnegative `C0`. The full algorithm physically writes `2^q`
pad units. Let `T` include all local enumeration and pad output. Then

```math
\log_2 T\ge\max(H_Q,V/r)
\ge\frac{H_Q+V/4}{1+r/4}
\ge\frac{4n}{r+4}.
\tag{11}
```

The middle inequality is just the fact that a weighted average is at
most the maximum of its two entries. It does not rely on a lower bound
for the source theorem's existential function `need`.

Writing `a_ell=log2 K_ell`, (9) gives
`a_(ell+1)=2a_ell+6ell+2x0>=2a_ell`. Consequently

```math
\log_2P\ge x_0 2^\ell\ge x_0r^2.
\tag{12}
```

For `P<=n`, equations (11)--(12) imply the exact lower bound

```math
\log_2T\ge\frac{4n}{4+\sqrt{\log_2n/x_0}}
=\Omega\!\left(\frac{n}{\sqrt{\log n}}\right).
\tag{13}
```

Choosing `P>n` does not evade the bound. A dense `PInst` row contains
`2^((P-1)lambda)` explicitly encoded integers. Since `lambda>=1` and
`P,n` are integers, writing just one such row then costs at least `2^n`.
This also explains why a succinct representation of mostly zero rows
does not suffice for the unchanged dense interface.

The source schedule and dense interface are at the pinned Appendix B
revision
[`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth),
in [LevelsStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/LevelsStatement.lean)
and [PaddedDescent/F/Data.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/PaddedDescent/F/Data.lean).
The local proof and actual encoding of (10) are in
[half-rate-factoring.md](half-rate-factoring.md) and
[half-evaluator-audit.md](half-evaluator-audit.md).

For any proposed score `alpha>1/2`, choose
`epsilon=(alpha-1/2)/2`. The exponent
`n/(log n)^(alpha-epsilon)` is then asymptotically smaller than the
right-hand side of (13). Thus this implementation class cannot satisfy
that score's eventual runtime bound.

## 5. What the obvious table shortcuts do and do not change

**Sparse table writing.** In the existing proof, table construction costs
`2^O(n/P)`, while the pad is bounded by `2^O(n/r)`. The maximal-level
choice gives `n<P^400` and therefore `P>n^(1/400)`. Hence `n/P` is already
lower order than `n/sqrt(log n)`. Reducing the cost of generating those
tables cannot improve the leading logarithmic exponent while their
width and the explicit pad remain unchanged. A sparse intermediate form
followed by decompression still pays for the dense target rows.

**Prefix restrictions.** Every edge minor in Section 2 is present for
every prefix query, since the `faQ` and `mid` settings and their weights
are unchanged. Fixing more of the `faP` input does not shrink these
separate edge tensors. Avoiding some prefix queries reduces the number
of calls, but even one call using the construction above has the stated
cost. Deriving and imposing new restrictions on `faQ` from global
multiplication constraints would change the hypothesis and requires a
separate algorithm and cost proof.

**Adaptive grouping and separate wire compression.** Equations (7)--(8)
hold for every partition, including one chosen using the input or its
prefix. They allow each wire its own best rational factorization. The
obstruction therefore survives all such choices that retain the local
enumeration and product-of-wire channel interface.

**Routes not ruled out.** A joint factorization can use correlations
among several deterministic messages without factoring every edge
separately. A local solver can count compatible settings without visiting
their whole Cartesian product. A new analysis might show that the
unchanged evaluator needs a smaller pad on a family closed under its
recursive transformations. No result here settles any of those options.
They require an additional structural theorem, beyond the 0/1 and
deterministic-message observations used by the completed constructor.

## 6. A finite check of the identity minor, including a composite modulus

The proof above applies to every `s`. As a limited source-definition
check, the following code evaluates the full `faOut` sum with the source's
zero extension. It checks every entry of the proposed identity minor on
every `faQ`--`mid` edge for `s=0,1`:

```sh
python3 - <<'PY'
for s in (0, 1):
    h = 1 << s
    A = 2 * h
    d = 2 * h * h
    modulus = (1 << (4 * d)) + 1
    values = 1 << d
    checked = 0
    for a in range(A):
        for b in range(A):
            messages = []
            for v in range(values):
                digits = [v] + [0] * (h - 1)
                total = sum(
                    (digits[i] if i < h else 0)
                    * pow(16, A * i * b, modulus)
                    for i in range(A)
                )
                messages.append(pow(16, a * b, modulus) * total % modulus)
            # mid second-vector coordinate a is messages[w].
            for v in range(values):
                for w in range(values):
                    assert (messages[v] == messages[w]) == (v == w)
                    checked += 1
    if s == 1:
        assert modulus == 641 * 6700417
    print(f's={s}, modulus={modulus}, wires={A*A}, '
          f'minor_order={values}, entries={checked}: PASS')
PY
```

Both cases passed: 64 entries across four edges at `s=0`, and 1,048,576
entries across sixteen edges at `s=1`. The latter modulus is the composite
integer `4294967297=641*6700417`. This finite check catches a misuse of
field assumptions or the source's digit layout in these cases; the
all-input rank claim follows from the odd-modulus argument, not from
these examples.
