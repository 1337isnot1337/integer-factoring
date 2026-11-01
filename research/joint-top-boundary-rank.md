# Joint ranks of the complete raw group-pair boundary

## Result and scope

Arbitrary joint factorization of all original edges between two groups does
not remove the direct top-constructor obstruction. This extends the
[separate-wire result](half-top-structure.md): the factors may combine the
input, middle, inverse-transform, and carry edges between each group pair.
No separate factor for the `faQ`--`mid` layer is required.

The statement still concerns a particular representation. Its group-pair
tensor is the product of the original intergroup equality constraints on
the full Cartesian setting domains. Vertex weights and internal equalities
remain in the group-local factors. The algorithm visits every original
local setting and retains the explicit dense clique interface and pad of
the [half-rate composition](half-rate-factoring.md). It may choose arbitrary
groups and arbitrary signed factors over any field.

For this class, the complete construction cost is

```math
2^{\Omega(n/\sqrt{\log n})}.
```

This is not a lower bound for factoring or for a boundary relation after
local constraints have been imposed. Solving a group's internal equations
first, absorbing its vertex weights into a new joint representation, or
using a different evaluator remains outside the theorem. No larger
factoring score is established.

## 1. Source definitions and composite-ring interpolation

Use Rudich and Rousseau's multiplication network at revision
[`c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1),
with its original parameters

```math
h=2^s,\qquad A=2h,\qquad d=2h^2,\qquad
B=2^d,\qquad M=2^{4d}+1,\qquad n_0=Ahd=4h^4\ge n.
\tag{1}
```

The source's `pL` is denoted by `d`. A `faQ a` setting consists of `h`
independent digits in `{0,...,B-1}`. A `mid b` setting has two independent
vectors in `(Z/MZ)^A`. The `faQ a` message to `mid b` is

```math
16^{ab}\sum_{i=0}^{h-1}u_i\zeta^{ib}\pmod M,
\qquad \zeta=16^A,\qquad\zeta^h=-1.
\tag{2}
```

Every difference of distinct powers `1,zeta,...,zeta^(A-1)` is a unit
modulo the composite integer `M`. Here is a direct proof. For every odd
prime `p` dividing `M`, the identity `zeta^h=-1` holds modulo `p`.
Since `A=2h` is a power of two and `-1` is not `1` modulo `p`, the
multiplicative order of `zeta` modulo `p` is exactly `A`. Thus no such
prime divides `zeta^i-zeta^j` for `0<=i<j<A`, so that difference is
coprime to `M`. The square Vandermonde determinant on any distinct
selected points is a product of units and is therefore a unit itself.

Consequently, for any set `J` of `t` recipient indices, varying the first
`min(t,h)` digits and setting the rest to zero gives at least

```math
B^{\min(t,h)}
\tag{3}
```

different message tuples. For `t<=h`, use the square Vandermonde matrix on
those recipients. For `t>=h`, use any `h` of them. Coordinate multipliers
`16^(ab)` are units and preserve injectivity.

The new step below embeds this identity minor in the complete raw
group-pair tensor.

## 2. Neutralize every other edge type at once

Let `U,V` be disjoint groups of original gates. Write `q_U,q_V` for their
numbers of `faQ` gates and `m_U,m_V` for their numbers of `mid` gates. Let
`E_(U,V)` be the zero-one matrix indexed by their full Cartesian setting
spaces, whose entry is the product of all original edge equalities joining
the two groups. Local weights and equalities wholly within a group are
not included.

**Theorem.** Over every coefficient field,

```math
\boxed{\quad
\mathrm{rank}E_{U,V}\ \ge\
2^{d\{q_U\min(m_V,h)+q_V\min(m_U,h)\}}.
\quad}
\tag{4}
```

The proof uses the original `mLink`, `mV`, and `mS` in
[`MulNet.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean),
and `midOut`, `ibOut` in
[`Transform.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean).

Restrict settings on both sides as follows:

- Every `faP` digit is zero.
- The first input vector of every middle gate is zero.
- Every inverse-transform gate setting is zero.
- Both coordinates of every carry gate setting are zero.

All these are permitted elements of the original typed setting domains.
The source formula for `midOut` is bilinear in its two incoming vectors.
Thus the zero first vector makes every outgoing middle message zero,
independently of its second vector. The inverse-transform messages are
then zero. A zero carry setting sends zero in both directions.

Every intergroup equality except `faQ`--`mid` now equals one: `faP` sends
zero, middle first coordinates are zero, middle-to-inverse messages are
zero, and inverse/carry messages are zero. Hence the restricted matrix is
exactly the product of the two directed `faQ`--`mid` interfaces.

For each `faQ` gate in `U`, choose the message family (3) on the `m_V`
recipients in `V`. The matching coordinates of the middle second vectors
are independent for different `faQ` gates. Selecting one column for each
resulting tuple gives an identity minor of size
`B^(q_U min(m_V,h))`. The reverse direction uses the `faQ` settings in `V`
and the second-vector coordinates of middle gates in `U`; these are
disjoint variables. It supplies the transpose of another identity minor.
Their tensor product is an identity matrix of the order in (4).

An identity minor has the same rank over every field. The ring `Z/MZ`
defines the finite message map; it is not being used as a field for the
rank calculation. This proves (4) without a prime-modulus assumption and
without relying on the rank of a factor before it is composed.

The zero `faP` restriction need not satisfy a divisor-prefix weight, and
the zero carry restriction need not satisfy the fixed target `N`. That
does not affect a rank statement about the raw edge matrix. It is exactly
why this proof does not apply once those weights have been absorbed.

## 3. An exact global tradeoff with no rank-estimate error

Partition all gates into `P>=2` groups. Put

```math
H=\max_u\log_2 |\mathcal S_u|,\qquad
C=\sum_{u<v}\log_2\mathrm{rank}E_{u,v},
\tag{5}
```

where `S_u` is the full Cartesian set of original settings in group `u`.
The same coefficient field is used for the pair ranks. Let `q_u,m_u`
count `faQ` and middle gates in group `u`, and put `q_max=max q_u`.
Both count sums are `A`.

Then

```math
\boxed{\quad H+C/2\ge n_0\ge n.\quad}
\tag{6}
```

If every `m_u<=h`, (4) yields

```math
C\ge d\sum_{u\ne v}q_um_v
  =d\left(A^2-\sum_uq_um_u\right)
  \ge d(A^2-Aq_{\max}).
\tag{7}
```

Also `H>=hd q_max`. Since `A=2h`, the dependence on `q_max` cancels:

```math
H+C/2\ge hdq_{\max}+\frac d2(A^2-Aq_{\max})
        =\frac{dA^2}{2}=Ahd=n_0.
\tag{8}
```

If some group contains more than `h` middle gates, its Cartesian domain
has at least `M^(2A(h+1))` settings. Therefore

```math
H\ge2A(h+1)\log_2M>8Ad(h+1)>8n_0,
\tag{9}
```

which proves (6) in the remaining case. Empty groups cause no problem.

This argument uses only the coarse exact interpolation minor, so there
is no packing loss, no asymptotic rank error, and no restriction such as
`P<=h`.

Suppose the algorithm factors each complete `E_(u,v)` through one channel
of at most `lambda` bits, with arbitrary field-valued factors. Matrix rank
is at most its channel cardinality. Thus

```math
C\le\binom P2\lambda\le P^2\lambda/2.
```

With `V=P^2 lambda`, equation (6) gives

```math
\boxed{\quad H+V/4\ge n.\quad}
\tag{10}
```

## 4. Cost under the unchanged dense interface and explicit pad

Visiting every original local setting costs at least `2^H` steps. The
half-rate composition uses `P=K_ell`, `r=2^floor(ell/2)` and physically
writes a sufficient pad with exponent `q_pad>V/r`. Hence its complete
construction time `T` satisfies

```math
\log_2T\ge\max(H,V/r)
 \ge\frac{H+V/4}{1+r/4}
 \ge\frac{4n}{r+4}.
\tag{11}
```

The source schedule gives `log2 P>=x0 r^2`, where the fixed integer
`x0>=18`. If `P<=n`,

```math
\log_2T\ge\frac{4n}{4+\sqrt{\log_2n/x_0}}
=\Omega\!\left(\frac n{\sqrt{\log n}}\right).
\tag{12}
```

If `P>n`, the unchanged dense interface with `lambda>=1` already writes
at least `2^((P-1)lambda)>=2^n` entries in one row. Thus the conclusion
holds without an additional upper restriction on `P`.

This is an output/enumeration bound for the stated implementation class.
It does not show that the sufficient pad is mathematically necessary, or
that a local solver must enumerate the original settings.

## 5. Replay and remaining route

Run the standard-library checker:

```sh
python3 tools/check_joint_top_boundary_rank.py \
  --output evidence/joint-top-boundary-rank.json
```

It reconstructs literal source messages for every edge type and checks
complete raw-boundary identity minors on finite partitions. It includes
both directions between the groups, nontrivial middle/inverse/carry edge
placements, and the composite modulus `2^32+1`. Separate exact integer
checks cover the global tradeoff for complete count compositions.
Negative controls impose a local input restriction and demonstrate that
the raw identity minor can disappear, preserving the theorem's scope.
The companion [evidence](../evidence/joint-top-boundary-rank.json) records
the exact counts and outcomes.

The remaining compression target is therefore more specific than joint
linear algebra on all raw edges of one group pair. It must use internal
compatibility or target-dependent weights before that interface is fixed,
preserve correlations across several group pairs, or replace the current
explicit evaluation cost. A constructor that solves local constraints
before encoding a boundary falls outside the stated hypotheses.

The transform and preceding separate-wire argument retain their source
attribution. During this work, local unpublished research also supplied
the fan-out interpolation observation; no public provenance is asserted
for that local material. All facts used here have self-contained proofs
above or exact pinned upstream definitions. This extension and its checker
were developed with OpenAI assistant support. No Lean source is changed.
