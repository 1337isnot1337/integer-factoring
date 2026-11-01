# An explicit top clique with polynomial-time table queries

## Theorem and scope

The original multiplication network has a stronger top-table property
than the local-setting enumeration bound uses. With a suitable grouping,
each prescribed boundary has **at most one valid internal setting**.
The setting can be recovered by interpolation and forward propagation.

Use the pinned source's parameters

```math
S=2^s,\quad A=2S,\quad h=S,\quad d=2S^2,\quad
B=2^d,\quad M=2^{4d}+1,\quad f=4d+1,
```

with `s` minimal such that `n<=4S⁴`. Let `P` be a power of
two with `2<=P<=A`, and put `L=A/P`.

There is a deterministic, explicitly grouped `P`-vertex clique with

```math
\lambda=(8L^2+1)f,\qquad
P^2\lambda=(8A^2+P^2)f\le9A^2f\le1296n.
\tag{1}
```

Its contraction is the original divisor-prefix count, or the weighted
bit-multivariate divisor polynomial. All metadata are generated in
polynomial time. **One table entry is computed in polynomial bit time
and polynomial space**, apart from the separately charged cost of
polynomially many operations on the target coefficient algebra.

The unweighted table values are zero or one. A weighted entry is zero
or the monomial weight of its unique valid local setting.
No grouping search or enumeration of local gate settings is required
for a query.

This is an oracle representation of the tables. Writing every entry
still requires time proportional to the number of entries, which is
generally exponential. The streamed evaluator already proved in
[rough-input-transfer.md](rough-input-transfer.md) can use these
procedures directly.

The theorem strengthens the initial constructor. It does not by itself
improve the logarithmic saving exponent beyond one half: the later grid
recursion still uses its existing crossing profile.

Primary source: Isaac Rudich and Louis-Martin Rousseau's
[factoring reference at `c91cc3f4947caca8977e81a5d14244f73b5ae2b1`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/tree/c91cc3f4947caca8977e81a5d14244f73b5ae2b1).
The exact gates and messages are in
[`MulNet.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/MulNet.lean),
and the arithmetic transform is in
[`Transform.lean`](https://github.com/IsaacRudich/SubexponentialDeterministicFactoring/blob/c91cc3f4947caca8977e81a5d14244f73b5ae2b1/FactoringAux/Transform.lean).
The constructor below is an ordinary argument, with no new Lean proof.

## 1. Partial Fourier outputs determine an input gate

Put `omega=16 modulo M` and `xi=omega^A`.
The source has

```math
\omega^{A^2/2}=-1,\qquad
\xi^{A/2}=-1.
\tag{2}
```

An input gate `FA[r]` contains `h=A/2` digits
`u_0,...,u_(h-1) in [0,B)`. Its output to `MID[k]` is

```math
y_k=\omega^{rk}\sum_{i=0}^{h-1}u_i\xi^{ik}\pmod M.
\tag{3}
```

Any `h` distinct outputs in (3) determine the digits uniquely.
The point is that this interpolation works over the **composite**
residue ring `Z/MZ`; it does not assume `M` prime.

Take distinct `k,l in [0,A)`. Write their nonzero difference modulo
`A` as `2^v o` with `o` odd. If `z=xi^(k-l)`, then

```math
z^{A/2^{v+1}}=-1.
```

The finite geometric identity gives

```math
(1-z)\sum_{j=0}^{A/2^{v+1}-1}z^j=2.
\tag{4}
```

Because `M` is odd, two is a unit modulo `M`. Thus `1-z`
is a unit, and so is `xi^k-xi^l`.
The Vandermonde determinant at any distinct selected points
`xi^(k_1),...,xi^(k_h)` is consequently a product of units.

Divide the selected `y_k` by the unit twiddle `omega^(rk)`.
Lagrange interpolation with these points recovers the unique
polynomial of degree less than `h`. Its coefficients are the
residues of the `u_i`. Read each coefficient in `[0,M)` and
reject if any coefficient is at least `B`. Since `B<M`, there is
at most one compatible digit vector.

All inverses required by Lagrange interpolation are known to be units.
The integer Euclidean algorithm computes them in polynomial bit time.
Alternatively, (4) gives an explicit geometric-sum inverse.
This avoids an unjustified field-style pivoting argument over `Z/MZ`.
The matrix order is `h=O(n^(1/4))` and its entries have
`O(sqrt(n))` bits, so interpolation costs `poly(n)`.

## 2. A unique-extension lemma for arbitrary gate clusters

Let `U` be any subset of the original gates containing at most
`h=A/2` middle gates. Prescribe the message on every edge leaving
`U`. Then there is at most one choice of the original gate settings
inside `U` that satisfies all internal equalities and these boundary
messages. Its existence and weight can be computed in polynomial time.

For each input gate inside `U`, at least `A-h=h` middle gates
are outside `U`. Its messages to those middle gates are prescribed
boundary values. Apply Section 1 to any `h` of them, then check all
remaining prescribed outputs, digit ranges, and input-prefix conditions.
Every contained input gate has now been determined or rejected.

Propagate the other gates in their directed order:

1. A middle setting is its two incoming vectors. Each entry is either
   the computed output of a contained input gate or a prescribed
   boundary value from an external input gate.
2. An inverse-gate setting is its incoming middle-output vector,
   determined in the same way.
3. A carry setting consists of its incoming coefficient and incoming
   carry. The coefficient comes from its inverse gate or the boundary.
   The incoming carry is zero at the initial position, a previously
   computed carry if the predecessor is inside `U`, or its prescribed
   boundary value otherwise.

Reject an incoming arithmetic value outside `[0,M)`, an illegal
carry value, any failed prefix or digit condition, or any mismatch on
an outgoing boundary message. Check the final zero-carry condition if
the final gate is contained in `U`.

The directed dependency graph is acyclic: input gates precede middle
gates, which precede inverse gates, which precede the ordered carry
chain. Thus the procedure determines every setting. Any valid setting
must follow exactly these interpolation and propagation steps, proving
uniqueness. Conversely, passing every local weight, internal equality,
and boundary check gives a valid setting.

In the unweighted network its contribution is one. For the weighted
network, multiply the input-bit weights present in the recovered
first-factor digits. This is exactly the contribution that the original
local-setting sum would have returned. A zero substituted weight is
allowed.

The condition on middle gates is substantive. If a cluster contains an
input gate and every middle gate, no outgoing Fourier value of that
input is prescribed. Taking the second-factor boundary values to be
zero makes all middle outputs zero regardless of those first-factor
digits. Such a boundary can have many valid local settings.

## 3. Explicit grouping and exact edge counts

Partition `[0,A)` into the `P` contiguous blocks

```math
I_a=\{aL,\ldots,(a+1)L-1\},\qquad0\le a<P.
```

Place each of `FA_P[r],FA_Q[r],MID[r],IB[r]` in group `a`
when `r in I_a`. For a carry gate with index

```math
j=qA+r,\qquad0\le q,r<A,
```

place it in the group containing `q`. Thus each group contains
`L` middle gates and `AL` carry gates. Since `P>=2`,
`L<=h`, so Section 2 applies to every group.

For distinct groups `a,b`, count the original edges crossing
between them:

- The two input families to the middle family contribute `4L²`.
- Middle to inverse contributes `2L²`.
- Inverse to carry contributes `2L²`. Indeed, `IB[r]` connects
  to the carry indexed `qA+r`, so each ordered pair of index blocks
  supplies `L²` such edges.
- The carry chain contributes one edge when the groups are consecutive,
  and zero otherwise. All other carry-chain edges stay inside a group.

Consequently the exact unordered boundary width is

```math
Y_{ab}=\bigl(8L^2+[|a-b|=1]\bigr)f.
\tag{5}
```

Concatenate these original messages in a common deterministic order
at both endpoints, then append zero padding to length
`lambda=(8L²+1)f`. Reject a query with nonzero padding.
The coarsening identity and the off-diagonal codec are exactly those
proved in [half-top-constructor.md](half-top-constructor.md),
Sections 2 and 4. Replacing its local enumeration by Section 2
preserves every table entry.

The grouping is explicit. There is no coloring search or existence
oracle. Formula (5) can be generated directly, and all original
boundary-edge lists have polynomial size.

## 4. Information volume and costs

Minimality of the source parameter gives `S⁴<=4n` for `n>=2`.
Since `S>=1`,

```math
A^2f=4S^2(8S^2+1)
=32S^4+4S^2\le36S^4\le144n.
```

Equation (1) follows immediately.
There are `A²+4A=O(sqrt(n))` original gates. Their settings,
messages, transforms, and ordinary integer arithmetic have polynomial
descriptions and computation costs. Recovering and storing all settings
for one group uses polynomial space. Multiplying its at most `n`
relevant bit weights uses polynomially many coefficient-algebra
operations.

Therefore metadata generation and every table query are polynomial
time and space. A written table has
`2^((P-1)lambda)` entries per node, so its explicit output still
costs exponential time. The query theorem does not erase that output
cost; it allows a streamed evaluator to avoid writing the table.

For the ordinary integer count the entry-height parameter can be taken
to be one, since the entries lie in `{0,1}`. This does not bound the
global count by one: many different shared edge labellings can satisfy
the resulting clique.

## 5. A larger top level fits the same streamed recursion

The schedule `K_j` used by
[subexponential-piece-transfer.md](subexponential-piece-transfer.md)
consists of powers of two. Choose the largest level `ell` with
`P=K_ell<=A`. Then `P` divides `A`, so the explicit grouping
above applies. As before, finitely many smaller inputs can use trial
division.

Now `A=Theta(n^(1/4))` and `K_(ell+1)<=P^10`, giving
`log P=Theta(log n)` and
`rho^ell=Theta((log n)^(-log2(1/rho)))`.
This permits top levels up to `Theta(n^(1/4))` instead of the old
sufficient restriction `P^40<=n`.

To retain the uniform label lower bound in the relaxed-family proof,
replace its volume cutoff `sqrt(n)` by `n^(3/4)`:
stop when `V^4<=n³`, using integer arithmetic.
Since `P²<=A²<=8sqrt(n)`, every node that performs a grid step has

```math
\lambda_j=\frac{V_j}{K_j^2}>
\frac{n^{3/4}}{P^2}\ge\frac18n^{1/4}.
\tag{6}
```

The rounding contribution has relative size

```math
\frac{K_{j-1}^2}{V_j}
\le\frac{P}{n^{3/4}}
=O(n^{-1/2}).
\tag{7}
```

Thus the same bounded product of width errors and the same geometric
entropy sum apply. The parent-query series
`sum_j V_j/K_(j-1)=O(n rho^ell)` is unchanged.
The grouping searches inside later levels cost
`exp(O(n^(1/8)(log n)^2))`, and an early leaf costs at most
`exp(O(n^(3/4)))`. Both are below every fixed logarithmic-power
bound considered there. The initial table-query cost is now polynomial.

This establishes a fully explicit starting instance for the streamed
recursion with the same rate and space conclusions. It does not
change the leading profile contraction `rho` or the square-growth
schedule, so it gives no new logarithmic saving exponent on its own.

## 6. Exact verification

The companion checker tests Vandermonde inversion without assuming a
prime modulus, recovers source input digits from arbitrary selected
half-size output sets, checks the exact grouped edge counts, and compares
interpolation-based cluster queries against an independent enumeration
of input-gate digits in small cases. It includes malformed boundaries
and a control with too many internal middle gates.

These checks validate the transcription. The ordinary proofs above
establish uniqueness and the uniform polynomial query bound.
