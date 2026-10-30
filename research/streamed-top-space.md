# Sublinear workspace for the original top-table queries

## Theorem and scope

Use the source parameters and the explicit grouping in
[polynomial-top-queries.md](polynomial-top-queries.md):

```math
A=2S,\qquad h=S,\qquad d=2S^2,\qquad
B=2^d,\qquad M=2^{4d}+1,\qquad f=4d+1.
```

Here `S` is the least power of two for which `n<=4S^4`.
Let `P` be a power of two with `2<=P<=A`, and write `L=A/P`.
An original gate belongs to group `a` exactly as specified in that
note: the four families indexed by `r` use the block containing `r`,
and carry `j=qA+r` uses the block containing `q`.

For each group and each prescribed, padded boundary, its table value
can be computed in polynomial time using

```math
O\!\left(\frac nP+Af+\log(n+2)+W\right)                 \tag{1}
```

work-tape bits. The boundary is included in this bound. The integer
`N`, and any specified prefix, are read-only inputs. The term `W` is
the workspace of an on-demand bit-weight routine. A constant number of
coefficient-ring elements is also required; in the half-rate application
the ring has sixteen elements, so those elements have constant size.
The running time is polynomial in `n` times the bit-weight query time.

Weights may be placed on **either or both factor-input families**.
In particular a character may depend on the first factor while the
single square-zero tag marks a bit of the second factor. This does not
change the uniqueness argument or the space bound.

The existing top query already has a unique valid local setting when
it is nonzero. The improvement here is a deterministic implementation
that reconstructs parts of that setting on demand. It stores neither
all original gate settings nor a matrix of interpolation coefficients.
No new crossing identity or time exponent is claimed.

The source pin, gate formulas, and composite-modulus unit proof are the
ones recorded in the top-query note. Everything below is an ordinary
algorithmic argument, without new Lean formalization.

## 1. Store the packed boundary; generate its addresses

An edge between two distinct groups carries at most

```math
\lambda=(8L^2+1)f
```

bits. Thus one queried group has a packed boundary of length

```math
(P-1)\lambda
\le (8A^2/P+P)f
\le 9A^2f/P
=O(n/P).                                                   \tag{2}
```

The last inequality uses `P<=A` and `A^2 f<=144n`.
Padding is checked by scanning its positions.

Read a requested bit of `N` or the prefix directly from its read-only
tape. Construct at most one `d`-bit input digit when a carry check needs
it. In particular, do not obtain a bit by copying or shifting an entire
`n`-bit input integer into a work register. The implementation exposes
this distinction with an immutable byte-tape accessor.

There is no stored list of original gates, boundary edges, or offsets.
Generate edges in four nested-loop families:

1. The two input families to the middle family.
2. Middle to inverse.
3. Inverse `r` to carry `qA+r`.
4. Consecutive carries.

The membership test is integer division of a gate index by `L`, after
division by `A` for a carry. To locate a requested boundary message,
scan this fixed edge order, retain just its ordinal among edges to the
specified other group, and compute its packed offset. The padded block
for each other group has the known length `lambda`. This uses
`O(log(n+2))` index bits and polynomial time. All endpoints use the same
order. Scanning instead of storing addresses is sufficient for the
claimed polynomial time bound.

## 2. One input gate uses `O(hf)` interpolation storage

Each contained input gate has at least `h` prescribed outputs to
external middle gates. Choose the first `h` of them. After removing
the source twiddle, the points and values are

```math
x_k=\xi^k,\qquad v_k=\omega^{-rk}y_k\pmod M,
\qquad \xi=\omega^A.
```

Every difference of two selected points is a unit modulo `M`, by the
geometric identity in the top-query note. Lagrange interpolation gives

```math
F(X)=\sum_i v_i
 \left(\prod_{j\ne i}(x_i-x_j)^{-1}\right)
 \prod_{j\ne i}(X-x_j).                                  \tag{3}
```

Maintain one result coefficient array of length `h` and one current
product coefficient array of length `h`. Multiply a product by
`X-x_j` in place, updating coefficients from high degree to low degree.
The denominator and its inverse use a constant number of `f`-bit
registers. Points and values may be stored in two further length-`h`
arrays or recomputed. In either implementation, storage is `O(hf)`.
There is no `h` by `h` matrix.

Euclidean inversion modulo `M`, modular multiplication, and modular
powering use `O(f+log(n+2))` bits with schoolbook arithmetic. Every
arithmetic operand here has `O(f)` bits. Interpolation remains
polynomial time, and all inversions are of proved units.

Read the recovered coefficients in `[0,M)` and reject a coefficient
outside `[0,B)`. Check the applicable input-prefix condition by
reading its corresponding bits. Evaluate the polynomial to check all
remaining prescribed outputs of this input gate.

The recovered digits can then be discarded. If a later middle-gate
query needs one of this input's outputs, repeat the interpolation.
Its result is the same because the prescribed boundary determines the
input uniquely.

## 3. Recompute arithmetic messages through a constant-depth graph

Use three procedures, each returning one arithmetic message:

| Procedure | Data retained while it runs |
| --- | --- |
| `input_output(kind,r,k)` | One input interpolation, then one evaluation |
| `middle_output(k,r)` | Two length-`A` input vectors and transform accumulators |
| `inverse_output(r,q)` | One length-`A` middle-output vector and an inverse-transform accumulator |

A missing predecessor inside the group is replaced by its prescribed
boundary message. Every incoming residue is checked to lie in
`[0,M)`. When a predecessor is inside, the procedure calls the
preceding row of the table. All transforms use the original source
formulas.

There are at most three live arithmetic procedure frames: inverse,
middle, and input. Consequently all their arrays together use
`O(Af)` bits. Repeated queries may recompute the same interpolation
and transform many times, but every loop has polynomial length and
the dependency depth is constant. Their total time remains polynomial
in `n`. No memoization is used.

Validate every contained input, middle, and inverse gate in separate
passes. In each pass recompute its outgoing messages to external
gates and compare them with the boundary. Incoming boundary messages
are used verbatim and checked for legal range. Internal messages are
generated by the same deterministic procedures at both endpoints, so
these checks are exactly the original local-equality checks.

## 4. Stream the contiguous carry segment

Group `a` contains precisely the consecutive carries

```math
aLA\le j<(a+1)LA.
```

Its initial carry is zero when `a=0`, and otherwise is the prescribed
incoming boundary carry. At position `j=qA+r`, obtain the coefficient
from `inverse_output(r,q)` when that inverse gate is inside the group;
otherwise read the boundary coefficient. Let `c` be the current carry.
Check the source carry condition

```math
(\text{coefficient}+c)\bmod B
   =\text{the }j\text{th base-}B\text{ digit of }N,
```

and replace `c` by the quotient by `B`. Check all carry ranges. At the
end compare with the outgoing boundary carry, or require zero if the
global final carry is inside this group.

Only one carry is retained. Its residue queries invoke the
constant-depth arithmetic procedures of Section 3. Thus the carry
chain adds time, but it adds only `O(f+log(n+2))` storage. In particular,
its length `AL` does not multiply the workspace.

This also explains why no recursive call is made to a predecessor
carry: the chain is scanned once in its natural order.

## 5. Accumulate either factor's weights once

Use one designated pass over the contained input gates to accumulate
the coefficient-ring weight. For each recovered input digit, scan its
bits and query the specified weight for its family and global bit
position. Multiply these weights into a single running accumulator.
Calls made later solely to reconstruct arithmetic messages do not
modify this accumulator.

This gives the monomial weight of the unique local setting exactly
once. It works for weights on the first factor, the second factor,
or both, including zero weights. If any local check fails, return zero.
If every check passes, return the accumulated weight.

For the fixed-ring isolation algorithm, computing a requested character
phase and one tag bit takes `O(k+log(n+2))` control bits. Hence `W` in
(1) can be that bound. No array of phases or witness bits is needed.

The uniqueness theorem proves correctness: any valid local setting
must equal every recovered input and every propagated message.
Conversely, the checks cover every local gate condition, boundary
condition, and internal equality. Therefore a passing setting is
valid and has exactly the returned weight.

## 6. Size at the scheduled top level

Since `A=Theta(n^(1/4))` and `f=Theta(n^(1/2))`,

```math
Af=O(n^{3/4}).                                             \tag{4}
```

Choose the largest scheduled top size `P=K_ell<=A`, with the schedule
from the relaxed transfer. Its proved growth bound
`K_(ell+1)<=K_ell^10` gives `P>A^(1/10)` once the schedule exists.
Thus `P=n^Omega(1)` and

```math
\frac nP+Af+\log(n+2)
=o\!\left(\frac n{\sqrt{\log(n+2)}}\right).               \tag{5}
```

This supplies the improved top-query component for a sublinear-space
half-rate evaluator. Equation (5) is only a top-query claim: recursive
metadata, local factors, and witness recovery require their own
storage accounting.

## 7. Exact checks

The implementation in
[streamed_top_queries.py](../tools/streamed_top_queries.py) contains
only scalar parameters, a packed boundary tape, and the constant-depth
procedures above. Its interpolation uses in-place polynomial products.
The checker compares it with the independently stored-setting
implementation from the original top-query checker. It covers valid
boundaries, malformed messages, nonzero padding, prefix rejection,
and separate weights on both factor-input families.

```sh
python3 tools/check_streamed_top_space.py \
  --output evidence/streamed-top-space.json
```

The finite checks test exact gate transcriptions and branch coverage.
The workspace and polynomial-time bounds follow from the explicit
array sizes and constant-depth recomputation argument, not from
runtime memory measurements on those finite examples.
