# Endpoint-only affine frames with the complete Gaussian tail

The literal full-tail reduction now has a checked full-rank constructor.
It retains both states of the lowest occupation mode and the entire
Gaussian tail. Its active matrix has order 256 and full generic rank
256, independently of the number of hot modes. The first complete
binary coefficient census finds no improvement over its Gaussian
control.

## 1. Declared family and exact calculation

Let \(M\ge4\) be a power of two. As in
[`tensor-literal-full-tail-core.md`](tensor-literal-full-tail-core.md),
the \(M\) hot modes lie in a reflection interval of length \(4M\).
On the lowest mode, choose arbitrary binary two-by-two matrices
\(A,B\) and set

\[
 L(t)=I+(1+t)A,\qquad P(t)=I+(1+t)B.
\]

Then \(L(1)P(1)=I\). The construction does not require
\(L(t)P(t)=I\) away from the endpoint. The physical alphabet has
dimension \(2^M\). Its lowest-mode heights are \((0,g)\), and
the remaining modes retain heights \(2j+1\), \(1\le j<M\).
The completed census uses \(g=1\).

The [constructor](../tools/tensor_literal_full_tail_affine_frames.py)
expands each weighted low endomorphism into the identity, the number
operator, the creator, and the annihilator. It transports these operators
using the literal eight-mode formulas, clears each local map, and applies
all four maps to the active Gaussian diagonal. Both places are computed
separately. Each receipt retains the total clearing power, the cleared
orders, and the actual Laurent orders after subtraction. At parameter
one, the active matrix is the identity, proving its full generic rank.

The active coefficient degree is at most

\[
 32M+\sum_{i=1}^4 d_i,
\]

where \(d_i\) is the degree of the corresponding cleared low map.
The constructor checks this bound. The native Smith calculation must
return all 256 orders before a profile is accepted. These profiles
have no capped or omitted tail orders.

The spectators and physical tail are retained by their exact products

\[
 \prod_{j=2}^{M-3}(1+z^{2j+1})(1+z^{2M+2j+1}),
 \qquad
 (1+z^g)\prod_{j=1}^{M-1}(1+z^{2j+1}),
\]

respectively. The eight-mode calculation reduces the cost of finding the
complete profile. It does not remove the physical alphabet dimension or
the polynomial-map degree from the eventual transfer cost.

## 2. Independent original-kernel checks

The [checker](../tools/check_tensor_literal_full_tail_affine_frames.py)
contracts the original raw Laurent Pfaffian tensor at \(M=4\).
Its four cases include the Gaussian control, a nonconstant diagonal
endpoint-only map, opposite off-diagonal maps, and a dense map with
assigned height \(g=2\). The last case also checks the placement of the
infinity height gauges.

The [retained receipt](../evidence/tensor-literal-full-tail-affine-controls.json)
records 524,288 complete raw Laurent coefficient comparisons and the
same number of endpoint checks. All 2,048 literal Smith orders agree
with the active constructor. A separate binary Toeplitz elimination,
which does not call the native Smith routine, verifies 556 jet ranks
through the largest order in every case.

Replay without replacing the retained evidence:

```sh
python tools/check_tensor_literal_full_tail_affine_frames.py \
  --output /tmp/tensor-literal-full-tail-affine-controls.json
```

## 3. Complete binary census at sixteen hot modes

There are 16 binary two-by-two matrices. The
[complete census](../evidence/tensor-literal-full-tail-affine-m16-census.json)
tests all 256 ordered pairs \((A,B)\) at \(M=16\), with natural
lowest-mode height \(g=1\). It contains 131,072 complete active orders,
covering both cuts of every pair. The physical alphabet dimension in
each case is 65,536. The entire spectator contribution is included in
every pressure fit.

The best observed fitted value is

\[
 \rho=0.7083157559380195,
 \qquad s=0.16261972089322593.
\]

This is the existing full Gaussian control \(A=B=0\). Eight other
pairs tie it to the reported numerical precision; none gives a smaller
observed fitted value. The preliminary
[16-case pilot](../evidence/tensor-literal-full-tail-affine-m16-pilot.json)
is retained separately.

The fitted parameter searches use a logarithmic grid followed by local
refinement. Their global minima are not certified. The exact part of
the census is the complete pair catalogue, the endpoint identities,
the polynomial matrices, the Laurent clearing costs, and the full
Smith profiles. The Gaussian control already lies below the separate
compact-alphabet benchmark; reproducing it here is not a new pressure
frontier or a stronger factoring result.

To repeat the complete census:

```sh
python tools/tensor_literal_full_tail_affine_frames.py \
  --modes 16 --workers 4 \
  --output /tmp/tensor-literal-full-tail-affine-m16-census.json
```
