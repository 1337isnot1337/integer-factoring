# A rectangular Fourier constructor with balanced carry relays

## Result

The [rough-length constructor](rough-length-fourier-constructor.md) can
retain its leading constant twelve over a larger parameter range.
A rectangular transform and deterministic equality relays give

```math
\lambda=\frac{6T}{P^2}f+4(r+2),\qquad
P^2\lambda\le\left(12+4\frac{P^2}{T}+o(1)\right)n,
\tag{1}
```

uniformly for every power of two `2<=P<=sqrt(T)`. Here `T=A^2` is the
transform length from the rough-length construction, `f` is the arithmetic
residue width, and `r` is the carry width defined below.

Thus `P^2 lambda<=(12+o(1))n` whenever `P=o(sqrt(T))`. In particular it
holds for every fixed range `n^epsilon<=P<=n^delta`, with
`0<epsilon<=delta<1/2`. The previous column grouping required `delta<1/4`.
At the largest permitted `P`, the uniform constant is at most `16+o(1)`.
The full table-query workspace is `O(n/P)` throughout this range.

For a uniform scalar clique evaluator with time `E(P,lambda)` and
workspace `W(P,lambda)`, the complete transfer remains

```math
T_{\rm fact}(n)\le(E(P,\lambda)+1)\operatorname{poly}(n)
                         \exp(O(n/P+P)),\qquad
S_{\rm fact}(n)\le W(P,\lambda)+O(n/P+P+\log n).
\tag{2}
```

Every preprocessing and reconstruction cost is included. This improves
the constructor parameters; it supplies no new local crossing evaluator.
The achieved logarithmic saving exponent and factoring score remain `1/2`.

## 1. Parameters and exact factor lengths

Use the uniformly constructed metadata from the rough-length note:

```math
\mu=\lceil\log_2(n+2)\rceil,\quad
g=\lceil\log_2(\mu+2)\rceil,\quad
A=\max\{2^j:(2^j)^2\mu g\le n\},\quad T=A^2.
```

Set `d=ceil((n+1)/T)`, `B=2^d`, and `H=T/2`. Small input lengths are
handled directly. Build the CRT arithmetic ring `Z/MZ` with `M>TB^2`
and a componentwise primitive `T`th root `omega`. Put

```math
f=\lceil\log_2M\rceil,\qquad r=\log_2H+d.
```

Every true carry is at most `H(B-1)<2^r`. The CRT construction proves

```math
Tf\le(2+o(1))n,\qquad T(r+2)\le(1+o(1))n,
\qquad r+2\le f.                                     \tag{3}
```

The second inequality uses `Td<=n+T` and
`T log n=O(n/g)=o(n)`. The error bounds are independent of `P`.

For a power of two `2<=P<=A`, remove all prime factors at most `2^s`,
where `s=Td/P`. For the resulting virtual rough cofactor `N`, enumerate
the exact bit lengths `k,l` of every proper pair `xy=N`. Their base-`B`
lengths satisfy

```math
H_x+H_y-1\le T,\qquad \min(H_x,H_y)\le H,\qquad
H_x,H_y\le T-T/P.                                   \tag{4}
```

The exact high-bit conditions give `xy<B^T`, even when the rounded digit
lengths sum to `T+1`. The trivial factor pairs are handled explicitly.
For `P=2`, roughness already leaves no proper branch; the network formulas
remain well-defined but are unnecessary for that case.

## 2. The rectangular transform

Choose

```math
R=T/P,\qquad S=P,\qquad T=RS,\qquad
\zeta=\omega^R,\quad\eta=\omega^S.
```

There are `R` input gates for each factor, `S` middle gates, and `R`
inverse gates. The input families send

```math
X_{a,b}=\omega^{ab}\sum_{i:a+Ri<H_x}x_{a+Ri}\zeta^{ib},
\qquad 0\le a<R,\quad0\le b<S,                       \tag{5}
```

and the analogous messages `Y`. Middle gate `b` computes

```math
\widehat x_{b,v}=\sum_{a<R}X_{a,b}\eta^{av},\qquad
\widehat y_{b,v}=\sum_{a<R}Y_{a,b}\eta^{av},\qquad
Z_{b,a}=R^{-1}\sum_{v<R}\widehat x_{b,v}\widehat y_{b,v}\eta^{-av}.
\tag{6}
```

Inverse gate `a` recovers

```math
c_{a+Ri}=S^{-1}\sum_{b<S}Z_{b,a}\omega^{-ab}\zeta^{-ib}. \tag{7}
```

The full transform at frequency `b+Sv` expands with digit `a+Ri` as
`omega^(ab) eta^(av) zeta^(ib)`. This proves (5)--(7) are exactly the
length-`T` DFT. The inverses of `R,S` exist modulo every CRT prime.
By (4), cyclic convolution is ordinary convolution, with each integer
coefficient below `HB^2<M`. Canonical residues recover these coefficients
exactly. The carry equations and telescoping proof are unchanged: each
factor pair has one arithmetic completion, and every valid completion
has product exactly `N`.

## 3. Grouping and unique interpolation

Put middle gate `b` in group `b`. Split the `R` input and inverse gates
into `P` consecutive intervals of size `R/P=T/P^2`. Place carry `a+Ri`
with inverse gate `a`. All divisions are integral because `P<=A` and
the parameters are powers of two.

Each input row has at most `P-1` coefficients, since (4) gives
`H_x,H_y<=R(P-1)`. Exactly `P-1` transform outputs go to other groups,
so the prescribed boundary determines the row polynomial. Distinct
powers of `zeta` have unit differences at every CRT prime; interpolation
over `Z/MZ` is valid. Remove the twiddle, interpolate, then check digit
ranges, high bits, all prescribed outputs, and requested bit pins.
An empty row is uniquely zero.

Once these rows are fixed, the local middle outputs and inverse
coefficients follow from (6)--(7). Each group has `S=P` disjoint carry
segments. A segment starts at its prescribed boundary carry or the
initial zero, and ends at a prescribed carry or the final zero.
Consequently every boundary assignment has at most one local setting.

The three complete bipartite arithmetic layers contribute exactly
`6RS/P^2=6T/P^2` words to each unordered pair of distinct groups.
All inverse-to-coefficient messages are internal.

## 4. Exactly balanced carry routing

Before rerouting, each carry row crosses the cycle of group boundaries
`0,1,...,P-1,0`; the final row omits its last wrap edge. More generally,
let there be `S` rows with `P|S`. For row `i`, choose relay group
`c=i mod P`. Replace each external edge from group `a` to its cyclic
neighbor `b` by `a -> c -> b`. If `c` is an endpoint, retain the direct
edge. Otherwise insert a degree-two equality vertex in group `c`.

An equality vertex has exactly one setting when its two words agree.
It therefore preserves the original contraction and factor-pair
bijection, without an additional multiplicity.

Fix two distinct groups `a,b`. They receive routed segments only from
rows with relay `a` or `b`. In either case, the other endpoint has two
cycle neighbors and contributes two segments per row. There are `S/P`
rows for each relay. The exact segment count is

```math
\frac{4S}{P}-[\{a,b\}=\{0,P-1\}].                    \tag{8}
```

The subtraction is the omitted final wrap, whose relay is `P-1`, so it
would have been direct. Formula (8) also holds for `P=2`, where the
two cycle boundaries are distinct directed occurrences.

For `S=P`, every clique pair carries at most four carry words. There
are exactly `S(P-2)<=T` new equality vertices and `P-2` per group.
Their local checks compare two prescribed boundary words. They introduce
no recursive reconstruction or additional choices.

## 5. Scalar phase controls and narrower carry words

The scalar recovery interface needs at most three phase states. Reserve
two bits above the `r` carry bits, giving width `r+2<=f`. Reject carries
above `H(B-1)` and invalid phase codes. Outside the phase path, require
both extra bits to be zero.

The phase follows the original row-zero chain of groups along its new
routed edges. Each group adds its first-factor input-bit contribution
exactly once at its original row-zero segment. Transit relays copy the
entire packed word and add nothing. The path revisits group zero through
distinct relay vertices; these visits are not repeated phase updates.
The final group applies the terminal scalar weight before the row-zero
wrap. This is a single deterministic control path for each factor pair.

The local query obtains its incoming state from the boundary or initial
pin and streams its contribution while regenerating input rows. It
therefore implements the same scalar phase count as the earlier
output-sensitive recovery construction.

## 6. Width, time, and workspace

Pack `6T/P^2` arithmetic words of width `f` and at most four carry words
of width `r+2` on each clique edge. Force unused padding bits to zero.
This proves the exact width in (1); equation (3) gives its uniform
asymptotic bound after multiplication by `P^2`.

The full boundary has `O(Tf/P+P(r+2))=O(n/P)` bits, because `P^2<=T`.
Interpolation uses `O(Pf)` scratch bits within the same bound. Regenerate
middle and inverse outputs by nested loops and running residues. A
middle output invokes only input-row interpolation; an inverse output
invokes only middle outputs and input rows. The depth is constant.
Process one carry segment at a time and stream the relay equality checks.

Every loop range is polynomial in `n`. Thus queries have polynomial bit
time relative to the virtual rough-cofactor and phase-input generators.
The CRT metadata have the charged polynomial cost from the companion
note. The largest middle gate has mass `O(Rf)=O(n/P)`; an input or inverse
gate has mass `O(Pf)<=O(n/P)`. Total information mass remains `O(Tf)=O(n)`.

Small-prime removal and recovery cost `exp(O(s+P+log n))`, with workspace
`O(s+P+log n)` and `s=Theta(n/P)`. Substituting their bit generators into
the fixed-depth query routines proves (2). At `P=A`, the added workspace
is `O(sqrt(n log n log log n))`. No evaluator at that point is assumed.

## 7. Exact controls and attribution

The [checker](../tools/check_routed_rectangular_fourier.py) compares the
rectangular transform with an independent full DFT and integer convolution.
It reconstructs every group from its boundary, including exact input
lengths, arithmetic outputs, carries, relays, and the complete three-state
phase path. Wrong targets and corrupted relays must be rejected. Separate
topology controls check (8), including `P=2` and unequal rectangles.

```sh
python3 tools/check_routed_rectangular_fourier.py \
  --output evidence/routed-rectangular-fourier-controls.json
```

The [receipt](../evidence/routed-rectangular-fourier-controls.json) records
54 routing parameter cases, 223,974 segments, 109,368 relay vertices,
36 transform pairs, 264 complete group reconstructions, and 72 rejection
controls. These finite controls store small arrays for comparison; the
regeneration proof establishes the asymptotic workspace bound.

The multiplication-network reduction is due to Isaac Rudich and
Louis-Martin Rousseau. This note builds on the repository's rough-length,
CRT Fourier, and output-sensitive recovery proofs. The rectangular
layout, balanced equality routing, and combined bounds were developed
with substantial OpenAI assistant support. No new Lean formalization
is used.
