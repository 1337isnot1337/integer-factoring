# Graded braid channels from the invertible three-strand object

The three-strand Jones--Wenzl projector gives a simple way to vary the
central vacuum rank throughout the height-three braid quotient. Replace
the central layer by disjoint f_3 triples and leave k strands free. Its
exact vacuum rank is

```math
D_k=\begin{cases}1,&k=0,\\F_{k-1},&k\ge1,\end{cases}            \tag{1}
```

where F_j is Fibonacci. With the same marked boundary seam as before,
the degree-b four-tensor has actual cut bounds

```math
\boxed{\quad r(T_b)\le h_{n,b}D_k,\qquad
c(T_b)\le\binom nb v_b,\quad}
\qquad h_{n,b}=\min(b+1,n-b+1).                 \tag{2}
```

Here v_b is the vacuum dimension on 2b interface strands. This extends
the [triplet middle channels](fusion-projector-seams.md): the row
exponent can now be set anywhere from zero to two, up to polynomial
prefactors, while retaining the same column-interface bound.

These are candidate rank factors. The currently implemented prefix
placements miss the five- and thirteen-label crossing even before a
cost filter. Neither (1) nor the joint envelope proves a complete
crossing decomposition or a score above `1/2`.

## 1. A local projector with invertible fusion charge

Use the order-five F16 parameter and
[negligible-wall quotient](quantum-braid-path-quotient.md) of the
earlier braid family. Let `delta=q+q^-1`. In characteristic two,

```math
f_3=I+e_0+e_1+\delta^{-1}(e_0e_1+e_1e_0).       \tag{3}
```

For example the local Jones--Wenzl recursion gives
`f_2=I+delta^-1 e_0` and `f_3=f_2+f_2 e_1 f_2`, since
`[2]_q=[3]_q=delta`. Expanding, using `e_0e_1e_0=e_0` and
`delta^-1+delta^-2=1`, gives (3). The same TL relations show
`f_3 e_i=e_i f_3=0` for i=0,1 and `f_3^2=f_3`.

The retained simple labels are0,1,2,3; **label0 is the tensor unit**.
The local fusion resolutions give

```math
0\otimes1=1,\quad 1\otimes1=0\oplus2,\quad
2\otimes1=1\oplus3,\quad 3\otimes1=2.
```

The discarded label4 is negligible. The retained resolutions have
unit Gram factors, so these are actual direct fusion channels in the
quotient. Associativity, or iterating the same resolutions, gives

```math
3\otimes3=0,\qquad3\otimes2=1,\qquad
2\otimes2=0\oplus2.
```

Thus write `g=3`, `tau=2`: g is invertible with `g^2=unit`,
`tau^2=unit+tau`, and label1 is `g tau`. These multiplicity identities
use integer channel counts; two copies are not canceled just because
the coefficient field has characteristic two.

## 2. Exact central rank

Suppose disjoint triples occupy3t physical strands and k strands are
left free, with `3t+k=2n`. The image of the product projector is the
fusion object

```math
g^t(g\tau)^k=g^{t+k}\tau^k=\tau^k.
```

The last equality holds because `t+k` is even. Iterating
`tau^2=unit+tau` gives

```math
\tau^k=F_{k-1}\,\mathrm{unit}\oplus F_k\,\tau
\quad(k\ge1).
```

The vacuum multiplicity is therefore D_k in (1). The product of
disjoint local idempotents projects onto exactly this image; unit
fusion norms make the multiplicity equal to its matrix rank on the
vacuum quotient. This is not a claim about the ordinary ambient
rank of the physical f_3 blocks.

The multiplicity does not depend on the placement of the disjoint
triples. The compiler uses one explicit placement: free prefixes of
lengths l,r in the two n-strand halves, followed by f_3 triples.
Require `0<=l,r<=n` and `n-l,n-r` divisible by three. Then k=l+r.
The whole central projector acts inside the two physical blocks.
Local block vacuum dimensions can vanish even when its global
vacuum rank is nonzero; the implementation does not assume otherwise.

## 3. Braid coefficients and both cuts

Use the anti-diagonal block-braid circuit. Replace its entire central
layer by this projector. On the n marked cross-boundary cells use
`qI+y e_(n-1)`; use `qI+e_i` on the other outer cells. Apply the
fixed scalar `q^(n-n^2/2)` and project to the usual paired vacuum
boundaries. This defines `T(y)=sum_b T_b y^b` over F16. No strict
endpoints are assigned to this auxiliary polynomial.

Each outer half has y-degree at most n/2. Taking coefficient b
therefore gives h_(n,b) outer coefficient pairs. Every row factor
passes through the central image of dimension D_k, proving the
first bound in (2).

Choose the b marked boundary gates that supply e instead of the
scalar qI. There are `binom(n,b)` choices. For each choice, cutting
between the two physical blocks opens2b strands. The
[vacuum-interface lemma](fusion-projector-seams.md#1-vacuum-interfaces-instead-of-physical-bit-channels)
factors the resulting column matrix through a space of dimension v_b.
The central f_3 blocks remain block-local, so that lemma applies
without another auxiliary factor. This proves the second bound.

Choose the same common exact binary alphabet at its least even
physical length n. Its size L obeys `L=Theta(phi^n)`, and cropping
cannot increase either rank. The path compiler, coefficient entries
and exact dense factors have allowed `2^O(lambda)` cost. The field
is fixed; the earlier whole-tensor descent handles its extension.

## 4. Programmable grades and an explicit cost filter

Let `s=n mod3`, `n=s+3m`. The allowed free counts are
`s,s+3,...,n`. Their two-block sums include every
`2s+3a`, `0<=a<=2m`. Thus for any fixed `0<=theta<=2` one can
choose l,r uniformly with `l+r=theta n+O(1)`. Equation (1) gives

```math
r(T_b)\le\operatorname{poly}(n)L^\theta.
```

If `b=eta n+O(1)`, the other exponent is bounded by

```math
\min\!\left(2,\frac{H(\eta)}{\ln\varphi}+2\eta\right),
\qquad H(u)=-u\ln u-(1-u)\ln(1-u).              \tag{4}
```

This is an upper profile, not a lower bound or a statement that every
coefficient is nonzero. It provides a wider collection of graded
candidates than a few isolated edge and middle regimes.

There are only `O(n^3)` prefix candidates indexed by l,r,b. For a
fixed rational `gamma=p/q`, retain those whose explicit bounds
R,C in (2) satisfy

```math
(RC)^q\le L^{2p}.                              \tag{5}
```

This filter uses integers of polynomial bit length. Include the cut
reflections interchanging Q,R, as well as the two cheap pairing
tensors. Reflection exchanges the cut ranks and preserves SWAP.
The retained collection has rank cost at most
`poly(n)L^gamma`. If `1<gamma<sqrt(2)` and SWAP belongs to its span
at every sufficiently large binary width, the existing strict
post-half transfer applies. Moving gamma toward one would approach
the score-one transfer endpoint, **conditional on that unproved
spanning statement**.

Exact linear solving for the span costs `2^O(lambda)` because the
candidate count is polynomial and tensor entries are effective.
Expose a solution for `SWAP+I` in distinct positive monomial degrees
with constant coefficient I. Alternatively append candidate pieces
by endpoint-vanishing terms before recombination. This satisfies the
strict crossing interface if a favorable spanning identity is found.

The prefix placements do not currently supply such an identity.
At full alphabets of5 and13 labels, even the **unfiltered** span
with both cut orientations and the two pairings misses SWAP.
Its exact dimensions are80 and146. No cost filter can repair those
finite failures. The two-label span contains SWAP, without passing
the universal finite rank-sum criterion.

A new placement/grouping identity or a fused kernel is needed. The
rank formula allows other disjoint placements, but enumerating all
tilings would introduce exponentially many pieces. A successful
grouping must retain polynomially many pieces and prove the ranks
of their sums; the small rank of each individual projector does not
give that automatically. Wider unchanged prefix scans are not the
remaining research target.

## 5. Reproduction

[The compiler](../tools/quantum_braid_graded_seam.py) and
[checker](../tools/check_quantum_braid_graded_seam.py) give exact
F16 controls. The [receipt](../evidence/invertible-fusion-graded-seams.json)
records global projector ranks and internal annihilations, complete
small coefficient profiles, independent scalar propagation and the
unfiltered span failures.

```sh
make verify-braid-graded-seam
```

NumPy is needed for the finite matrix replay. This extends the earlier
braid quotient, vacuum interfaces and fusion-projector channel work.
The strongest complete factoring score remains `1/2`, and no new
Lean formalization or independent mathematical acceptance is claimed.
