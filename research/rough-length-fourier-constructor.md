# A shorter Fourier constructor after small-prime removal

## Result and scope

The [small-word Fourier constructor](small-word-fourier-constructor.md)
can use a transform with half the earlier leading information volume.
The construction first removes small prime factors and then branches over
the exact bit lengths of the two remaining factors. Its carry-column
grouping satisfies

```math
P^2\lambda_{\rm col}\le(12+o(1))n
\qquad\text{when }P=o(\sqrt A).                 \tag{1}
```

The previous constant for that grouping was `24`. This is a constructor
constant, not a new logarithmic saving exponent. The complete achieved
factoring score remains `1/2`.

All preprocessing is charged. For a scalar clique evaluator with time
`E(P,lambda)` and workspace `W(P,lambda)`, the resulting complete factoring
bound, over a fixed coefficient field, is

```math
\begin{aligned}
T_{\rm fact}(n)
&\le(E(P,\lambda)+1)\operatorname{poly}(n)
       \exp(O(n/P+P)),\\
S_{\rm fact}(n)
&\le W(P,\lambda)+O(n/P+P+\log n).
\end{aligned}                                             \tag{2}
```

For the carry-column layout assume `P^2<=A`; then its additional workspace
is `O(n/P)`. A fixed polynomial range `n^epsilon<=P<=n^delta`, with
`0<epsilon<=delta<1/4`, makes the added logarithmic time smaller than
`n/sqrt(log n)`. The evaluator's own admissible parameters and costs must
still be checked when applying (2).

## 1. Parameters and the rough cofactor

Treat sufficiently small input lengths directly. For the remaining `n`,
put

```math
\mu=\lceil\log_2(n+2)\rceil,\qquad
g=\lceil\log_2(\mu+2)\rceil.
```

Let `A>=2` be the largest power of two with `A^2 mu g<=n`. Set

```math
T=A^2,\quad H=T/2,\quad
d=\left\lceil\frac{n+1}{T}\right\rceil,
\quad B=2^d.
\tag{3}
```

Choose a power of two `2<=P<=A`, write `L=A/P`, and set the integer
small-prime threshold exponent

```math
s=\frac{Td}{P}=ALd.                              \tag{4}
```

The useful elementary bounds are

```math
\frac{n}{4\mu g}<T\le\frac{n}{\mu g},\qquad
n+1\le Td\le n+T,\qquad
d=\Theta(\mu g),\qquad s=\Theta(n/P).             \tag{5}
```

Remove all prime powers supported on primes at most `2^s`. Let `R` be
the remaining cofactor, represented through its read-only bit generator.
The [output-sensitive transfer](output-sensitive-space-transfer.md)
constructs those bits in time `exp(O(s+log n))` and workspace
`O(s+log n)`, without storing an `n`-bit cofactor. It also produces the
removed prime factors. These costs are included in (2).

If `R>1`, every prime divisor of `R` exceeds `2^s`. Let `m=bits(R)<=n`.
For each proper factor pair `xy=R`, both factors exceed `2^s`. Their exact
bit lengths `k,l` therefore satisfy

```math
s+1\le k,l\le m-s,
\qquad k+l\in\{m,m+1\}.                          \tag{6}
```

The upper bound follows from the lower bound for the opposite factor and
`k+l<=m+1`. Enumerate the at most `2m` ordered length pairs satisfying
(6). Each proper factor pair occurs in exactly one branch. If the list
is empty, there is no proper factor pair; in particular, roughness proves
`R` prime whenever `R>1` and `bits(R)<=2s`.

The two trivial pairs `(1,R)` and `(R,1)` are handled explicitly. When
`R=1`, handle the single pair `(1,1)` instead. Thus weighted divisor
counts, prefixes, and bit pins retain their full original meaning.

## 2. The shorter transform is exact in every length branch

For one branch, let

```math
H_x=\lceil k/d\rceil,\qquad H_y=\lceil l/d\rceil.
```

The input digit arrays have these respective lengths. Require the highest
bit of each factor to be one, and the unused high bits of its last digit
to be zero. These are local conditions at the appropriate input gate.

For positive integers `k,l,d`,

```math
\lceil k/d\rceil+\lceil l/d\rceil-1
\le\lceil(k+l)/d\rceil.
```

Consequently (5)--(6) give

```math
H_x+H_y-1\le T,\qquad
\min(H_x,H_y)\le H=T/2.                          \tag{7}
```

A length-`T` cyclic convolution therefore equals the ordinary integer
convolution. Its coefficients satisfy

```math
c_j=\sum_{u+v=j}x_u y_v,
\qquad 0\le c_j\le H(B-1)^2.                     \tag{8}
```

Use the explicit CRT modulus construction from Section 9 of the companion
note: multiply successive primes `q_i=1 mod T` until their product `M`
exceeds `TB^2`, and combine exact `T`th roots componentwise. This gives
polynomial-time metadata construction and `f=ceil(log_2 M)`-bit arithmetic.
The original two-dimensional DFT identities remain valid. Only the number
of digits at an input gate changes: gate `FA_x[r]` sends

```math
X_{r,u}=\omega^{ru}\sum_{\substack{i\ge0\\r+Ai<H_x}}
                  x_{r+Ai}\xi^{iu},\qquad 0\le u<A,
\tag{9}
```

and `FA_y[r]` uses `H_y` in the same formula. Middle and inverse gates
are unchanged. Since `M>TB^2`, their canonical inverse-transform residues
are exactly (8).

The previous carry cap also remains valid:

```math
C_*=H(B-1),\qquad
c_j+v_j\le H(B-1)^2+C_*=BC_*.
```

Every candidate in a length branch has
`xy<2^(k+l)<=2^(n+1)<=B^T`, so its true final carry is zero. This argument
uses the bit-length restrictions; the rounded digit lengths alone can
sum to `T+1`. Requiring the usual digit equations against `R`, initial
carry zero, and final carry zero gives exactly `xy=R`, by telescoping.
There is one arithmetic extension of each pair and no extension of a
false pair. Summing the branch contractions therefore counts each proper
factor pair once.

## 3. Small-prime removal restores unique boundary interpolation

Shortening the transform without restricting factor lengths would break
the earlier uniqueness proof: an input row could have `A` coefficients
but only `A-L` outputs prescribed outside its group. This is why the
preprocessing in Section 1 is part of the construction.

Here, (4)--(6) instead give

```math
k,l\le n-s<Td-s=(T-T/P)d.
```

It follows that

```math
H_x,H_y\le T-T/P=A(A-L).                          \tag{10}
```

Each residue-class row `r+A i` therefore contains at most `A-L` input
digits. The group containing that input gate has exactly `L` middle gates,
so its boundary prescribes `A-L` transform outputs. Select as many of
those outputs as the actual row length. After removing the twiddle in
(9), they determine the row polynomial uniquely.

The interpolation points have unit differences modulo every CRT prime,
as proved in the companion note. Interpolation over `Z/MZ` is therefore
valid. A row of length zero is the unique empty input row. Reject a
recovered row if a digit range, high-bit condition, prescribed extra
output, or requested factor bit condition fails.

Once all local input rows are fixed, the unchanged directed reconstruction
determines the local middle, inverse, and carry settings. Thus every
prescribed group boundary has at most one valid setting. A table query
has polynomial time relative to bit access to `R`. It interpolates one
row at a time, with `O(Af)` scratch bits and a constant-depth chain of
regenerated transform outputs. Substituting the virtual cofactor bit
generator contributes the time and space already charged in (2).

## 4. Information volume and the constant twelve

The edge topology is unchanged: two complete input/middle layers, one
middle/inverse layer, the inverse/coefficient edges, and the carry chain.
The two input families together use `O(n)` setting bits; all other gate
and edge bounds remain `O(Tf)`. The largest gate has mass `O(n/A)`.

The CRT proof gives, for one fixed effective Linnik exponent `ell`,

```math
f\le2d+(\ell+1)\log_2T+\ell\log_2R_0+O(1),       \tag{11}
```

where `R_0` is the auxiliary power of three in that metadata proof.
For the parameters (3), the sufficient number of CRT primes and `R_0`
are both `O(g)`. Hence

```math
Tf\le2Td+O(T\log n)
\le2n+O(n/g)=(2+o(1))n.                          \tag{12}
```

Place each carry `iA+r` in the group containing inverse gate `IB[r]`.
The exact edge count and padding from the companion note give

```math
\lambda_{\rm col}=(6L^2+e_P)f,\qquad
e_P=\begin{cases}2A-1,&P=2,\\A,&P\ge4,\end{cases}
```

and therefore

```math
P^2\lambda_{\rm col}=6Tf+P^2e_Pf.                \tag{13}
```

For `P=o(sqrt(A))`, the second term is `o(Tf)`. Equations (12)--(13)
prove (1). The full boundary and query scratch occupy

```math
O(Tf/P+PAf+Af)=O(n/P)\quad\text{if }P^2\le A.    \tag{14}
```

For completeness, the original contiguous-carry grouping also works.
Its width is `(8L^2+1)f`, it has `O(n/P)` query space for every `P<=A`,
and it gives the constant `16+o(1)` when `P=o(A)`.

The carry cap has not changed. With `r=log_2H+d`, every carry is below
`2^r`, while `f-r>=d+2`. The scalar phase control fits the existing carry
slot on the row-zero path between consecutive groups. No additional
boundary width is needed for that control.

## 5. Complete preprocessing and recovery cost

The scalar oracle sums at most `2m` branch evaluations and adds the two
explicit trivial contributions. A branch index, the exact bit length `m`,
and all transform metadata use polynomial construction time and
`O(n/P)` workspace. The two trivial weights are evaluated by streaming
the bits of `R`; no full factor word is retained.

Roughness gives `log_2 tau(R)<=n/s<=P`. The uniform isolation/recovery
argument in the output-sensitive transfer therefore has control size
`O(P+log n)` and overhead `exp(O(s+P+log n))`. The cofactor bit generator
has overhead `exp(O(s+log n))`; substituting it into the constant-depth
table routines changes only the absolute constant in that exponent.
The polynomial number of length branches is also absorbed. Together
with `s=Theta(n/P)`, these facts prove (2), including complete prime
factorization and least-prime-factor output.

This does not assert that the new preprocessing is a practical improvement.
It can enumerate an enormous small-prime range. Its purpose here is the
fully charged asymptotic constructor bound and the exact interpolation
condition (10).

## Exact replay and provenance

Run

```sh
python3 tools/check_rough_length_fourier_constructor.py \
  --output evidence/rough-length-fourier-constructor.json
```

The checker covers exact branch capacities, variable-length transforms,
integer convolution and carries, external-output interpolation, rough
cofactor divisor coverage, and concrete CRT metadata. A negative control
exhibits distinct valid digit rows with identical external messages when
the new capacity bound is omitted. The finite controls supplement the
proof; they do not establish its asymptotic costs by benchmarking.

Factor-length splitting is also discussed in
[the p-adic carry note](padic-digit-carry-transfer.md). The added point
here is that charged small-prime removal makes the shorter transform
compatible with the existing unique-boundary interpolation and grouping.
The construction builds on the multiplication-network reduction of Isaac
Rudich and Louis-Martin Rousseau and the two companion transfer notes.
It used substantial OpenAI assistant support. No new Lean formalization
or factoring score improvement is claimed.
