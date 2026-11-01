# Exact syndrome ranks for restricted bilinear trace kernels

## Result and scope

The [binary trace kernels](trace-bilinear-conservation-kernels.md) have an
exact rank formula after restricting every external leg to an arbitrary
label subset. The formula counts realized syndromes and needs no Gaussian
elimination. It also handles separate affine controls on the two kernel
variables.

The available finite constant-weight profiles do not supply a crossing
decomposition. This note additionally identifies an exact obstruction to
one proposed synthesis step: ordinary products of two low-degree kernel
matrices vanish. Kernel-matrix multiplication is different from row
multiplication of the four-leg tensors, and no free cut-rank bound for the
former is asserted.

## 1. An exact block formula for every label subset

Let `G=F_2^t`, let `X` be any nonempty subset of `G`, and restrict all four
external labels of the tensor

```math
T_K(x,y,u,v)=[x+y=u+v]K(x+y,x+u)
```

to `X`. For a row-block label `S`, put

```math
L_S=X\cap(X+S).
```

Both free row and column coordinates of that block range over `L_S`.
Suppose the kernel at this `S` has the form

```math
K(S,w)=[B_S(w)=b_S],
```

where `B_S:G->F_2^k` is linear and `b_S` is fixed. Define the realized
syndrome set

```math
A_S=\{B_S(x):x\in L_S\}.
```

**Syndrome-rank theorem.** The rank of this row block is exactly

```math
\boxed{|A_S\cap(A_S+b_S)|.}
\tag{1}
```

**Proof.** Its entry at free labels `x,u` is
`[B_S(x)+B_S(u)=b_S]`. Group rows and columns by their syndromes. A row
syndrome `a` is joined only to the column syndrome `a+b_S`. When both are
present, the corresponding rectangular matrix is all ones and has rank
one. Distinct row syndromes use disjoint column groups. Exactly the
syndromes in (1) contribute. No independence assumption on the coordinates
of `B_S` is necessary. The proof works over any coefficient field because
it only ranks zero and all-one blocks.

If, for fixed `w`, the kernel likewise has a linear-syndrome description
in `S`, the identical argument computes each column block. Both full cut
ranks are sums of the corresponding block counts.

### The symmetric trace kernel

For

```math
K(S,w)=\prod_{j=0}^{k-1}\operatorname{Tr}(\alpha^j Sw),
```

take

```math
\sigma_S(x)=(\operatorname{Tr}(\alpha^j Sx))_{0\le j<k},
\quad b_S=\boldsymbol1.
```

The kernel and common label restriction are symmetric in `S,w`. Therefore

```math
\boxed{
r(T_K|_X)=c(T_K|_X)
 =\sum_{S\in G}|A_S\cap(A_S+\boldsymbol1)|,
\quad A_S=\sigma_S(L_S).
}
\tag{2}
```

The block at `S=0` contributes zero when `k>=1`. For the full label set,
assuming `k<=t` and `alpha` has degree `t`, every other syndrome map is
surjective, and (2) recovers `(2^t-1)2^k`.

Formula (2) can be strictly smaller than the elementary cap
`sum_S min(|L_S|,2^k)`. A syndrome present in the image contributes nothing
unless its prescribed partner is also present.

## 2. Affine controls and an exact finite algorithm

Multiply the kernel by separate binary affine controls `a(S)` and `b(w)`,
each a product of affine linear forms. If `a(S)=0`, discard the row block.
Otherwise append the linear parts of the forms defining `b` to the trace
syndrome. Their required target bits are one plus their affine constants.
The same formula (1) applies. For the column cut, exchange the roles of
the two controls. Thus this extension still needs only syndrome sets;
the controls need not have independent linear parts for this exact formula.

Enumerating ordered pairs `(x,y) in X^2` generates all row-block domains:
set `S=x+y` and insert `sigma_S(x)` in `A_S`. Each element of `L_S` is
generated once at its corresponding ordered pair. A bitset represents
the syndrome set when `k` is modest. This uses `O(|X|^2)` syndrome
evaluations, followed by finite set intersections.

The native implementation uses unordered pairs and inserts both endpoints.
It computes `Sx=x^2+xy` and `Sy=y^2+xy` with one field multiplication.
Field log tables are constructed and checked explicitly; they are an
implementation optimization, not a rank assumption. It stores up to six
syndrome bits in a 64-bit set. Smaller syndrome images are obtained by
projection, and translation by the all-one syndrome reverses the bitset.

For a constant-weight label set `X={x:|x|=m}`, a second exact count checks
the generated block sizes. If `|S|=2s`, then

```math
|L_S|=\binom{2s}{s}\binom{t-2s}{m-s},
\tag{3}
```

with value zero when the indices are inadmissible. Odd-weight `S` have
empty blocks. Equation (3) follows by choosing half of `S` and the
remaining `m-s` positions outside `S`.

## 3. The current constant-weight evidence

The optional screen uses the lexicographically first irreducible binary
polynomial at each width, its polynomial basis for Hamming weight, and
`alpha=X` in that field. The evidence records the exact defining polynomial.
It covers the six pairs

```text
(t,m)=(12,3),(12,4),(12,5),(14,4),(14,5),(14,6),
```

with every `k=1,...,6`. These are 36 individual tensor profiles.

| `t,m,k` | `|X|` | Exact rank in each cut | Block-size cap |
|---|---:|---:|---:|
| `12,4,4` | 495 | 30550 | 31696 |
| `14,5,5` | 2002 | 256964 | 259200 |
| `12,3,6` | 220 | 15236 | 46528 |

At `k` near the chosen label weight, the observed ranks remain close to
the cap. The third row shows a substantial syndrome-image reduction, but
its exact rank still exceeds `|X|^sqrt(2)`. Some low-`k` pieces are below
that individual threshold. None of these observations supplies the missing
sum-to-crossing identity or bounds the cost of such a synthesis.

## 4. Ordinary kernel compositions and their vanishing range

This section works over a coefficient field of characteristic two.
Write `K_f(S,w)=f(Sw)` for a full-label trace-product kernel. On nonzero
labels, expansion in powers gives

```math
K_f(S,w)=\sum_n c_n S^n w^n.
```

Multiplying two such **ordinary kernel matrices** pairs exponent `n` from
one kernel with exponent `-n` from the other. Their product therefore has
opposite multiplicative exponents, which are aligned with those in the
equality kernel. This explains why composition is a possible algebraic
operation to investigate. It does not make it a cheap tensor operation:
the intermediate kernel label has size `2^t`, and its contraction needs
an independently proved resource and cut-rank bound.

There is also an exact vanishing condition. If the two trace products
have respectively `k` and `l` factors and `k+l<t`, then

```math
\boxed{K_fK_g=0.}
\tag{4}
```

An elementary Boolean-polynomial proof avoids spectral assumptions. For
fixed outer labels, the summand in the intermediate variable `w` is a
product of `k+l` binary linear forms. Its degree is below `t`. The sum of
any such polynomial over `F_2^t` is zero: each monomial misses a coordinate
and hence occurs an even number of times. In particular, a bare trace
kernel with `2k<t` is square-zero as an ordinary matrix.

More generally, for a diagonal intermediate control `h(w)` of Boolean
degree at most `a`,

```math
K_f\,\operatorname{diag}(h)\,K_g=0
\quad\text{when }k+l+a<t.
\tag{5}
```

This follows from exactly the same degree argument. Separate controls
on the outer variables do not change the vanishing. Equations (4)--(5)
concern the full additive label cube. Restricting the intermediate sum
to a smaller subset can remove the cancellation, so no unchanged
restriction theorem is asserted for those products.

## Reproduction

```sh
python3 tools/check_trace_weight_restriction.py --output /tmp/trace-weight-restriction-verification.json
python3 tools/check_trace_weight_restriction.py --full --output /tmp/trace-weight-restriction-verification-full.json
```

The checker requires C++17. The small gate compares native syndrome counts
with independently constructed four-leg cuts, including a six-bit syndrome
to check the full bitset range. It also checks the ordinary-composition
identities directly in small fields. The optional full gate regenerates
all 36 profiles and verifies their block-size caps using (3).
