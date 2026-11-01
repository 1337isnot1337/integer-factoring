# Binary trace kernels with small crossing ranks and large ordinary rank

## Result and scope

There are explicit binary conservation tensors whose two crossing ranks are
small relative to the ordinary matrix rank of their kernels. They disprove
the proposed inequality

```math
r(T_K)c(T_K)\ge d\,\operatorname{rank}(K)^2
\tag{1}
```

even when the coefficient field is `F_2`.

A [smaller exact coefficient certificate](trace-kernel-small-certificate.md)
subsequently gives `t=19,k=7`, the first possible extension degree in this
trace-product family. The symbolic parameters proved here are `t=33`, `k=11`, and `d=2^33`.
The kernel is binary and has

```math
r(T_K)=c(T_K)=(2^{33}-1)2^{11},\qquad
\operatorname{rank}(K)\ge\binom{33}{11}=193536720.
\tag{2}
```

These bounds contradict (1) by an exact integer comparison. In particular,
an individual tensor can have both cut ranks below `d^(4/3)` while its
ordinary kernel rank is larger than `d^(5/6)`. A decomposition of the
crossing into such tensors has **not** been constructed. This note proves a
scope counterexample and supplies candidate pieces; it does not improve the
factoring score.

The stronger inequality with `d^2 rank(K)` was already excluded for general
kernels by the linear example in
[the diagonal-word note at the pinned PR #27 revision](https://github.com/1337isnot1337/integer-factoring/blob/44832e5316a4bd8bd6caebe962c31aa6897fb77d/research/diagonal-word-rank-sum-barrier.md).
The present example addresses the weaker expression (1), which would itself
have implied a general `d^(3/2)` rank-sum obstruction by subadditivity of
ordinary kernel rank. It shows why that extension is unavailable.

## 1. The binary construction and its two cuts

Identify the additive label group with the field `L=F_(2^t)`. Fix an element
`alpha` of degree `t` over `F_2`, and choose `1<=k<t`. Define

```math
f(z)=\prod_{j=0}^{k-1}\operatorname{Tr}_{L/\mathbf F_2}(\alpha^j z),
\qquad
K(S,w)=f(Sw)\in\mathbf F_2.
\tag{3}
```

Field multiplication appears only in the formula for `K`; conservation and
the two convolutions use the additive label group. The four-leg tensor is

```math
T_K(x,y,u,v)=[x+y=u+v]K(x+y,x+u).
\tag{4}
```

The row and column cuts split into additive-convolution blocks:

```math
r(T_K)=\sum_S\operatorname{rank}C_{K(S,\cdot)},\qquad
c(T_K)=\sum_w\operatorname{rank}C_{K(\cdot,w)},\qquad
C_g(a,b)=g(a+b).
\tag{5}
```

For `S!=0`, the `k` binary linear forms
`w -> Tr(alpha^j S w)` are linearly independent. Indeed,
`1,alpha,...,alpha^(k-1)` are independent, multiplication by `S` is
invertible, and the trace pairing is nondegenerate. For completeness, the
trace polynomial `X+X^2+...+X^(2^(t-1))` is nonzero and has degree less than
`2^t`, so it does not vanish on the whole field. Multiplication by any
nonzero element then proves nondegeneracy of the pairing.

Thus `K(S,·)` is the indicator of one affine coset of a subspace of
codimension `k`. Convolution by that indicator has rank exactly `2^k`:
after grouping rows and columns by the quotient cosets, its matrix consists
of `2^k` disjoint nonzero all-one blocks, up to a permutation. This rank
calculation holds even in characteristic two. The block at `S=0` vanishes.
Symmetry of (3) gives the identical column calculation. Therefore

```math
\boxed{r(T_K)=c(T_K)=(2^t-1)2^k.}
\tag{6}
```

The use of multiplication in `L` is a uniform, explicit rule for binary
tensor entries. Given an irreducible polynomial defining `L`, an entry can
be evaluated with polynomially many operations on `t` bits.

## 2. The ordinary kernel rank: a Vandermonde certificate

The zero row and zero column of `K` vanish. On `L^*`, write the function
`f` as its unique polynomial

```math
f(z)=\sum_{n=0}^{2^t-2}c_n z^n,
\quad\text{with exponents reduced modulo }2^t-1.
\tag{7}
```

The evaluation matrix `V=(S^n)_(S in L^*,0<=n<2^t-1)` is square and
invertible: it is a Vandermonde matrix on distinct nonzero field elements.
The nonzero-label submatrix of `K` is `V diag(c_n) V^T`. Consequently

```math
\operatorname{rank}(K)=|\{n:c_n\ne0\}|.
\tag{8}
```

This calculation takes place over `L`, but the matrix `K` has binary
entries. Its rank is unchanged by extending the coefficient field, so (8)
also computes its ordinary rank over `F_2`.

Expand the `k` trace factors in (3):

```math
\operatorname{Tr}(\alpha^j z)
  =\sum_{i=0}^{t-1}(\alpha^{2^i})^j z^{2^i}.
```

For a set `A` of `k` distinct positions, let `n_A=sum_(i in A)2^i`.
The coefficient at this exponent is

```math
c_{n_A}
 =\det\bigl((\alpha^{2^i})^j\bigr)_{
                 0\le j<k,\ i\in A}
 =\prod_{i<i',\ i,i'\in A}
                  (\alpha^{2^{i'}}-\alpha^{2^i})\ne0.
\tag{9}
```

The permanent from assigning the `k` distinct positions to the trace
factors equals the determinant in characteristic two. The conjugates
`alpha,alpha^2,...,alpha^(2^(t-1))` are pairwise distinct because `alpha`
has degree `t`. The usual Vandermonde determinant therefore does not vanish.

### Why repeated positions do not contribute to (9)

This point includes the reduction of exponents modulo `2^t-1`. Regard an
assignment to the trace factors as `k` tokens placed in `t` cyclic binary
positions. Two tokens at one position may be replaced by one token at the
next position; position `t-1` is followed by position zero. This replacement
preserves the exponent modulo `2^t-1` and decreases the number of tokens
by one. It must terminate with at most one token at each position.

Since the initial number is `k<t`, the terminal set is a proper subset of
the positions. Distinct proper subsets have distinct residues modulo
`2^t-1`, except for the empty/full-set pair; the full set is unavailable
here. If any position was repeated initially, at least one carry occurs,
so the terminal set has fewer than `k` positions. It cannot have residue
`n_A`, whose canonical proper subset has exactly `k` positions. Thus only
the bijective assignments contribute to (9).

All `binom(t,k)` coefficients in (9) are nonzero, proving

```math
\boxed{\operatorname{rank}(K)\ge\binom tk.}
\tag{10}
```

## 3. An explicit finite counterexample

Take `t=33,k=11` and define `L` using

```math
X^{33}+X^{10}+1.
```

The checker proves this polynomial irreducible by testing every monic
constant-one polynomial of degree at most 16 for divisibility. Any proper
factorization of a degree-33 polynomial would have such a factor. There
are exactly 65,535 candidate divisors. Set `alpha` to the residue of `X`.
This completely specifies the binary kernel (3), without storing its
exponentially large matrix.

The integer comparison is

```math
\binom{33}{11}^2
 =37456461988358400
 >36028797018963968
 =2^{33+2\cdot11}.
\tag{11}
```

Equations (6), (10), and (11) give

```math
r(T_K)c(T_K)<2^{2t+2k}
             <2^t\binom tk^2
             \le d\,\operatorname{rank}(K)^2,
```

contradicting (1). The checker also verifies the strict inequality using
the exact ranks `(2^t-1)2^k`, not just their power-of-two upper bounds.

More generally, along `k/t -> c` with `0<c<1/2`, (6) has logarithmic
size `(1+c)t+o(t)`, while (10) has size at least
`H_2(c)t-O(log t)`. The proposed inequality fails exponentially whenever
`H_2(c)>1/2+c`; `c=1/3` is one such choice. This asymptotic observation
concerns the counterexample family alone.

## 4. Controlled pieces retain the same product bound

The construction permits useful asymmetric variants. Let `a(S)` be the
product of `a` binary affine linear forms in `S` with independent linear
parts, and let `b(w)` be the analogous product of `b` forms in `w`.
Their supports have sizes `d/2^a` and `d/2^b`. Put

```math
K'(S,w)=a(S)b(w)K(S,w).
```

For each supported `S`, its row is a product of at most `k+b` affine
linear forms in `w`; if consistent it is a coset indicator, with
convolution rank at most `2^(k+b)`. Otherwise it is zero. The analogous
column statement uses `k+a` forms. Hence

```math
r(T_{K'})\le d\,2^{k+b-a},\qquad
c(T_{K'})\le d\,2^{k+a-b},\qquad
\sqrt{r(T_{K'})c(T_{K'})}\le d\,2^k.
\tag{12}
```

Ordinary coordinate monomials are a special case. The support discount
is essential in (12). Arbitrary polynomials of the same degrees need not
have either the stated support size or the stated convolution-rank bound.

These controlled pieces can have larger total Boolean degree than the
bare family while retaining (12). Whether they synthesize the crossing
with a competitive total rank sum remains open here.

## 5. The remaining synthesis problem

In convention (4), the crossing kernel is `[S=w]`. Its Boolean polynomial
in fixed additive coordinates is

```math
\prod_{i=1}^t(1+S_i+w_i),
```

which has total degree `t`. Bare kernels (3) have total degree at most
`2k`. Their affine substitutions retain this bound, so when `2k<t` their
span cannot contain the crossing. The controlled pieces in Section 4
are outside this elementary degree restriction.

There is a second useful way to see the alignment problem. On `L^*`, the
kernel `f(Sw)` has matched multiplicative exponents `S^n w^n`. Equality
uses opposite exponents `S^n w^(-n)`. Replacing `w` by its inverse would
align them, but inversion does not preserve the additive convolutions in
(5). Formula (6) therefore cannot be carried through that replacement.

As a diagnostic only, suppose a hypothetical synthesis used about
`d/binom(t,k)` pieces of cost `d 2^k`. Its exponent would be
`2+c-H_2(c)`, minimized at `c=1/3` to
`log2(8/3)=1.415037499...`. This is slightly larger than `sqrt(2)`.
No such synthesis is supplied, and the calculation is neither an upper
bound for an actual crossing decomposition nor a lower bound on better
syntheses. It indicates how small an additional improvement could matter
if an efficient synthesis were found.

## Reproduction

```sh
python3 tools/check_trace_bilinear_kernel.py --output /tmp/trace-bilinear-kernel.json
```

The standard-library checker verifies the degree-33 irreducibility and
integer certificate. For smaller fields it constructs the binary kernels,
computes ordinary ranks and both additive-convolution ranks exactly,
checks every top-weight Vandermonde coefficient, and compares the block
formula with direct four-leg cuts. It does not claim to construct or rank
the full `2^33`-label matrix explicitly. The proof above supplies the large
rank certificate. No Lean formalization is used.
