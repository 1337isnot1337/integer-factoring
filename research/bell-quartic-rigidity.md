# Self-reciprocal quartics require a particular Hadamard non-rigidity tradeoff

## Result and scope

There is a direct algebraic connection between self-reciprocal quartic
[Bell crossing paths](bell-diagonal-crossings.md) and matrix rigidity.
For a path with a rank-one constant row coefficient, any rational
order/degree pressure `Q` supplies an explicit sparse change of the
`d`-by-`d` Walsh-Hadamard matrix satisfying

```math
\boxed{\quad
\mathcal R_{H_d}\!\left(\left\lfloor\frac{2Q}{M}\right\rfloor\right)
 \le\lfloor\sqrt Q\rfloor,
\qquad M=(d^2-d)/2-1>0.
\quad}                                                    \tag{1}
```

Here `R_H(r)` is the fewest changed entries needed to reduce the matrix
rank to at most `r`. This is a necessary condition, not a converse.
It does not exclude quartics or improve the current factoring exponent.
It replaces an unstructured parameter hunt with a concrete sparse-matrix
construction question and its additional symmetry obligations.

The theorem holds over every field of characteristic different from two.
It uses ordinary algebra with substantial OpenAI Codex assistance.
Independent review and Lean formalization are not claimed.

## 1. A chirped reshape turns Bell duality into rank-preserving matrix algebra

Write `d=2^k`, `G=F_2^k`, and let

```math
H_{a,b}=(-1)^{a\cdot b},\qquad H^2=dI.
```

For an eigenvalue vector `f` with `d^2` coordinates define

```math
\mathcal C(f)_{a,b}=(-1)^{a\cdot b}f_{a,b}.
```

This reshape preserves the number of nonzero entries. The explicit Bell
transform gives the identity

```math
\mathcal C(Tf)=d^{-1}H\mathcal C(f)^T H.
\tag{2}
```

Expand `(-1)^((a+c) dot (b+e))` to verify (2) entry by entry.
In particular, this operation preserves matrix rank. A translation of
the vector index is, after this reshape, a row permutation, column
permutation and two invertible diagonal sign matrices. Thus the same
rank-preservation statement holds for `S=T P_g`, where `P_g` translates
the index by `g`. Such transforms are regular involutions as well.
This includes the translated quadratic coordinates in the folded and
projected constructions.

Notice that `C(1)=H`: the crossing endpoint becomes a full-rank
Hadamard matrix, rather than an unstructured all-ones matrix.

## 2. The midpoint coefficient is a sparse modification

Let

```math
f(z)=a_0+a_1z+a_2z^2+a_3z^3+a_4z^4,
\qquad Sf(z)=z^4 f(1/z),\qquad f(1)=\boldsymbol1.
```

Then `a_4=Sa_0`, `a_3=Sa_1`, and `a_2=Sa_2`. Put
`B=C(a_0+a_1)` and `C=C(a_2)`. The endpoint reads

```math
H-C=B+\mathcal C\bigl(S(a_0+a_1)\bigr).
```

By (2),

```math
\operatorname{rank}(H-C)\le2\operatorname{rank}B\le2R_1,
\qquad \#\operatorname{supp}C\le R_2,                    \tag{3}
```

where `R_j` counts eigenvalue polynomials whose order at zero is at
most `j`. The support of `a_0+a_1` is contained in that counted set;
matrix rank is at most the number of nonzero entries. Changing exactly
the nonzero entries of `C` converts `H` to the matrix `H-C`.

More generally the same identity applies to every self-reciprocal even
degree: its midpoint is the sparse modification and the sum of the
lower coefficients supplies the two rank terms. The quartic case below
gives a particularly useful quantitative restriction.

## 3. Turning a pressure witness into the rigidity parameters

Assume `a_0` has exactly one nonzero coordinate. Its image `a_4=Sa_0`
has entries `+alpha,-alpha`. As in the
[cubic endpoint argument](bell-cubic-outer-rank-floor.md), at least
`M=(d^2-d)/2-1` other coordinates must have a nonzero middle coefficient.
Their orders are at most three.

For a scale `x>=1`, self-reciprocity makes the two profile sums reversals.
Write

```math
F(x)=\sum_j x^{2-\ell_j},\qquad Q(x)=F(x)^2,
```

where `ell_j` are the row orders. Three elementary estimates are

```math
R_1\le F/x,\qquad R_2\le F,\qquad F\ge M/x.
```

The first two count only the terms of order at most one and two,
respectively. The third uses the `M` forced middle coordinates.
Consequently `x>=M/F`, and

```math
R_1\le F^2/M=Q/M,\qquad R_2\le\sqrt Q.
```

Combining with (3), and retaining integer ranks and support counts,
proves (1). There is no norm inequality or positivity assumption on
field elements, so the argument applies in fixed finite characteristic.

For a growing-alphabet family at the critical pressure
`Q<d^(2sqrt(2))`, the necessary asymptotic parameters are

```math
\operatorname{rank}(H-C)
 <(4+o(1))d^{2\sqrt2-2},\qquad
\#\operatorname{supp}C<d^{\sqrt2}.
\tag{4}
```

The rank exponent is approximately `0.828427` and the changed-entry
exponent approximately `1.414214`. At small finite alphabets the rank
bound in (1) can exceed `d` and be vacuous. Also, tensor-powering a fixed
quartic seed increases its degree; (4) is not a necessary growing-matrix
condition for every possible fixed-seed factoring improvement.

## 4. What existing non-rigidity results do and do not supply

[Alman and Williams, Theorem 1.1](https://arxiv.org/abs/1611.05558)
give an all-field Walsh-Hadamard non-rigidity tradeoff with changed-entry
exponent `1+epsilon` and rank exponent `1-Theta(epsilon^2/log(1/epsilon))`
for sufficiently small `epsilon`. This motivates examining sparse
perturbations, but the theorem as stated does not specify the numerical
tradeoff (4). No required matrix construction follows merely by citing it.

Even an appropriate low-rank sparse modification would not yet be a
quartic crossing. The midpoint must satisfy the fixed-transform equation,
and the remaining matrix must split into a **coordinate-sparse** half
and its transform, with a nonzero rank-one outer coefficient. Matrix
rank alone does not give that sparse half. Actual pressure, coefficient
construction and the full bit cost must still be checked.

There is an explicit stronger converse at the decomposition level.
Choose a free eigenvalue vector `b` and a nonzero point vector `a_0`, then set

```math
a_1=b-a_0,\quad a_2=\boldsymbol1-b-Sb,
\quad a_3=Sa_1,\quad a_4=Sa_0.
```

The midpoint is automatically fixed by `S`, and the endpoint is exact.
The constructive matrix question is whether both `C(b)` and
`C(1-b-Sb)` are sparse with a useful complete order distribution.
This uses one free vector rather than five independent coefficient
vectors. Its pressure must still be checked; matrix rank alone does
not give the needed coordinate sparsity.

The next constructive target is therefore a sparse self-dual modification
with this stronger decomposition, or a non-self-reciprocal/other-degree
family outside the reduction. The strongest proved factoring score stays
`1/2`.

## 5. Exact replay

```sh
python3 tools/check_bell_quartic_rigidity.py --output evidence/bell-quartic-rigidity.json
```

The [saved controls](../evidence/bell-quartic-rigidity.json) construct 60
literal self-reciprocal quartics in characteristics `3,5,7` at physical
dimensions `2,4,8,16`. They verify the chirped transform identity, its
rank preservation under translations, the endpoint matrix decomposition,
and 240 exact rational pressure-to-rigidity comparisons. These controls
support the displayed algebra; they are not a global quartic search.
