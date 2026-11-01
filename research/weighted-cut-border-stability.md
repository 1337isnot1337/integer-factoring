# Weighted cut budgets survive formal one-parameter degenerations

## Exact scope

Let `F` be any field, `R=F[[t]]`, and `K=F((t))`. Let `E_d` denote
the crossing tensor, whose row matrix is the swap on `F^d tensor F^d`.
Suppose a tensor `T(t)` has entries in `R`, specializes to `E_d` at
`t=0`, and has a finite decomposition

```math
T(t)=\sum_i A_i(t),\qquad A_i(t)\in K^{d\times d\times d\times d}.
                                                               \tag{1}
```

The individual pieces may have poles. Write `r_i,c_i` for their two
cut ranks **over the Laurent-series field `K`**. Then, for every
integer `0<=s<=d`,

```math
\sum_i\min(r_i,s^2c_i)\ge sd,\qquad
\sum_i\min(s^2r_i,c_i)\ge sd.                       \tag{2}
```

In addition, `sum r_i>=d^2` and `sum c_i>=d^2`. Consequently the
harmonic lower bound from the exact-decomposition argument also holds:

```math
\sum_i\sqrt{r_i c_i}
\ge d(2H_{2d}-H_d-1).                              \tag{3}
```

In particular, this formal degeneration cannot pass the strict finite
rank-sum criterion in dimensions `2<=d<=41`.

The statement concerns a supplied one-parameter formal degeneration.
It does not silently identify every abstract Zariski-closure point
over every base field with such a family. It also does **not** extend
the rank-one/Fitting or two-product complement refinements, or the
later cutoff of fifty-two, without separate arguments.

## 1. A saturated quotient removes poles

First write `T(t)=A(t)+B(t)`, and put `k=c(B(t))`. A column-rank
factorization over `K` gives, in the row arrangement,

```math
B(t)=\sum_{j=1}^k P_j(t)\otimes Q_j(t).              \tag{4}
```

Fix a constant `s`-dimensional subspace `U` of `F^d`, and define the
`K`-subspace

```math
L_K=\sum_j P_j(t)(U\otimes_FK),\qquad
\ell=\dim_K L_K\le ks.                             \tag{5}
```

Its lattice `L_R=L_K intersect R^d` is saturated: if a nonzero element
of `R` times `v in R^d` belongs to `L_R`, then `v` belongs to `L_K`
and hence to `L_R`. Therefore `R^d/L_R` is torsion-free. Since `R` is
a principal ideal domain and this quotient is finitely generated,
the quotient is free of rank `d-ell`. A basis gives a surjective
integral map

```math
q(t):R^d\longrightarrow R^{d-\ell},\qquad\ker q=L_R.
                                                               \tag{6}
```

The reduction `q(0)` is still surjective. After extending scalars to
`K`, the map kills `L_K`. Thus on `U tensor F^d`,

```math
(q(t)\otimes I)A(t)=(q(t)\otimes I)T(t).             \tag{7}
```

The right side is integral even though the pieces on the left of
(1) need not be. Reduce it at `t=0`. The swap sends `U tensor F^d`
onto `F^d tensor U`, so its projected reduction has rank

```math
rank_F((q(0)\otimes I)E_d|_{U\otimes F^d})
 =(d-\ell)s.                                       \tag{8}
```

The rank over `K` of an integral matrix is at least the rank of its
reduction: a nonzero reduced minor lifts to a nonzero minor over `K`.
Restriction and projection cannot increase rank. Hence

```math
r(A(t))\ge s(d-\ell)\ge sd-s^2k,
```

or equivalently

```math
\boxed{\ r(A(t))+s^2c(B(t))\ge sd.\ }              \tag{9}
```

If `ks>=d`, the last lower bound may be nonpositive; the argument
and inequality remain valid. The zero-rank case `k=0` is included.

Saturation matters here. Multiplying a generator by a power of `t`
does not change its generic span but can make its unsaturated
reduction vanish. The quotient in (6) retains the generic dimension
and a surjective reduction, which is exactly what (8) needs.

## 2. All subsets and the harmonic consequence

For a subset `J`, apply (9) to `B=sum_{i in J} A_i` and its
complement. Rank subadditivity yields

```math
\sum_{i\notin J}r_i+s^2\sum_{i\in J}c_i\ge sd.     \tag{10}
```

Minimizing over `J` chooses the smaller summand for each piece and
gives the first inequality in (2). Exchange the two cuts for the
second. Since both cut matrices of `T(0)` are invertible, both cut
ranks of `T(t)` are `d^2`, which gives the ordinary rank budgets.

The positive harmonic interpolant in
[the weighted-cut theorem](weighted-cut-rank-budget.md) now applies
without change: its proof uses only (2), the ordinary totals, and
the individual rank interval `1<=r_i,c_i<=d^2` for nonzero pieces.
This proves (3). The existing exact rational checker proves that
this harmonic bound is above `d^sqrt(2)` for `2<=d<=41`.

Thus permitting poles and taking one formal leading coefficient
does not evade the weighted quotient obstruction. Whether the
stronger complement theorems are stable under the same operation
remains a separate question.

## 3. Finite cancellation controls

```sh
python tools/check_weighted_cut_border_stability.py \
  --output evidence/weighted-cut-border-stability.json
```

The checker uses exact Laurent polynomials over finite fields. Its
test pieces have genuine poles, and their sum is an integral
perturbation of the crossing. It constructs a moving saturated image
space and integral quotient, checks annihilation and pole cancellation
after projection, and verifies the exact reduced rank in (8).
It also checks the harmonic cutoff with rational powers. These
controls support the proof; they do not replace its all-field
saturated-lattice argument.
The argument and checks were developed with substantial OpenAI
assistant support.
