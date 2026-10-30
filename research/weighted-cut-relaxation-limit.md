# An exact limitation of the scalar cut budgets

## Proposition

The weighted quotient budgets, ordinary rank totals, and the rank-one
and rank-two complement bounds do **not**, on their own, force the
finite rank-sum score above `d^sqrt(2)` in all dimensions. Their scalar
relaxation admits integer profiles with cost `O(d log d)`.

This is an abstract list of rank pairs. **No tensors realizing the
list are constructed or claimed to exist.** In particular, this note
does not improve a factoring algorithm or establish an admissible
crossing decomposition.

Let `d=2^m`, with `m>=1`. For every `0<=j<=m`, include both pairs

```math
(r_j,c_j)=(d2^j,d/2^j),\qquad(c_j,r_j).              \tag{1}
```

The two occurrences of `(d,d)` at `j=0` are separate pieces and are
both included. Every coordinate is a positive integer at most `d^2`.
Writing `R,C` for the two total ranks and `S` for the geometric-rank
cost gives the exact formulas

```math
R=C=2d^2+d-1,\qquad S=2d(m+1).                     \tag{2}
```

Thus the ordinary budgets `R,C>=d^2` hold.

## Every weighted quotient budget holds

Fix an integer `1<=s<d`, and let `t` be the largest power of two not
exceeding `s`. Both forward pairs corresponding to `t` and `2t`
occur in (1). Their combined contribution to the first weighted cut is

```math
\begin{aligned}
\min(dt,s^2d/t)+\min(2dt,s^2d/(2t))
 &=dt+\frac{ds^2}{2t}\\
 &\ge ds.
\end{aligned}                                                    \tag{3}
```

Indeed, `t<=s<2t`, so the first minimum chooses its first entry and
the second minimum chooses its second entry. Put `x=s/t` in `[1,2)`.
The combined contribution divided by `dt` is

```math
1+x^2/2\ge x.                                     \tag{4}
```

This follows from `1+x^2/2-x=((x-1)^2+1)/2>0`.

For `s=d`, the endpoint pair `(d^2,1)` alone contributes `d^2=sd`.
For `s=0`, the claim is immediate. Swapping the paired orientations
proves the other weighted cut. Therefore, for all `0<=s<=d`,

```math
\sum_i\min(r_i,s^2c_i)\ge sd,\qquad
\sum_i\min(s^2r_i,c_i)\ge sd.                       \tag{5}
```

These termwise-minimum inequalities imply every subset form of the
quotient budget. For any index subset `J`,

```math
\sum_{i\notin J}r_i+s^2\sum_{i\in J}c_i
\ge\sum_i\min(r_i,s^2c_i)\ge sd.                   \tag{6}
```

There is no need to check only selected subsets or to assume a
particular way of grouping the abstract pieces.

## Every forced small-cut complement bound holds

Suppose an index subset satisfies `sum_{i in J} c_i<=2`. Each
coordinate is positive, so such a subset has at most two members.
In fact it has at most one: there is exactly one piece with `c_i=1`,
and any other member has `c_i>=2`. This argument includes `d=2`, where
the duplicate center creates two distinct pieces with `c_i=2`.

The empty subset leaves row total `R>=d^2`. A nonempty eligible
subset removes at most `d^2` from `R`, so its complement has row total

```math
\sum_{i\notin J}r_i\ge d^2+d-1
                         \ge\frac{d(d-1)}2.        \tag{7}
```

This satisfies the [all-field rank-two complement theorem](two-product-complement.md). For the
only possible rank-one subset, the selected abstract pair is
`(d^2,1)`. Every product/Fitting lower bound
`d^2-h(h+1)/2`, for `0<=h<=d`, is at most `d^2`, whereas the remaining
row total is at least `d^2+d-1`. Thus even all those bounds are
satisfied. The opposite-cut statements follow by symmetry.

This is exactly the information that the scalar rank sums force.
For actual tensors, cancellation could make a large subset have
column rank one or two even when its sum of individual column ranks
exceeds two. The abstract profile provides no information about such
cancellation and makes no claim that those extra tensor constraints
can be met.

## Consequence and stop criterion

For powers of two, (2) gives

```math
S=2d(\log_2d+1)=O(d\log d)=o(d^{\sqrt2}).            \tag{8}
```

The previously proved harmonic budget supplies a lower bound of order
`d log d` for the same relaxation. Along powers of two, its asymptotic
order is `Theta(d log d)`, even after the rank-one and rank-two complement
constraints above are imposed. Improving constants or finite cutoffs
inside this relaxation cannot by itself prove a universal exclusion.

At the concrete dimension `d=4096`, the abstract cost is `106496`.
The checker proves the strict comparison using integers only:

```math
106496^{985}<4096^{1393},\qquad1393^2-2\cdot985^2=-1.
```

This gives `106496<4096^(1393/985)<4096^sqrt(2)`.

The next construction question is whether actual nested subspaces or
compatible operator pieces can achieve comparable ranks. A useful
new obstruction must use information beyond these individual cut
rank totals, such as compatible image spaces, cancellation, or ranks
of structured partial sums.

## Exact controls

```sh
python tools/check_weighted_cut_relaxation_limit.py \
  --output evidence/weighted-cut-relaxation-limit.json
```

The checker includes both center occurrences, every integer `s`, both
orientations, and every subset whose relevant summed cut rank is at
most two, for powers of two through `4096`. The proposition itself is
the general proof above; the finite calculation is an arithmetic
control and does not assert tensor realizability. The proof and checker
were developed with substantial OpenAI assistant support.
