# Where fast matrix multiplication can affect the present recursion

The current evaluator does not obtain a new logarithmic saving exponent by
replacing its fixed terminal contraction with fast matrix multiplication.
Faster grouped-table construction also leaves the leading contraction
coefficient unchanged. This is an audit of the supplied scalar recursion,
not an impossibility result for algorithms that use matrix multiplication.

## The source of the leading coefficient

All source links below refer to `anthropics/formal-math` at revision
`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`.

The general-pieces evaluator performs three distinct operations:

1. It enumerates the crossing-piece patterns.
2. For each pattern, it chooses a row or column grouping and constructs one
   smaller scalar clique instance.
3. It evaluates those child instances and adds their scalar values.

These operations are explicit in
[`Rec/Level.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Rec/Level.lean)
and
[`Rec/Descend.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Rec/Descend.lean).
The recursion interface supplies a scalar evaluation routine for one child.
It does not supply a matrix of related outputs or shared bilinear factors for
different children.

The count in
[`Rec/Count.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Rec/Count.lean)
separates the child evaluation bound `T` from the construction factors.
Using the source's notation, these latter factors include

\[
 s=2^{2N/K},\qquad u=2^{c\lambda},
 \qquad \Xi^{A+G+C_4},
\]

while the leading part of the child bound is

\[
 T\le 2^{e(\rho+3\zeta)N}
      \Xi^{A+K^2+(c+5)C_0K^2}.
\]

Here `N` is the parent volume. The factor multiplying the current exponent
`e` is obtained from the **volume of the chosen child**, rather than the
number of scalar multiplication operations in an otherwise fixed matrix
product. This distinction matters when proposing a per-level saving.

The grouped-table construction in
[`Group/Lemma.lean`](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/rates/GeneralPieces/Group/Lemma.lean)
enumerates groupings and writes their tables. Its table entries include sums
with at most `2^(2N/K)` terms. Sharing such sums through matrix multiplication
can reduce a construction exponent. In the existing growing-group schedule,
that is a summable additive error of order `1/K`; it is not the `rho e` term.

## A precise recurrence observation

Let `0 < rho < 1`. Suppose the exponent coefficients of a proposed revision
still satisfy the same scalar accounting recurrence

\[
 e_{j+1}=(\rho+\epsilon_j)e_j+\delta_j,
 \qquad e_0>0,\quad \epsilon_j,\delta_j\ge0,
\tag{1}
\]

where

\[
 \sum_{j\ge0}\epsilon_j<\infty,
 \qquad
 \sum_{j\ge0}\frac{\delta_j}{\rho^{j+1}}<\infty.
\tag{2}
\]

Then

\[
 e_j=\Theta(\rho^j).
\tag{3}
\]

Indeed, set `a_j=e_j/rho^j`. The recurrence becomes

\[
 a_{j+1}=(1+\epsilon_j/\rho)a_j
          +\delta_j/\rho^{j+1}.
\]

Thus `a_j>=e_0`. Its upper bound follows by unrolling the recurrence and using

\[
 \prod_{i\ge0}(1+\epsilon_i/\rho)
 \le\exp\left(\rho^{-1}\sum_i\epsilon_i\right)<\infty,
\]

together with the second series in (2). This proves (3).

For the square-root vertex schedule, `log K_j=Theta(2^j)`. Equation (3)
therefore gives the same logarithmic saving power

\[
 \alpha=\log_2(1/\rho).
\]

Changing a positive fixed terminal coefficient changes the implicit constant
in this analysis. Decreasing the nonnegative summable construction errors
also changes constants. Neither change alters this power. Equation (3)
describes the certified recurrence, rather than a runtime lower bound for
every implementation of the underlying mathematical function.

The existing terminal clique has a fixed number of vertices but an unbounded
label width. A faster fixed-clique solver can be a useful practical change.
It still supplies a fixed positive coefficient in the same recursion. Its
matrix-multiplication saving cannot be multiplied once again at every higher
level unless a corresponding new higher-level operation is specified.

## Why merely batching the list of scalar children is insufficient

The current piece-pattern and grouping counts are bounded by powers of
`Xi`. Their logarithms contribute to the polynomial-overhead part of the
recurrence. If a batch improvement only replaces a factor `a` in that count
by `a^beta`, with fixed `beta<1`, it removes
`(1-beta) log(a)` from this overhead. It does not change the width of any
child or the `e rho N` contribution of evaluating a child.

This observation does not preclude a joint evaluation algorithm whose saving
depends exponentially on the child's width. It specifies what such an
algorithm must establish: shared structure in the child inputs and a
multi-output contraction theorem, with its output count and table
construction costs included. The scalar interface does not provide that
structure automatically.

Likewise, a length-`R` inner product with independent inputs has bilinear rank
`R`: its coefficient matrix is the `R`-by-`R` identity, and each bilinear
multiplication contributes a rank-one matrix. A faster matrix multiplication
uses sharing among many inner products. It cannot replace an isolated inner
product by the same saving without identifying those shared inputs and
accounting for all requested outputs.

## Concrete ways a revised theorem could escape this audit

A valid improvement can change the child volume, replace scalar child calls
by a structured multi-output primitive, or use the complete row/column rank
tradeoff with a different recursion. The last possibility already has a
precise analysis for the present Gaussian pieces in
[`rectangular-grid-frontier.md`](rectangular-grid-frontier.md). Its
quarter-circle profile retains the one-half frontier even for unequal and
adaptive grid splits. That result concerns the declared profile model; a new
crossing family or a new contraction operation is outside its scope.

For a multi-output proposal, the useful next object is an explicit map from
one level's whole family of grouped tables to a collection of matrix products
or another fast algebraic primitive. The map must preserve exact scalar
recovery, expose the dimensions of every input and output, and improve the
coefficient of `eN` after construction and output costs are charged. No such
map is asserted in this checkpoint.

No Lean source changes or builds accompany this audit. The factoring bound
remains `exp(O(n/sqrt(log n)))`.
