# Polynomial-time values of the characteristic-two crossing pieces

## Result and scope

Every final piece in the [characteristic-two construction](past-half-pieces.md)
can be computed in polynomial time in its partition weight. A single
polynomial Pfaffian gives all `4h+1` pieces for four labels at once. The
algorithm enumerates neither the dual polynomial nor the intermediate
partition spaces.

This removes a computational restriction on future searches. Labels may
now lie at any fixed polynomial weight in the binary label width, including
weights much larger than `lambda^2`. Exact matrix elimination then computes
rank factorizations within the `2^O(lambda)` cost allowed by the recursive
evaluator. A uniform bound on the resulting ranks remains necessary.
The result does **not** establish a width coefficient below `1/sqrt(2)`
or a factoring exponent above `1/2`.

The original operators and crossing identity are Claude (Anthropic)'s
Appendix B, pinned to
[`anthropics/formal-math@16a5ce41304fcdaea4e0be7c132e1c361dc26da5`](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth).
Their exact definitions are in
[AppHalfPiecesStatement.lean](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/proof/statements/AppHalfPiecesStatement.lean).
This note proves a new characteristic-two evaluation formula, using the
stable quotient and normalization proved in [the modular note](past-half-pieces.md).
It is ordinary mathematics, with no new Lean formalization.

## 1. The explicit matrix

Let `x,v,u,y` be partitions of the same positive weight `h` into distinct
odd parts. The order of arguments below is chosen to distinguish the two
operator labels `x,v` from the output and input states `u,y`.
The crossing piece itself has the usual argument order `A_s(x,y,u,v)`.

Form a symmetric, zero-diagonal matrix `M(z)` over `F_2[z]`. Its vertices
are four blocks, in order:

```math
U=\{U_t:t\in u\},\qquad A=\{A_n:n\in x\},\qquad
B=\{B_m:m\in v\},\qquad Y=\{Y_t:t\in y\}.
```

Every entry within one block is zero. The six cross-block entries are

| Pair | Entry of `M(z)` |
| --- | --- |
| `U_t,A_n` | `z^(n OR t) (1+z)` |
| `A_n,Y_t` | `z^((n AND t)-1) (1+z)` |
| `A_n,B_m` | `[bits(n-1) subset bits(m-1)] z^n` |
| `U_t,B_m` | `[bits(m-1) subset bits(t-1)] (1+z)^m` |
| `B_m,Y_t` | `[m AND t = 1] z^t (1+z)^m` |
| `U_t,Y_r` | `[t=r] z^r` |

Here `AND` and `OR` are bitwise operations on nonnegative integers.
The subset condition means that every binary one in the first number
also occurs in the second. Since every part is odd, `n AND t>=1`,
so all the displayed exponents are nonnegative.

**Profile theorem.** The final crossing pieces satisfy

```math
\boxed{\operatorname{Pf}M(z)
       =z^{2h}\sum_{s=-2h}^{2h}A_s(x,y,u,v)z^s.}
\tag{1}
```

In characteristic two the Pfaffian is the sum, over all perfect matchings
of the matrix vertices, of the products of their entries. Its usual
matching signs equal one. The number of vertices is even because each
label's number of odd parts has parity `h`.

## 2. The normal-order matching identity

In the quotient

```math
\overline V_2=\mathbb F_2[t_1,t_3,t_5,\ldots]/(t_1^2,t_3^2,t_5^2,\ldots),
```

write `c_k` for multiplication by `t_k` and `d_k` for its derivative,
where `k` is odd. They satisfy

```math
c_k^2=d_k^2=0,\qquad
d_kc_l+c_ld_k=[k=l],
```

and different creators, or different derivatives, commute. These are the
creation and annihilation relations needed for the matching expansion.
The source's positive-mode factor `k` is one in `F_2` for retained modes.

For odd `j`, use the field convention

```math
\mathcal T_j(z)
=\sum_{k\ge1\atop k\ {\rm odd}}
 \binom{k-1}{j-1}c_k z^k
 +\sum_{k\ge1\atop k\ {\rm odd}}
 \binom{k+j-1}{j-1}d_k z^{-k}.
\tag{2}
```

The coefficients are reduced modulo two. The coefficient of `z^r` in
a normally ordered product of these fields is the source operator
`t^mu[r]`. Normal ordering puts all derivatives to the right, so they
act before the creators.

Represent the output coefficient by derivatives for the parts of `u`
on the left, and the input state by creators for the parts of `y` on
the right. Between them place the normally ordered first and second
operator fields, initially with independent formal mode variables.
Equality of the input and output weights forces their shifts to be
`r,-r`. Retain the first weight `z^r` and discard the second mode weight.
This is coefficientwise bookkeeping, not an analytic evaluation of an
infinite field at a number. The coefficient is a sum over perfect matchings:
an annihilator may pair with a creator to its right, with the pairing
coefficient supplied by (2). There are no contractions within either
already normally ordered operator block.

This matching rule follows directly by moving a derivative past each
creator using `d_k c_l=c_l d_k+[k=l]`. The unpaired derivative at the
far right kills the vacuum. Induction on the number of factors pairs
all the remaining vertices. In characteristic two no permutation sign
remains. The result is the Pfaffian of the matrix of allowed pairings.

Only six types of pairing remain: `U-A`, `A-Y`, `A-B`, `U-B`, `B-Y`,
and `U-Y`. We now fold the dual polynomials and the angle coefficients
into those entries.

## 3. Fold the dual and angle weights before expansion

For a homogeneous operator monomial `a=t^mu` of weight `j`, let

```math
F_a(z)=\sum_{r\in\mathbb Z}
 [t^u]a[r](t^v)[-r]t^y\,z^r.
```

The series is read as a formal Laurent series at infinity. Its exponents
are bounded above; each coefficient is a finite sum. This is convenient
because intermediate states can have arbitrarily large weight before
the source's signed finite difference is applied.

The definition of `Q_s` gives

```math
\sum_s[t^u]Q_s(a,t^v)t^y\,z^s
=z^{-j}(1+z)^{h+j}F_a(z).
\tag{3}
```

Indeed, put `r=s-h+i` and sum
`binom(h+j,i) z^(r+h-i)` over `i`. The signs disappear in
characteristic two. The source's support identity later makes the
result zero outside `-2h<=s<=2h`.

Set `w=(1+z)/z`. In the first operator block, replace the field attached
to a part `n` of `x` by

```math
\sum_{j\ge1\atop j\ {\rm odd}}\binom nj w^j\mathcal T_j(z).
\tag{4}
```

The sum ends at `j=n`. Formula (4) combines the original
`a_n=sum_j binom(n,j)t_j` with its weight factor `w^j`.
Terms with an even operator index vanish in the quotient. Products with
a repeated odd index vanish as well, as proved in the modular note.
Thus multilinearity of the normal-order matching formula includes the
whole dual polynomial without listing its monomials.

Finally, multiply the vertex `A_n` by `z^n`, `B_m` by `(1+z)^m`,
and `Y_t` by `z^t`. Scaling a vertex scales every matching by its factor.
The product of all these factors is

```math
z^{\sum_{n\in x}n+\sum_{t\in y}t}
 (1+z)^{\sum_{m\in v}m}=z^{2h}(1+z)^h.
```

The factor `(1+z)^h` is exactly the remaining common factor in (3).
It remains to compute the six resulting matrix blocks.

## 4. The six blocks

The elementary binary binomial rule is
`binom(a,b)=1 mod2` exactly when the bits of `b` form a subset of
the bits of `a`. It follows by factoring
`(1+T)^a=prod_(j:bit_j(a)=1)(1+T^(2^j))` over `F_2`.
The same factorization gives

```math
\sum_{j\ {\rm odd}\atop j\subset b}w^j
=w(1+w)^{b-1}\qquad(b\ {\rm odd}).
\tag{5}
```

### Output to first operator

For `U_t,A_n`, the unscaled creator coefficient in (4) is
`sum_j binom(n,j) binom(t-1,j-1) w^j`. Its allowed indices are
the odd subsets of `n AND t`. Apply (5), then use `1+w=1/z`
and the vertex and mode factors `z^(n+t)`:

```math
z^{n+t}w(1+w)^{(n\ {\rm AND}\ t)-1}
=z^{n\ {\rm OR}\ t}(1+z).
```

### First operator to input

For `A_n,Y_t`, the annihilator coefficient has allowed indices `j`
with the non-low bits of `j` contained in those of `n` and disjoint
from those of `t`. Thus the sum is
`w(1+w)^(n AND NOT t)`. The mode factor `z^(-t)` cancels the
`Y_t` scaling, leaving

```math
z^n w(1+w)^{n\ {\rm AND}\ ({\rm NOT}\ t)}
=z^{(n\ {\rm AND}\ t)-1}(1+z).
```

### Between the two operator blocks

After the `A_n` and `B_m` scalings, this block is

```math
(1+z)^{m+1}
\sum_{k\ {\rm odd}\atop m\subset k}
 z^{(n\ {\rm AND}\ k)-1-k}.
\tag{6}
```

If a bit belongs to `n` but not to `m`, toggling that bit of `k`
pairs the summands in (6). It preserves the condition `m subset k`
and the exponent, because the bit contributes to both `n AND k`
and `k`. Each Laurent coefficient has finitely many such terms,
so the cancellation is valid and (6) is zero.

Otherwise `n subset m`, and `n AND k=n` for every term. The elementary
formal power series identity

```math
\sum_{k\supset m}z^{-k}
=z^{-m}\prod_{j:\,\mathrm{bit}_j(m)=0}(1+z^{-2^j})
=\frac{z}{(1+z)^{m+1}}
```

then reduces (6) to `z^n`. All these `k` are odd because `m` is odd.
This proves the `A_n,B_m` entry without retaining an infinite series
in the final matrix.

### The three remaining blocks

The `U_t,B_m` creator coefficient is `binom(t-1,m-1)`, giving the
displayed subset condition and the `B_m` scaling `(1+z)^m`.
The `B_m,Y_t` annihilator coefficient is `binom(t+m-1,m-1)`.
For odd `m,t`, it is one precisely when `m AND t=1`; the vertex
scalings give `z^t(1+z)^m`. The unaltered input-output pairing is
`[t=r]`, multiplied by `z^r` from `Y_r`.

These calculations prove every entry of the matrix in Section 1.
The normal-order matching identity and (3) now prove (1).

## 5. Polynomial bit cost

A partition of `h` into `k` distinct odd parts has `k^2<=h`.
Consequently the matrix in (1) has at most `4sqrt(h)` vertices.
Give a vertex the weight of its part. Each displayed polynomial entry
has degree at most the sum of its two endpoint weights. For example,

```math
(n\ {\rm OR}\ t)+1
=n+t-(n\ {\rm AND}\ t)+1\le n+t.
```

Every matching therefore has degree at most the sum of all vertex
weights, namely `4h`. The same bound holds for every principal
Pfaffian minor. The entire output profile has at most `4h+1` bits.

One can compute it by the fraction-free Pfaffian elimination implemented
in [modular_pfaffian.py](../tools/modular_pfaffian.py). Here is why
its divisions are exact and its intermediate lengths remain polynomial.
After an even prefix of pivot vertices `I` has been eliminated, the
remaining entry at `i,j` is `Pf M[I union {i,j}]`. If the next pivot
pair is `a,b`, update each remaining pair `i,j` by

```math
\frac{
  \operatorname{Pf}(Iab)\operatorname{Pf}(Iij)
 +\operatorname{Pf}(Iai)\operatorname{Pf}(Ibj)
 +\operatorname{Pf}(Iaj)\operatorname{Pf}(Ibi)}
 {\operatorname{Pf}(I)}
=\operatorname{Pf}(Iabij).
\tag{7}
```

The empty-prefix denominator is one. Identity (7) is the Pfaffian
condensation identity in characteristic two. For a nonzero prefix
Pfaffian it follows by eliminating that prefix over `F_2(z)` and
expanding the Pfaffian of the remaining four-by-four alternating block.
Multiplying back gives (7) as a polynomial identity, proving exact
division. Simultaneous row and column swaps introduce no sign here.
If the remaining first row has no nonzero pivot, the corresponding
Schur-complement row is zero and the full Pfaffian is zero.

Thus all stored entries are principal Pfaffian minors of degree at
most `4h`. The numerator of an update has degree at most `8h`.
There are `O(h^(3/2))` updates, and schoolbook polynomial multiplication
and exact division cost `O(h^2)` bit operations each. Matrix construction
has smaller polynomial cost. In particular, `O(h^4)` is a uniform
bit-operation upper bound for the whole profile. Polynomial coefficients
are bits, so no hidden rational coefficient height is present.

Alternatively, `Pf(M)^2=det(M)` in characteristic two. A polynomial-time
determinant algorithm followed by extraction of its even-degree
coefficients gives the same profile. The implementation uses (7).
The variable `z` here indexes crossing pieces. It is separate from the
divisor variable in the coefficient algebra used for factoring. These
exact polynomial divisions compute prime-field constants before any
coefficient-algebra table is evaluated; they do not invert elements of
that algebra.

## 6. An effective actual-rank interface at polynomially larger weights

Fix a constant `C`. Suppose a uniform algorithm computes an integer
`h(lambda)<=C(lambda+1)^C` and injectively maps all `lambda`-bit strings
to distinct-odd-part partitions of that weight. Assume the weight and
each label can be computed in `2^O(lambda)` bit operations, with fixed
constants. The selection must be a uniform algorithm; an arbitrary
uncomputable family of subsets is not covered.

The profile theorem computes every entry of every final piece in
`poly(h)=poly(lambda)` time after its four labels have been obtained.
There are `D=2^lambda` labels, `D^4` tensor entries, and `4h+1`
pieces. Thus all the final matrices can be constructed in
`2^O(lambda)` bit operations, even if `h` is much larger than
`lambda^2`. No enumeration of the `exp(Theta(sqrt h))` ambient
partition space occurs.

Perform deterministic Gaussian elimination over `F_2` on each
`D^2` by `D^2` flattening. This also takes `2^O(lambda)` time.
For a matrix `A` of nonzero rank `r`, choose independent columns `J`
and rows `I` so that `B=A[I,J]` is invertible. Then

```math
A=A[:,J]\,B^{-1}A[I,:]
```

is an explicit rank-`r` factorization. Pad it with zero channels to
`2^ceil(log2 r)` channels. For a zero piece, use a zero factor with
one channel and width zero. All factor entries lie in `F_2`.
The widths, both factors, and their entries are therefore uniformly
computable within `2^O(lambda)` time. Their widths are at most
`2lambda`.

Consequently, if `1/2<=rho<1` and such a family has a proved uniform
actual-rank bound

```math
w(r_s^{\rm row})+w(r_s^{\rm col})
\le4\rho\lambda+c\log_2(2\lambda),
\qquad
w(r)=\begin{cases}0&r=0,\\\lceil\log_2r\rceil&r>0,
\end{cases}
\tag{8}
```

then it satisfies every effective modular-pieces hypothesis in
[Section 6 of the modular note](past-half-pieces.md#6-conditional-transfer-to-ring-arithmetic).
The polynomial piece count, unit normalization, extension to commutative
`F_2` coefficient algebras, and factor-computation costs are already paid.
The coefficient-algebra representation costs are covered by that transfer.
Its factoring exponent is
`log_2(1/rho)` by the complete transfer there.

This is a reduction of the remaining research problem to (8), not a
proof of (8) with a better `rho`. Gaussian elimination computes factors
but does not establish a uniform asymptotic rank saving. A search for a
label set whose own cost exceeds `2^O(lambda)` would also need a new
cost argument.

## 7. Exact checks and the tested high-weight labels

Run

```sh
python3 tools/check_modular_pfaffian.py
python3 tools/check_modular_pfaffian.py --full
```

The checker compares fraction-free elimination with independent
perfect-matching expansions on 80 small alternating polynomial matrices.
It verifies 992 exact product/division cases and 768 folded-block
identities against direct binomial sums. Removing a folded factor
provides a negative control.

The full suite compares the profile formula against the original
normal-order and intermediate-partition calculation at weights
8, 12, and 16. It checks 45,122 individual piece entries. The two
paths share no dual-monomial or intermediate-state enumeration in
the Pfaffian calculation.

The new evaluation algorithm also tests labels close to a filled set
of distinct odd modes. For a partition `epsilon` of an excitation
`E`, with at most `k` parts and padded with zeroes, set

```math
\mu_i=2(k-i)+1+2\epsilon_i\quad(1\le i\le k).
```

These are distinct odd parts with weight `h=k^2+2E`.
For `k>=E`, the label count is `p(E)`, while `k` can raise the ambient
weight far above the encoding's excitation scale. This is a concrete
test of the newly accessible high-weight regime.

The nine recorded cases do not improve the width coefficient. For example,
`(k,E)=(2,2),(4,2),(8,2),(16,2)` all have two labels, but their maximum
rank products are respectively `4,8,16,16`. The five-label cases
`(4,4)` and `(8,4)` both have maximum rank product 380.
These are exact characteristic-two ranks on the stated finite domains.
They give neither an asymptotic lower bound nor an improved factoring rate.

The full profiles, verification counts, and source hashes are recorded in
[modular-pfaffian.json](../evidence/modular-pfaffian.json).
