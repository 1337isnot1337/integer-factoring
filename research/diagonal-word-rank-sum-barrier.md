# A rank-sum barrier for all unsymmetrized diagonal words

## Theorem

Work over any field of characteristic two. Let `I` be the four-by-four
row identity for one binary crossing, and put `N=E_2+I`. For
`lambda>=1`, set `d=2^lambda` and consider the **entire** tensor-word
space

```math
T_f=\sum_{A\subseteq[\lambda]}c_A
       \bigotimes_{j=1}^{\lambda}
       \begin{cases}N,&j\in A,\\ I,&j\notin A.\end{cases}        \tag{1}
```

No symmetry between positions is assumed. Coefficients may vary
arbitrarily over the field. If `s(f)` counts the nonzero coefficients
`c_A`, then the two cut ranks satisfy

```math
\boxed{\quad r(T_f)c(T_f)\ge d^2s(f).\quad}           \tag{2}
```

Consequently every exact decomposition of the crossing into pieces
from this full word space satisfies

```math
\sum_i\sqrt{r(T_{f_i})c(T_{f_i})}\ge d^{3/2}.         \tag{3}
```

There is no restriction on the number of pieces. Since `3/2>sqrt(2)`,
such a decomposition never passes the strict finite rank-sum criterion.
With polynomially many pieces, (3) forces leading width coefficient
`rho>=3/4`. This is a scoped construction barrier, not a lower bound
for integer factoring or for arbitrary crossing tensors.

The space includes arbitrary weighted monomial and exterior-form
combinations. It is larger than a position-invariant coefficient span.
The matching of the same subset on the two sides is essential;
off-diagonal nonlinear selector interactions are outside the theorem.

## 1. The local row and column blocks

In characteristic two, `N` has rank one and satisfies `N^2=0`.
Its constant Jordan form consists of one two-dimensional nilpotent
block and two zero scalar blocks. Therefore the local row module is

```math
F[X]/(X^2)\ \oplus\ F\ \oplus\ F,
```

with `N` acting by multiplication by `X` in the first summand and
by zero on the others.

For the column cut, write `R=ww^T`, with
`w=e_0 tensor e_0+e_1 tensor e_1`. The column matrices of `I,N` are
respectively `R,E_2+R`. Here `R^2=0`, `E_2R=RE_2=R`, and
`J=E_2+R` satisfies `J^2=I` and `JR=R`. Multiplication by `J` turns
the local column pencil `aR+bJ` into `aR+bI`.

Both `R` and `N` are nonzero rank-one square-zero operators, so they
have the same constant Jordan form. Tensoring these changes of basis
shows that the column rank of (1) is the row rank of the
complement-coefficient polynomial

```math
f^*(X)=\sum_A c_A X_{[\lambda]\setminus A}.           \tag{4}
```

This is a rank statement over every characteristic-two field. It does
not rely on the coefficients being only zero or one.

## 2. A monomial gives a triangular minor

Let `J` be a subset of the variables, and let `f_J` denote the
restriction obtained by setting all variables outside `J` to zero.
On the square-zero algebra

```math
F[X_j:j\in J]/(X_j^2:j\in J),
```

suppose `c_A!=0` for some `A subset J`. Consider the multiplication
matrix with columns indexed by `X_B` and rows indexed by `X_(A union B)`,
where `B` ranges over subsets of `J minus A`.

The entry in the row for `B` and column for `C` is zero unless
`C subset B`. In the diagonal case it is exactly `c_A`. Ordering
subsets by cardinality gives a triangular square minor of size
`2^(|J|-|A|)` with nonzero diagonal. Thus

```math
rank(\text{multiplication by }f_J)\ge2^{|J|-|A|}.    \tag{5}
```

Choose `A` of smallest cardinality among supported subsets of `J`.
Define

```math
\delta(J)=\min\{|A|:c_A\ne0,\ A\subseteq J\},\qquad
u(J)=2^{-\delta(J)},                                \tag{6}
```

with `u(J)=0` if there is no such subset. The function `u` is increasing
on the Boolean lattice.

The active-variable block indexed by `J` occurs with multiplicity
`2^(lambda-|J|)` in the local direct-sum decomposition. Summing (5)
over every block gives

```math
r(T_f)\ge d\sum_Ju(J).                             \tag{7}
```

Apply the same argument to (4). Put

```math
\Delta(J)=\max\{|A|:c_A\ne0,\ J\subseteq A\},\qquad
v(J)=2^{\Delta(J)},                                \tag{8}
```

with value zero when the set is empty. The function `v` is decreasing.
Complementing the block index in (7) yields

```math
c(T_f)\ge\sum_Jv(J).                               \tag{9}
```

## 3. An elementary correlation inequality

For an increasing real function `u` and a decreasing real function `v`
on the uniform Boolean lattice,

```math
\left(\sum_Ju(J)\right)\left(\sum_Jv(J)\right)
\ge d\sum_Ju(J)v(J).                               \tag{10}
```

For completeness, prove this by induction on the number of bits.
Condition on the last bit and let `u_0,u_1` and `v_0,v_1` be the
conditional means. The induction hypothesis bounds the conditional
mean of the product by `u_b v_b`. Monotonicity gives
`u_0<=u_1` and `v_0>=v_1`, so

```math
\frac{u_0v_0+u_1v_1}{2}
-\frac{(u_0+u_1)(v_0+v_1)}4
=\frac{(u_1-u_0)(v_1-v_0)}4\le0.
```

This proves (10), including the zero-dimensional base case.

Whenever `c_A!=0`, the definitions give
`delta(A)<=|A|<=Delta(A)`, hence `u(A)v(A)>=1`. Combining
(7), (9), and (10) now gives

```math
r(T_f)c(T_f)
\ge d\left(\sum_Ju(J)\right)\left(\sum_Jv(J)\right)
\ge d^2\sum_Ju(J)v(J)
\ge d^2s(f),
```

which is (2).

## 4. From a piece to the whole crossing

The tensor words in (1) are linearly independent, since `I,N` are
linearly independent locally. The crossing is

```math
E_d=(I+N)^{\otimes\lambda}
    =\sum_{A\subseteq[\lambda]}N_A.                 \tag{11}
```

In any decomposition into this word space, every one of the `d`
word coefficients in (11) must appear in at least one piece. Therefore
`sum_i s(f_i)>=d`. Equation (2) implies

```math
\sum_i\sqrt{r_i c_i}
\ge d\sum_i\sqrt{s(f_i)}
\ge d\sqrt{\sum_i s(f_i)}
\ge d^{3/2}.
```

For `M` pieces, some rank product is at least `d^3/M^2`. If
`M` is polynomial in `lambda`, its logarithm is at least
`3lambda-O(log lambda)`, giving the stated `rho>=3/4`.

## 5. Why this is not a general kernel-rank inequality

In the additive-kernel representation, the space (1) consists of
matched monomials `sum_A c_A S_A N_A`. Its ordinary kernel-matrix rank
equals the number of nonzero `c_A`, because the two monomial-to-entry
changes of basis are invertible.

It would be false to replace this matched space by arbitrary kernels.
For `lambda=4`, the kernel `K(S,w)=S dot w` over `F_2` has ordinary
matrix rank four. Every nonzero row and column is a linear Boolean
form, whose convolution rank is two. Hence its two tensor-cut ranks
are both `30`, and

```math
30^2=900<16^2\cdot4=1024.
```

This kernel is outside (1). The checker retains it as a negative
control against enlarging the theorem's scope.

## Exact replay

```sh
python tools/check_diagonal_word_rank_sum.py \
  --output evidence/diagonal-word-rank-sum.json
```

The replay checks the local cut conversion, monomial minors and
correlation bounds on every binary word polynomial through four
positions, plus extension-field controls and the explicit general
kernel counterexample. The all-field, all-size statement is the
ordinary proof above; no new Lean formalization is asserted.

The construction, proof development, and replay used substantial OpenAI assistant support.
