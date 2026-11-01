# Degree-sensitive auxiliary paths in Smith crossing factorizations

## Result and scope

The local Smith-profile criterion can retain the auxiliary transition
orders of its factor matrices. This replaces the independent channel
sums by two finite nonnegative matrices. Exact rational positive-vector
inequalities for these matrices give a new sufficient crossing bound.

All auxiliary matrices stay in position order. The opposite factor in
each cut is scalar, so no cross-half commutation hypothesis is needed.
Finite jets always give valid order lower bounds. At the complete
polynomial precision proved below, raw Schur elimination also certifies
which entries are exactly zero.

No improving crossing is asserted here. The result strengthens the
search interface in
[smith-profile-crossing-transfer.md](smith-profile-crossing-transfer.md)
and is compatible with the independently telescoping gauges in
[auxiliary-gauge-smith-transfer.md](auxiliary-gauge-smith-transfer.md).

## 1. Fixed local data and regular factors

Let `K(z;x,y,u,v)` be a fixed degree-`L` polynomial `D`-by-`D` matrix
tensor over a fixed finite field. External labels range over an alphabet
of size `q`. Assume the exact endpoint and boundary identities

```math
f^Te=1,\qquad K(1;x,y,u,v)e=[x=v][y=u]e.
```

The ordered scalar product has coefficients

```math
f^T K(z)_1\cdots K(z)_t e=\sum_{k=0}^{Lt}z^k A_{t,k},              \tag{1}
```

which sum to the crossing on `q^t` labels.

For either cut, place both auxiliary indices on one physical half.
After any chosen rational auxiliary gauge, use a local Laurent-series
factorization at infinity of the form

```math
K_\infty(w)=\sum_{i=1}^{b}w^{-d_i}U_i(w;\text{left})
                                 v_i(w;\text{right}),            \tag{2}
```

where `b=q^2` for a full Smith factorization, `U_i` is a `D`-by-`D`
matrix, and `v_i` is scalar. Every outside series is regular at zero.
For the other cut, use

```math
K_0(z)=\sum_{i=1}^{b}z^{\nu_i}U'_i(z;\text{left})
                              v'_i(z;\text{right}).              \tag{3}
```

The side carrying the auxiliary indices can be chosen independently in
the two cuts. The integers `d_i,nu_i` may have either sign if gauges
introduce poles. Restore every scalar Laurent shift used in calculating
these integers.

The telescoping gauge proof supplies nonnegative fixed boundary shifts
`H_infinity,H_zero`: after multiplying each transformed boundary vector
by its specified power, both vectors are regular. Their total powers
are `H_infinity` and `H_zero`. With no gauges these shifts are zero.
Let `F_infinity,E_infinity` and `F_zero,E_zero` contain all coordinates
on which the corresponding regular boundary vectors can be nonzero.
Using all `D` coordinates is always safe.

## 2. Certified transition costs

For each infinity channel `i` and auxiliary edge `a -> b`, supply a
nonnegative integer lower bound `delta_(i,a,b)` for

```math
\min_{\text{left}}\operatorname{ord}_w (U_i)_{a,b}
+\min_{\text{right}}\operatorname{ord}_w v_i.                     \tag{4}
```

A channel with an identically zero scalar factor may be discarded.
Assigning infinity to an edge requires an exact zero certificate for
one of its factors. Otherwise a finite lower bound suffices. Define
independent costs `delta'_(i,a,b)` from (3).

Fix a rational `s>=1`, and form the nonnegative rational matrices

```math
(W_\infty(s))_{a,b}=\sum_i s^{d_i-\delta_{i,a,b}},\qquad
(W_0(s))_{a,b}=\sum_i s^{-\nu_i-\delta'_{i,a,b}},                  \tag{5}
```

where an exactly zero edge contributes zero. The row index is the
current auxiliary state and the column index the next state. Trimming
states outside every path from the boundary's initial support to its
final support is valid when the absent edges have been certified zero.

### How finite jets supply the costs

Normalize a local rational matrix by clearing denominators and a common
Laurent power, then apply complete raw Schur pivots. If the current
pivot is `z^theta u(z)`, with `u(0)!=0`, one valid rank-one extraction is

```math
L=R[:,j]/z^\theta,\qquad
V=R[i,:]/(z^\theta u),\qquad
R\longmapsto R-z^\theta LV.                                    \tag{6}
```

Starting from precision `P`, both `L` and `V` are known modulo
`z^(P-theta)`. A nonzero jet gives its exact order; a zero jet gives the
safe lower bound `P-theta`. It does not by itself give infinity.
Taking the minimum over all physical-half rows gives the corresponding
auxiliary-entry lower bound in (4). A scalar denominator that is a unit
at zero does not change these orders.

The same rule applies after the infinity substitution and its common
Laurent normalization. The loss is controlled by the normalized pivot
order `theta`, not by an unrestored degree parameter `d_i`.

### A complete-precision exact-zero certificate

Suppose the normalized matrix is polynomial, has `b` columns, and its
entry degrees are at most `J`. At the `j`th raw Schur step, every
remaining entry is, up to sign, a bordered `j`-by-`j` minor of the
original matrix divided by the selected `(j-1)`-by-`(j-1)` pivot minor.
The denominator is nonzero. This is the ordinary Schur-complement
determinant identity, valid after arbitrary pivot row and column
permutations.

A nonzero numerator has degree at most `jJ`; the denominator is a
polynomial of nonnegative valuation. Thus every nonzero residual entry
has valuation at most `jJ<=bJ`. Consequently

```math
P=bJ+1                                                       \tag{7}
```

certifies every zero residual jet as an exact rational zero. A zero
`L` or `V` jet at precision `P-theta` means that its numerator residual
jet vanishes modulo `z^P`, so it too is exactly zero. At (7), the finite
jets determine all the orders needed in (4), including exact zeros.

This argument requires tracked raw Schur updates (6), or an equivalent
algorithm with all scaling factors accounted for. A generic truncated
matrix with unknown polynomial degree does not have this zero test.
For rational matrices first clear a supplied common denominator;
restore its valuation in the channel exponents. Its remaining unit
part changes no order or zero pattern.

## 3. Count paths within the full-word degree budget

For an infinity channel word `i_1,...,i_t` and auxiliary path
`a_0,...,a_t`, put

```math
D(i)=\sum_{j=1}^t d_{i_j},\qquad
\Delta(i,a)=\sum_{j=1}^t\delta_{i_j,a_{j-1},a_j}.
```

After boundary regularization, this path cannot contribute to
coefficient `w^(-k)` unless

```math
D(i)-\Delta(i,a)\ge k-H_\infty.                              \tag{8}
```

For an admitted path the two whole-word outside factors require at
most `D(i)+H_infinity-k+1` coefficient splits. This is at most
`C_infinity(t+1)` for a fixed constant, for example
`C_infinity=1+H_infinity+max(0,max_i d_i)`. Each split is one ordinary
outer product across the cut: the first factor is the ordered product
of the chosen auxiliary entries and the boundary entries, and the
second is a product of scalar opposite-half factors.

There is no independent degree split at every site. Introducing one
would create an unnecessary exponential factor.

Let `1_F,1_E` be the indicator vectors of the two boundary supports.
The weighted sum over all channel words and auxiliary paths is exactly

```math
Q_\infty(t)=1_{F_\infty}^T W_\infty(s)^t1_{E_\infty}.
```

For every path satisfying (8), its weight is at least
`s^(k-H_infinity)`. Therefore

```math
r(A_{t,k})\le C_\infty(t+1)s^{H_\infty-k}Q_\infty(t).          \tag{9}
```

Counting several paths for one channel word is harmless; it either
overcounts the original word factorization or explicitly expands its
auxiliary contraction into these path factors.

At zero, a path contributes only if
`sum_j(nu_(i_j)+delta'_(i_j,a_(j-1),a_j))<=k+H_zero`.
The analogous bound is

```math
c(A_{t,k})\le C_0(t+1)s^{H_0+k}Q_0(t),\qquad
Q_0(t)=1_{F_0}^T W_0(s)^t1_{E_0},                              \tag{10}
```

where `C_0=1+H_0+L+max(0,-min_i nu_i)` is sufficient.
Noncommuting auxiliary matrices were only expanded in their original
order. No change in the scalar polynomial (1) has occurred.

## 4. Exact positive-potential certificate

Supply strictly positive rational vectors `h_infinity,h_zero` and
positive rationals `beta_infinity,beta_zero` satisfying

```math
W_\infty h_\infty\le\beta_\infty h_\infty,\qquad
W_0h_0\le\beta_0h_0                                         \tag{11}
```

entrywise. These are finite rational inequalities. For either cut put

```math
C_W=\left(\max_{b\in E}{1\over h_b}\right)\sum_{a\in F}h_a.
```

Nonnegativity gives `Q(t)<=C_W beta^t`. Equations (9)--(10) imply

```math
r(A_{t,k})c(A_{t,k})\le C(t+1)^2 B_{\rm path}^{,t},\qquad
B_{\rm path}=\beta_\infty\beta_0,                            \tag{12}
```

and hence a total rank-sum bound
`sum_k sqrt(r(A_(t,k))c(A_(t,k)))<=C(t+1)^2 B_path^(t/2)`.
All boundary shifts and support constants are absorbed into fixed `C`.
At `lambda=t log2(q)`, the resulting width is

```math
4\rho\lambda+O(\log(\lambda+2)),\qquad
\rho={\log B_{\rm path}\over4\log q}.                          \tag{13}
```

A verified certificate with `B_path<q^(2sqrt(2))` would improve the
factoring exponent through the existing effective crossing transfer.
The rational inequality `B_path^5<q^14` is a stronger convenient test.

The coefficient entries still come directly from (1). Enumerating the
fixed channel and state alphabets, and computing the required `O(t)`
series coefficients, constructs the factorizations in `2^O(t)` actual
bit operations. Restricting the least `q^t>=2^lambda` word alphabet to
binary labels costs only a fixed rounding factor. The field arithmetic,
table construction, and recursion remain charged as in the existing
transfer.

## 5. Relation to gauges and earlier channel walks

For a fixed factorization, an integral diagonal auxiliary gauge shifts
an edge order by a difference of state heights. Its weighted matrix
changes by diagonal similarity, so its spectral radius is unchanged.
Recomputing Smith channels after a gauge can change the bound and remains
a separate search choice. The rational vectors in (11) are certificates
for the resulting nonnegative matrices, not unverified changes of the
original field tensor.

The path test retains compatibility of consecutive auxiliary states.
Independent channel minima discard that information. No equivalence
with the optimization over all rational auxiliary gauges is proved or
assumed here.

The earlier [finite channel-walk transfer](https://github.com/1337isnot1337/integer-factoring/blob/bb3bf763ed054f6fdaf882649c037ad810e16b97/research/twisted-channel-walk-transfer.md)
uses certified zero products of matrix-valued channels and rational
positive potentials. The present criterion combines that counting idea
with local Smith degree budgets and entrywise auxiliary series orders.
Its scalar opposite-half factors avoid the earlier cross-half
commutation requirement. The bordered-minor argument supplies a finite
exact-zero test for this specific series factorization.

## Verification

The checker computes Schur factors as exact rational functions over
`F_2`, independently of truncated-series elimination. It expands their
finite jets, checks every order and exact-zero claim against the rational
functions, and compares explicit small-word crossing ranks with the
path counts and rational-potential upper bounds. It includes scalar,
ordered auxiliary, and nilpotent-loop examples, together with negative
controls for insufficient precision and omitted boundary shifts.

```sh
python3 tools/check_degree_sensitive_smith_paths.py \
  --output evidence/degree-sensitive-smith-paths.json
```

Finite examples do not prove an improving all-length family. Equations
(8)--(13), together with a valid local certificate, are the transfer
theorem. This work used substantial OpenAI assistant support.
