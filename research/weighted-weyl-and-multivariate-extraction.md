# Weighted Weyl crossings and multivariate coefficient extraction

## Result and status

There is an explicit polynomial crossing family on weighted occupation
labels. It retains ordinary occupation carries, has short coefficient
formulas, and admits exact blockwise rank computation. For binary mode
weights its labels have low energy growth and its piece count is small
enough for the factoring transfer. The direct factorization bounds lose
two state indices, however, and the tested raw coefficients do not
improve the established half exponent.

A separate multivariate extraction lemma extends the finite-jet
interface. Its coefficient cost must be included. An explicit nonunit
example below has one ring channel in each cut but returns the full
crossing ranks after coefficient extraction. This separates a useful
finite algebraic identity from an unproved rank improvement.

## 1. A polynomial crossing from weighted occupations

Fix positive integer weights `w_1,...,w_r`, an energy `h`, and labels

```math
\mathcal L_h=\{\alpha\in\mathbb N^r:w\cdot\alpha=h\},
\qquad D(h)=|\mathcal L_h|.
```

Use coordinatewise binomial coefficients, with out-of-range lower
indices giving zero. Define a polynomial tensor

```math
P(z)(\alpha,y,u,\beta)=
[\alpha+u=\beta+y]\,[t^\alpha]
\prod_{j=1}^{r}
 (1+t_j)^{\beta_j}
 \bigl(1+(1-z)^{w_j}t_j\bigr)^{y_j}.                \tag{1}
```

Its degree in `z` is at most `w dot y=h`. At `z=1`, its coefficient
is `binom(beta,alpha)`. Since both labels have energy `h`, the
coordinatewise inequalities `alpha<=beta` hold only when they are
equal. The conservation constraint then forces `u=y`. Hence, over
the integers and after reduction to every field,

```math
P(1)=E_{D(h)},\qquad
P(z)=\sum_{i=0}^{h}z^i A_i,
\qquad \sum_{i=0}^{h}A_i=E_{D(h)}.                 \tag{2}
```

Formula (1) comes from normal ordering weighted Weyl operators. It
can also be taken as a direct binomial construction; (2) needs no
operator-theoretic assumption, factorial division, or characteristic
restriction.

### A row factorization

Expand (1) in powers of `1-z`:

```math
P(z)=\sum_{e=0}^{h}(1-z)^{h-e}B_e.                 \tag{3}
```

A summand of `B_e` is described by occupations `k,l,v,w` with

```math
\alpha=k+v,\quad y=l+v,\quad
\beta=k+w,\quad u=l+w,
\qquad |k|_w=|l|_w=e,
```

and coefficient `binom(k+w,k) binom(l+v,v)`. Factoring through
the pair `(k,l)` proves

```math
\operatorname{rank}_{\rm row}B_e\le D(e)^2,
\qquad
\operatorname{rank}_{\rm row}A_i
 \le\sum_{\substack{0\le e\le h-i\\
                    \binom{h-e}{i}\ne0\text{ in }F}}
       D(e)^2.                                    \tag{4}
```

For clarity, the left factor at `(alpha,y;k,l)` is
`[alpha-k=y-l>=0] binom(y,l)`, and the right factor at
`(k,l;u,beta)` is `[beta-k=u-l>=0] binom(beta,k)`.

There is no corresponding column bound `D(h-e)^2` for `B_e`.
In particular `B_h=E` has full column rank `D(h)^2`. Fixing `(v,w)`
still leaves equalities involving `k,l` across the other cut, so
discarding those constraints would be invalid.

### A column factorization at prime-power weights

Suppose the coefficient field has characteristic `p` and every
`w_j` is a power of `p`. Frobenius gives
`(1-z)^(w_j)=1-z^(w_j)`. Put `sigma=alpha+u=beta+y`. Expanding (1)
therefore gives

```math
A_i(\alpha,y,u,\beta)=
[\alpha+u=\beta+y]
\sum_{w\cdot\gamma=i}
 (-1)^{|\gamma|}\binom y\gamma
 \binom{\sigma-\gamma}{\alpha-\gamma}.             \tag{5}
```

Use channels `(gamma,eta)`, where `eta=sigma-gamma` has energy
`2h-i`. The two factors are

```math
\begin{aligned}
L_{\alpha,u;\gamma,\eta}
 &=[\alpha+u=\gamma+\eta]
   \binom\eta{\alpha-\gamma},\\
R_{\gamma,\eta;y,\beta}
 &=[y+\beta=\gamma+\eta]
   (-1)^{|\gamma|}\binom y\gamma.
\end{aligned}
```

Consequently

```math
\operatorname{rank}_{\rm col}A_i\le D(i)D(2h-i).    \tag{6}
```

The two-index costs in (4) and (6) are real. Low state growth alone
does not replace them by one state index.

## 2. Binary weights: effective construction and actual finite ranks

Take characteristic two, weights `1,2,4,...,2^R`, and `h=2^R`.
The labels are ordinary partitions of `h` into powers of two.
Elementary counting gives

```math
\log_2D(h)=\tfrac12 R^2+O(R\log(R+2)).              \tag{7}
```

For an upper bound, choose each occupation independently up to
`floor(h/w_j)`. The logarithm of the resulting box size is
`R(R+1)/2+O(1)`. For a lower bound, choose the positive-weight
occupations independently up to `floor(h/((R+1)w_j))`, then fill
the remaining energy with weight-one parts. These choices fit below
`h`, are distinct, and give the stated lower estimate.

Thus, if `lambda` is comparable to `log2 D(h)`, then
`h=2^O(sqrt(lambda))`. The number `h+1` of raw coefficients has
entropy `O(sqrt(lambda))`, and complete coefficient tables and
their exact factorizations can be constructed in `2^O(lambda)`
time. For arbitrary binary label lengths, one may choose the least
available energy with at least `2^lambda` labels and restrict all
four coordinates. These statements concern effectiveness, not an
improved rank bound.

For exact rank computation, conservation in (1) splits the row
matrix into blocks with fixed `alpha-y`, and the column matrix
into blocks with fixed `alpha+u`. Both block decompositions are
disjoint. The supplied C++ checker computes each coefficient matrix
over `F_2` inside these blocks and adds its exact block ranks.
It retains occupation carries; it does not replace them by bitwise
conservation.

The finite results are:

| Energy | Labels | Nonzero raw pieces | Rank sum, decimal | Half threshold, decimal |
| ---: | ---: | ---: | ---: | ---: |
| 16 | 36 | 9 | 1112.3566 | 158.8342 |
| 31 | 166 | 16 | 18701.1819 | 1379.4532 |
| 32 | 202 | 17 | 31631.8308 | 1820.7871 |

The decimals summarize exact integer rank profiles in the JSON files.
At energy 32 the coefficient span has dimension 17; the coefficients
at degrees `0,2,...,32` are independent. Optional binary block dumps
support collective searches using any linear combination of these
coefficients without constructing a dense four-dimensional tensor.

For energy 32, failure of the raw coefficient certificate is already
exact without numerical logarithms: `A_0` has ranks

```math
R_0=40804=202^2,\qquad C_0=1821>202.
```

Its single contribution exceeds `202^(3/2)`, which exceeds
`202^sqrt(2)`. This calculation excludes the raw decomposition.
The separate [common row filtration](weighted-weyl-row-filtration.md)
excludes arbitrary recombinations in the original coefficient span
at every power-of-two energy at least four. The
[protected-minor witness](weighted-weyl-protected-minor.md) also
excludes the full external-symmetry span at energy 32.

From the repository root, reproduce the exact finite profiles with:

```sh
c++ -O3 -std=c++17 tools/explore_weighted_binomial_crossing.cpp -o /tmp/weighted-weyl
/tmp/weighted-weyl 16
/tmp/weighted-weyl 31
/tmp/weighted-weyl 32
python3 tools/check_weighted_weyl_interfaces.py
```

The committed profiles are in `evidence/weighted-binomial-h16.json`,
`evidence/weighted-binomial-h31.json`, and
`evidence/weighted-binomial-h32.json`. Native profiling requires a
C++17 compiler; the interface checker uses only the Python standard library.

Reproduce the main finite checkpoint with

```sh
g++ -O3 -std=c++17 tools/explore_weighted_binomial_crossing.cpp -o /tmp/weighted-binomial-crossing
/tmp/weighted-binomial-crossing 32 > /tmp/weighted-binomial-h32.json
```

The implementation uses the compiler's `unsigned __int128` support
and accepts energies at most 127. It checks the crossing identity
at every conservation-supported entry and uses exact binary rank
elimination. The energy-32 run checks 3,115,724 entries across the
two cuts. Dense block storage in this checker is not a claim about
the working space of the asymptotic factoring algorithm.

The independent Python checker `check_weighted_weyl_interfaces.py`
verifies 5,018 entries of each supplied factorization and coefficient
change, 994 crossing entries, and 56 exact rank bounds over
characteristics two, three, and five. It also checks the nonunit and
divided-power examples below. Its output is
`weighted-weyl-interfaces.json`.

## 3. Multivariate finite-jet extraction

The [border extraction lemma](border-and-nonunit-seed-transfer.md)
has a direct multivariate form. Suppose supplied row/column
factorizations hold over

```math
F[z_1,\ldots,z_s]/(z_1^{H_1},\ldots,z_s^{H_s}),
```

and their tensor sum is `z^t(E+G)`, where every monomial of `G`
has positive total degree. Write `S=sum sqrt(r_i c_i)` for the
supplied channel dimensions. For every `k` with
`k t_j<H_j` for all `j`, type-grouping `k` tensor copies and taking
coefficient `z^(kt)` gives an exact field seed with rank sum at most

```math
\left(\prod_{j=1}^{s}(kt_j+1)\right)S^k.            \tag{8}
```

The coefficient of a matrix factorization is a sum over pairs of
nonnegative multiindices adding to `kt`; there are exactly the
displayed number of possible first indices. This proves the two rank
bounds and hence (8). For any fixed finite parameter set, the loss
is polynomial in `k`, but the necessary precision still has to be
supplied.

A polynomial multivariate border with the displayed componentwise
leading form can be reduced to the univariate theorem by substituting
`z_j=z`. The leading coefficient remains `E`, and every larger minor
that was a zero polynomial stays zero after substitution. More
generally a specified positive-weight leading order can be preserved
by suitable integer substitution weights. Thus a multivariate
parameterization can help discover certificates, but does not remove
their rank or precision obligations.

## 4. A one-channel nonunit identity and its full extraction cost

Let `K_j>=max_(alpha in L_h) alpha_j`, and consider

```math
R=F[t_1,\ldots,t_r,s_1,\ldots,s_r]/
       (t_j^{K_j+1},s_j^{K_j+1}:1\le j\le r).
```

The fully separable tensor

```math
t^\alpha s^y t^{K-\beta}s^{K-u}
   =t^K s^K E(\alpha,y,u,\beta)                    \tag{9}
```

has one ring channel in both cuts. To verify (9), a coordinate with
`alpha_j>beta_j` exceeds the truncation and vanishes. If no such
coordinate exists, equal positive weighted energy forces
`alpha=beta`. The same argument applies to `y,u`.

However, top-coefficient extraction over this ring has bilinear rank

```math
\dim_F R=\prod_j(K_j+1)^2.                          \tag{10}
```

In the monomial basis, the bilinear form taking the top coefficient
of a product pairs each monomial with its unique complement; its
matrix is a permutation matrix. Since the occupation box contains
all labels, (10) is at least `D(h)^2`. The extracted scalar tensor
indeed has both ranks `D(h)^2`. Furthermore `(t^K s^K)^2=0` for
`h>0`, so ordinary tensor-power amplification cannot amortize the
loss inside this fixed ring. This ring is not a principal chain ring;
the chain-ring budget lemma does not claim to cover it.

Divided-power parameters need a separate precision argument as well.
In characteristic `p`, every positive-augmentation element of a
divided-power polynomial algebra has `p`th power zero. For a basis
monomial `gamma_alpha`, the coefficient in `gamma_alpha^p` is the
product of multinomial integers `(p alpha_j)!/(alpha_j!)^p`.
If `alpha` is nonzero, at least one such integer is divisible by
`p` (its `p`-adic valuation is the base-`p` digit sum of `alpha_j`).
Frobenius then proves the assertion for every finite linear
combination. An exact scalar commutation identity over that algebra
does not automatically survive arbitrarily many ordinary tensor
copies with a nonzero common multiplier.

These examples leave open a genuinely extendible border construction,
or a new collective factorization of the Weyl coefficient tensors.
They prevent a coefficient extraction or a nonunit scalar from being
counted as free in either proposal.
