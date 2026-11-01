# Explicit Pascal zero spectra from binary linear maps

## Constructive extension

The [filtered exponent classification](filtered-pascal-exponents.md)
keeps the original odd zero spectrum. An arbitrary invertible binary
linear exponent map has another explicit zero spectrum, obtained from
minimum-weight row bases of its projected matrices. This includes
linear maps outside the filtered class.

The result supplies actual polynomial kernels and their complete zero
profiles. The infinity profiles must still be computed from those same
kernels. Increasing a zero order also changes the restored physical
height cost; it is not a free pressure improvement. No seed below
`1/sqrt(2)` is supplied here, and the factoring score remains `1/2`.

## 1. A weighted projected rank function

Let `M` be an invertible `r` by `r` matrix over `F_2`, with rows
numbered `0,...,r-1`, and put `L=2^r`. For a column set `I`, define

\[
 g(I)=\min\left\{\sum_{j\in J}2^j:
      M[J,I]\text{ is an invertible square matrix}\right\}.       \tag{1}
\]

Here `|J|=|I|`, and `g(empty)=0`. Restricting an invertible matrix
to any column set preserves column independence, so a basis exists.
The minimum is found by scanning rows from low index to high index
and retaining each independent row. Ordinary weighted-basis exchange
proves this greedy rule. Equivalently, if `r_j(I)` is the rank of
the first `j+1` rows restricted to `I`, then

\[
 g(I)=\sum_{j=0}^{r-1}2^j\bigl(r_j(I)-r_{j-1}(I)\bigr)
     =2^{r-1}|I|-\sum_{j=0}^{r-2}2^j r_j(I).             \tag{2}
\]

The convention is `r_(-1)=0`. The second expression follows by
telescoping and full column rank at the last row.

Each `r_j` is a submodular rank function on column sets. Indeed the
span of the union is the sum of the two spans, while the span of the
intersection is contained in their intersection. The dimension formula
gives the required inequality. Equation (2) therefore shows that `g`
is supermodular:

\[
 g(I)+g(J)\le g(I\cup J)+g(I\cap J).                    \tag{3}
\]

It is also increasing: a row basis after restriction to a larger
column set contains a row subset spanning the smaller projection,
whose positive weight is no greater. Finally,

\[
 g(\{0,\ldots,r-1\})=L-1.                              \tag{4}
\]

## 2. Exact orders of the Boolean subcube sums

Let `pi(i)` be the integer whose binary vector is `M` times the
binary vector of `i`. For each bit mask `q`, let `I` be its
complementary column set. Then

\[
 \boxed{\quad
 v_{T+1}\left(\sum_{i\supseteq q}T^{\pi(i)}\right)=g(I).
 \quad}                                                \tag{5}
\]

To prove (5), expand at `T=1`. The coefficient of degree `k` is
the sum over the affine subcube `i superset q` of
`binom(pi(i),k) mod2`. By the binary binomial identity this is the
product of the output-coordinate functions indexed by the set bits
of `k`.

On the free coordinates `I`, these output coordinates are affine
linear functions. If their linear parts have rank less than `|I|`,
their product is constant on fibers of even size. Its sum over the
subcube is zero. A nonzero coefficient therefore requires a spanning
row set. Its integer mask is at least the minimum basis weight (1),
because any spanning row set contains a basis of no greater weight.

At the minimum basis mask, the selected affine functions form a
bijection from the free subcube to `F_2^|I|`. The sum of their
product is one: only the all-one output contributes. This proves
the exact first nonzero coefficient in (5). The argument also
allows any fixed affine translation of the output bits.

## 3. The polynomial kernel and its complete zero spectrum

Use the same Pascal matrix and block kernel as in the preceding
constructor:

\[
 P_{ij}=[j\subseteq i],\quad c=1+u,\quad
 E_\infty(i)=2\pi(i)+1,
\]

\[
 A=\mathrm{diag}(c^{E_\infty(i)})P^{\mathsf T},\quad
 B=P^{\mathsf T}\mathrm{diag}(c^{2L-E_\infty(i)})P,
\]

\[
 \mathcal D_M=\begin{pmatrix}A&u^{2L}P\\I&B\end{pmatrix}.
                                                               \tag{6}
\]

All exponents are positive and smaller than `2L`. Hence the
polynomial inverse identities, determinant one, and the endpoint
`[[0,P],[I,0]]` hold exactly. The independent matching-type endpoint
argument from the weighted Pascal note applies.

Define new zero energies by

\[
 E_0(i)=2L-1-2g(\overline i).                           \tag{7}
\]

The bar denotes the complementary bit mask in `r` bits. These
are positive odd integers. Applying (5) to the complemented output
map gives the exact lower-right entry orders

\[
 v_u(B_{ij})=2g(\overline{i\mathbin{|}j}).              \tag{8}
\]

Set zero heights, in boundary order `(B,X,V,Y)`, to

\[
 h_0(i)=(0,E_0(i)-2L,0,E_0(i)).                         \tag{9}
\]

The internal `A` block is regular. With `I=bar i` and `J=bar j`,
the weighted order of an internal `B` entry is

\[
 2\bigl(L-1-g(I)-g(J)+g(I\cap J)\bigr)\ge0,            \tag{10}
\]

by (3), monotonicity, and (4). Thus both internal Gaussian factors
are unimodular over the zero valuation ring. The remaining cross
block is the weighted `[[0,I],[u^(2L)P,0]]`. Its full one-particle
Smith list consists of **two copies of every `E_0(i)`**. The full
Fock profile consists of all subset sums of this list.

If the row flag of `M` is a coordinate flag, (7) just permutes the
ordinary odd energy list. Other flags can give repeated energies and
a different mean. For example,

\[
 M=\begin{pmatrix}1&1\\0&1\end{pmatrix},\qquad
 \pi=(0,1,3,2),\qquad
 E_0=(1,5,5,7),\qquad E_\infty=(1,3,7,5).              \tag{11}
\]

This kernel is outside the preceding fixed-zero-profile class.

## 4. The physical charge and the other cut remain essential

For the candidate infinity heights
`(E_infty-L,L,E_infty-L,-L)`, each of the two matching-pair
physical energies is

\[
 G(i)=E_0(i)+E_\infty(i)-2L.                           \tag{12}
\]

The exact zero partition function at scale `tau` is
`prod_i (1+exp(-tau E_0(i)))^2`. The infinity factor uses the
actual `2L` by `2L` matrix (6), not the exponent list by itself.
Both the matrix and its full Clifford evaluator remain explicit for
independent physical heights.

Equation (3) with complementary sets gives
`E_0(i)+E_0(bar i)>=2L`. The resulting larger average zero energy
is accompanied by the charge (12). Omitting this change would
create a false gain. Likewise, `E_infty` is not asserted to be an
infinity Smith spectrum.

The [constructor](../tools/linear_pascal_matroid.py) computes all
rank weights by binary elimination. The
[checker](../tools/check_linear_pascal_matroid.py) compares (5) with
literal polynomial sums, verifies supermodularity, and reconstructs
the complete zero Clifford profiles and finite original Fock cuts.

```sh
python3 tools/check_linear_pascal_matroid.py \
  --output evidence/linear-pascal-matroid-energies.json
```

No new Lean formalization or independent external acceptance is claimed.
