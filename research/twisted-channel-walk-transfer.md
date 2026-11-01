# Scalar twists and finite channel walks for matrix-valued seeds

## Constructive goal and scope

The [commuting-halves transfer](commuting-halves-matrix-seed-transfer.md)
turns a fixed matrix-valued crossing
identity into ordinary field seeds after ordered tensor powers. This
note extends that construction in two explicit ways. Scalar commutation
relations can replace exact cross-half commutation. Certified zero
products within either half can then reduce the number of surviving
factor channels, with a finite weighted adjacency matrix recording the
saving.

These are effective sufficient conditions for a better finite seed.
They do not exhibit an improving seed or prove a factoring score above
`1/2`. All matrix sizes, alphabets, and fields below are fixed parts of
one finite certificate.

## 1. Scalar cross-half commutation is enough

Let `A_1,...,A_M` be tensors with entries in `Mat_b(F)` and

```math
\sum_i A_i=E_d I_b.
```

Suppose the row factorization of `A_i` has channels `a=1,...,r_i`:

```math
A_i(x,y,u,v)=\sum_a L_{i,a}(x,y)R_{i,a}(u,v).
```

Write `s=(i,a)` for a row channel. Require scalars
`chi_(s,t) in F^*` such that, for every pair of external inputs,

```math
R_s(u,v)L_t(x,y)=\chi_{s,t}L_t(x,y)R_s(u,v).       \tag{1}
```

The scalar may depend on the two channels, but not on their external
inputs. Impose an independent condition on the column factorizations,
which have `c_i` channels and may use different scalars. No relation is
required between two left factors or between two right factors.

For an ordered channel word `s_1,...,s_k`, repeated use of (1) gives

```math
L_{s_1}R_{s_1}\cdots L_{s_k}R_{s_k}
=\left(\prod_{j<l}\chi_{s_j,s_l}\right)
 (L_{s_1}\cdots L_{s_k})(R_{s_1}\cdots R_{s_k}).    \tag{2}
```

Each right factor crosses each later left factor exactly once. The
order within each half is unchanged. The prefactor is a scalar
depending only on the channel word, so it can be absorbed into the
left factor. In particular, scalar phases can vary among words of
the same piece type without changing the channel bound.

For a piece word `i_1,...,i_k`, (2) gives at most `prod_j r_(i_j)`
row channels and `prod_j c_(i_j)` column channels. Summing words of
one type, then applying `ell(T)=T_11` once, gives exactly the previous
bound

```math
S_k\le b\left(\sum_i\sqrt{r_i c_i}\right)^k.       \tag{3}
```

Finite homogeneous gradings with scalar commutation factors are
examples. Equation (1), rather than a grading terminology, is the
actual finite hypothesis. If the scalar depends on external inputs,
the present argument does not apply: the gathered phase may couple
the two sides of the factorization.

## 2. Omit channels whose ordered products vanish

Use the row alphabet

```math
\mathcal R=\{(i,a):1\le i\le M,\ 1\le a\le r_i\}.
```

Supply a directed graph `G_R` on this alphabet. Every omitted edge
`s -> t` must have one of the following exactly verified certificates:

```math
L_s(x,y)L_t(x',y')=0\quad\hbox{for all }x,y,x',y',
```

or

```math
R_s(u,v)R_t(u',v')=0\quad\hbox{for all }u,v,u',v'. \tag{4}
```

Additional edges may be retained even when their products vanish.
Thus the graph need only give a sound upper bound on possible channel
words. Define a column graph `G_C` with the corresponding independent
factor tables. All conditions in (4) are finite matrix identities.

If a channel word contains an omitted adjacent edge, one of the two
ordered half-products in (2) is zero. It may be removed from the
factorization for every external input. Let `R_nu(k)` count row
walks with `k` vertices whose piece labels have type `nu`, and let
`C_nu(k)` be the corresponding column count. For the matrix-valued
type piece `B_nu`, the supplied factor dimensions are therefore at
most `R_nu(k)` and `C_nu(k)`. The scalar pieces `ell(B_nu)` satisfy

```math
S_k\le b\sum_{|\nu|=k}\sqrt{R_\nu(k)C_\nu(k)}.    \tag{5}
```

This statement uses adjacent zero products only. Longer-word zero
identities can also be incorporated by a finite automaton that
remembers a fixed number of previous channel states, provided every
omitted word is certified to vanish. The simpler directed graphs
already give the transfer below.

## 3. A weighted spectral certificate

Choose positive real numbers `z_1,...,z_M`. If `i(s)` is the piece
label of a row channel, define

```math
(T_R(z))_{s,t}=[s\to t\text{ in }G_R]z_{i(t)},
\qquad (w_R(z))_s=z_{i(s)}.
```

Then the weighted walk count is

```math
P_R(k,z)=w_R(z)^T T_R(z)^{k-1}\mathbf 1
        =\sum_{|\nu|=k}R_\nu(k)z^\nu.              \tag{6}
```

Define `T_C`, `w_C`, and `P_C` using the column alphabet. Cauchy--Schwarz
applied to (5), with reciprocal piece weights in the second factor,
gives the exact finite bound

```math
S_k\le b\sqrt{P_R(k,z)P_C(k,z^{-1})}.              \tag{7}
```

For a directly checkable certificate, supply positive vectors
`h_R,h_C` and positive numbers `beta_R,beta_C` such that

```math
T_R(z)h_R\le\beta_R h_R,\qquad
T_C(z^{-1})h_C\le\beta_C h_C                     \tag{8}
```

entrywise. Set

```math
C_R=\frac{w_R(z)^T h_R}{\min_s(h_R)_s},\qquad
C_C=\frac{w_C(z^{-1})^T h_C}{\min_s(h_C)_s}.
```

Nonnegativity and induction give

```math
P_R(k,z)\le C_R\beta_R^{k-1},\qquad
P_C(k,z^{-1})\le C_C\beta_C^{k-1}.
```

Consequently the asymptotic rank-sum base is at most

```math
L=\sqrt{\beta_R\beta_C},\qquad
S_k\le b\sqrt{C_R C_C}\,L^{k-1}.                  \tag{9}
```

All data in (8) can be positive rationals, so validation requires
only exact rational arithmetic alongside the finite-field matrix
checks. An improving certificate needs
`L<d^(sqrt(2))`. Its resulting quality is
`rho=log(L)/(2log(d))`, up to an arbitrarily small slack if the
displayed upper bound is approached by certificates.

Equivalently one can study the sufficient asymptotic quantity

```math
\inf_{z_i>0}\sqrt{
 \operatorname{spr}(T_R(z))
 \operatorname{spr}(T_C(z^{-1}))}.                 \tag{10}
```

Here `spr` denotes spectral radius. For every strict upper bound on
a nonnegative matrix's spectral radius, a positive vector satisfying
the corresponding inequality exists: for `beta>spr(T)`, the vector
`sum_(j>=0)(T/beta)^j 1` is positive and has `Th<beta h`.
Strictness permits positive rational approximation when the matrix
entries are rational. Thus (8) is the finite certificate interface;
(10) is an optimization description, not a substitute for a
verified bound.

The infimum in (10) need not itself have a positive supervector
certificate with a constant prefactor. A reducible matrix can have
polynomial growth multiplying its spectral-radius power. For example,
`[[1,1],[0,1]]` has spectral radius one and a linear off-diagonal
entry in its powers, but no strictly positive vector satisfies
`Th<=h`. Bounds from strict rational upper certificates retain the
constant-prefactor assertion (9). These distinctions matter when
comparing a limiting base with a separate rank-sum lower bound.

When both graphs are complete, one may take constant vectors in
(8), with

```math
\beta_R=\sum_i r_i z_i,\qquad
\beta_C=\sum_i c_i/z_i.
```

The infimum of their product is
`(sum_i sqrt(r_i c_i))^2`, attained at
`z_i=sqrt(c_i/r_i)` up to a common factor. The channel-walk bound
therefore includes the original rank-sum theorem exactly.

## 4. Uniform construction and what remains to find

The actual matrix-valued type entries are still coefficients of the
ordered product of `sum_i A_i X_i`, in central variables `X_i`.
There are at most `(k+1)^M` type pieces and their entries have a
polynomial-size dynamic program for fixed certificate data. One
does not need to approximate eigenvalues or enumerate walks to
evaluate a type entry.

If explicit factor tables are wanted, enumerate the allowed channel
walks and use (2). The alphabets are fixed, so this costs `2^O(k)`.
The same bound covers complete scalar tables, exact field ranks,
and factorizations after restricting `d^k` labels to any desired
binary width. For `lambda` input bits, choosing
`k=ceil(lambda/log_2 d)` gives total rank width

```math
4\rho\lambda+O(1),\qquad
\rho=\frac{\log L}{2\log d},                      \tag{11}
```

with polynomial piece count and a `2^O(lambda)` constructor.

An actual advance must supply both the crossing identity and the
factor tables with their commutation phases and omitted-edge
certificates. An arbitrary sparse graph with a small spectral radius
does not suffice. The reduction identifies a concrete possible use
of matrix coefficients: nilpotent ordered products can remove
channels without a factor lost at every tensor copy. No such
improving crossing certificate is currently provided here.

## 5. Exact examples and checks

There is a small actual crossing certificate on which the graph
strictly reduces the supplied channel bound. Over `F_2`, take

```math
N_1=\begin{pmatrix}0&1\\0&0\end{pmatrix},\quad
N_2=\begin{pmatrix}0&0\\1&0\end{pmatrix},\quad
N_3=\begin{pmatrix}1&1\\1&1\end{pmatrix}.
```

They satisfy `N_1+N_2+N_3=I_2`, each `N_i^2=0`, and every product
`N_iN_j` with `i!=j` is nonzero. Let the external label size be two
and set `A_i=E_2 N_i`. Give each cut its trivial four-channel
factorization, putting the matrix `N_i` in the left factor and
scalar identity matrices in the right factor. The global
cross-half commutation requirement holds, although the left
coefficient matrices do not commute internally.

In each cut, all transitions between two channels of the same
piece may be removed. All transitions between different pieces
are retained. There are exactly

```math
3\,2^{k-1}4^k=12\,8^{k-1}
```

allowed channel walks of length `k`. Constant weights and constant
positive supervectors certify `beta_R=beta_C=8`, whereas the
unrefined supplied rank-sum base is `3 sqrt(4*4)=12`.
This is a strict reduction caused by noncommuting square-zero
coefficients, with every omitted edge checked for all external
inputs. It remains worse than the direct scalar crossing
certificate, and is far above `2^(sqrt(2))`; it is not an improved
factoring construction.

The scalar-twist control works over `F_3`. Let

```math
X=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
Z=\begin{pmatrix}1&0\\0&-1\end{pmatrix},
\qquad ZX=-XZ.
```

For the tensor `T` supported only at external input `(0,0,0,0)`,
the three pieces `T XZ`, `-T XZ`, and `E_2 I_2` sum to the
crossing. The first two have one channel in each cut, with `X`
in the left factor and `Z` in the right. They supply a nontrivial
phase in (2). Dropping that phase changes an explicitly tested
two-copy factor entry.

`check_twisted_channel_walks.py` verifies 6,192 all-input scalar
commutation tests, 11,520 all-input zero-product tests, 1,728
ordered gathered factor identities, and eight type-walk counts.
It also checks the reciprocal weighting and Cauchy bound exactly
on the complete-graph profiles `(r_1,c_1)=(1,4)` and
`(r_2,c_2)=(4,1)`, where rational weights `(2,1/2)` recover base
four. An invalid positive-supervector certificate is rejected.
The recorded run in `twisted-channel-walks.json` completed in
about two seconds. The finite examples validate the interface,
with the limitations stated above.

```sh
python3 tools/check_twisted_channel_walks.py
```

The exact evidence is recorded in `evidence/twisted-channel-walks.json`.
