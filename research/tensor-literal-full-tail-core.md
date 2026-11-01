# A literal full-tail reduction for general low physical frames

## Result and scope

A fixed homogeneous physical code on the lowest `k` hot modes, with the
complete occupation alphabet on the other `M-k` hot modes, has an exact
profile calculation on **8k active Gaussian modes**. The code need not
be closed under deletion. Both within-cut Gaussian factors are transported
explicitly instead of being discarded.

Assume `k,M` are powers of two and `M>=4k`. Let the low code have rank
`r`, with constant height-homogeneous maps `P,L` satisfying `LP=I_r`.
The physical alphabet has

\[
q=r2^{M-k}.
\]

For each cut, there is an explicit polynomial active matrix with:

| Quantity | Bound or exact value |
|---|---:|
| Number of active Gaussian modes | `8k` |
| Ambient matrix order | `2^(8k)` |
| Generic rank | `r^2 2^(6k)` |
| Coefficient degree in `t` | At most `32kM` |
| Remaining factors | Explicit diagonal Gaussian spectators |

Consequently, for fixed `k`, both complete natural profiles can be
represented and calculated in time polynomial in `M`, with constants
exponential in `k`. This statement concerns profile calculation. It
does not remove the physical alphabet size, the tensor degree, or the
costs in a complete factoring transfer.

The [previous full-tail method](coherent-full-tail-core.md) used `4k`
active modes for the compound matrix, but required deletion invariance
when passing from the literal reflection tensor to that compound.
The present construction removes that requirement by retaining the
internal Gaussian factors. Sections 1-6 first prove the homogeneous
projector version. The canonical embedding in Section 7 also removes
height homogeneity: it handles arbitrary fixed Laurent low-mode frames,
with only an endpoint split and with every clearing shift retained.

The ambient reflection interval here is explicitly **4M**, with the cold
modes frozen to vacuum. Thus the retained M-mode literal kernel has zero
same-leg blocks. This is not a claim about the literal tensor on all M
physical modes of interval `2M`, which has additional same-leg factors.

## 1. The polynomial crossing and its literal Gaussian factorization

Put `t=w^2`, and give hot mode `a` the inherited height `2a+1`.
On M hot modes, define

\[
\begin{aligned}
J_{ab}&=[a=b]t^{2a+1},\\
D_{ab}&=(1+t)t^{2(a\mathbin{|}b)+1},\\
A_{ab}&=(1+t)t^{2(a\mathbin{\&}b)},\\
F_{ab}&=(1+t)t^{2(a\mathbin{\&}\neg b)}.
\end{aligned}
\tag{1}
\]

The [literal hot kernel](tensor-literal-low-cell-retractions.md) has
the block matrix

\[
\begin{pmatrix}
0&A&D&J\\ A&0&J&D\\ D&J&0&A\\ J&D&A&0
\end{pmatrix}
\]

in role order `x,y,u,v`. After reversing the two column labels to `v,u`,
its occupation matrix factors as

\[
T_0=E_A\,\Gamma(B')\,E_A^{T},\qquad
B'=\begin{pmatrix}J&D\\D&J\end{pmatrix}.
\tag{2}
\]

Here `Gamma` is the occupation action of a one-particle matrix, and

\[
E_A=\prod_{a,b}(1+A_{ab}c_{x_a}c_{y_b}).
\]

The product is finite in the occupation algebra. Its transpose acts by
annihilators on the right. Formula (2) follows by expanding the selected
Pfaffian into within-row pairings, cross pairings, and within-column
pairings. In characteristic two there are no signs. Both E factors are
regular polynomial units; their inverses are the same products.

As in the earlier compound calculation,

\[
B'=L_0\Delta R_0,\quad
L_0=\begin{pmatrix}I&0\\F&I\end{pmatrix},\quad
R_0=\begin{pmatrix}I&F^T\\0&I\end{pmatrix},\quad
\Delta=\begin{pmatrix}J&0\\0&t^{2M}J\end{pmatrix}.
\tag{3}
\]

Indeed, `FJ=JF^T=D` and `F^2=(1+t^(2M))I`. The latter identity follows
by multiplying the two-by-two binary-digit factors of F. Thus all four
factors surrounding Delta in (2)-(3) are regular units.

For interval `4M`, the original common scalar is `t^(4M^2)`. A raw
physical coefficient is

\[
t^{4M^2-h(x)-h(y)}T_0(t;x,y,u,v).
\]

Each occupation height is at most `M^2`. Every regular coefficient has
degree at most `h(x)+h(y)+h(u)+h(v)`. Therefore the displayed raw
coefficient is polynomial and has degree at most `6M^2` in t.
Constant physical compression does not increase that degree. At t=1
the kernel is exactly the crossing, and `LP=I` preserves its endpoint.
The opposite scalar shifts `+4M^2` and `-4M^2` cancel in the two-place
counting product. They are not dropped from one place alone.

## 2. The extra reflected-mode identity

Let R denote reversal of the M mode indices: `R_ab=[a+b=M-1]`.
The additional identity needed for the literal factors is

\[
FA=AF^T=(1+t^{2M})R.
\tag{4}
\]

For example, the `(d,a)` entry of FA is

\[
(1+t)^2\sum_c t^{2(d\mathbin{\&}\neg c+c\mathbin{\&}a)}.
\]

The sum factors over the binary digits of c. A digit factor is zero
when the corresponding digits of d and a agree, and is `1+t^(2^(j+1))`
when they differ. Thus only `d=M-1-a` survives. The product telescopes
to `1+t^(2M)`, proving (4).

Conjugation through the literal factors in (2) now gives the following
complete pre-concentration row images, for `a<k`:

\[
\begin{aligned}
c_{x_a}&\mapsto c_{x_a}+\sum_b F_{ba}c_{y_b},\\
c_{y_a}&\mapsto c_{y_a},\\
a_{x_a}&\mapsto a_{x_a}+\sum_b A_{ab}c_{y_b},\\
a_{y_a}&\mapsto a_{y_a}+\sum_b F_{ab}a_{x_b}
 +\sum_b A_{ab}c_{x_b}+(1+t^{2M})c_{y_{M-1-a}}.
\end{aligned}
\tag{5}
\]

The corresponding column images are

\[
\begin{aligned}
c_{v_a}&\mapsto c_{x_a}+\sum_b A_{ab}a_{y_b},\\
a_{v_a}&\mapsto a_{x_a}+\sum_b F_{ba}a_{y_b},\\
c_{u_a}&\mapsto c_{y_a}+\sum_b F_{ab}c_{x_b}
 +\sum_b A_{ab}a_{x_b}+(1+t^{2M})a_{y_{M-1-a}},\\
a_{u_a}&\mapsto a_{y_a}.
\end{aligned}
\tag{6}
\]

These formulas include the creation terms in the row annihilators and
the annihilation terms in the column creators. Omitting those terms
would return to the compound-only calculation.

## 3. Two-direction weighted concentration

Write `b=c+kB`, with `0<=c<k`, `0<=B<N=M/k`, and put `x=t^(2k)`.
For fixed low index a, the F coefficients in (5)-(6) have tail dependence
either 1 or `x^B`. The A coefficients have tail dependence 1. The extra
reflected indices lie in the last block `B=N-1`.

Keep the blocks `B=0,N-1` individually. On the middle blocks
`B=1,...,N-2`, write `j=B-1` and `ell=N-2`. Thus ell is even. Define

\[
S=\sum_{j=0}^{\ell-1}x^j,\qquad T=1+x^{\ell-1}.
\]

Both are units at zero. On each residue class c, use a middle-block
change of basis C with first and last columns

\[
C_{j,0}=x^j/S,\qquad C_{j,\ell-1}=1/S.
\]

For an interior column `0<j<ell-1`, use

\[
C_j=e_j+
 \frac{1+x^{\ell-1-j}}{T}e_0+
 \frac{x^{\ell-1-j}+x^{\ell-1}}{T}e_{\ell-1}.
\tag{7}
\]

Let `ell_0=(1,...,1)` and `ell_1=(x^(ell-1-j))_j`. Direct multiplication,
using the even parity of ell, gives

\[
\ell_0C=e_0^T,\qquad \ell_1C=e_{\ell-1}^T,
\qquad C^{-1}(x^j)_j=Se_0,\qquad C^{-1}{\bf1}=Se_{\ell-1}.
\tag{8}
\]

C is invertible: its endpoint columns have the displayed independent
pairings, and its interior columns have independent interior coordinates.
Also `C(0)` is upper triangular with unit diagonal.

Every entry obeys the weighted valuation bound

\[
\mathrm{val}_t(C_{ij})\ge\max(0,2k(i-j)).
\tag{9}
\]

For the first column this follows from `x^i`; for the last column it is
automatic; and for an interior column the only entry below the diagonal
is divisible by `x^(ell-1-j)`. The weighted conjugate of C is therefore
integral and has unit determinant. Its inverse is integral as well, so
`U=C^{-1}` obeys the same bound.

Extend these changes by the identity on `B=0,N-1` and on both mode
groups. Set `V=Delta^{-1} C Delta`. The matrices U,V and their inverses
are regular, and `U Delta V=Delta`. Their occupation actions preserve
Smith profiles. Applying (8) to (5)-(6) leaves only four blocks per
residue and group:

\[
B\in\{0,1,N-2,N-1\}.
\tag{10}
\]

This is the claimed 8k-mode support. The reverse-power row in (8) is
needed for the constant column annihilator coefficients after the
weighted conjugation by Delta. Concentrating only the earlier
all-ones/power-vector pair would not cover the literal operators.

## 4. Polynomial operator images and the endpoint

Although S and T can vanish at t=1, all denominators cancel from the
transported low operators. The implementation gives their complete
polynomial formulas. To state them compactly, label the four retained
positions by `0,L,H,E`, corresponding to (10). For either group z, put

\[
\begin{aligned}
C_z^F(c)&=c_{z,c,0}+xS c_{z,c,L}+x^{N-1}c_{z,c,E},\\
C_z^A(c)&=c_{z,c,0}+S c_{z,c,H}+c_{z,c,E},\\
A_z^F(c)&=a_{z,c,0}+a_{z,c,L}+a_{z,c,E},\\
\widehat C_z^F(c)&=c_{z,c,0}+S c_{z,c,L}+c_{z,c,E},\\
\widehat A_z^F(c)&=a_{z,c,0}+x a_{z,c,L}+x^{N-1}a_{z,c,E},\\
\widehat A_z^A(c)&=a_{z,c,0}+a_{z,c,H}+a_{z,c,E}.
\end{aligned}
\]

Write `f_ca=(1+t)t^(2(c&~a))`, `g_ac=(1+t)t^(2(a&c))`, and
`R_t=1+t^(2M)`. The active row images are

\[
\begin{aligned}
c_{x_a}&\mapsto c_{x,a,0}+\sum_c f_{ca}C_y^F(c),\\
c_{y_a}&\mapsto c_{y,a,0},\\
a_{x_a}&\mapsto a_{x,a,0}+\sum_c g_{ac}C_y^A(c),\\
a_{y_a}&\mapsto a_{y,a,0}+\sum_c f_{ac}A_x^F(c)
 +\sum_c g_{ac}C_x^A(c)+R_t c_{y,k-1-a,E}.
\end{aligned}
\tag{11}
\]

The active column images are

\[
\begin{aligned}
c_{v_a}&\mapsto c_{x,a,0}+\sum_c g_{ac}\widehat A_y^A(c),\\
a_{v_a}&\mapsto a_{x,a,0}+\sum_c f_{ca}\widehat A_y^F(c),\\
c_{u_a}&\mapsto c_{y,a,0}+\sum_c f_{ac}\widehat C_x^F(c)
 +\sum_c g_{ac}\widehat A_x^A(c)+R_t a_{y,k-1-a,E},\\
a_{u_a}&\mapsto a_{y,a,0}.
\end{aligned}
\tag{12}
\]

Every coefficient has degree at most `2M`. Every added term vanishes at
t=1, since each contains `1+t` or `R_t`. Thus (11)-(12) specialize to the
original low creation and annihilation operators. The endpoint is checked
on these polynomial formulas; no rational concentration matrix is
evaluated at its pole.

## 5. The active matrix, both cuts, and the complete cost

Let `Pi=PL`. It is a constant projector of rank r. Because the maps
preserve height, the natural physical gauges pass through P and L.
The literal compressed cut therefore has the same nonzero Smith factors
as the full matrix sandwiched by its physical projectors. Extending P by
a constant basis of `ker L` proves this assertion directly.

At zero, substitute (11)-(12) into four copies of Pi. At infinity, use
Pi on the x,v roles and `Pi^T` on the other paired roles, with the
corresponding role exchange of the literal kernel. The same operator
calculation applies. The projectors are supported only on the 8k modes
in (10), so the full sandwiched matrix factors as

\[
K_{M,k}\otimes\Gamma(\Delta_{\rm spectator}).
\tag{13}
\]

The active physical indices are

\[
\mathcal A=\{0,\ldots,2k-1\}\cup\{M-2k,\ldots,M-1\}.
\]

Their Gaussian diagonal orders are `2j+1` and `2M+2j+1` for
`j in A`. The spectator orders are the same two expressions for
`2k<=j<M-2k`.

At the endpoint, the active matrix is Pi on each of the two low mode
groups, or Pi and its transpose for the other cut, tensored with the
identity on 6k remaining active modes. Its rank is exactly
`r^2 2^(6k)`. The active row projector gives the same generic upper
bound, so this endpoint supplies the complete generic rank. Including
the spectators gives

\[
r^2 2^{6k}\,2^{2M-8k}=r^2 2^{2M-2k}=q^2,
\]

as required.

Every matrix unit on k occupation modes is a product of at most 2k
creators and annihilators: the single-mode diagonal units are `ca`
and `1+ca`, and the off-diagonal units are c and a. Thus each transported
physical projector has coefficient degree at most `4kM`. The four
projectors contribute at most `16kM`. The sum of all active diagonal
orders is exactly `16kM`. Therefore

\[
\deg_t K_{M,k}\le32kM.
\tag{14}
\]

The active order is constant for fixed k. Polynomial matrix arithmetic
therefore constructs its complete profiles in time polynomial in M.
For example, a nonzero maximal-rank minor has degree at most
`r^2 2^(6k) 32kM`, which is a finite sufficient precision bound for
resolving every nonzero Smith order. A capped or failed smaller
calculation must still be identified as such.

If the active orders at a cut are `nu_i`, its complete channel polynomial
has the compact representation

\[
Z_{M,k}(z)=\left(\sum_i z^{\nu_i}\right)
 \prod_{j=2k}^{M-2k-1}(1+z^{2j+1})(1+z^{2M+2j+1}).
\tag{15}
\]

If the low physical heights are `h_1,...,h_r`, the physical partition is

\[
A_{M,k}(z)=\left(\sum_{i=1}^r z^{h_i}\right)
 \prod_{j=k}^{M-1}(1+z^{2j+1}).
\tag{16}
\]

Both (15) and (16), the raw degree bound, and the actual finite field
must enter any pressure or uniform-transfer claim. The active reduction
alone is not an improved factoring exponent.

## 6. Independent checks and complete small cores

The [operator constructor](../tools/tensor_literal_full_tail_core.py)
provides (11)-(12) for arbitrary supported k,M. Its dense active-matrix
implementation is deliberately limited to k=1. It computes the vacant
and occupied lowest-mode codes, including the occupied code that fails
deletion closure.

The [checker](../tools/check_tensor_literal_full_tail_core.py) independently
constructs C and its inverse with a common denominator, performs the
full 2M-mode conjugations before restricting support, and compares the
result with (11)-(12). The retained
[receipt](../evidence/tensor-literal-full-tail-controls.json) covers:

- 5,920 entries in the two inverse products;
- 5,920 weighted-regularity entries;
- 5,920 entries in the two reflected-mode products;
- 10,624 complete transported operator coefficients;
- 1,128 canonical anticommutator pairings;
- 6,016 endpoint operator coefficients;
- 128 complete active Smith orders compared with independently
  constructed literal Pfaffian profiles at M=4.

The last comparison includes both the deletion-closed vacant code and
the noncoherent occupied code. The controls use zero external numerical
weights, so the native Smith pivot weights are exact integers. The
[completed k=1 archive](../evidence/tensor-literal-full-tail-k1-profiles.json)
contains all 12 vacant/occupied cores at M=4,8,16,32,64,128, totaling
768 complete active orders. The largest represented physical alphabet
has `2^127` letters; its full profile is retained through the 64-order
core and the explicit spectator product, not as an exponential list.

```sh
python3 tools/check_tensor_literal_full_tail_core.py
python3 tools/tensor_literal_full_tail_core.py --modes 4 8 16 32 64 128 \
  --output evidence/tensor-literal-full-tail-k1-profiles.json
```

No matrix of order `2^(8k)` for k>1 is claimed evaluated here. The
ordinary proof supplies the general reduction; the finite checks test
its algebra and the implemented small cores.

## 7. Canonical embeddings remove height homogeneity

The projector argument in Section 5 is sufficient but not necessary.
There is a direct way to retain the actual rectangular physical maps,
including unequal heights and parameter dependence.

Let `d=2^k`. Take arbitrary fixed Laurent-polynomial maps
`L(t):F^d -> F^r` and `P(t):F^r -> F^d`, with only

\[
L(1)P(1)=I_r.
\tag{17}
\]

Tensor them with the identity on the remaining M-k hot modes. Assign
the low virtual labels any fixed integer heights `g_1,...,g_r`.
The complete physical height is `g_i+h(tail)`. Gauges that couple a
low label to the tail configuration are outside this locality statement.
Write

\[
D_h=\mathrm{diag}(t^{h(S)})_{S\subseteq[k]},\qquad
H=\mathrm{diag}(t^{g_i}),\qquad Q^*(t)=Q(t^{-1}).
\]

The physical gauge convention is the unflattened weight `g(x)+g(y)`
at zero and `g(u)+g(v)` at infinity. It is the natural reflection
convention in
[`logical_heights`](../tools/reflection_dual_rail_gaussian.py), whose
infinity-active roles are u,v, and in the
[protected-height calculation](coherent-protected-height25.md). In the
infinity flattening, one of these weights is on a row role and the
other is on a column role.

The actual zero-place row map on a low physical leg is `H L D_h^{-1}`;
its column map is P. To derive the other place directly, write the raw
coefficient as `F=t^(C-h(x)-h(y)) T_0`. The reciprocal identity gives

\[
F(t^{-1})=t^{-C-h(u)-h(v)}T_\infty(t).
\]

After contracting the starred physical maps and multiplying by
`t^(g(u)+g(v))`, the two factors `H` therefore belong to u and v.
This gives exactly the infinity line of the table below. The tail
gauges cancel through the identity tail maps, so all changes remain
supported on the lowest k modes. For paired types, x,v share one type
and y,u share the other; the total mean charge of both cuts is
`2(E_X g+E_Y g)`, with every virtual height included. Moving the
infinity weight from v to x gives the same paired-type mean but is
not asserted to preserve the local Smith profile.

Fix a constant coordinate inclusion `E:F^r -> F^d`, with `E^T E=I_r`.
For any rectangular row and column maps R,C, the full sandwich

\[
(ER)\,T\,(CE^T)=E(RTC)E^T
\tag{18}
\]

is the desired compressed matrix with zero rows and columns added.
Its nonzero Smith factors are exactly those of RTC. This statement
does not require R,C to split over the local ring. The inclusion E is
a calculation device after the physical gauges have already been
applied; it does not reassign any physical height.

For the two-leg cuts, use the tensor product of two such inclusions.
The four low endomorphisms to substitute into the active transport are:

| Cut | Row x | Other row | Column v | Other column |
|---|---|---|---|---|
| Zero | `E H L D_h^-1` | `E H L D_h^-1` | `P E^T` | `P E^T` |
| Infinity | `E L*` | `E H (P*)^T D_h^-1` | `D_h^-1 P* H E^T` | `(L*)^T E^T` |

The second row and fourth column in the infinity line refer to the
original u and y roles after the usual role exchange. The star reverses
every parameter exponent. This reversal is essential for moving frames.

Every low endomorphism is a polynomial of degree at most 2k in the low
creators and annihilators, with Laurent-polynomial coefficients. Clear
the four endomorphisms at a particular cut by nonnegative monomial
powers `c_1,...,c_4`. Let their cleared coefficient degrees be
`d_1,...,d_4`, and set `C_cut=sum c_i`. Substitution of (11)-(12) gives
a polynomial active matrix with

\[
\begin{aligned}
\text{ambient order}&=2^{8k},\\
\text{generic rank}&=r^2 2^{6k},\\
\deg_t K^{\rm cleared}_{M,k}&\le32kM+\sum_i d_i.
\end{aligned}
\tag{19}
\]

If its nonzero Smith orders are `nu_i`, the actual weighted orders
are `nu_i-C_cut`, followed by the same spectator subset sums as in
(15). The two cuts can have different clearing powers. Each is
subtracted from its own spectrum; they are not assumed to cancel.

To prove the rank assertion, observe that each row endomorphism has
rank at most r on its low group. This gives the generic upper bound.
At t=1, the transported operators return to their original low modes,
and (17) makes the product on each low group `EE^T`. The endpoint rank
is therefore exactly `r^2 2^(6k)`. This argument uses no all-parameter
identity LP=I.

The physical crossing itself is still effective and polynomial after
one common monomial clearing. If `t^m_L L` and `t^m_P P` are polynomial,
applying those four cleared physical maps multiplies the Laurent
crossing by `t^(2m_L+2m_P)`. It preserves the endpoint and adds opposite
scalar shifts at zero and infinity. A conservative tensor-degree bound
is

\[
6M^2+2\deg(t^{m_L}L)+2\deg(t^{m_P}P).
\tag{20}
\]

This common physical clearing is distinct from the cut-specific
calculation clearings in (19). Both are explicitly accounted for.
For fixed k and fixed low maps and heights, (19) gives a complete
profile calculation polynomial in M. The physical partition remains
exactly (16), using the assigned g heights rather than the ambient
heights of mixed map components.

## 8. A practical nonhomogeneous one-mode implementation

For k=1 and r=1, the active matrix has order 256 and rank 64 for every
M. The [rank-one driver](../tools/tensor_literal_full_tail_rank_one.py)
implements all constant F4 pairs `LP=1`, including maps mixing the
vacuum and occupied lowest mode. There are 20 normalized pairs and
13 Frobenius representatives. The assigned low virtual height is zero;
changing it merely makes a common shift of every physical height and
the matching weighted spectra, leaving the pressure equation unchanged.

If `P=(p_0,p_1)^T` and `L=(l_0,l_1)`, the zero cut clears each weighted
row by `t^[l_1!=0]`. Its total calculation clearing is therefore
`2[l_1!=0]`. The infinity cut has total clearing `2[p_1!=0]`.
The driver records the cleared and actual orders separately and
computes both cuts directly. It does not infer equality from an
orthogonality test on an unequal-height map.

The independent
[checker](../tools/check_tensor_literal_full_tail_rank_one.py) constructs
the original raw Laurent reflection tensor at M=4, contracts the actual
physical maps, applies the assigned weights, and compares its complete
64-order cuts with the 256-order active calculation. It also checks
the literal matrices with independent binary Toeplitz ranks. This tests
the canonical embedding, both different clearing powers, and the
nonhomogeneous gauge accounting separately from the operator-support
identities.

Its [completed control receipt](../evidence/tensor-literal-full-tail-rank-one-controls.json)
covers all 13 F4 representatives: 106,496 original raw-Laurent tensor
entries and endpoints, 1,664 complete Smith-order comparisons, and
1,042 independent Toeplitz-rank checks. Both cuts are computed in each
case, including different zero and infinity clearing powers.

```sh
python3 tools/check_tensor_literal_full_tail_rank_one.py
python3 tools/tensor_literal_full_tail_rank_one.py --modes 4 8 16 32 64 128 \
  --output evidence/tensor-literal-full-tail-f4-rank-one-profiles.json
```

The completed [small-M archive](../evidence/tensor-literal-full-tail-f4-rank-one-small-profiles.json)
and [large-M archive](../evidence/tensor-literal-full-tail-f4-rank-one-large-profiles.json)
cover all 13 representatives at M=4,8,16,32,64,128. These 78 records
contain 9,984 complete active orders, with both cuts directly computed
and no capped tails. They cover 120 normalized field pairs across the
six M values. The spectator products represent their complete physical
profiles.

| M | Best numerical pressure fit |
|---:|---:|
| 4 | 0.838523759044 |
| 8 | 0.800657042415 |
| 16 | 0.768698103476 |
| 32 | 0.745635296991 |
| 64 | 0.730280856772 |
| 128 | 0.720621592460 |

The best retained case at each M is the vacant coordinate low mode.
None of these finite fits improves the pressure frontier. The searches
use a logarithmic scale grid and bounded refinements; they do not certify
global numerical minima or exclude other low frames. Their purpose here
is to apply the uniform nonhomogeneous reduction to complete large-tail
profiles without constructing the exponentially large physical matrices.

The constructor and the uniform lemma do not by themselves give a
stronger rho or factoring exponent. No Lean formalization is added.
