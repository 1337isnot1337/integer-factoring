# Exact weighted aggregation of complete folded-Gram channels

## Result

The [complete folded-Gram construction](folded-gram-channels.md) can be
evaluated without enumerating double-residue and single-residue patterns.
This note aggregates those patterns into the same small Frobenius blocks.
The aggregation retains the exact residue energy, so it supports shells far
below central energy.

The two cuts have different energy conditions:

\[
 W_R=kM+(2j-LM),\qquad W_C=kM.
\]

Lowering the shell energy \(j\) changes the first condition while leaving
the second fixed. The exact algorithm below preserves that tradeoff.
It also computes the shell label dimension, so the full rank-sum bound can
be compared with its actual finite seed threshold.

This is an evaluation improvement for the already proved complete channel
bounds. It does not assert that those bounds equal the actual piece ranks
or establish a better factoring exponent.

## 1. The exterior module behind one residue

For one low residue, a folded generator records whether it belongs to the
unchanged set, the creator set, both, or neither. Under the two-dimensional
matrix used in the Johnson-rank proof, these four possibilities form

\[
 \Lambda V=\mathbf1\oplus V\oplus\det V.
 \tag{1}
\]

The scalar summands correspond to zero and two occupied modes. The natural
two-dimensional summand distinguishes the two single-occupation choices.
For \(N=M/2\) residues the full module is

\[
 \Lambda(V\otimes\mathbb F_2^N).
 \tag{2}
\]

In exterior degree \(2k\), the middle torus weight consists exactly of
ordered pairs of \(k\)-subsets \((\eta,C)\).

The Frobenius block decomposition in the earlier note is constant and
preserves the number of particles at every residue. It therefore preserves
both total exterior degree and total residue energy.

## 2. Unweighted multiplicities

In exterior degree \(2k\), let \(r\) be the number of occupied coordinates
of the first species. There are

\[
 \binom Nr\binom N{2k-r}
\]

such basis vectors. The torus character is

\[
 \chi_{N,k}(Y)=
 \sum_{r=0}^{2k}\binom Nr\binom N{2k-r}Y^r.
 \tag{3}
\]

For the explicit Frobenius block \(B_m\), let

\[
 b_m(Y)=\prod_{w\in\mathcal W_m}(1+Y^w).
\]

Its determinant twist \(d\) has character \(Y^d b_m(Y)\).
Using the elementary tensor splitting already proved, the degree-\(2k\)
part of (2) has a direct-sum decomposition

\[
 \chi_{N,k}(Y)=
 \sum_{m+2d=2k}\nu_{N,k,m}\,Y^d b_m(Y),
 \qquad \nu_{N,k,m}\in\mathbb Z_{\geq0}.
 \tag{4}
\]

The multiplicities can be recovered by triangular subtraction:
process \(m=2k,2k-2,\ldots,0\); at \(d=k-m/2\), the remaining coefficient
of \(Y^d\) is \(\nu_{N,k,m}\). Subtract that many copies of \(Y^db_m(Y)\).
The leading coefficient of every \(b_m\) is one, so the procedure is unique.
Its nonnegativity follows from the constructive decomposition, not from
an unproved character ansatz.

As a consistency identity,

\[
 \sum_{m+2d=2k}\nu_{N,k,m}\,
       \dim(B_m)_{\text{middle}}
 =\binom Nk^2.
 \tag{5}
\]

## 3. The exact residue-energy character

Let

\[
 E_r(X)=
 \sum_{\substack{A\subseteq\mathcal R_M\\|A|=r}}
 X^{\sum_{a\in A}a}.
 \tag{6}
\]

The residue energy of a pair of species subsets is the sum of their two
energies. It counts a residue twice when both species occupy it, as required.
Thus the character restricted to residue energy \(W\) is

\[
 \chi_{N,k,W}(Y)=
 \sum_{r=0}^{2k}
 \bigl[X^W\bigr]\bigl(E_r(X)E_{2k-r}(X)\bigr)Y^r.
 \tag{7}
\]

Because the explicit decomposition preserves residue energy, there are
nonnegative integer multiplicities

\[
 \chi_{N,k,W}(Y)=
 \sum_{m+2d=2k}
 \nu_{N,k,m}(W)\,Y^d b_m(Y).
 \tag{8}
\]

They are recovered by the same triangular subtraction. Summing (8) over
\(W\) gives (4). At the middle weight, a second useful identity is

\[
 \sum_{m+2d=2k}\nu_{N,k,m}(W)\,
       \dim(B_m)_{\text{middle}}
 =\bigl[X^W\bigr]E_k(X)^2.
 \tag{9}
\]

The right-hand side is an ordinary integer convolution. It is not a
characteristic-two square: these coefficients count the number of residue
patterns multiplying a binary matrix rank.

## 4. Gaussian-binomial evaluation

For \(N\) low odd modes,

\[
 E_r(X)=X^{r^2}
 \begin{bmatrix}N\\r\end{bmatrix}_{X^2}.
 \tag{10}
\]

Indeed, writing the selected modes as \(2a+1\) and subtracting the minimum
indices \(0,1,\ldots,r-1\) gives the usual subset generating polynomial.

Write \(G_{N,r}(q)\) for the Gaussian binomial in (10). Its coefficients
are computed by the finite recurrence

\[
 G_{n,r}(q)=G_{n-1,r}(q)+q^{n-r}G_{n-1,r-1}(q),
 \quad G_{n,0}=G_{n,n}=1.
 \tag{11}
\]

Complement symmetry gives \(G_{N,r}=G_{N,N-r}\), so only half the
polynomials need to be stored.

For \(s=2k-r\), the coefficient needed in (7) is zero unless
\(W-r^2-s^2\) is a nonnegative even integer. Otherwise it is

\[
 [q^{(W-r^2-s^2)/2}]\,G_{N,r}(q)G_{N,s}(q).
 \tag{12}
\]

Only one convolution coefficient is needed. The symmetry in \(r,s\)
also halves the character work.

The dimension of the actual \(L\)-part shell at energy \(j\) is

\[
 d(M,L,j)=
 [q^{(j-L^2)/2}]\,G_{N,L}(q),
 \tag{13}
\]

with zero for a negative or nonintegral exponent.

## 5. Aggregated complete block ranks

For \(m+2d=2k\), use the weighted middle-subset indices
\(\alpha,\beta\) of \(B_m\). Let

\[
 w(\alpha,\beta)=
 \sum_{a\in\alpha\cap\beta}w_a.
\]

Define

\[
 \beta_{L,k,m}(v)=
 \operatorname{rank}_{\mathbb F_2}
 \left[
 \binom{L-d-w(\alpha,\beta)}{L+v-k-d}
 \right]_{\alpha,\beta}.
 \tag{14}
\]

All tops are nonnegative because
\(d+w(\alpha,\beta)\leq d+m/2=k\leq L\).

The number of doubled residues has disappeared from (14). In the
unaggregated formula, a component with \(i\) double residues contributes a
determinant shift \(i\), while its single-residue tensor block contributes
another shift \(d'\). Their sum is always

\[
 i+d'=k-m/2=d.
\]

Thus grouping equal blocks before elimination is exact.

With \(\delta=2j-LM\), the complete channels are

\[
 \boxed{
 \begin{aligned}
 R_{Mv}
 &\leq
 \sum_{k=0}^L\sum_{m+2d=2k}
 \nu_{N,k,m}(kM+\delta)\,\beta_{L,k,m}(v),\\
 C_{Mv}
 &\leq
 \sum_{k=0}^L\sum_{m+2d=2k}
 \nu_{N,k,m}(kM)\,\beta_{L,k,m}(-v).
 \end{aligned}}
 \tag{15}
\]

These are exactly the same numbers as the residue-pattern sums in the
earlier note. No angle sum is added: all source angles were already combined
in the complete Gram polynomial.

The unweighted formulas replace \(\nu(W)\) by \(\nu\) from (4).
All binary block ranks in (14) depend on \(L,k,m,v\), not on \(N\), \(M\),
or \(j\). They can therefore be reused while scanning different cutoffs and
shell energies.

For a finite dimension \(d(M,L,j)\), cap each cut by \(d^2\). A wholly
integer upper bound for the rank sum is

\[
 S_{\rm int}=
 \sum_v
 \left\lceil
 \sqrt{\min(d^2,R_{Mv})\min(d^2,C_{Mv})}
 \right\rceil.
 \tag{16}
\]

This is safe to compare with the seed threshold \(d^{\sqrt2}\).
It is an upper bound obtained from the construction, not a measured optimal
rank sum.

## 6. Implementation and checks

Files:

- [Fast aggregation module](../tools/folded_gram_aggregate.py).
- [Independent aggregation checker](../tools/check_folded_gram_aggregate.py).
- [Aggregation evidence](../evidence/folded-gram-aggregation.json).

The main interfaces are:

    complete_gram_channels(M, L, j)
    unweighted_complete_gram_channels(N, L)
    shell_dimension(M, L, j)
    gaussian_distributions(N)[L]

The final expression returns the coefficient list of \(G_{N,L}\).
Coefficient \(q\) corresponds to the shell energy \(j=L^2+2q\).

The checker compares multiplicities with direct sums over residue patterns
and the original Frobenius tensor blocks. It checks Gaussian-binomial
coefficients against explicit subsets, and compares the complete weighted
and unweighted channel profiles with the unaggregated implementation.
Large examples use exactly the same integer formulas.

The code makes dense and noncentral scans practical. For example, the
weighted \(M=256,L=64,j=8192\) calculation has
45374927594179648855944192519456556 labels. No tensor of that dimension is
constructed; the computation uses subset coefficient distributions,
nonnegative block multiplicities, and small binary matrices.

## 7. A coarse noncentral screen

The [noncentral screen](../tools/screen_folded_gram_noncentral.py) evaluates
13 energies for each of \(M=256,L=32\) and \(M=256,L=64\).
It uses (16), with exact integer square-root ceilings throughout.
The [complete evidence](../evidence/folded-gram-noncentral-M256.json) records
all uncapped channel bounds, exact shell dimensions, and integer scores.

For each length, the best sampled value of
\(\rho=\log_2(S_{\rm int})/(2\log_2d)\) occurs at central energy:

| \(M,L\) | Best sampled \(j\) | \(\rho\) from the upper bound | \(\log_2(d^{\sqrt2}/S_{\rm int})\) |
|---|---:|---:|---:|
| \(256,32\) | 4096 | 0.808467918 | -18.525280681 |
| \(256,64\) | 8192 | 0.894349991 | -43.113666791 |

The required value is \(\rho<1/\sqrt2\). No tested upper bound meets it.
This is not an obstruction to smaller actual ranks.

Reproduce:

    python3 tools/screen_folded_gram_noncentral.py \
      --output evidence/folded-gram-noncentral-M256.json

An additional [516-case energy screen](../evidence/weighted-gram-energy-screen.json)
uses exact square-root ceilings at scale `10^9`. The
[reproduction script](../tools/explore_weighted_gram_costs.py) also writes
every capped channel profile. Its six length/cutoff grids give:

| \(M,L\) | Energies checked | Best sampled energy | \(\rho\) from the upper bound |
|---|---:|---:|---:|
| \(64,8\) | 48 | 256 | 0.916739488 |
| \(64,12\) | 60 | 384 | 0.928530778 |
| \(64,16\) | 64 | 512 | 0.998457683 |
| \(128,16\) | 96 | 1024 | 0.846972846 |
| \(128,24\) | 120 | 1536 | 0.864493411 |
| \(128,32\) | 128 | 2048 | 0.932316569 |

Central energy is best in every sampled grid. All comparisons use actual
shell dimensions. The integer check for an improved seed uses
`upper^985 < d^1393 * scale^985`, where `1393/985 < sqrt(2)`;
none succeeds. These failures concern the displayed constructive upper
bounds, not the unknown optimal piece ranks.
