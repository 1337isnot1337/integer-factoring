# Arbitrary identity selections above a prime power

## Result and scope

Let \(p\) be prime, \(N=p^k\) with \(k\ge1\), and
\(n=N+d\), where \(1\le d\le N\). Select any \(r\) distinct
physical identity directions in the interval of length \(n\).
Every occupation-preserving operator on those directions has a
regular zero reduction to a square core with at most
**\(2(r+2d)\) actual Gaussian modes**. Its complete zero profile
can be obtained from a Weyl matrix of order at most **\(p^{r+2d}\)**.
No formal modes or Gaussian-factor division are required.
For a global crossing, the original selected operator must satisfy the
existing invertible energy-filtered endpoint condition. The local
reduction itself applies to every occupation-preserving operator.

For the last \(r\) directions, with \(d\le r\le N\), the
stronger count is \(2(r+d)\) actual modes and a Weyl matrix of
order \(p^{r+d}\). The binary \(d=1\) case is the separately
implemented [above-power operator family](above-power-identity-operator-reduction.md).

The additional identity directions are existing physical directions
whose operators remain the identity. They add no global polynomial
degree or endpoint condition. This note proves the reduction and
provides exact polynomial geometry. It does not claim a new completed
pressure search or an improved factoring exponent for its enlarged
family.

## 1. Complete the selected directions using actual identities

Put

\[
\tau=z^{2P},\qquad c=1-\tau,
\qquad f_j(i;c)=\sum_h\binom ih\binom jh c^h,
\tag{1}
\]

with the existing odd ramification \(P\). Let \(S\) be the actual
selected set and define

\[
L=\{0,\ldots,d-1\},\qquad
V=\{N,\ldots,N+d-1\},\qquad
J=S\cup L\cup V.
\tag{2}
\]

The letter \(V\) in (2) denotes an index set; the physical variable
\(V_i\) is written with its index below. Let \(O=J\cap[0,N)\)
and \(R=|J|=|O|+d\). Then \(R\le r+2d\). Retain auxiliary
indices

\[
T=\{N-1-j:j\in O\}\cup V.
\tag{3}
\]

There are exactly \(R\) retained auxiliary directions. The augmented
identity directions are also exactly \(J\). Every direction in
\(J\setminus S\) carries the unchanged identity operator.

## 2. Extend the old reflected basis

On the old interval of length \(N\), let
\(H_{ij}=f_j(i;c)\), let \(\mathcal R\) reverse the indices,
and put \(B=H\mathcal R\). The
[prime-power theorem](prime-power-joint-identity-filters.md) gives a
regular Gaussian-preserving change with

\[
B^{-1}H=\mathcal R,\qquad
H^{\mathsf T}D^{-1}BD=c^{N-1}\mathcal R,
\qquad D=\operatorname{diag}_{i<N}(\tau^{i+1}).
\tag{4}
\]

For completeness, the parent reversal identity gives

\[
\operatorname{ord}_\tau B_{ij}\ge\max(0,i-j).
\tag{5}
\]

The determinant of \(H\) is a power of \(c\), since
\(H=P_N\operatorname{diag}(c^h)P_N^{\mathsf T}\) with a unit
triangular Pascal matrix \(P_N\). Thus both \(B\) and
\(D^{-1}BD\), and both inverses, are regular. Extend this change
by the identity on the \(d\) new auxiliary coordinates. Its
Gaussian companions remain regular.

Lucas' theorem gives, for \(i,j<N\) and \(h,k<d\),

\[
\begin{aligned}
f_{N+h}(i;c)&=f_h(i;c),\\
f_j(N+k;c)&=f_j(k;c),\\
f_{N+h}(N+k;c)&=(1+c^N)f_h(k;c)
                    =(2-\tau^N)f_h(k;c).
\end{aligned}
\tag{6}
\]

These equalities hold in every characteristic, including the
characteristic-two cancellation of the constant two.

After (4), the old rows of the coefficient matrices satisfy

\[
\begin{aligned}
(F_X)_{ij}&=\mathbf1_{i=N-1-j},&
(F_V)_{ij}&=c^{N-1}\mathbf1_{i=N-1-j}
&& (j<N),\\
(F_X)_{i,N+h}&=\mathbf1_{i=N-1-h},&
(F_V)_{i,N+h}&=c^{N-1}\mathbf1_{i=N-1-h}.
\end{aligned}
\tag{7}
\]

The new row \(N+k\) in both matrices has entry \(f_j(k;c)\)
on an old column and \((2-\tau^N)f_h(k;c)\) on column \(N+h\).
All coefficient support is therefore contained in (3).

### Both square form matrices are units

Let \(H_d=(f_h(k;c))_{0\le k,h<d}\). The Pascal factorization gives

\[
\det H_d=c^{d(d-1)/2}.
\tag{8}
\]

Every old column in \(O\setminus L\) has its individual old-row
pivot. These pivot rows vanish on all columns in \(L\cup V\).
Eliminating their new-row tails therefore leaves, after harmless
row/column permutations, the special blocks

\[
\begin{pmatrix}I&I\\H_d&(1+c^N)H_d\end{pmatrix},
\qquad
\begin{pmatrix}c^{N-1}I&c^{N-1}I\\H_d&(1+c^N)H_d\end{pmatrix}.
\tag{9}
\]

Their determinants are \((c^N)^d\det H_d\) and
\(c^{d(N-1)}(c^N)^d\det H_d\), respectively. Both are zero units.
In the full matrices, up to a nonzero constant from permutations,

\[
\det F_X=c^{Nd+d(d-1)/2},\qquad
\det F_V=c^{|O|(N-1)}\det F_X.
\tag{10}
\]

No additional interpolation or rationality assumption is used.

## 3. Weighted Weyl reduction and exact physical accounting

At aggregate gauge one, the active Gaussian and forms are

\[
\begin{aligned}
G_{\rm act}
 &=\prod_{j\in J}E(\tau^jY_jU_j)
   \prod_{i\in T}E(\tau^{i+1}V_iX_i),\\
L_j&=\tau^jY_j+c\sum_{i\in T}(F_V)_{ij}V_i,\\
R_j&=U_j-\tau^{-j}\sum_{i\in T}(F_X)_{ij}X_i,
\qquad E(Q)=\sum_{h=0}^{p-1}\frac{Q^h}{h!}.
\end{aligned}
\tag{11}
\]

For each extra actual identity \(j\in J\setminus S\), the removed
factor \(E(-L_jR_j)\) cancels the restored factor \(E(L_jR_j)\).
This follows from \((L_jR_j)^p=0\) in the truncated algebra.
The resulting tensor is exactly the original physical operator
crossing, expressed in the regular local basis. There is no added
Gaussian to divide out.

Put \(A=\operatorname{diag}_{j\in J}(\tau^j)\) and
\(D_T=\operatorname{diag}_{i\in T}(\tau^{i+1})\). The inverse
Weyl parameter is

\[
C^{-1}=I+A\mathcal S,\qquad
\mathcal S=c^{-1}F_X^{-1}D_TF_V^{-\mathsf T}.
\tag{12}
\]

Every entry of \(\mathcal S\) has positive valuation, since all
inverse form matrices are regular and \(i+1\ge1\). Therefore both
\(C^{-1}=I+A\mathcal S\) and
\(A^{-1}C^{-1}A=I+\mathcal S A\) reduce to the identity at zero.
Their inverses are regular. These are exactly the weighted
hypotheses of the [Weyl theorem](identity-operator-weyl-reduction.md).

One complete polynomial numerator uses

\[
\Delta=c\det(F_X)\det(F_V),\qquad
\mathcal N=\Delta I+
A\operatorname{adj}(F_X)D_T\operatorname{adj}(F_V)^{\mathsf T}.
\tag{13}
\]

Clearing \(\Delta^{R(p-1)}\) in the Weyl matrix is a local zero-unit
operation. It is never applied to the global polynomial crossing.
The operator matrix has order \(p^R\); its auxiliary shifts
reconstruct all \(p^{2R}\) active zero orders.

At zero, the remaining identity spectators are precisely the indices
outside \(J\), and the auxiliary spectators are outside \(T\).
At infinity, only the original \(S\) is replaced by the filtered
operator. In particular, every direction in \(J\setminus S\)
remains an ordinary infinity spectator. Both full physical
paired-type partitions remain unchanged. General aggregate gauge
\(\eta\) adds \(P(1-\eta)\) times the occupation charge to each
active zero order, as in the existing operator construction.

The original endpoint repair concerns only the actual selected
operator. With common Laurent clearing and maximum resulting
coefficient degree \(D_F\), a sufficient ordinary degree bound is

\[
P\bigl(\deg G+2nr(p-1)\bigr)+D_F.
\tag{14}
\]

The count in (14) is the original \(r\), not the larger \(R\).

## 4. Terminal specialization and implementation

If \(S=\{N+d-r,\ldots,N+d-1\}\), with \(d\le r\le N\),
then \(S\) contains all new directions and is disjoint from \(L\).
Thus \(R=r+d\), and the retained indices are

\[
T=\{0,\ldots,r-d-1\}\cup\{N-d,\ldots,N+d-1\}.
\tag{15}
\]

The exact constructor is
[`surplus_identity_mode_geometry.py`](../tools/surplus_identity_mode_geometry.py).
`selected_geometry_data` and `selected_weyl_data` handle arbitrary
selected tuples; `geometry_data` and `weyl_data` provide the terminal
specialization. All coefficients are complete prime-field Laurent
polynomials, with no source-width truncation.

The [checker](../tools/check_surplus_identity_mode_geometry.py) and
[evidence](../evidence/surplus-identity-mode-geometry.json) cover
16 geometries in characteristics two, three, and five, with parent
lengths through 32 and surpluses through three. They verify 465 Lucas
block identities, 48 unit-determinant identities, 370 full-versus-
reduced Gram entries, 370 inverse-contraction identities, and 740
weighted-regularity conditions. Five complete geometry and Weyl
regressions agree with the independently implemented binary
one-surplus engine.

```bash
python3 tools/check_surplus_identity_mode_geometry.py --output /tmp/surplus-identity-geometry-replay.json
```

Literal original physical-cut controls and a completed pressure
catalogue currently come from that separate binary one-surplus
package. The new general-prime geometry files do not claim additional
full-tensor replays or search results.
