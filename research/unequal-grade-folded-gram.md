# Independent operator and state grades in the complete Gram construction

## Result

The complete folded-Gram construction permits different operator and state
exterior grades, with independent homogeneous energies. This removes the
restriction that the two label domains come from one fixed-length shell.

Let the operator grade be \(L\), its lower energy be \(j\), the state grade
be \(K\), and the state energy be \(E\). The exact merged Gram block is a
**noncentral weight slice** of the same elementary Frobenius modules used in
[the equal-grade construction](folded-gram-channels.md).
Both complete cuts remain computable with small binary matrices and
Gaussian-binomial coefficients.

The extension is proved below and checked against complete Pfaffian
polynomials in both cuts. The recorded dense-grade and independent-energy
screens do not improve the factoring exponent. This is an all-size channel
formula and a reproducible search result, not an exclusion of all
unequal-grade constructions.

## 1. Domains and the crossing endpoint

Let \(M\geq4\) be a power of two and \(N=M/2\). Write

\[
 \mathcal R_M=\{1,3,\ldots,M-1\},\qquad
 \mathcal S_{L,j}=\{a\subseteq\mathcal R_M:|a|=L,\ \sum a=j\}.
\]

Choose \(d\) distinct parameter labels \(a_x\) in \(\mathcal S_{L,j}\), and
\(d\) distinct state labels \(y\) in \(\mathcal S_{K,E}\). The choices may be
arbitrary subsets. Set

\[
 h=LM-j,\qquad
 b_v=\{M-r:r\in a_v\},\qquad j\leq h.
\]

The vacuum parameter pairing is the identity, by the same equal-sum
Boolean-containment argument as in
[complementary shells](complement-shell-crossings.md).
It does not refer to the state grade. The coefficient pairing on the
chosen monomial states is also the identity. Consequently the complete
source polynomial has crossing endpoint

\[
 F(1)_{x,y,u,v}=\delta_{xv}\delta_{yu}.
 \tag{1}
\]

In particular, the maximum automatically available common dimension is

\[
 d=\min\{|\mathcal S_{L,j}|,\ |\mathcal S_{K,E}|\}.
 \tag{2}
\]

The two flattenings throughout this note are
\((x,y)\mid(v,u)\) and \((x,u)\mid(v,y)\). Parameter labels occupy the
\(x,v\) legs; state labels occupy the \(y,u\) legs.

The local complementary-pair Pfaffian remains valid with these unequal
grades. Its six edge polynomials, in the vertex order \(U,A,V,Y\), are

\[
 UA=t+t^2,\quad AY=1+t,\quad AV=t,\quad
 UV=VY=UY=1,\qquad t=z^M,
 \tag{3}
\]

with the same equality or complementation conditions on mode indices.
At \(t=1\), every \(A\) vertex must match its corresponding \(V\) vertex.
The remaining \(U,Y\) vertices must match identically. This also proves
(1) directly for all \(L,K\).

The product of local tensors is \(t^L F\), of degree at most \(2L\).
Its coefficient at degree \(L+v\) is the source coefficient at shift
\(Mv\). The native Pfaffian normalization is instead \(z^{E+j}F\);
therefore a local polynomial monomial \(t^a\) has native exponent

\[
 E+j-LM+aM.
 \tag{4}
\]

Formula (4), including its nonnegative boundary, is checked independently.

## 2. Creator count and the exact residue-energy boundary

Put \(\Delta=K-L\). If the first operator uses \(c\) creators, it uses
\(L-c\) annihilators and leaves \(c+\Delta\) unchanged input modes.
Thus

\[
 \max(0,-\Delta)\leq c\leq L,\qquad
 \ell=2c+\Delta
 \tag{5}
\]

is the intermediate exterior degree.

Write a folded generator as \(v_{\eta,C,Q}\), with unchanged low modes
\(\eta\), created residues \(C\), and winding \(Q\). Now

\[
 |\eta|=c+\Delta,\quad |C|=c,\quad
 D=\eta\cap C,\quad S=\eta\triangle C.
\]

Let \(i=|D|\), \(p=c-i\), and index the component by the
**creator-only** \(p\)-subset \(B=C\setminus\eta\) of \(S\). Then

\[
 |S|=2p+\Delta,\quad
 \eta=D\cup(S\setminus B),\quad C=D\cup B.
 \tag{6}
\]

Its residue energy is

\[
 W(D,S)=2\sum_{r\in D}r+\sum_{r\in S}r.
\]

For an incoming \(L\)-part operator of energy \(H\) and a low-mode state
of energy \(E\), coefficient nonvanishing gives

\[
 H\leq W+(L-c)M-E.
 \tag{7}
\]

Indeed, each creator parameter is at most the created residue, and each
annihilator parameter is at most \(M\) minus the annihilated low mode.
Summing these inequalities gives (7).

For the outgoing operator of energy \(H'\), all surviving unchanged modes
are low. The corresponding inequality is

\[
 H'\leq E+cM-W.
 \tag{8}
\]

Since \(H+H'=LM\), the lower and upper bounds agree exactly.
They force

\[
 \boxed{
 W_R=E+j-LM+cM,\qquad W_C=E-j+cM .
 }
 \tag{9}
\]

Here \(R\) uses the exchanged order, whose incoming operator has energy
\(j\); \(C\) uses the original order, whose incoming operator has energy
\(h\).

Every individual inequality underlying (7) and (8) must be an equality.
Thus the same exact mode equalities and complementations hold as in the
equal-grade proof. Incoming and outgoing coordinates are adjoints on this
boundary. This proves the Gram factorization and retains all high-mode
windings.

## 3. The merged Johnson matrix

For two generator indices \(B,B'\), put \(e=|B\cap B'|\).
The common-support generating function is

\[
 \langle v_B,v_{B'}\rangle
 =[q^{Q-i}]\,(1+q)^{-i-e}.
 \tag{10}
\]

The \(i\) double residues must have positive height. A single residue can
have positive height on both sides exactly when it belongs to
\(B\cap B'\). This proves (10), also when \(\Delta\ne0\).

The low-residue determinant coordinates do not depend on the winding.
All source angles can therefore be merged before taking a rank. Since
\((1+z)^{LM}=(1+z^M)^L\), the result at shift \(Mv\) is

\[
 \Gamma_{L,K,c,i}(v)=
 \left[
 \binom{L-i-|B\cap B'|}{L+v-c-i}
 \right]_{B,B'\in\binom{[\,2(c-i)+\Delta\,]}{c-i}} .
 \tag{11}
\]

For the other cut, replace \(v\) by \(-v\). Equivalently, the row
polynomial in the component is

\[
 t^{c+i}(1+t)^{L-i-|B\cap B'|};
 \tag{12}
\]

the column polynomial is its degree-\(2L\) reversal. All exponents are
nonnegative: \(|B\cap B'|\leq c-i\leq L-i\).

The incoming-coordinate checker stores unchanged-only subsets
\(A=S\setminus B\), rather than \(B\). Its intersection conversion is

\[
 |B\cap B'|=|A\cap A'|-\Delta.
 \tag{13}
\]

This correction is necessary when the grades differ.

## 4. Exact small blocks in an arbitrary weight slice

Use the explicit modules \(B_m\) and binary changes of basis proved in
[the folded-Gram rank algorithm](folded-gram-channels.md). Their Frobenius
weights satisfy

\[
 \mathcal W_0=(),\quad
 \mathcal W_{2r+1}=(1,2\mathcal W_r),\quad
 \mathcal W_{2r+2}=(1,1,2\mathcal W_r).
\]

They sum to \(m\), and the full dimension of \(B_m\) is at most
\((m+1)^2\). These are explicit tensor products of two-dimensional modules.
The decomposition is weight preserving, so it applies to an arbitrary
torus weight, including odd exterior degrees.

In exterior degree \(\ell=2c+\Delta\), a summand has the form
\(\det^a B_m\), where

\[
 a=(\ell-m)/2,\qquad
 \tau=c-a=(m-\Delta)/2.
 \tag{14}
\]

The relevant basis consists of subsets \(I\) of the Frobenius weight list
whose total weight is \(\tau\). For two such subsets, write
\(w(I,J)\) for their weighted intersection.

The exact merged rank of this small block is

\[
 \boxed{
 \beta_{L,K,c,m}(v)=
 \operatorname{rank}_{\mathbb F_2}
 \left[
 \binom{L-a-w(I,J)}{L+v-c-a}
 \right]_{|I|_{\mathcal W}=|J|_{\mathcal W}=\tau}.
 }
 \tag{15}
\]

If the slice is empty or the lower binomial argument is negative, its
rank is zero. Every top in (15) is nonnegative, because
\(a+w(I,J)\leq a+\tau=c\leq L\).

For completeness, apply the generic matrix
\(g(T)=\left(\begin{smallmatrix}1&1\\1&T\end{smallmatrix}\right)\) to the
single residues. Substitute \(T=(1+t)^{-1}\).
The raw prefactor in (12) is \(t^{c+i}(1+t)^{L-i}\).
A determinant twist \(a-i\) contributes
\((t/(1+t))^{a-i}\). The remaining entry of \(B_m(g(T))\) is
\((1+t)^{-w(I,J)}\). The result is exactly

\[
 t^{c+a}(1+t)^{L-a-w(I,J)},
\]

whose coefficient gives (15). This proves that aggregation over different
double-residue counts is valid without losing their cancellations.

## 5. Multiplicities and both complete bounds

Define the low-mode subset polynomials

\[
 E_r(X)=\sum_{\substack{A\subseteq\mathcal R_M\\|A|=r}}X^{\sum A}
       =X^{r^2}{N\brack r}_{X^2}.
\]

At exterior degree \(\ell\) and residue energy \(W\), the torus character
of the residue module is

\[
 \chi_{\ell,W}(Y)=
 \sum_{r=0}^{\ell}[X^W]E_r(X)E_{\ell-r}(X)\,Y^r.
 \tag{16}
\]

The local module is \(1\oplus V\oplus\det V\). Decompose (16) in the
explicit basis of characters

\[
 \chi_{\ell,W}(Y)=
 \sum_{\substack{0\leq m\leq\ell\\m\equiv\ell\pmod2}}
 \nu_{N,\ell,m}(W)\,
 Y^{(\ell-m)/2}
 \prod_{w\in\mathcal W_m}(1+Y^w).
 \tag{17}
\]

The multiplicities are nonnegative integers. They are obtained by
triangular subtraction, just as in
[weighted exterior aggregation](folded-gram-aggregation.md).
The same proof works for odd \(\ell\).

Combining (9), (15), and (17) gives the complete bounds

\[
 \boxed{
 \begin{aligned}
 R_{Mv}&\leq
 \sum_{c=\max(0,-\Delta)}^L\ \sum_m
 \nu_{N,\,2c+\Delta,\,m}(E+j-LM+cM)\,
 \beta_{L,K,c,m}(v),\\
 C_{Mv}&\leq
 \sum_{c=\max(0,-\Delta)}^L\ \sum_m
 \nu_{N,\,2c+\Delta,\,m}(E-j+cM)\,
 \beta_{L,K,c,m}(-v).
 \end{aligned}}
 \tag{18}
\]

Each cut may also be capped by \(d^2\). If either cut bound is zero, the
whole coefficient tensor is zero. The physical source range reduces to

\[
 \max\{-L,-\lfloor(E+j)/M\rfloor\}
 \leq v\leq
 \min\{L,\lfloor(E+h)/M\rfloor\}.
 \tag{19}
\]

The common dimension and all multiplicities are evaluated exactly with
integer Gaussian-binomial coefficients. The implementation computes the
integer upper bound

\[
 S_{\rm int}=\sum_v
 \left\lceil\sqrt{\min(d^2,R_{Mv})\min(d^2,C_{Mv})}\right\rceil.
 \tag{20}
\]

It reports \(\rho=\log_2 S_{\rm int}/(2\log_2d)\). A sufficient seed needs
\(\rho<1/\sqrt2\). Floating-point logarithms are descriptive; every rank,
shell count, and summand in (20) is computed exactly.

## 6. Useful symmetries

Exchanging \((L,j)\) with \((K,E)\) preserves (18), whenever the new
operator energy also lies below its central shell. Replace \(c\) by
\(c+\Delta\). The exterior degree and determinant twist stay unchanged.
The new block slice is complementary to the old one; its weighted
intersection changes from \(w\) to \(w+\Delta\), leaving both binomial
arguments unchanged. The row residue energy stays unchanged. The column
residue energy becomes \(\ell M-W_C\), which has the same multiplicity
under low-mode complementation.

Reflecting the state energy \(E\) to \(KM-E\) exchanges the row and column
multiplicities. It therefore exchanges the two cuts and sends \(v\) to
\(-v\). These statements concern the effective profile after removing
pieces with a zero cut. Such zero pieces can have unequal unused ambient
upper bounds before the physical-range truncation (19).

Setting \(K=L\) and \(E=j\) recovers the previous weighted formula exactly.

## 7. Recorded checks and search

The independent checker verifies weighted multiplicities by enumerating
local residue occupancies, including odd grades. It compares the small
blocks with direct dense Johnson matrices in noncentral slices.

Six complete source controls include grades \(2/4\), \(4/2\), \(3/4\),
\(2/3\), \(3/2\), and equal grades with unequal energies. Across them,
2,401 full entry polynomials are reconstructed in **both** cuts. Independent
native Pfaffian samples verify normalization (4), and exact ranks of every
nonzero coefficient are bounded by (18). Equal-grade specializations and
the two symmetries are also checked.

The accompanying screen records 43 cases with central operator energy and
174 cases with independently varied operator and state energies. It covers
\(M=64,128,256\), dense state grades, and operator grades through \(64\).
For the joint-energy grids the best recorded cases are:

| \(M\) | \(L=K\) | \(j=E\) | Selected labels \(d\) | Upper-bound \(\rho\) |
|---:|---:|---:|---:|---:|
| 64 | 8 | 256 | 178870 | 0.9167394878 |
| 128 | 16 | 1024 | 2990752400778 | 0.8469728458 |
| 256 | 32 | 4096 | 3227693268242421225516534768 | 0.8084679176 |

None meets the sufficient threshold. The actual coefficient ranks may be
smaller than the ambient Gram bounds: incoming coordinate maps can omit
large parts of these modules. Computing those images, rather than scanning
more energies with the same upper bound, is the next structural target.

## Reproduce

~~~bash
python3 tools/check_unequal_grade_gram.py \
  --output evidence/unequal-grade-gram.json
python3 tools/screen_unequal_grade_gram.py \
  --output evidence/unequal-grade-screen.json
~~~

The source module is [unequal_grade_gram.py](../tools/unequal_grade_gram.py).
The check report is [unequal-grade-gram.json](../evidence/unequal-grade-gram.json),
and the bounded search is
[unequal-grade-screen.json](../evidence/unequal-grade-screen.json).
