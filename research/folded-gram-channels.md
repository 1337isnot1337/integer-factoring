# Complete folded-Gram channels and exact Johnson ranks

## Result

The [residue-folded space](residue-folded-channels.md) can be compressed
further in the complementary fixed-length construction.

First, the incoming and outgoing operators force the same residue-energy
boundary. On that boundary, one is the adjoint of the other. Each cut
therefore factors through a **Gram form** on folded vectors, rather than
through their whole span.

Second, the incoming coordinates do not depend on the winding. All source
angles can consequently be combined inside that Gram form before taking its
rank. The result is a finite Johnson matrix with entries

\[
 \binom{L-i-|A\cap B|}{L+v-k-i}\pmod2,
 \qquad |A|=|B|=p,\quad k=i+p.
 \tag{1}
\]

Its size is \(\binom{2p}{p}\), but its exact rank can be computed using small
blocks. This note gives an elementary, constructive decomposition proving
that algorithm. It uses only a fixed eight-dimensional change of basis and
Frobenius twists. No external modular-representation theorem is assumed.

The resulting formulas provide both complete channel bounds, including
dense lengths. They do not by themselves establish a better factoring
exponent. The companion checker independently reconstructs 5,408 complete
entry polynomials in both cuts.

## 1. Setup

Let \(M\) be a power of two, \(N=M/2\), and
\(\mathcal R_M=\{1,3,\ldots,M-1\}\). Parameter labels and state labels are
equally sized collections of \(L\)-subsets of \(\mathcal R_M\), all of weight
\(j\leq LM/2\). They may be different collections.

The partner parameter is \(\bar v=\{M-r:r\in v\}\), of weight
\(h=LM-j\). Put

\[
 \delta=2j-LM,\qquad t=z^M.
\]

The complete crossing polynomial multiplied by \(t^L\) has degree at most
\(2L\). Its coefficient of degree \(L+v\) is the original piece at shift
\(s=Mv\). See
[the exact complementary-shell construction](complement-shell-crossings.md)
for vacuum duality and the local Pfaffian formula.

Recall the folded generators

\[
 v_{\eta,C,Q}
 =t^\eta[z^Q]\prod_{r\in C}
      \left(\sum_{q\geq0}t_{r+qM}z^q\right),
 \qquad |\eta|=|C|=k.
\]

Write

\[
 D=\eta\cap C,\quad S=\eta\triangle C,\quad
 |D|=i,\quad |S|=2p,\quad k=i+p,
\]

and

\[
 W(D,S)=2\sum_{r\in D}r+\sum_{r\in S}r.
\]

For fixed \(D,S,Q\), generators are indexed by \(p\)-subsets \(A\subseteq S\),
with \(\eta=D\cup A\) and \(C=D\cup(S\setminus A)\).

## 2. The exact residue-energy boundary

For an \(L\)-part parameter of weight \(H\) acting on a low-mode input of
weight \(E\), a term with \(k\) creators has residue energy \(W\) satisfying

\[
 H\leq W+(L-k)M-E.
 \tag{2}
\]

Indeed, a creation coefficient can be nonzero only when its parameter part
is at most the created residue. An annihilation coefficient can be nonzero
only when its parameter part is at most \(M\) minus the annihilated low
mode. Sum these inequalities.

In the original operator order, the first parameter has weight \(h\) and
the input has weight \(j\). Equation (2) gives \(W\geq kM\).
The final \(j\)-weight parameter and \(j\)-weight output give the opposite
inequality \(W\leq kM\). Thus only

\[
 W=kM
 \tag{3}
\]

contributes to the column cut.

In the exchanged order, the first parameter and input both have weight
\(j\), while the outgoing parameter has weight \(h\). The corresponding
equalities give

\[
 W=kM+\delta.
 \tag{4}
\]

All inequalities used in either derivation must be equalities individually.
Consequently, every surviving creation or annihilation edge is an exact
mode equality or exact mode complementation.

This also proves the adjoint relation needed below. On the boundary (3),
the outgoing \(a_x\) coefficients are the transpose of the incoming
\(\bar x\) coefficients. On the boundary (4), the outgoing \(\bar v\)
coefficients are the transpose of the incoming \(a_v\) coefficients.
The unchanged low modes coincide on both sides. Thus each fixed
intermediate-energy contribution is a Gram matrix of incoming vectors.

## 3. The unmerged folded Gram matrix

Give the regular Fock basis its usual coefficient pairing. For two folded
generators in the same \(D,S,Q\) component, let

\[
 e=|A\cap B|.
\]

A double residue has one fixed low mode and one positive-height mode,
contributing \(z/(1+z)\) to the common-support series. A single residue can
have positive height in both generators precisely when it lies outside
\(A\cup B\). There are \(e\) such residues, each contributing \(1/(1+z)\).
Therefore

\[
 \langle v_A,v_B\rangle
 =[z^Q]\frac{z^i}{(1+z)^{i+e}}
 =[z^{Q-i}]\,(1+z)^{-i-e}.
 \tag{5}
\]

For \(Q>0\), this is

\[
 \binom{Q+e-1}{Q-i}\pmod2.
 \tag{6}
\]

For \(Q=0\), the coefficient interpretation (5) fixes the boundary
unambiguously. In particular, the matrix is all ones if \(i=0\), and zero
if \(i>0\).

Different residue-occupancy patterns \(D,S\) have disjoint Fock supports, so
their Gram forms form a direct sum.

The Gram rank may be smaller than the dimension of the folded span. For
example, with no double residues, positive even \(p\) and \(Q=1\), the span dimension
is \(2p-1\), while the Gram rank is \(2p-2\).

## 4. Combining every source angle before taking rank

The incoming coefficient of a folded generator is a low-residue determinant.
It does not depend on \(Q\). The angle sum can therefore use one common
generator index set.

The surviving source angles are \(Ma\), where \(a\) is a binary submask of
\(L\). Over \(\mathbb F_2\),

\[
 \sum_{a\subseteq L}z^a=(1+z)^L.
\]

For the column cut, the unshifted winding index is \(L-v-k\).
Combining (5) over all angles gives

\[
 [z^{L-v-k-i}]\,(1+z)^{L-i-e}
 =\binom{L-i-e}{L-v-k-i}.
 \tag{7}
\]

For the row cut, replace the angle submask \(a\) by \(L-a\), which is again a
submask. The result is

\[
 [z^{L+v-k-i}]\,(1+z)^{L-i-e}
 =\binom{L-i-e}{L+v-k-i}.
 \tag{8}
\]

The top \(L-i-e\) is nonnegative because \(e\leq p\) and \(i+p\leq L\).
Thus the complete matrix is a finite polynomial coefficient. No Laurent
tail or unspecified high-energy convention remains.

An equivalent useful statement is that a component contributes the row
polynomial

\[
 t^{k+i}(1+t)^{L-i-e}.
 \tag{9}
\]

The column polynomial is its degree-\(2L\) reversal:

\[
 t^{L-k+e}(1+t)^{L-i-e}.
 \tag{10}
\]

Both have nonnegative exponents and degree at most \(2L\).

## 5. Complete two-cut channel counts

Let

\[
 \gamma_{L,i,p}(v)=
 \operatorname{rank}_{\mathbb F_2}
 \left[
 \binom{L-i-|A\cap B|}{L+v-(i+p)-i}
 \right]_{A,B\in\binom{[2p]}p}.
 \tag{11}
\]

A negative lower binomial argument gives a zero matrix.
The exact rank algorithm in the next sections evaluates (11) without
constructing its exponentially large displayed matrix.

Define \(n_{M}(i,p,W)\) as the number of disjoint low-residue sets \(D,S\)
with \(|D|=i,\ |S|=2p,\ W(D,S)=W\). These counts are coefficients of

\[
 \prod_{r\in\mathcal R_M}(1+uX^{2r}+wX^r).
 \tag{12}
\]

Then valid complete channel bounds are

\[
 \boxed{
 \begin{aligned}
 R_{Mv}
 &\leq\sum_{i+p\leq L}
 n_M(i,p,(i+p)M+\delta)\,\gamma_{L,i,p}(v),\\
 C_{Mv}
 &\leq\sum_{i+p\leq L}
 n_M(i,p,(i+p)M)\,\gamma_{L,i,p}(-v).
 \end{aligned}}
 \tag{13}
\]

Every source angle has already been included in these formulas. Adding
separate angle-channel dimensions again would lose this cancellation.
The actual complete matrix can have smaller rank because the incoming
coordinates need not fill the entire folded component.

For a uniform bound independent of the selected energy shell, use

\[
 n_M(i,p,W)\leq\binom Ni\binom{N-i}{2p}.
\]

This gives

\[
 U_{N,L}(v)=
 \sum_{\substack{i+p\leq L\\i+2p\leq N}}
 \binom Ni\binom{N-i}{2p}\gamma_{L,i,p}(v),
 \quad
 R_{Mv}\leq U_{N,L}(v),\quad C_{Mv}\leq U_{N,L}(-v).
 \tag{14}
\]

One can also cap each cut by \(d^2\), where \(d\) is the selected label
dimension.

There are \(L(N-L)+1\) possible energies for an \(L\)-subset of the \(N\)
low odd modes. Hence some shell, chosen at or below central energy by
complement symmetry, has at least

\[
 \left\lceil\frac{\binom NL}{L(N-L)+1}\right\rceil
 \tag{15}
\]

labels. Equations (14)--(15) give an effective all-size optimization problem
that includes dense \(L\).

## 6. An elementary eight-dimensional decomposition

Let \(V\) be the natural two-dimensional module over \(\mathbb F_2\).
Write \(V^{[2]}\) for its Frobenius twist: a matrix acts by squaring each
entry. There is a constant, weight-preserving decomposition

\[
 V^{\otimes3}\cong
 \bigl(V\otimes V^{[2]}\bigr)
 \oplus(\det V\otimes V)
 \oplus(\det V\otimes V).
 \tag{16}
\]

Here is a direct proof. Multiplication maps \(V^{\otimes3}\) onto
\(\operatorname{Sym}^3V\). Its section sends the four monomials to

\[
 000,\quad
 001+010+100,\quad
 011+101+110,\quad
 111.
 \tag{17}
\]

The coefficients \(\binom3k\) are all one in \(\mathbb F_2\), so
multiplication composed with this section is the identity. The section
intertwines a generic two-by-two matrix, as can also be checked directly.
Multiplication \(v\otimes w^{[2]}\mapsto vw^2\) identifies
\(\operatorname{Sym}^3V\) with \(V\otimes V^{[2]}\).

The kernel is the direct sum of the two determinant copies obtained by
antisymmetrizing slots \(1,2\) and slots \(2,3\). Each copy has dimension two;
their four explicit vectors are independent. This proves (16).

The checker verifies all 64 entries of the generic polynomial intertwining
identity and the invertibility of its constant binary change of basis.
Thus (16) is not inferred merely from a character equality.

## 7. Constructive Frobenius blocks for every tensor power

Define explicit polynomial modules

\[
 B_0=\mathbf1,\qquad
 B_{2r+1}=V\otimes B_r^{[2]},\qquad
 B_{2r+2}=V\otimes V\otimes B_r^{[2]}.
 \tag{18}
\]

Each \(B_m\) is a tensor product of two-dimensional factors with
power-of-two Frobenius weights. Let their list be \(\mathcal W_m\):

\[
 \mathcal W_0=(),\quad
 \mathcal W_{2r+1}=(1,2\mathcal W_r),\quad
 \mathcal W_{2r+2}=(1,1,2\mathcal W_r).
 \tag{19}
\]

The weights sum to \(m\).

There is a recursive direct-sum decomposition of \(V\otimes B_m\).
For \(m=0\), it is \(B_1\). For odd \(m\), it is \(B_{m+1}\).
For \(m=2r+2\), apply (16):

\[
 V\otimes B_m
 \cong 2\,\det V\otimes B_{m-1}
 \oplus V\otimes(V\otimes B_r)^{[2]}.
 \tag{20}
\]

If the already constructed decomposition of \(V\otimes B_r\) contains
\((\det V)^d\otimes B_k\), the final term of (20) contains
\((\det V)^{2d}\otimes B_{2k+1}\), with the same multiplicity.
Since \(r<m\), this is an induction.

Iterating tensor multiplication by \(V\) yields a constructive decomposition

\[
 V^{\otimes n}\cong
 \bigoplus_{m+2d=n}
 \left((\det V)^d\otimes B_m\right)^{\oplus\mu_{n,m,d}}.
 \tag{21}
\]

All changes of basis are constant and preserve the number of occupied
natural-module coordinates, with the appropriate determinant shift.
The construction is valid for a generic matrix over every commutative
\(\mathbb F_2\)-algebra.

The factor count in \(B_m\) is at most logarithmic up to a constant factor.
More concretely,

\[
 \dim B_m\leq(m+1)^2.
 \tag{22}
\]

This follows immediately by induction from (18). Consequently every block
needed below has polynomial size.

## 8. Exact rank formula for the Johnson matrices

Let

\[
 H_p(T)_{A,B}=T^{|A\cap B|},
 \qquad A,B\in\binom{[2p]}p.
\]

This is the weight-\(p\) compression of

\[
 g(T)^{\otimes2p},
 \qquad
 g(T)=\begin{pmatrix}1&1\\1&T\end{pmatrix}.
 \tag{23}
\]

Because (21) preserves weights, it remains a direct-sum equivalence after
compressing both domain and codomain to weight \(p\).

For \(m+2d=2p\), index the relevant part of \(B_m\) by subsets
\(\alpha\) of its Frobenius factors with

\[
 \sum_{a\in\alpha}w_a=m/2.
\]

Write

\[
 w(\alpha,\beta)=\sum_{a\in\alpha\cap\beta}w_a.
\]

The compressed block entries are

\[
 (T+1)^d\,T^{w(\alpha,\beta)}.
 \tag{24}
\]

### Increasing intersection

For the matrix
\(\left[\binom{A+|X\cap Y|}{B}\right]\), substitute \(T=1+z\),
multiply by \((1+z)^A\), and take the coefficient of \(z^B\).
Its exact rank is

\[
 \sum_{m+2d=2p}\mu_{2p,m,d}\,
 \operatorname{rank}_{\mathbb F_2}
 \left[
 \binom{A+w(\alpha,\beta)}{B-d}
 \right]_{\alpha,\beta}.
 \tag{25}
\]

The unmerged Gram matrix uses \(A=Q-1,\ B=Q-i\).

### Decreasing intersection

For \(\left[\binom{A-|X\cap Y|}{B}\right]\), substitute
\(T=(1+z)^{-1}\). Since \(T+1=z/(1+z)\), the corresponding exact formula is

\[
 \boxed{
 \sum_{m+2d=2p}\mu_{2p,m,d}\,
 \operatorname{rank}_{\mathbb F_2}
 \left[
 \binom{A-d-w(\alpha,\beta)}{B-d}
 \right]_{\alpha,\beta}.}
 \tag{26}
\]

For (11), take

\[
 A=L-i,\qquad B=L+v-(i+p)-i.
\]

Every top in (26) is nonnegative in this application:
the central weighted subsets give
\(d+w(\alpha,\beta)\leq d+m/2=p\), and \(i+p\leq L\).

Equations (21), (25), and (26) constitute an exact rank algorithm polynomial
in \(p\) and the bit lengths of the binomial arguments. Binary containment
evaluates each entry; only polynomial-size binary matrices are eliminated.
No exponentially large Johnson matrix is constructed.

## 9. Verification and implementation

Files:

- [Reusable exact rank and channel module](../tools/folded_gram_ranks.py).
- [Independent checker](../tools/check_folded_gram_channels.py).
- [Evidence](../evidence/folded-gram-channels.json).

Run:

    python3 tools/check_folded_gram_channels.py \
      --output evidence/folded-gram-channels.json

The checker verifies the generic decomposition (16), tensor characters
through degree 64, direct dense Johnson ranks, and explicit high-mode Gram
matrices. It independently reconstructs every complete entry in the
2-, 6-, and 8-label shell cases from both factorizations (9)--(10).
It also checks (13) against every measured piece of the six full profiles
through 496 labels, and records dense-length examples.

Useful interfaces:

    folded_gram_rank(i, p, Q)
    johnson_binomial_rank(p, top_offset, order, intersection_sign)
    complete_component_rank(L, i, p, signed_shift)
    complete_gram_channels(M, L, j)
    unweighted_complete_gram_channels(N, L)

The final two functions return the complete channel bounds after all-angle
cancellation. They must not be followed by another angle sum.

For scale, when \(p=64\), the displayed Johnson matrix has dimension
\(\binom{128}{64}\). The algorithm replaces it by weighted central blocks;
the largest explicit matrix in that case has dimension 64.

This is a positive all-size channel construction and an exact rank
algorithm. The distinction between these upper bounds and the actual
complete-piece ranks remains essential when assessing the factoring score.
