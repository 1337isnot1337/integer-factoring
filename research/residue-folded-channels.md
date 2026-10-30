# Exact high-mode residue folding and two-cut channel dimensions

## Result and scope

For characteristic-two source operators whose parameter parts are below a
power of two \(M\), all creation and annihilation coefficients are periodic
in the odd mode modulo \(M\). This does **not** permit dropping high modes.
It does permit replacing each infinite residue class by an exact generating
vector with a specified total winding.

This note proves:

1. Every intermediate vector from an \(L\)-part low-mode state and an
   \(L\)-part low-mode parameter lies in an explicit folded space.
2. The dimension of that space has an exact finite binomial formula,
   including the characteristic-two rank loss at the energy boundary.
3. Both complete source factorizations can use those dimensions as channel
   counts.
4. Incoming and outgoing factors can be evaluated without enumerating high
   modes.
5. At sufficiently high energy, the folded dimensions are periodic rather
   than growing with the partition number.

The dimension formula is exact for the stated **ambient folded space**.
A particular operator or selected label domain can reach a proper subspace.
Thus the resulting channel counts are valid upper bounds for the complete
pieces, not assertions that their matrix ranks attain those counts.
The present bounds do not establish a better factoring exponent.

## 1. Periodic coefficients and distinct created residues

Work in the square-zero regular-mode quotient over \(\mathbb F_2\):
its state basis consists of products of distinct odd-mode variables \(t_k\).
Let \(M\geq2\) be a power of two, and let

\[
 \mathcal R_M=\{1,3,\ldots,M-1\},\qquad N=M/2.
\]

For an odd parameter part \(n=2a+1<M\) and odd mode \(k=2c+1\), the creation
and annihilation coefficients are

\[
 C_n(k)=[a\subseteq c],\qquad
 A_n(k)=[a\mathbin{\&}c=0].
\]

Here containment refers to binary bits. Adding \(M\) to \(k\) adds \(M/2\)
to \(c\), leaving every bit visible to \(a<M/2\) unchanged. Both coefficients
therefore depend only on \(k\bmod M\).

The normal-ordered operator coefficient is a determinant of creation and
annihilation columns. Two created modes with the same residue give identical
creation columns, so that determinant is zero.

Take an input state with \(L\) low modes and a parameter with \(L\) distinct
low parts. If \(k\) fields create and \(L-k\) annihilate, exactly \(k\) input
modes remain unchanged. Write their set as \(\eta\). The created residues
form a set \(C\) of size \(k\).

For each residue \(r\), introduce the formal vector series

\[
 T_r(z)=\sum_{q\geq0}t_{r+qM}z^q.
\]

The corresponding contribution at energy \(K\) is a scalar multiple of

\[
 v_{\eta,C,Q}
  =t^\eta[z^Q]\prod_{r\in C}T_r(z),
 \qquad
 Q=\frac{K-\sum_{r\in\eta}r-\sum_{r\in C}r}{M},
 \qquad |\eta|=|C|=k\leq L.
\]

Terms with nonintegral or negative \(Q\) are zero. Products are taken in the
square-zero quotient. In particular, if \(r\in\eta\cap C\), the creator at
residue \(r\) must have positive height. The forbidden zero-height term
vanishes automatically.

Define \(\mathcal V_{M,L}(K)\) as the span of all these vectors for
\(\eta,C\subseteq\mathcal R_M\), \(|\eta|=|C|\leq L\). Let
\(D_{M,L}(K)=\dim\mathcal V_{M,L}(K)\), with \(D_{M,L}(K)=0\) for \(K<0\).
Every incoming operator image described above lies in this space.

## 2. Splitting by residue occupancy

Put

\[
 D=\eta\cap C,\qquad S=\eta\mathbin{\triangle}C,
 \qquad |D|=i,\quad |S|=2p,\quad i+p\leq L.
\]

The letter \(D\) in this section denotes a residue set; the function
\(D_{M,L}\) denotes the dimension.
Every residue of \(D\) has two occupied modes: one low and one of positive
height. Every residue of \(S\) has one occupied mode. Other residues are
absent. Different pairs \((D,S)\) therefore have disjoint state supports.
Their dimensions add.

The low-residue energy and winding are

\[
 W(D,S)=2\sum_{r\in D}r+\sum_{r\in S}r,\qquad
 Q=\frac{K-W(D,S)}M.
\]

Fix \(D,S,Q\). Choosing \(\eta,C\) amounts to choosing a \(p\)-subset
\(A\subseteq S\):

\[
 \eta=D\cup A,\qquad C=D\cup(S\setminus A).
\]

For a state in this component, let \(T\subseteq S\) be the single residues
whose occupied modes have positive height. The coefficient of that state in
\(v_A\) is simply

\[
 [A\cap T=\varnothing].
\]

The possible positive-support sets are:

- If \(i>0\), all \(T\) with \(|T|\leq Q-i\), provided \(Q\geq i\).
  The positive double residues absorb any excess winding.
- If \(i=0,Q=0\), only \(T=\varnothing\).
- If \(i=0,Q>0\), all nonempty \(T\) with \(|T|\leq Q\).

Columns with \(|T|>p\) are zero because a \(p\)-subset \(A\) cannot avoid
them. Multiple height assignments with the same \(T\) give repeated columns
in this span calculation; they do not increase its dimension.

## 3. Exact component dimension

Write \(\ell(p)=2^{v_2(p)}\) for the lowest nonzero binary place of \(p>0\).
The rank of the component above is

\[
 \Phi(i,p,Q)=
 \begin{cases}
 0,&Q<i,\\
 \binom{2p}{\min(p,Q-i)},&i>0,\ Q\geq i,\\
 1,&i=0,\ Q=0,\\
 0,&i=p=0,\ Q>0,\\
 \binom{2p}{\min(p,Q)}-[Q<\ell(p)],
   &i=0,\ p>0,\ Q>0.
 \end{cases}
 \tag{1}
\]

### Cumulative inclusion rank

The matrix with rows indexed by \(p\)-subsets \(B\) of a \(2p\)-set and
columns indexed by subsets \(T\) of size at most \(r\), with entry
\([T\subseteq B]\), has rank

\[
 \binom{2p}{\min(p,r)}
 \tag{2}
\]

over every field.

For completeness, a more general statement has row subsets of size \(k\)
in an \(n\)-set and rank \(\binom{n}{\min(r,k,n-k)}\).
In the interior \(r<\min(k,n-k)\), restriction to the slice where the last
coordinate is zero gives the degree-\(\leq r\) space on \((n-1,k)\).
Its kernel contains the last coordinate times the degree-\(\leq r-1\)
space on \((n-1,k-1)\). Induction and Pascal's identity give the lower
bound \(\binom nr\).

Over the rationals, every lower-degree subset indicator is a scalar multiple
of a sum of degree-\(r\) indicators, so the rank is at most \(\binom nr\).
The matrix has integer entries; its rank after reduction to any
characteristic cannot exceed its rational rank. At the boundary, point
indicators of \(k\)-subsets give full rank when \(r=k\), and complementing
coordinates handles \(r=n-k\). This proves (2).

In the folded component, replace \(A\) by its complement \(B=S\setminus A\)
and use this inclusion matrix.

### Removing the constant column

The case \(i=0,Q>0\) omits the constant column \(T=\varnothing\).
On the \(p\)-subset slice, that column lies in the span of positive-degree
columns of degree at most \(r\) if and only if \(r\geq\ell(p)\).

For sufficiency, the sum of all degree-\(\ell(p)\) monomials is the constant
\(\binom p{\ell(p)}=1\) in \(\mathbb F_2\).

For necessity, let \(1\leq r<\ell=\ell(p)\). Choose a fixed
\((p+\ell-1)\)-element subset of the \(2p\)-set and sum evaluations over all
its \(p\)-subsets. The constant column has sum
\(\binom{p+\ell-1}{p}=1\pmod2\). Every positive monomial of degree
\(j\leq r\) has sum either zero or

\[
 \binom{p+\ell-1-j}{p-j}=0\pmod2.
\]

Both parity statements follow immediately from binary containment. The
constant column is therefore independent of those positive-degree columns.
This proves the subtraction term in (1).

## 4. Exact total dimension and an eventual period

Since residue-occupancy supports are disjoint, the total dimension is

\[
 D_{M,L}(K)=
 \sum_{\substack{D,S\subseteq\mathcal R_M,\ D\cap S=\varnothing\\
                  |S|\ {\rm even},\ |D|+|S|/2\leq L\\
                  K-W(D,S)\in M\mathbb Z_{\geq0}}}
 \Phi\left(|D|,\frac{|S|}{2},\frac{K-W(D,S)}M\right).
 \tag{3}
\]

This is a finite all-energy formula. It does not enumerate partitions of
\(K\). A direct evaluation has \(2^{O(M)}\) terms and uses only polynomial
work in \(\log(K+2)\) per term.

For simultaneous finite-energy evaluation, the counts of pairs \(D,S\) by
their sizes and residue energy are the coefficients of

\[
 \prod_{\substack{1\leq r<M\\r\ {\rm odd}}}
       (1+uX^{2r}+vX^r).
 \tag{4}
\]

The checker uses this finite dynamic program and then applies (1).

Assume \(1\leq L\leq N\). Every component with \(k=i+p>0\) has full rank
\(\binom{2p}{p}\) once \(Q\geq k\). Also

\[
 W(D,S)=\sum\eta+\sum C\leq2k(M-k).
\]

Consequently, for

\[
 K\geq K_0=3LM-2L^2,
\]

every contributing component has reached its full rank, and

\[
 D_{M,L}(K+M)=D_{M,L}(K).
 \tag{5}
\]

The stable dimensions, summed over one full residue period, satisfy

\[
 \sum_{K=K_0}^{K_0+M-1}D_{M,L}(K)
   =\sum_{k=1}^L\binom Nk^2.
 \tag{6}
\]

To see (6), the full component rank counts all \(p\)-subsets \(A\subseteq S\).
Across all \(D,S,A\), these are precisely the ordered pairs
\((\eta,C)\) of equally sized \(k\)-subsets. The empty pair contributes only
at energy zero and is absent from the stable period.

For every energy, whether stable or not,

\[
 D_{M,L}(K)\leq\sum_{k=0}^L\binom Nk^2.
 \tag{7}
\]

Thus high energy does not force an unbounded channel dimension once the
mode residue set and field length are fixed.

## 5. Both complete source cuts

Take low-mode parameter families \(a_x,b_v\), each of length \(L\), with
weights \(j\leq h\). Take low-mode input and output states \(y,u\), each of
length \(L\) and energy \(E\). In the characteristic-two source piece

\[
 A_s(x,y,u,v)=
 [t^u]\sum_{i\subseteq j+h}
       a_x[s-h+i]\,b_v[h-s-i]\,t^y,
 \tag{8}
\]

the notation \(i\subseteq j+h\) means that \(i\) is a binary submask, the
precise nonzero-angle condition.

In the displayed order, the intermediate energy is

\[
 K_C=E+h-s-i.
\]

The incoming vector is obtained by applying \(b_v\) to \(y\), so it lies in
\(\mathcal V_{M,L}(K_C)\). This gives the column cut
\((x,u)\mid(v,y)\).

The source exchange identity reverses the two operators. The intermediate
energy is then

\[
 K_R=E+s-h+i.
\]

Applying \(a_x\) to \(y\) gives a vector in
\(\mathcal V_{M,L}(K_R)\), producing the row cut
\((x,y)\mid(v,u)\).

Therefore valid complete channel counts are

\[
 \boxed{
 R_s=\sum_{i\subseteq j+h}D_{M,L}(E+s-h+i),\qquad
 C_s=\sum_{i\subseteq j+h}D_{M,L}(E+h-s-i).
 }
 \tag{9}
\]

These replace the corresponding regular-partition counts in the
[homogeneous construction](homogeneous-central-selection.md).
They retain every angle and use exact energy arguments, including negative
energies. An invertible parameter-basis normalization does not increase the
channel counts. A zero piece can be padded with a zero channel if the
factorization convention requires positive dimensions.

## 6. Exact factor evaluation without high-mode enumeration

### Incoming coordinates

For a parameter \(\theta\) and input label \(Y\), choose
\(\eta\subseteq Y\) and a distinct creator residue set \(C\), with
\(|\eta|=|C|\). The coefficient of \(v_{\eta,C,Q}\) is the determinant with
creation columns \(C\) and annihilation columns \(Y\setminus\eta\).
This determinant uses only low residues.

Within each \((D,S,Q)\) component, choose independent rows of the finite
incidence matrix from Section 2. Those rows select an explicit basis of
folded generators. Binary elimination expresses every incoming generator in
that basis. Matrix sizes are bounded by \(2^{O(L)}\), independently of \(K\).

### Outgoing coordinates

Fix a folded generator and an outgoing low-mode operator coefficient.
Let \(D=\eta\cap C\), and let

\[
 B=C\setminus\eta.
\]

Choose the subset \(T\subseteq B\) whose creator modes have positive height.
The positive-height creator residues are exactly \(D\cup T\), with
\(m=|D|+|T|\).

For a fixed \(T\), every assignment of positive heights summing to \(Q\)
has the same outgoing operator coefficient. The determinant depends only
on residues, while all unchanged output modes are low. The number of
assignments is

\[
 \binom{Q-1}{m-1}
\]

when \(m>0,Q\geq m\). For \(m=0\), there is one assignment if \(Q=0\) and
none otherwise.

Hence the outgoing factor is an exact finite sum over \(T\subseteq B\).
Evaluate one representative state for each \(T\): give every positive
residue height one and put the remaining \(Q-m\) units on a fixed
distinguished residue. Multiply by the binary parity of the displayed
binomial coefficient.

Lucas containment computes that parity directly from \(Q-1,m-1\). There
are at most \(2^L\) terms. This evaluation preserves the full high-mode
tail and its boundary; it is not a truncation.

## 7. Exact examples and verification

Some ambient dimensions at energy \(K=LM\) are:

| \(M,L,K\) | Folded dimension | Regular-partition dimension |
|---|---:|---:|
| \(8,2,16\) | 5 | 5 |
| \(16,3,48\) | 68 | 87 |
| \(16,4,64\) | 184 | 276 |
| \(32,3,96\) | 1678 | 2048 |
| \(32,4,128\) | 9557 | 11564 |
| \(32,8,256\) | 832040 | 2824974 |
| \(64,16,1024\) | 51016358091980 | 824068326214949 |

These values concern the folded ambient space, not an individual coefficient
matrix. The companion evidence also applies (9) to every measured piece of
the six complete complementary-shell profiles, through 496 labels, and
checks both rank inequalities.

Files:

- [Dimension and factor checker](../tools/check_residue_folded_channels.py).
- [Exact evidence](../evidence/residue-folded-channels.json).

Run:

    python3 tools/check_residue_folded_channels.py \
      --output evidence/residue-folded-channels.json

The checker compares:

- 280 component formulas with direct binary incidence-matrix ranks.
- 112 folded-space dimensions with explicit high-mode vector expansions.
- 612 complete incoming operator vectors with the folded span.
- 300 outgoing-factor evaluations with explicit winding expansions.
- 40 outgoing evaluations with the original integer normal-order recursion.
- Eventual periodicity and the stable residue-total identity in four cases.

The exact dimension formula, the incoming containment theorem, and the
outgoing composition formula are proved above for all sizes. These finite
checks independently test their implementation and boundary conventions.
