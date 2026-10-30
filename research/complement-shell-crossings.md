# Explicit complement shells and their finite local tensor

## Result and scope

Fixed-length, fixed-energy subsets of the odd modes below a power of
two give explicit homogeneous dual labels with vacuum Gram equal to
the identity. No Gram-minor search or matrix inversion is required.
The resulting complete crossing has an exact reduction to a universal
four-state local tensor, one tensor for each complementary mode pair.
For labels of length \(L\), there are at most \(2L+1\) nonzero
shift pieces, independently of the mode cutoff.

The reduction proves the previously observed reciprocal-pencil
formula for complement-paired labels at every size. More general
shells use all four local states and need not obey conservation of
the individual mode-occupation bits.

The construction is effective at every binary label width. For
fixed \(L=3\), it also gives an explicit family with six possible
pieces and exactly \(M(M-4)/32\) labels. The measured complete ranks
do not meet the finite seed criterion. The all-size rank behavior of
the larger four-state family remains a separate question; this note
does not assert a factoring exponent above \(1/2\).

The operators are Anthropic's Appendix B operators at
[16a5ce41304fcdaea4e0be7c132e1c361dc26da5](https://github.com/anthropics/formal-math/tree/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth),
with the characteristic-two quotient described in
[modular-crossing.md](modular-crossing.md). The homogeneous
complete-profile formula is in
[homogeneous-central-selection.md](homogeneous-central-selection.md).
The following construction, ordinary proofs, and computations used
OpenAI assistant support. No new Lean formalization is claimed.

## 1. A shell has an identity vacuum Gram

Fix a power of two \(M\ge4\), a positive length \(L\), and an
energy \(j\le LM/2\). Let \(\mathcal A\) be any collection of
distinct \(L\)-element subsets of
\(\{1,3,\ldots,M-1\}\), each with sum \(j\). Define

\[
a^x=t^{A_x},\qquad
b^x=t^{\overline{A_x}},\qquad
\overline A=\{M-k:k\in A\},\qquad h=LM-j.
\tag{1}
\]

Use the same \(a\)-domain for the input and output state labels.
Thus state energy is \(E=j\), partner energy is \(h\), and
parameter energy is \(j\le h\).

**Theorem 1.** The raw vacuum pairing is

\[
[1]a^x[-h]b^v=[x=v].
\tag{2}
\]

**Proof.** The vacuum coefficient is the determinant of the
one-part contraction matrix \([n\mathbin{\&}m=1]\) over
\(\mathbf F_2\). Write \(n=2a+1\), \(k=2c+1\), and
\(m=M-k\). Since \(M/2\) is a power of two,

\[
[n\mathbin{\&}m=1]=[a\subseteq c].
\tag{3}
\]

In a nonzero determinant matching, each half-index in \(A_x\)
is at most its matched half-index in \(A_v\). Their sums both
equal \((j-L)/2\), so all these inequalities are equalities.
This forces \(A_x=A_v\). When the two sets coincide, distinctness
of the parts leaves exactly the identical matching. Its coefficient
is one. \(\square\)

The proof applies to every subcollection of the shell. It permits
arbitrary ordering of the labels and rectangular partner and
parameter energies. In particular, it does not require the sets
to be invariant under mode complementation.

## 2. A direct construction for every label width

For binary label width \(\lambda\ge1\), put

\[
K=2\lambda+16,\qquad L=K/2,
\]

and choose \(M\) to be the least power of two at least \(4K\).
Consider the \(L\)-subsets of the \(K\) low odd modes
\(\{1,3,\ldots,2K-1\}\). There are at most \(K^2+1\)
possible energies. A largest energy shell therefore contains at
least

\[
\frac{\binom K{K/2}}{K^2+1}
\ge\frac{2^K}{(K+1)(K^2+1)}
\ge2^\lambda
\tag{4}
\]

labels. For the last inequality, check \(\lambda=1\); as
\(\lambda\) increases by one, \(2^{K-\lambda}\) doubles while
\((K+1)(K^2+1)\), for \(K\ge18\), increases by a factor less
than two.

Choose the lowest-energy largest shell, and take its first
\(2^\lambda\) subsets in lexicographic order. Every chosen
mode is below \(M/2\), so \(j<LM/2<h\). Moreover
\(M<8K\), \(j\le K^2\), and \(h<LM=O(K^2)\).
Theorem 1 gives the exact duals immediately.

A subset-sum dynamic program counts each fixed-length shell and
unranks its first \(2^\lambda\) subsets. Its state indices are
bounded by \(O(K)\) modes, \(O(K)\) lengths, and \(O(K^2)\)
energies; its counts have \(O(K)\) bits. The complete construction,
including outputting the labels, takes \(2^{O(\lambda)}\) bit
operations. All label coefficients are bits.

This proves effective duality, not a favorable complete-piece
rank bound. The chosen constant in the energy bound is not used
to infer a new contraction rate.

## 3. Equality forces every intermediate energy to be divisible by \(M\)

The complete crossing is

\[
F_{x,y,u,v}(z)
=\sum_s [t^{A_u}]
 \sum_{i\subseteq LM}
 a^x[s-h+i]b^v[h-s-i]t^{A_y}\,z^s,
\tag{5}
\]

where \(i\subseteq LM\) means binary-digit containment.
The source identity and (2) give

\[
F_{x,y,u,v}(1)=[x=v][y=u].
\tag{6}
\]

Every source angle \(i\) is divisible by \(M\), because
\(M\) is a power of two and \(h+j=LM\).

Consider one term of the literal normal-order matching expansion
through an intermediate regular monomial \(T\). Suppose \(b^v\)
creates \(c\) parts and annihilates \(L-c\) parts. Since the input
length is \(L\), the intermediate length is \(2c\). The following
operator \(a^x\) must annihilate \(c\) parts and create \(L-c\).

For an intermediate odd part \(k\), let \(r(k)\) be its odd
residue in \(\{1,3,\ldots,M-1\}\), and define

\[
w_M(T)=\sum_{k\in T}r(k).
\]

This is a sum with multiplicity of residues; different actual
parts may have the same residue.

For the \(b^v\) action, write each parameter part as \(M-p\)
with \(p\in A_v\). A nonzero annihilation of an input part \(q\)
requires \(q\le p\). A nonzero creation of a part \(k\)
requires \(M-r(k)\le p\). These are consequences of the
binary containments in the one-part coefficient. Summing over
the determinant matching, and using
\(\sum A_v=\sum A_y=j\), gives

\[
j\ge
\sum_{\text{annihilated }q}q+cM
 -\sum_{\text{created }k}r(k),
\qquad w_M(T)\ge cM.
\tag{7}
\]

For the \(a^x\) action, a creation into an output part \(q\)
requires its parameter part \(p\le q\). An annihilation of
an intermediate part \(k\) requires \(p\le M-r(k)\).
There are \(c\) annihilations. Using
\(\sum A_x=\sum A_u=j\), the corresponding sum gives

\[
j\le
\sum_{\text{created }q}q+cM
 -\sum_{\text{annihilated }k}r(k),
\qquad w_M(T)\le cM.
\tag{8}
\]

Thus \(w_M(T)=cM\). More strongly, equality holds in every
individual inequality used in (7) and (8). This argument applies
to each nonzero matching before cancellation, not just to the
sum of matching coefficients.

Consequently the actual intermediate energy is divisible by
\(M\). It equals \(j+h-s-i=LM-s-i\), so \(s\) is divisible
by \(M\) as well.

**Theorem 2.** Every complete shift in (5) is a multiple of \(M\).
Since the full support lies in \([-2j,LM]\), the number of
possible nonzero pieces is at most

\[
\left\lfloor\frac{2j}{M}\right\rfloor+L+1
\le2L+1.
\tag{9}
\]

In particular, the residual \(F(z)-F(1)\) has Taylor order at
least \(M\) at \(z=1\). This follows from \(t=z^M\) and
\(t-1=(z-1)^M\) in characteristic two. It does not alone
bound a complete cut rank.

## 4. The exact universal four-state tensor

The equality in (7)--(8) gives more than sparse support.
Every contributing external creation becomes equality of modes,
and every contributing external annihilation becomes
complementarity of modes. The only infinite internal contraction
between parameters \(a\) and \(b=M-a\) is

\[
\sum_{q\ge0}z^{-(b+qM)}
=\frac{z^a}{1+z^M}.
\tag{10}
\]

It is legitimate to replace the original one-part coefficients
by these equality rules on the four external domains: every
original nonzero path already satisfies them, and all retained
equalities were valid original containments.

Put \(t=z^M\). Index the two modes in each complementary pair
\(\{r,M-r\}\) by bits \(0,1\). Each of the four labels has a
local occupation state in \(\{0,1,2,3\}\), the four subsets
of these two modes.

For local occupations \(x,y,u,v\), use four blocks of vertices
\(U,A,V,Y\), one vertex for each occupied mode in that block.
The \(V\) block is indexed by the uncomplemented label \(A_v\).
Define an alternating matrix with zero within-block entries
and the following cross entries:

| Vertex blocks | Local polynomial entry |
| --- | --- |
| \(U,A\) | \([u=a]\,(t+t^2)\) |
| \(A,Y\) | \([a\ne y]\,(1+t)\) |
| \(A,V\) | \([a=v]t\) |
| \(U,V\) | \([u\ne v]\) |
| \(V,Y\) | \([v=y]\) |
| \(U,Y\) | \([u=y]\) |

Here equality and inequality refer to the two local mode bits.
Let \(H_t(x,y,u,v)\) be its Pfaffian; an odd number of vertices
gives zero.

**Theorem 3.** If \(P=M/4\), and the local occupations of the
four labels at pair \(r\) are \(x_r,y_r,u_r,v_r\), then

\[
t^L F_{x,y,u,v}(z)
=\prod_{r=1}^{P}H_t(x_r,y_r,u_r,v_r),
\qquad t=z^M.
\tag{11}
\]

The right side is a polynomial of degree at most \(2L\).

**Proof.** Recompute the ordinary Wick Pfaffian after the
equality reduction. Divide the raw \(A_a\) vertex by its
factor \(z^a\). The product of these factors is \(z^j\),
which cancels the \(z^{-j}\) in the source angle profile.
The source's remaining multiplier is
\((1+z)^{LM}=(1+t)^L\). Assign one factor \(1+t\) to
each of the \(L\) vertices in the \(A\) block.
Before clearing denominators the six entries are
\(1+t,(1+t)/t,1,1,1,1\), with the equality conditions in
the table. In particular, (10) becomes the \(A,V\) entry
one. Multiplying every \(A\) vertex by \(t\) gives exactly
the table and multiplies the full Pfaffian by \(t^L\).

No edge joins different complementary mode pairs. The Pfaffian
therefore factors as the product in (11). Each \(A\) vertex
contributes degree at most two to its incident edge, and
every other edge has degree zero. There are \(L\) such
vertices in total, proving the degree bound. \(\square\)

The local tensor by itself is **not** a crossing on the unrestricted
four-state label set. At \(t=1\), proper local containment of
the \(A\) occupation in the \(V\) occupation can survive.
The equal global lengths force equality at every pair. All
crossing and rank claims in this note retain the shell restriction.

The reduction conserves the parity of the total occupation in
each complementary pair. It need not conserve each individual
mode bit. At \(M=16,L=4,j=32\), one nonzero complete entry has
abstract partner labels

\[
x=y=(1,3,13,15),\quad
u=(3,5,9,15),\quad v=(1,7,11,13),
\]

and shift support \(\{-64,32\}\). The XOR of their individual
mode masks is nonzero. This only excludes the direct
individual-mode conservation rule; it does not exclude changes
of basis or other representations of the tensor.

## 5. Complement-paired labels reduce to the reciprocal pencil

Suppose every label selects \(R\) complete complementary pairs.
Then \(L=2R\), \(j=h=RM\), and there are
\(\binom PR\) labels, where \(P=M/4\). The local states are
only \(0\) and \(3\). Writing those states as binary \(0,1\),
the nonzero values in (11) are:

| \((x,y,u,v)\) | \(H_t\) |
| --- | --- |
| \(0000\) | \(1\) |
| \(0011,0101,0110\) | \(1\) |
| \(1001\) | \(t^2\) |
| \(1010\) | \(t^2+t^4\) |
| \(1100\) | \(1+t^2\) |
| \(1111\) | \(1+t^2+t^4\) |

Let \(R_0,C_0,E_0\) be the three binary equality-pairing
tensors,

\[
R_0=[x=y][u=v],\quad C_0=[x=u][y=v],
\quad E_0=[x=v][y=u],\quad H_0=R_0+C_0+E_0.
\]

Put \(q(t)=1+t+t^2\), and let \(\mathcal P_R(t)\) denote

\[
\bigotimes_{r=1}^{P}(R_0+tH_0+t^2C_0)
\]

restricted to four weight-\(R\) binary labels. Then

\[
t^{2R}F(z)=q(t)^{2R-P}\mathcal P_R(t),
\qquad t=z^M,
\tag{12}
\]

as an identity of rational functions, with polynomial left side.
Equivalently, without negative powers of \(q\),

\[
q(t)^P t^{2R}F(z)=q(t)^{2R}\mathcal P_R(t).
\tag{13}
\]

To prove this entry by entry, odd local patterns vanish.
Among the six two-one patterns, equal weights of the four
labels force opposite orientations to have the same counts.
Let these counts be \(r,c,e\) for the \(R_0,C_0,E_0\)
orientations, and let \(b\) count four-one sites. Then
\(b+r+c+e=R\). The local table gives

\[
t^{2R}F(z)
=\left(q(t)^b t^{c+e}(1+t)^{r+c}\right)^2.
\tag{14}
\]

The ordinary reciprocal pencil gives the same expression
times \(q(t)^{P-2R}\), proving (12).
This establishes the full-entry reduction at all sizes.
In particular, this paired subfamily has residual Taylor
order at least \(2M\), since (14) is a polynomial in \(t^2\).
An additional rank estimate is still needed to turn such
flatness into a contraction bound.

## 6. An explicit fixed-length family with six possible pieces

For a power of two \(M\ge8\), set

\[
L=3,\qquad j=3M/2-1,\qquad h=3M/2+1.
\tag{15}
\]

The number of low-mode parameter labels is exactly

\[
d_M=\frac{M(M-4)}{32}.
\tag{16}
\]

For a direct count, put \(N=M/2\) and replace each odd mode
by its half-index in \(\{0,\ldots,N-1\}\). Three indices
must sum to \(3N/2-2\). The number of ordered bounded triples,
before imposing distinctness, is

\[
\binom{3N/2}{2}-3\binom{N/2}{2}=\frac{3N^2}{4}.
\]

For each of the three possible equality pairs there are
\(N/2\) triples. No triple has all three indices equal,
because \(3N/2-2\) is not divisible by three.
Subtract \(3N/2\), and divide by six. The result is
\(N(N-2)/8=M(M-4)/32\).

By (9), the only possible shifts are

\[
s/M\in\{-2,-1,0,1,2,3\}.
\tag{17}
\]

The local tensor computes each entry using at most \(4L\)
active mode pairs and polynomials of degree at most \(2L\).
For fixed \(L\), this cost is independent of \(M\), apart
from reading and comparing the mode indices. Complete rank
computation may still be expensive: a cut matrix has \(d^2\)
rows and columns. Formula (11) is a structure theorem, not
a claim that every such matrix has a cheap rank algorithm.

## 7. Complete finite ranks and reproduction

Let \(r_s,c_s\) be the two complete cut ranks and
\(S=\sum_s\sqrt{r_sc_s}\). The finite seed criterion here
requires \(S<d^{\sqrt2}\).

| \(M,L,j,h\) | Labels \(d\) | \(S\) | \(d^{\sqrt2}\) |
| --- | ---: | ---: | ---: |
| \(16,4,32,32\) | 8 | 87.984 | 18.931 |
| \(16,3,23,25\) | 6 | 60.120 | 12.603 |
| \(32,3,47,49\) | 28 | 835.961 | 111.324 |

For the length-three family, the exact profiles are:

| \(s/M\) | \(M=16:\ (r_s,c_s)\) | \(M=32:\ (r_s,c_s)\) |
| ---: | --- | --- |
| \(-2\) | \((3,10)\) | \((7,346)\) |
| \(-1\) | \((10,25)\) | \((95,422)\) |
| \(0\) | \((7,20)\) | \((194,392)\) |
| \(1\) | \((18,11)\) | \((524,99)\) |
| \(2\) | \((12,4)\) | \((378,8)\) |
| \(3\) | \((36,1)\) | \((784,1)\) |

The central length-four shell has ranks
\((1,64),(8,30),(6,6),(29,29),(6,6),(30,8),(64,1)\)
at \(s/M=-4,-2,-1,0,1,2,4\), respectively.

Reproduce the evidence with:

    python3 tools/check_complement_shell_crossings.py \
      --output evidence/complement-shell-crossings.json

The checker verifies exact vacuum Grams, complete summed
crossings, every coefficient rank in both cuts, every local
tensor entry against the full native Pfaffian for the three
recorded shells, and independent Python Pfaffian samples.
It also checks the direct all-width construction through
\(\lambda=5\) and \(5{,}664\) entries of the paired
reciprocal reduction. The evidence pins the Python checker,
native Pfaffian implementation, and local-tensor implementation
by SHA-256 hashes.

The all-size assertions follow from the proofs above. The finite
ranks neither establish an improved seed nor rule out a useful
large-\(M\) rank asymptotic for the four-state shell family.
