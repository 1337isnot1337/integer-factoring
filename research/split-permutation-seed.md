# A six-label permutation crossing with rank product 150

## Result and scope

This note gives an exact characteristic-two crossing on six labels. Its
61 nonzero coefficient pieces have row-rank times column-rank at most 150.
The exact integer inequality

\[
 150^5=75\,937\,500\,000<78\,364\,164\,096=6^{14}
\]

therefore gives the finite coefficient

\[
 \frac{\log 150}{4\log 6}<\frac7{10}<\frac1{\sqrt2}.                 \tag{1}
\]

This is a finite seed, not an improved asymptotic factoring exponent.
Tensoring every piece separately creates exponentially many pieces. The
relaxed transfer in [subexponential-piece-transfer.md](subexponential-piece-transfer.md)
does not allow that cost. Summing tensor coefficients requires new rank
bounds. Moreover, the direct extension using all permutations on an
increasing number of points is obstructed: one of its pieces has both
flattening ranks \(d^{2-o(1)}\), where \(d\) is the number of labels.

The construction is the identity-kernel specialization of the direct
Pfaffian construction in
[intersection-matching-crossing.md](intersection-matching-crossing.md).
All algebra needed here is given below. No new formal proof is assumed.

## 1. The crossing

Fix \(b\ge2\). A label is a permutation \(A\in S_b\), regarded as the
bipartite matching \(\{(i,A(i)):0\le i<b\}\). Thus \(d=b!\).
There are four label groups, ordered \(U,B,A,Y\). Within each group the
alternating matrix has zero entries. Entries between distinct groups are
zero unless their two particles are the same edge \(e=(i,j)\).

For a matching edge set

\[
 z_e=X_i,\qquad t_e=1+Y_j .
\]

On equal edges the six upper-triangular group blocks are

| Group pair | Entry |
| --- | --- |
| \(U,B\) | \(t_e\) |
| \(U,A\) | \(z_et_e\) |
| \(U,Y\) | \(z_e\) |
| \(B,A\) | \(z_e\) |
| \(B,Y\) | \(z_et_e\) |
| \(A,Y\) | \(t_e\) |

Let \(F(A,B,U,Y;X,Y)\) be the Pfaffian of this \(4b\)-particle
alternating matrix over \(\mathbb F_2[X,Y]\). At every coordinate equal
to one, \(t_e=0\) and \(z_e=1\). Only the \(U,Y\) and \(B,A\) blocks
remain, and their matching matrices are identities exactly when the
corresponding labels agree. Consequently

\[
 F(A,B,U,Y;\mathbf1,\mathbf1)
       =[A=B][U=Y].                                               \tag{2}
\]

The Pfaffian splits across modes because unequal edges never pair. A
mode appearing an odd number of times gives zero. A mode in all four
groups contributes

\[
 t_e^2+z_e^2t_e^2+z_e^2
       =1+Y_j^2+X_i^2Y_j^2.
\]

The product also has only even exponents when modes appear twice:
the permutation condition pairs their row and column contributions.
The formula in the next section proves this directly. Thus \(F\) is
multilinear in

\[
 P_i=X_i^2,\qquad Q_j=Y_j^2.
\]

Write \(T_{S,T}\) for the coefficient of \(P_SQ_T\), where \(S,T\subseteq
\{0,\ldots,b-1\}\). Equation (2) says
\(\sum_{S,T}T_{S,T}=[A=B][U=Y]\), with scale one and at most \(4^b\)
pieces. The row flattening has indices \((A,Y)\mid(B,U)\); the column
flattening has indices \((A,U)\mid(B,Y)\).

## 2. A compact entry formula

Identify every label with its edge set. If

\[
 A\mathbin{\triangle}Y\ne B\mathbin{\triangle}U,
\]

then \(F=0\). Otherwise let \(D=A\triangle Y\), with active row set
\(R(D)\) and active column set \(C(D)\). The nonempty components of
\(D\) are alternating cycles. On each such cycle, \(B\) chooses either
the \(A\)-edges or the \(Y\)-edges, and \(U\) chooses the other half.

Start with \(P_{R(D)}\). A cycle on which \(B=A\) contributes no
additional factor. A cycle on which \(B=Y\) contributes

\[
 \prod_{j\in C(\text{cycle})}(1+Q_j).
\]

For each inactive row \(i\notin R(D)\), multiply by

\[
 \begin{cases}
  1+Q_{A(i)},& A(i)\ne B(i),\\
  1+Q_{A(i)}+P_iQ_{A(i)},& A(i)=B(i).
 \end{cases}                                                       \tag{3}
\]

All these factors use disjoint sets of \(P\)- and \(Q\)-coordinates.
This proves multilinearity and gives an entry algorithm using a
polynomial number of operations in \(b\) per coefficient query.

To verify the formula, an active cycle with \(B=A\) pairs its modes
through the \(B,A\) and \(U,Y\) blocks. Their product is the cycle's
part of \(P_{R(D)}\). A cycle with \(B=Y\) uses the \(U,A\) and \(B,Y\)
blocks, and contributes the additional squared \(t\)-factors. On the
inactive complement, both \(A=Y\) and \(B=U\) are matchings on the same
remaining row and column sets. Modes in two groups contribute paired
\(t\)-factors; modes common to all four use the displayed four-mode
Pfaffian. Their product is exactly (3).

## 3. Exact row ranks for every size

Let \(\mathcal D_b\) be the set of distinct edge sets \(A\triangle Y\)
as \(A,Y\) range over \(S_b\). For \(D\in\mathcal D_b\), define
\(c_T(D)\) as the number of nonempty alternating cycles of \(D\) whose
column set meets \(T\). Use the falling factorial convention
\((a)_{\underline s}=a!/(a-s)!\) when \(0\le s\le a\), and zero when
\(s>a\).

The row rank of \(T_{S,T}\) is exactly

\[
 \boxed{\displaystyle
 r(S,T)=
 \sum_{\substack{D\in\mathcal D_b\\R(D)\subseteq S}}
   2^{c_T(D)}
   \bigl(|T\setminus C(D)|\bigr)_{\underline{|S\setminus R(D)|}} .
 }                                                               \tag{4}
\]

Indeed, the row matrix is block diagonal by \(D\), since a nonzero
entry requires \(A\triangle Y=B\triangle U\). Within a cycle, the
orientation matrix is

\[
 \begin{pmatrix}1&\prod_{j\in C}(1+Q_j)\\
 \prod_{j\in C}(1+Q_j)&1\end{pmatrix}.
\]

Its coefficient at the chosen \(Q\)-subset has rank one if that subset
is empty, and rank two otherwise. This gives the factor
\(2^{c_T(D)}\). The obligatory active factor \(P_{R(D)}\) gives the
condition \(R(D)\subseteq S\).

On the inactive complement, the coefficient in (3) is one precisely
when the two common matchings agree on \(S\setminus R(D)\), and the
images of these selected rows lie in \(T\setminus C(D)\). The matrix
is a disjoint union of nonempty all-ones blocks indexed by those
partial injections. Its rank is their number, the falling factorial
in (4). The cycle and inactive-complement blocks form a tensor
product. Summing their ranks proves (4).

The argument includes \(D=\varnothing\), for which there are no cycle
factors, and cases with empty \(S\) or \(T\).

## 4. Exact finite profiles

The independent checker compares the compact entry formula with a
direct recursive Pfaffian. It then builds every coefficient matrix
and computes both ranks by exact Gaussian elimination over
\(\mathbb F_2\).

| \(b\) | Labels \(d\) | Nonzero pieces | Nonzero degree orbits | Maximum rank product | A worst pair of ranks |
| --- | ---: | ---: | ---: | ---: | --- |
| 2 | 2 | 14 | 8 | 4 | \(2,2\) |
| 3 | 6 | 61 | 15 | 150 | \(15,10\) |
| 4 | 24 | 252 | 24 | 34,560 | \(90,384\) |

Ranks depend only on \((|S|,|T|)\), by the independent row and column
permutation actions. The six-label worst orbit is \((3,1)\); the
24-label worst orbit is \((4,0)\). Complete subset profiles, including
zero pieces, are saved in the evidence file.

For six labels, exact Gaussian elimination also supplies
factorizations attaining these ranks. They remain valid over every
commutative \(\mathbb F_2\)-algebra by scalar extension. Thus (1) is a
finite algebraic certificate, rather than a bound obtained by counting
declared intermediate states.

## 5. Why increasing the full permutation domain fails

Consider the piece \(T_{[b],\varnothing}\). Formula (4) becomes

\[
 r([b],\varnothing)
   =\#\{D\in\mathcal D_b:R(D)=[b]\}.                               \tag{5}
\]

There are \(d\,!b\) ordered pairs of edge-disjoint permutation
matchings, where \(!b\) is the number of derangements. Each full-support
union graph has \(2^{c(D)}\) ordered alternating decompositions, and
\(c(D)\le b/2\). For \(b\ge2\), the elementary derangement bound
\(!b\ge b!/3\) gives

\[
 r([b],\varnothing)\ge\frac{d^2}{3\,2^{b/2}}.                      \tag{6}
\]

There is also an identity submatrix in the column flattening. Restrict
its row pairs \((A,U)\) and column pairs \((B,Y)\) to edge-disjoint
permutation pairs. A nonzero entry must satisfy
\(A\sqcup U=B\sqcup Y\) and \(A\cap Y=\varnothing\). Since \(U\)
already contains every edge in the union outside \(A\), this forces
\(Y=U\) and \(B=A\). The restricted matrix is therefore the identity
on \(d\,!b\) pairs, so

\[
 c([b],\varnothing)\ge d\,!b\ge d^2/3.                            \tag{7}
\]

Equations (6)–(7), together with \(\log(b!)=\Theta(b\log b)\), imply
that both ranks are \(d^{2-o(1)}\). This particular family and
coefficient decomposition therefore has normalized worst rank
product tending to one, not to the desired value below
\(1/\sqrt2\).

This obstruction leaves open tensor powers of the fixed six-label
seed with a different coefficient grouping, or a different
deformation of its crossing. It does not give a lower bound on
arbitrary crossing constructions or factoring algorithms.

## 6. Curve interface for amplification experiments

One may substitute arbitrary polynomials \(P_i(z),Q_j(z)\in
\mathbb F_2[z]\) with value one at \(z=1\). Equation (2) still holds,
so the resulting univariate coefficient tensors form a scale-one
crossing. More generally, Laurent substitutions are valid after
multiplication by a common monomial; rational substitutions are valid
after clearing a common denominator that has value one at \(z=1\).

The implementation exposes exact coefficient matrices for these
specializations. A curve with a favorable tensor-coefficient rank
envelope could amplify the finite seed without enumerating all piece
assignments. The finite rank-product improvement alone does not prove
such an envelope.

## Reproduction

Run

\[
 \texttt{python tools/check\_split\_permutation\_seed.py
 --output evidence/split-permutation-seed.json}.
\]

Files:

- [tools/split_permutation_seed.py](../tools/split_permutation_seed.py):
  compact coefficient construction and curve specialization.
- [tools/check_split_permutation_seed.py](../tools/check_split_permutation_seed.py):
  an independent Pfaffian, exact rank calculations, and negative controls.
- [evidence/split-permutation-seed.json](../evidence/split-permutation-seed.json):
  all finite profiles and verification counts.

The recorded run checks 1,568 independent Pfaffian entries, all 336
subset row-rank formulas for \(b=2,3,4\), and 48 nonlinear curve
crossing identities. It also checks the full coefficient crossing
identity at all label quartets for these three sizes.
