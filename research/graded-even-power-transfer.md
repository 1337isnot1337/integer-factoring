# Graded powers of even polynomial crossing seeds

## Result and scope

Even polynomial tensors admit a general sign-safe tensor-power transfer.
For each contiguous matrix cut, the globally graded power has the same
Smith orders as the ordinary Kronecker power of the **cut-ordered local
matrix**. Fixed physical types make the change between global boundary
orders a constant sign. This gives the usual coefficient-rank pressure
criterion without requiring a Gaussian seed.

The distinction between a local tensor and its cut-ordered matrix matters
in odd characteristic. For example, the signed residual tensor has ordinary
row rank 15 over \(\mathbb F_3\), but its contiguous Grassmann row matrix
has rank 16. The smaller rank cannot be used in this graded construction.

The theorem also gives a new admissible search family: arbitrary even
polynomial perturbations of the [all-field Pascal Gaussian](all-field-pascal-crossings.md)
which vanish at its crossing endpoint. Their typed endpoints and
effectivity are settled below. A pressure witness below \(1/\sqrt2\)
remains necessary; none is asserted here.

## 1. The local and global conventions

Let each physical space have a finite ordered homogeneous basis
\(\mathcal A\), with parity \(\epsilon(a)\in\{0,1\}\). Order the
four physical boundaries as

~~~math
(V,Y,U,X).                                               \tag{1}
~~~

Let \(T(w;v,y,u,x)\) be a polynomial tensor over a field \(K\).
Assume it is even: every nonzero entry has

~~~math
\epsilon(v)+\epsilon(y)+\epsilon(u)+\epsilon(x)=0\pmod2.    \tag{2}
~~~

For \(m\) copies, first multiply the tensors in copy order, using the
graded tensor product. Then group each physical boundary's \(m\)
letters together, keeping copy order inside that boundary. In coordinates,
the resulting tensor has entries

~~~math
T^{[m]}(w;\mathbf v,\mathbf y,\mathbf u,\mathbf x)
=(-1)^\sigma\prod_{i=1}^m T(w;v_i,y_i,u_i,x_i),
\quad
\sigma=\sum_{i<j}\sum_{a>b}\epsilon(a_i)\epsilon(b_j),      \tag{3}
~~~

where \(a,b\) range over the boundary positions in (1). The sign is
independent of \(w\), and can be computed directly. Formula (3) is the
definition; an unsigned local tensor power is not substituted for it.

The two cuts used for the crossing decomposition are

~~~math
(X,Y)\mid(V,U),\qquad (X,U)\mid(V,Y).                    \tag{4}
~~~

Before flattening locally, move the four graded boundary blocks into the
corresponding order in (4), retaining all Koszul signs. Call the resulting
local matrices \(M_0(w)\) and \(M_\infty(w)\). The same convention
defines the cut-ordered global matrices.

## 2. The even-power theorem

**Theorem.** For either cut, its cut-ordered global matrix differs from
\(M_e(w)^{\otimes m}\) only by row and column permutations and diagonal
matrices with entries \(\pm1\). Consequently its Smith orders at any
valuation place are all \(m\)-fold sums of the local Smith orders.
The same row and column equivalence applies separately to every polynomial
coefficient matrix.

**Proof.** Write a local copy in contiguous cut order as \(L_iR_i\).
Reordering
\(L_1R_1\cdots L_mR_m\) to
\(L_1\cdots L_mR_1\cdots R_m\) contributes the sign

~~~math
(-1)^{\sum_{i<j}|R_i||L_j|}.                            \tag{5}
~~~

On the support of the tensor product, evenness gives
\(|R_i|=|L_i|\pmod2\) for every \(i\). Thus (5) equals

~~~math
(-1)^{\sum_{i<j}|L_i||L_j|},                            \tag{6}
~~~

which depends only on the row label. Grouping the two physical boundary
blocks inside each cut side introduces a row sign and a column sign as
well. These changes are invertible constant matrices. This proves the
matrix equivalence, including coefficientwise equivalence.

Over a discrete valuation ring, tensoring Smith decompositions gives
all sums of their diagonal orders, with multiplicity. Zero diagonal
entries remain zero and can be omitted. This proves the profile claim.
The argument does not use a Gaussian representation. \(\square\)

For a fixed physical type, the number of occurrences of every letter
on each boundary is fixed. Its total parity is therefore fixed. Moving
the **global** physical boundary blocks from (1) into either order (4)
has one constant sign on the whole restricted tensor. Hence either
cut-ordered global matrix can be used to bound the ordinary flattening
of the same typed polynomial tensor. Taking a submatrix cannot increase
the coefficient ranks.

### Why the hypotheses cannot be omitted

If a local tensor has both odd and even entries, the sign in (5) need
not be a row or column factor. The checker gives a two-letter example
whose cut-ordered local matrix has rank one. Its ordinary square has
rank one, whereas its graded square has rank two over \(\mathbb F_3\).

Evenness also does not identify an arbitrary ungraded flattening with
the required cut-ordered local matrix. The change of local boundary
order can put a sign between variables on opposite sides of a cut,
which can change its rank.

## 3. Crossing endpoints and vanishing perturbations

First suppose the local endpoint itself is an ordinary crossing:

~~~math
T(1;v,y,u,x)=[v=x][y=u].                                \tag{7}
~~~

Restrict all four physical words to one common type. Let \(q\) be the
number of odd letters in that type. On the support of (7), the shuffle
exponent (3) simplifies to

~~~math
\sum_{i<j}\epsilon(x_i)\epsilon(x_j)
+\sum_{i<j}\epsilon(y_i)\epsilon(y_j)
=2\binom q2\pmod2.                                     \tag{8}
~~~

It is even. Therefore the globally graded, typed tensor has the exact
ordinary crossing endpoint. If the two crossing pairs have different
fixed types, the endpoint scalar is
\((-1)^{\binom{q_X}{2}+\binom{q_Y}{2}}\), a fixed unit which can be
removed. No input-dependent sign correction is needed.

The Pascal Gaussian supplies a second endpoint mechanism. Let
\(\mathcal G_S(w)\) be the even exterior element associated to the
one-particle kernel of [the all-field construction](all-field-pascal-crossings.md).
The local physical basis consists of subsets of the finite mode set
\(S\), with parity equal to the subset size modulo two. Its globally
graded power is the Gaussian associated to the block sum of the local
one-particle kernels. Regrouping copy and mode indices inside a physical
boundary changes only that physical basis by signs and permutations.
At a crossing, the matching input and output changes cancel.

Choose any finite collection of even exterior monomials \(e_I\), each
of positive even degree, and polynomial coefficients \(f_I(w)\) with
\(f_I(-1)=0\). Define the local seed

~~~math
\boxed{\quad
T_S(w)=\left(1+\sum_I f_I(w)e_I\right)\mathcal G_S(w).
\quad}                                                   \tag{9}
~~~

This is polynomial and even. At \(w=-1\), every perturbation vanishes,
so its graded powers have the same endpoint as the Gaussian. A common
physical word type fixes every mode occupation count. The endpoint
theorem for the Gaussian therefore proves an ordinary crossing for
every copy count and common type in (9).

Terms of degree four or higher permit non-Gaussian seeds. No claim that
every nonzero such term is non-Gaussian is needed. This family includes
arbitrary coefficients over any fixed finite field. Its admissibility
does not depend on characteristic-two sign cancellation, and it needs
no residual-code repair at the endpoint.

To put the endpoint at one, substitute \(w=-u^2\). Further integer
ramification permits rational height choices. Every parameter and every
local polynomial is fixed before the copy count grows.

## 4. Pressure from the actual local cut profiles

Let \(\eta_{e,D}(a)\) be a fixed monomial height for each endpoint
\(e\in\{0,\infty\}\), physical boundary \(D\), and physical
letter \(a\). Apply those diagonal monomials before computing the
Smith orders of the corresponding local matrix in (4). Write the
finite orders as \(\nu_{e,j}\), and put

~~~math
W_e(s)=\sum_j s^{-\nu_{e,j}},\qquad
G(a)=\sum_{e,D}\eta_{e,D}(a),\qquad s>1.                  \tag{10}
~~~

The even-power theorem makes the global profile sum exactly
\(W_e(s)^m\). For a common type \(\mathbf q=(q_a)\), the physical
dimension and the restored total height satisfy

~~~math
D_m=\frac{m!}{\prod_a(mq_a)!},\qquad
\log D_m=mH(\mathbf q)+O(\log(m+1)),
\qquad
g_0+g_\infty=m\sum_a q_aG(a).                            \tag{11}
~~~

Here the alphabet is fixed, \(\mathbf q\) is rational, and \(m\)
is divisible by its denominator. The local heights may depend on the
whole physical letter. They need not be additive over occupied modes.

For completeness, a matrix with Smith orders \(\nu_j\) has its
coefficient at degree \(r\) of rank at most
\(\sum_j\max(0,r-\nu_j+1)\). Expand the two regular unit factors
of a Smith decomposition to see this bound. For any fixed \(s>1\),
the bound is at most \(C_s s^r\sum_j s^{-\nu_j}\).
At zero, coefficient \(A_i\) of the ungauged typed tensor occurs
at degree \(i+g_0\). At infinity it occurs at degree
\(g_\infty-i\). Multiplying the two bounds cancels \(i\).
Consequently every coefficient has row/column rank product at most

~~~math
C_s^2\left(W_0(s)W_\infty(s)
             s^{\sum_a q_aG(a)}\right)^m.                \tag{12}
~~~

Thus, for \(R=4\rho\), the strict criterion is

~~~math
\log(W_0(s)W_\infty(s))
+(\log s)\sum_a q_aG(a)<R H(\mathbf q).                   \tag{13}
~~~

Optimizing over the type gives the equivalent open criterion

~~~math
\boxed{\quad W_0(s)W_\infty(s)
 <\left(\sum_{a\in\mathcal A}s^{-G(a)/R}\right)^R.
 \quad}                                                  \tag{14}
~~~

The maximizing real type is proportional to \(s^{-G(a)/R}\).
A strict inequality permits a nearby rational type. This is the full
physical-letter partition function; replacing it by independent mode
occupations would restrict the search.

For a fixed local polynomial tensor, the number of coefficients is
\(O(m)\). Every entry of the global tensor (3) is computable by polynomial
arithmetic and an explicit sign. Enumerating the restricted tables and
computing coefficient rank factorizations takes \(\exp(O(m))\) bit
operations over a fixed finite field. A positive-entropy rational type
and the strict margin in (13) permit embedding every binary label width,
with only logarithmic dimension and rounding losses. These are effective
crossing pieces at parameter \(\rho\), so the fixed-field evaluator
and [scalar recovery transfer](fixed-prime-transfer-source-audit.md) apply.
No inequality below the current threshold is proved by the transfer alone.

## 5. The signed residual control

The [signed residual note](signed-residual-codes.md) uses an exterior
operator with output order \((V,Y)\), input order \((X,U)\), and
ordinary entries \(R_{v+4y,x+4u}\). Passing to exterior coefficients
in boundary order (1) contributes physical input signs
\((-1)^{\binom{|x|}{2}+\binom{|u|}{2}}\), which do not change any
flattening rank. Moving next to the row cut in (4) contributes, on the
even support, the sign

~~~math
(-1)^{\epsilon(v)\epsilon(y)+\epsilon(x)}.                \tag{15}
~~~

The last factor is a row sign. The factor
\((-1)^{\epsilon(v)\epsilon(y)}\) connects opposite sides of the
cut and is not a physical diagonal correction. Over \(\mathbb F_3\),
the ordinary residual row matrix has rank 15 and the corrected matrix
has rank 16. The column rank remains 12. This illustrates why the
local matrix convention is part of the hypothesis of (14).

The residual endpoint codes themselves remain valid as stated in their
original ordinary convention. This note does not complete their separate
polynomial regularity obligations or assert a useful residual pressure.

## 6. Exact reproduction

```sh
python3 tools/check_graded_even_power.py \
  --output evidence/graded-even-power-transfer.json
```

The checker compares literal permutations of occupied graded variables
with the row/column sign factorization, coefficient by coefficient, for
two- and four-letter alphabets over \(\mathbb F_2,\mathbb F_3,\mathbb F_5\).
It checks common-type crossing endpoints with arbitrary even polynomial
perturbations vanishing at the endpoint. It also includes the odd-seed
negative control and the exact residual rank comparison in Section 5.
These finite controls support the sign bookkeeping; the uniform theorem
is the argument in Sections 1--4.
