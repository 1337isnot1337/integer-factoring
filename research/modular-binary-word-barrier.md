# A characteristic-two barrier for all binary position-invariant word tensors

## Result and scope

Let the labels be all binary words of length \(t\), so \(D=2^t\).
Work over any field of characteristic two. Suppose

\[
 E_D(x,y,u,v)=[x=v][y=u]=\sum_i F_i(x,y,u,v),
\]

and every summand is invariant under simultaneously permuting the \(t\)
positions in all four words. There is no restriction to pairing tensors,
Gaussian tensors, coefficient families, or a fixed number of pieces.
Arbitrary linear recombinations within the full invariant space are allowed.

**Theorem.** Some summand has both flattening ranks at least

\[
 \operatorname{rank}_{(x,y)\mid(u,v)}F_i,\quad
 \operatorname{rank}_{(x,u)\mid(y,v)}F_i
 \ \ge\ 2^{\,3t/2-O((\log(t+1))^2)}.                 \tag{1}
\]

The same summand satisfies both inequalities. Therefore its total binary
channel width is at least

\[
 3t-O((\log(t+1))^2).
\]

Under the existing width-to-factoring transfer, such a family has
\(\rho\ge3/4\), and hence rate parameter at most
\(\log_2(4/3)=0.41503749\ldots\). In particular it cannot attain the
half rate or improve it through that interface.

This is a restriction on these invariant crossing decompositions, not on
factoring algorithms. It does not cover symmetry-breaking pieces, arbitrary
restricted subsets of a larger alphabet, or general nonbinary word seeds.
It supplies no improved factoring exponent. The repository's factoring
exponent remains \(1/2\).

The argument is genuinely modular. It uses one explicit eight-monomial
polynomial, a simple 14-dimensional characteristic-two representation, and
a composition-series rank inequality. It does not import the semisimple
Schur-Weyl rank formula from the
[characteristic-zero word barrier](https://github.com/1337isnot1337/integer-factoring/blob/05d9c1b8aed4f8494ce230f0869f9048ef879474/research/word-symmetry-barrier.md).
This is an ordinary proof with an exact checker, without new Lean work.

## 1. One common polynomial for the two cuts

Write a local four-leg tensor as a \(4\)-by-\(4\) row matrix \(A\), with
\(A_{2x+y,\,2u+v}=A(x,y,u,v)\). All indices below run from zero.
Let

\[
 X=\operatorname{diag}(0,1,0,1),\qquad
 \ell(Y)=Y_{00}+Y_{11}.
\]

Define the homogeneous degree-four polynomial

\[
 P(A)=\det A+\ell\bigl(A X\operatorname{adj}(A)\bigr).       \tag{2}
\]

Its explicit expression in characteristic two is

\[
\begin{split}
P(A)={}&(A_{00}A_{12}+A_{02}A_{10})
         (A_{21}A_{33}+A_{23}A_{31})\\
 &+(A_{01}A_{13}+A_{03}A_{11})
         (A_{20}A_{32}+A_{22}A_{30}).                    \tag{3}
\end{split}
\]

To see (3), expand the determinant by permutations. In the trace term in
(2), a determinant monomial occurs once for each selected entry in rows
\(\{0,1\}\) and columns \(\{1,3\}\). Adding the determinant keeps exactly
the permutations having zero or two such entries. The zero case and
the two case give the two products in (3).

The column arrangement replaces the local entry index
\((x,y,u,v)\) by \((x,u,y,v)\). Each of the eight monomials in (3) is
individually unchanged by that replacement. Consequently \(P\) is
**the same polynomial in both arrangements**.

Let \(e(x,y,u,v)=[x=v][y=u]\) be the one-bit crossing. Its row matrix
has ones precisely at entries \(A_{00},A_{12},A_{21},A_{33}\).
Exactly one monomial of (3) survives, so

\[
 P(e)=1.                                                \tag{4}
\]

## 2. The two simple representations supporting this polynomial

Extend the field to an algebraic closure \(K\). This preserves all matrix
ranks and invariant-tensor identities, and permits polynomial arguments
over an infinite field.

Put

\[
 L=\mathfrak{sl}_4(K)/(K I).
\]

Since the characteristic is two, \(I\) has trace zero. Thus
\(\dim L=15-1=14\). The action

\[
 B:[Y]\longmapsto \det(B)[B Y B^{-1}]
 \quad(B\in\operatorname{GL}_4(K))                       \tag{5}
\]

is a polynomial representation of degree four: its coefficients are the
entries of \(B Y\operatorname{adj}(B)\), after passing to the quotient.
It has highest weight \((2,1,1,0)\).

Here is a direct irreducibility proof. Under the diagonal torus, the
off-diagonal classes \([E_{ij}]\) have distinct nonzero weights. If an
invariant subspace contains a nonzero off-diagonal weight component, it
contains some \([E_{ij}]\). Permutation conjugations then give all such
classes. Conjugating \(E_{ij}\) by \(I+sE_{ji}\) gives, in characteristic
two,

\[
 E_{ij}+s(E_{ii}+E_{jj})+s^2E_{ji}.
\]

Varying \(s\) gives every diagonal trace-zero class as well. If an
invariant subspace initially contains only diagonal classes, a nonzero
one has a noncentral diagonal representative. Conjugation by a suitable
\(I+sE_{ij}\) produces a nonzero off-diagonal class. Therefore every
nonzero invariant subspace is \(L\). Twisting by the determinant does
not affect irreducibility.

The highest-weight class is \([E_{03}]\): upper triangular unipotents
fix it, and a diagonal matrix scales it by
\(b_0^2b_1b_2\). This gives the asserted weight.

Both \(X\) and the functional \(\ell\) in Section 1 are legitimate on
this quotient: \(\operatorname{tr}X=0\) and \(\ell(I)=0\).
Thus the second term in (2) is a matrix coefficient of the simple
representation (5). The first term is the coefficient of the
one-dimensional determinant representation, of highest weight
\((1,1,1,1)\).

For \(k=2^j\), applying entrywise Frobenius to (5) preserves
irreducibility: Frobenius is surjective on \(\operatorname{GL}_4(K)\).
The two polynomials in

\[
 P(A)^k=(\det A)^k+
       \ell\bigl(A X\operatorname{adj}(A)\bigr)^k         \tag{6}
\]

are therefore coefficients of simple degree-\(4k\) representations with
highest weights

\[
 \lambda_{\rm a}=(2k,k,k,0),\qquad
 \lambda_{\rm d}=(k,k,k,k).                             \tag{7}
\]

No assertion about an arbitrary tensor product being semisimple is used.

## 3. A rank inequality that survives characteristic two

Let \(U=K^q\), and let \(W=U^{\otimes m}\). A position-invariant matrix
on \(W\) commutes with the permutation action of \(S_m\).
Such matrices are spanned by \(B^{\otimes m}\), with
\(B\in\operatorname{GL}_q(K)\), in every characteristic.

One direct proof is by polynomial duality. An invariant matrix has one
coordinate for each degree-\(m\) monomial in the \(q^2\) local matrix
entries. A linear functional vanishing on every \(B^{\otimes m}\)
would give a polynomial vanishing on all invertible \(B\). This is a
Zariski-dense open set over \(K\), so the polynomial and the functional
are zero.

Let \(T\) be such an invariant matrix. Every polynomial
\(\operatorname{GL}_q\)-submodule of \(W\) is preserved by \(T\),
because \(T\) is a linear combination of the group actions
\(B^{\otimes m}\).

**Composition-series rank lemma.** If a simple polynomial representation
\(L_\lambda\) occurs \(c_\lambda\) times in a composition series of
\(W\), and the induced action of \(T\) on \(L_\lambda\) is nonzero, then

\[
 \operatorname{rank}T\ge c_\lambda.                    \tag{8}
\]

Choose a basis adapted to a composition series. The resulting matrix is
block triangular. Its rank is at least the sum of the ranks of its
diagonal blocks: choose a nonzero maximal minor in each block and combine
these into a block-triangular minor. Every quotient isomorphic to
\(L_\lambda\) carries the same induced operator, up to change of basis.
Each of its \(c_\lambda\) copies therefore contributes at least one to
that lower bound.

This is an inequality, not an additive-rank identity for a filtration.
Off-diagonal extension maps can increase the total rank and do not
invalidate (8).

### An elementary multiplicity lower bound

For a partition \(\lambda\) of \(m\) with at most \(q\) rows, let
\(f^\lambda\) be its number of standard Young tableaux. Then

\[
 c_\lambda\ge f^\lambda.                               \tag{9}
\]

Here is a highest-weight-vector proof, which also specifies the integral
tensors involved. For every standard tableau \(S\), put \(e_r\) in each
position whose tableau entry lies in row \(r\), and antisymmetrize the
positions in every column. This gives a nonzero tensor \(v_S\).
Signs may all become one in characteristic two; the construction remains
valid because the basis labels in each column are distinct.

Each \(v_S\) is a highest-weight vector of weight \(\lambda\). Every
column is an alternating tensor on \(e_1,\ldots,e_h\), which upper
triangular unipotents fix. These vectors are independent over any field.
Indeed, order the tensor-word basis lexicographically. The original
row-label word is the least word in \(v_S\), with coefficient one:
the positions in each column increase downward. Distinct standard
tableaux have distinct row-label words. In any nonzero linear
combination, its least such word cannot cancel.

Let \(M\) be the group submodule generated by all these \(v_S\).
The lower-triangular/diagonal/upper-triangular open cell of
\(\operatorname{GL}_q\) shows that all weights of \(M\) are at most
\(\lambda\) in dominance order, and that its weight-\(\lambda\) space
is exactly the span of the \(v_S\). The open cell suffices because it is
dense and the representation is polynomial.

By the highest-weight classification of simple polynomial
\(\operatorname{GL}_q\)-modules, the simple factors of \(M\) have highest
weights \(\mu\le\lambda\). Only \(L_\lambda\) among them can contain
weight \(\lambda\), and that weight has dimension one in \(L_\lambda\).
Torus weight spaces are exact through a composition series. Therefore
\(M\) contains exactly \(f^\lambda\) factors \(L_\lambda\), proving (9)
for the ambient tensor power as well.

The highest-weight classification of simple polynomial modules is the
remaining standard representation-theory input in this proof. This use
does not presume complete reducibility. The explicit tableaux argument supplies the multiplicity
estimate otherwise often obtained from a Weyl filtration.

## 4. Polynomial functionals and a repeated rank in both cuts

For a homogeneous polynomial \(Q\) of degree \(m\) in the \(16\) local
entries, define a functional on the position-invariant \(m\)-word
tensors as follows. For each degree-\(m\) monomial \(a^\nu\), let
\(F_\nu\) be the common tensor value on that orbit, and put

\[
 \Phi_Q(F)=\sum_\nu [a^\nu]Q\,F_\nu.                    \tag{10}
\]

There are no multinomial coefficients or factorial denominators in this
definition. In particular,

\[
 \Phi_Q(A^{\otimes m})=Q(A).                            \tag{11}
\]

If \(Q\) is a linear combination of matrix coefficients of some simple
polynomial modules, nonzero \(\Phi_Q(T)\) implies that the induced
operator on at least one of those modules is nonzero. This follows by
writing \(T\) as a linear combination of pure matrix powers, as in
Section 3, and applying (11).

For \(m=4k\) with \(k\) a power of two, choose \(Q=P^k\). Equations
(6)--(9) imply

\[
 \Phi_{P^k}(F)\ne0\quad\Longrightarrow\quad
 \operatorname{rank}F\ge
 \min\{f^{(2k,k,k,0)},f^{(k,k,k,k)}\}.                  \tag{12}
\]

The polynomial \(P\), hence \(P^k\), is unchanged by the two
arrangements. The same functional therefore gives (12) for both
flattenings. Equation (4) also gives

\[
 \Phi_{P^k}(E_{2^{4k}})=1.                              \tag{13}
\]

Thus every exact invariant decomposition has a summand detected by this
one functional, and that summand has both rank bounds.

## 5. An explicit exponential lower bound

The hook formula in Vandermonde form, with a partition padded to four
parts, is

\[
 f^\lambda=
 \frac{(4k)!\prod_{i<j}(\lambda_i-\lambda_j+j-i)}
      {\prod_{i=1}^4(\lambda_i+4-i)!}.
\]

For either partition in (7), its numerator product is at least one.
Also

\[
 \prod_i(\lambda_i+4-i)!
 \le (4k+4)^6\prod_i\lambda_i!
 \le (4k+4)^6(2k)!(k!)^2.
\]

For the determinant partition the last step uses
\((k!)^2\le(2k)!\). Consequently both dimensions are at least

\[
\begin{split}
 B(k)&=\frac{\binom{4k}{2k}\binom{2k}{k}}{(4k+4)^6}\\
 &\ge\frac{2^{6k}}{(4k+1)(2k+1)(4k+4)^6}
 \ge\frac{2^{6k}}{(4k+4)^8}.                           \tag{14}
\end{split}
\]

The central binomial coefficients are maximal in their respective
rows, giving the elementary middle inequality. The estimate is valid
for every \(k\ge1\), including the small cases where the displayed
rational lower bound is less than one.

## 6. All word lengths from binary position classes

Given \(t\ge4\), keep
\(T=4\lfloor t/4\rfloor\) positions and fix the remaining at most three
positions to zero in all four labels. The crossing restricts to
\(E_{2^T}\), invariant summands remain invariant on the retained
positions, and neither flattening rank increases under restriction.

Write

\[
 \lfloor t/4\rfloor=\sum_{j\in J}2^j.
\]

Partition the \(T\) retained positions into classes of sizes
\(m_j=4\cdot2^j\). Every summand is invariant under independent
permutations within these classes. Apply the tensor product of the
functionals \(\Phi_{P^{2^j}}\), one per class. Its value on the crossing
is one, so one summand has nonzero value for this product functional.

To justify multiplying the rank bounds for that one summand, use the
product group \(\prod_j\operatorname{GL}_4\) and its tensor space
\(\bigotimes_j(K^4)^{\otimes m_j}\). Its invariant-matrix algebra is
spanned by products of local pure powers, by the same polynomial-duality
argument as before. The product functional is a sum of coefficients of
simple product-group modules, selecting one of the two simple types
(7) in each class. Such a product module has composition multiplicity
at least the product of the multiplicities in the classes: tensor their
composition series. Tensor products of simple modules for the
independent group factors are simple over \(K\).

The composition-series rank lemma therefore yields, in each cut,

\[
 \operatorname{rank}F_i\ge
 \prod_{j\in J}B(2^j)
 \ge
 \frac{2^{3T/2}}{\prod_{j\in J}(4\cdot2^j+4)^8}.         \tag{15}
\]

This argument uses one product functional on the complete summand.
It does not choose different summands for different classes or cuts.

Since \(|J|\le1+\log_2 t\), \(T\ge t-3\), and every factor in the
denominator is at most \((t+4)^8\), equation (15) implies the explicit
bound

\[
 \log_2\operatorname{rank}F_i
 \ge\frac32t-\frac92
       -8(1+\log_2t)\log_2(t+4),                       \tag{16}
\]

which proves (1). No effective estimate for a convergence modulus or
enumeration algorithm is needed: this is a rank obstruction.

## 7. Verification and limits

~~~
python3 tools/check_modular_binary_word_barrier.py
python3 tools/check_modular_binary_word_barrier.py \
  --output evidence/modular-binary-word-barrier.json
~~~

The checker verifies the eight-monomial identity, its individual
monomial invariance under reshuffling, its value on the crossing,
the 196 independent coefficient polynomials of the 14-dimensional
representation, exact small tableaux vectors, and the hook and
integer lower bounds. It passes 74 independent standard-tableau vectors,
444 exact upper-unipotent distribution checks, 256 independent hook
formula and integer-bound cases, seven Frobenius powers, and 16,381
binary position-class decompositions. Its finite calculations supplement the
all-length proof above; they do not establish it by extrapolation.

The larger exploratory calculation found an 18-dimensional common
coefficient space after combining the 196-dimensional simple-adjoint
coefficient space with the determinant line. The simple-adjoint space
alone has a 17-dimensional common intersection which annihilates the
crossing. The determinant line is essential in this witness. Only the
explicit polynomial (2)--(3), not the exploratory search, is needed
for the theorem.

The characteristic-two full-word family is therefore closed at a
parameter strictly worse than the known half-rate construction.
Finding better factoring bounds still requires a different crossing
construction or pieces that escape the stated position symmetry.
