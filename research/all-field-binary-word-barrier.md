# A uniform all-field barrier for binary word and balanced-layer tensors

## The result

The [characteristic-two binary word barrier](modular-binary-word-barrier.md)
extends to every field, with a remainder uniform in the characteristic.
The field may even change with the word length.

**All-field theorem.** Let \(D=2^t\), and decompose the crossing on all
binary words as \(E_D=\sum_iF_i\). Suppose every piece is invariant under
simultaneously permuting the \(t\) positions in all four labels. Then,
over any field, some single piece has both flattening ranks at least

\[
 2^{\,3t/2-O((\log(t+1))^2)},                            \tag{1}
\]

where the implicit constant is absolute and independent of the field.
There is no restriction on the number of pieces or on the local tensors
used to generate the invariant space.

The same conclusion holds when the label domain is just the middle
Hamming layer, and a product form holds for independently permuted
position classes. These versions are proved in Section 6.

Thus these families require \(\rho\ge3/4\) in the existing crossing
transfer and cannot reach its half-rate threshold \(1/\sqrt2\).
The corresponding ceiling \(\log_2(4/3)<1/2\) is a restriction on this
method, not a factoring lower bound. No improved factoring exponent is
claimed. The earlier characteristic-zero full-word barrier has a
stronger constant in that particular characteristic; this note adds
uniformity across all fields and the balanced-layer application.

## 1. An integral eight-monomial polynomial

Keep the row arrangement \(A_{2x+y,\,2u+v}=A(x,y,u,v)\), and put

\[
 X=\operatorname{diag}(0,1,0,1),\qquad \ell(Y)=Y_{00}+Y_{11}.
\]

This time work over the integers initially and set

\[
 P(A)=\ell(A X\operatorname{adj}(A))-\det A.              \tag{2}
\]

An ordinary signed Laplace expansion gives

\[
\begin{split}
P(A)={}&(A_{00}A_{12}-A_{02}A_{10})
         (A_{21}A_{33}-A_{23}A_{31})\\
 &-(A_{01}A_{13}-A_{03}A_{11})
         (A_{20}A_{32}-A_{22}A_{30}).                    \tag{3}
\end{split}
\]

For a determinant monomial, the first term in (2) counts how many
selected entries lie in rows \(\{0,1\}\) and columns \(\{1,3\}\).
Subtracting the determinant leaves coefficients \(-1,0,1\) for
counts zero, one, and two. The two surviving Laplace terms give (3).

Every individual monomial in (3) is fixed by the reshuffling
\((x,y,u,v)\mapsto(x,u,y,v)\). Hence this is one common polynomial for
both cuts over the integers, and therefore over every field. On the
one-bit crossing \(e\), exactly its first monomial survives:

\[
 P(e)=1.                                                \tag{4}
\]

Reduction modulo two recovers exactly the polynomial in the earlier
note.

## 2. Odd characteristic and characteristic zero

Suppose \(2\ne0\) in the field, and extend to an algebraic closure.
The matrix

\[
 X_0=X-\tfrac12 I
\]

is traceless. Since \(\ell(I)=2\), equation (2) becomes

\[
 P(A)=\ell(A X_0\operatorname{adj}(A)).                  \tag{5}
\]

Thus \(P\) is a matrix coefficient of the polynomial representation

\[
 V=\det\otimes\mathfrak{sl}_4,\qquad
 A:Y\longmapsto\det(A)AYA^{-1}.                         \tag{6}
\]

It has dimension 15 and highest weight
\(\lambda_0=(2,1,1,0)\). It is simple whenever the characteristic does
not divide four, including characteristic zero. The direct torus and
elementary-conjugation proof in the characteristic-two note applies
without a quotient: a scalar traceless matrix is zero when \(4\ne0\).
Off-diagonal weight vectors generate all of \(\mathfrak{sl}_4\), and
a nonzero diagonal traceless vector produces an off-diagonal vector.

In characteristic two the corresponding simple representation has
dimension 14, and the determinant summand is also needed. That case
continues to use the earlier proof; it is not obtained by dividing
(5) by two.

## 3. Semisimple blocks below the characteristic

We use the following elementary form of the Schur-algebra connection.
Over an algebraically closed field, homogeneous degree-\(m\) polynomial
representations of \(\operatorname{GL}_4\) are modules for

\[
 S(4,m)=\operatorname{End}_{K[S_m]}((K^4)^{\otimes m}).
\]

Indeed, the pure powers \(A^{\otimes m}\) span this algebra by the
polynomial-duality argument in the earlier note. Every coefficient
of a degree-\(m\) polynomial representation is a linear functional on
that same space. This defines its action on \(S(4,m)\); compatibility
with multiplication follows first on pure powers and then by bilinearity.

If the characteristic is zero or is a prime \(p>m\), Maschke's theorem
makes \(K[S_m]\) semisimple. Decompose its tensor-space module into
simple isotypic components. The commutant is a direct sum of full
matrix algebras, so \(S(4,m)\) and all its modules are semisimple.

Consequently, if \(4a<p\), the representation \(V^{\otimes a}\) is
semisimple. The same assertion holds for every \(a\) in characteristic
zero. The polynomial \(P^a\) is a matrix coefficient of this tensor
power. Every constituent has a highest weight \(\mu\) satisfying

\[
 |\mu|=4a,\qquad
 \mu\preceq(2a,a,a,0),                                  \tag{7}
\]

where \(\preceq\) denotes dominance order. Every weight of \(V\) is
dominated by its highest weight \((2,1,1,0)\); tensoring adds weights.

In positive characteristic, entrywise Frobenius preserves simplicity
and semisimplicity over the algebraic closure. For every \(j\ge0\),
the polynomial

\[
 (P^a)^{p^j}=P^{ap^j}
\]

is therefore a sum of coefficients of simple representations with
highest weights

\[
 \lambda=p^j\mu\preceq(2k,k,k,0),\qquad k=ap^j.          \tag{8}
\]

When \(p=3\), use only \(a=1\). The base representation (6) is already
simple, so this case needs no assertion that degree four is below
the characteristic.

## 4. A rank bound uniform over all supporting weights

The composition-series rank lemma and the standard-tableau
multiplicity proof from the earlier note imply that a nonzero
coefficient action on a simple highest-weight module \(L_\lambda\)
forces rank at least \(f^\lambda\) in the word tensor space.
This is a filtration inequality, not a complete-reducibility assumption
for the large degree \(4k\).

For any partition \(\lambda\preceq(2k,k,k,0)\) of \(4k\),
log-convexity of integer factorials gives

\[
 \prod_{i=1}^4\lambda_i!\le(2k)!(k!)^2.                 \tag{9}
\]

An integer proof uses the elementary majorization moves: moving one
unit from a part \(b\) to a part \(a\ge b\) multiplies the factorial
product by \((a+1)/b\ge1\). A sequence of such moves reaches the
dominating partition. Equivalently, the reverse balancing moves
decrease the product.

The hook formula, followed by the central-binomial estimates, now gives
the same bound for every weight in (8):

\[
 f^\lambda\ge
 \frac{\binom{4k}{2k}\binom{2k}{k}}{(4k+4)^6}
 \ge\frac{2^{6k}}{(4k+4)^8}.                            \tag{10}
\]

This is why a mixture of several simple constituents is harmless.
A nonzero value of the one common functional \(\Phi_{P^k}\) forces
some supporting constituent to act nontrivially in each cut; all of
them satisfy (10). The constituent need not be the same in the two
cuts, while the detected tensor summand is the same.

The characteristic-two two-constituent proof also satisfies (10).
For its determinant type \((k,k,k,k)\), dominance (9) holds as well.

## 5. At most logarithmically many blocks, uniformly in the field

Let \(k=\lfloor t/4\rfloor\). In characteristic zero use a single
block of length \(4k\). In characteristic two use the binary
position classes in the earlier note.

For odd characteristic \(p\), write \(k=\sum_jd_jp^j\), with
\(0\le d_j<p\).

For \(p=3\), split each digit into at most two parts of size one.
For \(p\ge5\), put \(a_0=\lfloor(p-1)/4\rfloor\) and split each
nonzero digit \(d_j\) into positive parts \(a\le a_0\).
There are at most six parts per digit. This is immediate for
\(p=5,7\); for \(p\ge11\), \(a_0\ge2\) and
\((p-1)/a_0\le4+3/a_0<6\).

For each part take a position class of size \(4ap^j\) and apply the
functional \(\Phi_{P^{ap^j}}\). The hypotheses of Sections 3--4 hold
on that class. The total number of classes is at most

\[
 b\le6(1+\log_2t),                                     \tag{11}
\]

independently of \(p\). This remains true when \(p\) grows with \(t\).

Fix the at most three unused positions to a common word in all four
labels. The product of the class functionals detects the restricted
crossing. Applying the product-group composition-series argument
from the earlier note to the same detected summand gives

\[
 \operatorname{rank}_{\rm row}F_i,\quad
 \operatorname{rank}_{\rm col}F_i
 \ge \frac{2^{3T/2}}{\prod_{\nu=1}^{b}(m_\nu+4)^8},
 \qquad T=4\lfloor t/4\rfloor.                         \tag{12}
\]

In particular, both logarithmic ranks are at least

\[
 \frac32t-\frac92
 -48(1+\log_2t)\log_2(t+4),                             \tag{13}
\]

proving the uniform all-field theorem. The constants do not depend
on the characteristic or on the extension degree of a finite field.

## 6. Middle layers and independently permuted position classes

Every monomial of the integral polynomial (3) has exactly two ones
in each of the four local binary legs. Hence every monomial of
\(P^k\) has exactly \(2k\) ones in each leg, across \(4k\) positions.
The functional \(\Phi_{P^k}\) therefore queries only entries whose
four labels all lie in that middle Hamming layer.

Extend any tensor on the middle layer by zero to the full word
domain. Its two ranks are unchanged: the new flattenings merely add
zero rows and columns. Position invariance is preserved as well.
Although the crossing extends to a supported crossing rather than
the full-word crossing, every queried entry is inside the support.
The functional value remains one.

For a class of arbitrary length \(m\), retain
\(T=4\lfloor m/4\rfloor\) positions. On the remaining \(r=m-T<4\)
positions, fix a common pattern of weight \(\lfloor r/2\rfloor\).
The retained balanced labels then have total weight

\[
 T/2+\lfloor r/2\rfloor=\lfloor m/2\rfloor.
\]

Thus the proof applies to the middle layer for every class length,
with the same logarithmic remainder.

More generally, suppose labels are the product of middle layers in
classes of lengths \(m_1,\ldots,m_b\), and every summand is invariant
under independent permutations inside each class. Write
\(M=\sum_jm_j\). The product common functional forces one summand
to have both ranks at least

\[
 2^{\,3M/2-O(\sum_j(\log(m_j+4))^2)},                   \tag{14}
\]

again uniformly over all fields. A class of length less than four
can be fixed completely and absorbed in the displayed error.

The number of labels in this product domain satisfies

\[
 M-\sum_j\log_2(m_j+1)
 \le\log_2\prod_j\binom{m_j}{\lfloor m_j/2\rfloor}
 \le M.
\]

Therefore, whenever
\(\sum_j(\log(m_j+4))^2=o(M)\), equation (14) gives the same
\(\rho\ge3/4\) obstruction in terms of the actual label-bit length.
For example, \(O(R)\) classes of lengths at most \(O(R)\), with total
length \(\Theta(R^2)\), have only \(O(R\log^2R)=o(R^2)\) loss.

This is a direct restriction to a supported tensor domain. No claim
that an arbitrary low-rank matrix retains its rank after a Hamming-layer
restriction is needed.

## Verification

~~~
python3 tools/check_all_field_binary_word_barrier.py
python3 tools/check_all_field_binary_word_barrier.py \
  --output evidence/all-field-binary-word-barrier.json
~~~

The checker verifies the signed integral polynomial identity and its
balanced support, the odd-characteristic simple-adjoint coefficient
certificate, uniform digit splitting, and exact hook/factorial bounds
for dominated partitions. It passes 3,572 dominated-partition bounds,
57,358 block and middle-layer padding cases, and all 32 monomial
leg-balance checks. The 225 coefficient polynomials are independently
verified over characteristics 3, 5, 7, 11, and 101. It supplements the characteristic-two
checker and the written representation-theory argument; finite
calculations are not used to infer an asymptotic theorem.
