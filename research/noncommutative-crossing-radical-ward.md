# Ward factors from commuting symplectic forms

This note proves the all-size factor observed in the symplectic radical
crossings, extends it to several overlapping forms, and gives an exact
dimension formula for the resulting homogeneous input codes.

Let \(\eta=q+r\). A code contained in the image of a product of \(k\)
commuting symplectic operators forces the **entire restricted Pfaffian
tensor** to be divisible by \(\eta^k\). This holds before choosing any
output projection. A mode-orientation selector then cancels both
unwanted vertices while preserving the crossing.

The result includes mixed local grades on disjoint blocks, as well as
overlapping binary forms. It also explains why repeatedly nesting the
same global-pivot construction adds only a scalar factor.

All statements are over a field of characteristic two.

## 1. The state operator and its Ward identity

Use four exterior-variable blocks \(U,A,V,Y\), with complementary mode
involution \(i\mapsto\bar i\). The universal Pfaffian generating kernel
has alternating quadratic form
\[
 \mathcal B=
 p\sum_iU_iA_i+
 q\sum_iA_iY_{\bar i}+
 r\sum_iA_iV_i+
 \sum_iU_iV_{\bar i}+
 \sum_iV_iY_i+
 \sum_iU_iY_i.
 \tag{1}
\]
Its Pfaffian exponential means the perfect-matching generating
polynomial; no division by factorials is used.

Regard this kernel as a state operator
\(\mathscr M(A,U)\), with matrix rows indexed by \(Y\) and columns by
\(V\). Its entries are polynomials in \(p,q,r\) with coefficients in the
external exterior algebra on \(A,U\). Extracting external labels
\(x,u\) and state labels \(y,v\) recovers the tensor entry.

Let \(a_i\) create state mode \(i\), and let \(b_i\) contract it. Define
\[
 c_i=a_i+b_{\bar i}.
\]
In characteristic two,
\[
 c_i^2=0,\qquad c_ic_j=c_jc_i,\qquad c_i^{\mathsf T}=c_{\bar i}.
 \tag{2}
\]
Left multiplication by \(c_i\) acts on the kernel as
\(Y_i+\partial_{Y_{\bar i}}\); right multiplication acts as
\(\partial_{V_i}+V_{\bar i}\).

The standard first-derivative Pfaffian identity gives
\[
 \partial_{Y_{\bar i}}\mathcal F
  =(qA_i+U_{\bar i}+V_{\bar i})\mathcal F,
 \qquad
 \partial_{V_i}\mathcal F
  =(rA_i+U_{\bar i}+Y_i)\mathcal F.
\]
Their sum proves the exact operator identity
\[
 \boxed{\quad [c_i,\mathscr M]=\eta A_i\mathscr M,\qquad
 \eta=q+r.\quad}
 \tag{3}
\]
Here the commutator is a sum in characteristic two. Identity (3) holds
before any exterior-grade or label-code restriction.

## 2. Commuting square-zero forms give powers of the factor

Choose one representative \(i\) from each complementary mode pair.
For a weight vector \(w=(w_1,\ldots,w_P)\), put
\[
 Q_w=\sum_{i=1}^P w_i c_i c_{\bar i}.
 \tag{4}
\]
Equation (2) shows that all such operators commute, are self-adjoint,
and satisfy
\[
 Q_w^2=0.
 \tag{5}
\]

From (3),
\[
 [Q_w,\mathscr M]
 =\eta\sum_iw_i
  \bigl(A_{\bar i}c_i\mathscr M+A_i\mathscr M c_{\bar i}\bigr).
 \tag{6}
\]
Every \(Q_w\) commutes with every \(c_i\) and with the external variables.
Thus each additional commutator with another \(Q\) supplies another
factor of \(\eta\):
\[
 \operatorname{ad}_{Q_1}\cdots
 \operatorname{ad}_{Q_k}(\mathscr M)
 \quad\text{is divisible by }\eta^k.
 \tag{7}
\]
One can also write the exact formula
\[
 \operatorname{ad}_{Q_1}\cdots\operatorname{ad}_{Q_k}(\mathscr M)
 =\eta^k\mathscr M
   \prod_{a=1}^k
   \left(
    \sum_iw_{ai}(A_{\bar i}c_i+A_i c_{\bar i})
    +\eta\sum_iw_{ai}A_iA_{\bar i}
   \right).
 \tag{8}
\]

**Theorem 1.** Let \(Q=Q_1\cdots Q_k\), with \(k\ge1\). Suppose every
input vector belongs to \(\operatorname{im}Q\). Contract both input legs
\(y,v\) with these vectors. The resulting tensor is divisible by
\((q+r)^k\), for arbitrary output functionals.

**Proof.** Write the input matrix as \(I=QW\). Since the \(Q_a\) are
commuting and self-adjoint,
\[
 I^{\mathsf T}\mathscr M I
 =W^{\mathsf T}Q\mathscr M QW.
\]
Repeatedly commute the left \(Q_a\) through \(\mathscr M\). Every term
except the full iterated commutator contains a repeated \(Q_a\), and
vanishes by \(Q_a^2=0\). Hence
\[
 Q\mathscr M Q
 =\bigl(\operatorname{ad}_{Q_1}\cdots
        \operatorname{ad}_{Q_k}\mathscr M\bigr)Q.
 \tag{9}
\]
Equation (7) proves the claim. \(\square\)

The proof does not require an output parity condition or a coordinate
projection. Homogeneous grade is needed only when using the fixed-grade
crossing endpoint later.

Also,
\[
 Q^{\mathsf T}Q=Q^2=0,
\]
so every product-image input code is totally isotropic:
\[
 I^{\mathsf T}I=0.
 \tag{10}
\]

## 3. Recover the ordinary radical factor

On one pair, \(c_ic_{\bar i}\) is zero on the two singly occupied states.
On the empty/fully occupied two-dimensional space, its matrix is
\[
 \begin{pmatrix}1&1\\1&1\end{pmatrix}.
\]
For a set \(J\) of pairs, let \(Q_J\) use its indicator weights. Then
\[
 Q_J=\Omega_J+\Lambda_J+|J|+N_J,
 \tag{11}
\]
where \(N_J\) is physical degree parity on that set.

On parity \(N_J\equiv|J|\), put \(D_J=\Omega_J+\Lambda_J\). There
\[
 Q_J=D_J,\qquad D_J^2=0,\qquad
 \Omega_J\Lambda_J=\Omega_JD_J=D_J\Omega_J.
 \tag{12}
\]
Thus the local Gram radical
\(\operatorname{im}(\Omega_J\Lambda_J)\) is contained in
\(\operatorname{im}Q_J\). Theorem 1 with \(k=1\) proves the all-size
\(q=r\) vanishing for the radical codes in
[noncommutative-crossing-symplectic-radicals.md](noncommutative-crossing-symplectic-radicals.md).

In particular, the vanishing persists if all ambient output labels are
retained. It is a property of the two input subspaces, not an artifact
of their chosen coordinate pivots.

## 4. Exact homogeneous dimensions for overlapping forms

Let \(A=(w_{ai})\) be a \(k\times P\) weight matrix of full row rank.
Define
\[
 \mathcal U_{A,L}
  =\operatorname{im}(Q_{w_1}\cdots Q_{w_k})\cap\Lambda^L V.
 \tag{13}
\]

Fix a syndrome \(\sigma\) of singly occupied pairs. Write
\[
 s=|\sigma|,\qquad E=\{i:\sigma_i=0\},\qquad
 t=|E|,\qquad h=(L-s)/2,
\]
and let \(A_E\) be the restriction of \(A\) to the even positions.

**Theorem 2.** For each fixed orientation of the single modes, the
dimension of (13) in this syndrome is
\[
 \boxed{
 \begin{cases}
 \displaystyle\binom{t-2k}{h-k},
    &\operatorname{rank}A_E=k\ \text{and}\ A_EA_E^{\mathsf T}=0,\\[6pt]
 0,&\text{otherwise}.
 \end{cases}}
 \tag{14}
\]

**Proof: the full product image.** On an even pair, write \(X_i\) for
the flip exchanging empty and full occupancy. Then
\[
 c_ic_{\bar i}=1+X_i.
\]
On all \(t\) even pairs, the \(2^t\)-dimensional space is the regular
module for
\[
 F[q_1,\ldots,q_t]/(q_1^2,\ldots,q_t^2),
 \qquad q_i=1+X_i.
 \tag{15}
\]
For example, the products \(q_S\) applied to the all-full monomial give
an independent triangular basis.

If the \(k\) rows of \(A_E\) are dependent, their product as linear
forms in the \(q_i\) is zero. If they are independent, extend them to an
invertible linear change of the \(q_i\). In the new coordinates,
\[
 \operatorname{im}(Q_1\cdots Q_k)
 =\bigcap_{a=1}^k\ker Q_a,
 \qquad
 \dim\operatorname{im}(Q_1\cdots Q_k)=2^{t-k}.
 \tag{16}
\]
There is no additional multiplicity factor in (16).

**Proof: intersect with one homogeneous grade.** On the syndrome,
\[
 Q_a=\Omega_a+\Lambda_a+\kappa_a I,\qquad
 \kappa_a=\sum_{i\in E}w_{ai}.
\]
For a homogeneous vector, the three terms have degrees \(h+1,h-1,h\)
in the even-pair variables. Thus \(Q_av=0\) is equivalent to
\[
 \Omega_av=0,\qquad\Lambda_av=0,\qquad\kappa_av=0.
 \tag{17}
\]
Moreover,
\[
 \Lambda_a\Omega_b+\Omega_b\Lambda_a
 =\sum_{i\in E}w_{ai}w_{bi}.
 \tag{18}
\]
If any entry of \(A_EA_E^{\mathsf T}\) is nonzero, (17)--(18) force
\(v=0\).

Conversely, suppose \(A_EA_E^{\mathsf T}=0\) and the rows are independent.
Its diagonal entries imply \(\kappa_a^2=0\), hence \(\kappa_a=0\).
Choose linear coordinates \(e_1,\ldots,e_k,f_1,\ldots,f_k\), together
with \(t-2k\) residual variables, so that \(\Omega_a\) is multiplication
by \(e_a\) and \(\Lambda_a=\partial_{f_a}\). This is possible because
the independent creation vectors lie in the common kernel of the
annihilation covectors; choose dual \(f_a\), then a residual complement.

The common kernel in (17) consists exactly of monomials containing all
\(e_a\), no \(f_a\), and \(h-k\) residual variables. Its dimension is
\(\binom{t-2k}{h-k}\). Together with (16), this proves (14).
\(\square\)

For binary \(A\), the zero-Gram condition is the explicit quadratic
charge equation
\[
 \sum_{i\in E}a_i a_i^{\mathsf T}=0,
 \tag{19}
\]
where \(a_i\) is column \(i\) of \(A\). This connects dimension counting
to the multiplicities of the column types. The full-rank condition in
(14) remains necessary.

## 5. Exact crossings and the orientation selector

Each admissible syndrome has \(2^s\) choices of single-mode orientation.
Therefore
\[
 D(A;P,L)=
 \sum_{\substack{\sigma\subseteq[P]\\
                  \operatorname{rank}A_E=k\\
                  A_EA_E^{\mathsf T}=0}}
 2^{|\sigma|}
 \binom{P-|\sigma|-2k}{(L-|\sigma|)/2-k}
 \tag{20}
\]
is the full homogeneous input dimension. Nonintegral or out-of-range
binomial arguments contribute zero.

The mode-complement involution \(J\) exchanges the orientations of every
nonempty syndrome and leaves its even-coordinate basis unchanged.
Retain the orientations whose least singly occupied pair uses its low
mode. Omit the zero-syndrome component. The selected dimension is
\[
 d=\frac{D-D_0}{2},
 \tag{21}
\]
where
\[
 D_0=
 \begin{cases}
 \binom{P-2k}{L/2-k},&L\text{ even and }AA^{\mathsf T}=0,\\
 0,&\text{otherwise}.
 \end{cases}
\]
Here \(D_0\) is the zero-syndrome sector on which \(J\) acts trivially,
not the full fixed subspace of \(J\).

Choose a basis \(I\) of the selected input space and any left inverse
\(\Pi\) supported on the same syndrome/orientation sectors. Such a left
inverse exists by ordinary linear algebra within each sector. Then
\[
 \Pi I=I_d,\qquad I^{\mathsf T}I=0,\qquad\Pi JI=0.
 \tag{22}
\]
The three Pfaffian vertices consequently become
\[
 p^L\longmapsto0,\qquad
 q^L\longmapsto0,\qquad
 r^L\longmapsto E_d.
 \tag{23}
\]
The entire restricted tensor is divisible by \((q+r)^k\).

This is an explicit finite algebraic construction. Bases and left
inverses can be found using the coordinate changes in Theorem 2.
Their coefficients are minors of finite matrices. The dimension and
factor results do not assert a new uniform low-space implementation or
an improved complete rank-sum exponent.

For example, with two forms and binary column types \(1,2,3\), the
zero-Gram condition says that the number of active even columns of
each type is even. The overlapping examples below are exact:

| Column-type multiplicities | \(P=L\) | Full \(D\) | Selected \(d\) | Guaranteed factor |
|---|---:|---:|---:|---|
| \((2,2,2)\) | 6 | 14 | 6 | \((q+r)^2\) |
| \((3,2,2)\) | 7 | 68 | 34 | \((q+r)^2\) |
| \((3,3,2)\) | 8 | 312 | 156 | \((q+r)^2\) |

The dimensions alone do not determine the complete cut ranks.

## 6. Disjoint blocks, with mixed local grades

Partition all \(P\) pairs into \(k\) blocks, and take their indicators
as the rows of \(A\). In a syndrome, the zero-Gram condition says that every block
contains an even number of even positions. On each block, this is
precisely the ordinary radical condition from Section 3.

One may allow every local exterior grade of the required parity while
fixing only the total grade \(L\). The tensor product of the local
radical spaces is contained in the product image. Hence every
coefficient is divisible by \((q+r)^k\).

If \(R_b(z)\) is the ordinary radical dimension polynomial for a block
of \(b\) pairs, the dimension before the orientation selector is
\[
 [z^L]\prod_{j=1}^k R_{b_j}(z).
 \tag{24}
\]
This also follows from (14): summing the products of local binomial
dimensions at fixed total degree is Vandermonde's identity. It does
not require choosing a separate fixed grade on each block.

## 7. Nested global pivots collapse to scalar padding

There is a useful limit to the preceding construction. Start with a
global-pivot embedding
\[
 F(v)=\Omega_E\Lambda_E(z_pv),
\]
where \(p,q\) are two new even pairs and \(v\) uses neither.
Suppose the inner vector is already in a radical:
\[
 \Omega_{\mathrm{inner}}v=\Lambda_{\mathrm{inner}}v=0.
\]
Then
\[
 \Lambda_E(z_pv)=v,\qquad
 F(v)=(z_p+z_q)v.
 \tag{25}
\]
Thus every additional outer layer is a fixed two-pair factor.

For that factor, take input \(z_p+z_q\) and output coordinate \(z_q\).
The exact one-label Pfaffian polynomial is
\[
 q^2+r^2=(q+r)^2.
 \tag{26}
\]
Indeed, the two diagonal input choices contribute \(p^2\) and
\(p^2+q^2+r^2\), while the mixed choices vanish.

Consequently, \(k\) nested layers equal
\[
 (q+r)^{2k-2}
 \times\text{the single-layer tensor on }
 (P-2k+2,\ L-2k+2).
 \tag{27}
\]
The same core labels and projections remain after this reduction.

Multiplication of a vector-valued polynomial by a nonzero scalar
polynomial does not change the span of its coefficient vectors. One
proof applies any linear functional annihilating the new coefficients:
the scalar ring is an integral domain, so it also annihilates the old
polynomial. Hence (27) adds no new coefficient-space rank advantage.
Independent or genuinely overlapping constraints are needed for a
different construction.

## 8. Verification

Run
\[
 \texttt{python tools/check\_noncommutative\_crossing\_radical\_ward.py}
\]
with optional output
\(\texttt{--output evidence/noncommutative-crossing-radical-ward.json}\).

The checker verifies:

- all 512 local generic-polynomial equations in (3), with the state
  creation/contraction conventions checked explicitly;
- polynomial sandwich divisibility for one-form, disjoint-form, and
  overlapping-form examples, including nonzero power-two sandwiches;
- the exact image rank \(2^{t-k}\), the homogeneous intersection
  dimensions, and the zero-Gram/full-rank conditions in every syndrome
  of ten binary weight matrices through \(P=8\);
- direct sparse-vector instances of (25) and the full polynomial
  calculation (26).

These checks audit the ordinary proofs. Complete coefficient-rank
bounds for the overlapping codes remain a separate research task.

