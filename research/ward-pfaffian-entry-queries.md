# Polynomial-time physical Ward entry queries

A physical occupation-basis entry of the combined Ward kernel can be
evaluated by one Pfaffian of size at most \(k+1\). This works for
arbitrary overlaps among the \(k\) forms. It uses
\(O(Pk^2+k^3)\) field operations and \(O(k^2)\) field elements of
workspace, including exact handling of zero local entries.

The statement holds over every field of characteristic two. It is
an entry-query result, not a coefficient-matrix rank bound. The
nonlocal stabilized basis is outside its automatic scope.

## 1. The scalar coefficient problem

Let \(F\) be a field of characteristic two, let
\(A=[a_1\ \cdots\ a_P]\in F^{k\times P}\), and let
\(m_i,b_i\in F\). Define
\[
 \mathcal F(A;m,b)
 =[z_1\cdots z_k]\prod_{i=1}^P
       \left(m_i+b_i\sum_{a=1}^k A_{ai}z_a\right).
 \tag{1}
\]
The coefficient can be computed in the ordinary commutative
polynomial ring or in its quotient by \(z_a^2=0\); the answers
coincide. Expanding the coefficient gives
\[
 \mathcal F(A;m,b)
 =\sum_{\substack{S\subseteq[P]\\|S|=k}}
   \det A_S\prod_{i\in S}b_i\prod_{i\notin S}m_i.
 \tag{2}
\]
Indeed, the sum over the assignments of the \(k\) forms to a chosen
set of \(k\) sites is a permanent. In characteristic two it equals
the determinant. This is the characteristic assumption used below;
the same formula does not hold for a general odd-characteristic
commutative coefficient problem.

For the physical Ward kernel, write
\[
 M_P=\bigotimes_{i=1}^P S_i,\qquad
 T_i(\eta)=R_i+\eta O_i,\qquad
 W(\eta)=\prod_{a=1}^k\left(\sum_i A_{ai}T_i(\eta)\right).
\]
The [local operator identities](ward-form-operator-transfer.md)
give \(R_i^2=O_i^2=R_iO_i=O_iR_i=0\). Operators at distinct
sites commute, so the product over forms has the site-minor
expansion in (2). For a fixed tensor-product source and target, take
\[
 m_i=(S_i)_{\mathrm{out}_i,\mathrm{in}_i},\qquad
 b_i=(T_i(\eta)S_i)_{\mathrm{out}_i,\mathrm{in}_i}.
 \tag{3}
\]
Then (1) is precisely that entry of \(W(\eta)M_P\). Local
particle-hole permutations preserve this description. A nonlocal
change of basis need not preserve the scalar site product (1).

## 2. No zero ordinary factors

First assume every \(m_i\ne0\). Put
\[
 \beta_i=b_i/m_i,\qquad
 B=A\operatorname{diag}(\beta_i),\qquad
 u=\sum_i\beta_i a_i,\qquad
 Q=\sum_i\beta_i^2a_i a_i^{\mathsf T}.
 \tag{4}
\]
Let \(J\) be the \(P\times P\) alternating matrix with every
off-diagonal entry equal to one. In characteristic two,
\[
 H=BJB^{\mathsf T}=uu^{\mathsf T}+Q.
 \tag{5}
\]
The diagonal is zero: squaring a sum equals the sum of its squares.

For even \(k\), the Pfaffian minor-summation identity gives
\[
 \operatorname{Pf}(BJB^{\mathsf T})
 =\sum_{|S|=k}\det(B_S)\operatorname{Pf}(J_S)
 =\sum_{|S|=k}\det(B_S).
 \tag{6}
\]
The last equality follows because \(J_S\) has one term for every
perfect matching, and their number \((k-1)!!\) is odd. The empty
Pfaffian is one.

For completeness, the first identity can be seen directly by
expanding the Pfaffian into matchings and then expanding each
matrix entry. Fix the assignment of a column of \(B\) to each row.
If two rows receive the same column, terms cancel by exchanging
their matched partners; the term that pairs those two rows is zero
because \(J\) has zero diagonal. If the assignment is injective,
its contribution is the Pfaffian on its selected column set.
Summing the bijections onto that set supplies the determinant in
characteristic two. This proves (6) without a division by a
factorial.

For odd \(k\), append a new row and one dummy column:
\[
 \widetilde B=
 \begin{pmatrix}B&0\\0&1\end{pmatrix}.
\]
Every nonzero maximal minor must contain the dummy column. Applying
(6) in dimension \(k+1\) gives
\[
 \sum_{|S|=k}\det(B_S)
 =\operatorname{Pf}
    \begin{pmatrix}H&u\\u^{\mathsf T}&0\end{pmatrix}.
 \tag{7}
\]
Thus multiply the Pfaffian in (6) or (7) by \(\prod_i m_i\) to
obtain (1).

## 3. Exact handling of zero ordinary factors

Let
\[
 Z=\{i:m_i=0\},\qquad N=[P]\setminus Z,\qquad r=|Z|.
\]
Every site in \(Z\) must be selected in (2). The answer is zero if
one of its \(b_i\) is zero, if \(r>k\), or if the columns \(A_Z\)
are dependent.

Otherwise choose an invertible row transformation \(R\) such that
\[
 RA_Z=\begin{pmatrix}I_r\\0\end{pmatrix}.
 \tag{8}
\]
Use a fixed order of the forced columns when choosing these pivots.
Let \(L\) be the bottom \(d=k-r\) rows of \(R\). For the sites in
\(N\), compute the quantities in (4), then project them:
\[
 u'=Lu,\qquad Q'=LQL^{\mathsf T},\qquad
 H'=u'u'^{\mathsf T}+Q'.
 \tag{9}
\]
There is no contribution from forced sites in \(u\) or \(Q\).
Define
\[
 \mathcal P_d(u',Q')=
 \begin{cases}
  \operatorname{Pf}(H'),&d\text{ even},\\[2mm]
  \operatorname{Pf}\begin{pmatrix}H'&u'\\u'^{\mathsf T}&0\end{pmatrix},
       &d\text{ odd}.
 \end{cases}
\]
Then the complete formula is
\[
 \boxed{\displaystyle
 \mathcal F(A;m,b)
 = (\det R)^{-1}
    \left(\prod_{i\in Z}b_i\right)
    \left(\prod_{i\in N}m_i\right)
    \mathcal P_d(u',Q').}
 \tag{10}
\]

To prove the determinant factor, (2) implies
\(\mathcal F(RA;m,b)=\det(R)\mathcal F(A;m,b)\).
After (8), the forced sites supply \(z_1\cdots z_r\). The
coefficient still needed uses only \(z_{r+1},\ldots,z_k\), so the
remaining columns are exactly their projections \(La_i\).
Equations (6)--(7) now apply in dimension \(d\). When \(d=0\),
the Pfaffian is one and (10) still holds. Row and column swaps
introduce no signs in characteristic two.

## 4. Streamed work and workspace

A single pass through the input sites retains \(u,Q\), one product
of scalars, and at most \(k\) forced columns. Each nonforced site
updates one vector and one symmetric outer-product sum, using
\(O(k^2)\) field operations. The forced-column matrix, its row
transformation, and the two projected moments use \(O(k^2)\) field
elements. Their row reduction and multiplication cost \(O(k^3)\).

The remaining Pfaffian has dimension at most \(k+1\). Paired
elimination computes it in cubic time. With a nonzero pivot
\(h_{01}\), the smaller alternating matrix has entries
\[
 h'_{ij}=h_{ij}+
    \frac{h_{0i}h_{1j}+h_{1i}h_{0j}}{h_{01}},
 \qquad i,j\ge2,
\]
and its Pfaffian is multiplied by \(h_{01}\). If no partner of the
current first row is nonzero, the Pfaffian is zero. A row swap
together with its column swap puts a nonzero partner in the pivot
position when one exists.

For \(k\ge1\), the resulting bounds are
\[
 \boxed{\quad O(Pk^2+k^3)\text{ field operations},\qquad
                O(k^2)\text{ field elements of workspace}.\quad}
 \tag{11}
\]
The case \(k=0\) is the direct product \(\prod_i m_i\). The input
matrix, source and target indices, and local-factor descriptions
are read-only inputs, not included in the workspace. The ordered
stream must make each site column and its local entry pair
available. Index and loop counters use additional logarithmic
bit storage. Field-element workspace is distinct from bit
complexity; the latter also includes the representation and cost
of each field operation.

This removes the exponential form-subset memory for physical scalar
queries even when the support intervals have large overlap. The
[interval-width theorem](ward-interval-operator-transfer.md)
retains its separate role for the nonlocal stabilized operator.

## 5. Exact formal coefficient queries over finite fields

Suppose the local entries in (3) are polynomials in a fixed number
\(v\) of parameters, and a total degree bound \(D\) for the requested
entry is supplied. Let the constant form weights lie in an
explicit finite field \(F_q\) of characteristic two. Work in a
finite extension \(K/F_q\) with more than \(D\) elements and choose
\(D+1\) distinct interpolation nodes.

Evaluate (10) on the \((D+1)^v\) product grid. Tensor-product
Lagrange interpolation gives every formal coefficient exactly.
For one requested coefficient, multiply each scalar query by the
product of its one-variable Lagrange coefficient weights and sum.
This avoids storing the entire coefficient array. Grid points
where some \(m_i=0\) require no exclusion: Section 3 evaluates
them exactly. In particular, the method is not an interpolation
of an unevaluated rational expression at its poles.

For fixed \(v\), this takes polynomially many field operations in
\(P,k,D\). In an explicit binary-field representation, the bit
work is polynomial in \(P,k,D,\log q\). An extension can be made
effective without an unbounded search: if \(q>D\), none is needed.
Otherwise take the least positive \(e\) with \(q^e>D\). Then
\(q^e\le qD\le D^2\). Enumerate monic degree-\(e\) polynomials over
\(F_q\), testing irreducibility by polynomial gcds with
\(X^{q^j}-X\) for \(1\le j\le\lfloor e/2\rfloor\).
A reducible polynomial has a factor of one of those degrees, and
an irreducible degree-\(e\) polynomial exists. Modular repeated
squaring computes the needed powers. This elementary construction
already has polynomial bit cost for the stated regime.

For the local physical Pfaffian, there are at most eight occupied
vertices and hence at most four edge factors. If the edge weights
are affine in \(p,q,r\), every \(m_i\) has total degree at most
four. Since \(R_i,O_i\) are parameter-independent and \(T_i\) is
linear in \(\eta\), an inserted local factor has total degree at
most five in \((p,q,r,\eta)\). Each contributing term chooses
exactly \(k\) inserted sites. A safe degree bound is therefore
\[
 D=4P+k.\tag{12}
\]
An additional external power of \(\eta\) must be added to this
bound, or handled as a known coefficient-index shift. More
generally, local ordinary degree at most \(d_0\) and inserted
degree at most \(d_0+1\) give \(D=d_0P+k\).

These are exact formal **entry** queries. They do not calculate the
rank of a coefficient matrix containing exponentially many such
entries, and they do not evaluate an arbitrary nonlocal stabilized
entry without further work.

## 6. Implementation and exact controls

The finite-field implementation is
[ward_pfaffian_queries.py](../tools/ward_pfaffian_queries.py).
It stores the nonforced moments before eliminating the forced
columns, so no earlier input columns need to be revisited. Its
paired elimination accepts singular alternating matrices and the
empty matrix.

Run

```sh
python3 tools/check_ward_pfaffian_queries.py --output /tmp/ward-pfaffian-queries.json
```

The [checker](../tools/check_ward_pfaffian_queries.py) compares the
streamed result with direct squarefree-polynomial multiplication:
5,520 exhaustive binary inputs and 2,160 seeded cases over
\(\mathbf F_2,\mathbf F_4,\mathbf F_8,\mathbf F_{16}\). It tests
dependent forced columns, too many forced sites, nontrivial row
scalings, and zero local factors. Another 960 controls compare
paired Pfaffian elimination with the explicit matching expansion.

The checker also reconstructs 27 two-variable polynomial queries
on 567 interpolation points, including 357 points with a zero
ordinary factor, and compares 108 literal one-, two-, and
three-site physical operator entries with direct global operator
action. The [evidence](../evidence/ward-pfaffian-queries.json)
records these counts and hashes every source used.
