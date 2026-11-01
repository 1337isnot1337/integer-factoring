# Fixed-image rank-one outputs: a universal all-power barrier

## Theorem and scope

Fix any vector \(u\in K^7\), where \(K\) is any field of characteristic two.
Allow every complementary output map
\[
 C=uv^T,\qquad v\in K^7,
\]
in the cubic construction with the fixed cyclic output selector.
Let \(\mathcal V_u\) be the full combined coefficient space of all these maps.
Then, for every positive integer \(n\), every decomposition of \(E^{\otimes n}\) into tensors from \(\mathcal V_u^{\otimes n}\) has rank cost
\[
 \boxed{\quad S_n\geq252^{n/2}>7^{\sqrt2 n}.\quad}
\]
The vector \(v\), its scalar multiples, and the coefficient combinations may all vary.
The output image vector \(u\) remains fixed within the local coefficient space.
The proof also allows different fixed image vectors at different tensor positions.

This closes the dense rank-one pair in [the mixed F8 search](cubic-mixed-F8-spaces.md), with substantially broader parameter and field scope.
It does not cover arbitrary mixtures of different image vectors within one local space.

## An enlarged tensor class

Let \(\mathcal B\) be the baseline coefficient space for \(C=0\).
Write \(e(B)\) for its crossing coefficient, so \(e(E)=1\).
The exact physical parent kernel has two useful properties:

1. All coefficient blocks with one low output and one high output vanish.
2. The high/high coefficients are invariant under simultaneously exchanging the two physical outputs and the two logical inputs.

Consequently, a fixed-image output perturbation contributes tensors of the form
\[
 U(x,y,v,w)=u_xu_w K(y,v),
\]
where \(K\) is symmetric. The symmetry follows after summing the physical high-output coefficients against the same covector on both output legs.

It is enough to prove the theorem for the larger linear space
\[
 \mathcal W_u
 =\{B+U:B\in\mathcal B,\ U(x,y,v,w)=u_xu_wK(y,v),\ K=K^T\}.
\]
The matrix \(K\) is arbitrary in this enlarged space; no realizability assumption is needed.
The parent identities and the baseline restrictions below are checked entrywise over F2, so they hold over every characteristic-two field.

## Exact baseline restrictions

Let
\[
 \mathcal P=\{(a,b):0\leq a<b\leq6,\ b-a\notin\{1,6\}\},
\]
the fourteen unordered nonedges of the seven-cycle, and let \(\mathcal O\) contain both orientations of these pairs.
For every baseline coefficient tensor \(B\), direct matrix identities give
\[
 R(B)\big|_{\mathcal P\times\mathcal P}=e(B)I_{14},
 \qquad
 C(B)\big|_{\mathcal O\times\mathcal O}=e(B)I_{28}.
\]
These are identities for the actual baseline family. They do not assume that an arbitrary change of logical subspace preserves the cyclic output selector.

## Nine-dimensional row compression

On the row restriction, an update has entries
\[
 R(U)\bigl[(x,y),(v,w)\bigr]=u_xu_wK(y,v).
\]
Group the fourteen row coordinates by their second label \(y\).
There are five groups, with sizes \(1,2,3,4,4\).
Within each group, every column of every update is proportional to the fixed vector whose entries are \(u_x\).
Thus all updates have their row image in a common subspace of dimension at most five.

Choose a rank-nine map \(L_u:K^{14}\to K^9\) annihilating this common subspace, and choose a right inverse \(P_u\).
These maps depend only on \(u\), not on \(K\) or on the baseline coefficient tensor.
Then
\[
 L_u\,R(T)\big|_{\mathcal P\times\mathcal P}\,P_u=e(T)I_9
 \qquad(T\in\mathcal W_u).
\]
This identity also shows that the functional \(e\) is well defined on the enlarged space, even if its baseline-plus-update representation is not unique.

For the recorded dense F8 vector, the checker constructs the eliminators explicitly: in each group, remove its first coordinate and subtract the appropriate multiple of that row from the others.
For image vectors with zero coordinates, the common image dimension can only decrease, so a rank-nine annihilator still exists.

## Twenty-eight-dimensional column triangular form

On the column restriction, every update is an outer product
\[
 C(U)\big|_{\mathcal O\times\mathcal O}=ab^T,
 \qquad
 a_{(x,w)}=u_xu_w,\quad b_{(v,y)}=K(y,v).
\]
The vector \(a\) is common to all updates. Symmetry of \(K\) and the presence of both orientations of each nonedge give
\[
 b^Ta
 =\sum_{\{i,j\}\in\mathcal P}u_iu_j\bigl(K(i,j)+K(j,i)\bigr)=0.
\]
If \(a=0\), all restricted updates vanish.
If \(a\ne0\), choose a basis beginning with \(a\).
Every update has image in the first basis line and annihilates that line. Its matrix therefore has possible entries only in the first row, with zero first column.
All updates are simultaneously strictly upper triangular in this constant basis.

Thus every restricted column matrix in \(\mathcal W_u\) is simultaneously upper triangular with all 28 diagonal entries equal to \(e(T)\).
The dimension remains 28: the common square-zero image does not require discarding a coordinate.

## Tensor-power conclusion

The row identity tensors to a scalar compression of dimension \(9^n\).
The column triangular form tensors to a triangular restriction of dimension \(28^n\).
Both statements survive arbitrary linear recombinations in \(\mathcal W_u^{\otimes n}\).
Therefore
\[
 e^{\otimes n}(T)\ne0
 \quad\Longrightarrow\quad
 \operatorname{rank}R(T)\geq9^n,
 \qquad\operatorname{rank}C(T)\geq28^n.
\]
Since \(e^{\otimes n}(E^{\otimes n})=1\), at least one piece of every crossing decomposition has nonzero value of this functional.
That piece alone has cost at least \(\sqrt{9^n28^n}\).
The exact inequalities
\[
 252^6>7^{17},\qquad17^2>2\cdot12^2
\]
show \(252^{1/2}>7^{17/12}>7^{\sqrt2}\).

## Certificate

The [checker](../tools/check_cubic_fixed_image_rank_one.py) verifies:

- all 4,900 baseline restriction coefficients;
- all 24,010 vanishing low/high parent entries and 12,005 high/high symmetry entries;
- the polynomial cancellation \(b^Ta=0\) for a basis of all 28 symmetric matrices, with formal image-vector variables;
- all 648 coefficients of the explicit nine-dimensional row compression for the recorded dense case;
- all 6,272 coefficients of the corresponding column factorizations, including every common-image annihilation identity.

The [evidence](../evidence/cubic-fixed-image-rank-one.json) records the row groups, explicit eliminators, column image vector, update covectors, and source hashes.

```bash
python tools/check_cubic_fixed_image_rank_one.py --output /tmp/cubic-fixed-image-rank-one.json
```
