# Complete cubic spaces from varying complementary outputs

## Results and scope

For the seven-label cyclic cubic construction, allowing every rank-one
complementary output map gives an exact **786-dimensional** coefficient
space over every field of characteristic two. Its column rank has a
closed reduction to a matrix of size at most 28, instead of 49.

Allowing complementary maps of rank at most two already gives the entire
**1,225-dimensional major-symmetric tensor space**. Allowing higher ranks
does not enlarge that mixed space.

These are descriptions of complete spaces in which to search for crossing
decompositions. They do not give a strict seed, improve the factoring exponent,
or prove optimality of the existing cost-36 decomposition. They concern
mixtures of coefficients from different output maps, not the smaller
coefficient space of one fixed map.

The prior [154-map construction](https://github.com/1337isnot1337/integer-factoring/blob/7154cec09b982b54a0fcc0bb7a303318fc1f90da/research/cubic-varying-image-affordable-escape.md)
and its [177-dimensional normal form](https://github.com/1337isnot1337/integer-factoring/blob/7154cec09b982b54a0fcc0bb7a303318fc1f90da/research/cubic-varying-image-normal-form-and-peels.md)
use only a specified subset of rank-one maps. The present results allow all
image vectors and all covectors. The physical cubic tensor, input, and cyclic
output selector are unchanged.

## 1. Conventions and the physical identities

Let \(F\) have characteristic two and put \(d=7\). Write tensor coordinates
as \(T(x,y,v,u)\). The two flattenings are

\[
R(T)_{(x,y),(v,u)}=T(x,y,v,u),\qquad
C(T)_{(x,u),(v,y)}=T(x,y,v,u).
\]

The crossing is \(E(x,y,v,u)=[x=v][y=u]\). Let

\[
N(x,y,v,u)=[x=u][y=v],\qquad
D(x,y,v,u)=E(x,y,v,u)[x-y=\pm1\pmod7].
\]

For a complementary output matrix \(C_0\), the output map is
\(\Pi=[I\ C_0]\) on the selected low/high cubic outputs. Use \(C_0\)
for this map to distinguish it from the column flattening \(C(T)\).

The prior [complementary-output identity](https://github.com/1337isnot1337/integer-factoring/blob/7154cec09b982b54a0fcc0bb7a303318fc1f90da/research/cubic-complementary-output-maps.md)
gives

\[
[q^3]G_{C_0}(x,y,v,u)=(C_0)_{xy}(C_0)_{uv}.                 \tag{1}
\]

At \(C_0=0\), the five raw coefficients at
\((r^3,pr^2,p^2r,qr^2,q^3)\) are

\[
(E,\ D+N,\ N+\Delta,\ D,\ 0),\qquad
\Delta(x,y,v,u)=[x=y=v=u].                                \tag{2}
\]

The physical parent has zero mixed low/high output blocks, and its high/high
blocks satisfy simultaneous output/input swap symmetry. Its whole parent has
that same major symmetry. These are identities over \(F_2\) of fixed finite
matrices. The checker regenerates them directly from the local Pfaffian,
before applying an output map; therefore they hold for arbitrary \(F\)-valued
output coefficients. No extrapolation in \(d\) is used here.

## 2. Every rank-one map: the exact 786-dimensional space

Let \(\mathrm{Sym}_d(F)\) mean the vector space of symmetric
\(d\)-by-\(d\) matrices. Define

\[
\mathcal S=\mathrm{span}\{T_{A,B}:A,B\in
\mathrm{Sym}_d(F)\},\qquad
T_{A,B}(x,y,v,u)=A_{xu}B_{yv}.                              \tag{3}
\]

**Theorem 1.** The full mixed coefficient space of all maps
\(C_0=ab^T\) is

\[
\boxed{\mathcal W=\mathcal S+\mathrm{span}\{E,D\}.}  \tag{4}
\]

**Upper inclusion.** The zero mixed parent blocks make the change from the
baseline quadratic in \(C_0\). For a rank-one map it has the form
\(a_xa_uK(y,v)\). Parent swap symmetry makes \(K\) symmetric, so every
change belongs to \(\mathcal S\). Equation (2) puts the baseline in the
right side of (4), since both \(N\) and \(\Delta\) belong to
\(\mathcal S\).

**Lower inclusion.** Equation (1) supplies
\(T_{aa^T,bb^T}\) as an actual raw coefficient. The binary vectors

\[
\mathcal U=\{e_i:0\le i<7\}\cup
\{e_i+e_j:0\le i<j<7\}                                  \tag{5}
\]

have outer squares forming a basis of \(\mathrm{Sym}_7(F)\).
Indeed, the diagonal squares give \(e_ie_i^T\), and

\[
(e_i+e_j)(e_i+e_j)^T+e_ie_i^T+e_je_j^T
=e_ie_j^T+e_je_i^T.
\]

Thus the \(28^2=784\) maps \(ab^T\), with \(a,b\in\mathcal U\),
already supply every direction of \(\mathcal S\). The zero map supplies
\(E,D\) by (2).

The dimension of \(\mathcal S\) is \(28^2=784\). No nonzero
combination of \(E,D\) belongs to it: every column matrix from
\(\mathcal S\) kills \(e_i\otimes e_j+e_j\otimes e_i\), while
\(eE+aD\) multiplies that vector by \(e+a[i-j=\pm1]\).
A cycle edge and a nonedge force \(a=e=0\). Therefore
\(\dim\mathcal W=786\), over every characteristic-two field.

## 3. Exact column ranks from a smaller matrix

Use the symmetric-matrix basis consisting of seven diagonal units followed by
the 21 alternating edge matrices. Let \(V\) have their vectorizations as its
columns. Every \(T\in\mathcal W\) has unique coordinates

\[
C(T)=\mathrm{diag}_{(i,j)}(g_{ij})+VQV^T,\qquad
g_{ii}=e,\quad g_{ij}=e+a[i-j=\pm1]\ (i\ne j),              \tag{6}
\]

where \(Q\) is an arbitrary \(28\)-by-\(28\) matrix. In particular,
\(Q\) is **not** required to be symmetric.

Let \(Z\) be the set of unordered off-diagonal pairs with \(g_{ij}=0\),
and let \(L\) contain the seven diagonal labels and the pairs in \(Z\).
Let \(I_{\mathrm{diag}}\) be diagonal on \(L\), equal to one on its seven
diagonal labels and zero elsewhere.

**Theorem 2.**

\[
\boxed{\mathrm{rank}C(T)
=2(21-|Z|)+\mathrm{rank}(Q_{L,L}+eI_{\mathrm{diag}}).} \tag{7}
\]

**Proof.** For each unordered pair choose the ordered basis

\[
a_{ij}=e_i\otimes e_j+e_j\otimes e_i,\qquad
b_{ij}=e_i\otimes e_j.
\]

Place all \(a_{ij}\) first, then the diagonal vectors, then all \(b_{ij}\).
The symmetric update \(VQV^T\) kills every \(a_{ij}\) and has image in
the span of the \(a_{ij}\) and diagonal vectors. The diagonal graph operator
acts by \(g_{ij}\) on both \(a_{ij},b_{ij}\).

For each nonzero \(g_{ij}\), its \(a_{ij}\) column and its \(b_{ij}\) row
provide two independent pivots. Removing them introduces no Schur-complement
correction: the first column has no other entries, and the last row has no
other entries. The remaining nonzero rectangle has rows indexed by the
zero-weight \(a_{ij}\) and diagonal vectors, and columns indexed by the
diagonal vectors and zero-weight \(b_{ij}\). After permuting its row/column
orders, it is precisely \(Q_{L,L}+eI_{\mathrm{diag}}\). This proves (7).

| Parameters | Forced column rank | Remaining matrix size |
| --- | ---: | ---: |
| \(e=a=0\) | 0 | 28 |
| \(e=0,\ a\ne0\) | 14 | 21 |
| \(e=a\ne0\) | 28 | 14 |
| \(e\ne0,\ e\ne a\) | 42 | 7 |

The last line includes \(a=0\). Formula (7) gives an exact rank, not just a
lower bound or a factor-channel count. Row ranks still require their own
calculation; the theorem does not infer them from column ranks.

## 4. Rank-two maps already give the entire major-symmetric space

Let \(J\) swap the two indices in a \(49\)-coordinate row-cut column.
Define the major-symmetric space by

\[
\mathcal M=\{T:R(T)J=(R(T)J)^T\}.                          \tag{8}
\]

**Theorem 3.** The full mixed coefficient space of arbitrary complementary
maps is \(\mathcal M\), of dimension \(49\cdot50/2=1225\). Maps of
rank at most two suffice to span it.

Parent major symmetry is preserved when the same output map is applied to
both output legs. Thus every raw coefficient lies in \(\mathcal M\).
Conversely, (1) gives

\[
R([q^3]G_{C_0})J=\mathrm{vec}(C_0)
                         \mathrm{vec}(C_0)^T.         \tag{9}
\]

Take \(C_0\) to be one matrix unit or the sum of two distinct matrix units.
Every such map has rank at most two. Their 1,225 outer squares form a basis
of the symmetric \(49\)-by-\(49\) matrices, by the same diagonal/off-diagonal
argument used in (5). This proves both claims over every characteristic-two
field.

This theorem identifies the precise limit of enlarging the complementary
map while keeping the two physical output maps equal. The resulting
optimization is already a general major-symmetric crossing problem. It
does not mean that any one map supplies all 1,225 directions, or that a cheap
decomposition exists in this space.

## 5. Verification and attribution

Run the standalone standard-library checker:

```bash
python tools/check_cubic_full_image_span.py --output /tmp/cubic-full-image-span.json
```

The [recorded result](../evidence/cubic-full-image-span.json) contains the
source hash and scope. The checker regenerates 3,925 raw matrices from the
784 spanning rank-one maps and the zero map, checks 18,847,850 entries in
both flattenings against (4), and proves the 784- and 786-dimensional span
claims by exact binary elimination. It checks the parent identities before
specialization, verifies the 1,225-dimensional outer-square basis, and
compares (7) against full matrix ranks over \(F_2,F_4\), with two pivot orders.

The physical Pfaffian/cubic construction and earlier restricted-map results
are from PR #29 at the revision linked above. The underlying modular crossing
operators come from Anthropic's Appendix B. This extension, proof, and checker
used substantial OpenAI assistant support. No new Lean proof is claimed.
