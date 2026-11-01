# An overlapping Ward family with interval width two

This construction keeps the number of active forms bounded while the
number of forms grows. Its label dimensions have an exact two-state
recurrence. The central-grade dimension has the same exponential rate
as the total dimension. These are construction and entry-evaluation
results; the first complete rectangular rank bound below does not give
a factoring-score improvement.

## 1. The form matrix

For \(k\ge1\), take the binary column types
\[
 e_1,\ e_1+e_2,\ e_2,\ e_2+e_3,\ldots,
 e_{k-1}+e_k,\ e_k.
 \tag{1}
\]
Repeat each type twice, keeping the displayed order. The resulting
matrix \(A_k\) has \(k\) rows and
\[
 P=2(2k-1)=4k-2
\]
physical pair sites. It has full row rank. Its total Gram matrix is
zero, since each column appears twice.

The \(2k-1\) distinct quadratic charges \(aa^{\mathsf T}\) are
linearly independent. The charges of the \(e_i\) form the diagonal
basis. The charge of \(e_i+e_{i+1}\) has an off-diagonal entry at
\((i,i+1)\) that appears in no other type. Consequently the quadratic
charge rank is exactly \(2k-1\).

For \(k\ge2\), at most two form-support intervals are active at any
site in this order. At a site of type \(e_i+e_{i+1}\), the inclusive
prefix has rank \(i+1\) and the inclusive suffix has rank \(k-i+1\).
Thus their rank sum minus \(k\) equals two. The
[fixed-order basis theorem](ward-interval-operator-transfer.md)
therefore shows that the optimum interval width is exactly two.
For \(k=1\) it is one.

The stabilized \(k\)-form Ward operator consequently has width at
most \(8\cdot2^2=32\), independent of \(k\). A fixed stabilized
operator entry takes \(O(P)\) field operations and a constant number
of field elements of workspace, under the ordered input-access
conditions of the interval theorem. Position and grade counters,
input storage, and the bit representation of field elements are
separate. This does not bound a coefficient-matrix rank.

## 2. Admissible subsets are connected spanning subgraphs

The product-image construction requires an even-position set \(E\)
with
\[
 \operatorname{rank}A_E=k,
 \qquad A_EA_E^{\mathsf T}=0.
 \tag{2}
\]
Independence of the distinct quadratic charges implies that a
zero-Gram set contains either both copies or neither copy of every
type. Identify \(e_i\) with an edge from a distinguished vertex \(0\)
to vertex \(i\), and \(e_i+e_{i+1}\) with the edge from \(i\) to
\(i+1\). The \(2k-1\) types are the edges of this fan graph: one path
on vertices \(1,\ldots,k\), with all its vertices adjacent to \(0\).

For any selected type set, the form columns are its reduced incidence
matrix. A vector in its left kernel is a vertex potential that is
constant on every selected-edge component, with potential zero at
vertex \(0\). Its column rank is therefore
\(k+1-\#\{\text{components}\}\). Condition (2) holds exactly when the
selected type edges form a connected spanning subgraph.

Let \(N_{k,t}\) count these subgraphs with \(t\) selected type edges.
These are ordinary integer counts. Their edge-counting generating
polynomial is not taken modulo two.

## 3. Exact two-state counting recurrence

Build the graph in the order \(1,\ldots,k\). Let \(C_j(z)\) count
partial edge sets in which every processed vertex is connected to
the distinguished vertex \(0\). Let \(U_j(z)\) count partial sets
whose only component not yet connected to \(0\) contains the most
recent vertex \(j\). An older component that no longer contains
the current vertex can never acquire another path edge and is
discarded.

When vertex \(j+1\) is introduced, its new spoke and its new path
edge give
\[
 \begin{pmatrix}C_{j+1}(z)\\U_{j+1}(z)\end{pmatrix}
 =
 \begin{pmatrix}2z+z^2&z^2\\1&z\end{pmatrix}
 \begin{pmatrix}C_j(z)\\U_j(z)\end{pmatrix},
 \qquad
 \begin{pmatrix}C_1(z)\\U_1(z)\end{pmatrix}
 =\begin{pmatrix}z\\1\end{pmatrix}.
 \tag{3}
\]
From a connected state, the new vertex joins the connected component
through its spoke, its path edge, or both, giving \(2z+z^2\). If
neither edge is selected, the new vertex starts the pending
component. From a pending state, the path edge is required; selecting
the spoke then closes that component. At the end, only the connected
state is allowed. Hence
\[
 C_k(z)=\sum_tN_{k,t}z^t.
\]

Write \(c_k=C_k(1)\). The numerical transition matrix has trace four
and determinant two, so
\[
 c_0=0,\qquad c_1=1,\qquad
 c_k=4c_{k-1}-2c_{k-2}.
 \tag{4}
\]
Equivalently,
\[
 c_k=\frac{(2+\sqrt2)^k-(2-\sqrt2)^k}{2\sqrt2}.
 \tag{5}
\]

## 4. Total and central-grade dimensions

The [homogeneous dimension formula](ward-code-charge-entropy.md)
assigns a type set with \(t\) edges the contribution
\[
 2^{P-2t}x^{P-2t+2k}(1+x^2)^{2t-2k}.
\]
Thus the complete dimension polynomial is
\[
 D_k(x)=\sum_{t=k}^{2k-1}
 N_{k,t}\,2^{P-2t}x^{P-2t+2k}(1+x^2)^{2t-2k}.
 \tag{6}
\]
Its total dimension and its central-grade dimension are different:
\[
 \begin{aligned}
 D_k^{\rm total}&=D_k(1)=2^{P-2k}c_k=2^{2k-2}c_k,\\
 D_k^{\rm central}&=[x^P]D_k(x)
   =\sum_{t=k}^{2k-1}N_{k,t}\,2^{P-2t}
     \binom{2t-2k}{t-k}.
 \end{aligned}
 \tag{7}
\]
The zero-syndrome central-grade summand is the set with every site
even, so its dimension is
\[
 D_k^0=\binom{P-2k}{(P-2k)/2}
      =\binom{2k-2}{k-1}.
 \tag{8}
\]
An orientation selector chooses one of each complementary pair of
nonempty single-site syndromes. Complementation has no fixed
orientation on a nonempty syndrome. The resulting number of central
labels is exactly
\[
 d_k=\frac{D_k^{\rm central}-D_k^0}{2}.
 \tag{9}
\]

| Forms \(k\) | Pair sites \(P\) | Connected type sets \(c_k\) | Total dimension | Central-grade dimension | Selected central labels \(d_k\) |
|---:|---:|---:|---:|---:|---:|
| 2 | 6 | 4 | 16 | 14 | 6 |
| 3 | 10 | 14 | 224 | 174 | 84 |
| 4 | 14 | 48 | 3,072 | 2,140 | 1,060 |
| 5 | 18 | 164 | 41,984 | 26,454 | 13,192 |
| 6 | 22 | 560 | 573,440 | 329,412 | 164,580 |
| 7 | 26 | 1,912 | 7,831,552 | 4,131,916 | 2,065,496 |

For \(0\le u\le k-1\),
\[
 \frac{4^u}{2u+1}\le\binom{2u}{u}\le4^u.
\]
Apply this to each summand of \(7\). It gives
\[
 \frac{D_k^{\rm total}}{2k-1}
 \le D_k^{\rm central}\le D_k^{\rm total}.
 \tag{10}
\]
Moreover \(D_k^0\le4^{k-1}\), which is exponentially smaller than
\(D_k^{\rm total}\). Therefore total dimension, central dimension,
and the positive sequence \(d_k\) all have the same asymptotic
dimension rate per physical pair:
\[
 \lim_{k\to\infty}\frac{\log_2 d_k}{4k-2}
 =\frac{2+\log_2(2+\sqrt2)}4
 =0.942888325790903\ldots.
 \tag{11}
\]
This is a label-dimension rate, not a factoring exponent or a rank
bound.

## 5. First complete rectangular coefficient profile

For \(k=3,P=L=10\), the total Gram rank is zero and the binary joint
vanishing order used by the current engine is six. The
[channel-free rectangular bound](noncommutative-crossing-ward-charge-filtration.md)
gives the following complete coefficient upper profiles:

| Central labels | Integer sum of geometric-mean rank bounds | Earlier one-sided bound | Computed upper-bound exponent |
|---:|---:|---:|---:|
| 174 | 801,804 | 983,460 | 1.3175493069460942 |
| 84 | 192,340 | 234,676 | 1.3729996601803454 |

The stored exponent is \(\log_2 S/(2\log_2 d)\). The largest
explicit stable block has dimension 21. These are upper bounds
obtained after enlarging the actual code to complete rectangles;
they are not actual restricted coefficient ranks. Their failure to
give a seed does not exclude a sharper computation that retains the
fan's form and code correlations.

## 6. Reproduction

Run

```sh
python3 tools/check_ward_fan_family.py --output /tmp/ward-fan-family.json
python3 tools/screen_ward_fan_bounds.py --forms 3 --output /tmp/ward-fan-first-profile.json
```

The [dimension checker](../tools/check_ward_fan_family.py) exhausts
all type subsets for \(k\le7\), checks connectivity against exact
binary rank, and verifies the recurrence, Gram rank, optimal interval
width, and all homogeneous dimensions through \(k=12\). Its
[evidence](../evidence/ward-fan-family.json) includes source hashes.
The [profile tool](../tools/screen_ward_fan_bounds.py) records both
central label counts and the complete rectangular rank profile in
[ward-fan-first-profile.json](../evidence/ward-fan-first-profile.json).
