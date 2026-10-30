# Exact singleton profiles for sparse Boolean holonomy

## Result and scope

Let \(S=\{0,1,2,4,\ldots,2^{b-1}\}\), \(T=2^b-1\), and work
in characteristic two. The matrix

\[
H_S(a,c)=z^{T-(a\mathbin{|}c)}
\]

has an explicit Smith profile at zero. List
\(2^{b-1},2^{b-2},\ldots,1,0\) in descending order and pair
consecutive entries. A pair \((a,c)\) contributes two orders
\(T-a-c\). If a zero remains unpaired, it contributes the order \(T\).

For example, at \(b=8\), the paired sums are \(192,48,12,3\),
with an unpaired zero. The orders are

\[
63,63,207,207,243,243,252,252,255.
\]

The proof below applies to arbitrary distinct positive integer weights,
not only powers of two. It also supplies an exact Gaussian encoding of
the larger Hamming-ball matrices. That encoding does not by itself prove
a Smith formula for balls of radius two or three. The singleton profile
is a control for the sparse-mode construction, not an improved factoring
exponent.

## 1. Weighted singleton theorem

Choose distinct positive integers \(w_1>\cdots>w_b>0\), and include
\(w_0=0\). Choose \(T\) at least the sum of the two largest weights;
when \(b=1\), it suffices that \(T\ge w_1\). Define a symmetric matrix
by

\[
\begin{aligned}
H_{0,0}&=z^T,\\
H_{i,i}=H_{0,i}=H_{i,0}&=z^{T-w_i}\quad(i>0),\\
H_{i,j}&=z^{T-w_i-w_j}\quad(i,j>0,\ i\ne j).
\end{aligned}                                                    \tag{1}
\]

Pair the decreasing list \((w_1,\ldots,w_b,0)\). Its zero Smith
orders are two copies of \(T-a-c\) for each pair \((a,c)\), plus
\(T\) for a final unpaired zero.

Write \(D=\operatorname{diag}(z^{-w_i})\), with the zero-weight
coordinate included. Then

\[
H=z^T D\left(J+\operatorname{diag}(d_i)\right)D,
\qquad d_0=0,\quad d_i=1+z^{w_i}\ (i>0),                   \tag{2}
\]

where \(J\) is the all-one matrix. Laurent factors in (2) are only
notation; every entry of (1) is polynomial under the assumption on \(T\).

## 2. Regular elimination of the two largest weights

The induction maintains the same matrix shape with \(J\) replaced by
\(uJ\), where \(u\in F[[z]]^\times\), and

\[
\operatorname{ord}(u+1)>a                                  \tag{3}
\]

for the largest remaining positive weight \(a\). Initially \(u=1\),
so (3) holds with infinite order.

Let \(a>c>0\) be the two largest remaining weights and put
\(\theta=T-a-c\). Their principal block, after division by
\(z^\theta\), has diagonal entries of orders \(c\) and \(a\),
and its two off-diagonal entries are the unit \(u\). Its reduction
at zero is therefore

\[
\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\]

The block has two Smith orders equal to \(\theta\). Every entry
coupling it to a remaining weight \(w<c\) has order at least
\(\theta\): the two orders are \(T-a-w\) and \(T-c-w\).
Consequently its block Schur elimination uses invertible matrices over
\(F[[z]]\), so it preserves Smith orders.

In the unweighted middle matrix of (2), the Schur complement again has
the form \(u'J+\operatorname{diag}(d_i)\). Direct inversion of the
two-by-two block gives the exact recurrence

\[
\boxed{\displaystyle
u'=\frac{u\,d_a d_c}{d_a d_c+u(d_a+d_c)},
\qquad d_a=1+z^a,\quad d_c=1+z^c.}                          \tag{4}
\]

The denominator is a unit: its constant term is one. Moreover

\[
u'+1=
\frac{(u+1)d_a d_c+u(z^a+z^c)}
     {d_a d_c+u(z^a+z^c)}.
\]

Condition (3) and \(a>c\) show that this numerator has order exactly
\(c\). Thus

\[
\operatorname{ord}(u'+1)=c.
\]

Every remaining positive weight is smaller than \(c\), so the
induction hypothesis is restored. The next paired sum is smaller than
\(a+c\), and hence the next Smith orders are larger than \(\theta\).

## 3. The final zero coordinate

If only the zero coordinate remains, its entry is \(z^T u\), with
order \(T\).

If a positive weight \(a\) remains together with zero, divide their
block by \(z^{T-a}\). In the order \((a,0)\), it is

\[
\begin{pmatrix}
z^{-a}(u+1+z^a)&u\\
u&z^a u
\end{pmatrix}.
\]

Condition (3) makes this a regular matrix whose reduction is
\(\left(\begin{smallmatrix}1&1\\1&0\end{smallmatrix}\right)\).
It is invertible over the power-series ring. Both final orders are
therefore \(T-a\), as the pairing rule asserts.

As a consistency check, the sum of all orders is

\[
(b+1)T-2\sum_{i=1}^b w_i.                                  \tag{5}
\]

Indeed, the middle determinant in (2) is
\(\prod_{i=1}^b(1+z^{w_i})\), a unit. This gives (5) independently.

## 4. An exact Gaussian encoding for every Hamming ball

Let now
\(S_k=\{a\in\{0,\ldots,2^b-1\}:\operatorname{popcount}(a)\le k\}\).
Identify a mask with its subset of the \(b\) weighted bit positions.
Append one additional mode of weight zero. To each subset \(A\) of
size at most \(k\), append this mode exactly when \(|A|\) has parity
different from \(k\). The resulting set \(\widetilde A\) has parity
\(k\), size at most \(k\), and the same total weight. This is a
bijection onto the subsets of \(b+1\) modes of that parity and size.

Let \(w_i\) denote the mode weights, including the added zero. Define
alternating matrices and a cross block by

\[
A_{ij}=z^{-w_i-w_j}\ (i\ne j),\quad A_{ii}=0,
\qquad C=A+\operatorname{diag}(z^{-w_i}).
\]

Then

\[
\boxed{\displaystyle
H_{S_k}(A,B)=z^T\operatorname{Pf}
 \begin{pmatrix}A&C\\C^T&A\end{pmatrix}
 [\widetilde A_+,\widetilde B_-].}                          \tag{6}
\]

The bracket means the principal submatrix on the indicated plus and
minus mode sets, not a matrix product. The total number of selected
modes is even, so the Pfaffian is defined. The empty Pfaffian is one.

Here is a direct characteristic-two proof. Before attaching the factors
\(z^{-w_i}\) to the selected modes, take the alternating matrix whose
off-diagonal entries are all one. Add \(1+z^{w_i}\) on the edge
joining the plus and minus copies of mode \(i\). The Pfaffian of every
even all-one alternating matrix is one in characteristic two, because
its number of perfect matchings is odd. Expanding in the added edges
therefore gives

\[
\prod_{i\in\widetilde A\cap\widetilde B}
\bigl(1+(1+z^{w_i})\bigr)
=z^{\operatorname{wt}(\widetilde A\cap\widetilde B)}.
\]

Multiplying by the factors for the selected modes and by \(z^T\)
produces \(z^{T-\operatorname{wt}(A\cup B)}\), which is the
required entry. The added mode has weight zero and contributes no change.

Equation (6) identifies the sparse matrix with a Gaussian kernel restricted
to a degree filtration. It suggests exterior-algebra methods for higher
radius. Removing an internal Gaussian block can require Laurent, rather
than regular, basis changes; such a step would not automatically preserve
Smith orders. No higher-radius profile formula is asserted here.

## 5. Exact controls

Run:

```sh
python3 tools/check_sparse_holonomy_singletons.py \
  --output evidence/sparse-holonomy-singletons.json
```

The checker compares the pairing rule with the exact adaptive Smith
engine for \(b=1,\ldots,8\), verifies selected independent block-Toeplitz
ranks, and checks nonbinary choices of distinct integer weights. It
checks the unit recurrence (4) over \(\mathbb F_2(z)\) and the Gaussian
identity by a separate recursive Pfaffian computation on small Hamming
balls. All larger-size conclusions above have direct proofs.

Files: [checker](../tools/check_sparse_holonomy_singletons.py),
[tests](../tests/test_sparse_holonomy_singletons.py), and
[evidence](../evidence/sparse-holonomy-singletons.json).
