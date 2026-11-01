# Label entropy for binary Ward codes

Adding annihilating forms gives a deeper Ward factorization, but it
also reduces the number of logical labels. This note makes that
tradeoff explicit for the homogeneous product-image construction.
It applies to binary forms. It does not bound arbitrary crossing
constructions or arbitrary factoring algorithms.

Let \(A\in\mathbf F_2^{k\times P}\) have rank \(k\), with columns
\(a_1,\ldots,a_P\). Let \(D_L\) be the dimension of the full
homogeneous grade-\(L\) product-image code from
[the overlapping-form theorem](noncommutative-crossing-radical-ward.md).
Define the quadratic-charge rank
\[
 q=\dim_{\mathbf F_2}\operatorname{span}
       \{a_i a_i^{\mathsf T}:1\le i\le P\}.
\]
An orientation selector or a smaller crossing code has at most
\(D_L\) labels, so all bounds below apply to it as well.

## 1. Exact total dimension

Call a subset \(E\subseteq[P]\) admissible when
\[
 \operatorname{rank}A_E=k,
 \qquad A_EA_E^{\mathsf T}=0.
\]
Write \(N_{\rm adm}\) for the number of such subsets. For a fixed
admissible \(E\), put \(t=|E|\) and \(s=P-t\). The dimension formula
gives its contribution to the generating polynomial as
\[
 2^s z^{s+2k}(1+z^2)^{t-2k}.
 \tag{1}
\]
Here \(t\ge2k\), since a self-orthogonal subspace of dimension
\(k\) in the nonsingular \(t\)-dimensional bilinear space has
dimension at most \(t/2\).

At \(z=1\), every admissible subset contributes exactly
\(2^{P-2k}\). Consequently
\[
 \boxed{\qquad
 \sum_L D_L=2^{P-2k}N_{\rm adm}
       \le2^{\,2P-2k-q}\le2^{\,2P-3k}.
 \qquad}
 \tag{2}
\]
The first inequality uses the linear map
\(\mathbf F_2^P\to\operatorname{Sym}_k(\mathbf F_2)\) sending
the indicator of \(E\) to \(\sum_{i\in E}a_i a_i^{\mathsf T}\).
Its kernel has exactly \(2^{P-q}\) elements. Some kernel elements
may fail the full-rank condition, so they need not be admissible.

For the second inequality, the diagonal projection sends
\(a_i a_i^{\mathsf T}\) to \(a_i\). Since the columns span
\(\mathbf F_2^k\), this proves \(q\ge k\).

## 2. Minimal charge rank means disjoint form supports

**Proposition.** The equality \(q=k\) holds exactly when an
invertible row change of \(A\) makes every nonzero column a
coordinate vector.

Choose \(k\) independent columns and apply an invertible row change
that sends them to the coordinate vectors. Their outer products
span the entire diagonal-matrix space. If \(q=k\), every other
column outer product is diagonal. Over \(\mathbf F_2\), a vector
with a diagonal outer product has at most one nonzero coordinate.
The converse is immediate.

Thus the forms have disjoint supports precisely in the minimal-charge
case. Genuine overlap requires \(q>k\). Each additional independent
quadratic charge costs a factor two in the upper bound (2).
This does not compare every fixed-grade profile, and it does not
exclude an overlapping code whose actual coefficient ranks are
smaller. It identifies the label cost that such a construction must
offset.

There is an exact refinement for disjoint supports. Let their sizes
be \(b_1,\ldots,b_k\), and allow additional zero columns. An
admissible even set must contain a positive even number of sites
from each support. Hence
\[
 \sum_LD_L
 =2^{2P-3k}\prod_{i=1}^k(1-2^{1-b_i}).
 \tag{3}
\]
If a support has size one, the code is zero. Among disjoint-support
matrices with fixed \(P,k\), the total dimension is maximized by
using no zero columns and taking support sizes as equal as possible.
Indeed, \(\log(1-2^{1-b})\) is increasing and concave for \(b>1\).

For equal supports of fixed size \(b\), the achievable total-dimension
rate is therefore
\[
 \frac1P\log_2\sum_LD_L
 =2-\frac3b+\frac1b\log_2(1-2^{1-b}).
 \tag{4}
\]
At least one homogeneous grade retains this exponential rate, since
there are at most \(2P+1\) grades.

## 3. A bound at each grade

For every real \(z>0\), summing (1) gives
\[
 D(z):=\sum_LD_Lz^L
 =z^{2k}(1+z^2)^{-2k}
   \sum_{E\ {m admissible}}
       (2z)^{P-|E|}(1+z^2)^{|E|}.
 \tag{5}
\]
Drop the full-rank condition and sum over the kernel of the charge
map. Choose \(q\) pivot coordinates of that linear system. Each
choice of its \(P-q\) free coordinates fixes the pivots. Since
\(1+z^2\ge2z\), the pivot contribution is at most
\((1+z^2)^q\), while the free-coordinate sum is
\((1+2z+z^2)^{P-q}\). Therefore
\[
 \boxed{\quad
 D_L\le\inf_{z>0}
 z^{2k-L}(1+z^2)^{q-2k}(1+z)^{2(P-q)}.
 \quad}
 \tag{6}
\]
This is a positive-parameter bound; no coefficientwise inequality
for the rational expression is asserted. For fixed \(z\), its
right side decreases as \(q\) increases. Replacing \(q\) by \(k\)
gives a universal upper bound for this class of binary Ward codes.

When \(k/P\to\kappa\), \(q/P\to\theta\), and \(L/P\to\lambda\),
(6) yields the entropy constraint
\[
 \limsup\frac{\log_2D_L}{P}
 \le\inf_{z>0}\left[
 (2\kappa-\lambda)\log_2z
 +(\theta-2\kappa)\log_2(1+z^2)
 +2(1-\theta)\log_2(1+z)
 \right].
 \tag{7}
\]
The rates satisfy \(\theta\ge\kappa\). At the central grade,
using \(z=1\) already gives the bound \(2-2\kappa-\theta\),
and hence \(2-3\kappa\) without information about the charge rank.

## 4. Consequence for further rank work

The new
[dual-incidence descent](noncommutative-crossing-ward-charge-filtration.md)
depends on \(P,L,k\) through complete rectangular Pfaffian ranks.
The binary joint-order refinement also uses the total Gram rank.
Taking many overlapping forms does not by itself avoid the label
loss in (2) or (6). A successful proportional-\(k\) construction
must balance that loss against the reduction of both cut ranks.

The remaining opportunity is concrete: preserve the actual code
and its form-column correlations inside the rank computation, or
find a fixed-grade profile that improves the comparison after the
charge penalty. Full rectangular bounds alone discard those
correlations. Failure of one such bound is not a lower bound on the
actual restricted tensor.

## 5. Exact verification

The checker
[check_ward_code_charge_entropy.py](../tools/check_ward_code_charge_entropy.py)
enumerates complete small ranges of full-rank column multisets,
checks the homogeneous dimension sums, the charge-kernel counts,
and the minimal-charge classification. It also verifies (5)--(6)
at seven exact rational positive parameters per case.
The output is
[ward-code-charge-entropy.json](../evidence/ward-code-charge-entropy.json).

