# Minimal representations make boundary module certificates complete

## Result and scope

The finite-module hierarchy can detect every strict asymptotic improvement
in the exact grouped-channel word bound of a fixed rational factorization,
including improvement caused by fixed scalar boundaries. First reduce the
auxiliary representation to its reachable and observable part over
\(F(\tau)\). Then apply the free-boundary module construction.

The reduction is necessary. Without it, a coordinate that is invisible to
the terminal vector can prevent every finite precision from detecting a
true exponential gain. A two-dimensional exact example below exhibits
this failure and its one-dimensional repair.

This extends the free-boundary completeness theorem in
[the cancellation-module note](cancellation-module-smith-automata.md).
It concerns one fixed local factorization, one fixed rational evaluation
point, and the sum of bounds for individual channel words. It does not
claim completeness for actual coefficient ranks, which can have additional
cancellation or linear dependence between distinct channel words.
No improved factoring exponent follows from the controls here.

## 1. Keep physical choices grouped by channel

Let \(F\) be a finite field and \(K=F(\tau)\). Fix a finite channel
alphabet \(I\), finite nonempty physical choice sets \(A_i\), matrices
\(U_i(a)\in K^{D\times D}\), and boundary vectors
\(f^T\in K^{1\times D}\), \(e\in K^{D\times1}\).
Let \(s>1\) be rational and let \(h_i\) be integer channel exponents.

For a channel word \(w=i_1\cdots i_t\), write

\[
U_w(\mathbf a)=U_{i_1}(a_1)\cdots U_{i_t}(a_t),
\qquad h(w)=\sum_{j=1}^t h_{i_j}.
\]

The scalar and free-boundary orders are

\[
\begin{aligned}
v_b(w)&=\min_{\mathbf a}\operatorname{ord}_\tau
                   \bigl(f^TU_w(\mathbf a)e\bigr),\\
v_*(w)&=\min_{\mathbf a,r,c}\operatorname{ord}_\tau
                   \bigl(U_w(\mathbf a)_{r,c}\bigr).
\end{aligned}                                                     \tag{1}
\]

Use \(\operatorname{ord}(0)=+\infty\) and \(s^{-\infty}=0\).
Define

\[
Q_b(t)=\sum_{w\in I^t}s^{h(w)-v_b(w)},\qquad
Q_*(t)=\sum_{w\in I^t}s^{h(w)-v_*(w)}.                         \tag{2}
\]

Each channel word is counted once. The physical choices enter through
the minimum in (1), not through an additional sum in (2). This distinction
is essential: two identical physical choices do not double a channel's
weight.

## 2. Exact reachable and observable reduction

Let \(\mathcal S\) be the finite alphabet of all individual matrices
\(U_i(a)\). First restrict the row representation to

\[
\mathcal R=\operatorname{span}_K\{f^TU_{\sigma_1}\cdots
U_{\sigma_m}:m\ge0,\ \sigma_j\in\mathcal S\}.
\]

This space is stable under right multiplication by every matrix in
\(\mathcal S\). It is computed by repeated linear-span closure, in
dimension at most \(D\). Within this restricted representation, take
the quotient by

\[
\mathcal N=\{v\in\mathcal R:
vU_{\sigma_1}\cdots U_{\sigma_m}e=0
\text{ for every }m\ge0\text{ and every choice of }\sigma_j\}.
\]

The space \(\mathcal N\) is also stable. Equivalently, form the span
of the observable columns \(U_{\sigma_1}\cdots U_{\sigma_m}e\) and
use the induced row quotient. Ordinary exact Gaussian elimination over
\(K\) computes both reductions.

Every scalar output \(f^TU_w(\mathbf a)e\) is preserved exactly.
The channel labels and their physical choice sets are preserved as well.
After reduction to dimension \(d>0\), the prefix rows span \(K^{1\times d}\)
and the suffix columns span \(K^{d\times1}\). We call this representation
minimal. If the dimension is zero, every scalar word is identically zero.

Both full spans have bases formed by literal words of length at most
\(d-1\). Indeed, start with the nonzero boundary vector and close its
span under the matrix alphabet. Unless the span is already full, its
dimension must increase at the next step; otherwise it is stable forever.

The exact checker verifies the two intertwining identities underlying
this construction. If \(R\) contains the reachable row basis, and \(C\)
contains the observable column basis in reachable coordinates, then

\[
R U_\sigma=\widehat U_\sigma R,
\qquad
\widehat U_\sigma C=C\widetilde U_\sigma.                     \tag{3}
\]

Together with the corresponding boundary identities, (3) proves equality
of all scalar words, not only the finitely many words used as test cases.

## 3. A finite window compares scalar and free rates

Assume now that the representation is minimal. Choose literal prefix
words \(p_1,\ldots,p_d\) with fixed physical choices whose rows
\(r_j=f^TU_{p_j}\) form a basis. Choose literal suffix words
\(q_1,\ldots,q_d\) with fixed physical choices whose columns
\(c_k=U_{q_k}e\) form a basis. Put these rows into an invertible matrix
\(R\), and these columns into an invertible matrix \(C\).

For every middle word and every one of its physical choices,

\[
U_w(\mathbf a)=R^{-1}\bigl(RU_w(\mathbf a)C\bigr)C^{-1}.     \tag{4}
\]

The \((j,k)\) entry inside parentheses is the scalar output for
\(p_jwq_k\), with its prefix and suffix physical choices fixed.
Let

\[
\eta=\min_{r,c}\operatorname{ord}(R^{-1}_{r,c})
     +\min_{r,c}\operatorname{ord}(C^{-1}_{r,c}),
\]

where zeros are ignored when taking a minimum. The valuation inequality
for a finite sum and (4) give

\[
v_*(w)\ge\eta+\min_{j,k}v_b(p_jwq_k).                       \tag{5}
\]

Allowing all physical choices in the completed word can only decrease
its minimum order. Thus the passage from fixed prefix/suffix choices to
the full \(v_b\) in (5) has the displayed direction.

For a fixed pair \((j,k)\), the map \(w\mapsto p_jwq_k\) is
injective on words of length \(t\). Multiplying (5) by the channel
weights and summing yields

\[
Q_*(t)\le
\sum_{j,k}s^{-\eta-h(p_j)-h(q_k)}
 Q_b\bigl(t+|p_j|+|q_k|\bigr).                              \tag{6}
\]

This proves the needed bounded-multiplicity comparison for grouped
channels. It does not count physical path words as independent terms.
All shifts on the right lie between zero and \(2d-2\).

In the other direction, if
\(\alpha=\min_r\operatorname{ord}(f_r)\) and
\(\gamma=\min_r\operatorname{ord}(e_r)\), then

\[
v_b(w)\ge\alpha+\gamma+v_*(w),\qquad
Q_b(t)\le s^{-\alpha-\gamma}Q_*(t).                         \tag{7}
\]

The free partition function is submultiplicative, since
\(v_*(uv)\ge v_*(u)+v_*(v)\). Therefore

\[
\beta_*=\lim_{t\to\infty}Q_*(t)^{1/t}
       =\inf_{t\ge1}Q_*(t)^{1/t}.                          \tag{8}
\]

The finite-window comparison (6) and the upper comparison (7) imply

\[
\boxed{\displaystyle
\limsup_{t\to\infty}Q_b(t)^{1/t}=\beta_* .}                 \tag{9}
\]

A fixed finite window and fixed positive multiplicative constants do not
change an exponential limsup rate. This remains true if the rate is zero.
One must use a limsup on the left: scalar boundary outputs can vanish
at every odd length, even in a minimal representation.

## 4. Regularization preserves the weighted rate

The minimal matrices are rational, and may have poles at zero. For each
channel choose an integer \(\delta_i\) such that every entry of
\(\tau^{-\delta_i}U_i(a)\) is regular. Replace

\[
U_i(a)\longmapsto\tau^{-\delta_i}U_i(a),
\qquad h_i\longmapsto h_i-\delta_i.                         \tag{10}
\]

The expression \(h(w)-v_b(w)\) is unchanged for every word. The same
holds for the free-boundary expression. Multiplying either boundary
vector by a fixed power of \(\tau\) changes \(Q_b(t)\) by a fixed
positive factor, independent of \(t\). Thus one can make all matrices
and boundaries regular without changing either rate in (9).

A channel whose reduced matrices are all zero can be deleted: every
positive-length word containing it has zero contribution in (2). If
all channels vanish, both positive-length rates are zero and no transfer
matrix is needed.

## 5. Finite rational certificates recover every strict boundary rate

Apply the free-boundary module construction at precision \(P\) to the
regularized minimal representation. Write its partition function as
\(Q_{P,*}(t)\), and its exponential rate as \(\beta_P\). The
[free-boundary theorem](cancellation-module-smith-automata.md) gives

\[
Q_{P,*}(t)\downarrow Q_*(t)\quad\text{for every fixed }t,
\qquad
\inf_{P\ge1}\beta_P=\beta_* .                              \tag{11}
\]

The second conclusion uses submultiplicativity and interchange of the
two infima over \(P\) and \(t\); it requires no uniform convergence
in the length. Combining (9) and (11) proves

\[
\boxed{\displaystyle
\inf_{P\ge1}\beta_P=
\limsup_{t\to\infty}Q_b(t)^{1/t}.}                         \tag{12}
\]

Each finite module graph has rational nonnegative weights at rational
\(s\). Every state is reachable and its free-boundary terminal weight
is one. Its rate is the spectral radius of its transition matrix
\(T_P\). For every rational \(b>\beta_P\), the positive rational
vector

\[
h=(I-T_P/b)^{-1}\mathbf1
\]

satisfies \(T_Ph=b(h-\mathbf1)<bh\). Consequently every strict
inequality above the exact right side of (12) has a finite certificate:
the reduced representation, a precision, and a positive rational
potential satisfying the usual transfer inequality.

For the two crossing cuts, carry out the reduction independently while
keeping the same evaluation point \(s\). If the product of the two
exact fixed-boundary word rates is strictly below a target, some finite
pair of module precisions and rational potentials is also below that
target. Boundary normalization contributes only fixed constants to the
coefficient bounds. The original crossing transfer applies unchanged.

This is a completeness statement for the specified wordwise sufficient
criterion. It is not a decision procedure for the best possible crossing
decomposition, and does not identify actual matrix ranks with word counts.

## 6. Why the preprocessing and limsup are needed

Consider one channel, one physical choice, exponent zero, and

\[
U=\begin{pmatrix}1&0\\0&\tau\end{pmatrix},\qquad
f^T=(1,1),\qquad e=(0,1)^T.
\]

The exact scalar output is \(\tau^t\), so its rate is \(1/s\).
In the unreduced boundary module graph, precision \(P\) gives
\(m_P(0^t)=\min(t,P)\). The resulting rate is one for every fixed
\(P\). The first coordinate remains present in the row module and
sets its normalization, although every suffix output ignores it.

The observable quotient removes that coordinate. The reduced matrix is
the one-dimensional \((\tau)\), and every module precision now gives
the exact increment one at every step. Its transfer weight is \(1/s\).

For the limsup issue, take

\[
U=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad f=e=(1,0)^T.
\]

This representation is minimal. Its scalar partition function is one
at even lengths and zero at odd lengths, whereas its free partition
function is always one. The limsup rate is one; the liminf is zero.

## 7. Reproducible exact checks

Run:

```sh
python3 tools/check_minimal_boundary_modules.py \
  --output evidence/minimal-boundary-modules.json
python3 -m unittest discover -s tests -p 'test_minimal_boundary_modules.py'
```

The checker implements exact reachable and observable reduction over
\(\mathbb F_2(\tau)\), verifies the all-length intertwining identities,
and compares literal scalar outputs independently. Eight rationally
gauged four-dimensional embeddings reduce to dimension two. These cases
check exact channel regularization and both finite-window inequalities.
The elementary, ordered-pairing, and mixed-loop crossing factorizations
provide 24 further reduction controls. The duplicate physical-choice,
parity-gap, zero-series, and invisible-coordinate examples test the
distinctions used in the theorem.

Files: [checker](../tools/check_minimal_boundary_modules.py),
[tests](../tests/test_minimal_boundary_modules.py), and
[exact evidence](../evidence/minimal-boundary-modules.json).
