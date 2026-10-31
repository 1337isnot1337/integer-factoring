# Fixed blocks strengthen the local Smith certificate

## Result and scope

Several consecutive positions of an ordered polynomial crossing can be
treated as one local block. The external alphabet grows, while the auxiliary
dimension stays fixed. Computing the block's actual Smith profile can
improve the certificate obtained by multiplying single-position profiles.

For fixed choices of cut, auxiliary side, rational gauge, and evaluation
point, both profile sums are submultiplicative in the block length. Hence
the normalized block bound converges to its infimum. This is a hierarchy
of sufficient bounds, not a claim of completeness for actual coefficient
ranks or for cancellations at particular scalar boundaries.

An exact three-position example improves the ordered-pairing control from
the single-position bound \(16\) to a normalized bound below \(16\).
It remains much weaker than the elementary scalar certificate and does
not improve the factoring exponent.

The local polynomial interface and the controls are from PR #30 at
[28105165edb8e9a0427468dc833f3c09bc0a7fb7](https://github.com/1337isnot1337/integer-factoring/tree/28105165edb8e9a0427468dc833f3c09bc0a7fb7),
especially its
[Smith-profile transfer](https://github.com/1337isnot1337/integer-factoring/blob/28105165edb8e9a0427468dc833f3c09bc0a7fb7/research/smith-profile-crossing-transfer.md).
The block construction uses that transfer unchanged.

## 1. Ordered block matrices

Let \(K(z;x,y,u,v)\) be a fixed polynomial \(D\)-by-\(D\) matrix
of degree at most \(L\), with physical alphabet size \(q\). Suppose

\[
f^Te=1,\qquad K(1;x,y,u,v)e=[x=v][y=u]e.
\]

For a block of length \(r\), define

\[
K_r(z;\mathbf x,\mathbf y,\mathbf u,\mathbf v)
=K(z;x_1,y_1,u_1,v_1)\cdots K(z;x_r,y_r,u_r,v_r).              \tag{1}
\]

Its physical alphabet has size \(q^r\), its degree is at most \(rL\),
and its auxiliary dimension is still \(D\). The endpoint condition
holds by successive application to \(e\).

In either cut, place both auxiliary indices on one physical half.
Let \(M_r(z)\) be the resulting rectangular matrix with \(q^{2r}\)
columns. Contracting its auxiliary indices by \(f_a e_b\) at \(z=1\)
gives the invertible flattening of the crossing. It therefore has full
generic column rank.

For consistent choices of cut and auxiliary side,

\[
M_{r+s}(z)=C_{r,s}\bigl(M_r(z)\otimes M_s(z)\bigr),             \tag{2}
\]

up to fixed row and column permutations. The constant matrix \(C_{r,s}\)
identifies the right auxiliary index of the first block with the left
auxiliary index of the second, and sums that common index. Both indices
are on the same physical side, so this is a left matrix multiplication
on the displayed flattening. No physical channel is summed or omitted.

## 2. Left contraction can only increase Smith orders

Work over \(O=F[[\tau]]\). Let \(M\) and \(CM\) have the same full
generic column rank \(n\), with \(C\) regular over \(O\). Write their
sorted Smith orders as \(a_1\le\cdots\le a_n\) and
\(b_1\le\cdots\le b_n\). Then

\[
b_i\ge a_i\quad(1\le i\le n).                                \tag{3}
\]

Here is a proof that gives coordinatewise inequalities, rather than
only inequalities between sums of initial segments. For every
positive integer \(\ell\), multiplication by \(C\) gives a surjection
from the image of \(M\) modulo \(\tau^\ell\) onto the image of \(CM\)
modulo \(\tau^\ell\). The first image needs exactly
\(\#\{i:a_i<\ell\}\) generators as an \(O/(\tau^\ell)\)-module;
this follows directly from its Smith decomposition. A quotient needs
no more generators, so

\[
\#\{i:b_i<\ell\}\le\#\{i:a_i<\ell\}.
\]

For all \(\ell\), these inequalities imply (3). The same statement
holds for Laurent matrices after a common scalar regularization and
restoration of the shift. If generic rank drops, append infinite
orders; that case is unnecessary for the crossing matrices here.

The Smith orders of a tensor product are the multiset of all pairwise
sums of the two order lists. This follows by tensoring their invertible
Smith basis changes.

## 3. Submultiplicative profile sums

Let \(\nu_{r,i}\) be the orders of \(M_r\) at zero. At infinity,
let \(\mu_{r,i}\) be the orders of
\(\tau^{rL}M_r(\tau^{-1})\), and put

\[
d_{r,i}=rL-\mu_{r,i}.
\]

Apply (3) to (2), both at zero and after the infinity substitution.
The sorted zero-order list for \(r+s\) is coordinatewise at least
the pairwise sums of the two zero-order lists. The sorted descending
degree list is coordinatewise at most the pairwise sums of the degree
lists.

For a fixed real \(a\ge1\), define

\[
U_r(a)=\sum_i a^{d_{r,i}},\qquad
V_r(a)=\sum_i a^{-\nu_{r,i}}.
\]

It follows that

\[
U_{r+s}(a)\le U_r(a)U_s(a),\qquad
V_{r+s}(a)\le V_r(a)V_s(a).                                  \tag{4}
\]

Use one cut at infinity and the other at zero, with each choice fixed
across block lengths. Then

\[
B_r(a)=U_r(a)V_r(a),\qquad
B_{r+s}(a)\le B_r(a)B_s(a).                                  \tag{5}
\]

In particular \(B_r(a)\le B_1(a)^r\), and subadditivity of
\(\log B_r(a)\) gives

\[
\lim_{r\to\infty}B_r(a)^{1/r}
=\inf_{r\ge1}B_r(a)^{1/r}.                                   \tag{6}
\]

The evaluation point is fixed in (4)--(6). Optimizing it separately at
each block length does not by itself prove submultiplicativity of that
new optimized sequence.

One fixed rational auxiliary gauge in each cut is also allowed. Products
of the conjugated blocks telescope, so (2) still uses constant auxiliary
contraction. Regularize each Laurent factor, apply (3), and remove its
shift. Choosing a different gauge independently at every block length is
outside the stated submultiplicativity claim.

## 4. Effective use of one fixed block

Apply the existing Smith transfer to \(K_r\) as one fixed local crossing.
At a fixed rational \(a\), its width parameter is

\[
\rho_r=\frac{\log B_r(a)}{4r\log q}.                          \tag{7}
\]

The block length, alphabet, polynomial degree, field, and auxiliary
dimension are fixed certificate data. Constructing the full input tables
at arbitrary eventual label width still takes \(2^{O(\lambda)}\) bit
operations; the fixed block changes its constant. Restricting the least
word alphabet large enough for \(2^\lambda\) labels incurs only a fixed
rounding factor.

Thus a block satisfying
\(B_r(a)<q^{\,2\sqrt2 r}\) would pass the half-rate threshold.
If the infimum in (6) is strictly below \(q^{2\sqrt2}\), some finite
block suffices. This does not assert that the infimum equals an optimal
rank rate, or that a passing block exists.

## 5. An exact three-position improvement

The elementary control is the usual scalar pencil tensored with the
auxiliary identity. Its block profiles are tensor products, and the
checker recovers equality:

\[
\min_{\text{sampled }a} B_r(a)=(256/27)^r,\qquad a=3.
\]

For the ordered-pairing control, write

\[
R=[x=y][u=v],\quad C=[x=u][y=v],\quad E=[x=v][y=u],\quad H=R+C+E
\]

over \(\mathbb F_2\). Use the fixed auxiliary matrices

\[
J_R=\begin{pmatrix}1&0\\0&0\end{pmatrix},\quad
J_H=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
J_C=\begin{pmatrix}1&0\\0&1\end{pmatrix}.
\]

The local matrix is

\[
K(z)=R(I+(1+z)J_R)+H(I+(1+z)J_H)
       +zC(I+(1+z)J_C).                                    \tag{8}
\]

It has degree two and endpoint \(K(1)=E I_2\). Both auxiliary indices
are placed on the first physical half.

| Block length | External alphabet | Best sampled \(B_r\) | Evaluation point |
| --- | ---: | ---: | ---: |
| 1 | 2 | \(16\) | \(1\) |
| 2 | 4 | \(256\) | \(1\) |
| 3 | 8 | \(1241625872629/320000000\) | \(5/4\) |

The last exact number is less than \(16^3=4096\). Its cube root is
about \(15.713765\), a strict improvement in this control's normalized
Smith certificate. The comparison is against the recorded rational
evaluation grid; no global optimization over every real point or gauge
is claimed. This remains worse than \(256/27\) for the elementary family
and does not improve the factoring bound.

## 6. Exact finite precision and reproduction

The checker uses minimum-valuation Schur elimination over truncated
series. It accepts a profile only after every column has a visible
nonzero pivot. This is an exact stopping condition, even if the chosen
precision is less than the general sufficient bound \(nL+1\).

At a pivot of order \(\theta<P\), all remaining entries have order at
least \(\theta\). Dividing the pivot by its power loses \(\theta\)
known coefficients in the unit factor. When this unit or a row multiplier
is multiplied by another pivot-row entry of order at least \(\theta\),
the lost tail starts at order \(P\). Thus every elimination step preserves
the exact residual matrix modulo \(\tau^P\). A visible pivot has its
true order, and an all-zero residual jet is never declared exact zero.
If not all pivots are visible, the algorithm increases the precision.
The degree bound \(nL+1\) guarantees termination for full generic rank.

Run:

    python3 tools/check_block_smith_crossings.py \
      --output evidence/block-smith-crossings.json

For both controls and every block through length three, the checker
independently reconstructs every Smith order from block-Toeplitz ranks.
It verifies the complete matrix endpoint, coordinatewise inequalities
in (3), exact rational versions of (4), and the displayed strict
three-position gain. A high-order diagonal entry is a negative control
against treating a truncated zero as an exact zero.

The output is [block-smith-crossings.json](../evidence/block-smith-crossings.json);
the implementation is
[check_block_smith_crossings.py](../tools/check_block_smith_crossings.py).
The finite calculations validate the profiles. The all-block
submultiplicativity theorem has the separate module-quotient proof above.

This work used substantial OpenAI assistant support. It adds ordinary
proofs and exact certificates, with no new Lean formalization or improved
factoring exponent.
